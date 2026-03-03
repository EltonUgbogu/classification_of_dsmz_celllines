#!/usr/bin/env python3
"""
hema_cohort_aggregator.py

Aggregate STAR ReadsPerGene.out.tab files for the HEMA (CLL/BL/MPN) cohort.

Outputs:
- final_hema_count_matrix.parquet        (samples × Ensembl_ID)
- final_hema_count_matrix_hgnc.parquet   (samples × HGNC_Symbol)
- final_hema_sample_metadata.csv
- final_hema_gene_annotation.csv
- final_hema_ensembl_to_hgnc.tsv
- final_hema_gene_list_ensembl.txt
- final_hema_gene_list_hgnc.txt
"""

import argparse
from pathlib import Path
from typing import List, Dict, Optional, Set

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate STAR count files into matrices.")
    p.add_argument(
        "--counts-dirs",
        nargs="+",
        required=True,
        help="Directories containing STAR valid_count/*.tab files.",
    )
    p.add_argument(
        "--gtf",
        required=True,
        help="GTF annotation file (Ensembl gene_ids etc.).",
    )
    p.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix, e.g. /home/.../final_hema",
    )
    p.add_argument(
        "--strand-col",
        type=int,
        choices=[2, 3, 4],
        default=None,
        help="STAR ReadsPerGene column (2=unstranded, 3=forward, 4=reverse). "
             "If omitted, auto-detect based on non-zero counts.",
    )
    return p.parse_args()


def list_count_files(counts_dirs: List[str]) -> Dict[str, Dict]:
    """Returns dict of sample -> {path, project}"""
    files: Dict[str, Dict] = {}
    for d in counts_dirs:
        dpath = Path(d)
        if not dpath.exists():
            print(f"[WARN] Counts dir does not exist: {d}")
            continue
        # Extract project ID from path (e.g., .../cll/GSE237907/results/star/valid_count -> GSE237907)
        project = None
        for part in dpath.parts:
            if part.startswith("GSE"):
                project = part
                break
        if project is None:
            project = "unknown"
        for f in sorted(dpath.glob("*.tab")):
            sample = f.stem
            if sample in files:
                print(f"[WARN] Duplicate sample {sample}: {f} (keeping first).")
                continue
            files[sample] = {"path": f, "project": project}
    if not files:
        raise SystemExit("No *.tab files found in any --counts-dirs")
    print(f"[INFO] Found {len(files)} STAR count files.")
    return files


def detect_stranded_col(df: pd.DataFrame) -> int:
    mask = ~df.iloc[:, 0].astype(str).str.startswith("N_")
    candidates = {}
    for col in [2, 3, 4]:
        if col >= df.shape[1]:
            continue
        counts = df.loc[mask, col]
        candidates[col] = counts.sum()
    if not candidates:
        raise ValueError("No suitable count columns (2,3,4) found in STAR counts.")
    best_col = max(candidates, key=candidates.get)
    print(f"[INFO] Auto-detected stranded column: {best_col}")
    return best_col


def load_star_counts_tsv(
    path: Path,
    sample_name: str,
    stranded_col: Optional[int] = None,
) -> pd.Series:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        engine="python",
        dtype={0: "string"},
        na_filter=False,
    )
    if df.shape[1] < 4:
        raise ValueError(f"{path} has <4 columns; expected STAR ReadsPerGene format.")

    col = stranded_col if stranded_col is not None else detect_stranded_col(df)
    label = {2: "unstranded", 3: "forward/+", 4: "reverse/-"}.get(col, "unknown")
    print(f"   [{sample_name}] using column {col} ({label}) from {path.name}")

    mask = ~df.iloc[:, 0].astype(str).str.startswith("N_")
    gene_ids = df.loc[mask, 0].astype(str)
    counts = df.loc[mask, col].astype("int64")
    s = pd.Series(counts.values, index=gene_ids.values, name=sample_name)
    return s


def build_count_matrix(files: Dict[str, Dict], strand_col: Optional[int]) -> pd.DataFrame:
    series_list = []
    for sample, info in files.items():
        path = info["path"]
        print(f"[LOAD] {sample} → {path}")
        s = load_star_counts_tsv(path, sample, stranded_col=strand_col)
        series_list.append(s)
    df = pd.concat(series_list, axis=1).T  # samples × genes
    df.index.name = "sample"
    print(f"[INFO] Count matrix shape: {df.shape[0]} samples × {df.shape[1]} genes")
    return df


def parse_gtf(gtf_path: Path, gene_ids: Set[str]) -> pd.DataFrame:
    rows = []
    print(f"[INFO] Parsing GTF for gene annotations: {gtf_path}")
    with gtf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attr = parts
            if feature != "gene":
                continue

            attrs = {}
            for field in attr.split(";"):
                field = field.strip()
                if not field or " " not in field:
                    continue
                key, value = field.split(" ", 1)
                attrs[key] = value.strip().strip('"')

            gid = attrs.get("gene_id")
            if gid not in gene_ids:
                continue

            row = {
                "Ensembl_ID": gid,
                "HGNC_Symbol": attrs.get("gene_name", ""),
                "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                "seqname": seqname,
                "strand": strand,
            }
            rows.append(row)

    annot = pd.DataFrame.from_records(rows).drop_duplicates("Ensembl_ID")
    print(f"[INFO] Parsed annotation rows: {annot.shape[0]} for {len(gene_ids)} genes")
    return annot


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_dirs = [str(Path(d)) for d in args.counts_dirs]
    files = list_count_files(counts_dirs)

    print("[STEP] Building Ensembl count matrix…")
    ens_mat = build_count_matrix(files, strand_col=args.strand_col)

    # Build metadata with sample and project columns
    meta_records = []
    for sample in ens_mat.index:
        project = files[sample]["project"] if sample in files else "unknown"
        meta_records.append({"sample": sample, "project": project})
    meta = pd.DataFrame.from_records(meta_records)
    meta.to_csv(out_dir / "final_hema_sample_metadata.csv", index=False)
    print(f"[INFO] Metadata with project column: {meta['project'].value_counts().to_dict()}")

    gene_ids = set(ens_mat.columns.astype(str))
    gene_annot = parse_gtf(Path(args.gtf), gene_ids)
    gene_annot.to_csv(out_dir / "final_hema_gene_annotation.csv", index=False)

    ens_path = out_dir / "final_hema_count_matrix.parquet"
    print(f"[SAVE] Ensembl matrix → {ens_path}")
    ens_mat.to_parquet(ens_path)

    if "HGNC_Symbol" in gene_annot.columns:
        mapping = gene_annot[["Ensembl_ID", "HGNC_Symbol"]].drop_duplicates()
        mapping = mapping[mapping["HGNC_Symbol"] != ""]
        mapping.to_csv(
            out_dir / "final_hema_ensembl_to_hgnc.tsv",
            sep="\t",
            index=False,
        )

        m = mapping.set_index("Ensembl_ID")["HGNC_Symbol"]
        common = ens_mat.columns.intersection(m.index)
        print(f"[INFO] Genes with HGNC symbol: {len(common)}")
        mat_sub = ens_mat[common]
        mat_sub.columns = m.loc[common].values
        hgnc_mat = mat_sub.groupby(axis=1, level=0).sum()

        hgnc_path = out_dir / "final_hema_count_matrix_hgnc.parquet"
        print(f"[SAVE] HGNC matrix → {hgnc_path}")
        hgnc_mat.to_parquet(hgnc_path)

        ens_mat.columns.to_series().to_csv(
            out_dir / "final_hema_gene_list_ensembl.txt",
            index=False,
            header=False,
        )
        hgnc_mat.columns.to_series().to_csv(
            out_dir / "final_hema_gene_list_hgnc.txt",
            index=False,
            header=False,
        )
    else:
        print("[WARN] HGNC_Symbol not found; skipping HGNC matrix.")
        ens_mat.columns.to_series().to_csv(
            out_dir / "final_hema_gene_list_ensembl.txt",
            index=False,
            header=False,
        )

    print("\n[DONE] HEMA cohort aggregation complete.")


if __name__ == "__main__":
    main()
