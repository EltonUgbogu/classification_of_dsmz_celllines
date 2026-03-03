#!/usr/bin/env python3
"""
nbl_cohort_aggregator.py
NBL COHORT PROCESSING PIPELINE
Merges STAR count files → two matrices:
    1. final_nbl_count_matrix.parquet        → samples × Ensembl_ID
    2. final_nbl_count_matrix_hgnc.parquet   → samples × HGNC_Symbol
Outputs:
    • final_nbl_count_matrix.parquet
    • final_nbl_count_matrix_hgnc.parquet
    • final_nbl_sample_metadata.csv
    • final_nbl_gene_annotation.csv
    • final_nbl_ensembl_to_hgnc.tsv
    • final_nbl_gene_list_ensembl.txt
    • final_nbl_gene_list_hgnc.txt
"""

import pandas as pd
from pathlib import Path
from typing import Optional, Set, Dict, List

# ==============================
# CONFIGURATION
# ==============================
ROOT_DIR = Path("/home/chu25/data/nbl")
COHORTS = ["GSE100148", "GSE189367", "GSE49711", "SRP409177"]

# INPUT MAPS
BIOTYPE_MAP_PATH = Path("/home/chu25/data/ensembl_biotype.tsv")
HGNC_MAP_PATH    = Path("/home/chu25/data/ensembl_to_hgnc_map.tsv")  # Ensembl_ID, HGNC_Symbol

# OUTPUTS
MATRIX_ENS       = ROOT_DIR / "final_nbl_count_matrix.parquet"
MATRIX_HGNC      = ROOT_DIR / "final_nbl_count_matrix_hgnc.parquet"
META_FILE        = ROOT_DIR / "final_nbl_sample_metadata.csv"
GENE_ANNOT_CSV   = ROOT_DIR / "final_nbl_gene_annotation.csv"
ENS_TO_HGNC_TSV  = ROOT_DIR / "final_nbl_ensembl_to_hgnc.tsv"
GENE_LIST_ENS    = ROOT_DIR / "final_nbl_gene_list_ensembl.txt"
GENE_LIST_HGNC   = ROOT_DIR / "final_nbl_gene_list_hgnc.txt"

# ==============================
# LOAD FILTERS & MAPPINGS
# ==============================
def load_protein_coding_ids() -> Optional[Set[str]]:
    if not BIOTYPE_MAP_PATH.exists():
        print("Biotype file missing – skipping protein-coding filter.")
        return None
    print(f"Loading biotype map from {BIOTYPE_MAP_PATH}")
    biotype = pd.read_csv(BIOTYPE_MAP_PATH, sep="\t")
    pc = set(biotype.loc[biotype["gene_biotype"] == "protein_coding", "Ensembl_ID"])
    print(f"   → {len(pc):,} protein-coding genes")
    return pc

def load_hgnc_map() -> Dict[str, str]:
    if not HGNC_MAP_PATH.exists():
        raise FileNotFoundError(f"HGNC map not found: {HGNC_MAP_PATH}")
    print(f"Loading Ensembl → HGNC map from {HGNC_MAP_PATH}")
    df = pd.read_csv(HGNC_MAP_PATH, sep="\t", dtype=str)
    if not {"Ensembl_ID", "HGNC_Symbol"}.issubset(df.columns):
        raise ValueError("HGNC map must have columns: Ensembl_ID, HGNC_Symbol")
    df["Ensembl_ID"] = df["Ensembl_ID"].str.replace(r"\.\d+$", "", regex=True)
    hgnc_dict = dict(zip(df["Ensembl_ID"], df["HGNC_Symbol"].fillna("")))
    print(f"   → {len(hgnc_dict):,} mappings loaded")
    return hgnc_dict

# ==============================
# STRANDEDNESS DETECTION
# ==============================
def detect_stranded_col(df_raw: pd.DataFrame) -> int:
    df = df_raw[~df_raw.iloc[:, 0].astype(str).str.startswith("N_", na=False)].copy()
    if df.empty or df.shape[1] < 4:
        return 2
    s2 = pd.to_numeric(df.iloc[:, 1], errors="coerce").fillna(0).sum()
    s3 = pd.to_numeric(df.iloc[:, 2], errors="coerce").fillna(0).sum()
    s4 = pd.to_numeric(df.iloc[:, 3], errors="coerce").fillna(0).sum()
    if s3 >= 1.5 * s4 and s3 >= 0.6 * s2:
        return 3
    if s4 >= 1.5 * s3 and s4 >= 0.6 * s2:
        return 4
    return 2

# ==============================
# LOAD ONE STAR .tab FILE
# ==============================
def load_star_counts_tsv(
    path: Path,
    sample_name: str,
    stranded_col: Optional[int] = None,
    protein_coding_ids: Optional[Set[str]] = None,
) -> pd.Series:
    df = pd.read_csv(
        path, sep="\t", header=None, engine="python",
        dtype={0: "string"}, na_filter=False
    )
    if df.shape[1] < 4:
        raise ValueError(f"{path} has <4 columns")
    col = stranded_col if stranded_col is not None else detect_stranded_col(df)
    strand_label = {2: "unstranded", 3: "forward/+", 4: "reverse/-"}.get(col, "unknown")
    print(f"   [{sample_name}] Using column {col} ({strand_label})")
    sub = df.iloc[:, [0, col - 1]].copy()
    sub.columns = ["gene_id_ver", sample_name]
    sub["gene_id_ver"] = sub["gene_id_ver"].str.strip()
    sub = sub[~sub["gene_id_ver"].str.startswith("N_")]
    sub[sample_name] = pd.to_numeric(sub[sample_name], errors="coerce").fillna(0).astype("int32")
    sub["Ensembl_ID"] = sub["gene_id_ver"].str.replace(r"\.\d+$", "", regex=True)
    if protein_coding_ids is not None:
        before = len(sub)
        sub = sub[sub["Ensembl_ID"].isin(protein_coding_ids)]
        print(f"   → protein-coding: {before} → {len(sub)} genes")
    return sub[["Ensembl_ID", sample_name]].set_index("Ensembl_ID")[sample_name], strand_label

# ==============================
# MAIN
# ==============================
def main():
    pc_ids = load_protein_coding_ids()
    hgnc_map = load_hgnc_map()

    cohort_matrices: List[pd.DataFrame] = []
    meta_rows: List[Dict[str, str]] = []
    all_ensembl_ids: Set[str] = set()

    print("Starting NBL cohort aggregation...\n")

    for cohort in COHORTS:
        cohort_path = ROOT_DIR / cohort / "results" / "valid_count"
        tab_files = sorted(cohort_path.glob("*.tab"))
        if not tab_files:
            print(f"WARNING: No .tab files in {cohort_path}")
            continue
        print(f"Processing cohort: {cohort} ({len(tab_files)} samples)")

        series_list: List[pd.Series] = []
        for fp in tab_files:
            sample_name = f"{cohort}_{fp.stem}"
            try:
                counts_series, strand = load_star_counts_tsv(
                    path=fp,
                    sample_name=sample_name,
                    stranded_col=None,
                    protein_coding_ids=pc_ids,
                )
                series_list.append(counts_series)
                all_ensembl_ids.update(counts_series.index)
                meta_rows.append({
                    "sample": sample_name,
                    "cohort": cohort,
                    "file": str(fp),
                    "strand_info": strand
                })
            except Exception as e:
                print(f"   ERROR {fp.name}: {e}")
                continue

        if not series_list:
            print(f"   No valid samples for {cohort}. Skipping.")
            continue

        cohort_df = pd.concat(series_list, axis=1).fillna(0).astype("int32")
        print(f"   Cohort matrix: {cohort_df.shape[1]} samples × {cohort_df.shape[0]} genes")
        cohort_matrices.append(cohort_df)

    if not cohort_matrices:
        raise RuntimeError("No data processed.")

    print("\nMerging all cohorts...")
    final_ens = pd.concat(cohort_matrices, axis=1, join="outer").fillna(0).astype("int32")
    final_ens = final_ens.loc[~final_ens.index.duplicated(keep="first")].sort_index()

    # -------------------------------
    # 1. Save Ensembl matrix
    # -------------------------------
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(META_FILE, index=False)
    print(f"Metadata saved → {META_FILE}")

    final_t_ens = final_ens.T
    final_t_ens = final_t_ens.sort_index(axis=1)
    final_t_ens = final_t_ens.reindex(meta_df["sample"])

    print(f"\nSaving Ensembl matrix → {MATRIX_ENS}")
    final_t_ens.to_parquet(MATRIX_ENS, index=True, engine="pyarrow", compression="zstd")
    print(f"   → {final_t_ens.shape[0]} samples × {final_t_ens.shape[1]} Ensembl genes")

    # -------------------------------
    # 2. Map to HGNC
    # -------------------------------
    print("\nMapping Ensembl → HGNC...")
    hgnc_series = final_ens.index.to_series().map(hgnc_map)
    hgnc_series = hgnc_series.fillna(final_ens.index.to_series())  # fallback to Ensembl

    final_hgnc = final_ens.copy()
    final_hgnc.index = hgnc_series
    final_hgnc = final_hgnc.groupby(level=0).sum().astype("int32")
    final_hgnc = final_hgnc.loc[~final_hgnc.index.duplicated()].sort_index()

    final_t_hgnc = final_hgnc.T
    final_t_hgnc = final_t_hgnc.sort_index(axis=1)
    final_t_hgnc = final_t_hgnc.reindex(meta_df["sample"])

    print(f"Saving HGNC matrix → {MATRIX_HGNC}")
    final_t_hgnc.to_parquet(MATRIX_HGNC, index=True, engine="pyarrow", compression="zstd")
    print(f"   → {final_t_hgnc.shape[0]} samples × {final_t_hgnc.shape[1]} HGNC genes")

    # -------------------------------
    # 3. Save annotations
    # -------------------------------
    print(f"\nSaving gene annotation → {GENE_ANNOT_CSV}")
    gene_annot = pd.DataFrame({
        "Ensembl_ID": final_ens.index,
        "HGNC_Symbol": hgnc_series.values
    })
    gene_annot.to_csv(GENE_ANNOT_CSV, index=False)

    print(f"Saving Ensembl→HGNC map → {ENS_TO_HGNC_TSV}")
    gene_annot.to_csv(ENS_TO_HGNC_TSV, sep="\t", index=False)

    print(f"Saving gene list (Ensembl) → {GENE_LIST_ENS}")
    final_t_ens.columns.to_series().to_csv(GENE_LIST_ENS, index=False, header=False)

    print(f"Saving gene list (HGNC) → {GENE_LIST_HGNC}")
    final_t_hgnc.columns.to_series().to_csv(GENE_LIST_HGNC, index=False, header=False)

    print("\nPIPELINE COMPLETE!")
    print("   • Two matrices: Ensembl + HGNC")
    print("   • Sample names preserved (index=True)")
    print("   • All annotations saved")
    print("   • Ready for purity, DE, survival, WGCNA")

if __name__ == "__main__":
    main()