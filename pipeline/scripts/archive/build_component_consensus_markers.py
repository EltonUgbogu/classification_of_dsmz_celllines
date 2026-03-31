#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path
import pandas as pd
import numpy as np

ANCHOR_RE = re.compile(r"^anchor_(.+?)_vs_outside_component_(\d+)\.tsv$")
ISOLATE_RE = re.compile(r"^isolate_(.+?)_vs_rest\.tsv$")

def strip_ensg_version(gene_id: str) -> str:
    if pd.isna(gene_id):
        return gene_id
    return str(gene_id).split(".", 1)[0]

def read_deseq2_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # Standardize expected cols
    required = ["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")
    df = df.copy()
    df["gene_id_raw"] = df["gene_id"]
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)

    # Coerce numeric columns
    for c in ["baseMean", "log2FoldChange", "pvalue", "padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Clean obvious junk:
    # - baseMean == 0 rows are typically all-zero or non-tested
    # - drop rows with NA log2FC
    df = df[df["baseMean"].fillna(0) > 0].copy()
    df = df[~df["log2FoldChange"].isna()].copy()

    return df

def split_up_down(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    up = df[df["log2FoldChange"] > 0].copy()
    down = df[df["log2FoldChange"] < 0].copy()
    return up, down

def filter_sig(df: pd.DataFrame, padj_max: float, baseMean_min: float, lfc_min: float) -> pd.DataFrame:
    # Treat NA padj as not significant
    return df[
        (df["padj"].notna()) &
        (df["padj"] <= padj_max) &
        (df["baseMean"] >= baseMean_min) &
        (df["log2FoldChange"].abs() >= lfc_min)
    ].copy()

def aggregate_consensus(anchor_dfs: list[pd.DataFrame], direction: str) -> pd.DataFrame:
    """
    direction: 'UP' or 'DOWN'
    anchor_dfs already filtered to the correct sign + significance threshold.
    """
    if not anchor_dfs:
        return pd.DataFrame()

    # For each anchor, keep key stats
    long = []
    for i, df in enumerate(anchor_dfs, start=1):
        tmp = df[["gene_id", "baseMean", "log2FoldChange", "padj"]].copy()
        tmp["anchor_idx"] = i
        long.append(tmp)
    long = pd.concat(long, ignore_index=True)

    # Support = # anchors where gene appears
    g = long.groupby("gene_id", as_index=False).agg(
        support_n=("anchor_idx", "nunique"),
        median_log2FC=("log2FoldChange", "median"),
        median_abs_log2FC=("log2FoldChange", lambda x: float(np.median(np.abs(x)))),
        median_baseMean=("baseMean", "median"),
        best_padj=("padj", "min"),
        n_obs=("log2FoldChange", "size"),
    )

    # Rank: prefer genes supported by multiple anchors and strong effect
    g["rank_score"] = g["support_n"] * g["median_abs_log2FC"]

    # Add direction label (for downstream joins)
    g["direction"] = direction
    g = g.sort_values(["support_n", "rank_score", "best_padj"], ascending=[False, False, True]).reset_index(drop=True)
    return g

def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

def main():
    ap = argparse.ArgumentParser(description="Build component consensus markers from DESeq2 anchor tables.")
    ap.add_argument("--root", required=True, help="Pipeline results root, e.g. results/unsupervised")
    ap.add_argument("--outdir", required=True, help="Output directory for consensus results")
    ap.add_argument("--padj", type=float, default=0.05, help="padj threshold")
    ap.add_argument("--basemean", type=float, default=10.0, help="baseMean threshold (mean normalized counts)")
    ap.add_argument("--lfc", type=float, default=0.5, help="abs(log2FC) threshold")
    ap.add_argument("--core_support_min", type=int, default=2, help="Minimum anchor support for core markers")
    args = ap.parse_args()

    root = Path(args.root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Find disease dirs by existence of deseq2_markers
    diseases = []
    for d in root.iterdir():
        if d.is_dir() and (d / "deseq2_markers").exists():
            diseases.append(d.name)
    if not diseases:
        raise SystemExit(f"No disease directories found under {root} containing deseq2_markers/")

    qc_rows = []

    for disease in sorted(diseases):
        tables_dir = root / disease / "deseq2_markers" / "tables"
        if not tables_dir.exists():
            print(f"[WARN] Missing tables dir for {disease}: {tables_dir}")
            continue

        # Group anchors by component
        anchors_by_comp: dict[str, list[tuple[str, Path]]] = {}
        for f in tables_dir.glob("anchor_*_vs_outside_component_*.tsv"):
            m = ANCHOR_RE.match(f.name)
            if not m:
                continue
            anchor_name, comp = m.group(1), m.group(2)
            anchors_by_comp.setdefault(comp, []).append((anchor_name, f))

        # Process each component
        for comp, anchor_files in sorted(anchors_by_comp.items(), key=lambda x: int(x[0])):
            comp_out = outdir / disease / f"component_{comp}"
            comp_out.mkdir(parents=True, exist_ok=True)

            # Load + filter each anchor, keep UP/DOWN sets for QC overlap
            up_sets = {}
            down_sets = {}
            up_filtered = []
            down_filtered = []

            for anchor_name, path in sorted(anchor_files, key=lambda x: x[0]):
                df = read_deseq2_table(path)
                up, down = split_up_down(df)

                up_sig = filter_sig(up, args.padj, args.basemean, args.lfc)
                down_sig = filter_sig(down, args.padj, args.basemean, args.lfc)

                up_sets[anchor_name] = set(up_sig["gene_id"].tolist())
                down_sets[anchor_name] = set(down_sig["gene_id"].tolist())

                up_sig = up_sig.assign(anchor=anchor_name, component=comp, disease=disease)
                down_sig = down_sig.assign(anchor=anchor_name, component=comp, disease=disease)

                up_filtered.append(up_sig)
                down_filtered.append(down_sig)

                # Save cleaned per-anchor significant tables (optional but useful)
                up_sig.to_csv(comp_out / f"{anchor_name}.UP.sig.tsv", sep="\t", index=False)
                down_sig.to_csv(comp_out / f"{anchor_name}.DOWN.sig.tsv", sep="\t", index=False)

            # Consensus aggregation
            cons_up = aggregate_consensus(up_filtered, direction="UP")
            cons_down = aggregate_consensus(down_filtered, direction="DOWN")

            # Core vs extended
            core_up = cons_up[cons_up["support_n"] >= args.core_support_min].copy()
            core_down = cons_down[cons_down["support_n"] >= args.core_support_min].copy()
            ext_up = cons_up[cons_up["support_n"] < args.core_support_min].copy()
            ext_down = cons_down[cons_down["support_n"] < args.core_support_min].copy()

            cons_up.to_csv(comp_out / "component_consensus.UP.all.tsv", sep="\t", index=False)
            cons_down.to_csv(comp_out / "component_consensus.DOWN.all.tsv", sep="\t", index=False)
            core_up.to_csv(comp_out / "component_consensus.UP.core.tsv", sep="\t", index=False)
            core_down.to_csv(comp_out / "component_consensus.DOWN.core.tsv", sep="\t", index=False)
            ext_up.to_csv(comp_out / "component_consensus.UP.extended.tsv", sep="\t", index=False)
            ext_down.to_csv(comp_out / "component_consensus.DOWN.extended.tsv", sep="\t", index=False)

            # QC: Jaccard overlap between anchors (UP and DOWN)
            anchors = sorted(up_sets.keys())
            for i in range(len(anchors)):
                for j in range(i + 1, len(anchors)):
                    a, b = anchors[i], anchors[j]
                    qc_rows.append({
                        "disease": disease,
                        "component": comp,
                        "anchor_a": a,
                        "anchor_b": b,
                        "jaccard_up": jaccard(up_sets[a], up_sets[b]),
                        "jaccard_down": jaccard(down_sets[a], down_sets[b]),
                        "n_up_a": len(up_sets[a]),
                        "n_up_b": len(up_sets[b]),
                        "n_down_a": len(down_sets[a]),
                        "n_down_b": len(down_sets[b]),
                    })

    qc = pd.DataFrame(qc_rows)
    qc_out = outdir / "anchor_agreement_qc.tsv"
    qc.to_csv(qc_out, sep="\t", index=False)
    print(f"[OK] Wrote QC: {qc_out}")
    print(f"[OK] Consensus outputs in: {outdir}")

if __name__ == "__main__":
    main()
