#!/usr/bin/env python3
"""
build_component_consensus_markers_v2.py

Within-profile aggregation of anchor/component DESeq2 tables into signed
component consensus markers and a disease-level signed catalog.

Accepts either:
  --tables-dir  Direct path to DESeq2 component_vs_rest tables (preferred)
  --root        Legacy multi-disease layout: {root}/{disease}/deseq2_markers/tables/

Outputs per component:
  {anchor}.UP.sig.tsv / {anchor}.DOWN.sig.tsv
  component_consensus.UP.all.tsv / .core.tsv / .extended.tsv
  component_consensus.DOWN.all.tsv / .core.tsv / .extended.tsv

Global outputs:
  anchor_agreement_qc.tsv
  all_anchor_sig_long.tsv
  component_consensus_summary.tsv   (one row per component)
  disease_consensus_catalog.tsv     (signed per-gene rows across components)
"""

import argparse
import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ANCHOR_RE = re.compile(r"^anchor_(.+?)_vs_outside_component_(\d+)\.tsv$")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def strip_ensg_version(gene_id: str) -> str:
    """Remove version suffix from Ensembl gene identifiers."""
    if pd.isna(gene_id):
        return gene_id
    return str(gene_id).split(".", 1)[0]


def read_deseq2_table(path: Path) -> pd.DataFrame:
    """Load and preprocess a DESeq2 results table."""
    df = pd.read_csv(path, sep="\t")
    required = ["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")

    df = df.copy()
    df["gene_id_raw"] = df["gene_id"]
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)

    for c in ["baseMean", "log2FoldChange", "pvalue", "padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df[df["baseMean"].fillna(0) > 0].copy()
    df = df[~df["log2FoldChange"].isna()].copy()
    return df


def split_up_down(df: pd.DataFrame):
    up = df[df["log2FoldChange"] > 0].copy()
    down = df[df["log2FoldChange"] < 0].copy()
    return up, down


def filter_sig(df: pd.DataFrame, padj_max: float, baseMean_min: float,
               lfc_min: float) -> pd.DataFrame:
    return df[
        (df["padj"].notna()) &
        (df["padj"] <= padj_max) &
        (df["baseMean"] >= baseMean_min) &
        (df["log2FoldChange"].abs() >= lfc_min)
    ].copy()


def topn_fallback(df: pd.DataFrame, n: int) -> pd.DataFrame:
    df2 = df[df["padj"].notna()].copy()
    if df2.empty:
        return df2
    df2["abs_lfc"] = df2["log2FoldChange"].abs()
    df2 = df2.sort_values(
        ["padj", "abs_lfc"], ascending=[True, False]
    ).head(n).drop(columns=["abs_lfc"])
    return df2


def aggregate_consensus(anchor_dfs: list, direction: str, n_anchors: int) -> pd.DataFrame:
    """Aggregate DE results across anchors, including support_frac."""
    if not anchor_dfs:
        return pd.DataFrame()

    long = []
    for i, df in enumerate(anchor_dfs, start=1):
        tmp = df[["gene_id", "baseMean", "log2FoldChange", "padj"]].copy()
        tmp["anchor_idx"] = i
        long.append(tmp)
    long = pd.concat(long, ignore_index=True)

    g = long.groupby("gene_id", as_index=False).agg(
        support_n=("anchor_idx", "nunique"),
        median_log2FC=("log2FoldChange", "median"),
        median_abs_log2FC=("log2FoldChange", lambda x: float(np.median(np.abs(x)))),
        median_baseMean=("baseMean", "median"),
        best_padj=("padj", "min"),
        n_obs=("log2FoldChange", "size"),
    )
    g["n_anchors"] = n_anchors
    g["support_frac"] = g["support_n"] / max(n_anchors, 1)
    g["rank_score"] = g["support_n"] * g["median_abs_log2FC"]
    g["direction"] = direction
    g = g.sort_values(
        ["support_n", "rank_score", "best_padj"],
        ascending=[False, False, True]
    ).reset_index(drop=True)
    return g


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def comp_sort_key(c: str):
    try:
        return int(c)
    except ValueError:
        return 10**9


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build component consensus markers from DESeq2 anchor tables."
    )
    # Input modes
    ap.add_argument("--tables-dir",
                    help="Direct path to component_vs_rest DESeq2 tables (preferred)")
    ap.add_argument("--root",
                    help="Legacy pipeline root with {root}/{disease}/deseq2_markers/tables/")
    ap.add_argument("--disease",
                    help="Disease name for single-profile mode (used with --tables-dir)")

    # Output
    ap.add_argument("--outdir", required=True,
                    help="Output directory for consensus results")

    # Thresholds
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--basemean", type=float, default=1.0)
    ap.add_argument("--lfc", type=float, default=0.5)
    ap.add_argument("--core-support-min", type=int, default=2,
                    help="Minimum anchor count for core markers")
    ap.add_argument("--core-support-frac", type=float, default=0.30,
                    help="Minimum support fraction for core markers")
    ap.add_argument("--min-anchor-genes", type=int, default=25)
    ap.add_argument("--fallback-topn", type=int, default=200)
    ap.add_argument("--flip-policy", choices=["remove", "majority"], default="remove",
                    help="How to handle genes appearing UP in some anchors, DOWN in others")

    args = ap.parse_args()

    if not args.tables_dir and not args.root:
        ap.error("Provide either --tables-dir or --root")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Resolve input: build a dict  disease -> tables_dir_path
    # -----------------------------------------------------------------------
    disease_tables: dict[str, Path] = {}

    if args.tables_dir:
        tables_path = Path(args.tables_dir)
        if not tables_path.is_dir():
            sys.exit(f"[ERROR] --tables-dir does not exist: {tables_path}")
        disease_name = args.disease or "unknown"
        disease_tables[disease_name] = tables_path
    else:
        root = Path(args.root)
        if not root.is_dir():
            sys.exit(f"[ERROR] --root does not exist: {root}")
        for d in sorted(root.iterdir()):
            marker_dir = d / "deseq2_markers" / "tables"
            if d.is_dir() and marker_dir.exists():
                disease_tables[d.name] = marker_dir
        if not disease_tables:
            sys.exit(f"[ERROR] No disease directories with deseq2_markers/tables/ under {root}")

    # -----------------------------------------------------------------------
    # Accumulators
    # -----------------------------------------------------------------------
    qc_rows: list[dict] = []
    all_sig_long: list[pd.DataFrame] = []
    comp_summary_rows: list[dict] = []
    catalog_frames: list[pd.DataFrame] = []

    # -----------------------------------------------------------------------
    # Process each disease
    # -----------------------------------------------------------------------
    for disease, tables_dir in sorted(disease_tables.items()):
        print(f"[INFO] Processing disease={disease}  tables_dir={tables_dir}")

        # Discover anchor files grouped by component
        anchors_by_comp: dict[str, list[tuple[str, Path]]] = {}
        for f in tables_dir.glob("anchor_*_vs_outside_component_*.tsv"):
            m = ANCHOR_RE.match(f.name)
            if not m:
                continue
            anchor_name, comp = m.group(1), m.group(2)
            anchors_by_comp.setdefault(comp, []).append((anchor_name, f))

        if not anchors_by_comp:
            print(f"[WARN] No anchor files found in {tables_dir}")
            continue

        # -------------------------------------------------------------------
        # Process each component
        # -------------------------------------------------------------------
        for comp, anchor_files in sorted(anchors_by_comp.items(),
                                          key=lambda kv: comp_sort_key(kv[0])):
            comp_out = outdir / disease / f"component_{comp}"
            comp_out.mkdir(parents=True, exist_ok=True)

            n_anchors = len(anchor_files)
            up_sets: dict[str, set] = {}
            down_sets: dict[str, set] = {}
            up_filtered: list[pd.DataFrame] = []
            down_filtered: list[pd.DataFrame] = []
            # Track which anchors used fallback, per direction
            fallback_up_anchors: list[str] = []
            fallback_down_anchors: list[str] = []

            for anchor_name, path in sorted(anchor_files, key=lambda x: x[0]):
                df = read_deseq2_table(path)
                up, down = split_up_down(df)

                up_sig = filter_sig(up, args.padj, args.basemean, args.lfc)
                down_sig = filter_sig(down, args.padj, args.basemean, args.lfc)

                fb_up = False
                fb_down = False
                if len(up_sig) < args.min_anchor_genes:
                    up_sig = topn_fallback(up, args.fallback_topn)
                    fb_up = True
                    fallback_up_anchors.append(anchor_name)
                if len(down_sig) < args.min_anchor_genes:
                    down_sig = topn_fallback(down, args.fallback_topn)
                    fb_down = True
                    fallback_down_anchors.append(anchor_name)

                up_sets[anchor_name] = set(up_sig["gene_id"].tolist())
                down_sets[anchor_name] = set(down_sig["gene_id"].tolist())

                up_sig = up_sig.assign(
                    anchor=anchor_name, component=comp, disease=disease,
                    direction="UP", fallback_used=fb_up
                )
                down_sig = down_sig.assign(
                    anchor=anchor_name, component=comp, disease=disease,
                    direction="DOWN", fallback_used=fb_down
                )

                up_filtered.append(up_sig)
                down_filtered.append(down_sig)

                cols = ["disease", "component", "anchor", "direction", "gene_id",
                        "gene_id_raw", "baseMean", "log2FoldChange", "padj", "fallback_used"]
                all_sig_long.append(up_sig[cols])
                all_sig_long.append(down_sig[cols])

                up_sig.to_csv(comp_out / f"{anchor_name}.UP.sig.tsv", sep="\t", index=False)
                down_sig.to_csv(comp_out / f"{anchor_name}.DOWN.sig.tsv", sep="\t", index=False)

            # ---------------------------------------------------------------
            # Handle flip genes
            # ---------------------------------------------------------------
            all_up_genes = set().union(*up_sets.values()) if up_sets else set()
            all_down_genes = set().union(*down_sets.values()) if down_sets else set()
            flip_genes = all_up_genes & all_down_genes
            n_flip = len(flip_genes)

            if flip_genes:
                if args.flip_policy == "remove":
                    for k in up_sets:
                        up_sets[k] -= flip_genes
                    for k in down_sets:
                        down_sets[k] -= flip_genes
                    up_filtered = [df[~df["gene_id"].isin(flip_genes)].copy()
                                   for df in up_filtered]
                    down_filtered = [df[~df["gene_id"].isin(flip_genes)].copy()
                                     for df in down_filtered]
                elif args.flip_policy == "majority":
                    # Keep gene in the direction supported by more anchors
                    for gene in flip_genes:
                        n_up = sum(1 for s in up_sets.values() if gene in s)
                        n_down = sum(1 for s in down_sets.values() if gene in s)
                        if n_up >= n_down:
                            # Remove from DOWN
                            for k in down_sets:
                                down_sets[k].discard(gene)
                            down_filtered = [
                                df[~((df["gene_id"] == gene))].copy()
                                if gene in set(df["gene_id"]) else df
                                for df in down_filtered
                            ]
                        else:
                            for k in up_sets:
                                up_sets[k].discard(gene)
                            up_filtered = [
                                df[~((df["gene_id"] == gene))].copy()
                                if gene in set(df["gene_id"]) else df
                                for df in up_filtered
                            ]

            # ---------------------------------------------------------------
            # Aggregate consensus
            # ---------------------------------------------------------------
            cons_up = aggregate_consensus(up_filtered, "UP", n_anchors)
            cons_down = aggregate_consensus(down_filtered, "DOWN", n_anchors)

            # ---------------------------------------------------------------
            # Hybrid core logic: support_n >= min AND support_frac >= frac
            # ---------------------------------------------------------------
            def _is_core(row):
                return (row["support_n"] >= args.core_support_min and
                        row["support_frac"] >= args.core_support_frac)

            for cons, label in [(cons_up, "UP"), (cons_down, "DOWN")]:
                if cons.empty:
                    cons["core_status"] = pd.Series(dtype=str)
                else:
                    cons["core_status"] = cons.apply(
                        lambda r: "core" if _is_core(r) else "extended", axis=1
                    )

            core_up = cons_up[cons_up["core_status"] == "core"].copy() if not cons_up.empty else cons_up
            core_down = cons_down[cons_down["core_status"] == "core"].copy() if not cons_down.empty else cons_down
            ext_up = cons_up[cons_up["core_status"] == "extended"].copy() if not cons_up.empty else cons_up
            ext_down = cons_down[cons_down["core_status"] == "extended"].copy() if not cons_down.empty else cons_down

            cons_up.to_csv(comp_out / "component_consensus.UP.all.tsv", sep="\t", index=False)
            cons_down.to_csv(comp_out / "component_consensus.DOWN.all.tsv", sep="\t", index=False)
            core_up.to_csv(comp_out / "component_consensus.UP.core.tsv", sep="\t", index=False)
            core_down.to_csv(comp_out / "component_consensus.DOWN.core.tsv", sep="\t", index=False)
            ext_up.to_csv(comp_out / "component_consensus.UP.extended.tsv", sep="\t", index=False)
            ext_down.to_csv(comp_out / "component_consensus.DOWN.extended.tsv", sep="\t", index=False)

            # ---------------------------------------------------------------
            # Fallback support tracking
            # ---------------------------------------------------------------
            def _fallback_support(cons_df, fallback_anchors, filtered_list):
                """Count how many fallback anchors contributed to each gene."""
                if cons_df.empty or not fallback_anchors:
                    if cons_df.empty:
                        return pd.Series(dtype=int), pd.Series(dtype=float)
                    return (pd.Series(0, index=cons_df.index),
                            pd.Series(0.0, index=cons_df.index))
                fb_set = set(fallback_anchors)
                fb_genes: dict[str, int] = {}
                for df in filtered_list:
                    if df.empty:
                        continue
                    fb_rows = df[df["anchor"].isin(fb_set)]
                    for g in fb_rows["gene_id"].unique():
                        fb_genes[g] = fb_genes.get(g, 0) + 1
                fb_n = cons_df["gene_id"].map(lambda g: fb_genes.get(g, 0))
                fb_frac = fb_n / max(n_anchors, 1)
                return fb_n, fb_frac

            fb_up_n, fb_up_frac = _fallback_support(cons_up, fallback_up_anchors, up_filtered)
            fb_down_n, fb_down_frac = _fallback_support(cons_down, fallback_down_anchors, down_filtered)

            if not cons_up.empty:
                cons_up["fallback_support_n"] = fb_up_n.values
                cons_up["fallback_support_frac"] = fb_up_frac.values
            if not cons_down.empty:
                cons_down["fallback_support_n"] = fb_down_n.values
                cons_down["fallback_support_frac"] = fb_down_frac.values

            # ---------------------------------------------------------------
            # Jaccard QC
            # ---------------------------------------------------------------
            anchors = sorted(up_sets.keys())
            jacc_up_vals = []
            jacc_down_vals = []
            for i in range(len(anchors)):
                for j in range(i + 1, len(anchors)):
                    a, b = anchors[i], anchors[j]
                    j_up = jaccard(up_sets[a], up_sets[b])
                    j_down = jaccard(down_sets[a], down_sets[b])
                    jacc_up_vals.append(j_up)
                    jacc_down_vals.append(j_down)
                    qc_rows.append({
                        "disease": disease, "component": comp,
                        "anchor_a": a, "anchor_b": b,
                        "jaccard_up": j_up, "jaccard_down": j_down,
                        "n_up_a": len(up_sets[a]), "n_up_b": len(up_sets[b]),
                        "n_down_a": len(down_sets[a]), "n_down_b": len(down_sets[b]),
                        "n_flip_genes": n_flip,
                    })

            # ---------------------------------------------------------------
            # Component summary row
            # ---------------------------------------------------------------
            comp_summary_rows.append({
                "disease": disease,
                "component": comp,
                "n_anchors": n_anchors,
                "n_core_up": len(core_up),
                "n_core_down": len(core_down),
                "n_extended_up": len(ext_up),
                "n_extended_down": len(ext_down),
                "n_flip_genes": n_flip,
                "n_fallback_up_anchors": len(fallback_up_anchors),
                "n_fallback_down_anchors": len(fallback_down_anchors),
                "mean_jaccard_up": float(np.mean(jacc_up_vals)) if jacc_up_vals else float("nan"),
                "mean_jaccard_down": float(np.mean(jacc_down_vals)) if jacc_down_vals else float("nan"),
            })

            # ---------------------------------------------------------------
            # Accumulate for disease catalog
            # ---------------------------------------------------------------
            for cons_df in [cons_up, cons_down]:
                if cons_df.empty:
                    continue
                cat = cons_df.copy()
                cat["disease"] = disease
                cat["component"] = comp
                catalog_frames.append(cat)

    # -----------------------------------------------------------------------
    # Write global outputs
    # -----------------------------------------------------------------------
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(outdir / "anchor_agreement_qc.tsv", sep="\t", index=False)

    if all_sig_long:
        pd.concat(all_sig_long, ignore_index=True).to_csv(
            outdir / "all_anchor_sig_long.tsv", sep="\t", index=False
        )

    # Component summary
    pd.DataFrame(comp_summary_rows).to_csv(
        outdir / "component_consensus_summary.tsv", sep="\t", index=False
    )

    # Disease consensus catalog
    if catalog_frames:
        catalog = pd.concat(catalog_frames, ignore_index=True)
        # Ensure required columns exist with defaults for any missing
        for col, default in [("fallback_support_n", 0), ("fallback_support_frac", 0.0)]:
            if col not in catalog.columns:
                catalog[col] = default
        # Reorder columns to match spec
        col_order = [
            "disease", "component", "gene_id", "direction",
            "support_n", "support_frac", "n_anchors",
            "median_log2FC", "median_abs_log2FC", "median_baseMean",
            "best_padj", "rank_score",
            "fallback_support_n", "fallback_support_frac", "core_status",
        ]
        extra = [c for c in catalog.columns if c not in col_order]
        catalog = catalog[[c for c in col_order if c in catalog.columns] + extra]
        catalog.to_csv(outdir / "disease_consensus_catalog.tsv", sep="\t", index=False)
    else:
        # Write empty catalog with headers
        pd.DataFrame(columns=[
            "disease", "component", "gene_id", "direction",
            "support_n", "support_frac", "n_anchors",
            "median_log2FC", "median_abs_log2FC", "median_baseMean",
            "best_padj", "rank_score",
            "fallback_support_n", "fallback_support_frac", "core_status",
        ]).to_csv(outdir / "disease_consensus_catalog.tsv", sep="\t", index=False)

    print(f"[OK] Consensus outputs in: {outdir}")
    print(f"[OK] component_consensus_summary.tsv: {len(comp_summary_rows)} components")
    print(f"[OK] disease_consensus_catalog.tsv: "
          f"{sum(len(df) for df in catalog_frames)} gene-rows")


if __name__ == "__main__":
    main()
