#!/usr/bin/env python3
"""
build_component_consensus_markers_v2.py

Within-profile aggregation of anchor/component DESeq2 tables into signed
component consensus markers plus a summary export layer that feeds
build_pan_cancer_features.py.

Primary CLI (used by Snakemake):
  --tables-dir   Path to component_vs_rest or isolate-vs-rest tables for a single profile
  --profile-name Logical profile name (e.g. brca, nbl, rbl)
  --outdir       Output directory for consensus results

Legacy multi-disease layout via --root is retained for backwards
compatibility but deprecated.

Outputs per component:
  {anchor}.UP.sig.tsv / {anchor}.DOWN.sig.tsv
  component_consensus.UP.all/core/extended.tsv
  component_consensus.DOWN.all/core/extended.tsv

Profile-level outputs under {outdir}/:
  anchor_agreement_qc.tsv
  all_anchor_sig_long.tsv
  component_consensus_summary.tsv
  disease_consensus_catalog.tsv (legacy contract)

New summary exports under {outdir}/summary/:
  consensus_feature_catalog.tsv
  gene_recurrence_across_components.tsv
  unique_feature_set_recurrence_ge_2.txt (+ UP/DOWN splits)
  consensus_export_done.txt
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
ISOLATE_RE = re.compile(r"^isolate_(.+?)_vs_rest\.tsv$")


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
    if not anchor_dfs or n_anchors <= 0:
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
    g["n_anchors_in_component"] = n_anchors
    g["support_frac"] = g["support_n"] / n_anchors
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
        return (0, int(c))
    except ValueError:
        return (1, str(c))


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def discover_contrast_tables(tables_dir: Path) -> tuple[dict, str]:
    """Discover supported DESeq2 table layouts."""
    anchors_by_comp = {}

    for f in sorted(tables_dir.glob("anchor_*_vs_outside_component_*.tsv")):
        m = ANCHOR_RE.match(f.name)
        if not m:
            continue
        anchor_name, comp = m.group(1), m.group(2)
        anchors_by_comp.setdefault(comp, []).append((anchor_name, f))

    if anchors_by_comp:
        return anchors_by_comp, "anchor_vs_outside_component"

    for f in sorted(tables_dir.glob("isolate_*_vs_rest.tsv")):
        m = ISOLATE_RE.match(f.name)
        if not m:
            continue
        isolate_name = m.group(1)
        anchors_by_comp.setdefault(isolate_name, []).append((isolate_name, f))

    return anchors_by_comp, "isolate_vs_rest"


def write_summary_exports(profile: str,
                          catalog: pd.DataFrame,
                          summary_dir: Path,
                          n_components: int,
                          recurrence_k: int) -> None:
    """Write summary-layer exports consumed by downstream scripts."""
    if catalog.empty:
        raise SystemExit(f"[ERROR] No consensus genes for profile={profile}; cannot export summary")

    ensure_dir(summary_dir)

    feature_df = catalog.copy()
    feature_df["profile"] = profile
    feature_df["marker_class"] = feature_df["core_status"].map(
        lambda x: "core" if str(x).lower() == "core" else "extended"
    )

    feature_columns = [
        "profile", "component", "gene_id", "direction", "marker_class",
        "support_n", "support_frac", "median_log2FC", "median_abs_log2FC",
        "median_baseMean", "best_padj", "rank_score",
    ]
    # Ensure all required columns exist even if empty
    for col in feature_columns:
        if col not in feature_df.columns:
            feature_df[col] = []
    feature_df = feature_df[feature_columns]
    feature_catalog_path = summary_dir / "consensus_feature_catalog.tsv"
    feature_df.to_csv(feature_catalog_path, sep="\t", index=False)

    recurrence = feature_df.groupby(["gene_id", "direction"], as_index=False).agg(
        freq=("component", "nunique"),
        mean_support_n=("support_n", "mean"),
        mean_support_frac=("support_frac", "mean"),
        mean_abs_log2FC=("median_abs_log2FC", "mean"),
        max_rank_score=("rank_score", "max"),
    )
    recurrence["profile"] = profile
    recurrence["n_components"] = n_components
    recurrence = recurrence[[
        "profile", "gene_id", "direction", "freq", "n_components",
        "mean_support_n", "mean_support_frac", "mean_abs_log2FC",
        "max_rank_score"
    ]]
    recurrence_path = summary_dir / "gene_recurrence_across_components.tsv"
    recurrence.to_csv(recurrence_path, sep="\t", index=False)

    threshold = max(int(recurrence_k), 1)
    genes_any = sorted({
        row.gene_id for row in recurrence.itertuples()
        if row.freq >= threshold
    })
    genes_up = sorted({
        row.gene_id for row in recurrence.itertuples()
        if row.freq >= threshold and row.direction == "UP"
    })
    genes_down = sorted({
        row.gene_id for row in recurrence.itertuples()
        if row.freq >= threshold and row.direction == "DOWN"
    })

    def _write_gene_list(path: Path, genes) -> None:
        with open(path, "w") as handle:
            for gene in genes:
                handle.write(f"{gene}\n")

    _write_gene_list(summary_dir / "unique_feature_set_recurrence_ge_2.txt", genes_any)
    _write_gene_list(summary_dir / "unique_feature_set_recurrence_ge_2.UP.txt", genes_up)
    _write_gene_list(summary_dir / "unique_feature_set_recurrence_ge_2.DOWN.txt", genes_down)

    with open(summary_dir / "consensus_export_done.txt", "w") as handle:
        handle.write(f"profile={profile}\n")


def process_profile(profile: str, tables_dir: Path, outdir: Path, args: argparse.Namespace) -> None:
    profile = str(profile)
    tables_dir = Path(tables_dir)
    outdir = Path(outdir)
    ensure_dir(outdir)

    print(f"[INFO] Processing profile={profile} tables_dir={tables_dir}")

    anchors_by_comp, table_mode = discover_contrast_tables(tables_dir)

    if not anchors_by_comp:
        raise SystemExit(
            f"[ERROR] No supported DESeq2 tables found for profile={profile} in {tables_dir}. "
            "Expected anchor_*_vs_outside_component_*.tsv or isolate_*_vs_rest.tsv"
        )

    print(f"[INFO] Detected DESeq2 table layout: {table_mode}")

    qc_rows = []
    all_sig_long = []
    comp_summary_rows = []
    catalog_frames = []
    processed_components = 0

    summary_dir = outdir / "summary"

    for comp, anchor_files in sorted(anchors_by_comp.items(), key=lambda kv: comp_sort_key(kv[0])):
        processed_components += 1
        comp_out = ensure_dir(outdir / f"component_{comp}")

        n_anchors = len(anchor_files)
        up_sets = {}
        down_sets = {}
        up_filtered = []
        down_filtered = []
        fallback_up_anchors = []
        fallback_down_anchors = []

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
                profile=profile,
                disease=profile,
                component=comp,
                direction="UP",
                anchor=anchor_name,
                fallback_used=fb_up,
            )
            down_sig = down_sig.assign(
                profile=profile,
                disease=profile,
                component=comp,
                direction="DOWN",
                anchor=anchor_name,
                fallback_used=fb_down,
            )

            up_filtered.append(up_sig)
            down_filtered.append(down_sig)

            cols = [
                "profile", "disease", "component", "anchor", "direction",
                "gene_id", "gene_id_raw", "baseMean", "log2FoldChange",
                "padj", "fallback_used"
            ]
            all_sig_long.append(up_sig[cols])
            all_sig_long.append(down_sig[cols])

            up_sig.to_csv(comp_out / f"{anchor_name}.UP.sig.tsv", sep="\t", index=False)
            down_sig.to_csv(comp_out / f"{anchor_name}.DOWN.sig.tsv", sep="\t", index=False)

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
                up_filtered = [df[~df["gene_id"].isin(flip_genes)].copy() for df in up_filtered]
                down_filtered = [df[~df["gene_id"].isin(flip_genes)].copy() for df in down_filtered]
            elif args.flip_policy == "majority":
                # Keep gene in whichever direction has majority support
                for gene in list(flip_genes):
                    up_count = sum(gene in up_sets[a] for a in up_sets)
                    down_count = sum(gene in down_sets[a] for a in down_sets)
                    if up_count >= down_count:
                        for k in down_sets:
                            down_sets[k].discard(gene)
                        down_filtered = [df[df["gene_id"] != gene].copy() for df in down_filtered]
                    else:
                        for k in up_sets:
                            up_sets[k].discard(gene)
                        up_filtered = [df[df["gene_id"] != gene].copy() for df in up_filtered]

        cons_up = aggregate_consensus(up_filtered, "UP", n_anchors)
        cons_down = aggregate_consensus(down_filtered, "DOWN", n_anchors)

        def _is_core(row):
            return (
                row["support_n"] >= args.core_support_min and
                row["support_frac"] >= args.core_support_frac
            )

        for cons in [cons_up, cons_down]:
            if not cons.empty:
                cons["core_status"] = cons.apply(lambda r: "core" if _is_core(r) else "extended", axis=1)
            else:
                cons["core_status"] = pd.Series(dtype=str)

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

        def _fallback_support(cons_df, fallback_anchors, filtered_list):
            if cons_df.empty or not fallback_anchors:
                if cons_df.empty:
                    return pd.Series(dtype=int), pd.Series(dtype=float)
                return pd.Series(0, index=cons_df.index), pd.Series(0.0, index=cons_df.index)
            fb_set = set(fallback_anchors)
            fb_genes = {}
            for df in filtered_list:
                if df.empty:
                    continue
                fb_rows = df[df["anchor"].isin(fb_set)]
                for g in fb_rows["gene_id"].unique():
                    fb_genes[g] = fb_genes.get(g, 0) + 1
            fb_n = cons_df["gene_id"].map(lambda g: fb_genes.get(g, 0))
            fb_frac = fb_n / n_anchors
            return fb_n, fb_frac

        fb_up_n, fb_up_frac = _fallback_support(cons_up, fallback_up_anchors, up_filtered)
        fb_down_n, fb_down_frac = _fallback_support(cons_down, fallback_down_anchors, down_filtered)

        if not cons_up.empty:
            cons_up["fallback_support_n"] = fb_up_n.values
            cons_up["fallback_support_frac"] = fb_up_frac.values
        if not cons_down.empty:
            cons_down["fallback_support_n"] = fb_down_n.values
            cons_down["fallback_support_frac"] = fb_down_frac.values

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
                    "profile": profile,
                    "component": comp,
                    "anchor_a": a,
                    "anchor_b": b,
                    "jaccard_up": j_up,
                    "jaccard_down": j_down,
                    "n_up_a": len(up_sets[a]),
                    "n_up_b": len(up_sets[b]),
                    "n_down_a": len(down_sets[a]),
                    "n_down_b": len(down_sets[b]),
                    "n_flip_genes": n_flip,
                })

        comp_summary_rows.append({
            "profile": profile,
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

        for cons_df in [cons_up, cons_down]:
            if cons_df.empty:
                continue
            cat = cons_df.copy()
            cat["profile"] = profile
            cat["disease"] = profile
            cat["component"] = comp
            catalog_frames.append(cat)

    if not catalog_frames:
        raise SystemExit(f"[ERROR] No valid consensus genes for profile={profile}")

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(outdir / "anchor_agreement_qc.tsv", sep="\t", index=False)

    if all_sig_long:
        pd.concat(all_sig_long, ignore_index=True).to_csv(
            outdir / "all_anchor_sig_long.tsv", sep="\t", index=False
        )

    pd.DataFrame(comp_summary_rows).to_csv(
        outdir / "component_consensus_summary.tsv", sep="\t", index=False
    )

    catalog = pd.concat(catalog_frames, ignore_index=True)
    for col, default in [("fallback_support_n", 0), ("fallback_support_frac", 0.0)]:
        if col not in catalog.columns:
            catalog[col] = default
    col_order = [
        "profile", "disease", "component", "gene_id", "direction",
        "support_n", "support_frac", "n_anchors_in_component",
        "median_log2FC", "median_abs_log2FC", "median_baseMean",
        "best_padj", "rank_score",
        "fallback_support_n", "fallback_support_frac", "core_status",
    ]
    extra = [c for c in catalog.columns if c not in col_order]
    catalog = catalog[[c for c in col_order if c in catalog.columns] + extra]
    catalog_path = outdir / "disease_consensus_catalog.tsv"
    catalog.to_csv(catalog_path, sep="\t", index=False)

    write_summary_exports(profile, catalog, summary_dir, processed_components, args.recurrence_k)

    print(f"[OK] Consensus outputs for profile={profile} in {outdir}")
    print(f"[OK] component_consensus_summary.tsv: {len(comp_summary_rows)} components")
    print(f"[OK] disease_consensus_catalog.tsv: {len(catalog)} gene-rows")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build component consensus markers from DESeq2 anchor tables."
    )
    ap.add_argument("--tables-dir",
                    help="Path to component_vs_rest tables for the active profile")
    ap.add_argument("--profile-name",
                    help="Logical profile/disease name (e.g. brca, nbl)")
    ap.add_argument("--root",
                    help="[Deprecated] Legacy root with {root}/{disease}/deseq2_markers/tables/")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for profile-specific consensus results")
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
    ap.add_argument("--recurrence-k", type=int, default=2,
                    help="Minimum component count for recurrence gene summaries")

    args = ap.parse_args()

    if not args.tables_dir and not args.root:
        ap.error("Provide --tables-dir for per-profile mode or --root for legacy multi-profile mode")

    if args.tables_dir and not args.profile_name:
        ap.error("--profile-name is required when --tables-dir is provided")

    jobs = []

    if args.tables_dir:
        tables_path = Path(args.tables_dir)
        if not tables_path.is_dir():
            sys.exit(f"[ERROR] --tables-dir does not exist: {tables_path}")
        jobs.append((args.profile_name, tables_path, Path(args.outdir)))
    elif args.root:
        root = Path(args.root)
        if not root.is_dir():
            sys.exit(f"[ERROR] --root does not exist: {root}")
        for d in sorted(root.iterdir()):
            marker_dir = d / "deseq2_markers" / "tables"
            if d.is_dir() and marker_dir.exists():
                jobs.append((d.name, marker_dir, Path(args.outdir) / d.name))
        if not jobs:
            sys.exit(f"[ERROR] No profile directories with deseq2_markers/tables/ under {root}")

    for profile, tables_dir, profile_outdir in jobs:
        process_profile(profile, tables_dir, profile_outdir, args)


if __name__ == "__main__":
    main()
