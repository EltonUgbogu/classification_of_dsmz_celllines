#!/usr/bin/env python3
"""
build_pan_cancer_features.py

Cross-profile merge step: reads disease-level signed consensus catalogs
(produced by build_component_consensus_markers_v2.py) and combines them
with optional isolate markers into a pan-cancer feature panel.

Selection hierarchy:
  1. PAN_CORE      -- genes in >= min_cross_disease diseases
  2. DISEASE_SPECIFIC -- per-disease additions (capped)
  3. ISOLATE_SPECIFIC -- per-disease isolate additions (capped)

Outputs:
  pan_cancer_features.tsv          -- full provenance table
  pan_cancer_features.UP.txt       -- UP gene IDs only
  pan_cancer_features.DOWN.txt     -- DOWN gene IDs only
  pan_cancer_features_clean.txt    -- all gene IDs (no comments)
  pan_cancer_features_summary.tsv  -- summary stats
"""

import argparse
import os
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd


# ---------------------------------------------------------------------------
# Ensembl ID normalisation (same logic as consensus script)
# ---------------------------------------------------------------------------

def strip_ensg_version(gene_id: str) -> str:
    if pd.isna(gene_id):
        return gene_id
    return str(gene_id).split(".", 1)[0]


def normalize_gene_ids(df: pd.DataFrame, col: str = "gene_id") -> pd.DataFrame:
    """Strip Ensembl version from a gene_id column."""
    df = df.copy()
    df[col] = df[col].map(strip_ensg_version)
    return df


# ---------------------------------------------------------------------------
# Gene annotation loading for ribo/MT filtering
# ---------------------------------------------------------------------------

def load_gene_annotation(path: str) -> dict:
    """Load a TSV with at least gene_id and gene_type/gene_name columns.

    Returns dict mapping gene_id -> dict of annotation fields.
    """
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        print(f"[WARN] gene_annotation_tsv missing 'gene_id' column: {path}",
              file=sys.stderr)
        return {}
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)
    records = {}
    for _, row in df.iterrows():
        records[row["gene_id"]] = row.to_dict()
    return records


def is_ribo_or_mt(gene_id: str, annotation: dict) -> bool:
    """Check if gene is ribosomal or mitochondrial using annotation."""
    if not annotation:
        return False
    info = annotation.get(gene_id, {})
    gene_name = str(info.get("gene_name", "")).upper()
    gene_type = str(info.get("gene_type", info.get("gene_biotype", ""))).lower()

    # Ribosomal protein genes
    if gene_name.startswith(("RPL", "RPS", "MRPL", "MRPS")):
        return True
    # Mitochondrial genes
    if gene_name.startswith(("MT-", "MT_")):
        return True
    # rRNA biotype
    if "rrna" in gene_type:
        return True
    return False


# ---------------------------------------------------------------------------
# Catalog loading
# ---------------------------------------------------------------------------

def load_catalog(path: str) -> pd.DataFrame:
    """Load a disease_consensus_catalog.tsv."""
    df = pd.read_csv(path, sep="\t")
    required = ["gene_id", "direction", "component"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] Catalog {path} missing columns: {missing}")
    return normalize_gene_ids(df)


def load_isolate_dir(dirpath: str) -> set:
    """Load gene IDs from isolate marker files in a directory."""
    genes = set()
    if not dirpath or not os.path.isdir(dirpath):
        return genes
    p = Path(dirpath)
    for f in sorted(p.glob("*.txt")):
        with open(f) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith("#"):
                    genes.add(strip_ensg_version(g))
    for f in sorted(p.glob("*.tsv")):
        try:
            df = pd.read_csv(f, sep="\t")
            if "gene_id" in df.columns:
                for g in df["gene_id"]:
                    genes.add(strip_ensg_version(str(g)))
        except Exception:
            pass
    return genes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build pan-cancer feature set from disease consensus catalogs."
    )
    ap.add_argument(
        "--catalog", action="append", default=[], metavar="DISEASE=PATH",
        help="Repeatable. Map disease name to its disease_consensus_catalog.tsv. "
             "Example: --catalog BRCA=/path/to/catalog.tsv --catalog NBL=/path/to/catalog.tsv"
    )
    ap.add_argument(
        "--isolate-dir", action="append", default=[], metavar="DISEASE=PATH",
        help="Repeatable. Map disease name to its isolate markers directory. "
             "Example: --isolate-dir BRCA=/path/to/isolate_degs/markers"
    )
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--min-cross-disease", type=int, default=2)
    ap.add_argument("--cap-specific", type=int, default=300,
                    help="Max disease-specific genes per disease")
    ap.add_argument("--cap-isolate", type=int, default=50,
                    help="Max isolate genes per disease")
    ap.add_argument("--remove-ribo-mt", action="store_true", default=False,
                    help="Filter ribosomal/mitochondrial genes (requires --gene-annotation-tsv)")
    ap.add_argument("--gene-annotation-tsv", default="",
                    help="TSV with gene_id, gene_name, gene_type for ribo/MT filtering")

    args = ap.parse_args()

    # Parse --catalog entries
    catalogs: dict[str, str] = {}
    for entry in args.catalog:
        if "=" not in entry:
            sys.exit(f"[ERROR] --catalog must be DISEASE=PATH, got: {entry}")
        disease, path = entry.split("=", 1)
        catalogs[disease.upper()] = path

    if not catalogs:
        sys.exit("[ERROR] At least one --catalog DISEASE=PATH is required")

    # Parse --isolate-dir entries
    isolate_dirs: dict[str, str] = {}
    for entry in args.isolate_dir:
        if "=" not in entry:
            sys.exit(f"[ERROR] --isolate-dir must be DISEASE=PATH, got: {entry}")
        disease, path = entry.split("=", 1)
        isolate_dirs[disease.upper()] = path

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load gene annotation if filtering requested
    annotation = {}
    filtering_applied = False
    if args.remove_ribo_mt:
        if args.gene_annotation_tsv and os.path.exists(args.gene_annotation_tsv):
            annotation = load_gene_annotation(args.gene_annotation_tsv)
            if annotation:
                filtering_applied = True
                print(f"[INFO] Loaded {len(annotation)} gene annotations for ribo/MT filtering",
                      file=sys.stderr)
            else:
                print("[WARN] Gene annotation file loaded but empty; skipping ribo/MT filtering",
                      file=sys.stderr)
        else:
            print("[WARN] --remove-ribo-mt set but no valid --gene-annotation-tsv provided; "
                  "skipping ribo/MT filtering", file=sys.stderr)

    # -----------------------------------------------------------------------
    # Load all catalogs
    # -----------------------------------------------------------------------
    disease_catalogs: dict[str, pd.DataFrame] = {}
    for disease, path in sorted(catalogs.items()):
        print(f"[INFO] Loading catalog for {disease}: {path}", file=sys.stderr)
        cat = load_catalog(path)
        if "disease" not in cat.columns:
            cat["disease"] = disease
        disease_catalogs[disease] = cat
        print(f"  {len(cat)} gene-rows, {cat['component'].nunique()} components",
              file=sys.stderr)

    # -----------------------------------------------------------------------
    # Step 1: PAN_CORE -- genes appearing in >= min_cross_disease diseases
    # -----------------------------------------------------------------------
    gene_diseases: dict[str, set] = defaultdict(set)
    gene_best_row: dict[str, dict] = {}
    gene_components: dict[str, set] = defaultdict(set)

    for disease, cat in disease_catalogs.items():
        for _, row in cat.iterrows():
            gid = row["gene_id"]
            gene_diseases[gid].add(disease)
            gene_components[gid].add(f"{disease}:{row['component']}")
            # Keep the row with the best rank_score
            rs = row.get("rank_score", 0) or 0
            prev = gene_best_row.get(gid)
            if prev is None or rs > (prev.get("rank_score", 0) or 0):
                gene_best_row[gid] = row.to_dict()

    pan_core_genes = {
        g for g, diseases in gene_diseases.items()
        if len(diseases) >= args.min_cross_disease
    }
    print(f"[INFO] PAN_CORE: {len(pan_core_genes)} genes", file=sys.stderr)

    # -----------------------------------------------------------------------
    # Step 2: DISEASE_SPECIFIC additions
    # -----------------------------------------------------------------------
    selected = set(pan_core_genes)
    disease_specific: dict[str, set] = {}

    for disease, cat in sorted(disease_catalogs.items()):
        candidates = cat[~cat["gene_id"].isin(selected)].copy()
        if "rank_score" in candidates.columns:
            candidates = candidates.sort_values("rank_score", ascending=False)
        chosen = []
        seen = set()
        for _, row in candidates.iterrows():
            gid = row["gene_id"]
            if gid in seen or gid in selected:
                continue
            seen.add(gid)
            chosen.append(gid)
            if len(chosen) >= args.cap_specific:
                break
        disease_specific[disease] = set(chosen)
        selected.update(chosen)
        print(f"[INFO] {disease}-specific: {len(chosen)} genes", file=sys.stderr)

    # -----------------------------------------------------------------------
    # Step 3: ISOLATE_SPECIFIC additions
    # -----------------------------------------------------------------------
    isolate_specific: dict[str, set] = {}
    for disease in sorted(catalogs.keys()):
        iso_dir = isolate_dirs.get(disease)
        iso_genes = load_isolate_dir(iso_dir) if iso_dir else set()
        candidates = iso_genes - selected
        chosen = sorted(candidates)[:args.cap_isolate]
        isolate_specific[disease] = set(chosen)
        selected.update(chosen)
        if iso_genes:
            print(f"[INFO] {disease}-isolates: {len(chosen)} genes "
                  f"(from {len(iso_genes)} candidates)", file=sys.stderr)

    # -----------------------------------------------------------------------
    # Apply ribo/MT filtering (if configured and annotation available)
    # -----------------------------------------------------------------------
    removed_ribo_mt = set()
    if filtering_applied:
        for g in list(selected):
            if is_ribo_or_mt(g, annotation):
                removed_ribo_mt.add(g)
        selected -= removed_ribo_mt
        pan_core_genes -= removed_ribo_mt
        for d in disease_specific:
            disease_specific[d] -= removed_ribo_mt
        for d in isolate_specific:
            isolate_specific[d] -= removed_ribo_mt
        print(f"[INFO] Removed {len(removed_ribo_mt)} ribosomal/mitochondrial genes",
              file=sys.stderr)

    # -----------------------------------------------------------------------
    # Build output table with provenance
    # -----------------------------------------------------------------------
    rows = []
    rank_counter = 0

    def _add_rows(gene_set, feature_class):
        nonlocal rank_counter
        for gid in sorted(gene_set):
            rank_counter += 1
            info = gene_best_row.get(gid, {})
            rows.append({
                "gene_id": gid,
                "direction": info.get("direction", "UNKNOWN"),
                "feature_class": feature_class,
                "source_diseases": ",".join(sorted(gene_diseases.get(gid, set()))),
                "source_components": ",".join(sorted(gene_components.get(gid, set()))),
                "disease_count": len(gene_diseases.get(gid, set())),
                "max_support_n": info.get("support_n", 0),
                "max_support_frac": info.get("support_frac", 0.0),
                "selection_rank": rank_counter,
            })

    _add_rows(pan_core_genes, "PAN_CORE")
    for disease in sorted(disease_specific.keys()):
        _add_rows(disease_specific[disease], "DISEASE_SPECIFIC")
    for disease in sorted(isolate_specific.keys()):
        _add_rows(isolate_specific[disease], "ISOLATE_SPECIFIC")

    features_df = pd.DataFrame(rows)

    # -----------------------------------------------------------------------
    # Write outputs
    # -----------------------------------------------------------------------
    features_df.to_csv(outdir / "pan_cancer_features.tsv", sep="\t", index=False)

    up_genes = set(features_df[features_df["direction"] == "UP"]["gene_id"])
    down_genes = set(features_df[features_df["direction"] == "DOWN"]["gene_id"])

    with open(outdir / "pan_cancer_features.UP.txt", "w") as f:
        for g in sorted(up_genes):
            f.write(f"{g}\n")

    with open(outdir / "pan_cancer_features.DOWN.txt", "w") as f:
        for g in sorted(down_genes):
            f.write(f"{g}\n")

    with open(outdir / "pan_cancer_features_clean.txt", "w") as f:
        for g in sorted(selected):
            f.write(f"{g}\n")

    # Summary
    summary_rows = []
    summary_rows.append({
        "category": "PAN_CORE", "n_genes": len(pan_core_genes),
        "n_up": len(pan_core_genes & up_genes),
        "n_down": len(pan_core_genes & down_genes),
    })
    for disease in sorted(disease_specific.keys()):
        gs = disease_specific[disease]
        summary_rows.append({
            "category": f"{disease}_SPECIFIC", "n_genes": len(gs),
            "n_up": len(gs & up_genes), "n_down": len(gs & down_genes),
        })
    for disease in sorted(isolate_specific.keys()):
        gs = isolate_specific[disease]
        summary_rows.append({
            "category": f"{disease}_ISOLATE", "n_genes": len(gs),
            "n_up": len(gs & up_genes), "n_down": len(gs & down_genes),
        })
    summary_rows.append({
        "category": "TOTAL", "n_genes": len(selected),
        "n_up": len(up_genes), "n_down": len(down_genes),
    })
    if removed_ribo_mt:
        summary_rows.append({
            "category": "REMOVED_RIBO_MT", "n_genes": len(removed_ribo_mt),
            "n_up": 0, "n_down": 0,
        })

    pd.DataFrame(summary_rows).to_csv(
        outdir / "pan_cancer_features_summary.tsv", sep="\t", index=False
    )

    print(f"\n[OK] Wrote {len(selected)} features to {outdir}/", file=sys.stderr)
    print(f"  pan_cancer_features.tsv       ({len(features_df)} rows)", file=sys.stderr)
    print(f"  pan_cancer_features.UP.txt    ({len(up_genes)} genes)", file=sys.stderr)
    print(f"  pan_cancer_features.DOWN.txt  ({len(down_genes)} genes)", file=sys.stderr)
    print(f"  pan_cancer_features_clean.txt ({len(selected)} genes)", file=sys.stderr)
    if filtering_applied:
        print(f"  Ribo/MT filtering: {len(removed_ribo_mt)} genes removed", file=sys.stderr)
    elif args.remove_ribo_mt:
        print("  Ribo/MT filtering: requested but no annotation available", file=sys.stderr)


if __name__ == "__main__":
    main()
