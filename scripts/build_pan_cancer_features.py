#!/usr/bin/env python3
"""
build_pan_cancer_features.py

Config-driven merge step. Consumes the per-profile consensus summary
exports produced by build_component_consensus_markers_v2.py (one
directory per profile containing gene_recurrence_across_components.tsv)
and combines them into a marker-derived pan-cancer feature panel.

Selection hierarchy:
  1. PAN_CORE          -- genes present in >= min_cross_disease profiles
  2. PROFILE_SPECIFIC  -- per-profile additions (capped per direction)
  3. ISOLATE_RESCUE    -- optional per-profile isolate-vs-rest marker rescue
                          genes, capped and ranked deterministically

Outputs:
  pan_cancer_features.tsv
  pan_cancer_features.UP.txt
  pan_cancer_features.DOWN.txt
  pan_cancer_features_clean.txt
  pan_cancer_feature_build_summary.tsv
  pan_cancer_feature_build_report.md
  pan_cancer_features_done.txt
"""

import argparse
import os
import re
import sys
from pathlib import Path

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
    """Load a TSV with gene_id and gene_name columns for ribo/MT filtering."""
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        print(f"[WARN] gene_annotation_tsv missing 'gene_id' column: {path}", file=sys.stderr)
        return {}
    name_col = None
    for candidate in ["gene_name", "symbol", "gene_symbol"]:
        if candidate in df.columns:
            name_col = candidate
            break
    if name_col is None:
        print(f"[WARN] gene_annotation_tsv missing gene symbol column: {path}", file=sys.stderr)
        return {}
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)
    return dict(zip(df["gene_id"], df[name_col].astype(str)))


def is_ribo_or_mt(gene_id: str, annotation: dict) -> bool:
    if not annotation:
        return False
    gene_symbol = annotation.get(gene_id, "").upper()
    if not gene_symbol:
        return False
    return bool(re.match(r"^(RPL|RPS|MRPL|MRPS|MT-|HIST)", gene_symbol))


# ---------------------------------------------------------------------------
# Catalog loading
# ---------------------------------------------------------------------------

def load_recurrence_table(profile: str, summary_dir: Path) -> pd.DataFrame:
    path = Path(summary_dir) / "gene_recurrence_across_components.tsv"
    if not path.exists():
        sys.exit(f"[ERROR] Missing recurrence table for profile {profile}: {path}")
    df = pd.read_csv(path, sep="\t")
    required = [
        "gene_id", "direction", "freq", "n_components",
        "mean_support_n", "mean_support_frac", "mean_abs_log2FC", "max_rank_score"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {path} missing required columns: {missing}")
    df = normalize_gene_ids(df)
    df["direction"] = df["direction"].str.upper()
    df["profile"] = profile
    return df


def load_gene_list(path: Path) -> list[str]:
    """Load a one-gene-per-line marker file, preserving file order."""
    genes = []
    with open(path) as handle:
        for line in handle:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.append(strip_ensg_version(gene))
    return genes


def load_contrast_table(path: Path) -> pd.DataFrame:
    """Load a DESeq2 contrast table and normalize gene IDs when available."""
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        return pd.DataFrame()
    return normalize_gene_ids(df)


def resolve_marker_path(marker_dir: Path, manifest_path: Path | None, row: pd.Series) -> Path | None:
    """Resolve a marker file from a manifest row without hardcoding topN."""
    if "marker_file" in row.index and pd.notna(row["marker_file"]):
        marker_file = Path(str(row["marker_file"]))
        if marker_file.is_absolute():
            return marker_file
        base = manifest_path.parent.parent if manifest_path is not None else marker_dir.parent
        return base / marker_file

    contrast = str(row["contrast"])
    matches = sorted(marker_dir.glob(f"{contrast}_markers_top*.txt"))
    return matches[0] if matches else None


def load_isolate_rescue_candidates(
    profile: str,
    marker_dir: Path,
    manifest_path: Path | None = None
) -> pd.DataFrame:
    """
    Load isolate-vs-rest marker candidates for one profile.

    Candidate ranking is intentionally deterministic.  Genes are ordered by:
      1. number of isolate contrasts supporting the gene (descending)
      2. best adjusted p value across isolate contrasts (ascending)
      3. max absolute log2 fold change (descending)
      4. best within-file marker rank (ascending)
      5. gene_id (ascending)

    This avoids the non-reproducible set-order tie behavior in the legacy
    isolate rescue script while preserving the same high-level criterion.
    """
    marker_dir = Path(marker_dir)
    if not marker_dir.exists():
        print(f"[WARN] Isolate marker directory missing for {profile}: {marker_dir}",
              file=sys.stderr)
        return pd.DataFrame()

    records = []
    table_dir = marker_dir.parent / "tables"

    if manifest_path is not None and manifest_path.exists():
        manifest = pd.read_csv(manifest_path, sep="\t")
        if "contrast" not in manifest.columns:
            sys.exit(f"[ERROR] Marker manifest missing contrast column: {manifest_path}")
        manifest = manifest[manifest["contrast"].astype(str).str.startswith("isolate_")].copy()
        manifest_rows = [row for _, row in manifest.sort_values("contrast").iterrows()]
    else:
        manifest_rows = []
        for marker_path in sorted(marker_dir.glob("isolate_*_markers_top*.txt")):
            suffix_match = re.search(r"_markers_top\d+\.txt$", marker_path.name)
            contrast = marker_path.name[:suffix_match.start()] if suffix_match else marker_path.stem
            manifest_rows.append(pd.Series({"contrast": contrast, "marker_file": str(marker_path)}))

    for row in manifest_rows:
        marker_path = resolve_marker_path(marker_dir, manifest_path, row)
        if marker_path is None or not marker_path.exists():
            print(f"[WARN] Missing isolate marker file for {profile}: {row['contrast']}",
                  file=sys.stderr)
            continue
        suffix_match = re.search(r"_markers_top\d+\.txt$", marker_path.name)
        contrast = marker_path.name[:suffix_match.start()] if suffix_match else marker_path.stem
        table = load_contrast_table(table_dir / f"{contrast}.tsv")
        table_by_gene = {}
        if not table.empty:
            table_by_gene = {str(row["gene_id"]): row for _, row in table.iterrows()}

        for marker_rank, gene_id in enumerate(load_gene_list(marker_path), start=1):
            row = table_by_gene.get(gene_id)
            padj = pd.NA
            abs_lfc = pd.NA
            direction = "UNKNOWN"
            if row is not None:
                if "padj" in row.index and pd.notna(row["padj"]):
                    padj = float(row["padj"])
                if "log2FoldChange" in row.index and pd.notna(row["log2FoldChange"]):
                    lfc = float(row["log2FoldChange"])
                    abs_lfc = abs(lfc)
                    direction = "UP" if lfc >= 0 else "DOWN"
            records.append({
                "profile": profile,
                "gene_id": gene_id,
                "source_contrast": contrast,
                "source_file": str(marker_path),
                "marker_rank": marker_rank,
                "padj": padj,
                "abs_log2FC": abs_lfc,
                "direction": direction,
            })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df["padj_sort"] = pd.to_numeric(df["padj"], errors="coerce").fillna(float("inf"))
    df["abs_lfc_sort"] = pd.to_numeric(df["abs_log2FC"], errors="coerce").fillna(0.0)

    def collapse_direction(values: pd.Series) -> str:
        known = [v for v in values if v in {"UP", "DOWN"}]
        if not known:
            return "UNKNOWN"
        counts = pd.Series(known).value_counts()
        if len(counts) > 1 and counts.iloc[0] == counts.iloc[1]:
            best = df.loc[values.index].sort_values(
                ["abs_lfc_sort", "padj_sort", "gene_id"],
                ascending=[False, True, True]
            ).iloc[0]
            return best["direction"]
        return counts.index[0]

    grouped = df.groupby("gene_id", sort=True).agg(
        n_isolate_contrasts=("source_contrast", "nunique"),
        best_padj=("padj_sort", "min"),
        max_abs_log2FC=("abs_lfc_sort", "max"),
        best_marker_rank=("marker_rank", "min"),
        source_contrast=("source_contrast", lambda x: ";".join(sorted(set(x)))),
        source_file=("source_file", lambda x: ";".join(sorted(set(x)))),
        direction=("direction", collapse_direction),
    ).reset_index()

    grouped["profile"] = profile
    return grouped


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build pan-cancer feature set from disease consensus catalogs."
    )
    ap.add_argument(
        "--profile-dir", action="append", default=[], metavar="PROFILE=PATH",
        help="Repeatable. Map profile name to its consensus summary directory."
    )
    ap.add_argument(
        "--profile-marker-dir", action="append", default=[], metavar="PROFILE=PATH",
        help="Repeatable. Map profile name to its DESeq2 markers directory for isolate rescue."
    )
    ap.add_argument(
        "--profile-marker-manifest", action="append", default=[], metavar="PROFILE=PATH",
        help="Repeatable. Map profile name to marker_sets_manifest.tsv for isolate rescue."
    )
    ap.add_argument("--output-dir", required=True, help="Output directory")
    ap.add_argument("--min-cross-disease", type=int, default=2,
                    help="Minimum number of profiles required for PAN_CORE status")
    ap.add_argument("--cap-specific", type=int, default=300,
                    help="Max PROFILE_SPECIFIC genes per profile per direction")
    ap.add_argument("--cap-isolate", type=int, default=50,
                    help="Max ISOLATE_RESCUE genes per profile. 0 disables isolate rescue.")
    ap.add_argument("--remove-ribo-mt", action="store_true", default=False,
                    help="Drop RPL/RPS/MRPL/MRPS/MT-/HIST genes when annotation provided")
    ap.add_argument("--gene-annotation-tsv", default="",
                    help="TSV with gene_id + gene_symbol columns for ribo/MT filtering")

    args = ap.parse_args()

    profile_dirs = {}
    for entry in args.profile_dir:
        if "=" not in entry:
            sys.exit(f"[ERROR] --profile-dir must be PROFILE=PATH, got: {entry}")
        profile, path = entry.split("=", 1)
        profile = profile.strip()
        if not profile:
            sys.exit(f"[ERROR] Invalid profile name in --profile-dir: {entry}")
        profile_dirs[profile] = path.strip()

    if not profile_dirs:
        sys.exit("[ERROR] At least one --profile-dir PROFILE=PATH is required")

    profile_marker_dirs = {}
    for entry in args.profile_marker_dir:
        if "=" not in entry:
            sys.exit(f"[ERROR] --profile-marker-dir must be PROFILE=PATH, got: {entry}")
        profile, path = entry.split("=", 1)
        profile = profile.strip()
        if not profile:
            sys.exit(f"[ERROR] Invalid profile name in --profile-marker-dir: {entry}")
        profile_marker_dirs[profile] = path.strip()

    profile_marker_manifests = {}
    for entry in args.profile_marker_manifest:
        if "=" not in entry:
            sys.exit(f"[ERROR] --profile-marker-manifest must be PROFILE=PATH, got: {entry}")
        profile, path = entry.split("=", 1)
        profile = profile.strip()
        if not profile:
            sys.exit(f"[ERROR] Invalid profile name in --profile-marker-manifest: {entry}")
        profile_marker_manifests[profile] = path.strip()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    recurrence_frames = []
    for profile, summary_dir in sorted(profile_dirs.items()):
        print(f"[INFO] Loading recurrence summary for {profile}: {summary_dir}", file=sys.stderr)
        recurrence_frames.append(load_recurrence_table(profile, summary_dir))

    if not recurrence_frames:
        sys.exit("[ERROR] No recurrence tables loaded")

    data = pd.concat(recurrence_frames, ignore_index=True)
    if data.empty:
        sys.exit("[ERROR] Recurrence tables contain no genes")

    min_cross = max(int(args.min_cross_disease), 1)
    cap_specific = max(int(args.cap_specific), 0)
    cap_isolate = max(int(args.cap_isolate), 0)

    agg = data.groupby(["direction", "gene_id"]).agg(
        n_profiles=("profile", "nunique"),
        mean_abs_log2FC=("mean_abs_log2FC", "mean"),
        mean_support_frac=("mean_support_frac", "mean"),
        profiles_present=("profile", lambda x: sorted(set(x)))
    ).reset_index()
    agg["profiles_present"] = agg["profiles_present"].apply(lambda lst: ",".join(lst))
    pan_core_df = agg[agg["n_profiles"] >= min_cross].copy()

    pan_core_records = []
    pan_core_genes_by_dir = {"UP": set(), "DOWN": set()}
    for direction in sorted(pan_core_df["direction"].unique()):
        dir_df = pan_core_df[pan_core_df["direction"] == direction].copy()
        if dir_df.empty:
            continue
        dir_df = dir_df.sort_values(
            ["n_profiles", "mean_support_frac", "mean_abs_log2FC", "gene_id"],
            ascending=[False, False, False, True]
        )
        pan_core_genes_by_dir[direction] = set(dir_df["gene_id"].tolist())
        for _, row in dir_df.iterrows():
            pan_core_records.append({
                "gene_id": row["gene_id"],
                "direction": direction,
                "feature_class": "PAN_CORE",
                "owner_profile": "PAN",
                "profiles_present": row["profiles_present"],
                "n_profiles": int(row["n_profiles"]),
                "mean_abs_log2FC": float(row["mean_abs_log2FC"]),
                "mean_support_frac": float(row["mean_support_frac"]),
                "source_contrast": "gene_recurrence_across_components",
                "source_file": "",
                "source_rank": pd.NA,
            })

    specific_records = []
    for direction in ["UP", "DOWN"]:
        dir_core_genes = pan_core_genes_by_dir.get(direction, set())
        for profile in sorted(profile_dirs.keys()):
            subset = data[(data["profile"] == profile) & (data["direction"] == direction)].copy()
            if subset.empty:
                continue
            subset = subset[~subset["gene_id"].isin(dir_core_genes)]
            if subset.empty:
                continue
            subset = subset.sort_values(
                ["freq", "mean_support_frac", "mean_abs_log2FC", "gene_id"],
                ascending=[False, False, False, True]
            )
            if cap_specific > 0:
                subset = subset.head(cap_specific)
            for _, row in subset.iterrows():
                specific_records.append({
                    "gene_id": row["gene_id"],
                    "direction": direction,
                    "feature_class": "PROFILE_SPECIFIC",
                    "owner_profile": profile,
                    "profiles_present": profile,
                    "n_profiles": 1,
                    "mean_abs_log2FC": float(row["mean_abs_log2FC"]),
                    "mean_support_frac": float(row["mean_support_frac"]),
                    "source_contrast": "gene_recurrence_across_components",
                    "source_file": "",
                    "source_rank": pd.NA,
                })

    feature_rows = pan_core_records + specific_records
    if not feature_rows:
        sys.exit("[ERROR] No pan-cancer features selected")

    features_df = pd.DataFrame(feature_rows)
    selected_genes = set(features_df["gene_id"])

    isolate_records = []
    isolate_audit_records = []
    if cap_isolate > 0:
        print(f"[INFO] isolate rescue cap = {cap_isolate} per profile",
              file=sys.stderr)
        for profile in sorted(profile_dirs.keys()):
            marker_dir = profile_marker_dirs.get(profile)
            if not marker_dir:
                print(f"[WARN] No --profile-marker-dir provided for {profile}; "
                      "skipping isolate rescue for this profile",
                      file=sys.stderr)
                continue
            manifest_path = profile_marker_manifests.get(profile)
            candidates = load_isolate_rescue_candidates(
                profile,
                Path(marker_dir),
                Path(manifest_path) if manifest_path else None
            )
            if candidates.empty:
                print(f"[INFO] {profile}: no isolate rescue candidates", file=sys.stderr)
                continue
            candidates = candidates[~candidates["gene_id"].isin(selected_genes)].copy()
            candidates = candidates.sort_values(
                [
                    "n_isolate_contrasts", "best_padj", "max_abs_log2FC",
                    "best_marker_rank", "gene_id"
                ],
                ascending=[False, True, False, True, True],
            ).head(cap_isolate)

            for source_rank, (_, row) in enumerate(candidates.iterrows(), start=1):
                isolate_audit_records.append({
                    "profile": profile,
                    "gene_id": row["gene_id"],
                    "direction": row["direction"],
                    "n_isolate_contrasts": int(row["n_isolate_contrasts"]),
                    "best_padj": float(row["best_padj"]),
                    "max_abs_log2FC": float(row["max_abs_log2FC"]),
                    "best_marker_rank": int(row["best_marker_rank"]),
                    "source_rank": source_rank,
                    "source_contrast": row["source_contrast"],
                    "source_file": row["source_file"],
                })
                isolate_records.append({
                    "gene_id": row["gene_id"],
                    "direction": row["direction"],
                    "feature_class": "ISOLATE_RESCUE",
                    "owner_profile": profile,
                    "profiles_present": profile,
                    "n_profiles": 1,
                    "mean_abs_log2FC": float(row["max_abs_log2FC"]),
                    "mean_support_frac": float(row["n_isolate_contrasts"]),
                    "source_contrast": row["source_contrast"],
                    "source_file": row["source_file"],
                    "source_rank": source_rank,
                })
            selected_genes.update(candidates["gene_id"])
            print(f"[INFO] {profile}: selected {len(candidates)} isolate rescue genes",
                  file=sys.stderr)

    if isolate_records:
        features_df = pd.concat([features_df, pd.DataFrame(isolate_records)],
                                ignore_index=True)

    isolate_audit_path = outdir / "isolate_rescue_audit.tsv"
    isolate_audit_columns = [
        "profile", "gene_id", "direction", "n_isolate_contrasts",
        "best_padj", "max_abs_log2FC", "best_marker_rank", "source_rank",
        "source_contrast", "source_file",
    ]
    pd.DataFrame(isolate_audit_records, columns=isolate_audit_columns).to_csv(
        isolate_audit_path, sep="\t", index=False
    )

    annotation = {}
    filtering_applied = False
    removed_ribo_mt = set()
    if args.remove_ribo_mt:
        if args.gene_annotation_tsv and os.path.exists(args.gene_annotation_tsv):
            annotation = load_gene_annotation(args.gene_annotation_tsv)
            if annotation:
                filtering_applied = True
                mask = features_df["gene_id"].map(lambda g: is_ribo_or_mt(g, annotation))
                removed_ribo_mt = set(features_df.loc[mask, "gene_id"])
                features_df = features_df[~mask].reset_index(drop=True)
                print(f"[INFO] Removed {len(removed_ribo_mt)} ribosomal/mitochondrial genes",
                      file=sys.stderr)
            else:
                print("[WARN] Annotation file provided but gene symbols missing; skipping ribo/MT filtering",
                      file=sys.stderr)
        else:
            print("[WARN] --remove-ribo-mt requested but --gene-annotation-tsv not provided; skipping",
                  file=sys.stderr)

    if features_df.empty:
        sys.exit("[ERROR] No features remain after filtering")

    features_df["selection_rank"] = range(1, len(features_df) + 1)
    features_df = features_df[[
        "gene_id", "direction", "feature_class", "owner_profile",
        "profiles_present", "n_profiles", "mean_abs_log2FC",
        "mean_support_frac", "source_contrast", "source_file",
        "source_rank", "selection_rank"
    ]]

    features_path = outdir / "pan_cancer_features.tsv"
    features_df.to_csv(features_path, sep="\t", index=False)

    def _write_gene_list(path: Path, genes):
        with open(path, "w") as handle:
            for gene in genes:
                handle.write(f"{gene}\n")

    up_genes_ordered = features_df[features_df["direction"] == "UP"]["gene_id"].tolist()
    down_genes_ordered = features_df[features_df["direction"] == "DOWN"]["gene_id"].tolist()
    clean_genes = list(dict.fromkeys(features_df["gene_id"].tolist()))

    _write_gene_list(outdir / "pan_cancer_features.UP.txt", up_genes_ordered)
    _write_gene_list(outdir / "pan_cancer_features.DOWN.txt", down_genes_ordered)
    _write_gene_list(outdir / "pan_cancer_features_clean.txt", clean_genes)

    summary_rows = []
    for direction in ["UP", "DOWN"]:
        dir_df = features_df[features_df["direction"] == direction]
        pan_n = len(dir_df[dir_df["feature_class"] == "PAN_CORE"])
        if pan_n:
            summary_rows.append({
                "direction": direction,
                "feature_class": "PAN_CORE",
                "owner_profile": "PAN",
                "n_genes": pan_n,
            })
        for profile in sorted(profile_dirs.keys()):
            prof_n = len(dir_df[(dir_df["feature_class"] == "PROFILE_SPECIFIC") &
                                (dir_df["owner_profile"] == profile)])
            if prof_n:
                summary_rows.append({
                    "direction": direction,
                    "feature_class": "PROFILE_SPECIFIC",
                    "owner_profile": profile,
                    "n_genes": prof_n,
                })
            iso_n = len(dir_df[(dir_df["feature_class"] == "ISOLATE_RESCUE") &
                               (dir_df["owner_profile"] == profile)])
            if iso_n:
                summary_rows.append({
                    "direction": direction,
                    "feature_class": "ISOLATE_RESCUE",
                    "owner_profile": profile,
                    "n_genes": iso_n,
                })
        summary_rows.append({
            "direction": direction,
            "feature_class": "TOTAL",
            "owner_profile": "ALL",
            "n_genes": len(dir_df),
        })

    summary_rows.append({
        "direction": "ALL",
        "feature_class": "TOTAL",
        "owner_profile": "ALL",
        "n_genes": len(features_df),
    })
    if filtering_applied:
        summary_rows.append({
            "direction": "ALL",
            "feature_class": "REMOVED_RIBO_MT",
            "owner_profile": "PAN",
            "n_genes": len(removed_ribo_mt),
        })

    summary_path = outdir / "pan_cancer_feature_build_summary.tsv"
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_path, sep="\t", index=False)

    report_path = outdir / "pan_cancer_feature_build_report.md"
    by_class = (
        features_df.groupby(["feature_class", "owner_profile"], sort=True)
        .size()
        .reset_index(name="n_genes")
        .sort_values(["feature_class", "owner_profile"])
    )
    with open(report_path, "w") as handle:
        handle.write("# Pan-cancer feature build report\n\n")
        handle.write("## Inputs\n\n")
        handle.write("| Profile | Consensus summary directory | Isolate marker manifest |\n")
        handle.write("|---|---|---|\n")
        for profile in sorted(profile_dirs.keys()):
            manifest = profile_marker_manifests.get(profile, "")
            handle.write(f"| {profile} | `{profile_dirs[profile]}` | `{manifest}` |\n")

        handle.write("\n## Selection rule\n\n")
        handle.write(
            f"- PAN_CORE: genes recurrent in at least {min_cross} disease profiles "
            "using the consensus recurrence tables.\n"
        )
        handle.write(
            f"- PROFILE_SPECIFIC: disease-specific recurrent genes not already selected "
            f"as PAN_CORE, capped at {cap_specific} genes per profile per direction.\n"
        )
        if cap_isolate > 0:
            handle.write(
                f"- ISOLATE_RESCUE: graph-derived isolate-vs-rest marker genes from "
                f"`marker_sets_manifest.tsv`, excluding already selected genes and capped "
                f"at {cap_isolate} genes per profile.\n"
            )
            handle.write(
                "- Isolate rescue ranking: number of supporting isolate contrasts desc, "
                "best adjusted P value asc, max absolute shrunken log2FC desc, marker "
                "rank asc, Ensembl gene ID asc.\n"
            )
        else:
            handle.write("- ISOLATE_RESCUE: disabled by `cap_isolate = 0`.\n")
        handle.write(
            "- The final clean list preserves first selection order and removes duplicate "
            "gene IDs after Ensembl version stripping.\n"
        )

        handle.write("\n## Build summary\n\n")
        handle.write(f"- Feature rows: {len(features_df)}\n")
        handle.write(f"- Unique clean genes: {len(clean_genes)}\n")
        handle.write(f"- UP genes: {len(up_genes_ordered)}\n")
        handle.write(f"- DOWN genes: {len(down_genes_ordered)}\n")
        if filtering_applied:
            handle.write(f"- Ribosomal/mitochondrial/histone genes removed: {len(removed_ribo_mt)}\n")
        elif args.remove_ribo_mt:
            handle.write("- Ribosomal/mitochondrial/histone filtering requested but not applied because annotation was unavailable.\n")
        else:
            handle.write("- Ribosomal/mitochondrial/histone filtering: not requested.\n")

        handle.write("\n## Counts by feature class\n\n")
        handle.write("| Feature class | Owner profile | Genes |\n")
        handle.write("|---|---:|---:|\n")
        for _, row in by_class.iterrows():
            handle.write(f"| {row['feature_class']} | {row['owner_profile']} | {int(row['n_genes'])} |\n")

        handle.write("\n## Output files\n\n")
        handle.write(f"- `{features_path}`: full feature manifest with provenance columns.\n")
        handle.write(f"- `{outdir / 'pan_cancer_features_clean.txt'}`: ordered unique gene list for expression extraction.\n")
        handle.write(f"- `{summary_path}`: compact machine-readable count summary.\n")
        handle.write(f"- `{isolate_audit_path}`: selected isolate-rescue candidates and deterministic ranking fields.\n")

    done_path = outdir / "pan_cancer_features_done.txt"
    with open(done_path, "w") as handle:
        handle.write("done\n")

    print(f"\n[OK] Wrote {len(features_df)} features to {outdir}/", file=sys.stderr)
    print(f"  pan_cancer_features.tsv       ({len(features_df)} rows)", file=sys.stderr)
    print(f"  pan_cancer_features.UP.txt    ({len(up_genes_ordered)} genes)", file=sys.stderr)
    print(f"  pan_cancer_features.DOWN.txt  ({len(down_genes_ordered)} genes)", file=sys.stderr)
    print(f"  pan_cancer_features_clean.txt ({len(clean_genes)} genes)", file=sys.stderr)
    print(f"  isolate_rescue_audit.tsv      ({len(isolate_audit_records)} rows)", file=sys.stderr)
    print(f"  pan_cancer_feature_build_report.md", file=sys.stderr)
    if filtering_applied:
        print(f"  Ribo/MT filtering: {len(removed_ribo_mt)} genes removed", file=sys.stderr)
    elif args.remove_ribo_mt:
        print("  Ribo/MT filtering: requested but annotation unavailable", file=sys.stderr)


if __name__ == "__main__":
    main()
