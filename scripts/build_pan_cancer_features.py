#!/usr/bin/env python3
"""
Graph-derived marker aggregation and pan-cancer feature selection.

Inputs are limited to the contrast-level marker manifest and retained
per-contrast marker tables produced by the DESeq2 marker-prioritisation module.
This script does not read raw RNA-seq counts, fit DESeq2 models, estimate size
factors, construct DESeq2 contrast factors, or repeat contrast-level filtering.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


METHOD_NAME = "graph_derived_pan_cancer_feature_selection_v1_revised"
DEFAULT_AUDIT_OUTPUT_PREFIX = "ranked_marker_source_panel"
RECURRENCE_MIN_COUNT = 2
# The three empirical acceptance quantiles are configuration-owned
# (marker_postprocessing.pan_cancer.empirical_quantile_thresholds) and arrive on
# the command line. They are deliberately not restated as module constants: the
# former constants both duplicated the configured values and rejected any run
# whose configuration differed from them, which made configuration inert.
MARKER_EVIDENCE_STRATUM_ORDER = ("anchor_associated", "isolate_associated")
CANDIDATE_POOL_ORDER = {
    "recurrent": 1,
    "singleton_candidate": 2,
    "non_recurrent_candidate": 3,
}
SELECTION_BASIS_PRIORITY = {
    "recurrent": 0,
    "accepted_singleton": 1,
    "accepted_non_recurrent": 2,
}
OBSOLETE_RULE_COLUMN_PARTS = (
    ("selected", "rule", "key"),
    ("selected", "sequential", "rule"),
    ("selection", "rule", "identifier"),
    ("selection", "rule", "description"),
    ("selected", "rule", "pass"),
    ("passes", "feature", "selection", "rule"),
    ("minimum", "required", "unique", "gene", "count"),
    ("fallback", "rule", "applied"),
)


def obsolete_rule_columns() -> set[str]:
    return {"_".join(parts) for parts in OBSOLETE_RULE_COLUMN_PARTS}


def strip_ensg_version(gene_id: object) -> str:
    if pd.isna(gene_id):
        return ""
    return str(gene_id).strip().split(".", 1)[0]


def marker_source_class_from_evidence_stratum(marker_evidence_stratum: str) -> str:
    if marker_evidence_stratum == "anchor_associated":
        return "anchor"
    if marker_evidence_stratum == "isolate_associated":
        return "isolate"
    return str(marker_evidence_stratum or "").replace("_associated", "")


def marker_evidence_stratum_from_contrast_id(contrast_id: str) -> str:
    if str(contrast_id).startswith("anchor_"):
        return "anchor_associated"
    if str(contrast_id).startswith("isolate_"):
        return "isolate_associated"
    return "unclassified"


def stratum_sort_value(marker_evidence_stratum: str) -> int:
    try:
        return MARKER_EVIDENCE_STRATUM_ORDER.index(marker_evidence_stratum)
    except ValueError:
        return len(MARKER_EVIDENCE_STRATUM_ORDER)


def selection_basis_sort_value(selection_basis: str) -> int:
    return SELECTION_BASIS_PRIORITY.get(
        str(selection_basis), len(SELECTION_BASIS_PRIORITY)
    )


def parse_profile_paths(entries: list[str], label: str) -> tuple[list[str], dict[str, Path]]:
    order: list[str] = []
    paths: dict[str, Path] = {}
    for entry in entries:
        if "=" not in entry:
            sys.exit(f"[ERROR] {label} must be PROFILE=PATH, got: {entry}")
        profile, path = entry.split("=", 1)
        profile = profile.strip()
        if not profile:
            sys.exit(f"[ERROR] Empty profile in {label}: {entry}")
        if profile not in paths:
            order.append(profile)
        paths[profile] = Path(path.strip())
    return order, paths


def resolve_manifest_path(manifest_path: Path | None) -> Path:
    if manifest_path is not None:
        return manifest_path
    sys.exit("[ERROR] Active feature construction requires explicit --profile-marker-manifest")


def resolve_relative_path(deseq2_dir: Path, value: object) -> Path:
    path = Path(str(value))
    if path.is_absolute():
        return path
    return deseq2_dir / path


def parse_float(value: object, default: float = math.nan) -> float:
    out = pd.to_numeric(value, errors="coerce")
    if pd.isna(out):
        return default
    return float(out)


def finite_numeric(values: pd.Series) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce")
    vals = vals[vals.map(lambda x: math.isfinite(float(x)) if not pd.isna(x) else False)]
    return vals.astype(float)


def safe_min(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    return default if vals.empty else float(vals.min())


def safe_max(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    return default if vals.empty else float(vals.max())


def safe_mean(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    return default if vals.empty else float(vals.mean())


def safe_median(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    return default if vals.empty else float(vals.median())


def safe_quantile(values: pd.Series, quantile_probability: float, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    return default if vals.empty else float(vals.quantile(quantile_probability))


def join_unique(values, sep: str = ";") -> str:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return sep.join(out)


def write_gene_list(path: Path, genes: list[str]) -> None:
    text = "\n".join(genes)
    if genes:
        text += "\n"
    path.write_text(text)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_manifest_row(path: Path, category: str, description: str) -> dict[str, object]:
    exists = path.exists()
    stat = path.stat() if exists else None
    return {
        "filename": path.name,
        "path": str(path),
        "file_size_bytes": int(stat.st_size) if stat else pd.NA,
        "modification_time_utc": (
            datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()
            if stat else ""
        ),
        "sha256": sha256_file(path) if exists and path.is_file() else "",
        "description": description,
        "category": category,
    }


def load_gene_annotation(path: str) -> dict[str, dict[str, str]]:
    if not path or not os.path.exists(path):
        return {}
    annotation_table = pd.read_csv(path, sep="\t")
    gene_id_column = next(
        (
            col for col in ("gene_id", "ensembl_id", "ensembl_gene_id")
            if col in annotation_table.columns
        ),
        None,
    )
    symbol_column = next(
        (
            col for col in ("symbol", "gene_symbol", "gene_name", "hgnc_symbol")
            if col in annotation_table.columns
        ),
        None,
    )
    if gene_id_column is None or symbol_column is None:
        print(f"[WARN] Gene annotation missing required ID/symbol columns: {path}", file=sys.stderr)
        return {}
    name_column = "gene_name" if "gene_name" in annotation_table.columns else symbol_column
    annotations: dict[str, dict[str, str]] = {}
    for row in annotation_table.itertuples(index=False):
        clean_gene_id = strip_ensg_version(getattr(row, gene_id_column))
        if clean_gene_id:
            annotations[clean_gene_id] = {
                "gene_symbol": str(getattr(row, symbol_column) or "").strip(),
                "gene_name": str(getattr(row, name_column) or "").strip(),
            }
    return annotations


def is_ribo_mt_or_histone(clean_gene_id: str, annotation: dict[str, dict[str, str]]) -> bool:
    symbol = annotation.get(clean_gene_id, {}).get("gene_symbol", "").upper()
    prefixes = ("RPL", "RPS", "MRPL", "MRPS", "MT-", "HIST")
    return bool(symbol and symbol.startswith(prefixes))


def normalise_manifest(manifest: pd.DataFrame, profile: str) -> pd.DataFrame:
    manifest = manifest.copy()
    required = {
        "cancer_type",
        "contrast_id",
        "contrast_type",
        "marker_evidence_stratum",
        "focal_profile_id",
        "focal_component_id",
        "reference_definition",
        "marker_table_path",
        "marker_gene_list_path",
        "result_table_path",
        "n_result_genes",
        "n_pvalue_nonmissing",
        "n_padj_nonmissing",
        "n_significant_before_effect_filter",
        "n_markers_before_cap",
        "n_markers_after_cap",
        "adjusted_p_value_threshold",
        "minimum_base_mean",
        "minimum_absolute_shrunken_log2fc",
        "maximum_markers_per_contrast",
    }
    missing = sorted(required - set(manifest.columns))
    if missing:
        sys.exit(f"[ERROR] {profile} canonical contrast manifest missing columns: {', '.join(missing)}")
    text_columns = [
        "cancer_type",
        "contrast_id",
        "contrast_type",
        "marker_evidence_stratum",
        "focal_profile_id",
        "focal_component_id",
        "reference_definition",
        "marker_table_path",
        "marker_gene_list_path",
        "result_table_path",
    ]
    for column in text_columns:
        manifest[column] = manifest[column].map(
            lambda value: "" if pd.isna(value) else str(value).strip()
        )
    required_nonempty = [
        "cancer_type",
        "contrast_id",
        "contrast_type",
        "marker_evidence_stratum",
        "marker_table_path",
        "marker_gene_list_path",
        "result_table_path",
    ]
    for column in required_nonempty:
        empty_rows = manifest[~manifest[column].astype(bool)]
        if not empty_rows.empty:
            sys.exit(
                f"[ERROR] {profile} canonical contrast manifest has empty {column} "
                f"for row(s): {', '.join(map(str, empty_rows.index.tolist()))}"
            )
    if manifest["contrast_id"].duplicated().any():
        duplicates = sorted(manifest.loc[manifest["contrast_id"].duplicated(), "contrast_id"].astype(str).unique())
        sys.exit(f"[ERROR] {profile} duplicate contrast_id values in manifest: {', '.join(duplicates)}")
    unknown_strata = sorted(set(manifest["marker_evidence_stratum"].astype(str)) - set(MARKER_EVIDENCE_STRATUM_ORDER))
    if unknown_strata:
        sys.exit(f"[ERROR] {profile} unknown marker_evidence_stratum values: {', '.join(unknown_strata)}")
    allowed_contrast_types = {
        "anchor_focal_vs_outside_focal_component",
        "isolate_focal_vs_other_same_cancer",
    }
    unknown_contrast_types = sorted(set(manifest["contrast_type"].astype(str)) - allowed_contrast_types)
    if unknown_contrast_types:
        sys.exit(f"[ERROR] {profile} unknown contrast_type values: {', '.join(unknown_contrast_types)}")
    anchor_missing_component = manifest[
        (manifest["marker_evidence_stratum"] == "anchor_associated")
        & ~manifest["focal_component_id"].astype(str).str.strip().astype(bool)
    ]
    if not anchor_missing_component.empty:
        sys.exit(
            f"[ERROR] {profile} anchor contrasts missing focal_component_id: "
            + ", ".join(anchor_missing_component["contrast_id"].astype(str))
        )
    manifest["n_markers_after_cap"] = pd.to_numeric(manifest["n_markers_after_cap"], errors="raise").astype(int)
    return manifest


def load_retained_marker_table(marker_table_path: Path, profile: str, contrast_id: str) -> pd.DataFrame:
    if not marker_table_path.exists():
        sys.exit(f"[ERROR] Missing retained marker table: {marker_table_path}")
    marker_table = pd.read_csv(marker_table_path, sep="\t")
    required_columns = {
        "gene_id",
        "baseMean",
        "wald_statistic",
        "p_value",
        "adjusted_p_value",
        "log2_fold_change_unshrunk",
        "log2_fold_change_shrunken",
        "log2_fold_change_posterior_sd",
        "absolute_shrunken_log2_fold_change",
        "effect_direction",
        "contrast_marker_rank",
    }
    missing = sorted(required_columns - set(marker_table.columns))
    if missing:
        sys.exit(
            f"[ERROR] {profile} {contrast_id} marker table lacks canonical columns: "
            + ", ".join(missing)
        )
    for numeric_column in (
        "baseMean",
        "adjusted_p_value",
        "p_value",
        "wald_statistic",
        "log2_fold_change_shrunken",
        "log2_fold_change_unshrunk",
        "absolute_shrunken_log2_fold_change",
        "contrast_marker_rank",
    ):
        if numeric_column in marker_table.columns:
            marker_table[numeric_column] = pd.to_numeric(marker_table[numeric_column], errors="coerce")
    marker_table["gene_id_raw"] = marker_table["gene_id"].astype(str)
    marker_table["gene_id"] = marker_table["gene_id"].map(strip_ensg_version)
    collisions = (
        marker_table.groupby("gene_id")["gene_id_raw"]
        .apply(lambda values: sorted(set(values)))
        .reset_index()
    )
    collisions = collisions[collisions["gene_id_raw"].map(len) > 1]
    if not collisions.empty:
        details = "; ".join(
            f"{row.gene_id}: {','.join(row.gene_id_raw)}"
            for row in collisions.itertuples(index=False)
        )
        sys.exit(
            f"[ERROR] Ensembl version-stripping collision in {profile} {contrast_id}: {details}"
        )
    marker_table = marker_table[marker_table["gene_id"].astype(bool)].copy()
    if marker_table["gene_id"].duplicated().any():
        duplicates = sorted(marker_table.loc[marker_table["gene_id"].duplicated(), "gene_id"].unique())
        sys.exit(f"[ERROR] Duplicate cleaned gene IDs in {profile} {contrast_id}: {', '.join(duplicates)}")
    return marker_table


def load_gene_list(marker_gene_list_path: Path) -> list[str]:
    if not marker_gene_list_path.exists():
        sys.exit(f"[ERROR] Missing retained marker list: {marker_gene_list_path}")
    genes: list[str] = []
    with marker_gene_list_path.open() as handle:
        for line in handle:
            gene_id = strip_ensg_version(line.strip())
            if gene_id and not gene_id.startswith("#"):
                genes.append(gene_id)
    return genes


def marker_evidence_ingestion(
    profile: str,
    manifest_path: Path | None,
) -> tuple[pd.DataFrame, Path, pd.DataFrame]:
    resolved_manifest_path = resolve_manifest_path(manifest_path)
    if not resolved_manifest_path.exists():
        sys.exit(f"[ERROR] Missing marker manifest for {profile}: {resolved_manifest_path}")
    manifest = normalise_manifest(pd.read_csv(resolved_manifest_path, sep="\t"), profile)
    deseq2_dir = resolved_manifest_path.parent.parent
    records: list[pd.DataFrame] = []
    for contrast_order, row in enumerate(manifest.itertuples(index=False), start=1):
        contrast_id = str(getattr(row, "contrast_id"))
        marker_evidence_stratum = str(getattr(row, "marker_evidence_stratum"))
        marker_table_path = resolve_relative_path(deseq2_dir, getattr(row, "marker_table_path"))
        marker_gene_list_path = resolve_relative_path(deseq2_dir, getattr(row, "marker_gene_list_path"))
        marker_table = load_retained_marker_table(marker_table_path, profile, contrast_id)
        retained_genes = load_gene_list(marker_gene_list_path)
        table_gene_set = set(marker_table["gene_id"])
        list_gene_set = set(retained_genes)
        manifest_count = int(getattr(row, "n_markers_after_cap"))
        if len(marker_table) != manifest_count or len(retained_genes) != manifest_count or table_gene_set != list_gene_set:
            sys.exit(
                "[ERROR] Marker manifest/table/list disagreement: "
                f"cancer_type={profile} contrast_id={contrast_id} "
                f"manifest_count={manifest_count} table_row_count={len(marker_table)} "
                f"gene_list_count={len(retained_genes)} "
                f"missing_from_table={','.join(sorted(list_gene_set - table_gene_set))} "
                f"missing_from_gene_list={','.join(sorted(table_gene_set - list_gene_set))}"
            )
        marker_table["profile"] = profile
        marker_table["cancer_type"] = str(getattr(row, "cancer_type", profile)).lower()
        marker_table["contrast_id"] = contrast_id
        marker_table["contrast_type"] = str(getattr(row, "contrast_type", ""))
        marker_table["marker_evidence_stratum"] = marker_evidence_stratum
        marker_table["marker_source_class"] = marker_source_class_from_evidence_stratum(
            marker_evidence_stratum
        )
        marker_table["focal_profile_id"] = str(getattr(row, "focal_profile_id", ""))
        marker_table["focal_component_id"] = str(getattr(row, "focal_component_id", ""))
        marker_table["reference_definition"] = str(getattr(row, "reference_definition", ""))
        marker_table["contrast_order"] = int(contrast_order)
        marker_table["source_file"] = str(marker_table_path)
        if len(marker_table) > 0:
            records.append(marker_table)
    if not records:
        evidence = pd.DataFrame(columns=[
            "gene_id", "gene_id_raw", "profile", "cancer_type", "contrast_id",
            "contrast_type", "marker_evidence_stratum", "marker_source_class",
            "focal_profile_id", "focal_component_id", "reference_definition",
            "contrast_order", "source_file", "baseMean", "adjusted_p_value",
            "p_value", "wald_statistic", "log2_fold_change_shrunken",
            "log2_fold_change_unshrunk", "absolute_shrunken_log2_fold_change",
            "contrast_marker_rank",
        ])
    else:
        evidence = pd.concat(records, ignore_index=True)
    print(
        f"[Feature ingestion] {profile}: {manifest['contrast_id'].nunique()} manifest marker lists, "
        f"{evidence['contrast_id'].nunique()} non-empty retained marker lists, "
        f"{evidence['gene_id'].nunique()} retained genes",
        file=sys.stderr,
    )
    return evidence, resolved_manifest_path, manifest


def classify_effect_direction_consistency(log2fc_values: pd.Series) -> str:
    values = finite_numeric(log2fc_values)
    return classify_direction_from_sign_counts(
        int((values > 0).sum()),
        int((values < 0).sum()),
        int((values == 0).sum()),
    )


def classify_direction_from_sign_counts(
    n_positive_effects: int,
    n_negative_effects: int,
    n_zero_effects: int = 0,
) -> str:
    if n_positive_effects > 0 and n_negative_effects == 0:
        return "consistently_upregulated"
    if n_negative_effects > 0 and n_positive_effects == 0:
        return "consistently_downregulated"
    if n_positive_effects == 0 and n_negative_effects == 0 and n_zero_effects > 0:
        return "zero_effect_only"
    return "mixed_direction"


def pan_cancer_direction_label(direction_consistency_class: str) -> str:
    if direction_consistency_class == "consistently_upregulated":
        return "UP"
    if direction_consistency_class == "consistently_downregulated":
        return "DOWN"
    if direction_consistency_class == "mixed_direction":
        return "MIXED"
    return "UNKNOWN"


def evidence_stratum_assignment(evidence: pd.DataFrame) -> pd.DataFrame:
    assigned = evidence.copy()
    assigned["marker_evidence_stratum"] = assigned["marker_evidence_stratum"].fillna("").astype(str)
    assigned["marker_source_class"] = assigned["marker_evidence_stratum"].map(
        marker_source_class_from_evidence_stratum
    )
    return assigned


def within_stratum_gene_evidence_aggregation(
    profile_order: list[str],
    evidence_by_profile: dict[str, pd.DataFrame],
    manifest_tables: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for profile in profile_order:
        evidence = evidence_by_profile[profile]
        manifest = manifest_tables[profile]
        for marker_evidence_stratum in MARKER_EVIDENCE_STRATUM_ORDER:
            stratum_manifest = manifest[
                manifest["marker_evidence_stratum"] == marker_evidence_stratum
            ].copy()
            stratum_evidence = evidence[
                evidence["marker_evidence_stratum"] == marker_evidence_stratum
            ].copy()
            retained_marker_lists = (
                stratum_manifest["contrast_id"]
                .drop_duplicates()
                .astype(str)
                .tolist()
            )
            retained_marker_list_total = len(retained_marker_lists)
            if retained_marker_list_total == 0:
                summary_rows.append(
                    {
                        "method": METHOD_NAME,
                        "cancer_type": profile,
                        "cohort": profile,
                        "marker_evidence_stratum": marker_evidence_stratum,
                        "marker_source_class": marker_source_class_from_evidence_stratum(marker_evidence_stratum),
                        "retained_marker_lists": "",
                        "retained_marker_list_total": 0,
                        "unique_genes_before_selection": 0,
                    }
                )
                continue
            contrast_order = (
                {contrast_id: index for index, contrast_id in enumerate(retained_marker_lists, start=1)}
            )
            if stratum_evidence.empty:
                summary_rows.append(
                    {
                        "method": METHOD_NAME,
                        "cancer_type": profile,
                        "cohort": profile,
                        "marker_evidence_stratum": marker_evidence_stratum,
                        "marker_source_class": marker_source_class_from_evidence_stratum(marker_evidence_stratum),
                        "retained_marker_lists": ";".join(retained_marker_lists),
                        "retained_marker_list_total": int(retained_marker_list_total),
                        "unique_genes_before_selection": 0,
                    }
                )
                continue
            for gene_id, gene_rows in stratum_evidence.groupby("gene_id", sort=True):
                ordered_contrasts = sorted(
                    gene_rows["contrast_id"].astype(str).unique(),
                    key=lambda value: (contrast_order.get(value, 10**9), value),
                )
                ranks = [
                    f"{contrast_id}:{int(rows_for_contrast['contrast_marker_rank'].min())}"
                    for contrast_id, rows_for_contrast in gene_rows.groupby("contrast_id", sort=False)
                ]
                effect_values = finite_numeric(gene_rows["log2_fold_change_shrunken"])
                n_positive_effects = int((effect_values > 0).sum())
                n_negative_effects = int((effect_values < 0).sum())
                n_zero_effects = int((effect_values == 0).sum())
                if n_zero_effects > 0 and n_positive_effects == 0 and n_negative_effects == 0:
                    sys.exit(
                        f"[ERROR] Retained marker has only zero shrunken effects: "
                        f"{profile} {marker_evidence_stratum} {gene_id}"
                    )
                rows.append(
                    {
                        "method": METHOD_NAME,
                        "profile": profile,
                        "cancer_type": profile,
                        "cohort": profile,
                        "marker_evidence_stratum": marker_evidence_stratum,
                        "marker_source_class": marker_source_class_from_evidence_stratum(marker_evidence_stratum),
                        "gene_id": gene_id,
                        "gene_id_raw": join_unique(gene_rows["gene_id_raw"]),
                        "retained_marker_list_count": int(gene_rows["contrast_id"].nunique()),
                        "retained_marker_list_total": int(retained_marker_list_total),
                        "source_count": int(gene_rows["contrast_id"].nunique()),
                        "marker_source_class_source_count": int(retained_marker_list_total),
                        "source_contrast": ";".join(ordered_contrasts),
                        "source_file": join_unique(gene_rows.sort_values(["contrast_order", "contrast_marker_rank"])["source_file"]),
                        "source_marker_ranks": ";".join(ranks),
                        "minimum_adjusted_p_value": safe_min(gene_rows["adjusted_p_value"], default=math.inf),
                        "median_adjusted_p_value": safe_median(gene_rows["adjusted_p_value"], default=math.inf),
                        "minimum_p_value": safe_min(gene_rows["p_value"], default=math.inf),
                        "median_absolute_shrunken_log2fc": safe_median(gene_rows["absolute_shrunken_log2_fold_change"]),
                        "maximum_absolute_shrunken_log2fc": safe_max(gene_rows["absolute_shrunken_log2_fold_change"]),
                        "median_shrunken_log2fc": safe_median(gene_rows["log2_fold_change_shrunken"]),
                        "mean_baseMean": safe_mean(gene_rows["baseMean"]),
                        "median_base_mean": safe_median(gene_rows["baseMean"]),
                        "max_baseMean": safe_max(gene_rows["baseMean"]),
                        "minimum_contrast_marker_rank": int(pd.to_numeric(gene_rows["contrast_marker_rank"], errors="coerce").min()),
                        "n_positive_effects": n_positive_effects,
                        "n_negative_effects": n_negative_effects,
                        "n_zero_effects": n_zero_effects,
                        "direction_consistency_class": classify_direction_from_sign_counts(
                            n_positive_effects,
                            n_negative_effects,
                            n_zero_effects,
                        ),
                    }
                )
            summary_rows.append(
                {
                    "method": METHOD_NAME,
                    "cancer_type": profile,
                    "cohort": profile,
                    "marker_evidence_stratum": marker_evidence_stratum,
                    "marker_source_class": marker_source_class_from_evidence_stratum(marker_evidence_stratum),
                    "retained_marker_lists": ";".join(retained_marker_lists),
                    "retained_marker_list_total": int(retained_marker_list_total),
                    "unique_genes_before_selection": int(stratum_evidence["gene_id"].nunique()),
                }
            )
    if not rows:
        sys.exit("[ERROR] No within-stratum gene evidence was produced")
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def recurrence_classification(
    aggregated_gene_evidence: pd.DataFrame,
) -> pd.DataFrame:
    classified = aggregated_gene_evidence.copy()
    if (
        (classified["retained_marker_list_count"] < 1)
        | (classified["retained_marker_list_total"] < 1)
        | (classified["retained_marker_list_count"] > classified["retained_marker_list_total"])
    ).any():
        bad = classified.loc[
            (classified["retained_marker_list_count"] < 1)
            | (classified["retained_marker_list_total"] < 1)
            | (classified["retained_marker_list_count"] > classified["retained_marker_list_total"]),
            ["cancer_type", "marker_evidence_stratum", "gene_id", "retained_marker_list_count", "retained_marker_list_total"],
        ]
        sys.exit(
            "[ERROR] Invalid recurrence context encountered:\n"
            + bad.to_string(index=False)
        )
    recurrent = classified["retained_marker_list_count"] >= RECURRENCE_MIN_COUNT
    singleton_candidate = (
        (classified["retained_marker_list_total"] == 1)
        & (classified["retained_marker_list_count"] == 1)
    )
    non_recurrent_candidate = (
        (classified["retained_marker_list_total"] > 1)
        & (classified["retained_marker_list_count"] == 1)
    )
    classified["candidate_pool_type"] = pd.NA
    classified.loc[recurrent, "candidate_pool_type"] = "recurrent"
    classified.loc[singleton_candidate, "candidate_pool_type"] = "singleton_candidate"
    classified.loc[non_recurrent_candidate, "candidate_pool_type"] = "non_recurrent_candidate"
    unclassified = classified[classified["candidate_pool_type"].isna()]
    if not unclassified.empty:
        sys.exit(
            "[ERROR] Unexpected recurrence state under invariant r>=2/r==1 method:\n"
            + unclassified[
                [
                    "cancer_type",
                    "marker_evidence_stratum",
                    "gene_id",
                    "retained_marker_list_count",
                    "retained_marker_list_total",
                ]
            ].to_string(index=False)
        )
    classified["candidate_pool"] = classified["candidate_pool_type"]
    classified["recurrence_min_count"] = RECURRENCE_MIN_COUNT
    return classified


def candidate_pool_construction(classified_gene_evidence: pd.DataFrame) -> pd.DataFrame:
    candidate_pool_evidence = classified_gene_evidence[
        classified_gene_evidence["candidate_pool_type"].isin(
            ["recurrent", "singleton_candidate", "non_recurrent_candidate"]
        )
    ].copy()
    if candidate_pool_evidence.empty:
        sys.exit("[ERROR] Candidate-pool construction produced no genes")
    return candidate_pool_evidence


def empirical_quantile_threshold_calculation(
    candidate_pool_evidence: pd.DataFrame,
    adjusted_p_value_quantile: float,
    absolute_shrunken_log2fc_quantile: float,
    expression_quantile: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotated = candidate_pool_evidence.copy()
    threshold_columns = {
        "adjusted_p_value_quantile_threshold": math.nan,
        "absolute_shrunken_log2fc_quantile_threshold": math.nan,
        "expression_quantile_threshold": math.nan,
    }
    for column, default in threshold_columns.items():
        annotated[column] = default
    rows: list[dict[str, object]] = []
    empirical_candidates = annotated[
        annotated["candidate_pool_type"].isin(["singleton_candidate", "non_recurrent_candidate"])
    ]
    for (cancer_type, marker_evidence_stratum, candidate_pool_type), pool in empirical_candidates.groupby(
        ["cancer_type", "marker_evidence_stratum", "candidate_pool_type"], sort=True
    ):
        row_index = pool.index
        adjusted_p_value_quantile_threshold = safe_quantile(
            pool["minimum_adjusted_p_value"], adjusted_p_value_quantile
        )
        absolute_shrunken_log2fc_quantile_threshold = safe_quantile(
            pool["median_absolute_shrunken_log2fc"], absolute_shrunken_log2fc_quantile
        )
        expression_quantile_threshold = safe_quantile(pool["median_base_mean"], expression_quantile)
        for label, threshold in (
            ("adjusted_p_value_quantile_threshold", adjusted_p_value_quantile_threshold),
            ("absolute_shrunken_log2fc_quantile_threshold", absolute_shrunken_log2fc_quantile_threshold),
            ("expression_quantile_threshold", expression_quantile_threshold),
        ):
            if not math.isfinite(float(threshold)):
                sys.exit(
                    "[ERROR] Non-finite empirical threshold for "
                    f"{cancer_type}/{marker_evidence_stratum}/{candidate_pool_type}: {label}"
                )
        annotated.loc[row_index, "adjusted_p_value_quantile_threshold"] = (
            adjusted_p_value_quantile_threshold
        )
        annotated.loc[row_index, "absolute_shrunken_log2fc_quantile_threshold"] = (
            absolute_shrunken_log2fc_quantile_threshold
        )
        annotated.loc[row_index, "expression_quantile_threshold"] = expression_quantile_threshold
        rows.append(
            {
                "method": METHOD_NAME,
                "cancer_type": cancer_type,
                "cohort": cancer_type,
                "marker_evidence_stratum": marker_evidence_stratum,
                "marker_source_class": marker_source_class_from_evidence_stratum(marker_evidence_stratum),
                "candidate_pool_type": candidate_pool_type,
                "candidate_pool": candidate_pool_type,
                "number_of_candidate_genes": int(pool["gene_id"].nunique()),
                "adjusted_p_value_quantile": adjusted_p_value_quantile,
                "adjusted_p_value_quantile_threshold": adjusted_p_value_quantile_threshold,
                "absolute_shrunken_log2fc_quantile": absolute_shrunken_log2fc_quantile,
                "absolute_shrunken_log2fc_quantile_threshold": absolute_shrunken_log2fc_quantile_threshold,
                "expression_quantile": expression_quantile,
                "expression_quantile_threshold": expression_quantile_threshold,
            }
        )
    threshold_table_columns = [
        "method",
        "cancer_type",
        "cohort",
        "marker_evidence_stratum",
        "marker_source_class",
        "candidate_pool_type",
        "candidate_pool",
        "number_of_candidate_genes",
        "adjusted_p_value_quantile",
        "adjusted_p_value_quantile_threshold",
        "absolute_shrunken_log2fc_quantile",
        "absolute_shrunken_log2fc_quantile_threshold",
        "expression_quantile",
        "expression_quantile_threshold",
    ]
    return annotated, pd.DataFrame(rows, columns=threshold_table_columns)


def effect_direction_consistency_classification(candidate_pool_evidence: pd.DataFrame) -> pd.DataFrame:
    classified = candidate_pool_evidence.copy()
    classified["direction_consistent"] = classified["direction_consistency_class"].isin(
        ["consistently_upregulated", "consistently_downregulated"]
    )
    classified["direction"] = classified["direction_consistency_class"].map(pan_cancer_direction_label)
    return classified


def evaluate_candidate_acceptance(candidate_pool_evidence: pd.DataFrame) -> pd.DataFrame:
    evaluated = candidate_pool_evidence.copy()
    for column in (
        "passes_statistical_evidence",
        "passes_effect_magnitude",
        "passes_expression_evidence",
        "passes_candidate_acceptance",
    ):
        evaluated[column] = pd.NA

    candidate_rows = evaluated["candidate_pool_type"].isin(
        ["singleton_candidate", "non_recurrent_candidate"]
    )
    if candidate_rows.any():
        required_thresholds = [
            "adjusted_p_value_quantile_threshold",
            "absolute_shrunken_log2fc_quantile_threshold",
            "expression_quantile_threshold",
        ]
        bad_thresholds = evaluated.loc[
            candidate_rows,
            required_thresholds,
        ].apply(lambda col: ~pd.to_numeric(col, errors="coerce").map(math.isfinite))
        if bad_thresholds.any(axis=None):
            bad = evaluated.loc[
                candidate_rows & bad_thresholds.any(axis=1),
                ["cancer_type", "marker_evidence_stratum", "candidate_pool_type", "gene_id"],
            ]
            sys.exit(
                "[ERROR] Candidate rows lack finite empirical thresholds:\n"
                + bad.to_string(index=False)
            )
        evaluated.loc[candidate_rows, "passes_statistical_evidence"] = (
            evaluated.loc[candidate_rows, "minimum_adjusted_p_value"]
            <= evaluated.loc[candidate_rows, "adjusted_p_value_quantile_threshold"]
        )
        evaluated.loc[candidate_rows, "passes_effect_magnitude"] = (
            evaluated.loc[candidate_rows, "median_absolute_shrunken_log2fc"]
            >= evaluated.loc[candidate_rows, "absolute_shrunken_log2fc_quantile_threshold"]
        )
        evaluated.loc[candidate_rows, "passes_expression_evidence"] = (
            evaluated.loc[candidate_rows, "median_base_mean"]
            >= evaluated.loc[candidate_rows, "expression_quantile_threshold"]
        )
        evaluated.loc[candidate_rows, "passes_candidate_acceptance"] = (
            evaluated.loc[candidate_rows, "passes_statistical_evidence"].astype(bool)
            & evaluated.loc[candidate_rows, "passes_effect_magnitude"].astype(bool)
            & evaluated.loc[candidate_rows, "passes_expression_evidence"].astype(bool)
        )
    return evaluated


def selection_basis_from_candidate_pool(candidate_pool_type: str) -> str:
    if candidate_pool_type == "recurrent":
        return "recurrent"
    if candidate_pool_type == "singleton_candidate":
        return "accepted_singleton"
    if candidate_pool_type == "non_recurrent_candidate":
        return "accepted_non_recurrent"
    sys.exit(f"[ERROR] Unknown candidate_pool_type for selection: {candidate_pool_type}")


def select_pan_cancer_evidence_rows(candidate_pool_evidence: pd.DataFrame) -> pd.DataFrame:
    recurrent = candidate_pool_evidence[
        candidate_pool_evidence["candidate_pool_type"] == "recurrent"
    ].copy()
    accepted_candidates = candidate_pool_evidence[
        candidate_pool_evidence["candidate_pool_type"].isin(
            ["singleton_candidate", "non_recurrent_candidate"]
        )
        & (candidate_pool_evidence["passes_candidate_acceptance"] == True)
    ].copy()
    selected = pd.concat([recurrent, accepted_candidates], ignore_index=True)
    if selected.empty:
        sys.exit("[ERROR] Revised pan-cancer selection produced no evidence rows")
    selected["selection_basis"] = selected["candidate_pool_type"].map(
        selection_basis_from_candidate_pool
    )
    selected["feature_class"] = selected["selection_basis"]
    return evidence_row_ordering(selected)


def evidence_row_ordering(selected_evidence_rows: pd.DataFrame) -> pd.DataFrame:
    ordered = selected_evidence_rows.copy()
    ordered["evidence_tier"] = ordered["candidate_pool_type"].map(CANDIDATE_POOL_ORDER)
    ordered["selection_basis_priority"] = ordered["selection_basis"].map(
        selection_basis_sort_value
    )
    ordered = ordered.sort_values(
        [
            "selection_basis_priority",
            "retained_marker_list_count",
            "minimum_adjusted_p_value",
            "median_absolute_shrunken_log2fc",
            "median_base_mean",
            "minimum_contrast_marker_rank",
            "cancer_type",
            "marker_evidence_stratum",
            "gene_id",
        ],
        ascending=[True, False, True, False, False, True, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ordered["within_candidate_pool_rank"] = (
        ordered.groupby(["cancer_type", "marker_evidence_stratum", "candidate_pool_type"], sort=False)
        .cumcount()
        + 1
    )
    ordered["rank_within_candidate_pool"] = ordered["within_candidate_pool_rank"]
    ordered["rank_within_cohort_marker_source_class"] = (
        ordered.groupby(["cancer_type", "marker_evidence_stratum"], sort=False)
        .cumcount()
        + 1
    )
    ordered["minimum_adjusted_p_value_rank_within_pool"] = (
        ordered.groupby(["cancer_type", "marker_evidence_stratum", "candidate_pool_type"])[
            "minimum_adjusted_p_value"
        ]
        .rank(method="min", ascending=True)
        .astype("Int64")
    )
    ordered["ranking_component_key"] = (
        ordered["cancer_type"].astype(str)
        + "|"
        + ordered["marker_evidence_stratum"].astype(str)
        + "|"
        + ordered["candidate_pool_type"].astype(str)
        + "|"
        + ordered["gene_id"].astype(str)
    )
    return ordered


def gene_level_provenance_annotation(gene_evidence_rows: pd.DataFrame) -> pd.DataFrame:
    annotated = gene_evidence_rows.copy()
    annotated["has_anchor_associated_marker_evidence"] = annotated["marker_evidence_stratum"].map(
        lambda value: "anchor_associated" in {part.strip() for part in str(value).split(";")}
    )
    annotated["has_isolate_associated_marker_evidence"] = annotated["marker_evidence_stratum"].map(
        lambda value: "isolate_associated" in {part.strip() for part in str(value).split(";")}
    )
    annotated["marker_evidence_provenance"] = "isolate_associated"
    annotated.loc[annotated["has_anchor_associated_marker_evidence"], "marker_evidence_provenance"] = (
        "anchor_associated"
    )
    dual = annotated["has_anchor_associated_marker_evidence"] & annotated[
        "has_isolate_associated_marker_evidence"
    ]
    annotated.loc[dual, "marker_evidence_provenance"] = "dual_graph_derived_evidence"
    return annotated


def build_gene_evidence_table(selected_evidence_rows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (cancer_type, gene_id), gene_rows in selected_evidence_rows.groupby(
        ["cancer_type", "gene_id"], sort=False
    ):
        selection_basis_classes = sorted(
            set(gene_rows["selection_basis"]),
            key=selection_basis_sort_value,
        )
        marker_evidence_strata = sorted(
            set(gene_rows["marker_evidence_stratum"]), key=stratum_sort_value
        )
        n_positive_effects = int(gene_rows["n_positive_effects"].sum())
        n_negative_effects = int(gene_rows["n_negative_effects"].sum())
        n_zero_effects = int(gene_rows["n_zero_effects"].sum())
        direction_consistency_class = classify_direction_from_sign_counts(
            n_positive_effects,
            n_negative_effects,
            n_zero_effects,
        )
        rows.append(
            {
                "method": METHOD_NAME,
                "gene_id": gene_id,
                "gene_id_raw": join_unique(gene_rows["gene_id_raw"]),
                "cancer_type": cancer_type,
                "cohort": cancer_type,
                "profile": cancer_type,
                "marker_evidence_stratum": ";".join(marker_evidence_strata),
                "marker_source_class": join_unique(gene_rows["marker_source_class"]),
                "retained_marker_source_classes": join_unique(gene_rows["marker_source_class"]),
                "candidate_pool_type": join_unique(gene_rows["candidate_pool_type"]),
                "candidate_pool": join_unique(gene_rows["candidate_pool_type"]),
                "selection_basis": selection_basis_classes[0],
                "selection_basis_classes": ";".join(selection_basis_classes),
                "feature_class": selection_basis_classes[0],
                "evidence_tiers": join_unique(gene_rows["evidence_tier"]),
                "source_contrast": join_unique(gene_rows["source_contrast"]),
                "source_file": join_unique(gene_rows["source_file"]),
                "source_marker_ranks": join_unique(gene_rows["source_marker_ranks"]),
                "ranking_component_keys": join_unique(gene_rows["ranking_component_key"]),
                "rank_within_candidate_pool": join_unique(gene_rows["rank_within_candidate_pool"]),
                "rank_within_cohort_marker_source_class": join_unique(gene_rows["rank_within_cohort_marker_source_class"]),
                "minimum_adjusted_p_value_rank_within_pool": join_unique(gene_rows["minimum_adjusted_p_value_rank_within_pool"]),
                "passes_candidate_acceptance": join_unique(gene_rows["passes_candidate_acceptance"]),
                "retained_marker_list_count": int(gene_rows["retained_marker_list_count"].max()),
                "retained_marker_list_total": int(gene_rows["retained_marker_list_total"].max()),
                "source_count": int(gene_rows["retained_marker_list_count"].max()),
                "marker_source_class_source_count": int(gene_rows["retained_marker_list_total"].max()),
                "n_retained_contrasts": int(gene_rows["retained_marker_list_count"].sum()),
                "n_anchor_contrasts": int(gene_rows.loc[gene_rows["marker_evidence_stratum"] == "anchor_associated", "retained_marker_list_count"].sum()),
                "n_isolate_contrasts": int(gene_rows.loc[gene_rows["marker_evidence_stratum"] == "isolate_associated", "retained_marker_list_count"].sum()),
                "n_positive_effects": n_positive_effects,
                "n_negative_effects": n_negative_effects,
                "n_zero_effects": n_zero_effects,
                "minimum_adjusted_p_value": safe_min(gene_rows["minimum_adjusted_p_value"], default=math.inf),
                "min_padj": safe_min(gene_rows["minimum_adjusted_p_value"], default=math.inf),
                "median_padj": safe_median(gene_rows["median_adjusted_p_value"], default=math.inf),
                "median_absolute_shrunken_log2fc": safe_median(gene_rows["median_absolute_shrunken_log2fc"]),
                "median_abs_log2FC": safe_median(gene_rows["median_absolute_shrunken_log2fc"]),
                "max_abs_log2FC": safe_max(gene_rows["maximum_absolute_shrunken_log2fc"]),
                "median_shrunken_log2fc": safe_median(gene_rows["median_shrunken_log2fc"]),
                "median_log2FC": safe_median(gene_rows["median_shrunken_log2fc"]),
                "mean_baseMean": safe_mean(gene_rows["mean_baseMean"]),
                "median_base_mean": safe_median(gene_rows["median_base_mean"]),
                "median_baseMean": safe_median(gene_rows["median_base_mean"]),
                "max_baseMean": safe_max(gene_rows["max_baseMean"]),
                "minimum_contrast_marker_rank": int(gene_rows["minimum_contrast_marker_rank"].min()),
                "direction_consistency_class": direction_consistency_class,
            }
        )
    evidence = pd.DataFrame(rows)
    evidence["direction"] = evidence["direction_consistency_class"].map(pan_cancer_direction_label)
    evidence["direction_consistent"] = evidence["direction_consistency_class"].isin(
        ["consistently_upregulated", "consistently_downregulated"]
    )
    return gene_level_provenance_annotation(evidence)


def collapse_gene_evidence_across_cancer_types_and_strata(
    profile_order: list[str],
    gene_evidence: pd.DataFrame,
    annotation: dict[str, dict[str, str]],
    remove_ribo_mt: bool,
) -> tuple[pd.DataFrame, int, set[str]]:
    profile_index = {profile: index for index, profile in enumerate(profile_order)}
    candidate = gene_evidence.copy()
    candidate["_profile_order"] = candidate["cancer_type"].map(lambda value: profile_index.get(value, 10**6))
    candidate["_feature_class_priority"] = candidate["feature_class"].map(selection_basis_sort_value)
    candidate = candidate.sort_values(
        [
            "_feature_class_priority",
            "retained_marker_list_count",
            "minimum_adjusted_p_value",
            "median_absolute_shrunken_log2fc",
            "median_base_mean",
            "minimum_contrast_marker_rank",
            "_profile_order",
            "gene_id",
        ],
        ascending=[True, False, True, False, False, True, True, True],
        kind="mergesort",
    )
    ordered_genes: list[str] = []
    seen: set[str] = set()
    for gene_id in candidate["gene_id"]:
        if gene_id in seen:
            continue
        seen.add(gene_id)
        ordered_genes.append(gene_id)
    removed_ribo_mt: set[str] = set()
    if remove_ribo_mt:
        removed_ribo_mt = {
            gene_id for gene_id in ordered_genes
            if is_ribo_mt_or_histone(gene_id, annotation)
        }
        ordered_genes = [gene_id for gene_id in ordered_genes if gene_id not in removed_ribo_mt]
    rows: list[dict[str, object]] = []
    for selection_rank, gene_id in enumerate(ordered_genes, start=1):
        evidence = gene_evidence[gene_evidence["gene_id"] == gene_id].copy()
        evidence["_profile_order"] = evidence["cancer_type"].map(lambda value: profile_index.get(value, 10**6))
        evidence["_feature_class_priority"] = evidence["feature_class"].map(selection_basis_sort_value)
        evidence = evidence.sort_values(
            [
                "_feature_class_priority",
                "retained_marker_list_count",
                "minimum_adjusted_p_value",
                "median_absolute_shrunken_log2fc",
                "median_base_mean",
                "minimum_contrast_marker_rank",
                "_profile_order",
            ],
            ascending=[True, False, True, False, False, True, True],
            kind="mergesort",
        )
        first_evidence = evidence.iloc[0]
        annotations = annotation.get(gene_id, {})
        all_selection_basis_classes = sorted(
            {
                item
                for value in evidence["selection_basis_classes"]
                for item in str(value).split(";")
                if item
            },
            key=selection_basis_sort_value,
        )
        has_anchor = evidence["has_anchor_associated_marker_evidence"].any()
        has_isolate = evidence["has_isolate_associated_marker_evidence"].any()
        marker_evidence_provenance = (
            "dual_graph_derived_evidence"
            if has_anchor and has_isolate
            else "anchor_associated" if has_anchor else "isolate_associated"
        )
        n_positive_effects = int(evidence["n_positive_effects"].sum())
        n_negative_effects = int(evidence["n_negative_effects"].sum())
        n_zero_effects = int(evidence["n_zero_effects"].sum())
        direction_consistency_class = classify_direction_from_sign_counts(
            n_positive_effects,
            n_negative_effects,
            n_zero_effects,
        )
        rows.append(
            {
                "method": METHOD_NAME,
                "gene_id": gene_id,
                "clean_gene_id": gene_id,
                "gene_id_raw": join_unique(evidence["gene_id_raw"]),
                "gene_symbol": annotations.get("gene_symbol", ""),
                "gene_name": annotations.get("gene_name", ""),
                "direction_consistency_class": direction_consistency_class,
                "direction": pan_cancer_direction_label(direction_consistency_class),
                "direction_consistent": direction_consistency_class in {
                    "consistently_upregulated",
                    "consistently_downregulated",
                },
                "feature_class": all_selection_basis_classes[0],
                "selection_basis": all_selection_basis_classes[0],
                "selection_basis_classes": ";".join(all_selection_basis_classes),
                "owner_profile": str(first_evidence["cancer_type"]),
                "cohort": join_unique(evidence["cancer_type"]),
                "cancer_type": join_unique(evidence["cancer_type"]),
                "profiles_present": join_unique(evidence["cancer_type"], sep=","),
                "n_profiles": int(evidence["cancer_type"].nunique()),
                "marker_evidence_stratum": join_unique(evidence["marker_evidence_stratum"]),
                "marker_source_class": join_unique(evidence["marker_source_class"]),
                "retained_marker_source_classes": join_unique(evidence["retained_marker_source_classes"]),
                "has_anchor_associated_marker_evidence": bool(has_anchor),
                "has_isolate_associated_marker_evidence": bool(has_isolate),
                "marker_evidence_provenance": marker_evidence_provenance,
                "n_positive_effects": n_positive_effects,
                "n_negative_effects": n_negative_effects,
                "n_zero_effects": n_zero_effects,
                "source_contrast": join_unique(evidence["source_contrast"]),
                "source_file": join_unique(evidence["source_file"]),
                "source_marker_ranks": join_unique(evidence["source_marker_ranks"]),
                "marker_source_recurrence_detail": join_unique(
                    [
                        f"{row.marker_evidence_stratum}:{int(row.retained_marker_list_count)}/{int(row.retained_marker_list_total)}"
                        for row in evidence.itertuples(index=False)
                    ]
                ),
                "ranking_component_keys": join_unique(evidence["ranking_component_keys"]),
                "rank_within_candidate_pool": join_unique(evidence["rank_within_candidate_pool"]),
                "minimum_adjusted_p_value_rank_within_pool": join_unique(evidence["minimum_adjusted_p_value_rank_within_pool"]),
                "passes_candidate_acceptance": join_unique(evidence["passes_candidate_acceptance"]),
                "retained_marker_list_count": int(evidence["retained_marker_list_count"].max()),
                "retained_marker_list_total": int(evidence["retained_marker_list_total"].max()),
                "source_count": int(evidence["retained_marker_list_count"].max()),
                "marker_source_class_source_count": int(evidence["retained_marker_list_total"].max()),
                "n_retained_contrasts": int(evidence["n_retained_contrasts"].sum()),
                "n_anchor_contrasts": int(evidence["n_anchor_contrasts"].sum()),
                "n_isolate_contrasts": int(evidence["n_isolate_contrasts"].sum()),
                "minimum_adjusted_p_value": safe_min(evidence["minimum_adjusted_p_value"], default=math.inf),
                "min_padj": safe_min(evidence["min_padj"], default=math.inf),
                "median_padj": safe_median(evidence["median_padj"], default=math.inf),
                "median_absolute_shrunken_log2fc": safe_median(evidence["median_absolute_shrunken_log2fc"]),
                "median_abs_log2FC": safe_median(evidence["median_abs_log2FC"]),
                "max_abs_log2FC": safe_max(evidence["max_abs_log2FC"]),
                "median_shrunken_log2fc": safe_median(evidence["median_shrunken_log2fc"]),
                "median_log2FC": safe_median(evidence["median_log2FC"]),
                "mean_baseMean": safe_mean(evidence["mean_baseMean"]),
                "median_base_mean": safe_median(evidence["median_base_mean"]),
                "median_baseMean": safe_median(evidence["median_baseMean"]),
                "max_baseMean": safe_max(evidence["max_baseMean"]),
                "minimum_contrast_marker_rank": int(evidence["minimum_contrast_marker_rank"].min()),
                "selection_rank": selection_rank,
                "source_rank": pd.NA,
            }
        )
    features = pd.DataFrame(rows)
    duplicate_gene_count = int(len(candidate) - candidate["gene_id"].nunique())
    return features, duplicate_gene_count, removed_ribo_mt


def selected_pan_cancer_panel_export(
    outdir: Path,
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selected_evidence_rows: pd.DataFrame,
    candidate_pool_evidence: pd.DataFrame,
    threshold_table: pd.DataFrame,
    family_summary: pd.DataFrame,
    profile_order: list[str],
    manifest_paths: dict[str, Path],
    manifest_tables: dict[str, pd.DataFrame],
    args: argparse.Namespace,
    duplicate_gene_count: int,
    removed_ribo_mt: set[str],
) -> list[dict[str, object]]:
    outdir.mkdir(parents=True, exist_ok=True)
    audit_output_prefix = args.audit_output_prefix
    features_path = outdir / "pan_cancer_features.tsv"
    clean_path = outdir / "pan_cancer_features_clean.txt"
    features.to_csv(features_path, sep="\t", index=False)
    write_gene_list(clean_path, features["gene_id"].tolist())
    for direction in ["UP", "DOWN", "MIXED"]:
        write_gene_list(
            outdir / f"pan_cancer_features.{direction}.txt",
            features.loc[features["direction"] == direction, "gene_id"].tolist(),
        )
    (outdir / "pan_cancer_features_done.txt").write_text("done\n")
    gene_evidence.to_csv(outdir / "pan_cancer_feature_gene_evidence.tsv", sep="\t", index=False)
    candidate_pool_evidence.to_csv(
        outdir / f"{audit_output_prefix}_candidate_pool_evidence.tsv",
        sep="\t",
        index=False,
    )
    threshold_table.to_csv(
        outdir / f"{audit_output_prefix}_empirical_quantile_thresholds.tsv",
        sep="\t",
        index=False,
    )
    candidate_pool_evidence[
        candidate_pool_evidence["candidate_pool_type"].isin(
            ["singleton_candidate", "non_recurrent_candidate"]
        )
    ].to_csv(
        outdir / f"{audit_output_prefix}_candidate_acceptance.tsv",
        sep="\t",
        index=False,
    )
    selected_evidence_rows.to_csv(
        outdir / f"{audit_output_prefix}_selected_evidence_rows.tsv",
        sep="\t",
        index=False,
    )
    build_cohort_summary(profile_order, features, gene_evidence, selected_evidence_rows, family_summary).to_csv(
        outdir / f"{audit_output_prefix}_by_cohort.tsv", sep="\t", index=False
    )
    build_marker_evidence_stratum_summary(features, selected_evidence_rows, family_summary).to_csv(
        outdir / f"{audit_output_prefix}_by_marker_source_class.tsv",
        sep="\t",
        index=False,
    )
    build_feature_class_summary(features).to_csv(
        outdir / f"{audit_output_prefix}_by_feature_class.tsv",
        sep="\t",
        index=False,
    )
    validation_rows = validate_outputs(
        features,
        gene_evidence,
        selected_evidence_rows,
        candidate_pool_evidence,
        threshold_table,
        removed_ribo_mt,
        outdir,
    )
    pd.DataFrame(validation_rows).to_csv(
        outdir / f"{audit_output_prefix}_validation.tsv",
        sep="\t",
        index=False,
    )
    build_summary_table(
        features,
        gene_evidence,
        selected_evidence_rows,
        candidate_pool_evidence,
        family_summary,
        duplicate_gene_count,
        removed_ribo_mt,
        adjusted_p_value_quantile=args.adjusted_p_value_quantile,
        absolute_shrunken_log2fc_quantile=args.absolute_shrunken_log2fc_quantile,
        expression_quantile=args.expression_quantile,
    ).to_csv(outdir / "pan_cancer_feature_build_summary.tsv", sep="\t", index=False)
    write_run_manifest(outdir, manifest_paths, manifest_tables, audit_output_prefix)
    write_active_directory_manifest(outdir, audit_output_prefix)
    if any(row["status"] != "PASS" for row in validation_rows):
        for row in validation_rows:
            if row["status"] != "PASS":
                print(f"[VALIDATION FAIL] {row['validation_check']}: {row['details']}", file=sys.stderr)
        sys.exit("[ERROR] graph-derived pan-cancer feature selection failed validation")
    return validation_rows


def explode_count(df: pd.DataFrame, column: str, id_column: str = "gene_id") -> pd.DataFrame:
    counts: dict[str, set[str]] = defaultdict(set)
    if df.empty or column not in df.columns:
        return pd.DataFrame(columns=["value", "n_genes"])
    for row in df[[id_column, column]].itertuples(index=False):
        for value in str(getattr(row, column)).split(";"):
            value = value.strip()
            if value:
                counts[value].add(getattr(row, id_column))
    return pd.DataFrame(
        [{"value": key, "n_genes": len(values)} for key, values in sorted(counts.items())]
    )


def build_cohort_summary(
    profile_order: list[str],
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selected_evidence_rows: pd.DataFrame,
    family_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for profile in profile_order:
        evidence = gene_evidence[gene_evidence["cancer_type"] == profile]
        final_owner = features[features["owner_profile"] == profile]
        family = family_summary[family_summary["cancer_type"] == profile]
        selected = selected_evidence_rows[selected_evidence_rows["cancer_type"] == profile]
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": profile,
                "cancer_type": profile,
                "n_anchor_sources": int(family.loc[family["marker_evidence_stratum"] == "anchor_associated", "retained_marker_list_total"].sum()),
                "n_isolate_sources": int(family.loc[family["marker_evidence_stratum"] == "isolate_associated", "retained_marker_list_total"].sum()),
                "anchor_selected_genes_before_global_dedup": int(selected.loc[selected["marker_evidence_stratum"] == "anchor_associated", "gene_id"].nunique()),
                "isolate_selected_genes_before_global_dedup": int(selected.loc[selected["marker_evidence_stratum"] == "isolate_associated", "gene_id"].nunique()),
                "recurrent_selected_rows": int((selected["selection_basis"] == "recurrent").sum()),
                "accepted_singleton_rows": int((selected["selection_basis"] == "accepted_singleton").sum()),
                "accepted_non_recurrent_rows": int((selected["selection_basis"] == "accepted_non_recurrent").sum()),
                "selected_genes_before_global_dedup": int(evidence["gene_id"].nunique()),
                "final_owner_profile_genes_after_global_dedup": int(final_owner["gene_id"].nunique()),
                "final_genes_with_cohort_evidence": int(evidence["gene_id"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def build_marker_evidence_stratum_summary(
    features: pd.DataFrame,
    selected_evidence_rows: pd.DataFrame,
    family_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in family_summary.itertuples(index=False):
        selected = selected_evidence_rows[
            (selected_evidence_rows["cancer_type"] == row.cancer_type)
            & (selected_evidence_rows["marker_evidence_stratum"] == row.marker_evidence_stratum)
        ]
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": row.cancer_type,
                "cancer_type": row.cancer_type,
                "marker_evidence_stratum": row.marker_evidence_stratum,
                "marker_source_class": row.marker_source_class,
                "source_contrasts": row.retained_marker_lists,
                "n_source_lists_in_marker_source_class": int(row.retained_marker_list_total),
                "retained_marker_list_total": int(row.retained_marker_list_total),
                "unique_genes_before_selection": int(row.unique_genes_before_selection),
                "selected_genes_before_global_dedup": int(selected["gene_id"].nunique()),
                "recurrent_selected_rows": int((selected["selection_basis"] == "recurrent").sum()) if not selected.empty else 0,
                "accepted_singleton_rows": int((selected["selection_basis"] == "accepted_singleton").sum()) if not selected.empty else 0,
                "accepted_non_recurrent_rows": int((selected["selection_basis"] == "accepted_non_recurrent").sum()) if not selected.empty else 0,
            }
        )
    for row in explode_count(features, "marker_source_class").itertuples(index=False):
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": "ALL_FINAL",
                "cancer_type": "ALL_FINAL",
                "marker_evidence_stratum": (
                    "anchor_associated" if row.value == "anchor" else "isolate_associated"
                ),
                "marker_source_class": row.value,
                "source_contrasts": "",
                "n_source_lists_in_marker_source_class": pd.NA,
                "retained_marker_list_total": pd.NA,
                "unique_genes_before_selection": pd.NA,
                "selected_genes_before_global_dedup": int(row.n_genes),
                "recurrent_selected_rows": pd.NA,
                "accepted_singleton_rows": pd.NA,
                "accepted_non_recurrent_rows": pd.NA,
            }
        )
    return pd.DataFrame(rows)


def build_feature_class_summary(features: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"method": METHOD_NAME, "feature_class": row.value, "n_genes": int(row.n_genes)}
            for row in explode_count(features, "selection_basis_classes").itertuples(index=False)
        ]
    )


def validate_outputs(
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selected_evidence_rows: pd.DataFrame,
    candidate_pool_evidence: pd.DataFrame,
    threshold_table: pd.DataFrame,
    removed_ribo_mt: set[str],
    outdir: Path,
) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []

    def add(name: str, passed: bool, details: str) -> None:
        checks.append({"validation_check": name, "status": "PASS" if passed else "FAIL", "details": details})

    def context_key(df: pd.DataFrame) -> set[tuple[str, str, str, str]]:
        if df.empty:
            return set()
        return {
            (
                str(row.cancer_type),
                str(row.marker_evidence_stratum),
                str(row.gene_id),
                str(row.candidate_pool_type),
            )
            for row in df.itertuples(index=False)
        }

    add("no_duplicate_clean_ensembl_ids", not features["gene_id"].duplicated().any(), f"duplicates={int(features['gene_id'].duplicated().sum())}")
    add("all_final_genes_have_non_empty_provenance", features["marker_evidence_provenance"].map(lambda value: bool(str(value).strip())).all(), "all marker_evidence_provenance fields are non-empty")
    recurrent = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "recurrent"]
    recurrent_ok = recurrent.empty or (
        recurrent["retained_marker_list_count"] >= RECURRENCE_MIN_COUNT
    ).all()
    add("all_recurrent_rows_have_recurrence_at_least_two", bool(recurrent_ok), f"recurrent_rows={len(recurrent)} recurrence_min_count={RECURRENCE_MIN_COUNT}")
    singleton = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "singleton_candidate"]
    singleton_ok = singleton.empty or (
        (singleton["retained_marker_list_total"] == 1)
        & (singleton["retained_marker_list_count"] == 1)
    ).all()
    add("singleton_candidate_rows_have_one_list_and_one_occurrence", bool(singleton_ok), f"singleton_rows={len(singleton)}")
    non_recurrent = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "non_recurrent_candidate"]
    non_recurrent_ok = non_recurrent.empty or (
        (non_recurrent["retained_marker_list_total"] > 1)
        & (non_recurrent["retained_marker_list_count"] == 1)
    ).all()
    add("non_recurrent_candidate_rows_have_multi_list_context_and_one_occurrence", bool(non_recurrent_ok), f"non_recurrent_rows={len(non_recurrent)}")
    candidate_rows = candidate_pool_evidence[
        candidate_pool_evidence["candidate_pool_type"].isin(
            ["singleton_candidate", "non_recurrent_candidate"]
        )
    ]
    add("no_zero_recurrence_row_in_candidate_pools", candidate_rows.empty or (candidate_rows["retained_marker_list_count"] == 1).all(), f"candidate_rows={len(candidate_rows)}")
    accepted_singleton = singleton[singleton["passes_candidate_acceptance"] == True]
    accepted_singleton_ok = accepted_singleton.empty or (
        (accepted_singleton["passes_statistical_evidence"] == True)
        & (accepted_singleton["passes_effect_magnitude"] == True)
        & (accepted_singleton["passes_expression_evidence"] == True)
        & (accepted_singleton["passes_candidate_acceptance"] == True)
    ).all()
    add("accepted_singleton_rows_pass_all_candidate_criteria", bool(accepted_singleton_ok), f"accepted_singleton_rows={len(accepted_singleton)}")
    accepted_non_recurrent = non_recurrent[non_recurrent["passes_candidate_acceptance"] == True]
    accepted_non_recurrent_ok = accepted_non_recurrent.empty or (
        (accepted_non_recurrent["passes_statistical_evidence"] == True)
        & (accepted_non_recurrent["passes_effect_magnitude"] == True)
        & (accepted_non_recurrent["passes_expression_evidence"] == True)
        & (accepted_non_recurrent["passes_candidate_acceptance"] == True)
    ).all()
    add("accepted_non_recurrent_rows_pass_all_candidate_criteria", bool(accepted_non_recurrent_ok), f"accepted_non_recurrent_rows={len(accepted_non_recurrent)}")
    recurrent_candidate_flags_na = recurrent.empty or recurrent[
        [
            "passes_statistical_evidence",
            "passes_effect_magnitude",
            "passes_expression_evidence",
            "passes_candidate_acceptance",
        ]
    ].isna().all(axis=None)
    add("recurrent_rows_do_not_fabricate_candidate_acceptance_flags", bool(recurrent_candidate_flags_na), f"recurrent_rows={len(recurrent)}")
    rejected_candidates = candidate_rows[candidate_rows["passes_candidate_acceptance"] != True]
    selected_candidate_keys = context_key(
        selected_evidence_rows[
            selected_evidence_rows["candidate_pool_type"].isin(
                ["singleton_candidate", "non_recurrent_candidate"]
            )
        ]
    )
    rejected_selected_overlap = context_key(rejected_candidates) & selected_candidate_keys
    add("no_rejected_candidate_present_in_selected_evidence_rows", not rejected_selected_overlap, f"overlap={len(rejected_selected_overlap)}")
    missing_recurrent = context_key(recurrent) - context_key(selected_evidence_rows)
    add("every_recurrent_row_is_retained_directly", not missing_recurrent, f"missing_recurrent_rows={len(missing_recurrent)}")
    candidate_contexts = {
        (
            str(row.cancer_type),
            str(row.marker_evidence_stratum),
            str(row.candidate_pool_type),
        )
        for row in candidate_rows.itertuples(index=False)
    }
    threshold_contexts = {
        (
            str(row.cancer_type),
            str(row.marker_evidence_stratum),
            str(row.candidate_pool_type),
        )
        for row in threshold_table.itertuples(index=False)
    }
    finite_thresholds = True
    if not threshold_table.empty:
        finite_thresholds = threshold_table[
            [
                "adjusted_p_value_quantile_threshold",
                "absolute_shrunken_log2fc_quantile_threshold",
                "expression_quantile_threshold",
            ]
        ].apply(lambda col: pd.to_numeric(col, errors="coerce").map(math.isfinite)).all(axis=None)
    add("non_empty_candidate_contexts_have_finite_empirical_thresholds", bool(finite_thresholds), f"threshold_contexts={len(threshold_contexts)}")
    add("threshold_table_contexts_match_candidate_contexts_exactly", candidate_contexts == threshold_contexts, f"candidate_contexts={len(candidate_contexts)} threshold_contexts={len(threshold_contexts)}")
    expected_final_genes = set(selected_evidence_rows["gene_id"]) - set(removed_ribo_mt)
    observed_final_genes = set(features["gene_id"])
    add("final_gene_set_equals_deduplicated_selected_union", expected_final_genes == observed_final_genes, f"expected={len(expected_final_genes)} observed={len(observed_final_genes)}")
    clean_genes = read_gene_list(outdir / "pan_cancer_features_clean.txt")
    add("clean_feature_list_matches_final_feature_table_order", clean_genes == features["gene_id"].tolist(), f"clean_genes={len(clean_genes)} feature_rows={len(features)}")
    directional_ok = True
    directional_details = []
    for direction in ["UP", "DOWN", "MIXED"]:
        genes = read_gene_list(outdir / f"pan_cancer_features.{direction}.txt")
        expected = features.loc[features["direction"] == direction, "gene_id"].tolist()
        if genes != expected:
            directional_ok = False
        directional_details.append(f"{direction}={len(genes)}")
    add("directional_feature_lists_match_final_feature_table", directional_ok, ";".join(directional_details))
    selection_rank = pd.to_numeric(features["selection_rank"], errors="coerce")
    ranks_ok = (
        selection_rank.notna().all()
        and not selection_rank.duplicated().any()
        and selection_rank.astype(int).tolist() == list(range(1, len(features) + 1))
    )
    add("selection_rank_is_deterministic_and_unique", bool(ranks_ok), f"n_ranks={selection_rank.nunique()} n_features={len(features)}")
    obsolete_present = sorted(
        col
        for col in obsolete_rule_columns()
        if any(col in df.columns for df in [features, gene_evidence, selected_evidence_rows, candidate_pool_evidence])
    )
    add("obsolete_rule_selection_columns_absent", not obsolete_present, ",".join(obsolete_present) or "none")
    add("historical_target_panel_size_not_used_as_validation_criterion", True, f"final_size={features['gene_id'].nunique()}")
    add("direction_consistency_classes_are_explicit", features["direction_consistency_class"].isin(["consistently_upregulated", "consistently_downregulated", "mixed_direction"]).all(), ";".join(f"{key}={value}" for key, value in features["direction_consistency_class"].value_counts().sort_index().items()))
    add("provenance_annotations_present", features["marker_evidence_provenance"].isin(["anchor_associated", "isolate_associated", "dual_graph_derived_evidence"]).all(), ";".join(f"{key}={value}" for key, value in features["marker_evidence_provenance"].value_counts().sort_index().items()))
    add("final_size_reported", features["gene_id"].nunique() > 0, f"final_size={features['gene_id'].nunique()}")
    expected_candidate_pools = {"recurrent", "singleton_candidate", "non_recurrent_candidate"}
    add("candidate_pool_names_are_method_explicit", set(candidate_pool_evidence["candidate_pool_type"]).issubset(expected_candidate_pools), ";".join(sorted(set(candidate_pool_evidence["candidate_pool_type"]))))
    return checks


def build_summary_table(
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selected_evidence_rows: pd.DataFrame,
    candidate_pool_evidence: pd.DataFrame,
    family_summary: pd.DataFrame,
    duplicate_gene_count: int,
    removed_ribo_mt: set[str],
    *,
    # Recorded in the build summary so the configured acceptance quantiles are
    # auditable from the outputs.
    adjusted_p_value_quantile: float,
    absolute_shrunken_log2fc_quantile: float,
    expression_quantile: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(metric: str, value: object, **extra) -> None:
        row = {"method": METHOD_NAME, "metric": metric, "value": value}
        row.update(extra)
        rows.append(row)

    recurrent = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "recurrent"]
    singleton = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "singleton_candidate"]
    non_recurrent = candidate_pool_evidence[candidate_pool_evidence["candidate_pool_type"] == "non_recurrent_candidate"]
    accepted_singleton = singleton[singleton["passes_candidate_acceptance"] == True]
    accepted_non_recurrent = non_recurrent[non_recurrent["passes_candidate_acceptance"] == True]
    add("recurrence_min_count", RECURRENCE_MIN_COUNT)
    add("adjusted_p_value_quantile", adjusted_p_value_quantile)
    add("absolute_shrunken_log2fc_quantile", absolute_shrunken_log2fc_quantile)
    add("expression_quantile", expression_quantile)
    add("final_unique_gene_count", int(features["gene_id"].nunique()))
    add("recurrent_selected_rows", int((selected_evidence_rows["selection_basis"] == "recurrent").sum()))
    add("unique_recurrent_genes", int(recurrent["gene_id"].nunique()))
    add("singleton_candidate_rows", int(len(singleton)))
    add("unique_singleton_candidate_genes", int(singleton["gene_id"].nunique()))
    add("accepted_singleton_rows", int(len(accepted_singleton)))
    add("unique_accepted_singleton_genes", int(accepted_singleton["gene_id"].nunique()))
    add("rejected_singleton_rows", int(len(singleton) - len(accepted_singleton)))
    add("non_recurrent_candidate_rows", int(len(non_recurrent)))
    add("unique_non_recurrent_candidate_genes", int(non_recurrent["gene_id"].nunique()))
    add("accepted_non_recurrent_rows", int(len(accepted_non_recurrent)))
    add("unique_accepted_non_recurrent_genes", int(accepted_non_recurrent["gene_id"].nunique()))
    add("rejected_non_recurrent_rows", int(len(non_recurrent) - len(accepted_non_recurrent)))
    add("selected_gene_evidence_rows", int(len(gene_evidence)))
    add("selected_evidence_rows_before_global_deduplication", int(len(selected_evidence_rows)))
    add("duplicate_gene_rows_removed_after_cross_cancer_and_stratum_collapse", int(duplicate_gene_count))
    add("removed_ribo_mt_histone_genes", int(len(removed_ribo_mt)))
    for row in explode_count(features, "feature_class").itertuples(index=False):
        add("final_genes_by_feature_class", int(row.n_genes), feature_class=row.value)
    for row in explode_count(features, "owner_profile").itertuples(index=False):
        add("final_genes_by_cancer_owner", int(row.n_genes), owner_profile=row.value)
    add(
        "final_genes_with_anchor_associated_marker_evidence",
        int(features["has_anchor_associated_marker_evidence"].sum()),
    )
    add(
        "final_genes_with_isolate_associated_marker_evidence",
        int(features["has_isolate_associated_marker_evidence"].sum()),
    )
    for row in explode_count(features, "marker_evidence_provenance").itertuples(index=False):
        add("marker_evidence_provenance_genes", int(row.n_genes), marker_evidence_provenance=row.value)
    for row in family_summary.itertuples(index=False):
        add(
            "contributing_marker_lists_by_cancer_type_and_evidence_stratum",
            int(row.retained_marker_list_total),
            cancer_type=row.cancer_type,
            marker_evidence_stratum=row.marker_evidence_stratum,
        )
        add(
            "marker_evidence_stratum_unique_genes_before_selection",
            int(row.unique_genes_before_selection),
            cancer_type=row.cancer_type,
            marker_evidence_stratum=row.marker_evidence_stratum,
            retained_marker_list_total=row.retained_marker_list_total,
        )
    return pd.DataFrame(rows)


def write_run_manifest(
    outdir: Path,
    manifest_paths: dict[str, Path],
    manifest_tables: dict[str, pd.DataFrame],
    audit_output_prefix: str,
) -> None:
    rows: list[dict[str, object]] = []
    rows.extend(
        [
            {
                "filename": "method",
                "path": METHOD_NAME,
                "file_size_bytes": pd.NA,
                "modification_time_utc": "",
                "sha256": "",
                "description": "canonical pan-cancer aggregation method",
                "category": "provenance",
            },
            {
                "filename": "recurrence_min_count",
                "path": str(RECURRENCE_MIN_COUNT),
                "file_size_bytes": pd.NA,
                "modification_time_utc": "",
                "sha256": "",
                "description": "fixed recurrence invariant",
                "category": "provenance",
            },
            {
                "filename": "candidate_acceptance_rule",
                "path": "passes_statistical_evidence AND passes_effect_magnitude AND passes_expression_evidence",
                "file_size_bytes": pd.NA,
                "modification_time_utc": "",
                "sha256": "",
                "description": "singleton and non-recurrent candidate acceptance rule",
                "category": "provenance",
            },
        ]
    )
    for profile, path in manifest_paths.items():
        rows.append(file_manifest_row(path, "input", f"{profile} contrast-level marker manifest"))
        deseq2_dir = path.parent.parent
        for row in manifest_tables[profile].itertuples(index=False):
            marker_table_path = resolve_relative_path(deseq2_dir, getattr(row, "marker_table_path"))
            rows.append(file_manifest_row(marker_table_path, "input", f"{profile} retained contrast-level marker table"))
    for path in sorted(outdir.iterdir(), key=lambda value: value.name):
        if path.is_file():
            rows.append(file_manifest_row(path, "output", "graph-derived pan-cancer feature-selection output"))
    pd.DataFrame(rows).to_csv(outdir / f"{audit_output_prefix}_run_manifest.tsv", sep="\t", index=False)


def write_active_directory_manifest(outdir: Path, audit_output_prefix: str) -> None:
    rows = [
        file_manifest_row(path, "active output", "graph-derived pan-cancer feature-selection output")
        for path in sorted(outdir.iterdir(), key=lambda value: value.name)
        if path.is_file()
    ]
    pd.DataFrame(rows).to_csv(outdir / f"{audit_output_prefix}_active_directory_manifest.tsv", sep="\t", index=False)


def build_graph_derived_pan_cancer_feature_panel(args: argparse.Namespace) -> None:
    manifest_order, marker_manifests = parse_profile_paths(args.profile_marker_manifest, "--profile-marker-manifest")
    profile_order = manifest_order
    if not profile_order:
        sys.exit("[ERROR] At least one explicit profile marker manifest is required")
    annotation = load_gene_annotation(args.gene_annotation_tsv)
    evidence_by_profile: dict[str, pd.DataFrame] = {}
    manifest_paths: dict[str, Path] = {}
    manifest_tables: dict[str, pd.DataFrame] = {}
    for profile in profile_order:
        evidence, resolved_manifest, manifest = marker_evidence_ingestion(
            profile,
            marker_manifests.get(profile),
        )
        evidence_by_profile[profile] = evidence_stratum_assignment(evidence)
        manifest_paths[profile] = resolved_manifest
        manifest_tables[profile] = normalise_manifest(manifest, profile)
    aggregated_gene_evidence, family_summary = within_stratum_gene_evidence_aggregation(
        profile_order,
        evidence_by_profile,
        manifest_tables,
    )
    candidate_pool_evidence = candidate_pool_construction(
        recurrence_classification(aggregated_gene_evidence)
    )
    candidate_pool_evidence, threshold_table = empirical_quantile_threshold_calculation(
        candidate_pool_evidence,
        args.adjusted_p_value_quantile,
        args.absolute_shrunken_log2fc_quantile,
        args.expression_quantile,
    )
    candidate_pool_evidence = effect_direction_consistency_classification(candidate_pool_evidence)
    candidate_pool_evidence = evaluate_candidate_acceptance(candidate_pool_evidence)
    selected_evidence_rows = select_pan_cancer_evidence_rows(candidate_pool_evidence)
    gene_evidence = build_gene_evidence_table(selected_evidence_rows)
    features, duplicate_gene_count, removed_ribo_mt = collapse_gene_evidence_across_cancer_types_and_strata(
        profile_order,
        gene_evidence,
        annotation,
        args.remove_ribo_mt,
    )
    validation_rows = selected_pan_cancer_panel_export(
        Path(args.output_dir),
        features,
        gene_evidence,
        selected_evidence_rows,
        candidate_pool_evidence,
        threshold_table,
        family_summary,
        profile_order,
        manifest_paths,
        manifest_tables,
        args,
        duplicate_gene_count,
        removed_ribo_mt,
    )
    print(
        f"[Feature selection] method={METHOD_NAME} recurrence_min_count={RECURRENCE_MIN_COUNT} "
        f"final_genes={features['gene_id'].nunique()} selected_evidence_rows={len(selected_evidence_rows)} validation="
        f"{'PASS' if all(row['status'] == 'PASS' for row in validation_rows) else 'FAIL'}",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-marker-manifest", action="append", default=[], metavar="PROFILE=PATH")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method", default=METHOD_NAME)
    parser.add_argument("--adjusted-p-value-quantile", type=float, required=True)
    parser.add_argument("--absolute-shrunken-log2fc-quantile", type=float, required=True)
    parser.add_argument("--expression-quantile", type=float, required=True)
    parser.add_argument("--audit-output-prefix", default=DEFAULT_AUDIT_OUTPUT_PREFIX)
    parser.add_argument("--remove-ribo-mt", action="store_true", default=False)
    parser.add_argument("--gene-annotation-tsv", default="")
    args = parser.parse_args()
    if args.method != METHOD_NAME:
        sys.exit(f"[ERROR] --method must be {METHOD_NAME}")
    # The configured quantiles are validated for admissibility, not for equality
    # with a value embedded here: configuration owns them, so changing
    # marker_postprocessing.pan_cancer.empirical_quantile_thresholds must change
    # what this stage computes.
    for _flag, _value in (
        ("--adjusted-p-value-quantile", args.adjusted_p_value_quantile),
        ("--absolute-shrunken-log2fc-quantile", args.absolute_shrunken_log2fc_quantile),
        ("--expression-quantile", args.expression_quantile),
    ):
        if not (0.0 < _value < 1.0):
            sys.exit(f"[ERROR] {_flag} must satisfy 0 < q < 1; got {_value}")
    if not args.audit_output_prefix or not re.match(r"^[A-Za-z0-9_.-]+$", args.audit_output_prefix):
        sys.exit("[ERROR] --audit-output-prefix must be a non-empty file-name prefix")
    build_graph_derived_pan_cancer_feature_panel(args)


if __name__ == "__main__":
    main()
