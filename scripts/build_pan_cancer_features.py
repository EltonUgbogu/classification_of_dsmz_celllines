#!/usr/bin/env python3
"""
Build the ranked marker-source pan-cancer panel from retained marker outputs.

This script intentionally consumes only the already generated retained marker
lists, marker manifests, and full DESeq2 tables. It does not alter the upstream
per-contrast anchor or isolate marker-selection thresholds/caps.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import re
import shlex
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


METHOD_NAME = "ranked_marker_source_pan_cancer_panel"
OUTPUT_PREFIX = "ranked_marker_source_panel"
MARKER_SOURCE_CLASS_ORDER = ("anchor", "isolate")
RECURRENT_CLASSES = {"anchor_source_recurrent", "isolate_source_recurrent"}
EVIDENCE_CLASS_PRIORITY = (
    "multi_source_recurrent",
    "multi_cohort_recurrent",
    "anchor_source_recurrent",
    "isolate_source_recurrent",
    "anchor_singleton_ranked",
    "isolate_singleton_ranked",
    "ranked_nonrecurrent_marker",
)
OLD_LOGIC_RESCUE_PATTERN = r"isolate_extension|isolate_rescue|rescued"
RULE_SPECS = (
    ("strict_iqr_median_baseMean", "strict_iqr", "median"),
    ("strict_iqr_q1_baseMean", "strict_iqr", "q1"),
    ("relaxed_iqr_median_baseMean", "relaxed_iqr", "median"),
)


def strip_ensg_version(gene_id: object) -> str:
    if pd.isna(gene_id):
        return ""
    return str(gene_id).split(".", 1)[0]


def marker_source_class_from_contrast(contrast_id: str) -> str:
    if contrast_id.startswith("anchor_"):
        return "anchor"
    if contrast_id.startswith("isolate_"):
        return "isolate"
    return "other"


def marker_source_class_sort_value(marker_source_class: str) -> int:
    try:
        return MARKER_SOURCE_CLASS_ORDER.index(marker_source_class)
    except ValueError:
        return len(MARKER_SOURCE_CLASS_ORDER)


def evidence_class_sort_value(evidence_class: str) -> int:
    try:
        return EVIDENCE_CLASS_PRIORITY.index(evidence_class)
    except ValueError:
        return len(EVIDENCE_CLASS_PRIORITY)


def parse_profile_paths(entries: list[str], label: str) -> tuple[list[str], dict[str, Path]]:
    order: list[str] = []
    paths: dict[str, Path] = {}
    for entry in entries:
        if "=" not in entry:
            sys.exit(f"[ERROR] {label} must be PROFILE=PATH, got: {entry}")
        profile, path = entry.split("=", 1)
        profile = profile.strip()
        if not profile:
            sys.exit(f"[ERROR] Invalid profile name in {label}: {entry}")
        if profile not in paths:
            order.append(profile)
        paths[profile] = Path(path.strip())
    return order, paths


def parse_labeled_paths(entries: list[str]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for entry in entries:
        if "=" not in entry:
            sys.exit(f"[ERROR] --previous-feature-set must be LABEL=PATH, got: {entry}")
        label, path = entry.split("=", 1)
        label = label.strip()
        if not label:
            sys.exit(f"[ERROR] Empty previous feature-set label in: {entry}")
        paths[label] = Path(path.strip())
    return paths


def resolve_manifest_path(marker_dir: Path | None, manifest_path: Path | None) -> Path:
    if manifest_path is not None:
        return manifest_path
    if marker_dir is None:
        sys.exit("[ERROR] Need either --profile-marker-manifest or --profile-marker-dir")
    return marker_dir / "marker_sets_manifest.tsv"


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
    if vals.empty:
        return default
    return float(vals.min())


def safe_max(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    if vals.empty:
        return default
    return float(vals.max())


def safe_mean(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    if vals.empty:
        return default
    return float(vals.mean())


def safe_median(values: pd.Series, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    if vals.empty:
        return default
    return float(vals.median())


def safe_quantile(values: pd.Series, q: float, default: float = math.nan) -> float:
    vals = finite_numeric(values)
    if vals.empty:
        return default
    return float(vals.quantile(q))


def join_unique(values, sep: str = ";") -> str:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = str(value or "").strip()
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
            if stat
            else ""
        ),
        "sha256": sha256_file(path) if exists and path.is_file() else "",
        "description": description,
        "category": category,
    }


def load_gene_annotation(path: str) -> dict[str, dict[str, str]]:
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    id_col = next(
        (col for col in ("gene_id", "ensembl_id", "ensembl_gene_id") if col in df.columns),
        None,
    )
    if id_col is None:
        print(f"[WARN] Gene annotation missing Ensembl ID column: {path}", file=sys.stderr)
        return {}
    symbol_col = next(
        (col for col in ("symbol", "gene_symbol", "gene_name", "hgnc_symbol") if col in df.columns),
        None,
    )
    if symbol_col is None:
        print(f"[WARN] Gene annotation missing gene symbol column: {path}", file=sys.stderr)
        return {}
    name_col = "gene_name" if "gene_name" in df.columns else symbol_col
    annotations: dict[str, dict[str, str]] = {}
    for row in df.itertuples(index=False):
        clean_id = strip_ensg_version(getattr(row, id_col))
        if not clean_id:
            continue
        annotations[clean_id] = {
            "gene_symbol": str(getattr(row, symbol_col) or "").strip(),
            "gene_name": str(getattr(row, name_col) or "").strip(),
        }
    return annotations


def is_ribo_mt_or_histone(clean_gene_id: str, annotation: dict[str, dict[str, str]]) -> bool:
    symbol = annotation.get(clean_gene_id, {}).get("gene_symbol", "").upper()
    return bool(symbol and re.match(r"^(RPL|RPS|MRPL|MRPS|MT-|HIST)", symbol))


def load_gene_list(path: Path) -> list[str]:
    genes: list[str] = []
    with path.open() as handle:
        for line in handle:
            gene = strip_ensg_version(line.strip())
            if gene and not gene.startswith("#"):
                genes.append(gene)
    return genes


def load_contrast_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"[ERROR] Missing DESeq2 result table: {path}")
    df = pd.read_csv(path, sep="\t")
    required = ["gene_id", "log2FoldChange", "padj"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {path} missing required columns: {missing}")
    df = df.copy()
    df["gene_id_raw"] = df["gene_id"].astype(str)
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)
    for col in ("baseMean", "log2FoldChange", "padj", "normalised_count_in_test_sample"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        else:
            df[col] = math.nan
    return df


def load_profile_marker_evidence(
    profile: str,
    marker_dir: Path | None,
    manifest_path: Path | None,
) -> tuple[pd.DataFrame, Path, pd.DataFrame]:
    manifest_path = resolve_manifest_path(marker_dir, manifest_path)
    if not manifest_path.exists():
        sys.exit(f"[ERROR] Missing marker manifest for {profile}: {manifest_path}")
    deseq2_dir = manifest_path.parent.parent
    manifest = pd.read_csv(manifest_path, sep="\t")
    required = ["contrast", "marker_file", "table_file", "n_markers"]
    missing = [col for col in required if col not in manifest.columns]
    if missing:
        sys.exit(f"[ERROR] {manifest_path} missing required columns: {missing}")

    records: list[dict[str, object]] = []
    for contrast_order, row in enumerate(manifest.itertuples(index=False), start=1):
        contrast_id = str(getattr(row, "contrast"))
        marker_source_class = marker_source_class_from_contrast(contrast_id)
        if marker_source_class not in MARKER_SOURCE_CLASS_ORDER:
            continue
        marker_path = resolve_relative_path(deseq2_dir, getattr(row, "marker_file"))
        table_path = resolve_relative_path(deseq2_dir, getattr(row, "table_file"))
        if not marker_path.exists():
            sys.exit(f"[ERROR] Missing marker list for {profile} {contrast_id}: {marker_path}")
        table = load_contrast_table(table_path)
        table_by_gene = {
            str(table_row["gene_id"]): table_row
            for _, table_row in table.drop_duplicates("gene_id", keep="first").iterrows()
        }
        for marker_rank, gene_id in enumerate(load_gene_list(marker_path), start=1):
            table_row = table_by_gene.get(gene_id)
            if table_row is None:
                sys.exit(
                    f"[ERROR] Marker gene {gene_id} from {marker_path} not found in {table_path}"
                )
            log2fc = parse_float(table_row["log2FoldChange"])
            records.append(
                {
                    "profile": profile,
                    "contrast": contrast_id,
                    "marker_source_class": marker_source_class,
                    "contrast_order": int(contrast_order),
                    "gene_id": gene_id,
                    "gene_id_raw": str(table_row.get("gene_id_raw", gene_id)),
                    "marker_rank": int(marker_rank),
                    "padj": parse_float(table_row["padj"], default=math.inf),
                    "log2FoldChange": log2fc,
                    "abs_log2FC": abs(log2fc) if not math.isnan(log2fc) else math.nan,
                    "baseMean": parse_float(table_row.get("baseMean", math.nan)),
                    "normalised_count_in_test_sample": parse_float(
                        table_row.get("normalised_count_in_test_sample", math.nan)
                    ),
                    "source_file": str(marker_path),
                    "table_file": str(table_path),
                }
            )
    if not records:
        sys.exit(f"[ERROR] No retained marker evidence loaded for {profile}")
    return pd.DataFrame(records), manifest_path, manifest


def collapse_direction(evidence: pd.DataFrame) -> str:
    if evidence.empty:
        return "UNKNOWN"
    has_up = (pd.to_numeric(evidence["log2FoldChange"], errors="coerce") > 0).any()
    has_down = (pd.to_numeric(evidence["log2FoldChange"], errors="coerce") < 0).any()
    if has_up and has_down:
        return "MIXED"
    if has_up:
        return "UP"
    if has_down:
        return "DOWN"
    return "UNKNOWN"


def collapse_direction_from_values(values: pd.Series) -> str:
    directions = {str(v).upper() for v in values if str(v).strip()}
    if "MIXED" in directions:
        return "MIXED"
    has_up = "UP" in directions
    has_down = "DOWN" in directions
    if has_up and has_down:
        return "MIXED"
    if has_up:
        return "UP"
    if has_down:
        return "DOWN"
    return "UNKNOWN"


def aggregate_family_gene_metrics(
    profile: str,
    marker_source_class: str,
    family_evidence: pd.DataFrame,
    n_sources: int,
    recurrence_k: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if family_evidence.empty:
        return pd.DataFrame(rows)
    contrast_order = (
        family_evidence[["contrast", "contrast_order"]]
        .drop_duplicates()
        .set_index("contrast")["contrast_order"]
        .to_dict()
    )
    for gene_id, gene_rows in family_evidence.groupby("gene_id", sort=True):
        source_contrasts = sorted(
            gene_rows["contrast"].astype(str).unique(),
            key=lambda c: (contrast_order.get(c, 10**9), c),
        )
        source_files = [
            str(x)
            for x in gene_rows.sort_values(["contrast_order", "marker_rank"])["source_file"].unique()
        ]
        source_count = int(gene_rows["contrast"].nunique())
        direction = collapse_direction(gene_rows)
        if n_sources == 1:
            candidate_pool = "singleton_source"
        elif source_count >= recurrence_k:
            candidate_pool = "marker_source_recurrent"
        else:
            candidate_pool = "nonrecurrent_marker"
        ranks = []
        for contrast, contrast_rows in gene_rows.groupby("contrast", sort=False):
            ranks.append(f"{contrast}:{int(contrast_rows['marker_rank'].min())}")
        rows.append(
            {
                "cohort": profile,
                "marker_source_class": marker_source_class,
                "gene_id": gene_id,
                "gene_id_raw": join_unique(gene_rows["gene_id_raw"]),
                "source_count": source_count,
                "marker_source_class_source_count": int(n_sources),
                "candidate_pool": candidate_pool,
                "source_contrast": ";".join(source_contrasts),
                "source_file": ";".join(source_files),
                "source_marker_ranks": ";".join(ranks),
                "min_padj": safe_min(gene_rows["padj"], default=math.inf),
                "median_padj": safe_median(gene_rows["padj"], default=math.inf),
                "max_abs_log2FC": safe_max(gene_rows["abs_log2FC"]),
                "median_abs_log2FC": safe_median(gene_rows["abs_log2FC"]),
                "median_log2FC": safe_median(gene_rows["log2FoldChange"]),
                "mean_baseMean": safe_mean(gene_rows["baseMean"]),
                "median_baseMean": safe_median(gene_rows["baseMean"]),
                "max_baseMean": safe_max(gene_rows["baseMean"]),
                "best_marker_rank": int(gene_rows["marker_rank"].min()),
                "direction": direction,
                "direction_consistent": direction != "MIXED",
            }
        )
    return pd.DataFrame(rows)


def sort_rank_components(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    ranked = df.copy()
    ranked["_direction_consistency_order"] = ranked["direction_consistent"].map(bool).astype(int)
    ranked = ranked.sort_values(
        [
            "cohort",
            "marker_source_class",
            "candidate_pool",
            "source_count",
            "min_padj",
            "median_abs_log2FC",
            "median_baseMean",
            "_direction_consistency_order",
            "best_marker_rank",
            "gene_id",
        ],
        ascending=[True, True, True, False, True, False, False, False, True, True],
        kind="mergesort",
    )
    ranked["rank_within_candidate_pool"] = (
        ranked.groupby(["cohort", "marker_source_class", "candidate_pool"], sort=False).cumcount() + 1
    )
    ranked["rank_within_cohort_marker_source_class"] = (
        ranked.groupby(["cohort", "marker_source_class"], sort=False).cumcount() + 1
    )
    ranked["best_padj_rank_within_pool"] = (
        ranked.groupby(["cohort", "marker_source_class", "candidate_pool"])["min_padj"]
        .rank(method="min", ascending=True)
        .astype("Int64")
    )
    ranked["ranking_component_key"] = (
        ranked["cohort"].astype(str)
        + "|"
        + ranked["marker_source_class"].astype(str)
        + "|"
        + ranked["candidate_pool"].astype(str)
        + "|"
        + ranked["gene_id"].astype(str)
    )
    return ranked.drop(columns=["_direction_consistency_order"])


def annotate_empirical_rule_flags(
    metrics: pd.DataFrame,
    padj_quantile: float,
    abs_log2fc_quantile: float,
    relaxed_minimum_criteria: int,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    if metrics.empty:
        return metrics, []
    annotated = metrics.copy()
    for col in (
        "padj_q1_threshold",
        "abs_log2fc_q3_threshold",
        "baseMean_median_threshold",
        "baseMean_q1_threshold",
    ):
        annotated[col] = math.nan
    for col in (
        "padj_q1_pass",
        "abs_log2fc_q3_pass",
        "baseMean_median_pass",
        "baseMean_q1_pass",
        "strict_iqr_median_pass",
        "strict_iqr_q1_pass",
        "relaxed_iqr_median_pass",
    ):
        annotated[col] = False
    annotated["relaxed_iqr_median_criteria_met"] = 0

    threshold_rows: list[dict[str, object]] = []
    empirical_pools = annotated[
        annotated["candidate_pool"].isin(["singleton_source", "nonrecurrent_marker"])
    ]
    for (cohort, marker_source_class, candidate_pool), pool in empirical_pools.groupby(
        ["cohort", "marker_source_class", "candidate_pool"], sort=True
    ):
        idx = pool.index
        padj_threshold = safe_quantile(pool["min_padj"], padj_quantile)
        abs_threshold = safe_quantile(pool["median_abs_log2FC"], abs_log2fc_quantile)
        base_median = safe_median(pool["median_baseMean"])
        base_q1 = safe_quantile(pool["median_baseMean"], 0.25)

        padj_pass = pool["min_padj"] <= padj_threshold
        effect_pass = pool["median_abs_log2FC"] >= abs_threshold
        base_median_pass = pool["median_baseMean"] >= base_median
        base_q1_pass = pool["median_baseMean"] >= base_q1
        relaxed_criteria = (
            padj_pass.astype(int) + effect_pass.astype(int) + base_median_pass.astype(int)
        )
        direction_ok_for_relaxed = pool["direction_consistent"] | (relaxed_criteria >= 3)

        annotated.loc[idx, "padj_q1_threshold"] = padj_threshold
        annotated.loc[idx, "abs_log2fc_q3_threshold"] = abs_threshold
        annotated.loc[idx, "baseMean_median_threshold"] = base_median
        annotated.loc[idx, "baseMean_q1_threshold"] = base_q1
        annotated.loc[idx, "padj_q1_pass"] = padj_pass.values
        annotated.loc[idx, "abs_log2fc_q3_pass"] = effect_pass.values
        annotated.loc[idx, "baseMean_median_pass"] = base_median_pass.values
        annotated.loc[idx, "baseMean_q1_pass"] = base_q1_pass.values
        annotated.loc[idx, "strict_iqr_median_pass"] = (
            padj_pass & effect_pass & base_median_pass
        ).values
        annotated.loc[idx, "strict_iqr_q1_pass"] = (
            padj_pass & effect_pass & base_q1_pass
        ).values
        annotated.loc[idx, "relaxed_iqr_median_criteria_met"] = relaxed_criteria.values
        annotated.loc[idx, "relaxed_iqr_median_pass"] = (
            (relaxed_criteria >= relaxed_minimum_criteria) & direction_ok_for_relaxed
        ).values

        for rule_key, rule_name, base_mode in RULE_SPECS:
            if rule_key == "strict_iqr_median_baseMean":
                pass_col = "strict_iqr_median_pass"
                base_threshold = base_median
                n_pass = int((padj_pass & effect_pass & base_median_pass).sum())
            elif rule_key == "strict_iqr_q1_baseMean":
                pass_col = "strict_iqr_q1_pass"
                base_threshold = base_q1
                n_pass = int((padj_pass & effect_pass & base_q1_pass).sum())
            else:
                pass_col = "relaxed_iqr_median_pass"
                base_threshold = base_median
                n_pass = int(((relaxed_criteria >= relaxed_minimum_criteria) & direction_ok_for_relaxed).sum())
            threshold_rows.append(
                {
                    "method": METHOD_NAME,
                    "rule_key": rule_key,
                    "rule_name": rule_name,
                    "cohort": cohort,
                    "marker_source_class": marker_source_class,
                    "candidate_pool": candidate_pool,
                    "n_candidate_genes": int(pool["gene_id"].nunique()),
                    "padj_quantile": padj_quantile,
                    "padj_q1_threshold": padj_threshold,
                    "abs_log2fc_quantile": abs_log2fc_quantile,
                    "abs_log2fc_q3_threshold": abs_threshold,
                    "baseMean_threshold_mode": base_mode,
                    "baseMean_threshold": base_threshold,
                    "relaxed_minimum_criteria": relaxed_minimum_criteria if rule_name == "relaxed_iqr" else pd.NA,
                    "n_padj_q1_pass": int(padj_pass.sum()),
                    "n_abs_log2fc_q3_pass": int(effect_pass.sum()),
                    "n_baseMean_median_pass": int(base_median_pass.sum()),
                    "n_baseMean_q1_pass": int(base_q1_pass.sum()),
                    "n_rule_pass": n_pass,
                    "pass_column": pass_col,
                }
            )
    return annotated, threshold_rows


def build_ranking_components(
    profile_order: list[str],
    evidence_by_profile: dict[str, pd.DataFrame],
    recurrence_k: int,
    padj_quantile: float,
    abs_log2fc_quantile: float,
    relaxed_minimum_criteria: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_metrics: list[pd.DataFrame] = []
    family_summary_rows: list[dict[str, object]] = []
    for profile in profile_order:
        evidence = evidence_by_profile[profile]
        for marker_source_class in MARKER_SOURCE_CLASS_ORDER:
            family_evidence = evidence[evidence["marker_source_class"] == marker_source_class].copy()
            source_contrasts = list(
                family_evidence.sort_values("contrast_order")["contrast"].drop_duplicates()
            )
            n_sources = len(source_contrasts)
            if n_sources == 0:
                family_summary_rows.append(
                    {
                        "method": METHOD_NAME,
                        "cohort": profile,
                        "marker_source_class": marker_source_class,
                        "source_contrasts": "",
                        "n_source_lists_in_marker_source_class": 0,
                        "recurrence_testable": False,
                        "recurrence_k": recurrence_k,
                        "unique_genes_before_selection": 0,
                        "recurrent_candidate_genes": 0,
                        "singleton_candidate_genes": 0,
                        "nonrecurrent_candidate_genes": 0,
                    }
                )
                continue
            metrics = aggregate_family_gene_metrics(
                profile, marker_source_class, family_evidence, n_sources, recurrence_k
            )
            all_metrics.append(metrics)
            family_summary_rows.append(
                {
                    "method": METHOD_NAME,
                    "cohort": profile,
                    "marker_source_class": marker_source_class,
                    "source_contrasts": ";".join(map(str, source_contrasts)),
                    "n_source_lists_in_marker_source_class": n_sources,
                    "recurrence_testable": n_sources >= 2,
                    "recurrence_k": recurrence_k,
                    "unique_genes_before_selection": int(metrics["gene_id"].nunique()),
                    "recurrent_candidate_genes": int(
                        (metrics["candidate_pool"] == "marker_source_recurrent").sum()
                    ),
                    "singleton_candidate_genes": int(
                        (metrics["candidate_pool"] == "singleton_source").sum()
                    ),
                    "nonrecurrent_candidate_genes": int(
                        (metrics["candidate_pool"] == "nonrecurrent_marker").sum()
                    ),
                }
            )

    if not all_metrics:
        sys.exit("[ERROR] No marker-source metrics were produced")
    ranking_components = pd.concat(all_metrics, ignore_index=True)
    ranking_components, threshold_rows = annotate_empirical_rule_flags(
        ranking_components,
        padj_quantile=padj_quantile,
        abs_log2fc_quantile=abs_log2fc_quantile,
        relaxed_minimum_criteria=relaxed_minimum_criteria,
    )
    ranking_components = sort_rank_components(ranking_components).reset_index(drop=True)
    family_summary = pd.DataFrame(family_summary_rows)
    threshold_df = pd.DataFrame(threshold_rows)
    return ranking_components, family_summary, threshold_df


def rule_pass_column(rule_key: str) -> str:
    if rule_key == "strict_iqr_median_baseMean":
        return "strict_iqr_median_pass"
    if rule_key == "strict_iqr_q1_baseMean":
        return "strict_iqr_q1_pass"
    if rule_key == "relaxed_iqr_median_baseMean":
        return "relaxed_iqr_median_pass"
    sys.exit(f"[ERROR] Unknown rule key: {rule_key}")


def evidence_class_for_candidate(row: pd.Series) -> str:
    if row["candidate_pool"] == "marker_source_recurrent":
        return f"{row['marker_source_class']}_source_recurrent"
    if row["candidate_pool"] == "singleton_source":
        return f"{row['marker_source_class']}_singleton_ranked"
    return "ranked_nonrecurrent_marker"


def candidate_tier(row: pd.Series) -> int:
    if row["candidate_pool"] == "marker_source_recurrent":
        return 1
    if row["candidate_pool"] == "singleton_source":
        return 2
    return 3


def estimate_unique_feature_count(
    selections: pd.DataFrame,
    annotation: dict[str, dict[str, str]],
    remove_ribo_mt: bool,
) -> int:
    genes = list(dict.fromkeys(selections["gene_id"].astype(str).tolist()))
    if remove_ribo_mt:
        genes = [gene for gene in genes if not is_ribo_mt_or_histone(gene, annotation)]
    return len(set(genes))


def build_selection_rows_for_rule(
    ranking_components: pd.DataFrame,
    rule_key: str,
    desired_min_size: int,
    force_exact_size: bool,
    annotation: dict[str, dict[str, str]],
    remove_ribo_mt: bool,
) -> tuple[pd.DataFrame, bool]:
    if force_exact_size:
        sys.exit("[ERROR] force_exact_size=true is not supported for this evidence-ranked method")
    pass_col = rule_pass_column(rule_key)
    base = ranking_components[
        (ranking_components["candidate_pool"] == "marker_source_recurrent")
        | (
            (ranking_components["candidate_pool"] == "singleton_source")
            & ranking_components[pass_col].map(bool)
        )
    ].copy()
    topup_activated = estimate_unique_feature_count(base, annotation, remove_ribo_mt) < desired_min_size
    if topup_activated:
        topup = ranking_components[
            (ranking_components["candidate_pool"] == "nonrecurrent_marker")
            & ranking_components[pass_col].map(bool)
        ].copy()
        selected = pd.concat([base, topup], ignore_index=True)
    else:
        selected = base
    if selected.empty:
        sys.exit(f"[ERROR] Rule {rule_key} selected no candidate rows")
    selected = selected.copy()
    selected["selected_rule_key"] = rule_key
    selected["selected_rule_pass_column"] = pass_col
    selected["selected_rule_pass"] = True
    selected["evidence_class"] = selected.apply(evidence_class_for_candidate, axis=1)
    selected["evidence_tier"] = selected.apply(candidate_tier, axis=1)
    selected["selection_class_priority"] = selected["evidence_class"].map(evidence_class_sort_value)
    selected["selection_marker_source_class_priority"] = selected["marker_source_class"].map(marker_source_class_sort_value)
    selected = selected.sort_values(
        [
            "evidence_tier",
            "selection_class_priority",
            "source_count",
            "min_padj",
            "median_abs_log2FC",
            "median_baseMean",
            "direction_consistent",
            "rank_within_candidate_pool",
            "gene_id",
        ],
        ascending=[True, True, False, True, False, False, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    selected["selection_row_rank"] = range(1, len(selected) + 1)
    return selected, bool(topup_activated)


def primary_evidence_class(classes: set[str]) -> str:
    if not classes:
        return ""
    return sorted(classes, key=evidence_class_sort_value)[0]


def build_gene_evidence_table(
    selections: pd.DataFrame,
    evidence_by_profile: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    selected_pairs = selections[["cohort", "gene_id"]].drop_duplicates()
    class_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    selected_family_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    selection_rows_by_pair: dict[tuple[str, str], pd.DataFrame] = {}
    for (cohort, gene_id), group in selections.groupby(["cohort", "gene_id"], sort=False):
        key = (cohort, gene_id)
        selection_rows_by_pair[key] = group.copy()
        for row in group.itertuples(index=False):
            class_by_pair[key].add(row.evidence_class)
            selected_family_by_pair[key].add(row.marker_source_class)

    rows: list[dict[str, object]] = []
    for pair in selected_pairs.itertuples(index=False):
        cohort = pair.cohort
        gene_id = pair.gene_id
        key = (cohort, gene_id)
        profile_evidence = evidence_by_profile[cohort]
        retained_evidence = profile_evidence[profile_evidence["gene_id"] == gene_id].copy()
        selected = selection_rows_by_pair[key]
        classes = set(class_by_pair[key])
        recurrent_families = set(
            selected.loc[selected["evidence_class"].isin(RECURRENT_CLASSES), "marker_source_class"]
        )
        if len(recurrent_families) > 1:
            classes.add("multi_source_recurrent")
        ordered_classes = sorted(classes, key=evidence_class_sort_value)
        selected_families = sorted(
            selected_family_by_pair[key], key=marker_source_class_sort_value
        )
        retained_families = sorted(
            set(retained_evidence["marker_source_class"].astype(str)), key=marker_source_class_sort_value
        )
        recurrence_details = [
            f"{row.marker_source_class}:{int(row.source_count)}/{int(row.marker_source_class_source_count)}"
            for row in selected.itertuples(index=False)
        ]
        rows.append(
            {
                "method": METHOD_NAME,
                "gene_id": gene_id,
                "gene_id_raw": join_unique(retained_evidence["gene_id_raw"]),
                "cohort": cohort,
                "marker_source_class": ";".join(selected_families),
                "retained_marker_source_classes": ";".join(retained_families),
                "evidence_classes": ";".join(ordered_classes),
                "primary_evidence_class": primary_evidence_class(classes),
                "evidence_tiers": join_unique(selected["evidence_tier"]),
                "marker_source_recurrence_detail": ";".join(recurrence_details),
                "ranking_component_keys": join_unique(selected["ranking_component_key"]),
                "selected_rule_key": join_unique(selected["selected_rule_key"]),
                "rank_within_candidate_pool": join_unique(selected["rank_within_candidate_pool"]),
                "rank_within_cohort_marker_source_class": join_unique(selected["rank_within_cohort_marker_source_class"]),
                "best_padj_rank_within_pool": join_unique(selected["best_padj_rank_within_pool"]),
                "padj_q1_threshold": safe_min(selected["padj_q1_threshold"]),
                "abs_log2fc_q3_threshold": safe_min(selected["abs_log2fc_q3_threshold"]),
                "baseMean_median_threshold": safe_min(selected["baseMean_median_threshold"]),
                "baseMean_q1_threshold": safe_min(selected["baseMean_q1_threshold"]),
                "padj_q1_pass": bool(selected["padj_q1_pass"].any()),
                "abs_log2fc_q3_pass": bool(selected["abs_log2fc_q3_pass"].any()),
                "baseMean_median_pass": bool(selected["baseMean_median_pass"].any()),
                "baseMean_q1_pass": bool(selected["baseMean_q1_pass"].any()),
                "strict_iqr_median_pass": bool(selected["strict_iqr_median_pass"].any()),
                "strict_iqr_q1_pass": bool(selected["strict_iqr_q1_pass"].any()),
                "relaxed_iqr_median_pass": bool(selected["relaxed_iqr_median_pass"].any()),
                "relaxed_iqr_median_criteria_met": int(
                    pd.to_numeric(selected["relaxed_iqr_median_criteria_met"], errors="coerce")
                    .fillna(0)
                    .max()
                ),
                "source_contrast": join_unique(
                    retained_evidence.sort_values(["contrast_order", "marker_rank"])["contrast"]
                ),
                "source_file": join_unique(
                    retained_evidence.sort_values(["contrast_order", "marker_rank"])["source_file"]
                ),
                "source_marker_ranks": join_unique(
                    [
                        f"{row.contrast}:{int(row.marker_rank)}"
                        for row in retained_evidence.sort_values(
                            ["contrast_order", "marker_rank"]
                        ).itertuples(index=False)
                    ]
                ),
                "source_count": int(selected["source_count"].max()),
                "marker_source_class_source_count": int(selected["marker_source_class_source_count"].max()),
                "n_retained_contrasts": int(retained_evidence["contrast"].nunique()),
                "n_anchor_contrasts": int(
                    retained_evidence.query("marker_source_class == 'anchor'")["contrast"].nunique()
                ),
                "n_isolate_contrasts": int(
                    retained_evidence.query("marker_source_class == 'isolate'")["contrast"].nunique()
                ),
                "min_padj": safe_min(retained_evidence["padj"], default=math.inf),
                "median_padj": safe_median(retained_evidence["padj"], default=math.inf),
                "max_abs_log2FC": safe_max(retained_evidence["abs_log2FC"]),
                "median_abs_log2FC": safe_median(retained_evidence["abs_log2FC"]),
                "median_log2FC": safe_median(retained_evidence["log2FoldChange"]),
                "mean_baseMean": safe_mean(retained_evidence["baseMean"]),
                "median_baseMean": safe_median(retained_evidence["baseMean"]),
                "max_baseMean": safe_max(retained_evidence["baseMean"]),
                "best_marker_rank": int(retained_evidence["marker_rank"].min()),
                "direction": collapse_direction(retained_evidence),
                "direction_consistent": collapse_direction(retained_evidence) != "MIXED",
            }
        )
    gene_evidence = pd.DataFrame(rows)
    if gene_evidence.empty:
        sys.exit("[ERROR] No selected gene evidence rows were produced")
    return gene_evidence


def build_final_features(
    profile_order: list[str],
    gene_evidence: pd.DataFrame,
    annotation: dict[str, dict[str, str]],
    remove_ribo_mt: bool,
) -> tuple[pd.DataFrame, int, set[str]]:
    profile_index = {profile: i for i, profile in enumerate(profile_order)}
    candidate = gene_evidence.copy()
    candidate["_profile_order"] = candidate["cohort"].map(lambda p: profile_index.get(p, 10**6))
    candidate["_class_order"] = candidate["primary_evidence_class"].map(evidence_class_sort_value)
    candidate["_direction_consistency_order"] = candidate["direction_consistent"].map(bool).astype(int)
    candidate = candidate.sort_values(
        [
            "_class_order",
            "source_count",
            "min_padj",
            "median_abs_log2FC",
            "median_baseMean",
            "_direction_consistency_order",
            "_profile_order",
            "best_marker_rank",
            "gene_id",
        ],
        ascending=[True, False, True, False, False, False, True, True, True],
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
        removed_ribo_mt = {gene for gene in ordered_genes if is_ribo_mt_or_histone(gene, annotation)}
        ordered_genes = [gene for gene in ordered_genes if gene not in removed_ribo_mt]

    rows: list[dict[str, object]] = []
    evidence_by_gene = {
        gene_id: df.copy()
        for gene_id, df in gene_evidence[gene_evidence["gene_id"].isin(ordered_genes)].groupby(
            "gene_id", sort=False
        )
    }
    for selection_rank, gene_id in enumerate(ordered_genes, start=1):
        evidence = evidence_by_gene[gene_id].copy()
        evidence["_profile_order"] = evidence["cohort"].map(lambda p: profile_index.get(p, 10**6))
        evidence["_class_order"] = evidence["primary_evidence_class"].map(evidence_class_sort_value)
        evidence = evidence.sort_values(
            [
                "_class_order",
                "source_count",
                "min_padj",
                "median_abs_log2FC",
                "median_baseMean",
                "direction_consistent",
                "_profile_order",
                "best_marker_rank",
            ],
            ascending=[True, False, True, False, False, False, True, True],
            kind="mergesort",
        )
        owner_profile = str(evidence.iloc[0]["cohort"])
        all_classes: set[str] = set()
        recurrent_cohorts: set[str] = set()
        recurrent_families: set[str] = set()
        for row in evidence.itertuples(index=False):
            classes = {x for x in str(row.evidence_classes).split(";") if x}
            all_classes.update(classes)
            if classes & RECURRENT_CLASSES:
                recurrent_cohorts.add(str(row.cohort))
                for family in str(row.marker_source_class).split(";"):
                    if family:
                        recurrent_families.add(family)
        if len(recurrent_families) > 1:
            all_classes.add("multi_source_recurrent")
        if len(recurrent_cohorts) > 1:
            all_classes.add("multi_cohort_recurrent")
        ordered_classes = sorted(all_classes, key=evidence_class_sort_value)
        annotations = annotation.get(gene_id, {})
        rows.append(
            {
                "method": METHOD_NAME,
                "gene_id": gene_id,
                "clean_gene_id": gene_id,
                "gene_id_raw": join_unique(evidence["gene_id_raw"]),
                "gene_symbol": annotations.get("gene_symbol", ""),
                "gene_name": annotations.get("gene_name", ""),
                "direction": collapse_direction_from_values(evidence["direction"]),
                "direction_consistent": bool(evidence["direction_consistent"].all()),
                "feature_class": primary_evidence_class(all_classes),
                "evidence_class": primary_evidence_class(all_classes),
                "evidence_classes": ";".join(ordered_classes),
                "owner_profile": owner_profile,
                "cohort": join_unique(evidence["cohort"]),
                "profiles_present": join_unique(evidence["cohort"], sep=","),
                "n_profiles": int(evidence["cohort"].nunique()),
                "marker_source_class": join_unique(evidence["marker_source_class"]),
                "retained_marker_source_classes": join_unique(evidence["retained_marker_source_classes"]),
                "source_contrast": join_unique(evidence["source_contrast"]),
                "source_file": join_unique(evidence["source_file"]),
                "source_marker_ranks": join_unique(evidence["source_marker_ranks"]),
                "marker_source_recurrence_detail": join_unique(evidence["marker_source_recurrence_detail"]),
                "ranking_component_keys": join_unique(evidence["ranking_component_keys"]),
                "selected_rule_key": join_unique(evidence["selected_rule_key"]),
                "rank_within_candidate_pool": join_unique(evidence["rank_within_candidate_pool"]),
                "best_padj_rank_within_pool": join_unique(evidence["best_padj_rank_within_pool"]),
                "source_count": int(evidence["source_count"].max()),
                "marker_source_class_source_count": int(evidence["marker_source_class_source_count"].max()),
                "n_retained_contrasts": int(evidence["n_retained_contrasts"].sum()),
                "n_anchor_contrasts": int(evidence["n_anchor_contrasts"].sum()),
                "n_isolate_contrasts": int(evidence["n_isolate_contrasts"].sum()),
                "min_padj": safe_min(evidence["min_padj"], default=math.inf),
                "median_padj": safe_median(evidence["median_padj"], default=math.inf),
                "max_abs_log2FC": safe_max(evidence["max_abs_log2FC"]),
                "median_abs_log2FC": safe_median(evidence["median_abs_log2FC"]),
                "median_log2FC": safe_median(evidence["median_log2FC"]),
                "mean_baseMean": safe_mean(evidence["mean_baseMean"]),
                "median_baseMean": safe_median(evidence["median_baseMean"]),
                "max_baseMean": safe_max(evidence["max_baseMean"]),
                "best_marker_rank": int(evidence["best_marker_rank"].min()),
                "source_rank": pd.NA,
                "selection_rank": selection_rank,
            }
        )
    features = pd.DataFrame(rows)
    duplicate_candidate_count = int(candidate["gene_id"].nunique() - features["gene_id"].nunique())
    return features, duplicate_candidate_count, removed_ribo_mt


def build_panel_for_rule(
    profile_order: list[str],
    ranking_components: pd.DataFrame,
    evidence_by_profile: dict[str, pd.DataFrame],
    annotation: dict[str, dict[str, str]],
    remove_ribo_mt: bool,
    desired_min_size: int,
    force_exact_size: bool,
    rule_key: str,
) -> dict[str, object]:
    selections, topup_activated = build_selection_rows_for_rule(
        ranking_components,
        rule_key,
        desired_min_size,
        force_exact_size,
        annotation,
        remove_ribo_mt,
    )
    gene_evidence = build_gene_evidence_table(selections, evidence_by_profile)
    features, duplicate_candidate_count, removed_ribo_mt = build_final_features(
        profile_order, gene_evidence, annotation, remove_ribo_mt
    )
    class_counts = explode_count(features, "evidence_classes")
    return {
        "rule_key": rule_key,
        "selections": selections,
        "gene_evidence": gene_evidence,
        "features": features,
        "duplicate_candidate_count": duplicate_candidate_count,
        "removed_ribo_mt": removed_ribo_mt,
        "topup_activated": topup_activated,
        "final_size": int(features["gene_id"].nunique()),
        "tier1_rows": int((selections["evidence_tier"] == 1).sum()),
        "tier2_rows": int((selections["evidence_tier"] == 2).sum()),
        "tier3_rows": int((selections["evidence_tier"] == 3).sum()),
        "evidence_class_counts": ";".join(
            f"{row.value}={row.n_genes}" for row in class_counts.itertuples(index=False)
        ),
    }


def choose_primary_build(builds: dict[str, dict[str, object]], args: argparse.Namespace) -> str:
    desired = int(args.desired_min_size)
    if args.primary_rule == "strict_iqr":
        return (
            "strict_iqr_q1_baseMean"
            if args.strict_basemean_threshold == "q1"
            else "strict_iqr_median_baseMean"
        )
    if args.primary_rule == "relaxed_iqr":
        return "relaxed_iqr_median_baseMean"
    if int(builds["strict_iqr_median_baseMean"]["final_size"]) >= desired:
        return "strict_iqr_median_baseMean"
    if int(builds["strict_iqr_q1_baseMean"]["final_size"]) >= desired:
        return "strict_iqr_q1_baseMean"
    if int(builds["relaxed_iqr_median_baseMean"]["final_size"]) >= desired:
        return "relaxed_iqr_median_baseMean"
    return "strict_iqr_median_baseMean"


def load_previous_features(path: Path) -> pd.DataFrame:
    if not path or not path.exists():
        return pd.DataFrame(columns=["gene_id"])
    df = pd.read_csv(path, sep="\t")
    id_col = "gene_id" if "gene_id" in df.columns else "clean_gene_id" if "clean_gene_id" in df.columns else None
    if id_col is None:
        return pd.DataFrame(columns=["gene_id"])
    df = df.copy()
    df["gene_id"] = df[id_col].map(strip_ensg_version)
    return df.drop_duplicates("gene_id", keep="first")


def build_previous_comparison(
    label: str,
    previous_path: Path,
    features: pd.DataFrame,
    outdir: Path,
) -> dict[str, int]:
    previous = load_previous_features(previous_path)
    new_by_gene = features.set_index("gene_id", drop=False)
    if previous.empty:
        removed = pd.DataFrame(columns=["gene_id"])
        added = features.copy()
        overlap = pd.DataFrame(columns=["gene_id"])
        counts = {
            f"{label}_previous_final_size": 0,
            f"{label}_retained_from_previous": 0,
            f"{label}_removed_from_previous": 0,
            f"{label}_newly_added": int(features["gene_id"].nunique()),
            f"{label}_class_changes": 0,
        }
    else:
        prev_by_gene = previous.set_index("gene_id", drop=False)
        previous_genes = set(prev_by_gene.index)
        new_genes = set(new_by_gene.index)
        removed_genes = sorted(previous_genes - new_genes)
        added_genes = sorted(new_genes - previous_genes)
        overlap_genes = sorted(previous_genes & new_genes)
        removed = prev_by_gene.loc[removed_genes].reset_index(drop=True) if removed_genes else pd.DataFrame(columns=previous.columns)
        added = new_by_gene.loc[added_genes].reset_index(drop=True) if added_genes else pd.DataFrame(columns=features.columns)
        overlap_rows: list[dict[str, object]] = []
        for gene_id in overlap_genes:
            prev_row = prev_by_gene.loc[gene_id]
            new_row = new_by_gene.loc[gene_id]
            prev_class = str(prev_row.get("feature_class", prev_row.get("evidence_class", "")))
            new_class = str(new_row.get("feature_class", new_row.get("evidence_class", "")))
            overlap_rows.append(
                {
                    "gene_id": gene_id,
                    "previous_feature_class": prev_class,
                    "new_feature_class": new_class,
                    "new_evidence_classes": new_row.get("evidence_classes", ""),
                    "class_changed": prev_class != new_class,
                    "previous_owner_profile": prev_row.get("owner_profile", ""),
                    "new_owner_profile": new_row.get("owner_profile", ""),
                    "owner_profile_changed": prev_row.get("owner_profile", "") != new_row.get("owner_profile", ""),
                    "previous_direction": prev_row.get("direction", ""),
                    "new_direction": new_row.get("direction", ""),
                    "direction_changed": prev_row.get("direction", "") != new_row.get("direction", ""),
                    "previous_source_contrast": prev_row.get("source_contrast", ""),
                    "new_source_contrast": new_row.get("source_contrast", ""),
                }
            )
        overlap = pd.DataFrame(overlap_rows)
        counts = {
            f"{label}_previous_final_size": int(len(previous_genes)),
            f"{label}_retained_from_previous": int(len(overlap_genes)),
            f"{label}_removed_from_previous": int(len(removed_genes)),
            f"{label}_newly_added": int(len(added_genes)),
            f"{label}_class_changes": int(overlap["class_changed"].sum()) if not overlap.empty else 0,
        }
    removed.to_csv(outdir / f"{OUTPUT_PREFIX}_removed_vs_{label}.tsv", sep="\t", index=False)
    added.to_csv(outdir / f"{OUTPUT_PREFIX}_added_vs_{label}.tsv", sep="\t", index=False)
    overlap.to_csv(outdir / f"{OUTPUT_PREFIX}_overlap_vs_{label}.tsv", sep="\t", index=False)
    return counts


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


def build_summary(
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selections: pd.DataFrame,
    family_summary: pd.DataFrame,
    comparison_counts: dict[str, int],
    duplicate_candidate_count: int,
    removed_ribo_mt: set[str],
    selected_rule_key: str,
    desired_min_size: int,
    topup_activated: bool,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(metric: str, value: object, **extra) -> None:
        row = {
            "method": METHOD_NAME,
            "selected_rule_key": selected_rule_key,
            "metric": metric,
            "value": value,
            "owner_profile": extra.pop("owner_profile", ""),
            "marker_source_class": extra.pop("marker_source_class", ""),
            "evidence_class": extra.pop("evidence_class", ""),
            "direction": extra.pop("direction", ""),
        }
        row.update(extra)
        rows.append(row)

    add("final_clean_genes", int(features["gene_id"].nunique()))
    add("desired_min_size", int(desired_min_size))
    add("final_size_meets_desired_minimum", bool(features["gene_id"].nunique() >= desired_min_size))
    add("topup_tier_activated", bool(topup_activated))
    add("selected_gene_evidence_rows", int(len(gene_evidence)))
    add("selected_marker_source_class_rows", int(len(selections)))
    add("duplicate_candidate_rows_removed_after_global_clean_id_deduplication", int(duplicate_candidate_count))
    add("removed_ribo_mt_histone_genes", int(len(removed_ribo_mt)))
    for key, value in comparison_counts.items():
        add(key, value)
    for direction, count in features["direction"].value_counts().sort_index().items():
        add("direction_genes", int(count), direction=str(direction))
    for cohort, count in features["owner_profile"].value_counts().sort_index().items():
        add("final_owner_profile_genes", int(count), owner_profile=str(cohort))
    for cohort, group in gene_evidence.groupby("cohort", sort=True):
        add("cohort_selected_genes_before_global_deduplication", int(group["gene_id"].nunique()), owner_profile=str(cohort))
    for row in explode_count(features, "evidence_classes").itertuples(index=False):
        add("evidence_class_genes", int(row.n_genes), evidence_class=row.value)
    for row in explode_count(features, "marker_source_class").itertuples(index=False):
        add("marker_source_class_genes", int(row.n_genes), marker_source_class=row.value)
    for row in family_summary.itertuples(index=False):
        add(
            "marker_source_class_candidate_genes",
            int(row.unique_genes_before_selection),
            owner_profile=row.cohort,
            marker_source_class=row.marker_source_class,
            recurrence_testable=row.recurrence_testable,
            n_source_lists_in_marker_source_class=row.n_source_lists_in_marker_source_class,
            recurrent_candidate_genes=row.recurrent_candidate_genes,
            singleton_candidate_genes=row.singleton_candidate_genes,
            nonrecurrent_candidate_genes=row.nonrecurrent_candidate_genes,
        )
    return pd.DataFrame(rows)


def build_cohort_summary(
    profile_order: list[str],
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selections: pd.DataFrame,
    family_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for profile in profile_order:
        evidence = gene_evidence[gene_evidence["cohort"] == profile]
        final_owner = features[features["owner_profile"] == profile]
        fam = family_summary[family_summary["cohort"] == profile]
        selected = selections[selections["cohort"] == profile]
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": profile,
                "n_anchor_sources": int(fam.loc[fam["marker_source_class"] == "anchor", "n_source_lists_in_marker_source_class"].sum()),
                "n_isolate_sources": int(fam.loc[fam["marker_source_class"] == "isolate", "n_source_lists_in_marker_source_class"].sum()),
                "anchor_selected_genes_before_global_dedup": int(selected.loc[selected["marker_source_class"] == "anchor", "gene_id"].nunique()),
                "isolate_selected_genes_before_global_dedup": int(selected.loc[selected["marker_source_class"] == "isolate", "gene_id"].nunique()),
                "tier1_recurrent_rows": int((selected["evidence_tier"] == 1).sum()),
                "tier2_singleton_ranked_rows": int((selected["evidence_tier"] == 2).sum()),
                "tier3_ranked_nonrecurrent_rows": int((selected["evidence_tier"] == 3).sum()),
                "selected_genes_before_global_dedup": int(evidence["gene_id"].nunique()),
                "final_owner_profile_genes_after_global_dedup": int(final_owner["gene_id"].nunique()),
                "final_genes_with_cohort_evidence": int(evidence["gene_id"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def build_marker_source_class_summary_output(
    features: pd.DataFrame,
    selections: pd.DataFrame,
    family_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in family_summary.itertuples(index=False):
        selected = selections[
            (selections["cohort"] == row.cohort)
            & (selections["marker_source_class"] == row.marker_source_class)
        ]
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": row.cohort,
                "marker_source_class": row.marker_source_class,
                "source_contrasts": row.source_contrasts,
                "n_source_lists_in_marker_source_class": int(row.n_source_lists_in_marker_source_class),
                "recurrence_testable": bool(row.recurrence_testable),
                "recurrence_k": int(row.recurrence_k),
                "unique_genes_before_selection": int(row.unique_genes_before_selection),
                "recurrent_candidate_genes": int(row.recurrent_candidate_genes),
                "singleton_candidate_genes": int(row.singleton_candidate_genes),
                "nonrecurrent_candidate_genes": int(row.nonrecurrent_candidate_genes),
                "selected_genes_before_global_dedup": int(selected["gene_id"].nunique()),
                "tier1_recurrent_rows": int((selected["evidence_tier"] == 1).sum()) if not selected.empty else 0,
                "tier2_singleton_ranked_rows": int((selected["evidence_tier"] == 2).sum()) if not selected.empty else 0,
                "tier3_ranked_nonrecurrent_rows": int((selected["evidence_tier"] == 3).sum()) if not selected.empty else 0,
            }
        )
    for row in explode_count(features, "marker_source_class").itertuples(index=False):
        rows.append(
            {
                "method": METHOD_NAME,
                "cohort": "ALL_FINAL",
                "marker_source_class": row.value,
                "source_contrasts": "",
                "n_source_lists_in_marker_source_class": pd.NA,
                "recurrence_testable": pd.NA,
                "recurrence_k": pd.NA,
                "unique_genes_before_selection": pd.NA,
                "recurrent_candidate_genes": pd.NA,
                "singleton_candidate_genes": pd.NA,
                "nonrecurrent_candidate_genes": pd.NA,
                "selected_genes_before_global_dedup": int(row.n_genes),
                "tier1_recurrent_rows": pd.NA,
                "tier2_singleton_ranked_rows": pd.NA,
                "tier3_ranked_nonrecurrent_rows": pd.NA,
            }
        )
    return pd.DataFrame(rows)


def build_evidence_class_summary(features: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"method": METHOD_NAME, "evidence_class": row.value, "n_genes": int(row.n_genes)}
            for row in explode_count(features, "evidence_classes").itertuples(index=False)
        ]
    )


def build_sensitivity_summary(
    builds: dict[str, dict[str, object]],
    selected_rule_key: str,
    desired_min_size: int,
) -> pd.DataFrame:
    rows = []
    for rule_key, _, _ in RULE_SPECS:
        build = builds[rule_key]
        rows.append(
            {
                "method": METHOD_NAME,
                "rule_key": rule_key,
                "selected_primary_rule": rule_key == selected_rule_key,
                "desired_min_size": desired_min_size,
                "final_size": int(build["final_size"]),
                "meets_desired_min_size": int(build["final_size"]) >= desired_min_size,
                "topup_tier_activated": bool(build["topup_activated"]),
                "tier1_selected_rows": int(build["tier1_rows"]),
                "tier2_selected_rows": int(build["tier2_rows"]),
                "tier3_selected_rows": int(build["tier3_rows"]),
                "evidence_class_counts": build["evidence_class_counts"],
            }
        )
    return pd.DataFrame(rows)


def validate_outputs(
    features: pd.DataFrame,
    clean_genes: list[str],
    gene_evidence: pd.DataFrame,
    selections: pd.DataFrame,
    ranking_components: pd.DataFrame,
    recurrence_k: int,
    selected_rule_key: str,
    desired_min_size: int,
) -> tuple[bool, list[dict[str, object]]]:
    checks: list[dict[str, object]] = []
    pass_col = rule_pass_column(selected_rule_key)

    def add(name: str, passed: bool, details: str) -> None:
        checks.append({"validation_check": name, "status": "PASS" if passed else "FAIL", "details": details})

    add("no_duplicate_clean_ensembl_ids", not features["gene_id"].duplicated().any(), f"duplicates={int(features['gene_id'].duplicated().sum())}")
    add("clean_gene_list_count_matches_feature_table_unique_ids", len(clean_genes) == features["gene_id"].nunique(), f"clean_list={len(clean_genes)} unique_table={features['gene_id'].nunique()}")
    add("every_retained_gene_has_evidence_class", features["evidence_classes"].map(lambda x: bool(str(x).strip())).all(), "all evidence_classes fields are non-empty")
    add("every_retained_gene_has_rank_metadata", features["ranking_component_keys"].map(lambda x: bool(str(x).strip())).all(), "all selected genes have ranking_component_keys")
    selected_keys = set(";".join(features["ranking_component_keys"].astype(str)).split(";"))
    ranked_keys = set(ranking_components["ranking_component_key"].astype(str))
    missing_keys = {key for key in selected_keys if key and key not in ranked_keys}
    add("rank_metadata_keys_exist", not missing_keys, f"missing_keys={len(missing_keys)}")

    recurrent = selections[selections["evidence_class"].isin(RECURRENT_CLASSES)]
    recurrent_ok = recurrent.empty or (
        (recurrent["source_count"] >= recurrence_k) & (recurrent["marker_source_class_source_count"] >= 2)
    ).all()
    add("every_marker_source_recurrent_gene_has_recurrence_ge_k", bool(recurrent_ok), f"recurrent_rows={len(recurrent)} recurrence_k={recurrence_k}")

    singleton = selections[selections["evidence_class"].isin(["anchor_singleton_ranked", "isolate_singleton_ranked"])]
    singleton_ok = singleton.empty or ((singleton["marker_source_class_source_count"] == 1) & singleton[pass_col].map(bool)).all()
    add("every_singleton_source_gene_passes_selected_rank_rule", bool(singleton_ok), f"singleton_rows={len(singleton)} selected_rule={selected_rule_key}")

    topup = selections[selections["evidence_class"] == "ranked_nonrecurrent_marker"]
    topup_ok = topup.empty or (
        (topup["marker_source_class_source_count"] >= 2)
        & (topup["source_count"] < recurrence_k)
        & topup[pass_col].map(bool)
    ).all()
    add("every_ranked_nonrecurrent_marker_passes_selected_rank_rule", bool(topup_ok), f"topup_rows={len(topup)} selected_rule={selected_rule_key}")

    no_old_rescue = not features["evidence_classes"].str.contains(
        OLD_LOGIC_RESCUE_PATTERN, case=False, regex=True
    ).any()
    add("no_old_isolate_extension_or_rescue_class", bool(no_old_rescue), "searched evidence_classes")

    rbl_anchor_all = ranking_components[
        (ranking_components["cohort"] == "rbl") & (ranking_components["marker_source_class"] == "anchor")
    ]
    rbl_anchor_selected = selections[
        (selections["cohort"] == "rbl") & (selections["marker_source_class"] == "anchor")
    ]
    rbl_anchor_ok = rbl_anchor_selected.empty or (
        (rbl_anchor_selected["marker_source_class_source_count"] == 1)
        & (rbl_anchor_selected["evidence_class"] == "anchor_singleton_ranked")
    ).all()
    add("rbl_anchor_genes_are_ranked_singleton_not_automatic", bool(rbl_anchor_ok), f"raw={rbl_anchor_all['gene_id'].nunique()} selected={rbl_anchor_selected['gene_id'].nunique()}")

    rbl_iso_recurrent = selections[
        (selections["cohort"] == "rbl") & (selections["evidence_class"] == "isolate_source_recurrent")
    ]
    rbl_iso_ok = rbl_iso_recurrent.empty or (rbl_iso_recurrent["source_count"] >= recurrence_k).all()
    add("rbl_isolate_source_recurrent_genes_have_recurrence_ge_k", bool(rbl_iso_ok), f"rbl_isolate_recurrent_rows={len(rbl_iso_recurrent)} recurrence_k={recurrence_k}")

    add("direction_counts_recomputed", features["direction"].isin(["UP", "DOWN", "MIXED", "UNKNOWN"]).all(), ";".join(f"{k}={v}" for k, v in features["direction"].value_counts().sort_index().items()))
    add("cohort_contribution_counts_recomputed", features["owner_profile"].notna().all(), ";".join(f"{k}={v}" for k, v in features["owner_profile"].value_counts().sort_index().items()))
    add("marker_source_class_contribution_counts_recomputed", features["marker_source_class"].map(lambda x: bool(str(x).strip())).all(), ";".join(f"{row.value}={row.n_genes}" for row in explode_count(features, "marker_source_class").itertuples(index=False)))
    add("final_size_reported", features["gene_id"].nunique() > 0, f"final_size={features['gene_id'].nunique()} desired_min_size={desired_min_size}")
    return all(row["status"] == "PASS" for row in checks), checks


def format_count_table(counts: pd.Series) -> list[str]:
    lines = ["| Value | Genes |", "|---|---:|"]
    for value, count in counts.items():
        lines.append(f"| {value} | {int(count)} |")
    return lines


def rbl_result_counts(
    ranking_components: pd.DataFrame,
    selections: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    features: pd.DataFrame,
) -> dict[str, int]:
    rbl_features = features[features["cohort"].astype(str).str.contains("rbl", regex=False)]
    return {
        "raw_rbl_anchor_singleton_candidates": int(
            ranking_components[
                (ranking_components["cohort"] == "rbl")
                & (ranking_components["marker_source_class"] == "anchor")
                & (ranking_components["candidate_pool"] == "singleton_source")
            ]["gene_id"].nunique()
        ),
        "rbl_anchor_singleton_ranked_selected": int(
            selections[
                (selections["cohort"] == "rbl")
                & (selections["evidence_class"] == "anchor_singleton_ranked")
            ]["gene_id"].nunique()
        ),
        "rbl_isolate_source_recurrent_selected": int(
            selections[
                (selections["cohort"] == "rbl")
                & (selections["evidence_class"] == "isolate_source_recurrent")
            ]["gene_id"].nunique()
        ),
        "rbl_ranked_nonrecurrent_selected": int(
            selections[
                (selections["cohort"] == "rbl")
                & (selections["evidence_class"] == "ranked_nonrecurrent_marker")
            ]["gene_id"].nunique()
        ),
        "final_rbl_contribution": int(rbl_features["gene_id"].nunique()),
        "rbl_gene_evidence_rows": int(gene_evidence[gene_evidence["cohort"] == "rbl"]["gene_id"].nunique()),
    }


def write_run_report(
    outdir: Path,
    args: argparse.Namespace,
    profile_order: list[str],
    manifest_paths: dict[str, Path],
    manifest_tables: dict[str, pd.DataFrame],
    family_summary: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    marker_source_class_summary: pd.DataFrame,
    features: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    selections: pd.DataFrame,
    ranking_components: pd.DataFrame,
    threshold_df: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    comparison_counts: dict[str, int],
    validation_rows: list[dict[str, object]],
    selected_rule_key: str,
) -> None:
    report_path = outdir / f"{OUTPUT_PREFIX}_run_report.md"
    validation_ok = all(row["status"] == "PASS" for row in validation_rows)
    evidence_class_counts = explode_count(features, "evidence_classes")
    marker_source_class_counts = explode_count(features, "marker_source_class")
    direction_counts = features["direction"].value_counts().sort_index()
    rbl_counts = rbl_result_counts(ranking_components, selections, gene_evidence, features)
    final_size = int(features["gene_id"].nunique())
    command = " ".join(shlex.quote(arg) for arg in sys.argv)

    lines: list[str] = []
    lines.append(f"# {METHOD_NAME}")
    lines.append("")
    lines.append("## Why The Method Changed")
    lines.append("")
    lines.append("The current panel is a ranked marker-source pan-cancer panel built from marker-source recurrent markers, empirically ranked singleton-source markers, and ranked non-recurrent marker candidates selected from the existing threshold-passing marker lists.")
    lines.append("")
    lines.append("Historically, the earlier feature-panel result was isolate-rescue dominated, while a later backup result reduced that issue but became RBL singleton-anchor dominated. The current method keeps the original per-contrast anchor and isolate marker-selection rules unchanged and revises only the final feature-construction layer.")
    lines.append("")
    lines.append("## Backup Details")
    lines.append("")
    lines.append(f"- Current active backup directory: `{args.backup_dir or 'not provided'}`")
    lines.append(f"- Preserved earlier backup directory: `{args.previous177_backup_dir or 'not provided'}`")
    lines.append("- Backup checksums were verified before the active feature-space directory was cleaned for this run.")
    lines.append("")
    lines.append("## Upstream Marker Rules")
    lines.append("")
    lines.append("The upstream anchor and isolate marker-generation scripts, marker thresholds, source lists, and DEG tables were not changed by this run. The inputs remain the retained anchor marker lists, retained isolate marker lists, marker manifests, full DESeq2 tables, and the gene annotation table.")
    lines.append("")
    lines.append("## Rebuild Command")
    lines.append("")
    lines.append("```bash")
    lines.append(command)
    lines.append("```")
    lines.append("")
    lines.append("## Input Files")
    lines.append("")
    lines.append("| Cohort | Marker manifest |")
    lines.append("|---|---|")
    for profile in profile_order:
        lines.append(f"| {profile.upper()} | `{manifest_paths[profile]}` |")
    lines.append("")
    lines.append(f"- Gene annotation table: `{args.gene_annotation_tsv or 'not provided'}`")
    lines.append("")
    lines.append("## Contrast Design")
    lines.append("")
    lines.append("| Cohort | Marker-source class | Source lists | Recurrence testable | Retained unique genes | Recurrent candidates | Singleton candidates | Nonrecurrent candidates |")
    lines.append("|---|---|---:|---|---:|---:|---:|---:|")
    for row in family_summary.itertuples(index=False):
        lines.append(
            f"| {row.cohort.upper()} | {row.marker_source_class} | {row.n_source_lists_in_marker_source_class} | "
            f"{row.recurrence_testable} | {row.unique_genes_before_selection} | "
            f"{row.recurrent_candidate_genes} | {row.singleton_candidate_genes} | "
            f"{row.nonrecurrent_candidate_genes} |"
        )
    lines.append("")
    lines.append("### Retained Marker Counts Per Source")
    lines.append("")
    lines.append("| Cohort | Contrast | Marker-source class | Retained markers |")
    lines.append("|---|---|---|---:|")
    for profile in profile_order:
        for row in manifest_tables[profile].itertuples(index=False):
            contrast = str(getattr(row, "contrast"))
            lines.append(
                f"| {profile.upper()} | {contrast} | {marker_source_class_from_contrast(contrast)} | {int(getattr(row, 'n_markers'))} |"
            )
    lines.append("")
    lines.append("## Selection Logic")
    lines.append("")
    lines.append("- Tier 1: marker-source recurrent markers. Anchor or isolate recurrence is testable only when the marker-source class has at least two retained marker source lists, and a gene must appear in at least `recurrence_k` source lists.")
    lines.append("- Tier 2: ranked singleton-source markers. If a cohort/marker-source has exactly one source list, recurrence is not testable; genes enter only if they pass the selected quantile/rank rule.")
    lines.append("- Tier 3: ranked nonrecurrent marker candidates. This tier is activated only when the Tier 1 and Tier 2 panel is below the desired minimum size, and candidates must pass the same selected quantile/rank rule.")
    lines.append("- The previous unbounded isolate-rescue extension is not used.")
    lines.append("")
    lines.append("## Empirical Thresholds")
    lines.append("")
    lines.append("The rule uses cohort/marker-source candidate-pool thresholds: Q1 adjusted p-value, Q3 median absolute log2 fold-change, and a baseMean support threshold. The selected rule for this run is `" + selected_rule_key + "`.")
    lines.append("")
    lines.append("| Rule | Cohort | Marker-source class | Candidate pool | Genes | Q1 padj | Q3 abs(log2FC) | baseMean threshold | Passing genes |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|---:|")
    for row in threshold_df.itertuples(index=False):
        if row.rule_key != selected_rule_key:
            continue
        lines.append(
            f"| {row.rule_key} | {str(row.cohort).upper()} | {row.marker_source_class} | "
            f"{row.candidate_pool} | {row.n_candidate_genes} | {row.padj_q1_threshold:.4g} | "
            f"{row.abs_log2fc_q3_threshold:.4g} | {row.baseMean_threshold:.4g} | {row.n_rule_pass} |"
        )
    lines.append("")
    lines.append("## Ranking Documentation")
    lines.append("")
    lines.append("Genes are ordered by higher source count, lower adjusted p-value, higher median absolute log2 fold-change, higher median baseMean, and directionally consistent evidence before mixed-direction evidence. Adjusted p-value ranks statistical evidence after multiple-testing correction; absolute log2 fold-change ranks marker effect size; baseMean ranks expression support; source count ranks recurrence/support breadth; and direction consistency improves interpretability without automatically removing mixed-direction genes.")
    lines.append("")
    lines.append("## Final Panel")
    lines.append("")
    lines.append(f"- Final gene count: {final_size}")
    lines.append(f"- Desired minimum size: {args.desired_min_size}")
    if final_size < args.desired_min_size:
        lines.append("- The final panel remains below the desired minimum because no additional candidates passed the selected evidence-ranked rule.")
    elif final_size > args.desired_min_size:
        lines.append("- The final panel is above the desired minimum because all selected genes satisfied the documented evidence rules and were not artificially capped.")
    else:
        lines.append("- The final panel reaches the desired minimum without imposing an exact-size target.")
    lines.append("")
    lines.append("### Cohort Counts")
    lines.append("")
    lines.append("| Cohort | Final owner-profile genes | Genes with cohort evidence | Anchor selected | Isolate selected | Tier 1 rows | Tier 2 rows | Tier 3 rows |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in cohort_summary.itertuples(index=False):
        lines.append(
            f"| {row.cohort.upper()} | {row.final_owner_profile_genes_after_global_dedup} | "
            f"{row.final_genes_with_cohort_evidence} | {row.anchor_selected_genes_before_global_dedup} | "
            f"{row.isolate_selected_genes_before_global_dedup} | {row.tier1_recurrent_rows} | "
            f"{row.tier2_singleton_ranked_rows} | {row.tier3_ranked_nonrecurrent_rows} |"
        )
    lines.append("")
    lines.append("### Marker-Source-Class Counts")
    lines.append("")
    lines.extend(format_count_table(marker_source_class_counts.set_index("value")["n_genes"] if not marker_source_class_counts.empty else pd.Series(dtype=int)))
    lines.append("")
    lines.append("### Evidence-Class Counts")
    lines.append("")
    lines.append("| Evidence class | Genes |")
    lines.append("|---|---:|")
    for row in evidence_class_counts.itertuples(index=False):
        lines.append(f"| {row.value} | {row.n_genes} |")
    lines.append("")
    lines.append("### Direction Counts")
    lines.append("")
    lines.extend(format_count_table(direction_counts))
    lines.append("")
    lines.append("## RBL-Specific Result")
    lines.append("")
    lines.append(f"- Raw RBL anchor singleton candidates: {rbl_counts['raw_rbl_anchor_singleton_candidates']}")
    lines.append(f"- RBL anchor singleton-ranked genes retained: {rbl_counts['rbl_anchor_singleton_ranked_selected']}")
    lines.append(f"- RBL isolate marker-source recurrent genes retained: {rbl_counts['rbl_isolate_source_recurrent_selected']}")
    lines.append(f"- RBL ranked nonrecurrent marker candidates retained: {rbl_counts['rbl_ranked_nonrecurrent_selected']}")
    lines.append(f"- Final RBL contribution: {rbl_counts['final_rbl_contribution']}")
    lines.append("")
    lines.append("## Comparison To Previous Results")
    lines.append("")
    for label in sorted({key.split('_previous_final_size')[0] for key in comparison_counts if key.endswith("_previous_final_size")}):
        lines.append(f"### {label}")
        lines.append("")
        lines.append(f"- Previous size: {comparison_counts.get(label + '_previous_final_size', 0)}")
        lines.append(f"- Overlap retained: {comparison_counts.get(label + '_retained_from_previous', 0)}")
        lines.append(f"- Removed from previous: {comparison_counts.get(label + '_removed_from_previous', 0)}")
        lines.append(f"- Added compared with previous: {comparison_counts.get(label + '_newly_added', 0)}")
        lines.append(f"- Class changes among retained genes: {comparison_counts.get(label + '_class_changes', 0)}")
        lines.append("")
    lines.append("## Sensitivity")
    lines.append("")
    lines.append("| Rule | Selected | Final genes | Meets desired minimum | Tier 1 rows | Tier 2 rows | Tier 3 rows |")
    lines.append("|---|---|---:|---|---:|---:|---:|")
    for row in sensitivity_summary.itertuples(index=False):
        lines.append(
            f"| {row.rule_key} | {row.selected_primary_rule} | {row.final_size} | "
            f"{row.meets_desired_min_size} | {row.tier1_selected_rows} | "
            f"{row.tier2_selected_rows} | {row.tier3_selected_rows} |"
        )
    lines.append("")
    lines.append("The primary rule was selected to prioritise the strict IQR rule when it produced an adequate panel, then the strict rule with Q1 expression support if the median expression filter was too restrictive, then the relaxed IQR rule if needed for adequacy.")
    lines.append("")
    lines.append("## Validation Status")
    lines.append("")
    lines.append("| Check | Status | Details |")
    lines.append("|---|---|---|")
    for row in validation_rows:
        lines.append(f"| {row['validation_check']} | {row['status']} | {row['details']} |")
    lines.append("")
    lines.append("## Downstream Readiness")
    lines.append("")
    if validation_ok:
        lines.append("The ranked marker-source panel passed validation and is ready for downstream dry-runs and reruns.")
    else:
        lines.append("The ranked marker-source panel did not pass validation and should not be used downstream until failures are resolved.")
    lines.append("")
    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text)
    (outdir / "pan_cancer_feature_build_report.md").write_text(report_text)


def write_downstream_plan(outdir: Path) -> None:
    plan_path = outdir / "downstream_rerun_plan_from_ranked_marker_source_panel.md"
    targets_path = outdir / "downstream_rerun_targets.tsv"
    dryrun_path = outdir / "downstream_rerun_dryrun.log"
    dependency_tsv = "results/unsupervised/pan_cancer/feature_space/pan_cancer_features.tsv"
    dependency_clean = "results/unsupervised/pan_cancer/feature_space/pan_cancer_features_clean.txt"
    target_rows = [
        ("marker_framework_query_preparation", "results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/query_manifest.tsv", dependency_tsv, "Snakemake dry-run"),
        ("gprofiler_by_cohort_query_sets", "results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/gprofiler_by_cohort_all.gmt", dependency_tsv, "Snakemake dry-run"),
        ("gprofiler_enrichment", "manual_or_API_gprofiler_submission_from_query_manifests", "results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/query_manifest.tsv", "manual/API after query prep"),
        ("background_validation_tables", "results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/background_validation.tsv", dependency_clean, "Snakemake dry-run"),
        ("feature_set_composition_summaries", f"results/unsupervised/pan_cancer/feature_space/{OUTPUT_PREFIX}_by_evidence_class.tsv", dependency_tsv, "already generated by panel build"),
        ("feature_set_plots_tables", "results/unsupervised/pan_cancer/feature_space/pan_cancer_features.tsv", dependency_tsv, "already generated by panel build"),
        ("evidence_class_contribution_outputs", f"results/unsupervised/pan_cancer/feature_space/{OUTPUT_PREFIX}_by_cohort.tsv", dependency_tsv, "already generated by panel build"),
        ("thesis_export_tables", "pending_after_gprofiler_and_plot_reruns", dependency_tsv, "downstream after accepted reruns"),
    ]
    pd.DataFrame(
        [
            {
                "analysis_step": step,
                "target_path": target,
                "depends_on": depends,
                "run_type": run_type,
                "notes": "Prepared for rerun after ranked_marker_source_pan_cancer_panel acceptance",
            }
            for step, target, depends, run_type in target_rows
        ]
    ).to_csv(targets_path, sep="\t", index=False)
    lines = [
        f"# Downstream Rerun Plan From {METHOD_NAME}",
        "",
        "The downstream rerun should start from the current active feature files:",
        "",
        f"- `{dependency_tsv}`",
        f"- `{dependency_clean}`",
        "",
        "Prepared downstream areas:",
        "",
        "1. marker-framework enrichment query preparation",
        "2. g:Profiler enrichment for the new panel after query preparation; no pan-cancer g:Profiler execution rule is currently exposed in this Snakefile",
        "3. background-validation tables",
        "4. feature-set composition summaries",
        "5. final feature-set plots/tables",
        "6. evidence-class/cohort/marker-source contribution plots",
        "7. thesis-export tables generated from the pan-cancer feature set",
        "",
        "Run dry-runs before full reruns for Snakemake-backed targets. Candidate targets and manual/API steps are listed in `downstream_rerun_targets.tsv`.",
        "",
        "Suggested dry-run command:",
        "",
        "```bash",
        "snakemake -n -p --config pipeline_profile=multicohort_cancer -- results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/gprofiler_by_cohort_all.gmt results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/query_manifest.tsv results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/background_validation.tsv",
        "```",
        "",
    ]
    plan_path.write_text("\n".join(lines))
    if not dryrun_path.exists():
        dryrun_path.write_text(
            "Dry-run log placeholder created by build script. A Snakemake dry-run can overwrite this file.\n"
        )


def write_thesis_notes(outdir: Path, selected_rule_key: str, final_size: int) -> None:
    notes_path = outdir / f"{OUTPUT_PREFIX}_thesis_update_notes.md"
    lines = [
        f"# Thesis Update Notes For {METHOD_NAME}",
        "",
        "Do not edit thesis files until the ranked marker-source panel is accepted.",
        "",
        "## Methods",
        "",
        "- Clarify that original per-contrast anchor and isolate marker-selection rules were unchanged.",
        "- Update final pan-cancer feature construction to describe the ranked marker-source panel.",
        "- Define marker-source recurrence with `recurrence_k = 2`.",
        "- Define singleton-source handling as empirically ranked singleton-source marker selection.",
        "- Define ranked nonrecurrent marker candidates as a conditional Tier 3 used only when the Tier 1/Tier 2 panel is below the desired minimum.",
        "- Define quantile/IQR criteria using Q1 adjusted p-value, Q3 absolute log2 fold-change, and baseMean expression support.",
        "- Define rank ordering by source count, adjusted p-value, absolute log2 fold-change, baseMean, and direction consistency.",
        "",
        "## Supplementary Methods",
        "",
        "- Include full contrast design and marker source counts per cohort.",
        "- Include empirical threshold tables per cohort/marker-source candidate pool.",
        "- Include sensitivity comparison between strict and relaxed IQR rules.",
        "",
        "## Algorithms",
        "",
        "- Update pan-cancer feature-construction pseudocode with inputs, recurrence threshold, evidence tiers, rank ordering, empirical quantile rules, and outputs.",
        "- Keep the algorithm concise and avoid script-specific implementation details.",
        "",
        "## Results",
        "",
        f"- Update final panel size to {final_size}.",
        "- Update cohort, marker-source, evidence-class, and direction counts from the current run outputs.",
        "- Update RBL interpretation using the ranked singleton-source and isolate marker-source recurrent counts.",
        "- Update enrichment results only after downstream reruns are completed.",
        "",
        "## Discussion",
        "",
        "- Explain why final construction logic was revised after the historical feature-panel results.",
        "- Explain that the panel balances marker-source recurrence evidence and empirically ranked marker evidence.",
        "- Avoid claiming that all features are recurrent.",
        "- Discuss caveats of singleton-source markers and ranked nonrecurrent marker candidates.",
        "",
        "## Current Build",
        "",
        f"- Method: `{METHOD_NAME}`",
        f"- Selected rule: `{selected_rule_key}`",
        f"- Current feature count: {final_size}",
        "",
    ]
    notes_path.write_text("\n".join(lines))


def write_active_directory_manifest(outdir: Path, descriptions: dict[str, tuple[str, str]]) -> None:
    manifest_path = outdir / f"{OUTPUT_PREFIX}_active_directory_manifest.tsv"
    rows: list[dict[str, object]] = []
    for path in sorted(outdir.iterdir(), key=lambda p: p.name):
        if not path.is_file():
            continue
        category, description = descriptions.get(
            path.name,
            ("current active output", "current active output from ranked marker-source panel"),
        )
        rows.append(file_manifest_row(path, category, description))
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    rows = []
    for path in sorted(outdir.iterdir(), key=lambda p: p.name):
        if path.is_file():
            category, description = descriptions.get(
                path.name,
                ("current active output", "current active output from ranked marker-source panel"),
            )
            rows.append(file_manifest_row(path, category, description))
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)


def write_run_manifest(
    outdir: Path,
    args: argparse.Namespace,
    profile_order: list[str],
    manifest_paths: dict[str, Path],
    manifest_tables: dict[str, pd.DataFrame],
    previous_feature_paths: dict[str, Path],
    output_descriptions: dict[str, tuple[str, str]],
) -> None:
    rows: list[dict[str, object]] = []
    for profile in profile_order:
        rows.append(file_manifest_row(manifest_paths[profile], "input", f"{profile} marker manifest"))
        deseq2_dir = manifest_paths[profile].parent.parent
        for row in manifest_tables[profile].itertuples(index=False):
            marker_path = resolve_relative_path(deseq2_dir, getattr(row, "marker_file"))
            table_path = resolve_relative_path(deseq2_dir, getattr(row, "table_file"))
            rows.append(file_manifest_row(marker_path, "input", f"{profile} retained marker list"))
            rows.append(file_manifest_row(table_path, "input", f"{profile} full DESeq2 table"))
    if args.gene_annotation_tsv:
        rows.append(file_manifest_row(Path(args.gene_annotation_tsv), "input", "gene annotation table"))
    for label, path in previous_feature_paths.items():
        rows.append(file_manifest_row(path, "input", f"{label} previous feature table for comparison"))
    if args.backup_dir:
        rows.append(file_manifest_row(Path(args.backup_dir) / "BACKUP_MANIFEST.tsv", "backup", "active backup manifest"))
    for name, (category, description) in output_descriptions.items():
        path = outdir / name
        if path.exists():
            rows.append(file_manifest_row(path, category, description))
    pd.DataFrame(rows).to_csv(outdir / f"{OUTPUT_PREFIX}_run_manifest.tsv", sep="\t", index=False)


def write_outputs(
    outdir: Path,
    args: argparse.Namespace,
    profile_order: list[str],
    manifest_paths: dict[str, Path],
    manifest_tables: dict[str, pd.DataFrame],
    evidence_by_profile: dict[str, pd.DataFrame],
    ranking_components: pd.DataFrame,
    family_summary: pd.DataFrame,
    threshold_df: pd.DataFrame,
    builds: dict[str, dict[str, object]],
    selected_rule_key: str,
    previous_feature_paths: dict[str, Path],
    annotation_path: str,
) -> tuple[dict[str, int], list[dict[str, object]]]:
    outdir.mkdir(parents=True, exist_ok=True)
    selected_build = builds[selected_rule_key]
    selections = selected_build["selections"]
    gene_evidence = selected_build["gene_evidence"]
    features = selected_build["features"]
    duplicate_candidate_count = int(selected_build["duplicate_candidate_count"])
    removed_ribo_mt = selected_build["removed_ribo_mt"]
    topup_activated = bool(selected_build["topup_activated"])

    features_path = outdir / "pan_cancer_features.tsv"
    clean_path = outdir / "pan_cancer_features_clean.txt"
    features.to_csv(features_path, sep="\t", index=False)
    clean_genes = features["gene_id"].tolist()
    write_gene_list(clean_path, clean_genes)
    for direction in ["UP", "DOWN", "MIXED"]:
        write_gene_list(
            outdir / f"pan_cancer_features.{direction}.txt",
            features.loc[features["direction"] == direction, "gene_id"].tolist(),
        )
    (outdir / "pan_cancer_features_done.txt").write_text("done\n")

    gene_evidence.to_csv(outdir / "pan_cancer_feature_gene_evidence.tsv", sep="\t", index=False)
    ranking_components.to_csv(outdir / f"{OUTPUT_PREFIX}_ranking_components.tsv", sep="\t", index=False)
    threshold_df.to_csv(outdir / f"{OUTPUT_PREFIX}_quantile_thresholds.tsv", sep="\t", index=False)
    sensitivity_summary = build_sensitivity_summary(builds, selected_rule_key, args.desired_min_size)
    sensitivity_summary.to_csv(outdir / f"{OUTPUT_PREFIX}_sensitivity_summary.tsv", sep="\t", index=False)
    selections.to_csv(outdir / f"{OUTPUT_PREFIX}_selected_marker_source_class_rows.tsv", sep="\t", index=False)

    cohort_summary = build_cohort_summary(
        profile_order, features, gene_evidence, selections, family_summary
    )
    cohort_summary.to_csv(outdir / f"{OUTPUT_PREFIX}_by_cohort.tsv", sep="\t", index=False)
    marker_source_class_summary = build_marker_source_class_summary_output(features, selections, family_summary)
    marker_source_class_summary.to_csv(outdir / f"{OUTPUT_PREFIX}_by_marker_source_class.tsv", sep="\t", index=False)
    evidence_class_summary = build_evidence_class_summary(features)
    evidence_class_summary.to_csv(outdir / f"{OUTPUT_PREFIX}_by_evidence_class.tsv", sep="\t", index=False)

    comparison_counts: dict[str, int] = {}
    for label, path in previous_feature_paths.items():
        comparison_counts.update(build_previous_comparison(label, path, features, outdir))

    validation_ok, validation_rows = validate_outputs(
        features,
        clean_genes,
        gene_evidence,
        selections,
        ranking_components,
        args.recurrence_k,
        selected_rule_key,
        args.desired_min_size,
    )
    pd.DataFrame(validation_rows).to_csv(outdir / f"{OUTPUT_PREFIX}_validation.tsv", sep="\t", index=False)

    summary = build_summary(
        features,
        gene_evidence,
        selections,
        family_summary,
        comparison_counts,
        duplicate_candidate_count,
        removed_ribo_mt,
        selected_rule_key,
        args.desired_min_size,
        topup_activated,
    )
    summary.to_csv(outdir / "pan_cancer_feature_build_summary.tsv", sep="\t", index=False)

    write_run_report(
        outdir,
        args,
        profile_order,
        manifest_paths,
        manifest_tables,
        family_summary,
        cohort_summary,
        marker_source_class_summary,
        features,
        gene_evidence,
        selections,
        ranking_components,
        threshold_df,
        sensitivity_summary,
        comparison_counts,
        validation_rows,
        selected_rule_key,
    )
    write_downstream_plan(outdir)
    write_thesis_notes(outdir, selected_rule_key, int(features["gene_id"].nunique()))

    output_descriptions = {
        "pan_cancer_features.tsv": ("current active output", "current ranked marker-source feature table"),
        "pan_cancer_features_clean.txt": ("compatibility output", "ordered clean Ensembl gene list"),
        "pan_cancer_features.UP.txt": ("compatibility output", "UP-direction feature subset"),
        "pan_cancer_features.DOWN.txt": ("compatibility output", "DOWN-direction feature subset"),
        "pan_cancer_features.MIXED.txt": ("compatibility output", "MIXED-direction feature subset"),
        "pan_cancer_features_done.txt": ("compatibility output", "Snakemake completion sentinel"),
        "pan_cancer_feature_build_summary.tsv": ("compatibility output", "machine-readable build summary"),
        "pan_cancer_feature_gene_evidence.tsv": ("current active output", "per-gene selected evidence table"),
        "pan_cancer_feature_build_report.md": ("report", "compatibility copy of ranked marker-source run report"),
        f"{OUTPUT_PREFIX}_by_cohort.tsv": ("current active output", "cohort contribution summary"),
        f"{OUTPUT_PREFIX}_by_marker_source_class.tsv": ("current active output", "marker-source contribution summary"),
        f"{OUTPUT_PREFIX}_by_evidence_class.tsv": ("current active output", "evidence-class contribution summary"),
        f"{OUTPUT_PREFIX}_ranking_components.tsv": ("current active output", "rank and quantile components for all candidates"),
        f"{OUTPUT_PREFIX}_quantile_thresholds.tsv": ("current active output", "empirical quantile thresholds"),
        f"{OUTPUT_PREFIX}_sensitivity_summary.tsv": ("current active output", "strict and relaxed IQR sensitivity results"),
        f"{OUTPUT_PREFIX}_selected_marker_source_class_rows.tsv": ("current active output", "selected marker-source candidate rows"),
        f"{OUTPUT_PREFIX}_validation.tsv": ("current active output", "validation check table"),
        f"{OUTPUT_PREFIX}_run_report.md": ("report", "ranked marker-source run report"),
        f"{OUTPUT_PREFIX}_run_manifest.tsv": ("current active output", "input/output run manifest"),
        f"{OUTPUT_PREFIX}_active_directory_manifest.tsv": ("current active output", "active directory file manifest"),
        "downstream_rerun_plan_from_ranked_marker_source_panel.md": ("downstream plan", "downstream rerun plan"),
        "downstream_rerun_targets.tsv": ("downstream plan", "downstream rerun targets"),
        "downstream_rerun_dryrun.log": ("downstream plan", "downstream dry-run log"),
        f"{OUTPUT_PREFIX}_thesis_update_notes.md": ("downstream plan", "thesis update notes for later editing"),
    }
    for label in previous_feature_paths:
        output_descriptions[f"{OUTPUT_PREFIX}_removed_vs_{label}.tsv"] = ("current active output", f"genes removed compared with {label}")
        output_descriptions[f"{OUTPUT_PREFIX}_added_vs_{label}.tsv"] = ("current active output", f"genes added compared with {label}")
        output_descriptions[f"{OUTPUT_PREFIX}_overlap_vs_{label}.tsv"] = ("current active output", f"overlap compared with {label}")

    write_run_manifest(
        outdir,
        args,
        profile_order,
        manifest_paths,
        manifest_tables,
        previous_feature_paths,
        output_descriptions,
    )
    write_active_directory_manifest(outdir, output_descriptions)

    if not validation_ok:
        for row in validation_rows:
            if row["status"] != "PASS":
                print(f"[VALIDATION FAIL] {row['validation_check']}: {row['details']}", file=sys.stderr)
        sys.exit("[ERROR] ranked marker-source panel failed validation")
    return comparison_counts, validation_rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=f"Build {METHOD_NAME} from retained DESeq2 marker lists."
    )
    parser.add_argument("--profile-marker-dir", action="append", default=[], metavar="PROFILE=PATH")
    parser.add_argument("--profile-marker-manifest", action="append", default=[], metavar="PROFILE=PATH")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method", default=METHOD_NAME)
    parser.add_argument("--recurrence-k", type=int, default=2)
    parser.add_argument("--desired-min-size", type=int, default=200)
    parser.add_argument("--force-exact-size", action="store_true", default=False)
    parser.add_argument("--marker-source-recurrence", action="store_true", default=False)
    parser.add_argument(
        "--singleton-source-policy",
        choices=["ranked_quantile", "empirical_ranked", "keep_as_single_source"],
        default="ranked_quantile",
    )
    parser.add_argument("--disable-old-isolate-rescue", action="store_true", default=False)
    parser.add_argument(
        "--cap-isolate",
        type=int,
        default=0,
        help="Deprecated compatibility option. Ignored because old isolate rescue is disabled.",
    )
    parser.add_argument("--primary-rule", choices=["auto", "strict_iqr", "relaxed_iqr"], default="auto")
    parser.add_argument("--strict-basemean-threshold", choices=["median", "q1"], default="median")
    parser.add_argument("--padj-quantile", type=float, default=0.25)
    parser.add_argument("--abs-log2fc-quantile", type=float, default=0.75)
    parser.add_argument("--relaxed-minimum-criteria", type=int, default=2)
    parser.add_argument("--remove-ribo-mt", action="store_true", default=False)
    parser.add_argument("--gene-annotation-tsv", default="")
    parser.add_argument("--previous-feature-set", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--previous-features", default="", help="Deprecated alias for previous177=PATH")
    parser.add_argument("--backup-dir", default="")
    parser.add_argument("--previous177-backup-dir", default="")
    args = parser.parse_args()

    if args.method != METHOD_NAME:
        sys.exit(f"[ERROR] --method must be {METHOD_NAME}")
    if not args.marker_source_recurrence:
        sys.exit("[ERROR] Refusing to run old logic: pass --marker-source-recurrence")
    if args.singleton_source_policy not in {"ranked_quantile", "empirical_ranked"}:
        sys.exit("[ERROR] This build requires --singleton-source-policy ranked_quantile")
    if not args.disable_old_isolate_rescue:
        sys.exit("[ERROR] This build requires --disable-old-isolate-rescue")
    if args.force_exact_size:
        sys.exit("[ERROR] This method does not support forced exact-size panels")
    if args.recurrence_k < 1:
        sys.exit("[ERROR] --recurrence-k must be >= 1")
    if args.desired_min_size < 1:
        sys.exit("[ERROR] --desired-min-size must be >= 1")
    if not (0 < args.padj_quantile < 1):
        sys.exit("[ERROR] --padj-quantile must be between 0 and 1")
    if not (0 < args.abs_log2fc_quantile < 1):
        sys.exit("[ERROR] --abs-log2fc-quantile must be between 0 and 1")
    if args.relaxed_minimum_criteria < 1:
        sys.exit("[ERROR] --relaxed-minimum-criteria must be >= 1")

    marker_order, marker_dirs = parse_profile_paths(args.profile_marker_dir, "--profile-marker-dir")
    manifest_order, marker_manifests = parse_profile_paths(args.profile_marker_manifest, "--profile-marker-manifest")
    profile_order = manifest_order or marker_order
    if not profile_order:
        sys.exit("[ERROR] At least one profile marker manifest or marker directory is required")

    outdir = Path(args.output_dir)
    annotation = load_gene_annotation(args.gene_annotation_tsv)
    previous_feature_paths = parse_labeled_paths(args.previous_feature_set)
    if args.previous_features and "previous177" not in previous_feature_paths:
        previous_feature_paths["previous177"] = Path(args.previous_features)
    previous_feature_paths.setdefault("previous177", Path("__missing_previous177__"))
    previous_feature_paths.setdefault("previous125", Path("__missing_previous125__"))

    evidence_by_profile: dict[str, pd.DataFrame] = {}
    manifest_paths: dict[str, Path] = {}
    manifest_tables: dict[str, pd.DataFrame] = {}
    for profile in profile_order:
        marker_dir = marker_dirs.get(profile)
        manifest_path = marker_manifests.get(profile)
        print(f"[INFO] Loading retained marker evidence for {profile}", file=sys.stderr)
        evidence, resolved_manifest, manifest = load_profile_marker_evidence(
            profile, marker_dir, manifest_path
        )
        evidence_by_profile[profile] = evidence
        manifest_paths[profile] = resolved_manifest
        manifest_tables[profile] = manifest

    ranking_components, family_summary, threshold_df = build_ranking_components(
        profile_order,
        evidence_by_profile,
        args.recurrence_k,
        args.padj_quantile,
        args.abs_log2fc_quantile,
        args.relaxed_minimum_criteria,
    )

    builds = {
        rule_key: build_panel_for_rule(
            profile_order,
            ranking_components,
            evidence_by_profile,
            annotation,
            args.remove_ribo_mt,
            args.desired_min_size,
            args.force_exact_size,
            rule_key,
        )
        for rule_key, _, _ in RULE_SPECS
    }
    selected_rule_key = choose_primary_build(builds, args)
    comparison_counts, validation_rows = write_outputs(
        outdir,
        args,
        profile_order,
        manifest_paths,
        manifest_tables,
        evidence_by_profile,
        ranking_components,
        family_summary,
        threshold_df,
        builds,
        selected_rule_key,
        previous_feature_paths,
        args.gene_annotation_tsv,
    )
    selected_size = builds[selected_rule_key]["final_size"]
    print(
        f"[OK] Wrote {METHOD_NAME} to {outdir} with {selected_size} genes "
        f"using {selected_rule_key}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
