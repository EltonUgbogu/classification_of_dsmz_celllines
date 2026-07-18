#!/usr/bin/env python3
"""Construct canonical graph-derived functional-enrichment query inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


PRIMARY_ANALYSIS_METHOD = "unordered_overrepresentation"
PRIMARY_RANKING_USED_FOR_ENRICHMENT = False
PRIMARY_ORDERED = False
RANKING_PROVENANCE_ROLE = "analytical_prioritisation_provenance"
PAN_CANCER_BACKGROUND_ID = "pan_cancer_feature_selection_background_genes"
ALLOWED_CONTRAST_TYPES = {
    "isolate_focal_vs_other_same_cancer",
    "anchor_focal_vs_outside_focal_component",
}
ALLOWED_MARKER_EVIDENCE_STRATA = {
    "isolate_associated",
    "anchor_associated",
}
MANIFEST_COLUMNS = [
    "analysis_method",
    "query_id",
    "query_output_directory_id",
    "query_family",
    "cancer_type",
    "attributed_cancer_type",
    "marker_evidence_stratum",
    "marker_evidence_source_class",
    "evidence_class",
    "contrast_id",
    "contrast_type",
    "direction",
    "direction_consistency_class",
    "gene_count",
    "gene_list_path",
    "background_id",
    "background_type",
    "background_count",
    "background_path",
    "background_sha256",
    "background_definition",
    "ordered",
    "ranking_used_for_enrichment",
    "ranked_genes_path",
    "rank_source",
    "ranking_role",
    "minimum_query_gene_count",
    "minimum_base_mean",
    "query_execution_status",
    "query_exclusion_reason",
]


@dataclass(frozen=True)
class BackgroundRecord:
    background_id: str
    background_type: str
    background_definition: str
    genes: tuple[str, ...]
    stage_path: Path
    final_path: Path
    sha256: str


@dataclass(frozen=True)
class ContrastRecord:
    cancer_type: str
    contrast_id: str
    contrast_type: str
    marker_evidence_stratum: str
    marker_genes: tuple[str, ...]
    background_genes: tuple[str, ...]
    ranks: pd.DataFrame
    minimum_base_mean: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Construct canonical graph-derived functional-enrichment queries."
    )
    parser.add_argument("--pipeline-root", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--selected-feature-genes", required=True)
    parser.add_argument("--eligible-gene-background", required=True)
    parser.add_argument(
        "--contrast-manifest",
        action="append",
        default=[],
        help="Cancer-type-labelled canonical manifest, formatted cancer_type=path.",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--minimum-query-gene-count", required=True, type=int)
    parser.add_argument("--enabled-query-families", required=True)
    parser.add_argument("--recurrence-minimum-retained-contrast-count", required=True, type=int)
    return parser.parse_args()


def strip_ensembl_version_suffix(value: object) -> str:
    return re.sub(r"\.[0-9]+$", "", str(value).strip())


def validate_canonical_gene_identifier_uniqueness(
    original_gene_ids: Iterable[object],
    source_file: Path,
    query_id: str = "",
) -> list[str]:
    canonical_gene_ids = canonicalize_gene_ids_preserve_rows(original_gene_ids, source_file, query_id)
    return list(dict.fromkeys(canonical_gene_ids))


def canonicalize_gene_ids_preserve_rows(
    original_gene_ids: Iterable[object],
    source_file: Path,
    query_id: str = "",
) -> list[str]:
    mapping: dict[str, set[str]] = defaultdict(set)
    canonical_gene_ids: list[str] = []
    for original in original_gene_ids:
        original_text = str(original).strip()
        if not original_text or original_text.lower() == "nan":
            continue
        canonical_gene_id = strip_ensembl_version_suffix(original_text)
        if not canonical_gene_id:
            continue
        mapping[canonical_gene_id].add(original_text)
        canonical_gene_ids.append(canonical_gene_id)
    collisions = {
        canonical: sorted(originals)
        for canonical, originals in mapping.items()
        if len(originals) > 1
    }
    if collisions:
        details = []
        for canonical, originals in sorted(collisions.items()):
            details.append(
                f"source_file={source_file} query_id_or_contrast_id={query_id or ''} "
                f"canonical_gene_id={canonical} original_identifiers={','.join(originals)}"
            )
        raise ValueError(
            "Ensembl version suffix removal produced gene identifier collisions:\n"
            + "\n".join(details)
        )
    return canonical_gene_ids


def read_gene_file(path: Path, source_label: str = "") -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing gene list {source_label}: {path}")
    values: list[str] = []
    with path.open() as handle:
        for line in handle:
            value = line.strip().split("\t", 1)[0]
            if value and value != "gene_id":
                values.append(value)
    return validate_canonical_gene_identifier_uniqueness(values, path, source_label)


def write_gene_file(path: Path, genes: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    unique_genes = list(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
    path.write_text(("\n".join(unique_genes) + "\n") if unique_genes else "")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sanitize_identifier(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    text = re.sub(r"_+", "_", text).strip("_")
    if not text:
        raise ValueError("Cannot create identifier from an empty value")
    return text


def resolve_path(path: object, base_dir: Path) -> Path:
    value = str(path).strip()
    if not value or value.lower() == "nan":
        raise ValueError(f"Missing required path relative to {base_dir}")
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return base_dir / candidate


def path_relative_to_root(path: Path, pipeline_root: Path) -> str:
    absolute_path = path.resolve()
    try:
        return str(absolute_path.relative_to(pipeline_root.resolve()))
    except ValueError:
        return os.path.relpath(absolute_path, pipeline_root.resolve())


def parse_enabled_query_families(value: str) -> set[str]:
    enabled = {item.strip() for item in value.split(",") if item.strip()}
    valid = {
        "complete_pan_cancer_feature_panel",
        "cancer_type_attributed_feature_set",
        "marker_evidence_stratum",
        "evidence_class",
        "direction_consistency_class",
        "contrast_level_marker_set",
        "enrichment_contrast_marker_recurrence",
    }
    unknown = sorted(enabled - valid)
    if unknown:
        raise ValueError(f"Unknown enabled query families: {', '.join(unknown)}")
    return enabled


def read_table(path: Path, required_columns: Iterable[str], table_name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {table_name}: {path}")
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [col for col in required_columns if col not in table.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns {missing}: {path}")
    return table


def split_multivalue(value: object, source_column: str, row_label: str) -> list[str]:
    parts = [part.strip() for part in str(value).split(";")]
    parts = [part for part in parts if part]
    if not parts:
        raise ValueError(f"{source_column} has no observed value for {row_label}")
    return parts


def load_feature_table(path: Path) -> pd.DataFrame:
    required = [
        "gene_id",
        "attributed_cancer_type",
        "marker_evidence_stratum",
        "evidence_classes",
        "direction_consistency_class",
        "direction",
    ]
    features = read_table(path, required, "pan-cancer feature table")
    features = features.copy()
    features["gene_id"] = canonicalize_gene_ids_preserve_rows(
        features["gene_id"], path, "pan_cancer_features"
    )
    if features["gene_id"].duplicated().any():
        duplicates = sorted(features.loc[features["gene_id"].duplicated(), "gene_id"].unique())
        raise ValueError(
            "pan_cancer_features.tsv contains duplicate canonical genes: "
            + ",".join(duplicates[:50])
        )
    for column in ["attributed_cancer_type", "direction_consistency_class", "direction"]:
        bad = features[column].astype(str).str.strip() == ""
        if bad.any():
            raise ValueError(f"pan_cancer_features.tsv has empty {column} values")
        features[column] = features[column].astype(str).str.strip()
    return features


def rank_rows_from_features(features: pd.DataFrame, genes: Iterable[str]) -> pd.DataFrame | None:
    if "selection_rank" not in features.columns:
        return None
    subset = features[features["gene_id"].isin(set(genes))][["gene_id", "selection_rank"]].copy()
    if subset.empty:
        return None
    subset["rank_stat"] = pd.to_numeric(subset["selection_rank"], errors="coerce")
    subset = subset.dropna(subset=["rank_stat"])
    if subset.empty:
        return None
    return subset[["gene_id", "rank_stat"]].drop_duplicates("gene_id")


def load_contrast_records(
    manifest_specs: list[str],
    pipeline_root: Path,
) -> list[ContrastRecord]:
    records: list[ContrastRecord] = []
    for spec in manifest_specs:
        if "=" not in spec:
            raise ValueError(f"--contrast-manifest must be cancer_type=path, got: {spec}")
        cancer_type, raw_manifest_path = spec.split("=", 1)
        cancer_type = cancer_type.strip()
        if not cancer_type:
            raise ValueError(f"Empty cancer type in --contrast-manifest {spec}")
        manifest_path = resolve_path(raw_manifest_path, pipeline_root)
        manifest = read_table(
            manifest_path,
            [
                "cancer_type",
                "contrast_id",
                "contrast_type",
                "marker_evidence_stratum",
                "marker_table_path",
                "marker_gene_list_path",
                "result_table_path",
                "n_markers_after_cap",
                "minimum_base_mean",
            ],
            "canonical contrast-level marker manifest",
        )
        manifest_dir = manifest_path.parent.parent
        unknown_contrast_types = sorted(set(manifest["contrast_type"]) - ALLOWED_CONTRAST_TYPES)
        if unknown_contrast_types:
            raise ValueError(
                f"{manifest_path} has unrecognized contrast_type values: "
                + ",".join(unknown_contrast_types)
            )
        unknown_strata = sorted(set(manifest["marker_evidence_stratum"]) - ALLOWED_MARKER_EVIDENCE_STRATA)
        if unknown_strata:
            raise ValueError(
                f"{manifest_path} has unrecognized marker_evidence_stratum values: "
                + ",".join(unknown_strata)
            )
        for row in manifest.itertuples(index=False):
            contrast_id = str(row.contrast_id)
            marker_gene_list_path = resolve_path(row.marker_gene_list_path, manifest_dir)
            marker_table_path = resolve_path(row.marker_table_path, manifest_dir)
            result_table_path = resolve_path(row.result_table_path, manifest_dir)
            marker_genes = tuple(read_gene_file(marker_gene_list_path, contrast_id))
            marker_table = read_table(
                marker_table_path,
                [
                    "gene_id",
                    "log2_fold_change_shrunken",
                    "absolute_shrunken_log2_fold_change",
                ],
                f"retained marker table for {contrast_id}",
            )
            marker_table = marker_table.copy()
            marker_table["gene_id"] = canonicalize_gene_ids_preserve_rows(
                marker_table["gene_id"], marker_table_path, contrast_id
            )
            marker_table["rank_stat"] = pd.to_numeric(
                marker_table["absolute_shrunken_log2_fold_change"], errors="coerce"
            )
            ranks = (
                marker_table[marker_table["gene_id"].isin(marker_genes)][["gene_id", "rank_stat"]]
                .dropna(subset=["rank_stat"])
                .drop_duplicates("gene_id")
            )
            result_table = read_table(
                result_table_path,
                [
                    "gene_id",
                    "adjusted_p_value",
                    "baseMean",
                    "log2_fold_change_shrunken",
                    "absolute_shrunken_log2_fold_change",
                ],
                f"post-prefilter contrast result table for {contrast_id}",
            )
            result_table = result_table.copy()
            result_table["gene_id"] = canonicalize_gene_ids_preserve_rows(
                result_table["gene_id"], result_table_path, contrast_id
            )
            minimum_base_mean = float(row.minimum_base_mean)
            adjusted = pd.to_numeric(result_table["adjusted_p_value"], errors="coerce")
            base_mean = pd.to_numeric(result_table["baseMean"], errors="coerce")
            background_genes = tuple(
                result_table.loc[
                    adjusted.notna() & base_mean.notna() & (base_mean >= minimum_base_mean),
                    "gene_id",
                ].drop_duplicates()
            )
            absent = sorted(set(marker_genes) - set(background_genes))
            if absent:
                raise ValueError(
                    "Contrast marker query is not a subset of its analytical background: "
                    f"cancer_type={cancer_type} contrast_id={contrast_id} "
                    f"query_size={len(marker_genes)} background_size={len(background_genes)} "
                    f"absent_genes={','.join(absent[:100])}"
                )
            records.append(
                ContrastRecord(
                    cancer_type=cancer_type,
                    contrast_id=contrast_id,
                    contrast_type=str(row.contrast_type),
                    marker_evidence_stratum=str(row.marker_evidence_stratum),
                    marker_genes=marker_genes,
                    background_genes=background_genes,
                    ranks=ranks,
                    minimum_base_mean=minimum_base_mean,
                )
            )
    return records


def prepare_staging(outdir: Path) -> Path:
    staging = outdir.parent / f".{outdir.name}.staging.{os.getpid()}"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    for subdir in ["queries", "backgrounds", "rankings"]:
        (staging / subdir).mkdir(parents=True, exist_ok=True)
    return staging


def publish_staging(staging: Path, outdir: Path) -> None:
    if outdir.exists():
        shutil.rmtree(outdir)
    staging.replace(outdir)


def main() -> int:
    args = parse_args()
    if args.minimum_query_gene_count < 1:
        raise ValueError("--minimum-query-gene-count must be at least 1")
    if args.recurrence_minimum_retained_contrast_count < 1:
        raise ValueError("--recurrence-minimum-retained-contrast-count must be at least 1")

    pipeline_root = Path(args.pipeline_root).resolve()
    outdir = resolve_path(args.outdir, pipeline_root)
    enabled_families = parse_enabled_query_families(args.enabled_query_families)

    features_path = resolve_path(args.features, pipeline_root)
    selected_gene_path = resolve_path(args.selected_feature_genes, pipeline_root)
    eligible_background_path = resolve_path(args.eligible_gene_background, pipeline_root)

    features = load_feature_table(features_path)
    selected_feature_genes = read_gene_file(selected_gene_path, "pan_cancer_selected_feature_gene_ids")
    eligible_background_genes = read_gene_file(
        eligible_background_path,
        "pan_cancer_feature_selection_eligible_gene_background",
    )
    feature_genes = list(features["gene_id"])

    missing_features_from_selected = sorted(set(feature_genes) - set(selected_feature_genes))
    missing_selected_from_features = sorted(set(selected_feature_genes) - set(feature_genes))
    if len(feature_genes) != len(selected_feature_genes) or missing_features_from_selected or missing_selected_from_features:
        raise ValueError(
            "Pan-cancer selected feature identity validation failed before query construction: "
            f"feature_table_count={len(feature_genes)} "
            f"selected_feature_gene_file_count={len(selected_feature_genes)} "
            f"missing_from_selected_file={','.join(missing_features_from_selected[:100])} "
            f"missing_from_feature_table={','.join(missing_selected_from_features[:100])}"
        )
    absent_from_background = sorted(set(feature_genes) - set(eligible_background_genes))
    if absent_from_background:
        raise ValueError(
            "Selected pan-cancer features are absent from the explicit eligible-gene background: "
            + ",".join(absent_from_background[:100])
        )

    contrast_records = load_contrast_records(args.contrast_manifest, pipeline_root)
    staging = prepare_staging(outdir)

    background_registry: dict[str, BackgroundRecord] = {}
    manifest_rows: list[dict[str, object]] = []
    integrity_rows: list[dict[str, object]] = []
    recurrence_rows: list[dict[str, object]] = []

    def stage_path_for(final_path: Path) -> Path:
        return staging / final_path.relative_to(outdir)

    def register_background(
        background_id: str,
        background_type: str,
        background_definition: str,
        genes: Iterable[str],
    ) -> BackgroundRecord:
        background_id = sanitize_identifier(background_id)
        genes_tuple = tuple(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
        final_path = outdir / "backgrounds" / f"{background_id}.tsv"
        stage_path = stage_path_for(final_path)
        if background_id in background_registry:
            previous = background_registry[background_id]
            if previous.genes != genes_tuple:
                raise ValueError(
                    f"background_id={background_id} refers to different gene sets across queries"
                )
            return previous
        write_gene_file(stage_path, genes_tuple)
        record = BackgroundRecord(
            background_id=background_id,
            background_type=background_type,
            background_definition=background_definition,
            genes=genes_tuple,
            stage_path=stage_path,
            final_path=final_path,
            sha256=sha256_file(stage_path),
        )
        background_registry[background_id] = record
        return record

    pan_background = register_background(
        PAN_CANCER_BACKGROUND_ID,
        "pan_cancer_feature_selection_eligible_gene_background",
        "Unique canonical genes eligible for graph-derived pan-cancer feature selection.",
        eligible_background_genes,
    )

    def add_integrity(check_name: str, status: str, query_id: str = "", details: str = "") -> None:
        integrity_rows.append(
            {
                "check_name": check_name,
                "status": status,
                "query_id": query_id,
                "details": details,
            }
        )

    def add_query(
        query_id: str,
        query_family: str,
        genes: Iterable[str],
        background: BackgroundRecord,
        cancer_type: str = "",
        attributed_cancer_type: str = "",
        marker_evidence_stratum: str = "",
        marker_evidence_source_class: str = "",
        evidence_class: str = "",
        contrast_id: str = "",
        contrast_type: str = "",
        direction: str = "",
        direction_consistency_class: str = "",
        ranks: pd.DataFrame | None = None,
        rank_source: str = "",
        minimum_base_mean: object = "",
    ) -> None:
        query_id = sanitize_identifier(query_id)
        query_output_directory_id = sanitize_identifier(query_id)
        genes_tuple = tuple(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
        absent = sorted(set(genes_tuple) - set(background.genes))
        if absent:
            raise ValueError(
                "Query genes are absent from analytical background: "
                f"cancer_type={cancer_type} contrast_id={contrast_id} query_id={query_id} "
                f"query_size={len(genes_tuple)} background_size={len(background.genes)} "
                f"absent_genes={','.join(absent[:100])}"
            )
        query_final_path = outdir / "queries" / query_output_directory_id / "genes.tsv"
        query_stage_path = stage_path_for(query_final_path)
        write_gene_file(query_stage_path, genes_tuple)
        ranked_final_path = Path("")
        ranked_stage_path = Path("")
        if ranks is not None and not ranks.empty:
            rank_table = ranks[ranks["gene_id"].isin(genes_tuple)].copy()
            rank_table = rank_table.drop_duplicates("gene_id")
            if set(rank_table["gene_id"]) != set(genes_tuple):
                missing_rank = sorted(set(genes_tuple) - set(rank_table["gene_id"]))
                raise ValueError(
                    f"Ranking provenance for {query_id} is missing query genes: "
                    + ",".join(missing_rank[:100])
                )
            ranked_final_path = outdir / "rankings" / query_output_directory_id / "ranked_genes.tsv"
            ranked_stage_path = stage_path_for(ranked_final_path)
            ranked_stage_path.parent.mkdir(parents=True, exist_ok=True)
            rank_table[["gene_id", "rank_stat"]].to_csv(ranked_stage_path, sep="\t", index=False)
        gene_count = len(genes_tuple)
        if gene_count == 0:
            status = "excluded_zero_gene_query"
            reason = "zero_gene_query"
        elif gene_count < args.minimum_query_gene_count:
            status = "excluded_below_minimum_query_gene_count"
            reason = "below_minimum_query_gene_count"
        else:
            status = "runnable"
            reason = ""
        row = {
            "analysis_method": PRIMARY_ANALYSIS_METHOD,
            "query_id": query_id,
            "query_output_directory_id": query_output_directory_id,
            "query_family": query_family,
            "cancer_type": cancer_type,
            "attributed_cancer_type": attributed_cancer_type,
            "marker_evidence_stratum": marker_evidence_stratum,
            "marker_evidence_source_class": marker_evidence_source_class,
            "evidence_class": evidence_class,
            "contrast_id": contrast_id,
            "contrast_type": contrast_type,
            "direction": direction,
            "direction_consistency_class": direction_consistency_class,
            "gene_count": gene_count,
            "gene_list_path": path_relative_to_root(query_final_path, pipeline_root),
            "background_id": background.background_id,
            "background_type": background.background_type,
            "background_count": len(background.genes),
            "background_path": path_relative_to_root(background.final_path, pipeline_root),
            "background_sha256": background.sha256,
            "background_definition": background.background_definition,
            "ordered": str(PRIMARY_ORDERED).upper(),
            "ranking_used_for_enrichment": str(PRIMARY_RANKING_USED_FOR_ENRICHMENT).upper(),
            "ranked_genes_path": path_relative_to_root(ranked_final_path, pipeline_root) if ranked_final_path else "",
            "rank_source": rank_source if ranked_final_path else "",
            "ranking_role": RANKING_PROVENANCE_ROLE if ranked_final_path else "",
            "minimum_query_gene_count": args.minimum_query_gene_count,
            "minimum_base_mean": minimum_base_mean,
            "query_execution_status": status,
            "query_exclusion_reason": reason,
        }
        manifest_rows.append(row)
        add_integrity("actual_query_gene_count_equals_manifest", "PASS", query_id, str(gene_count))
        add_integrity("query_genes_subset_of_background", "PASS", query_id, background.background_id)
        if status == "runnable":
            add_integrity("runnable_query_gene_count_at_least_minimum", "PASS", query_id, str(gene_count))
        else:
            add_integrity("query_exclusion_status_matches_minimum_gene_count", "PASS", query_id, status)

    if "complete_pan_cancer_feature_panel" in enabled_families:
        add_query(
            "complete_pan_cancer_selected_feature_panel",
            "complete_pan_cancer_feature_panel",
            feature_genes,
            pan_background,
            cancer_type="pan_cancer",
            ranks=rank_rows_from_features(features, feature_genes),
            rank_source="selection_rank" if "selection_rank" in features.columns else "",
        )

    if "cancer_type_attributed_feature_set" in enabled_families:
        for attributed_cancer_type, group in features.groupby("attributed_cancer_type", sort=True):
            add_query(
                f"cancer_type_attributed_feature_set__{attributed_cancer_type}",
                "cancer_type_attributed_feature_set",
                group["gene_id"],
                pan_background,
                cancer_type="pan_cancer",
                attributed_cancer_type=attributed_cancer_type,
                ranks=rank_rows_from_features(features, group["gene_id"]),
                rank_source="selection_rank" if "selection_rank" in features.columns else "",
            )

    if "marker_evidence_stratum" in enabled_families:
        stratum_to_genes: dict[str, list[str]] = defaultdict(list)
        for row in features.itertuples(index=False):
            for stratum in split_multivalue(row.marker_evidence_stratum, "marker_evidence_stratum", row.gene_id):
                stratum_to_genes[stratum].append(row.gene_id)
        for stratum, genes in sorted(stratum_to_genes.items()):
            add_query(
                f"marker_evidence_stratum__{stratum}",
                "marker_evidence_stratum",
                genes,
                pan_background,
                cancer_type="pan_cancer",
                marker_evidence_stratum=stratum,
                marker_evidence_source_class=stratum,
                ranks=rank_rows_from_features(features, genes),
                rank_source="selection_rank" if "selection_rank" in features.columns else "",
            )

    if "evidence_class" in enabled_families:
        evidence_class_to_genes: dict[str, list[str]] = defaultdict(list)
        for row in features.itertuples(index=False):
            for evidence_class in split_multivalue(row.evidence_classes, "evidence_classes", row.gene_id):
                evidence_class_to_genes[evidence_class].append(row.gene_id)
        for evidence_class, genes in sorted(evidence_class_to_genes.items()):
            add_query(
                f"evidence_class__{evidence_class}",
                "evidence_class",
                genes,
                pan_background,
                cancer_type="pan_cancer",
                evidence_class=evidence_class,
                ranks=rank_rows_from_features(features, genes),
                rank_source="selection_rank" if "selection_rank" in features.columns else "",
            )

    if "direction_consistency_class" in enabled_families:
        for direction_consistency_class, group in features.groupby("direction_consistency_class", sort=True):
            direction_values = sorted(set(group["direction"]))
            direction = direction_values[0] if len(direction_values) == 1 else ""
            add_query(
                f"direction_consistency_class__{direction_consistency_class}",
                "direction_consistency_class",
                group["gene_id"],
                pan_background,
                cancer_type="pan_cancer",
                direction=direction,
                direction_consistency_class=direction_consistency_class,
                ranks=rank_rows_from_features(features, group["gene_id"]),
                rank_source="selection_rank" if "selection_rank" in features.columns else "",
            )

    contrast_backgrounds: dict[tuple[str, str], BackgroundRecord] = {}
    if "contrast_level_marker_set" in enabled_families:
        for contrast in contrast_records:
            background = register_background(
                f"contrast_background__{contrast.cancer_type}__{contrast.contrast_id}",
                "contrast_specific_deseq2_tested_gene_background",
                (
                    "Genes with non-missing adjusted_p_value and baseMean >= "
                    f"{contrast.minimum_base_mean} in the canonical post-prefilter contrast table."
                ),
                contrast.background_genes,
            )
            contrast_backgrounds[(contrast.cancer_type, contrast.contrast_id)] = background
            add_query(
                (
                    "contrast_level_marker_set__"
                    f"{contrast.cancer_type}__{contrast.marker_evidence_stratum}__{contrast.contrast_id}"
                ),
                "contrast_level_marker_set",
                contrast.marker_genes,
                background,
                cancer_type=contrast.cancer_type,
                marker_evidence_stratum=contrast.marker_evidence_stratum,
                marker_evidence_source_class=contrast.marker_evidence_stratum,
                contrast_id=contrast.contrast_id,
                contrast_type=contrast.contrast_type,
                ranks=contrast.ranks,
                rank_source="absolute_shrunken_log2_fold_change",
                minimum_base_mean=contrast.minimum_base_mean,
            )

    if "enrichment_contrast_marker_recurrence" in enabled_families:
        grouped: dict[tuple[str, str], list[ContrastRecord]] = defaultdict(list)
        for contrast in contrast_records:
            grouped[(contrast.cancer_type, contrast.marker_evidence_stratum)].append(contrast)
        for (cancer_type, stratum), contrasts in sorted(grouped.items()):
            counts: Counter[str] = Counter()
            background_genes: list[str] = []
            rank_rows: list[pd.DataFrame] = []
            for contrast in contrasts:
                counts.update(set(contrast.marker_genes))
                background_genes.extend(contrast.background_genes)
                if not contrast.ranks.empty:
                    rank_rows.append(contrast.ranks)
            recurrent_genes = [
                gene for gene, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
                if count >= args.recurrence_minimum_retained_contrast_count
            ]
            for gene, count in sorted(counts.items(), key=lambda item: (item[0])):
                recurrence_rows.append(
                    {
                        "cancer_type": cancer_type,
                        "marker_evidence_stratum": stratum,
                        "gene_id": gene,
                        "retained_contrast_marker_set_count": count,
                        "retained_contrast_marker_set_total": len(contrasts),
                        "recurrence_minimum_retained_contrast_count": args.recurrence_minimum_retained_contrast_count,
                    }
                )
            background = register_background(
                f"enrichment_contrast_marker_recurrence_background__{cancer_type}__{stratum}",
                "enrichment_contrast_group_deseq2_tested_gene_background_union",
                (
                    "Union of contrast-specific backgrounds for cancer_type x "
                    "marker_evidence_stratum enrichment recurrence grouping."
                ),
                background_genes,
            )
            ranks = None
            if rank_rows:
                ranks = pd.concat(rank_rows, ignore_index=True)
                ranks = ranks[ranks["gene_id"].isin(recurrent_genes)]
                ranks = (
                    ranks.groupby("gene_id", as_index=False)["rank_stat"]
                    .median()
                    .sort_values(["rank_stat", "gene_id"], ascending=[False, True])
                )
            add_query(
                f"enrichment_contrast_marker_recurrence__{cancer_type}__{stratum}",
                "enrichment_contrast_marker_recurrence",
                recurrent_genes,
                background,
                cancer_type=cancer_type,
                marker_evidence_stratum=stratum,
                marker_evidence_source_class=stratum,
                ranks=ranks,
                rank_source="absolute_shrunken_log2_fold_change" if ranks is not None and not ranks.empty else "",
            )

    if not manifest_rows:
        raise ValueError("No enrichment query rows were constructed")

    manifest = pd.DataFrame(manifest_rows, columns=MANIFEST_COLUMNS)
    if manifest["query_id"].duplicated().any():
        duplicates = sorted(manifest.loc[manifest["query_id"].duplicated(), "query_id"].unique())
        raise ValueError("Duplicate query_id values in enrichment manifest: " + ",".join(duplicates))
    if manifest["query_output_directory_id"].duplicated().any():
        duplicates = sorted(
            manifest.loc[manifest["query_output_directory_id"].duplicated(), "query_output_directory_id"].unique()
        )
        raise ValueError(
            "Duplicate normalized query output directory identifiers: "
            + ",".join(duplicates)
        )

    complete = manifest[manifest["query_family"] == "complete_pan_cancer_feature_panel"]
    if len(complete) == 1:
        complete_final_path = pipeline_root / str(complete.iloc[0]["gene_list_path"])
        complete_query_genes = read_gene_file(
            stage_path_for(complete_final_path),
            "complete_pan_cancer_feature_panel",
        )
        missing_feature_from_query = sorted(set(feature_genes) - set(complete_query_genes))
        missing_query_from_feature = sorted(set(complete_query_genes) - set(feature_genes))
        count_equal = len(feature_genes) == len(selected_feature_genes) == len(complete_query_genes)
        set_equal = not missing_feature_from_query and not missing_query_from_feature
        if not count_equal or not set_equal:
            raise ValueError(
                "Complete-panel gene identity validation failed: "
                f"feature_table_count={len(feature_genes)} "
                f"selected_feature_gene_file_count={len(selected_feature_genes)} "
                f"enrichment_query_count={len(complete_query_genes)} "
                f"feature_genes_missing_from_query={','.join(missing_feature_from_query[:100])} "
                f"query_genes_missing_from_feature_table={','.join(missing_query_from_feature[:100])}"
            )
        add_integrity("complete_panel_gene_identity", "PASS", complete.iloc[0]["query_id"], str(len(feature_genes)))
    elif "complete_pan_cancer_feature_panel" in enabled_families:
        raise ValueError("Expected exactly one complete_pan_cancer_feature_panel query")

    runnable = manifest[manifest["query_execution_status"] == "runnable"]
    if (runnable["gene_count"].astype(int) < args.minimum_query_gene_count).any():
        raise ValueError("A runnable query has fewer genes than minimum_query_gene_count")
    if (runnable["gene_count"].astype(int) == 0).any():
        raise ValueError("A zero-gene query is marked runnable")
    if (runnable["ordered"].astype(str).str.upper() != "FALSE").any():
        raise ValueError("Primary enrichment query construction produced an ordered runnable query")
    if (runnable["ranking_used_for_enrichment"].astype(str).str.upper() != "FALSE").any():
        raise ValueError("Primary enrichment query construction marked ranking_used_for_enrichment TRUE")

    for background_id, group in manifest.groupby("background_id", sort=True):
        observed_hashes = set(group["background_sha256"])
        if len(observed_hashes) != 1:
            raise ValueError(f"background_id={background_id} has conflicting SHA256 values")
        add_integrity("repeated_background_id_has_identical_sha256", "PASS", "", background_id)

    manifest.to_csv(staging / "enrichment_query_manifest.tsv", sep="\t", index=False)
    manifest[manifest["query_execution_status"] != "runnable"].to_csv(
        staging / "skipped_queries.tsv",
        sep="\t",
        index=False,
    )
    summary = (
        manifest.groupby(["query_family", "query_execution_status"], dropna=False)
        .size()
        .reset_index(name="n_queries")
        .sort_values(["query_family", "query_execution_status"])
    )
    summary.to_csv(staging / "query_summary.tsv", sep="\t", index=False)
    pd.DataFrame(integrity_rows).to_csv(
        staging / "enrichment_query_background_integrity.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(recurrence_rows).to_csv(
        staging / "enrichment_contrast_marker_recurrence.tsv",
        sep="\t",
        index=False,
    )
    (staging / ".done").write_text("done\n")

    publish_staging(staging, outdir)
    print("[Enrichment DAG] Canonical graph-derived enrichment query set published.")
    print(
        "[Pan-cancer background] Feature-selection-eligible gene background consumed as "
        "explicit analytical background."
    )
    print(
        "[Primary enrichment] Primary analysis manifest uses unordered ORA with "
        "custom analytical backgrounds."
    )
    print(
        f"[Query validation] {len(manifest)} queries; "
        f"{len(runnable)} runnable; {len(manifest) - len(runnable)} excluded below interface gates."
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
