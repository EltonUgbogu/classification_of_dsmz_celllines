#!/usr/bin/env python3
"""Build pan-cancer marker-framework enrichment query/background files.

This helper prepares local inputs for manual g:Profiler upload or parser and
plotter testing. It uses the current graph-derived pan-cancer feature classes:
recurrent, accepted_singleton, and accepted_non_recurrent. It does not call
g:Profiler and does not rerun DESeq2 or feature construction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

COHORTS = ("BRCA", "NBL", "RBL")
DIRECTIONS = ("UP", "DOWN", "MIXED")
FEATURE_CLASS_QUERY_ORDER = (
    "RECURRENT",
    "ACCEPTED_SINGLETON",
    "ACCEPTED_NON_RECURRENT",
)
STALE_INPUT_MARKERS = (
    "stale_old_gprofiler_query_sets_20260610_173456/quarantine/gene_sets",
    "results/unsupervised/rbl/enrichment/query_sets/gene_sets",
    "_stale_query_sets_legacy_labels_",
)
FEATURE_REQUIRED_COLUMNS = (
    "gene_id",
    "direction",
    "feature_class",
    "owner_profile",
    "source_contrast",
    "source_file",
    "selection_rank",
)
MANIFEST_COLUMNS = (
    "query_name",
    "query_id",
    "category",
    "owner_profile",
    "direction",
    "source_contrast",
    "gene_count",
    "background_count",
    "genes_path",
    "background_path",
    "ranked_genes_path",
    "metadata_path",
    "custom_background_available",
    "background_strategy",
    "ranked_genes_available",
    "rank_source",
    "source_table",
    "source_table_sha256",
    "observed_feature_count",
    "notes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate pan-cancer marker-framework enrichment query sets."
    )
    parser.add_argument("--features", required=True, help="Current pan_cancer_features.tsv")
    parser.add_argument("--clean-features", required=True, help="Current pan_cancer_features_clean.txt")
    parser.add_argument(
        "--retained-source-contrasts",
        default="",
        help="Optional retained_gene_source_contrasts.tsv provenance cross-check table",
    )
    parser.add_argument(
        "--cohort-query-dir",
        required=True,
        help="Current Snakefile-managed gprofiler_by_cohort directory",
    )
    parser.add_argument("--outdir", required=True, help="Output marker-framework directory")
    parser.add_argument(
        "--profiles",
        default="brca,nbl,rbl",
        help="Comma-separated disease profiles to inspect for current marker manifests",
    )
    parser.add_argument(
        "--background-min-normalised-count",
        type=float,
        default=10.0,
        help="Minimum test-sample normalised count for DESeq2-tested background genes",
    )
    return parser.parse_args()


def fail(errors: list[str]) -> int:
    print("ERROR: marker-framework query-set build failed", file=sys.stderr)
    for err in errors:
        print(f"- {err}", file=sys.stderr)
    return 1


def clean_gene(value: object) -> str:
    return re.sub(r"\.[0-9]+$", "", str(value or "").strip())


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def parse_float(value: object) -> float | None:
    text = str(value or "").strip()
    if text == "" or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = reader.fieldnames or []
    return columns, rows


def write_tsv(path: Path, columns: tuple[str, ...] | list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def unique_ordered(values) -> list[str]:
    seen = set()
    out = []
    for value in values:
        gene = clean_gene(value)
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return out


def read_gene_list(path: Path) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    genes = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        first = True
        for line in handle:
            text = line.rstrip("\n")
            if not text.strip():
                continue
            field = text.split("\t", 1)[0].strip()
            if first and field.lower() in {"gene", "gene_id", "ensembl_gene_id"}:
                first = False
                continue
            first = False
            genes.append(field)
    return unique_ordered(genes)


def resolve_manifest_path(base_dir: Path, value: str) -> Path:
    path = Path(str(value or "").strip())
    if path.is_absolute():
        return path
    return base_dir / path


def check_stale_paths(label: str, value: str, errors: list[str]) -> None:
    norm = str(value or "")
    for marker in STALE_INPUT_MARKERS:
        if marker in norm:
            errors.append(f"{label} points at stale/quarantined enrichment tree: {norm}")


def rank_key(row: dict[str, str]) -> tuple[float, str]:
    rank = parse_float(row.get("selection_rank"))
    if rank is None:
        rank = float("inf")
    return (rank, clean_gene(row.get("gene_id")))


def feature_rank_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out = []
    for row in rows:
        gene = clean_gene(row.get("gene_id"))
        if not gene:
            continue
        rank = parse_float(row.get("max_abs_log2FC"))
        source = "max_abs_log2FC"
        if rank is None:
            selection = parse_float(row.get("selection_rank"))
            rank = 1.0 / (selection + 1.0) if selection is not None else None
            source = "inverse_selection_rank" if rank is not None else ""
        if rank is not None:
            out.append({"gene_id": gene, "rank_stat": rank, "rank_source": source})
    out.sort(key=lambda x: (-float(x["rank_stat"]), x["gene_id"]))
    return out


def read_deseq_table(path: Path, min_normalised_count: float) -> tuple[dict[str, dict[str, float | None]], list[str]]:
    if not path.exists() or path.stat().st_size == 0:
        return {}, []
    columns, rows = read_tsv(path)
    if not columns:
        return {}, []
    gene_col = "gene_id" if "gene_id" in columns else columns[0]
    norm_col = ""
    for candidate in (
        "normalised_count_in_test_sample",
        "normalized_count_in_test_sample",
        "test_sample_normalised_count",
        "test_sample_normalized_count",
    ):
        if candidate in columns:
            norm_col = candidate
            break
    info = {}
    background = []
    for row in rows:
        gene = clean_gene(row.get(gene_col))
        if not gene:
            continue
        lfc = parse_float(row.get("log2FoldChange"))
        padj = parse_float(row.get("padj"))
        norm = parse_float(row.get(norm_col)) if norm_col else None
        info[gene] = {"log2FoldChange": lfc, "padj": padj, "normalised_count_in_test_sample": norm}
        if padj is not None and norm is not None and norm >= min_normalised_count:
            background.append(gene)
    return info, unique_ordered(background)


def marker_rank_rows(genes: list[str], table_info: dict[str, dict[str, float | None]]) -> list[dict[str, object]]:
    out = []
    for gene in genes:
        lfc = table_info.get(gene, {}).get("log2FoldChange")
        if lfc is not None:
            out.append({"gene_id": gene, "rank_stat": abs(float(lfc)), "rank_source": "abs_log2FoldChange"})
    out.sort(key=lambda x: (-float(x["rank_stat"]), x["gene_id"]))
    return out


def read_marker_manifests(profiles: list[str], min_normalised_count: float):
    contrasts = []
    profile_backgrounds = defaultdict(list)
    validation_rows = []
    availability_notes = []

    for profile in profiles:
        p = profile.lower()
        manifest_path = Path("results") / "unsupervised" / p / "deseq2_markers" / "markers" / "marker_sets_manifest.tsv"
        if not manifest_path.exists():
            availability_notes.append(f"{p}: marker_sets_manifest.tsv missing")
            validation_rows.append({
                "scope": "marker_manifest",
                "item": str(manifest_path),
                "status": "missing",
                "observed": "0",
                "expected": "present",
                "notes": "per-contrast query sets unavailable for this profile",
            })
            continue
        base_dir = manifest_path.parent.parent
        columns, rows = read_tsv(manifest_path)
        required = {"contrast", "marker_file", "table_file"}
        missing = sorted(required.difference(columns))
        if missing:
            availability_notes.append(f"{p}: manifest missing columns " + ",".join(missing))
            validation_rows.append({
                "scope": "marker_manifest",
                "item": str(manifest_path),
                "status": "failed",
                "observed": ",".join(columns),
                "expected": ",".join(sorted(required)),
                "notes": "per-contrast query sets unavailable for this profile",
            })
            continue
        for row in rows:
            contrast = str(row.get("contrast") or "").strip()
            marker_path = resolve_manifest_path(base_dir, row.get("marker_file", ""))
            table_path = resolve_manifest_path(base_dir, row.get("table_file", ""))
            genes = read_gene_list(marker_path)
            table_info, background = read_deseq_table(table_path, min_normalised_count)
            profile_backgrounds[p.upper()].extend(background)
            up = [gene for gene in genes if (table_info.get(gene, {}).get("log2FoldChange") or 0) > 0]
            down = [gene for gene in genes if (table_info.get(gene, {}).get("log2FoldChange") or 0) < 0]
            contrasts.append({
                "profile": p,
                "owner_profile": p.upper(),
                "contrast": contrast,
                "marker_path": marker_path,
                "table_path": table_path,
                "genes": genes,
                "up": unique_ordered(up),
                "down": unique_ordered(down),
                "background": background,
                "table_info": table_info,
            })
            status = "pass" if genes and background else "failed"
            validation_rows.append({
                "scope": "per_contrast_background",
                "item": f"{p}:{contrast}",
                "status": status,
                "observed": str(len(background)),
                "expected": "non_empty",
                "notes": f"marker_genes={len(genes)}; marker_file={marker_path}; table_file={table_path}",
            })
    profile_backgrounds = {profile: unique_ordered(genes) for profile, genes in profile_backgrounds.items()}
    return contrasts, profile_backgrounds, validation_rows, availability_notes


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    for label, value in (
        ("features", args.features),
        ("clean_features", args.clean_features),
        ("retained_source_contrasts", args.retained_source_contrasts),
        ("cohort_query_dir", args.cohort_query_dir),
        ("outdir", args.outdir),
    ):
        check_stale_paths(label, value, errors)

    feature_path = Path(args.features)
    clean_feature_path = Path(args.clean_features)
    outdir = Path(args.outdir)
    gene_root = outdir / "gene_sets"
    if not feature_path.exists():
        errors.append(f"feature table missing: {feature_path}")
    if not clean_feature_path.exists():
        errors.append(f"clean feature list missing: {clean_feature_path}")
    if errors:
        return fail(errors)

    feature_columns, feature_rows_raw = read_tsv(feature_path)
    missing_feature_cols = sorted(set(FEATURE_REQUIRED_COLUMNS).difference(feature_columns))
    if missing_feature_cols:
        errors.append("feature table missing required columns: " + ",".join(missing_feature_cols))
    if errors:
        return fail(errors)

    feature_rows = []
    for row in feature_rows_raw:
        normalised = dict(row)
        normalised["gene_id"] = clean_gene(row.get("gene_id"))
        normalised["direction"] = str(row.get("direction") or "").strip().upper()
        normalised["feature_class"] = str(row.get("feature_class") or "").strip().upper()
        normalised["owner_profile"] = str(row.get("owner_profile") or "").strip().upper()
        feature_rows.append(normalised)
    feature_rows.sort(key=rank_key)

    source_table_rows = len(feature_rows)
    feature_genes = unique_ordered(row["gene_id"] for row in feature_rows)
    observed_feature_count = len(feature_genes)
    if source_table_rows != observed_feature_count:
        errors.append(
            f"feature table has {source_table_rows} rows but {observed_feature_count} unique clean gene IDs"
        )
    if any(not row["gene_id"] for row in feature_rows):
        errors.append("feature table contains blank gene IDs")
    bad_owners = sorted({row["owner_profile"] for row in feature_rows}.difference(COHORTS))
    if bad_owners:
        errors.append("unexpected owner_profile values: " + ",".join(bad_owners))
    bad_directions = sorted({row["direction"] for row in feature_rows}.difference(DIRECTIONS))
    if bad_directions:
        errors.append("unexpected direction values: " + ",".join(bad_directions))

    clean_genes = read_gene_list(clean_feature_path)
    validation_rows: list[dict[str, object]] = []
    validation_rows.append({
        "scope": "feature_table",
        "item": str(feature_path),
        "status": "pass" if source_table_rows == observed_feature_count else "failed",
        "observed": str(observed_feature_count),
        "expected": str(source_table_rows),
        "notes": "row count and unique clean gene count",
    })
    validation_rows.append({
        "scope": "clean_feature_list",
        "item": str(clean_feature_path),
        "status": "pass" if clean_genes == feature_genes else "failed",
        "observed": str(len(clean_genes)),
        "expected": str(observed_feature_count),
        "notes": "clean feature list matches feature table order" if clean_genes == feature_genes else "clean feature list differs from feature table",
    })
    if clean_genes != feature_genes:
        errors.append("clean feature list does not match current feature table")

    retained_path = Path(args.retained_source_contrasts) if args.retained_source_contrasts else None
    if retained_path is not None and retained_path.exists():
        retained_columns, retained_rows_raw = read_tsv(retained_path)
        retained_required = {"gene_id", "direction", "feature_class", "owner_profile", "source_contrast", "source_file", "selection_rank"}
        retained_missing = sorted(retained_required.difference(retained_columns))
        if retained_missing:
            validation_rows.append({
                "scope": "retained_source_contrasts",
                "item": str(retained_path),
                "status": "failed",
                "observed": ",".join(retained_columns),
                "expected": ",".join(sorted(retained_required)),
                "notes": "retained provenance table missing required columns",
            })
            errors.append("retained source contrast table missing required columns")
        else:
            def retained_key(row):
                return (
                    clean_gene(row.get("gene_id")),
                    str(row.get("direction") or "").strip().upper(),
                    str(row.get("feature_class") or "").strip().upper(),
                    str(row.get("owner_profile") or "").strip().upper(),
                    str(row.get("source_contrast") or "").strip(),
                    str(row.get("source_file") or "").strip(),
                    str(row.get("selection_rank") or "").strip(),
                )
            feature_keys = {retained_key(row) for row in feature_rows}
            retained_keys = {retained_key(row) for row in retained_rows_raw}
            status = "pass" if feature_keys == retained_keys and len(retained_rows_raw) == source_table_rows else "failed"
            validation_rows.append({
                "scope": "retained_source_contrasts",
                "item": str(retained_path),
                "status": status,
                "observed": str(len(retained_rows_raw)),
                "expected": str(source_table_rows),
                "notes": "retained provenance keys match feature table" if status == "pass" else "retained provenance keys differ from feature table",
            })
            if status != "pass":
                errors.append("retained source contrast table does not match current feature table")
    elif args.retained_source_contrasts:
        validation_rows.append({
            "scope": "retained_source_contrasts",
            "item": args.retained_source_contrasts,
            "status": "missing",
            "observed": "0",
            "expected": "present_if_configured",
            "notes": "optional retained provenance cross-check skipped because file is missing",
        })

    feature_sha = sha256_file(feature_path)
    profiles = [p.strip().lower() for p in args.profiles.split(",") if p.strip()]
    contrasts, profile_backgrounds, contrast_validation, marker_notes = read_marker_manifests(
        profiles, args.background_min_normalised_count
    )
    validation_rows.extend(contrast_validation)
    marker_universe = unique_ordered(gene for genes in profile_backgrounds.values() for gene in genes)

    outdir.mkdir(parents=True, exist_ok=True)
    gene_root.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    category_rows: list[dict[str, object]] = []

    def add_query(
        query_name: str,
        category: str,
        owner_profile: str,
        direction: str,
        source_contrast: str,
        genes: list[str],
        background: list[str],
        rank_rows: list[dict[str, object]],
        rank_source: str,
        source_table: Path,
        background_strategy: str,
        notes: str,
    ) -> None:
        query_name_safe = safe_name(query_name)
        qdir = gene_root / query_name_safe
        qdir.mkdir(parents=True, exist_ok=True)
        genes = unique_ordered(genes)
        background = unique_ordered(background)
        genes_path = qdir / "genes.tsv"
        background_path = qdir / "background.tsv"
        ranked_path = qdir / "ranked_genes.tsv"
        metadata_path = qdir / "metadata.tsv"
        write_tsv(genes_path, ("gene_id",), [{"gene_id": gene} for gene in genes])
        write_tsv(background_path, ("gene_id",), [{"gene_id": gene} for gene in background])
        rank_rows = [row for row in rank_rows if clean_gene(row.get("gene_id")) in set(genes)]
        write_tsv(ranked_path, ("gene_id", "rank_stat", "rank_source"), rank_rows)
        row = {
            "query_name": query_name_safe,
            "query_id": query_name_safe,
            "category": category,
            "owner_profile": owner_profile,
            "direction": direction,
            "source_contrast": source_contrast,
            "gene_count": len(genes),
            "background_count": len(background),
            "genes_path": str(genes_path),
            "background_path": str(background_path),
            "ranked_genes_path": str(ranked_path),
            "metadata_path": str(metadata_path),
            "custom_background_available": "TRUE" if background else "FALSE",
            "background_strategy": background_strategy,
            "ranked_genes_available": "TRUE" if rank_rows else "FALSE",
            "rank_source": rank_source if rank_rows else "",
            "source_table": str(source_table),
            "source_table_sha256": sha256_file(source_table) if source_table.exists() else "",
            "observed_feature_count": observed_feature_count,
            "notes": notes,
        }
        write_tsv(metadata_path, MANIFEST_COLUMNS, [row])
        manifest_rows.append(row)

    def rows_for(owner: str | None = None, direction: str | None = None, classes: set[str] | None = None):
        out = feature_rows
        if owner is not None:
            out = [row for row in out if row["owner_profile"] == owner]
        if direction is not None:
            out = [row for row in out if row["direction"] == direction]
        if classes is not None:
            out = [row for row in out if row["feature_class"] in classes]
        return sorted(out, key=rank_key)

    def feature_genes_from(rows: list[dict[str, str]]) -> list[str]:
        return unique_ordered(row["gene_id"] for row in rows)

    # A. Cohort/direction sets, matched against current Snakefile-managed outputs.
    cohort_query_dir = Path(args.cohort_query_dir)
    for cohort in COHORTS:
        cohort_background = profile_backgrounds.get(cohort, [])
        for direction in ("all",) + DIRECTIONS:
            selected_rows = rows_for(owner=cohort, direction=None if direction == "all" else direction)
            query_name = f"{cohort}_{direction}"
            genes = feature_genes_from(selected_rows)
            add_query(
                query_name=query_name,
                category="cohort_direction",
                owner_profile=cohort,
                direction=direction,
                source_contrast="",
                genes=genes,
                background=cohort_background,
                rank_rows=feature_rank_rows(selected_rows),
                rank_source="max_abs_log2FC_or_inverse_selection_rank",
                source_table=feature_path,
                background_strategy="cohort-specific union of current DESeq2-tested marker-selection backgrounds",
                notes="Current cohort/direction feature subset from pan_cancer_features.tsv",
            )
            current_path = cohort_query_dir / f"{query_name}.txt"
            if current_path.exists():
                current_genes = read_gene_list(current_path)
                status = "pass" if current_genes == genes else "failed"
                validation_rows.append({
                    "scope": "cohort_query_match",
                    "item": query_name,
                    "status": status,
                    "observed": str(len(genes)),
                    "expected": str(len(current_genes)),
                    "notes": f"compared with {current_path}",
                })
                if status != "pass":
                    errors.append(f"{query_name} does not match current gprofiler_by_cohort output")
            else:
                validation_rows.append({
                    "scope": "cohort_query_match",
                    "item": query_name,
                    "status": "missing",
                    "observed": str(len(genes)),
                    "expected": "current cohort file present",
                    "notes": f"missing {current_path}",
                })

    # B. Per-contrast marker sets from current marker manifests.
    for item in contrasts:
        base = item["profile"] + "__per_contrast__" + safe_name(item["contrast"])
        for suffix, genes in (("all", item["genes"]), ("up", item["up"]), ("down", item["down"])):
            add_query(
                query_name=f"{base}__{suffix}",
                category="per_contrast",
                owner_profile=item["owner_profile"],
                direction=suffix,
                source_contrast=item["contrast"],
                genes=genes,
                background=item["background"],
                rank_rows=marker_rank_rows(genes, item["table_info"]),
                rank_source="abs_log2FoldChange",
                source_table=item["marker_path"],
                background_strategy="contrast-specific DESeq2-tested genes with non-missing padj and sufficient test-sample expression",
                notes="table_file=" + str(item["table_path"]),
            )

    # C. Revised feature-class query sets.
    for feature_class in FEATURE_CLASS_QUERY_ORDER:
        class_rows_all = rows_for(classes={feature_class})
        if class_rows_all:
            label = feature_class.lower()
            for suffix, selected_rows in (
                ("all", class_rows_all),
                ("up", rows_for(direction="UP", classes={feature_class})),
                ("down", rows_for(direction="DOWN", classes={feature_class})),
                ("mixed", rows_for(direction="MIXED", classes={feature_class})),
            ):
                add_query(
                    query_name=f"feature_class_{label}_{suffix}",
                    category="feature_class",
                    owner_profile="pan_cancer",
                    direction=suffix,
                    source_contrast="",
                    genes=feature_genes_from(selected_rows),
                    background=marker_universe,
                    rank_rows=feature_rank_rows(selected_rows),
                    rank_source="max_abs_log2FC_or_inverse_selection_rank",
                    source_table=feature_path,
                    background_strategy="combined union of current cohort marker-selection backgrounds",
                    notes=f"Current feature_class subset: {label}",
                )

    # D. Current final pan-cancer feature set.
    for suffix, selected_rows in (
        ("all", feature_rows),
        ("up", rows_for(direction="UP")),
        ("down", rows_for(direction="DOWN")),
        ("mixed", rows_for(direction="MIXED")),
    ):
        add_query(
            query_name=f"final_pan_cancer_feature_set_{suffix}",
            category="final_pan_cancer_feature_set",
            owner_profile="pan_cancer",
            direction=suffix,
            source_contrast="",
            genes=feature_genes_from(selected_rows),
            background=marker_universe,
            rank_rows=feature_rank_rows(selected_rows),
            rank_source="max_abs_log2FC_or_inverse_selection_rank",
            source_table=feature_path,
            background_strategy="combined union of current cohort marker-selection backgrounds",
            notes="Current full feature table subset; count inferred dynamically",
        )

    # Validate final feature query contains the complete current feature set.
    final_rows = [row for row in manifest_rows if row["query_name"] == "final_pan_cancer_feature_set_all"]
    if final_rows:
        final_genes = read_gene_list(Path(final_rows[0]["genes_path"]))
        status = "pass" if set(final_genes) == set(feature_genes) and len(final_genes) == observed_feature_count else "failed"
        validation_rows.append({
            "scope": "final_feature_set_query",
            "item": "final_pan_cancer_feature_set_all",
            "status": status,
            "observed": str(len(final_genes)),
            "expected": str(observed_feature_count),
            "notes": "query contains all current feature-table genes",
        })
        if status != "pass":
            errors.append("final_pan_cancer_feature_set_all does not contain the full current feature set")

    categories = {
        "cohort_direction": "Current feature-table owner_profile/direction subsets",
        "per_contrast": "Current marker manifests under results/unsupervised/<profile>/deseq2_markers/markers",
        "feature_class": "Current revised feature_class subsets",
        "final_pan_cancer_feature_set": "Current full feature set",
    }
    for category, reason in categories.items():
        rows = [row for row in manifest_rows if row["category"] == category]
        status = "available" if rows else "unavailable"
        details = reason
        if category == "per_contrast" and marker_notes:
            details += "; " + "; ".join(marker_notes)
        if category == "feature_class" and not rows:
            details += "; no current feature_class values matched revised classes"
        category_rows.append({
            "category": category,
            "status": status,
            "generated_query_count": len(rows),
            "total_gene_count_across_queries": sum(int(row["gene_count"]) for row in rows),
            "notes": details,
        })

    manifest_rows.sort(key=lambda r: (str(r["category"]), str(r["owner_profile"]), str(r["source_contrast"]), str(r["direction"]), str(r["query_name"])))
    write_tsv(outdir / "query_manifest.tsv", MANIFEST_COLUMNS, manifest_rows)
    write_tsv(
        outdir / "query_set_counts.tsv",
        (
            "query_name",
            "category",
            "owner_profile",
            "direction",
            "source_contrast",
            "gene_count",
            "background_count",
            "custom_background_available",
            "observed_feature_count",
            "notes",
        ),
        manifest_rows,
    )
    write_tsv(
        outdir / "background_validation.tsv",
        (
            "query_name",
            "category",
            "owner_profile",
            "direction",
            "background_count",
            "custom_background_available",
            "background_strategy",
            "status",
            "notes",
        ),
        [
            {
                **row,
                "status": "pass" if int(row["background_count"]) > 0 else "missing_background",
            }
            for row in manifest_rows
        ],
    )
    write_tsv(
        outdir / "category_availability.tsv",
        ("category", "status", "generated_query_count", "total_gene_count_across_queries", "notes"),
        category_rows,
    )
    write_tsv(
        outdir / "validation_summary.tsv",
        ("scope", "item", "status", "observed", "expected", "notes"),
        validation_rows,
    )
    write_tsv(
        outdir / "source_feature_summary.tsv",
        ("metric", "value"),
        [
            {"metric": "feature_table", "value": str(feature_path)},
            {"metric": "feature_table_sha256", "value": feature_sha},
            {"metric": "feature_table_rows", "value": source_table_rows},
            {"metric": "feature_table_unique_genes", "value": observed_feature_count},
            {"metric": "clean_feature_list", "value": str(clean_feature_path)},
            {"metric": "marker_selection_universe_count", "value": len(marker_universe)},
            {"metric": "retained_source_contrasts", "value": args.retained_source_contrasts},
        ],
    )
    readme = outdir / "README_marker_framework_queries.md"
    readme.write_text(
        "# Pan-cancer marker-framework enrichment query sets\n\n"
        "Generated local query/background inputs for manual g:Profiler use or parser/plotter testing. "
        "Feature-class query sets use the revised classes `recurrent`, `accepted_singleton`, and `accepted_non_recurrent`. "
        "This build does not run g:Profiler, DESeq2, or pan-cancer feature construction.\n\n"
        f"Observed current feature count: {observed_feature_count}\n\n"
        "Each `gene_sets/<query_name>/` directory contains `genes.tsv`, `background.tsv`, "
        "`ranked_genes.tsv`, and `metadata.tsv`. Backgrounds are custom backgrounds when available. "
        "Per-contrast backgrounds use the matching current DESeq2 table; pan-cancer feature-class "
        "queries use the union of current cohort marker-selection backgrounds.\n",
        encoding="utf-8",
    )

    failed_validation = [row for row in validation_rows if row.get("status") == "failed"]
    if failed_validation:
        errors.extend("validation failed: " + str(row.get("scope", "")) + " " + str(row.get("item", "")) for row in failed_validation)
    if errors:
        return fail(errors)

    print(f"features={feature_path}")
    print(f"source_table_sha256={feature_sha}")
    print(f"observed_feature_count={observed_feature_count}")
    print("query_manifest=" + str(outdir / "query_manifest.tsv"))
    print(f"query_count={len(manifest_rows)}")
    print("gprofiler_run=not_run")
    print("deseq2_run=not_run")
    print("pan_cancer_feature_construction=not_run")
    print("validation_status=pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
