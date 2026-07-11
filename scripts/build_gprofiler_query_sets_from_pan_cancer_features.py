#!/usr/bin/env python3
"""Build cohort-specific g:Profiler query lists from pan-cancer features.

This helper prepares cohort-level query-list inputs and validation manifests
without running g:Profiler. Stale 171-gene references are used only as guard
patterns for archived paths, not as current feature-set descriptions.
"""

import argparse
import csv
import hashlib
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

COHORTS = ("BRCA", "NBL", "RBL")
DIRECTIONS = ("UP", "DOWN", "MIXED")
QUERY_DIRECTIONS = ("all",) + DIRECTIONS
REQUIRED_COLUMNS = (
    "gene_id",
    "direction",
    "owner_profile",
    "source_contrast",
    "source_file",
    "selection_rank",
)
STALE_171_SENTINEL = "final_pan_cancer_feature_set__all/genes.tsv"
STALE_RBL_QUERY_DIR = "results/unsupervised/rbl/enrichment/query_sets/gene_sets"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate validated cohort-specific g:Profiler query lists."
    )
    parser.add_argument("--features", required=True, help="pan_cancer_features.tsv")
    parser.add_argument("--outdir", required=True, help="Output query-list directory")
    parser.add_argument(
        "--retained_source_contrasts",
        "--retained-source-contrasts",
        default="",
        help="Optional retained_gene_source_contrasts.tsv cross-check table",
    )
    return parser.parse_args()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = reader.fieldnames or []
    return columns, rows


def read_gene_list(path):
    if os.path.getsize(path) == 0:
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return [line.rstrip("\n") for line in handle]


def normalise_row(row):
    out = dict(row)
    out["gene_id"] = (out.get("gene_id") or "").strip()
    out["direction"] = (out.get("direction") or "").strip().upper()
    out["owner_profile"] = (out.get("owner_profile") or "").strip().upper()
    out["source_contrast"] = (out.get("source_contrast") or "").strip()
    out["source_file"] = (out.get("source_file") or "").strip()
    out["selection_rank"] = (out.get("selection_rank") or "").strip()
    return out


def rank_key(row):
    raw = row.get("selection_rank", "")
    try:
        rank = float(raw)
    except ValueError:
        rank = float("inf")
    return (rank, row.get("gene_id", ""))


def query_id(cohort, direction):
    return f"{cohort}_{direction}"


def iter_query_ids():
    for cohort in COHORTS:
        for direction in QUERY_DIRECTIONS:
            yield cohort, direction, query_id(cohort, direction)


def count_query_rows(rows):
    counts = Counter()
    for row in rows:
        counts[query_id(row["owner_profile"], "all")] += 1
        counts[query_id(row["owner_profile"], row["direction"])] += 1
    return counts


def write_gene_list(path, genes):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(genes)
    if genes:
        text += "\n"
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_gmt(path, entries):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        for name, description, genes in entries:
            fields = [name, description] + list(genes)
            handle.write("\t".join(fields) + "\n")


def assert_nonempty_text_files_have_newline(paths):
    errors = []
    for path in paths:
        if os.path.getsize(path) == 0:
            continue
        with open(path, "rb") as handle:
            handle.seek(-1, os.SEEK_END)
            if handle.read(1) != b"\n":
                errors.append(f"non-empty text file lacks trailing newline: {path}")
    return errors


def key_set(rows):
    return {
        (
            row["gene_id"],
            row["direction"],
            row["owner_profile"],
            row["source_contrast"],
            row["source_file"],
            row["selection_rank"],
        )
        for row in rows
    }


def fail(errors):
    print("ERROR: g:Profiler query-set validation failed", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    return 1


def validate_output_lists(output_genes, query_genes, source_counts, unique_genes):
    errors = []

    for cohort, direction, name in iter_query_ids():
        observed = output_genes.get(name, [])
        expected = query_genes.get(name, [])
        if observed != expected:
            errors.append(f"{name}.txt does not match genes computed from the source table")
        if len(observed) != source_counts.get(name, 0):
            errors.append(
                f"{name}.txt has {len(observed)} genes; source table count is {source_counts.get(name, 0)}"
            )
        if not observed and source_counts.get(name, 0) != 0:
            errors.append(f"{name}.txt is empty but the source table count is {source_counts.get(name, 0)}")

    all_union = set()
    for cohort in COHORTS:
        cohort_all = set(output_genes[query_id(cohort, "all")])
        direction_union = set()
        for direction in DIRECTIONS:
            direction_union.update(output_genes[query_id(cohort, direction)])
        if cohort_all != direction_union:
            errors.append(f"{cohort}_all does not equal the union of {cohort}_UP/{cohort}_DOWN/{cohort}_MIXED")
        overlap = all_union.intersection(cohort_all)
        if overlap:
            examples = ",".join(sorted(overlap)[:10])
            errors.append(f"{cohort}_all is not disjoint from earlier cohort query sets; examples: {examples}")
        all_union.update(cohort_all)

    if all_union != unique_genes:
        missing = sorted(unique_genes - all_union)[:10]
        extra = sorted(all_union - unique_genes)[:10]
        errors.append(
            "union of cohort _all lists does not equal the full source gene set"
            f"; missing examples: {','.join(missing) or 'none'}"
            f"; extra examples: {','.join(extra) or 'none'}"
        )

    return errors


def main():
    args = parse_args()
    errors = []

    feature_path = os.path.normpath(args.features)
    retained_path = os.path.normpath(args.retained_source_contrasts) if args.retained_source_contrasts else ""
    outdir = os.path.normpath(args.outdir)

    paths_to_check = [("features", feature_path)]
    if retained_path:
        paths_to_check.append(("retained_source_contrasts", retained_path))
    for label, path in paths_to_check:
        if STALE_171_SENTINEL in path:
            errors.append(f"{label} points at a stale archived feature-set path: {path}")
    if STALE_RBL_QUERY_DIR in outdir:
        errors.append(f"outdir points at stale RBL query-set tree: {outdir}")

    if not os.path.exists(feature_path):
        errors.append(f"pan-cancer feature table does not exist: {feature_path}")
    if errors:
        return fail(errors)

    feature_columns, raw_feature_rows = read_tsv(feature_path)
    missing_feature_cols = [col for col in REQUIRED_COLUMNS if col not in feature_columns]
    if missing_feature_cols:
        errors.append("pan-cancer feature table missing required columns: " + ",".join(missing_feature_cols))
    if errors:
        return fail(errors)

    feature_rows = [normalise_row(row) for row in raw_feature_rows]
    source_table_rows = len(feature_rows)
    unique_genes = {row["gene_id"] for row in feature_rows}
    observed_feature_count = len(unique_genes)

    if source_table_rows != observed_feature_count:
        errors.append(
            f"pan-cancer feature table has {source_table_rows} rows but {observed_feature_count} unique gene_id values"
        )
    if any(not row["gene_id"] for row in feature_rows):
        errors.append("pan-cancer feature table contains blank gene_id values")

    owner_values = {row["owner_profile"] for row in feature_rows}
    direction_values = {row["direction"] for row in feature_rows}
    if owner_values - set(COHORTS):
        errors.append("owner_profile values outside BRCA/NBL/RBL: " + ",".join(sorted(owner_values - set(COHORTS))))
    if direction_values - set(DIRECTIONS):
        errors.append("direction values outside UP/DOWN/MIXED: " + ",".join(sorted(direction_values - set(DIRECTIONS))))
    if any(STALE_171_SENTINEL in row.get("source_file", "") for row in feature_rows):
        errors.append("feature source_file column references stale final_pan_cancer_feature_set__all/genes.tsv")
    if errors:
        return fail(errors)

    by_cohort = defaultdict(list)
    by_cohort_direction = defaultdict(list)
    for row in feature_rows:
        by_cohort[row["owner_profile"]].append(row)
        by_cohort_direction[(row["owner_profile"], row["direction"])].append(row)

    query_genes = {}
    for cohort in COHORTS:
        all_rows = sorted(by_cohort[cohort], key=rank_key)
        query_genes[query_id(cohort, "all")] = [row["gene_id"] for row in all_rows]
        for direction in DIRECTIONS:
            rows = sorted(by_cohort_direction[(cohort, direction)], key=rank_key)
            query_genes[query_id(cohort, direction)] = [row["gene_id"] for row in rows]

    source_counts = count_query_rows(feature_rows)

    retained_sha = ""
    retained_rows = []
    retained_row_count = ""
    retained_unique_gene_count = ""
    retained_crosscheck_status = "not_provided"
    if retained_path:
        if os.path.exists(retained_path):
            retained_columns, raw_retained_rows = read_tsv(retained_path)
            missing_retained_cols = [col for col in REQUIRED_COLUMNS if col not in retained_columns]
            if missing_retained_cols:
                errors.append(
                    "retained source contrast table missing required columns: "
                    + ",".join(missing_retained_cols)
                )
            if errors:
                return fail(errors)

            retained_rows = [normalise_row(row) for row in raw_retained_rows]
            retained_unique_genes = {row["gene_id"] for row in retained_rows}
            retained_row_count = str(len(retained_rows))
            retained_unique_gene_count = str(len(retained_unique_genes))
            retained_sha = sha256_file(retained_path)
            retained_crosscheck_status = "pass"

            if len(retained_rows) != source_table_rows:
                errors.append(
                    f"retained source contrast table has {len(retained_rows)} rows; source table has {source_table_rows}"
                )
            if len(retained_unique_genes) != observed_feature_count:
                errors.append(
                    "retained source contrast table has "
                    f"{len(retained_unique_genes)} unique gene_id values; source table has {observed_feature_count}"
                )
            if key_set(feature_rows) != key_set(retained_rows):
                errors.append("feature table rows do not exactly match retained source contrast cross-check keys")

            retained_counts = count_query_rows(retained_rows)
            for cohort, direction, name in iter_query_ids():
                if retained_counts.get(name, 0) != source_counts.get(name, 0):
                    errors.append(
                        f"{name} source count is {source_counts.get(name, 0)} "
                        f"but retained cross-check count is {retained_counts.get(name, 0)}"
                    )
        else:
            retained_crosscheck_status = "skipped_missing"

    if errors:
        return fail(errors)

    Path(outdir).mkdir(parents=True, exist_ok=True)
    output_paths = []
    gene_list_paths = {}
    for cohort, direction, name in iter_query_ids():
        path = os.path.join(outdir, f"{name}.txt")
        write_gene_list(path, query_genes[name])
        gene_list_paths[name] = path
        output_paths.append(path)

    output_genes = {name: read_gene_list(path) for name, path in gene_list_paths.items()}
    errors.extend(validate_output_lists(output_genes, query_genes, source_counts, unique_genes))
    for name, path in gene_list_paths.items():
        if os.path.getsize(path) == 0 and source_counts.get(name, 0) != 0:
            errors.append(f"{name}.txt is empty but the source table count is {source_counts.get(name, 0)}")
    if errors:
        return fail(errors)

    feature_sha = sha256_file(feature_path)
    manifest_path = os.path.join(outdir, "query_manifest.tsv")
    manifest_rows = []
    for cohort, direction, name in iter_query_ids():
        manifest_rows.append(
            {
                "query_id": name,
                "cohort": cohort,
                "direction": direction.upper() if direction != "all" else "all",
                "gene_count": str(len(query_genes[name])),
                "source_query_count": str(source_counts.get(name, 0)),
                "observed_feature_count": str(observed_feature_count),
                "source_table_rows": str(source_table_rows),
                "source_table_unique_genes": str(observed_feature_count),
                "source_table": feature_path,
                "source_table_sha256": feature_sha,
                "retained_source_contrasts": retained_path,
                "retained_source_contrasts_sha256": retained_sha,
                "retained_source_contrasts_rows": retained_row_count,
                "retained_source_contrasts_unique_genes": retained_unique_gene_count,
                "retained_source_crosscheck_status": retained_crosscheck_status,
                "path": os.path.join(outdir, f"{name}.txt"),
                "validation_status": "pass",
            }
        )
    write_tsv(
        manifest_path,
        [
            "query_id",
            "cohort",
            "direction",
            "gene_count",
            "source_query_count",
            "observed_feature_count",
            "source_table_rows",
            "source_table_unique_genes",
            "source_table",
            "source_table_sha256",
            "retained_source_contrasts",
            "retained_source_contrasts_sha256",
            "retained_source_contrasts_rows",
            "retained_source_contrasts_unique_genes",
            "retained_source_crosscheck_status",
            "path",
            "validation_status",
        ],
        manifest_rows,
    )
    output_paths.append(manifest_path)

    gmt_all_path = os.path.join(outdir, "gprofiler_by_cohort_all.gmt")
    write_gmt(
        gmt_all_path,
        [
            (
                query_id(cohort, "all"),
                f"PanCancerFeatureSet_owner_profile_{cohort};observed_feature_count={observed_feature_count}",
                query_genes[query_id(cohort, "all")],
            )
            for cohort in COHORTS
        ],
    )
    output_paths.append(gmt_all_path)

    gmt_directional_path = os.path.join(outdir, "gprofiler_by_cohort_directional.gmt")
    write_gmt(
        gmt_directional_path,
        [
            (
                query_id(cohort, direction),
                f"PanCancerFeatureSet_owner_profile_{cohort}_{direction};observed_feature_count={observed_feature_count}",
                query_genes[query_id(cohort, direction)],
            )
            for cohort in COHORTS
            for direction in DIRECTIONS
        ],
    )
    output_paths.append(gmt_directional_path)

    newline_errors = assert_nonempty_text_files_have_newline(output_paths)
    if newline_errors:
        return fail(newline_errors)

    print(f"features={feature_path}")
    print(f"source_table_sha256={feature_sha}")
    print(f"source_table_rows={source_table_rows}")
    print(f"source_table_unique_genes={observed_feature_count}")
    print(f"observed_feature_count={observed_feature_count}")
    print(f"retained_source_contrasts={retained_path}")
    print(f"retained_source_contrasts_sha256={retained_sha}")
    print(f"retained_source_crosscheck_status={retained_crosscheck_status}")
    print(f"output_dir={outdir}")
    for cohort, direction, name in iter_query_ids():
        print(f"{name}\t{len(query_genes[name])}")
    print("validation_status=pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
