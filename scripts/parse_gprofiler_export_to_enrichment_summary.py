#!/usr/bin/env python3
"""Convert a manual g:Profiler export into a heatmap-ready ranked marker-source-panel top-term summary.

This fallback parser does not run g:Profiler; it attaches functional-enrichment
query-manifest metadata to an externally downloaded term-level export.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

QUERY_COLS = ("query_id", "query_name", "query", "name")
TERM_NAME_COLS = ("term_name", "name", "term", "description")
TERM_ID_COLS = ("term_id", "native", "term")
P_VALUE_COLS = ("p_value", "p.value", "pvalue", "adjusted_p_value")
SOURCE_COLS = ("source", "database")
INTERSECTION_COLS = ("intersection", "intersections", "intersection_genes")

OUTPUT_COLUMNS = (
    "query_id",
    "query_name",
    "category",
    "query_family",
    "cohort",
    "contrast",
    "profile",
    "disease",
    "group_id",
    "direction",
    "source_marker_table",
    "gene_list_path",
    "background_path",
    "ordered",
    "rank_source",
    "source",
    "term_id",
    "term_name",
    "p_value",
    "intersection_size",
    "query_size",
    "gene_count",
    "background_count",
    "term_size",
    "effective_domain_size",
    "precision",
    "recall",
    "enrichment_ratio",
    "intersection",
    "top_intersection_genes",
    "rank_within_query",
    "rank_within_source",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Manual g:Profiler CSV/TSV export")
    parser.add_argument(
        "--query-manifest",
        required=True,
        help="functional-enrichment query manifest, such as ranked_marker_source_panel_enrichment_query_manifest.tsv",
    )
    parser.add_argument("--output", required=True, help="Heatmap-ready summary TSV")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="replace")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = sniff_delimiter(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        columns = reader.fieldnames or []
    return columns, rows


def pick(columns: list[str], candidates: tuple[str, ...]) -> str:
    lower_to_actual = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_actual:
            return lower_to_actual[candidate.lower()]
    return ""


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def main() -> int:
    args = parse_args()
    export_path = Path(args.input)
    manifest_path = Path(args.query_manifest)
    output_path = Path(args.output)
    if not export_path.exists():
        print(f"ERROR: g:Profiler export missing: {export_path}", file=sys.stderr)
        return 1
    if not manifest_path.exists():
        print(f"ERROR: query manifest missing: {manifest_path}", file=sys.stderr)
        return 1

    manifest_columns, manifest_rows = read_table(manifest_path)
    manifest_by_query = {}
    for row in manifest_rows:
        q = row.get("query_name") or row.get("query_id")
        if q:
            manifest_by_query[q] = row

    export_columns, export_rows = read_table(export_path)
    query_col = pick(export_columns, QUERY_COLS)
    source_col = pick(export_columns, SOURCE_COLS)
    term_id_col = pick(export_columns, TERM_ID_COLS)
    term_name_col = pick(export_columns, TERM_NAME_COLS)
    p_value_col = pick(export_columns, P_VALUE_COLS)
    intersection_col = pick(export_columns, INTERSECTION_COLS)
    required = {
        "query": query_col,
        "source": source_col,
        "term_name": term_name_col,
        "p_value": p_value_col,
    }
    missing = [name for name, col in required.items() if not col]
    if missing:
        print(
            "ERROR: export missing required interpretable columns: " + ",".join(missing),
            file=sys.stderr,
        )
        print("Observed columns: " + ",".join(export_columns), file=sys.stderr)
        return 1

    def value(row: dict[str, str], *names: str) -> str:
        for name in names:
            if name in row and str(row.get(name) or "") != "":
                return str(row.get(name) or "")
        return ""

    summary_rows = []
    for row in export_rows:
        q = str(row.get(query_col) or "").strip()
        if not q:
            continue
        meta = manifest_by_query.get(q, {})
        category = meta.get("category", meta.get("query_family", ""))
        cohort = meta.get("cohort", meta.get("owner_profile", meta.get("profile", meta.get("disease", ""))))
        contrast = meta.get("contrast", meta.get("source_contrast", meta.get("contrast_id", "")))
        source_marker_table = meta.get(
            "source_marker_table_path",
            meta.get("source_marker_table", meta.get("source_table", meta.get("marker_file", ""))),
        )
        gene_list_path = meta.get("gene_list_path", meta.get("genes_path", meta.get("genes_tsv", "")))
        background_path = meta.get("background_path", meta.get("background_tsv", ""))
        source = str(row.get(source_col) or "").strip()
        term_name = str(row.get(term_name_col) or "").strip()
        p_value = str(row.get(p_value_col) or "").strip()
        if not source or not term_name or not p_value:
            continue
        intersection = str(row.get(intersection_col) or "").strip() if intersection_col else ""
        summary_rows.append({
            "query_id": q,
            "query_name": q,
            "category": category,
            "query_family": meta.get("query_family", category),
            "cohort": cohort,
            "contrast": contrast,
            "profile": meta.get("profile", cohort),
            "disease": meta.get("disease", cohort),
            "group_id": meta.get("group_id", category),
            "direction": meta.get("direction", ""),
            "source_marker_table": source_marker_table,
            "gene_list_path": gene_list_path,
            "background_path": background_path,
            "ordered": "FALSE",
            "rank_source": meta.get("rank_source", ""),
            "source": source,
            "term_id": str(row.get(term_id_col) or "").strip() if term_id_col else "",
            "term_name": term_name,
            "p_value": p_value,
            "intersection_size": value(row, "intersection_size", "intersection.count", "term_intersection_size"),
            "query_size": value(row, "query_size", "query.size"),
            "gene_count": meta.get("gene_count", value(row, "query_size", "query.size")),
            "background_count": meta.get("background_count", ""),
            "term_size": value(row, "term_size", "term.size"),
            "effective_domain_size": value(row, "effective_domain_size", "effective.domain.size"),
            "precision": value(row, "precision"),
            "recall": value(row, "recall"),
            "enrichment_ratio": value(row, "enrichment_ratio"),
            "intersection": intersection,
            "top_intersection_genes": intersection,
            "rank_within_query": "",
            "rank_within_source": "",
        })

    if not summary_rows:
        print("ERROR: no usable rows parsed from g:Profiler export", file=sys.stderr)
        return 1
    write_tsv(output_path, OUTPUT_COLUMNS, summary_rows)
    provenance_path = output_path.with_suffix(output_path.suffix + ".provenance.tsv")
    # This provenance describes the parser path only. Live g:Profiler run
    # provenance is recorded by the runner output directory, not by this
    # manual-export fallback.
    write_tsv(
        provenance_path,
        ("metric", "value"),
        [
            {"metric": "manual_gprofiler_export", "value": str(export_path)},
            {"metric": "manual_gprofiler_export_sha256", "value": sha256_file(export_path)},
            {"metric": "query_manifest", "value": str(manifest_path)},
            {"metric": "query_manifest_sha256", "value": sha256_file(manifest_path)},
            {"metric": "parsed_rows", "value": len(summary_rows)},
            {"metric": "gprofiler_run", "value": "not_run_by_pipeline"},
        ],
    )
    print(f"manual_gprofiler_export={export_path}")
    print(f"query_manifest={manifest_path}")
    print(f"summary={output_path}")
    print(f"rows={len(summary_rows)}")
    print("gprofiler_run=not_run_by_pipeline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
