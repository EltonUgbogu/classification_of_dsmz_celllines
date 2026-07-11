#!/usr/bin/env python3
"""Write current graph-derived pan-cancer panel package reports.

Reports are generated from current feature, validation, threshold, enrichment,
graph, ranking, and multi-cohort package files. This script does not fabricate
validation records and does not encode a selected empirical rule.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
from pathlib import Path
import csv
import hashlib


ROOT = Path(__file__).resolve().parents[1]
PAN = ROOT / "results/unsupervised/pan_cancer"
MULTI = ROOT / "results/unsupervised/multicohort_cancer"
FS = PAN / "feature_space"
ENR = PAN / "enrichment/ranked_marker_source_panel"
DSMZ = PAN / "dsmz_similarity_graph/ranked_marker_source_panel"
RANK = PAN / "ranking/ranked_marker_source_panel"
EMBED = PAN / "embedding/ranked_marker_source_panel"
MOUT = MULTI / "ranked_marker_source_panel"
METHOD = "graph_derived_pan_cancer_feature_selection_v1_revised"
PREFIX = "ranked_marker_source_panel"


def read_tsv(path: Path, required: bool = False) -> list[dict[str, str]]:
    if not path.exists():
        if required:
            raise SystemExit(f"Required file is missing: {path}")
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def metrics(path: Path) -> dict[str, str]:
    return {row["metric"]: row["value"] for row in read_tsv(path)}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def count_field(rows: list[dict[str, str]], field: str) -> Counter:
    out: Counter = Counter()
    for row in rows:
        for value in str(row.get(field, "")).split(";"):
            value = value.strip()
            if value:
                out[value] += 1
    return out


def fmt_metrics(values: dict[str, str], keys: list[str]) -> str:
    lines = []
    for key in keys:
        lines.append(f"- {key}: {values.get(key, 'not_available')}")
    return "\n".join(lines)


feature = read_tsv(FS / "pan_cancer_features.tsv", required=True)
thresholds = read_tsv(FS / f"{PREFIX}_empirical_quantile_thresholds.tsv")
validation = read_tsv(FS / f"{PREFIX}_validation.tsv")
query = read_tsv(ENR / f"{PREFIX}_enrichment_query_manifest.tsv")
skipped = read_tsv(ENR / f"{PREFIX}_enrichment_skipped_queries.tsv")
dsmz = metrics(DSMZ / f"{PREFIX}_dsmz_graph_metrics.tsv")
multi = metrics(MOUT / f"{PREFIX}_multicohort_summary_metrics.tsv")
patient = metrics(RANK / f"{PREFIX}_patient_ranking_metrics.tsv")
balanced = metrics(RANK / f"{PREFIX}_retrieval_balanced_accuracy.tsv")
same_cancer_type_ranking = read_tsv(RANK / f"{PREFIX}_retrieval_accuracy.tsv")

owners = Counter(row.get("owner_profile", "") for row in feature)
classes = Counter(row.get("feature_class", "") for row in feature)
directions = Counter(row.get("direction", "") for row in feature)
selection_classes = count_field(feature, "selection_basis_classes")
marker_sources = count_field(feature, "marker_source_class")
validation_status = Counter(row.get("status", "") for row in validation)
final_gene_count = len({row.get("gene_id", "") for row in feature if row.get("gene_id", "")})
all_validation_pass = validation and all(row.get("status") == "PASS" for row in validation)

threshold_lines = "\n".join(
    "- {cancer_type}/{marker_evidence_stratum}/{candidate_pool_type}: "
    "p<={adjusted_p_value_quantile_threshold}, "
    "abs_shrunken_lfc>={absolute_shrunken_log2fc_quantile_threshold}, "
    "expression>={expression_quantile_threshold}".format(**row)
    for row in thresholds
) or "- not_available"

write(
    ENR / f"{PREFIX}_gprofiler_manual_execution_instructions.md",
    f"""# Pan-Cancer Panel g:Profiler Manual Execution

Use `{PREFIX}_enrichment_query_manifest.tsv` as the execution manifest. The
queries are derived from the current graph-derived pan-cancer feature panel and
the retained-marker evidence backgrounds prepared for each query. This script
does not assert that g:Profiler has completed unless current output files show
that status.

- Queries prepared: {len(query)}
- Queries meeting size threshold: {len(query) - len(skipped)}
- Queries skipped by size: {len(skipped)}
""",
)

write(
    ENR / f"{PREFIX}_enrichment_report.md",
    f"""# Pan-Cancer Panel Enrichment Report

The enrichment query manifest was generated from the current feature table.
Query categories should reflect cohort ownership, marker-source provenance,
feature class, direction, final panel membership, and retained per-contrast
marker sets.

- Method: `{METHOD}`
- Feature genes: {final_gene_count}
- Queries prepared: {len(query)}
- Runnable queries by size: {len(query) - len(skipped)}
- Skipped queries by size: {len(skipped)}

Functional enrichment should be interpreted as contextual annotation for the
selected marker panel, not as causal evidence of pathway regulation.
""",
)

write(
    DSMZ / f"{PREFIX}_dsmz_similarity_graph_report.md",
    f"""# Pan-Cancer Panel DSMZ Similarity Graph Report

The DSMZ graph package is summarised from current package metrics.

{fmt_metrics(dsmz, [
    "method",
    "feature_genes_requested",
    "feature_genes_used",
    "nodes",
    "edges",
    "connected_components",
    "communities",
    "modularity",
    "cancer_type_assortativity",
])}

Communities describe structure in a marker-derived expression space and should
not be described as formal molecular subtypes.
""",
)

overall = next((row for row in same_cancer_type_ranking if row.get("group") == "overall"), {})
write(
    RANK / f"{PREFIX}_ranking_report.md",
    f"""# Pan-Cancer Panel Ranking Report

Similarity ranking was summarised from current package metrics.

{fmt_metrics(patient, [
    "top1_accuracy",
    "top10_accuracy",
    "mrr",
])}

- Cell-line-centred top-1 agreement: {overall.get("accuracy", "not_available")}
- Balanced top-1 agreement: {balanced.get("balanced_accuracy", "not_available")}

The rankings define model-prioritisation evidence in the selected feature
space. They do not establish drug response, perturbation response, xenograft
behaviour, essentiality, clonal composition, microenvironmental fidelity, or
full biological equivalence.
""",
)

write(
    EMBED / f"{PREFIX}_embedding_report.md",
    f"""# Pan-Cancer Panel Embedding Report

The embedding package contains UMAP outputs generated from the current
graph-informed DESeq2 marker-derived pan-cancer feature panel. The embedding is
a diagnostic visualisation and is not used to define formal molecular subtypes.
""",
)

write(
    MOUT / f"{PREFIX}_multicohort_report.md",
    f"""# Pan-Cancer Panel Multi-Cohort Report

The multi-cohort package is summarised from current package metrics.

{fmt_metrics(multi, [
    "method",
    "cohorts_included",
    "samples_profiles",
    "feature_genes_requested",
    "feature_genes_found",
    "feature_genes_missing",
    "final_genes_used",
    "graph_nodes",
    "graph_edges",
    "connected_components",
    "communities",
    "modularity",
])}

These outputs support feature-space comparison and model prioritisation, not
claims of universal cancer biology or complete model fidelity.
""",
)

write(
    PAN / f"{PREFIX}_manuscript_update_notes.md",
    f"""# Pan-Cancer Panel Manuscript Update Notes

No manuscript or LaTeX file was edited by this report finaliser.

## Methods

The active method is `{METHOD}`. Contrast-level DESeq2 marker selection is
upstream of this panel builder. The final panel is assembled from recurrent
genes retained directly, accepted singleton candidates, and accepted
non-recurrent candidates, where candidate acceptance requires statistical
evidence, effect-magnitude evidence, and expression evidence.

## Current Panel

- Final genes: {final_gene_count}
- Cohort owners: {dict(sorted(owners.items()))}
- Directions: {dict(sorted(directions.items()))}
- Feature classes: {dict(sorted(classes.items()))}
- Selection-basis memberships: {dict(sorted(selection_classes.items()))}
- Marker-source memberships: {dict(sorted(marker_sources.items()))}
- Validation status counts: {dict(sorted(validation_status.items()))}

## Empirical Threshold Contexts

{threshold_lines}

Interpretation should remain proportional: the panel is a graph-informed,
contrast-derived marker feature space for downstream expression, ranking,
network, community, and enrichment analyses.
""",
)

validation_rows = [
    {
        "check": "feature_table_present",
        "status": "PASS" if feature else "FAIL",
        "detail": f"rows={len(feature)} unique_genes={final_gene_count}",
    },
    {
        "check": "feature_validation_loaded",
        "status": "PASS" if validation else "FAIL",
        "detail": f"validation_rows={len(validation)} status_counts={dict(validation_status)}",
    },
    {
        "check": "feature_validation_all_passed",
        "status": "PASS" if all_validation_pass else "FAIL",
        "detail": f"all_validation_pass={all_validation_pass}",
    },
    {
        "check": "threshold_table_loaded",
        "status": "PASS" if thresholds else "FAIL",
        "detail": f"threshold_contexts={len(thresholds)}",
    },
    {
        "check": "canonical_method",
        "status": "PASS" if {row.get("method") for row in feature} == {METHOD} else "FAIL",
        "detail": ";".join(sorted({row.get("method", "") for row in feature})),
    },
]
write_tsv(PAN / f"{PREFIX}_downstream_validation_checks.tsv", validation_rows)

final_report = f"""# Pan-Cancer Panel Downstream Report

## Summary

- Method: `{METHOD}`
- Final genes: {final_gene_count}
- Feature classes: {dict(sorted(classes.items()))}
- Validation all passed: {all_validation_pass}
- Queries prepared: {len(query)}
- DSMZ graph metrics available: {bool(dsmz)}
- Ranking metrics available: {bool(patient)}
- Multi-cohort metrics available: {bool(multi)}

## Method

The active feature panel follows one route: canonical retained marker lists,
recurrence counting within cancer type and marker-evidence stratum, recurrent
direct retention, all-three empirical acceptance for singleton and
non-recurrent candidates, final union, and deterministic provenance-preserving
deduplication.

## Downstream Interpretation

Expression, ranking, UMAP, graph, community, and enrichment outputs consume the
current clean feature list. They should be described as analyses in a
graph-informed DESeq2 marker-derived pan-cancer feature space.
"""
write(PAN / f"{PREFIX}_downstream_rerun_report.md", final_report)

manifest_roots = [ENR, DSMZ, RANK, EMBED, MOUT]
root_files = [
    PAN / f"{PREFIX}_downstream_rerun_report.md",
    PAN / f"{PREFIX}_downstream_validation_checks.tsv",
    PAN / f"{PREFIX}_manuscript_update_notes.md",
]
manifest_rows = []
seen: set[Path] = set()
for base in manifest_roots:
    for path in sorted(base.rglob("*")):
        if path.is_file() and path.name != f"{PREFIX}_downstream_rerun_manifest.tsv":
            seen.add(path)
for path in root_files:
    if path.exists():
        seen.add(path)
for path in sorted(seen):
    manifest_rows.append(
        {
            "path": str(path.relative_to(ROOT)),
            "size_bytes": path.stat().st_size,
            "modification_time": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(),
            "sha256": sha(path),
        }
    )
write_tsv(
    PAN / f"{PREFIX}_downstream_rerun_manifest.tsv",
    manifest_rows,
    ["path", "size_bytes", "modification_time", "sha256"],
)
print(f"reports and manifest complete: {len(manifest_rows)} files")
