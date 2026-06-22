# Supplementary data

This directory contains curated supplementary outputs for the active pan-cancer analysis, including the current 379-gene ranked marker-source pan-cancer feature panel and associated ranking, graph, enrichment, data-description, and manifest tables.

| Directory | Contents |
| --- | --- |
| `feature_space/` | Current pan-cancer feature panel and marker-source annotations. |
| `ranking/` | Bidirectional tumour-to-cell-line and cell-line-to-tumour ranking metrics. |
| `tumour_neighbourhood/` | Patient-neighbourhood and graph-resolution outputs. |
| `pan_cancer_similarity/` | DSMZ cell-line graph, network, and community outputs. |
| `model_prioritisation/` | Model-prioritisation scores and metadata. |
| `component_composition/` | Top-50 tumour-component composition summaries for ranking outputs. |
| `enrichment/` | Enrichment query/background manifests and endpoint-status outputs. |
| `data_description/` | Input descriptions, retained sample identifiers, exclusions, and processing notes. |
| `manifests/` | File manifests, checksums, and figure provenance. |
| `archive/` | Historical or provenance-only outputs not used as the active feature space. |

## Interpretation notes

- `precision@k` reports, for each cell-line group, the fraction of the top `k` retrieved patient tumour samples that share the cell line's cancer lineage. In `ranking/precision_at_k.tsv`, `k` is a top-N tumour-sample cutoff after ranking all tumour samples for that cell-line group; `k_type` is `absolute` for fixed cutoffs and `percentile` for cutoffs derived from 1%, 5%, or 10% of the tumour sample universe.
- `model_prioritisation/` compares all cell-line groups against patient tumour samples from the target cancer type named in `cancer_type`. Same-lineage rows are within that cancer type; outside-lineage rows are comparators.
- `tumour_neighbourhood/` includes adaptive tumour-neighbourhood exports: patient-neighbourhood membership, resolved cell-line neighbours, component assignments, and consensus component assignments.

## Path and checksum policy

Paths are repository-relative unless a table explicitly stores source identifiers. Entries in cohort-level `tmp` subdirectories below `data/` are repository-relative workflow intermediates, not operating-system temporary-directory paths.

Checksums in `manifests/checksums_sha256.tsv` refer only to curated supplementary files listed in that manifest. They do not certify raw, controlled-access, private, ignored, or unlisted source files.

Empty table values are encoded as `NA` when the value is not applicable or unavailable.
