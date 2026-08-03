# Pipeline Profiles And Output Namespaces

This note defines the boundary between selectable Snakemake
`pipeline_profile` values and output directory namespaces.

| Name | Type | Meaning |
|---|---|---|
| `multicohort_cancer` | `pipeline_profile` | Multicohort patient-referenced analysis context for joint BRCA/NBL/RBL inputs, tumour-neighbourhood representations, patient-referenced graph products, and post-hoc cancer-type interpretation. |
| `pan_cancer` | `pipeline_profile` and output namespace | Marker-derived pan-cancer analysis layer for feature-space, expression-matrix, tumour/cell-line mapping, ranking-diagnostic, enrichment-query, graph, and community outputs under `results/unsupervised/pan_cancer/`. |
| `brca` / `nbl` / `rbl` | `pipeline_profile` | Single-cohort analysis profiles for disease-specific clustering, graphs, marker calls, and supporting artefacts. |
| `heme` | `pipeline_profile` | Haematological benchmark profile present in the configuration. It is not part of the current solid-tumour pan-cancer feature panel. |

`multicohort_cancer` and `pan_cancer` are not interchangeable. The former is a
patient-referenced multicohort context; the latter is the pan-cancer
analysis/output layer.

## Current Routing

The selectable profile keys in `config/config.yaml` are now `multicohort_cancer`,
`pan_cancer`, `brca`, `nbl`, `rbl`, and `heme`.

Use `pipeline_profile=pan_cancer` for the marker-derived pan-cancer target family:

- disease-level marker outputs from `brca`, `nbl`, and `rbl`;
- the pan-cancer feature set in `results/unsupervised/pan_cancer/feature_space/`;
- pan-cancer feature-limited expression matrices in `results/unsupervised/pan_cancer/inputs/`;
- tumour/cell-line mapping and graph products in `results/unsupervised/pan_cancer/`;
- bidirectional ranking diagnostics in `results/unsupervised/pan_cancer/ranking/diagnostics/`;
- pan-cancer enrichment query sets in `results/unsupervised/pan_cancer/enrichment/query_sets/`.

For backward compatibility, the Snakefile may still declare pan-cancer rules
when `pipeline_profile=multicohort_cancer` is used with explicit pan-cancer
targets. That compatibility path is routing only. It does not make
`multicohort_cancer` the conceptual profile for pan-cancer analysis, and
pan-cancer targets are not part of the multicohort default target list.

## Multicohort Graph Interpretation Boundary

The `multicohort_cancer` graph is unsupervised with respect to cancer type.
Cancer-type labels are not used to construct the graph; they are applied after
graph construction to interpret connected components and cell-line
neighbourhoods. Retinoblastoma-labelled cell lines should therefore be described
as emerging as a connected group in the unsupervised multicohort similarity
space, not as a component caused by cancer-type information supplied to the
model.

The resolved-neighbour graph uses a global/local direction-intersection rule:
each cell line retains neighbours present in both the globally best-performing
representation and that cell line's winning representation or representations.
The majority edge-recurrence consensus graph is a separate product that filters
edges by recurrence across feature-distance representations.

## Dry-Run Examples

Dry-run a pan-cancer ranking diagnostic through the pan-cancer profile:

```bash
snakemake \
  -n -p \
  results/unsupervised/pan_cancer/ranking/diagnostics/ranking_diagnostic_metric_crosscheck.tsv \
  --config pipeline_profile=pan_cancer
```

Dry-run the multicohort patient-referenced target set through the multicohort
profile:

```bash
snakemake \
  --config pipeline_profile=multicohort_cancer \
  --rerun-triggers mtime \
  -n -p
```

## Future Cleanup

- Retire the `multicohort_cancer` compatibility route for explicit pan-cancer
  targets once active command examples and automation use `pipeline_profile=pan_cancer`.
- Avoid renaming current result directories until released outputs and
  provenance records have been updated or retired.
