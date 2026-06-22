# Pipeline Profiles And Output Namespaces

This note records the current boundary between selectable Snakemake
`pipeline_profile` values and output directory namespaces.

| Name | Type | Meaning |
|---|---|---|
| `multicohort_cancer` | `pipeline_profile` | Patient-referenced multicohort workflow context. It builds the joint BRCA/NBL/RBL inputs, patient-referenced multicohort graph products, and currently loads downstream cross-cohort/pan-cancer rules. |
| `pan_cancer` | output namespace | Downstream pan-cancer feature-space, expression-matrix, ranking-diagnostic, tumour/cell-line mapping, graph, and enrichment-query outputs under `results/unsupervised/pan_cancer/`. |
| `brca` / `nbl` / `rbl` | `pipeline_profile` | Single-cohort analysis profiles for disease-specific clustering, graphs, marker calls, and supporting artefacts. |
| `heme` | `pipeline_profile` | Single-cohort haematological profile present in the configuration. It is not part of the current solid tumour pan-cancer feature set. |

Do not interpret `pipeline_profile=multicohort_cancer` as meaning that every
target it can build is the multicohort patient-referenced graph. Some downstream
outputs under `results/unsupervised/pan_cancer/` are currently routed through this
profile because they depend on cross-cohort artefacts, but `pan_cancer` remains a
separate output namespace and analysis layer.

## Current Routing

There is currently no selectable `pipeline_profile=pan_cancer` in
`config/config.yaml`. The configured profile keys are `multicohort_cancer`,
`brca`, `nbl`, `rbl`, and `heme`.

In the current `Snakefile`, the marker post-processing branch is enabled only
when `profile_name == "multicohort_cancer"`. That branch builds or consumes:

- disease-level marker outputs from `brca`, `nbl`, and `rbl`;
- the pan-cancer feature set in `results/unsupervised/pan_cancer/feature_space/`;
- the pan-cancer feature-restricted expression matrices in `results/unsupervised/pan_cancer/inputs/`;
- tumour/cell-line mapping and graph products in `results/unsupervised/pan_cancer/`;
- bidirectional ranking diagnostics in `results/unsupervised/pan_cancer/ranking/diagnostics/`;
- pan-cancer enrichment query sets in `results/unsupervised/pan_cancer/enrichment/query_sets/`.

This makes `multicohort_cancer` the current workflow entrypoint for several
downstream cross-cohort products. It does not make `pan_cancer` and
`multicohort_cancer` interchangeable names.

## Multicohort Graph Interpretation Boundary

The `MULTICOHORT_CANCER` graph is unsupervised with respect to cancer type.
Cancer-type labels are not used to construct the graph; they are applied after
graph construction to interpret the resulting connected components and cell-line
neighbourhoods. Retinoblastoma-labelled cell lines should therefore be described
as emerging as a connected group in the unsupervised multicohort similarity
space, not as a component caused by cancer-type information supplied to the
model.

The resolved-neighbour graph is a stability-based edge filter. Unlike the union
graph, it retains only edges that meet the predefined support/intersection
criteria across feature–distance representations. This reduces the risk of
overfitting to representation-specific edges, but does not prove that overfitting
is absent. It should not be described as optimal unless it has been formally
optimised against a stated criterion.

## Dry-Run Example

Use the multicohort profile when dry-running the current bidirectional ranking
diagnostics target:

```bash
snakemake \
  --config pipeline_profile=multicohort_cancer \
  --rerun-triggers mtime \
  -n -p \
  results/unsupervised/pan_cancer/ranking/diagnostics/ranking_diagnostic_metric_crosscheck.tsv
```

This is a routing and entrypoint detail, not an identity between
`multicohort_cancer` and `pan_cancer`.

## TODO

- Consider introducing a clearer `pipeline_profile=pan_cancer` or a dedicated documented wrapper/entrypoint for downstream pan-cancer feature-space, ranking, and enrichment targets.
- If adding a `pan_cancer` profile would be too disruptive, document official command examples for pan-cancer targets that currently route through `pipeline_profile=multicohort_cancer`.
- Avoid renaming current result directories until thesis outputs and audit trails are final.
