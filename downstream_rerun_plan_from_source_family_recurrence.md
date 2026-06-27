# Downstream Rerun Plan From Source-Family Recurrence

## Context

The official pan-cancer feature set has been replaced by the source-family recurrence build in:

`results/unsupervised/pan_cancer/feature_space/`

The new active feature set contains 125 unique clean Ensembl IDs. Any downstream output that consumed either `pan_cancer_features.tsv` or `pan_cancer_features_clean.txt` from the previous 177-gene result is stale until rerun.

## Dry-Run Status

Dry-run log:

`/work/ugbogu/pipeline/downstream_rerun_dryrun.log`

Command used:

```bash
cd /work/ugbogu/pipeline && /work/ugbogu/.conda/envs/smk/bin/snakemake -n --cores 1 --config pipeline_profile=multicohort_cancer --rerun-triggers mtime -- /work/ugbogu/pipeline/results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds /work/ugbogu/pipeline/results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/query_manifest.tsv results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/query_manifest.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap/summary_pan_cancer_tumour_cell_line_alignment_umap.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_graph_edges.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_communities.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv
```

Dry-run result: PASS. Snakemake identified 9 local jobs to rerun.

## Local Targets Safe To Rerun

These are feature-set-dependent and do not require internet/API access:

- `build_pan_cancer_feature_expression_matrix`
- `build_pan_cancer_feature_expression_matrix_cell_lines_only`
- `build_gprofiler_by_cohort_query_sets`
- `build_pan_cancer_enrichment_marker_framework_query_sets`
- `plot_pan_cancer_tumour_cell_line_alignment_umap`
- `build_pan_cancer_cell_line_similarity_graph`
- `compute_pan_cancer_cell_line_communities`
- `compute_pan_cancer_cell_line_leiden_resolution_sweep`
- `compute_pan_cancer_cell_line_louvain_resolution_sweep`

Recommended combined local rerun command after backing up existing downstream outputs:

```bash
cd /work/ugbogu/pipeline && /work/ugbogu/.conda/envs/smk/bin/snakemake --cores 8 --config pipeline_profile=multicohort_cancer --rerun-triggers mtime -- /work/ugbogu/pipeline/results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds /work/ugbogu/pipeline/results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/query_manifest.tsv results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/query_manifest.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap/summary_pan_cancer_tumour_cell_line_alignment_umap.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_graph_edges.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_communities.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv /work/ugbogu/pipeline/results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv
```

## g:Profiler And Enrichment

Do not reuse g:Profiler outputs from the old 177-gene set.

The local marker-framework query/background preparation is safe to rerun now. After that, fresh g:Profiler results are required for the regenerated query sets. This step requires internet/API or manual web upload and should not be run automatically without deciding the submission method and backing up existing enrichment outputs.

After the fresh export is placed at:

`results/unsupervised/pan_cancer/enrichment/gprofiler_exports/marker_framework_gprofiler_export.csv`

rerun:

```bash
cd /work/ugbogu/pipeline && /work/ugbogu/.conda/envs/smk/bin/snakemake --cores 4 --config pipeline_profile=multicohort_cancer --rerun-triggers mtime -- results/unsupervised/pan_cancer/enrichment/figures/marker_framework/marker_framework_enrichment_summary_top_terms.tsv
```

## Backup Recommendation

Before executing the combined local rerun, back up these downstream directories if the current outputs need to remain reproducible:

- `results/unsupervised/pan_cancer/inputs/`
- `results/unsupervised/pan_cancer/enrichment/query_sets/gprofiler_by_cohort/`
- `results/unsupervised/pan_cancer/enrichment/query_sets/marker_framework/`
- `results/unsupervised/pan_cancer/enrichment/gprofiler_exports/`
- `results/unsupervised/pan_cancer/enrichment/figures/marker_framework/`
- `results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap/`
- `results/unsupervised/pan_cancer/cell_line_similarity/`

## Target Table

Machine-readable target details are in:

`/work/ugbogu/pipeline/downstream_rerun_targets.tsv`

