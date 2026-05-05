# Enrichment Stage Implementation Report

## Scope

The pipeline now includes a manifest-driven Stage 7B for gene set enrichment after DESeq2 differential expression analysis and before marker post-processing. The implementation avoids hardcoded diseases, contrasts, components, and cell lines by discovering query sets from DESeq2 marker manifests, component marker outputs, recurrence outputs, and configured pan-cancer feature files.

## Implemented query families

The enrichment builder emits separate query records for:

- Per-contrast marker sets.
- Component or isolate aggregated marker sets.
- Disease recurrence-filtered marker sets at the configured recurrence threshold.
- The strict expression-filtered union.
- The operative pan-cancer feature set when its configured file is present.
- The operative-minus-strict comparison set.
- Ordered strict union up/down queries.
- Ordered pan-cancer up/down queries only when explicitly enabled and rank statistics are supplied.

Unranked queries are split into `all`, `up`, and `down` directions. Ordered queries are split only into `up_ordered` and `down_ordered`.

## Background and marker eligibility

The DESeq2 isolate and component outputs now write `normalised_count_in_test_sample`. Enrichment backgrounds are drawn from genes with non-missing adjusted P values and `normalised_count_in_test_sample >= 10`, matching the expression-eligible marker universe used for foreground selection.

## Configuration

The default configuration now enables enrichment, pins the g:Profiler organism and correction settings, defines the primary source list, requires an archive URL, keeps short-link storage disabled, enables IEA sensitivity analysis, and treats both the strict union and operative pan-cancer feature set as primary analyses.

The strict union, operative feature set, and strict-vs-operative comparison are controlled by semantic configuration flags (`strict_union_primary`, `operative_feature_set_primary`, and `compare_strict_vs_operative`) rather than by names that encode observed gene counts. Current defaults enable all three.

The R conda environment now requests `r-gprofiler2=0.2.4`.

## Outputs

The new stage writes:

- `query_sets/query_manifest.tsv`
- `query_sets/skipped_queries.tsv`
- `query_sets/gene_sets/{query_id}/genes.tsv`
- `query_sets/gene_sets/{query_id}/background.tsv`
- `query_sets/gene_sets/{query_id}/ranked_genes.tsv`
- `gprofiler/corpus_manifest.tsv`
- `gprofiler/top_terms.tsv`
- `gprofiler/iea_sensitivity_summary.tsv`
- `gprofiler/gprofiler_version.tsv`
- Per-query primary and IEA-sensitivity g:Profiler full, significant, and raw RDS outputs.

Skipped query rows remain in the manifest and are also copied to `skipped_queries.tsv` with explicit reasons such as `fewer_than_min_genes`, `missing_rank_statistics`, `missing_background`, and `empty_direction_after_split`.

## g:Profiler execution

The runner pins the configured archive with `set_base_url()`, uses a custom background, `g_SCS` correction, `evcodes = TRUE`, `significant = FALSE`, and `as_short_link = FALSE`. It writes a primary IEA-included result set and, when configured, an IEA-excluded sensitivity result set. Top terms are ranked independently within `query_id + source + iea_mode`.

## Files changed

- `Snakefile`
- `config/config.yaml`
- `envs/tcga-r-env.yaml`
- `scripts/deseq2_isolate_degs.R`
- `scripts/deseq2_component_vs_rest.R`
- `scripts/build_enrichment_query_sets.R`
- `scripts/run_gprofiler_from_manifest.R`

## Validation performed

Static R parse checks were run locally for the new and edited R scripts before installation on the remote host. A remote syntax and Snakemake dry-run check should be used for final assessment after the files are copied into place.

## Operational notes

Existing DESeq2 tables that predate this change will not contain `normalised_count_in_test_sample`. To obtain the corrected enrichment backgrounds, the DESeq2 isolate and component outputs should be regenerated before running Stage 7B.
