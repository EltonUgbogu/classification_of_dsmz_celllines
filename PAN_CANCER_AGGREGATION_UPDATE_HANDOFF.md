# Pan-Cancer Aggregation Update Handoff

Last updated: 2026-07-11 Europe/Berlin, after core builder/config/Snakefile refactor phase.

## Repository State

- Repository root: `/Users/eltonugbogu/classification_of_dsmz_celllines`
- Current branch: `version2-probability-graph`
- Current HEAD: `13592660ac3eea3b62d5b31b053fbd9eb3decee8`
- Complete git status at start:

```text
## version2-probability-graph...origin/version2-probability-graph
 M .gitignore
 M README.md
 M Snakefile
 M config/config.yaml
 M docs/version1/PIPELINE_DOCUMENTATION.md
 M scripts/build_enrichment_query_sets.R
 M scripts/build_enrichment_summary_top_terms.R
 M scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
 M scripts/build_pan_cancer_cell_line_similarity_graph.R
 M scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
 M scripts/build_pan_cancer_features.py
 M scripts/build_pan_cancer_graph.R
 M scripts/build_ranked_marker_source_panel_enrichment_queries.py
 M scripts/compute_cell_line_louvain_resolution_sweep.R
 M scripts/compute_pan_cancer_cell_line_communities.R
 M scripts/compute_pan_cancer_cell_line_validation.R
 M scripts/compute_pan_cancer_communities.R
 M scripts/deseq2_component_vs_rest.R
 M scripts/deseq2_isolate_degs.R
 M scripts/finalize_ranked_marker_source_panel_reports.py
 M scripts/package_ranked_marker_source_panel_outputs.py
 M scripts/parse_gprofiler_export_to_enrichment_summary.py
 M scripts/plot_deg_eval_pan_cancer.R
 M scripts/plot_enrichment_top_terms_heatmap.R
 M scripts/plot_pan_cancer_two_panel.R
 M scripts/prepare_deseq2_inputs.R
 M scripts/run_gprofiler_from_manifest.R
 M scripts/select_gprofiler_terms.py
 M scripts/validate_deseq2_staged_inputs.R
 M scripts/write_deseq2_directional_marker_tables.R
?? .cursor/
?? .github/
?? AGENT.md
?? PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
?? SRP409177.runinfo.csv
?? archive/
?? docs/version1/DGE_REPORT.md
?? docs/version1/ecdf.md
?? docs/version2/version2_implementation_status.md
?? nbl_21_fastq_fast_check.tsv
?? nbl_missing_valid_count_samples.txt
?? preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt
?? scripts/derive_resolved_graph_node_statistics.py
?? supplementary/
?? sync_logs/
?? test_logs/
```

- Unrelated modifications already exist: yes. The working tree was dirty before this update began, including multiple files that are in scope for this mission. Do not reset or discard these changes.
- Initial diff stat at start: 30 tracked files changed, 4653 insertions, 4733 deletions.

## Files Inspected So Far

- `Snakefile`
- `config/config.yaml`
- `scripts/build_pan_cancer_features.py`
- `README.md` through search results only
- `rules/`, `scripts/`, `docs/`, `supplementary/` through filename and stale-term searches
- `rules/pan_cancer_enrichment_marker_framework.smk`
- `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`
- `scripts/build_gprofiler_query_sets_from_pan_cancer_features.py`
- `scripts/build_ranked_marker_source_panel_enrichment_queries.py`
- `scripts/package_ranked_marker_source_panel_outputs.py`
- `scripts/finalize_ranked_marker_source_panel_reports.py`
- `scripts/build_pan_cancer_expression_matrix.R`
- `scripts/score_tumour_cellline_mapping.R`

## Active Dependency Chain, Initial Verification

Current active wiring still includes the legacy pan-cancer feature construction rule:

1. Canonical disease-profile contrast marker manifests feed `construct_pan_cancer_feature_panel`.
2. `construct_pan_cancer_feature_panel` calls `scripts/build_pan_cancer_features.py`.
3. The rule currently declares final feature outputs, gene evidence, legacy ranked/audit outputs, sensitivity summary, selected marker-source-class rows, validation, previous177/previous125 comparisons, run manifest, and active-directory manifest.
4. `PAN_FEATURES_CLEAN` is consumed by expression matrix construction, UMAP, tumour/cell-line mapping, pan-cancer graph/community workflows, enrichment query builders, package outputs, and report finalisers.
5. `rules/pan_cancer_enrichment_marker_framework.smk` remains included from `Snakefile` and must be inspected before final schema changes.

Refined status: the pan-cancer feature-panel rule has been rewired to the revised builder outputs. `PAN_FEATURES_CLEAN` still drives expression, UMAP, mapping, graph/community, enrichment, package, and report products. The marker-framework enrichment route remains active and still requires rewrite. The method-specific enrichment builder and final report/package scripts remain active stale consumers and still require rewrite.

## Locked Revised Method

The implementation target is fixed:

1. Use canonical retained marker manifests and retained per-contrast marker tables/lists for BRCA, NBL, and RBL.
2. Count recurrence separately within `cancer_type x marker_evidence_stratum`.
3. Use invariant `RECURRENCE_MIN_COUNT = 2`.
4. Classify within-context rows as:
   - `recurrent`: `retained_marker_list_count >= 2`
   - `singleton_candidate`: `retained_marker_list_total == 1` and `retained_marker_list_count == 1`
   - `non_recurrent_candidate`: `retained_marker_list_total > 1` and `retained_marker_list_count == 1`
5. Do not admit zero-recurrence rows to candidate pools.
6. Compute empirical thresholds only for singleton and non-recurrent candidates within `cancer_type x marker_evidence_stratum x candidate_pool_type`.
7. Use quantiles: adjusted p-value 0.25, absolute shrunken log2FC 0.75, expression 0.50.
8. Candidate acceptance is exactly `P_g AND E_g AND X_g`.
9. Recurrent rows are retained directly and do not receive fabricated TRUE candidate criteria.
10. Final panel is `F = R union S union N`, deduplicated by gene with deterministic provenance-preserving ordering.
11. Remove minimum-panel-size selection, sequential R1/R2/R3-style rule selection, fallback rules, exact-size forcing, and hard-coded historical panel comparisons from the active canonical builder.

## Schema and Output Naming Plan

Stable final-table columns to preserve:

- `gene_id`
- `direction`
- `owner_profile`
- `source_contrast`
- `source_file`
- `selection_rank`

Revised final/evidence fields:

- `feature_class`: `recurrent`, `accepted_singleton`, `accepted_non_recurrent`
- `selection_basis`
- `selection_basis_classes`
- `candidate_pool_type`
- `marker_evidence_stratum`
- `marker_source_class`
- `marker_evidence_provenance`
- `retained_marker_list_count`
- `retained_marker_list_total`
- `minimum_adjusted_p_value`
- `median_absolute_shrunken_log2fc`
- `median_base_mean`
- `passes_statistical_evidence`
- `passes_effect_magnitude`
- `passes_expression_evidence`
- `passes_candidate_acceptance`

Obsolete fields to remove from active schema:

- `selected_rule_key`
- `selected_sequential_rule`
- any minimum-panel-size, fallback-rule, or selected-rule metadata

Recommended active audit outputs:

- `ranked_marker_source_panel_candidate_pool_evidence.tsv`
- `ranked_marker_source_panel_empirical_quantile_thresholds.tsv`
- `ranked_marker_source_panel_candidate_acceptance.tsv`
- `ranked_marker_source_panel_selected_evidence_rows.tsv`
- `ranked_marker_source_panel_by_cohort.tsv`
- `ranked_marker_source_panel_by_marker_source_class.tsv`
- `ranked_marker_source_panel_by_feature_class.tsv`
- `ranked_marker_source_panel_validation.tsv`
- `ranked_marker_source_panel_run_manifest.tsv`
- `ranked_marker_source_panel_active_directory_manifest.tsv`

Legacy active outputs to remove or replace:

- `*_sensitivity_summary.tsv`
- `*_selected_marker_source_class_rows.tsv`
- `*_removed_vs_previous177.tsv`
- `*_added_vs_previous177.tsv`
- `*_overlap_vs_previous177.tsv`
- `*_removed_vs_previous125.tsv`
- `*_added_vs_previous125.tsv`
- `*_overlap_vs_previous125.tsv`

## Files Expected To Change

- `PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md`
- `scripts/build_pan_cancer_features.py`
- `config/config.yaml`
- `Snakefile`
- `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`
- `rules/pan_cancer_enrichment_marker_framework.smk`
- `scripts/build_ranked_marker_source_panel_enrichment_queries.py`
- `scripts/package_ranked_marker_source_panel_outputs.py`
- `scripts/finalize_ranked_marker_source_panel_reports.py`
- `scripts/build_gprofiler_query_sets_from_pan_cancer_features.py` if stale comments/path sentinels are present
- `scripts/build_pan_cancer_expression_matrix.R` if stale method comments are present
- `scripts/score_tumour_cellline_mapping.R` if stale method comments are present
- `scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R` if stale labels/comments require updates
- `README.md` and active docs if they describe the old active method
- focused tests or new fixtures for pan-cancer feature construction

## Implementation Completed So Far

- `scripts/build_pan_cancer_features.py`
  - Removed active predefined rule parsing/evaluation, sequential rule selection, fallback rule, minimum unique gene count, forced exact-size panel, previous-panel comparison inputs, and obsolete rule-selection export products.
  - Added fixed `RECURRENCE_MIN_COUNT = 2`.
  - Revised recurrence classes to `recurrent`, `singleton_candidate`, and `non_recurrent_candidate` with non-recurrent defined by `retained_marker_list_count == 1`.
  - Added contextual empirical thresholds for singleton/non-recurrent candidates using adjusted p-value 0.25, shrunken LFC 0.75, and expression 0.50.
  - Added all-three candidate acceptance fields and direct recurrent retention.
  - Reworked selected rows, gene evidence, final feature collapse, build summary, validation, run manifest, and export names around `selection_basis` and `feature_class`.
- `config/config.yaml`
  - Removed active recurrence threshold, predefined rules, panel selection, forced exact size, and previous-panel path keys from `marker_postprocessing.pan_cancer`.
  - Replaced `alternative_expression_quantile` with `expression_quantile: 0.50`.
- `Snakefile`
  - Removed active pan-cancer rule-selection config wiring and obsolete builder CLI arguments.
  - Declared revised audit outputs.
  - Added archive-before-rsync logic for the active feature-space directory.
  - Removed the dead config reader for the old pan-cancer recurrence key while preserving the separate enrichment-sidecar recurrence threshold at 2.

## Checks Run So Far

- `PYTHONPYCACHEPREFIX=/private/tmp/codex_pycache python -m py_compile scripts/build_pan_cancer_features.py`: PASS.
- `python -m py_compile scripts/build_pan_cancer_features.py` without pycache redirect: failed because the sandbox could not write to `scripts/__pycache__`.
- YAML parse using `python -c 'import yaml'`: blocked because `PyYAML` is not installed in the current Python.
- `snakemake --dry-run ...`: blocked because `snakemake` is not on the current PATH.
- Stale-term search still reports active matches in README, `scripts/finalize_ranked_marker_source_panel_reports.py`, `scripts/build_ranked_marker_source_panel_enrichment_queries.py`, `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`, and `scripts/build_enrichment_query_sets.R`; archive-only matches also exist under `scripts/archive/`.

## Files Expected Only To Be Rerun

- Canonical retained DESeq2 marker manifests and retained marker tables should be inspected but not modified unless an interface defect is found.
- Expression matrix, UMAP, ranking, graph, community, enrichment, package, and report outputs should be regenerated after staged panel validation. Generated scientific outputs must not be hand-edited.

## Unresolved Questions

- Whether existing modified tracked files contain intentional user edits that conflict with the requested redesign.
- Whether `rules/pan_cancer_enrichment_marker_framework.smk` should be kept and rewritten or retired from the active DAG. Current evidence favours keeping and rewriting because it is included in the active Snakefile and depends on the current feature table.
- Whether `ranked_marker_source_panel` remains acceptable as an output prefix/package name while `graph_derived_pan_cancer_feature_selection_v1_revised` is the canonical revised Version 1 method identifier.
- Whether any active consumer requires old `DEG_SET` filenames as external stable artifacts.
- Whether active generated output directories can be safely archived and replaced in this sandbox/session.
- Which environment should be used for YAML parsing and Snakemake dry-run, since the current PATH lacks `snakemake` and current Python lacks `yaml`.

## Exact Next Action

Rewrite active downstream stale consumers in this order: `scripts/build_ranked_marker_source_panel_enrichment_queries.py`, `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`, `scripts/package_ranked_marker_source_panel_outputs.py`, `scripts/finalize_ranked_marker_source_panel_reports.py`, then README/stale comments and focused tests.

## Final Session Update - 2026-07-11T00:44:29Z

### Work completed

- Core builder updated to the revised method:
  - recurrence invariant is `RECURRENCE_MIN_COUNT = 2`;
  - recurrence is counted within `cancer_type x marker_evidence_stratum`;
  - active within-context classes are `recurrent`, `singleton_candidate`, and `non_recurrent_candidate`;
  - non-recurrent candidates require `retained_marker_list_count == 1`;
  - singleton candidates require `retained_marker_list_total == 1` and `retained_marker_list_count == 1`;
  - recurrent rows are retained directly;
  - singleton and non-recurrent candidates are accepted only when statistical, effect-magnitude, and expression criteria all pass;
  - final feature table is the deterministic one-row-per-gene collapse of selected `R union S union N`.
- Active pan-cancer config simplified to the revised method and quantiles:
  - adjusted p-value quantile 0.25;
  - absolute shrunken log2FC quantile 0.75;
  - expression quantile 0.50.
- Snakefile pan-cancer feature rule updated:
  - obsolete rule-selection config readers and CLI wiring removed;
  - revised candidate-pool, threshold, acceptance, selected-evidence, class-summary, validation, run-manifest, and active-directory-manifest outputs declared;
  - archive-before-rsync logic added for the active feature-space directory.
- Downstream active consumers updated or retired:
  - enrichment query builders now use revised classes `recurrent`, `accepted_singleton`, `accepted_non_recurrent`, plus final panel;
  - package and report finalisers now derive feature counts/classes/thresholds/validation from current outputs rather than constants;
  - UMAP script terminology changed from DEG-set naming to pan-cancer marker-panel naming;
  - legacy `scripts/plot_pan_cancer_umap_deg_set.R` retired to `scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set`.
- Focused method tests added in `tests/test_pan_cancer_features.py`.
- Active stale-code search completed after implementation.

### Files modified by this task

- `PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md`
- `README.md`
- `Snakefile`
- `config/config.yaml`
- `scripts/build_enrichment_query_sets.R`
- `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`
- `scripts/build_pan_cancer_features.py`
- `scripts/build_ranked_marker_source_panel_enrichment_queries.py`
- `scripts/export_ranked_marker_source_panel_matrices.R`
- `scripts/finalize_ranked_marker_source_panel_reports.py`
- `scripts/package_ranked_marker_source_panel_outputs.py`
- `scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R`
- `tests/test_pan_cancer_features.py`

### Files deleted, moved, or retired

- Deleted from active scripts tree:
  - `scripts/plot_pan_cancer_umap_deg_set.R`
- Retired archive copy:
  - `scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set`

### Outputs archived or replaced

- No active scientific output tree was replaced.
- No active feature-space directory was changed.
- No archive was created, because the real staged panel build could not run without the canonical contrast-level marker manifests.

### Validation and test results

- Python syntax compilation: PASS with pycache redirected to `/private/tmp/codex_pycache`.
- R parse checks: PASS for modified R scripts.
- YAML parse validation: PASS using Ruby YAML, because the current Python environment lacks `PyYAML`.
- Focused fixture tests: PASS when invoked directly through Python, because `pytest` is not installed.
- Targeted Snakemake dry-run for `construct_pan_cancer_feature_panel`: PASS using:
  - `/Users/eltonugbogu/miniforge3/envs/snakemake/bin/snakemake --runtime-source-cache-path /private/tmp/snakemake_runtime_cache --dry-run --cores 1 construct_pan_cancer_feature_panel --config pipeline_profile=pan_cancer`
- Builder staged real-data run: BLOCKED before producing outputs because the canonical manifests are absent:
  - `results/unsupervised/brca/deseq2_markers/markers/contrast_level_marker_manifest.tsv`
  - `results/unsupervised/nbl/deseq2_markers/markers/contrast_level_marker_manifest.tsv`
  - `results/unsupervised/rbl/deseq2_markers/markers/contrast_level_marker_manifest.tsv`
- Downstream dependency dry-run: BLOCKED by missing VST RDS inputs:
  - `/data/brca/brca_vst_joint.rds`
  - `/data/nbl/nbl_vst_joint.rds`
  - `/data/rbl/rbl_vst_joint.rds`
  - `/data/heme/heme_joint_expr_vst_batch_corrected.rds`
- Active stale-term search excluding `scripts/archive/**` and `archive/**`: PASS, no matches for the obsolete method terms requested in the mission.
- Active `379` search excluding archives: no panel-size constants. Remaining active hits are unrelated library identifiers or checksum fragments in supplementary data-description files.
- Archive-only stale matches remain in:
  - `scripts/archive/build_pan_cancer_features.py.bak_20260503_063432`
  - `scripts/archive/build_pan_cancer_features.py.bak_20260503_064528`
  - `scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set`

### Current blockers

- The checkout does not contain the canonical BRCA/NBL/RBL `contrast_level_marker_manifest.tsv` files required for the real staged pan-cancer builder run.
- The checkout does not contain the `/data/...` VST RDS files required for downstream expression/UMAP/ranking/graph dry-run expansion and reruns.
- Because of those missing inputs, final panel size, R/S/N counts, active output archival, active output replacement, and downstream regenerated outputs could not be completed in this session.

### Exact next action

Generate or restore the canonical profile marker manifests, then run the staged builder without touching the active feature-space directory:

```bash
/Users/eltonugbogu/miniforge3/envs/snakemake/bin/snakemake \
  --runtime-source-cache-path /private/tmp/snakemake_runtime_cache \
  --cores 1 \
  construct_pan_cancer_feature_panel \
  --config pipeline_profile=pan_cancer
```

If the intent is to avoid replacing the active feature-space tree on that first run, invoke `scripts/build_pan_cancer_features.py` directly into a staging directory using the three restored manifest paths and the revised CLI.

### Workflow jobs launched

- No non-dry-run Snakemake jobs were launched.
- No upstream DESeq2 jobs were rerun.
- No downstream scientific outputs were regenerated.

### Current branch and HEAD

- Branch: `version2-probability-graph`
- HEAD: `13592660ac3eea3b62d5b31b053fbd9eb3decee8`

### Current git status

```text
## version2-probability-graph...origin/version2-probability-graph
 M .gitignore
 M README.md
 M Snakefile
 M config/config.yaml
 M docs/version1/PIPELINE_DOCUMENTATION.md
 M scripts/build_enrichment_query_sets.R
 M scripts/build_enrichment_summary_top_terms.R
 M scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
 M scripts/build_pan_cancer_cell_line_similarity_graph.R
 M scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
 M scripts/build_pan_cancer_features.py
 M scripts/build_pan_cancer_graph.R
 M scripts/build_ranked_marker_source_panel_enrichment_queries.py
 M scripts/compute_cell_line_louvain_resolution_sweep.R
 M scripts/compute_pan_cancer_cell_line_communities.R
 M scripts/compute_pan_cancer_cell_line_validation.R
 M scripts/compute_pan_cancer_communities.R
 M scripts/deseq2_component_vs_rest.R
 M scripts/deseq2_isolate_degs.R
 M scripts/export_ranked_marker_source_panel_matrices.R
 M scripts/finalize_ranked_marker_source_panel_reports.py
 M scripts/package_ranked_marker_source_panel_outputs.py
 M scripts/parse_gprofiler_export_to_enrichment_summary.py
 M scripts/plot_deg_eval_pan_cancer.R
 M scripts/plot_enrichment_top_terms_heatmap.R
 M scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
 M scripts/plot_pan_cancer_two_panel.R
 D scripts/plot_pan_cancer_umap_deg_set.R
 M scripts/prepare_deseq2_inputs.R
 M scripts/run_gprofiler_from_manifest.R
 M scripts/select_gprofiler_terms.py
 M scripts/validate_deseq2_staged_inputs.R
 M scripts/write_deseq2_directional_marker_tables.R
?? .cursor/
?? .github/
?? AGENT.md
?? PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
?? PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
?? SRP409177.runinfo.csv
?? archive/
?? docs/version1/DGE_REPORT.md
?? docs/version1/ecdf.md
?? docs/version2/version2_implementation_status.md
?? nbl_21_fastq_fast_check.tsv
?? nbl_missing_valid_count_samples.txt
?? preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt
?? scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
?? scripts/derive_resolved_graph_node_statistics.py
?? supplementary/
?? sync_logs/
?? test_logs/
?? tests/
```

### Current git diff --stat

```text
 .gitignore                                         |   19 -
 README.md                                          |   31 +-
 Snakefile                                          | 1318 ++++-----
 config/config.yaml                                 |  203 +-
 docs/version1/PIPELINE_DOCUMENTATION.md            |   73 +-
 scripts/build_enrichment_query_sets.R              |   32 +-
 scripts/build_enrichment_summary_top_terms.R       |   14 +-
 ...profiler_query_sets_from_pan_cancer_features.py |    7 +-
 .../build_pan_cancer_cell_line_similarity_graph.R  |  123 +-
 ...ancer_enrichment_marker_framework_query_sets.py |  117 +-
 scripts/build_pan_cancer_features.py               | 3009 +++++++++-----------
 scripts/build_pan_cancer_graph.R                   |  190 +-
 ...anked_marker_source_panel_enrichment_queries.py |   44 +-
 .../compute_cell_line_louvain_resolution_sweep.R   |   21 +-
 scripts/compute_pan_cancer_cell_line_communities.R |   85 +-
 scripts/compute_pan_cancer_cell_line_validation.R  |   50 +-
 scripts/compute_pan_cancer_communities.R           |   35 +-
 scripts/deseq2_component_vs_rest.R                 |  761 ++---
 scripts/deseq2_isolate_degs.R                      | 1259 ++++----
 .../export_ranked_marker_source_panel_matrices.R   |   10 +-
 .../finalize_ranked_marker_source_panel_reports.py |  690 ++---
 .../package_ranked_marker_source_panel_outputs.py  |   60 +-
 ...parse_gprofiler_export_to_enrichment_summary.py |   15 +-
 scripts/plot_deg_eval_pan_cancer.R                 |   87 +-
 scripts/plot_enrichment_top_terms_heatmap.R        |  140 +-
 ...ot_pan_cancer_tumour_cell_line_alignment_umap.R |   74 +-
 scripts/plot_pan_cancer_two_panel.R                |   58 +-
 scripts/plot_pan_cancer_umap_deg_set.R             |  776 -----
 scripts/prepare_deseq2_inputs.R                    |  483 +++-
 scripts/run_gprofiler_from_manifest.R              |  475 +--
 scripts/select_gprofiler_terms.py                  |   14 +-
 scripts/validate_deseq2_staged_inputs.R            |  511 ++--
 scripts/write_deseq2_directional_marker_tables.R   |  401 +--
 33 files changed, 5091 insertions(+), 6094 deletions(-)
```

### Notes on unrelated modifications

The repository was dirty before this implementation began. Several modified files shown in the current git status and diff stat pre-date this task or were not part of this implementation. They were not reset or discarded.

## Method Label Correction - 2026-07-11T07:57:20Z

- The implemented aggregation method is labelled as the revised Version 1 method, not Version 2.
- Canonical active method identifier: `graph_derived_pan_cancer_feature_selection_v1_revised`.
- The current Git branch name `version2-probability-graph` is branch metadata only. It must not be used as evidence that the pan-cancer marker aggregation implemented here is a Version 2 or probability-based method.
- Updated active config, Snakefile default, builder/package/report constants, README text, Version 1 docs, and this handoff to use the revised Version 1 label.
- Verified Version 2 wording separately: remaining active `version2`, `Version 2`, `probability graph`, and `v2` wording belongs to explicit Version 2 development workflow files such as `Snakefile.v2`, `config/config_v2.yaml`, `docs/version2/`, and Version 2 helper scripts, or to the branch-name/status text in this handoff.
- Validation after relabelling:
  - Python syntax compilation for touched Python scripts: PASS.
  - R parse check for the touched matrix-export script: PASS.
  - Targeted Snakemake dry-run for `construct_pan_cancer_feature_panel`: PASS.

## Dirty-State Resolution Snapshot - 2026-07-11T08:05:21Z

### `git branch --show-current`

```text
version2-probability-graph
```

### `git rev-parse HEAD`

```text
13592660ac3eea3b62d5b31b053fbd9eb3decee8
```

### `git status --short`

```text
 M .gitignore
 M README.md
 M Snakefile
 M config/config.yaml
 M docs/version1/PIPELINE_DOCUMENTATION.md
 M scripts/build_enrichment_query_sets.R
 M scripts/build_enrichment_summary_top_terms.R
 M scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
 M scripts/build_pan_cancer_cell_line_similarity_graph.R
 M scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
 M scripts/build_pan_cancer_features.py
 M scripts/build_pan_cancer_graph.R
 M scripts/build_ranked_marker_source_panel_enrichment_queries.py
 M scripts/compute_cell_line_louvain_resolution_sweep.R
 M scripts/compute_pan_cancer_cell_line_communities.R
 M scripts/compute_pan_cancer_cell_line_validation.R
 M scripts/compute_pan_cancer_communities.R
 M scripts/deseq2_component_vs_rest.R
 M scripts/deseq2_isolate_degs.R
 M scripts/export_ranked_marker_source_panel_matrices.R
 M scripts/finalize_ranked_marker_source_panel_reports.py
 M scripts/package_ranked_marker_source_panel_outputs.py
 M scripts/parse_gprofiler_export_to_enrichment_summary.py
 M scripts/plot_deg_eval_pan_cancer.R
 M scripts/plot_enrichment_top_terms_heatmap.R
 M scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
 M scripts/plot_pan_cancer_two_panel.R
 D scripts/plot_pan_cancer_umap_deg_set.R
 M scripts/prepare_deseq2_inputs.R
 M scripts/run_gprofiler_from_manifest.R
 M scripts/select_gprofiler_terms.py
 M scripts/validate_deseq2_staged_inputs.R
 M scripts/write_deseq2_directional_marker_tables.R
?? .cursor/
?? .github/
?? AGENT.md
?? PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
?? PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
?? SRP409177.runinfo.csv
?? archive/
?? docs/version1/DGE_REPORT.md
?? docs/version1/ecdf.md
?? docs/version2/version2_implementation_status.md
?? nbl_21_fastq_fast_check.tsv
?? nbl_missing_valid_count_samples.txt
?? preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt
?? scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
?? scripts/derive_resolved_graph_node_statistics.py
?? supplementary/
?? sync_logs/
?? test_logs/
?? tests/
```

### `git diff --stat`

```text
 .gitignore                                         |   19 -
 README.md                                          |   31 +-
 Snakefile                                          | 1318 ++++-----
 config/config.yaml                                 |  203 +-
 docs/version1/PIPELINE_DOCUMENTATION.md            |   73 +-
 scripts/build_enrichment_query_sets.R              |   32 +-
 scripts/build_enrichment_summary_top_terms.R       |   14 +-
 ...profiler_query_sets_from_pan_cancer_features.py |    7 +-
 .../build_pan_cancer_cell_line_similarity_graph.R  |  123 +-
 ...ancer_enrichment_marker_framework_query_sets.py |  117 +-
 scripts/build_pan_cancer_features.py               | 3009 +++++++++-----------
 scripts/build_pan_cancer_graph.R                   |  190 +-
 ...anked_marker_source_panel_enrichment_queries.py |   44 +-
 .../compute_cell_line_louvain_resolution_sweep.R   |   21 +-
 scripts/compute_pan_cancer_cell_line_communities.R |   85 +-
 scripts/compute_pan_cancer_cell_line_validation.R  |   50 +-
 scripts/compute_pan_cancer_communities.R           |   35 +-
 scripts/deseq2_component_vs_rest.R                 |  761 ++---
 scripts/deseq2_isolate_degs.R                      | 1259 ++++----
 .../export_ranked_marker_source_panel_matrices.R   |   10 +-
 .../finalize_ranked_marker_source_panel_reports.py |  690 ++---
 .../package_ranked_marker_source_panel_outputs.py  |   60 +-
 ...parse_gprofiler_export_to_enrichment_summary.py |   15 +-
 scripts/plot_deg_eval_pan_cancer.R                 |   87 +-
 scripts/plot_enrichment_top_terms_heatmap.R        |  140 +-
 ...ot_pan_cancer_tumour_cell_line_alignment_umap.R |   74 +-
 scripts/plot_pan_cancer_two_panel.R                |   58 +-
 scripts/plot_pan_cancer_umap_deg_set.R             |  776 -----
 scripts/prepare_deseq2_inputs.R                    |  483 +++-
 scripts/run_gprofiler_from_manifest.R              |  475 +--
 scripts/select_gprofiler_terms.py                  |   14 +-
 scripts/validate_deseq2_staged_inputs.R            |  511 ++--
 scripts/write_deseq2_directional_marker_tables.R   |  401 +--
 33 files changed, 5091 insertions(+), 6094 deletions(-)
```

### `git diff --name-status`

```text
M	.gitignore
M	README.md
M	Snakefile
M	config/config.yaml
M	docs/version1/PIPELINE_DOCUMENTATION.md
M	scripts/build_enrichment_query_sets.R
M	scripts/build_enrichment_summary_top_terms.R
M	scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
M	scripts/build_pan_cancer_cell_line_similarity_graph.R
M	scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
M	scripts/build_pan_cancer_features.py
M	scripts/build_pan_cancer_graph.R
M	scripts/build_ranked_marker_source_panel_enrichment_queries.py
M	scripts/compute_cell_line_louvain_resolution_sweep.R
M	scripts/compute_pan_cancer_cell_line_communities.R
M	scripts/compute_pan_cancer_cell_line_validation.R
M	scripts/compute_pan_cancer_communities.R
M	scripts/deseq2_component_vs_rest.R
M	scripts/deseq2_isolate_degs.R
M	scripts/export_ranked_marker_source_panel_matrices.R
M	scripts/finalize_ranked_marker_source_panel_reports.py
M	scripts/package_ranked_marker_source_panel_outputs.py
M	scripts/parse_gprofiler_export_to_enrichment_summary.py
M	scripts/plot_deg_eval_pan_cancer.R
M	scripts/plot_enrichment_top_terms_heatmap.R
M	scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
M	scripts/plot_pan_cancer_two_panel.R
D	scripts/plot_pan_cancer_umap_deg_set.R
M	scripts/prepare_deseq2_inputs.R
M	scripts/run_gprofiler_from_manifest.R
M	scripts/select_gprofiler_terms.py
M	scripts/validate_deseq2_staged_inputs.R
M	scripts/write_deseq2_directional_marker_tables.R
```

### `git ls-files --others --exclude-standard`

```text
.cursor/projects/1769865452490/agent-transcripts/1a994437-e884-4296-ba37-5e96a18f85c4.txt
.cursor/projects/1769865452490/agent-transcripts/3bc656c0-2e16-456b-a993-f049bdb87c87.txt
.cursor/projects/1769865452490/agent-transcripts/66fc9769-aa44-47d5-bd59-603265e20cf1.txt
.cursor/projects/1770100137482/agent-transcripts/9c868c70-5045-4e9a-9d36-b8b66ba80616.txt
.cursor/projects/home-chu25-TCGA-hema-pipeline/terminals/3.txt
.cursor/projects/home-chu25-dsmz-bcc-analysis-unsupervised-main-scripts-consensusclust/terminals/1.txt
.cursor/projects/home-chu25-dsmz-bcc-analysis-unsupervised/terminals/1.txt
.cursor/projects/home-chu25-dsmz-dsmz-rbl-analysis-rbl-processing/terminals/1.txt
.cursor/projects/home-chu25-dsmz-leu-lym-dsmz-analysis/terminals/1.txt
.cursor/projects/home-chu25-dsmz-pipeline-pipeline-unsupervised/mcps/cursor-browser-extension/INSTRUCTIONS.md
.cursor/projects/home-chu25-dsmz-pipeline-pipeline-unsupervised/mcps/cursor-ide-browser/INSTRUCTIONS.md
.cursor/projects/home-chu25-dsmz-pipeline/terminals/2.txt
.cursor/projects/home-chu25-dsmz-pipeline/terminals/6.txt
.cursor/projects/home-chu25-dsmz-pipeline/terminals/7.txt
.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/agent-tools/53317590-b9b2-4477-ac97-b07a16a453bc.txt
.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/mcps/cursor-browser-extension/INSTRUCTIONS.md
.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/mcps/cursor-ide-browser/INSTRUCTIONS.md
.cursor/projects/tmp-c77eb14d-8aa5-48e8-ab49-fb918d645372/mcps/cursor-browser-extension/INSTRUCTIONS.md
.cursor/projects/tmp-c77eb14d-8aa5-48e8-ab49-fb918d645372/mcps/cursor-ide-browser/INSTRUCTIONS.md
.cursor/skills-cursor/create-rule/SKILL.md
.cursor/skills-cursor/create-skill/SKILL.md
.cursor/skills-cursor/create-subagent/SKILL.md
.cursor/skills-cursor/migrate-to-skills/SKILL.md
.cursor/skills-cursor/update-cursor-settings/SKILL.md
.github/workflows/deploy-pages.yml
AGENT.md
PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
SRP409177.runinfo.csv
archive/Snakefile.bak_20260503_062313
archive/Snakefile.bak_20260503_063432
archive/Snakefile.bak_20260503_064528
archive/Snakefile.bak_20260503_064759
archive/agent_passes_202605/NAMING_PASS_REPORT.md
archive/agent_passes_202605/NAMING_PASS_diff_Snakefile.patch
archive/agent_passes_202605/NAMING_PASS_diff_multicohort_communities.patch
archive/agent_passes_202605/NAMING_PASS_diff_plot_consensus.patch
archive/agent_passes_202605/NAMING_PASS_diff_publication.patch
archive/agent_passes_202605/NAMING_PASS_diff_scan_ambiguous_rbl.patch
archive/agent_passes_202605/NAMING_PASS_diff_visualize_resolved.patch
archive/agent_passes_202605/PATIENT_REF_PASS_REPORT.md
archive/agent_passes_202605/PATIENT_REF_PASS_diff_Snakefile.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_add_component_to_metadata.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_compute_and_plot_multicohort_cancer_communities.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_plot_consensus_graph.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_plot_publication_cell_line_similarity_and_resolved_networks.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_prepare_deseq2_inputs.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_scan_ambiguous_rbl_cellline_ids.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_verify_brca_pipeline_results.patch
archive/agent_passes_202605/PATIENT_REF_PASS_diff_visualize_resolved_dsmz_graph.patch
archive/agent_passes_202605/PLOTTING_PASS_REPORT.md
archive/agent_passes_202605/PLOTTING_PASS_diff_consensus.patch
archive/agent_passes_202605/PLOTTING_PASS_diff_resolved.patch
archive/agent_passes_202605/REPRODUCIBILITY_PASS_REPORT.md
archive/agent_passes_202605/REPRODUCIBILITY_PASS_diff_Snakefile.patch
archive/agent_passes_202605/REPRODUCIBILITY_PASS_diff_smk_yaml.patch
archive/agent_passes_202605/semantic_correctness_audit_report.md
archive/dsmz_cellline_graph_node_stats.tsv
archive/snakefile_backups_202605/Snakefile.bak_anchor_audit_20260527_1442
archive/snakefile_backups_202605/Snakefile.bak_anchor_schema_20260527_1513
archive/snakefile_backups_202605/Snakefile.bak_naming_20260527_1629
archive/snakefile_backups_202605/Snakefile.bak_naming_20260528
archive/snakefile_backups_202605/Snakefile.bak_patient_ref_20260528
archive/snakefile_backups_202605/Snakefile.bak_rename_20260528
archive/snakefile_backups_202605/Snakefile.bak_reproducibility_20260530
archive/validation/01_permutation_test_neighbourhood.R
archive/validation/02_random_baseline_comparison.R
archive/validation/03_silhouette_report.R
archive/validation/04_model_selection_summary.R
docs/version1/DGE_REPORT.md
docs/version1/ecdf.md
docs/version2/version2_implementation_status.md
nbl_21_fastq_fast_check.tsv
nbl_missing_valid_count_samples.txt
preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt
scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
scripts/derive_resolved_graph_node_statistics.py
supplementary/data_description/README.md
supplementary/data_description/analysis_input_files.tsv
supplementary/data_description/cell_line_sample_ids_by_cohort.tsv
supplementary/data_description/data_availability_statement.md
supplementary/data_description/dsmz_167_cell_line_ids.tsv
supplementary/data_description/dsmz_167_cell_line_metadata.tsv
supplementary/data_description/dsmz_rbl_cell_line_subset.tsv
supplementary/data_description/excluded_or_not_available_sources.tsv
supplementary/data_description/processing_summary.md
supplementary/data_description/reproducibility_checksums.tsv
supplementary/data_description/scripts/export_cell_line_ids_by_cohort.R
supplementary/data_description/scripts/export_dsmz_cell_line_ids.R
supplementary/data_description/scripts/export_tumour_purity_passed_ids.R
supplementary/data_description/tumour_purity_passed_sample_ids.tsv
sync_logs/jarvis_hpc_mac_sync_20260627_125015.txt
sync_logs/jarvis_hpc_mac_sync_20260627_224654.txt
sync_logs/jarvis_hpc_mac_sync_20260628_034822.txt
sync_logs/jarvis_hpc_mac_sync_20260628_042608.txt
sync_logs/jarvis_hpc_mac_sync_20260704_170747.txt
sync_logs/jarvis_hpc_mac_sync_20260705_061842.txt
sync_logs/jarvis_hpc_mac_sync_20260706_233722.txt
test_logs/profile_refactor_test_20260623_141241.txt
tests/test_pan_cancer_features.py
```

## Dirty-State Classification - 2026-07-11

Basis for classification:

- Pre-task dirty state was read from the `Complete git status at start` block in this handoff.
- Current state was captured above using the requested git commands.
- Tracked files were inspected with `git diff --numstat` and targeted file diffs.
- Untracked single files were inspected with line counts, file type checks, and header reads.
- Untracked directories were inspected with `du` and `find` inventories.
- `git status --ignored --short` was checked before making any cache/generated-file decision.
- No files were deleted, moved, staged, committed, reset, cleaned, or discarded during this classification pass.

| file | status | category | reason | keep/commit/archive/ignore/remove | action taken |
|---|---:|---|---|---|---|
| `.gitignore` | M | pre_existing_unrelated_change | Dirty before task; direct diff shows removal of ignore patterns for local editor state, logs, archives, runinfo, and supplementary scratch folders. This explains much of the current untracked noise. | user review; do not commit with pan-cancer implementation unless approved | none |
| `README.md` | M | documentation_change | Dirty before task; current diff includes revised V1 method label and selection-route text plus pre-existing workflow wording edits. | hunk-stage revised V1 documentation only after approval | none |
| `Snakefile` | M | snakefile_change | Dirty before task; current diff is mixed, including intended pan-cancer feature-rule rewiring and pre-existing graph/community/output changes. | hunk-stage revised V1 pan-cancer rule/config/output sections only after approval | none |
| `config/config.yaml` | M | config_change | Dirty before task; current diff is mixed, including intended revised V1 pan-cancer config and pre-existing DESeq2/profile/graph config changes. | hunk-stage revised V1 pan-cancer config only after approval | none |
| `docs/version1/PIPELINE_DOCUMENTATION.md` | M | documentation_change | Dirty before task; inspected diff includes revised Version 1 method wording plus pre-existing DGE documentation expansion. | hunk-stage revised V1 wording only after approval | none |
| `scripts/build_enrichment_query_sets.R` | M | downstream_schema_change | Dirty before task; inspected diff includes removal of obsolete strict-vs-operative/isolate-rescued query route and background/report wording updates. | commit with downstream schema changes if approved | none |
| `scripts/build_enrichment_summary_top_terms.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff preserves query-manifest metadata and appears to pre-date this task. | leave uncommitted pending user review | none |
| `scripts/build_gprofiler_query_sets_from_pan_cancer_features.py` | M | pre_existing_unrelated_change | Dirty before task; inspected diff is comment/docstring clarification only, not required by revised aggregation schema. | leave uncommitted pending user review | none |
| `scripts/build_pan_cancer_cell_line_similarity_graph.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff changes graph validation and cancer_type metadata handling, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py` | M | downstream_schema_change | Dirty before task; inspected diff replaces obsolete recurrent/core and isolate-extension/rescued classes with revised feature classes. | commit with downstream schema changes if approved | none |
| `scripts/build_pan_cancer_features.py` | M | core_aggregation_change | Dirty before task; inspected diff implements revised Version 1 recurrence, candidate thresholding, all-three acceptance, selected evidence, validation, and schema. | commit as core implementation after review; use care because file was already dirty | none |
| `scripts/build_pan_cancer_graph.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff changes simplified kNN graph validation/component handling, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/build_ranked_marker_source_panel_enrichment_queries.py` | M | downstream_schema_change | Dirty before task; inspected diff updates method label, removes 379 assertion, uses candidate-pool evidence, and revised classes. | commit with downstream schema changes if approved | none |
| `scripts/compute_cell_line_louvain_resolution_sweep.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff changes lineage/cancer_type metadata handling, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/compute_pan_cancer_cell_line_communities.R` | M | pre_existing_unrelated_change | Dirty before task; current diff belongs to cell-line community handling, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/compute_pan_cancer_cell_line_validation.R` | M | pre_existing_unrelated_change | Dirty before task; current diff belongs to cell-line graph/community validation, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/compute_pan_cancer_communities.R` | M | pre_existing_unrelated_change | Dirty before task; current diff belongs to graph/community calculations, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/deseq2_component_vs_rest.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff rewrites upstream DESeq2 component marker production. User explicitly said not to alter upstream producers unless needed; this should not be included in aggregation commit. | leave uncommitted pending user review | none |
| `scripts/deseq2_isolate_degs.R` | M | pre_existing_unrelated_change | Dirty before task; current diff is upstream DESeq2 marker producer work, not part of aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/export_ranked_marker_source_panel_matrices.R` | M | downstream_schema_change | Not dirty in pre-task status; inspected diff removes hard-coded 379 and updates method label. | commit with downstream schema changes if approved | none |
| `scripts/finalize_ranked_marker_source_panel_reports.py` | M | downstream_schema_change | Dirty before task; inspected diff replaces stale hard-coded report generation with dynamic revised V1 report generation. File mode changed from executable to non-executable and should be reviewed before commit. | commit with downstream schema changes if approved; review executable bit | none |
| `scripts/package_ranked_marker_source_panel_outputs.py` | M | downstream_schema_change | Dirty before task; inspected diff derives panel counts dynamically, updates method label, and updates UMAP filenames. | commit with downstream schema changes if approved | none |
| `scripts/parse_gprofiler_export_to_enrichment_summary.py` | M | pre_existing_unrelated_change | Dirty before task; inspected diff is parser/reporting comment and manifest help cleanup, not essential to aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/plot_deg_eval_pan_cancer.R` | M | pre_existing_unrelated_change | Dirty before task; current diff relates DE evaluation plotting and shrunken LFC fields, outside the aggregation builder commit. | leave uncommitted pending user review | none |
| `scripts/plot_enrichment_top_terms_heatmap.R` | M | downstream_schema_change | Dirty before task; inspected diff removes obsolete isolate-extension/recurrent-core heatmap categories and maps current marker-source-panel query metadata. | commit with downstream schema changes only if approved; otherwise user review because pre-task dirty | none |
| `scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R` | M | downstream_schema_change | Not dirty in pre-task status; inspected diff removes DEG_SET/--deg_set active terminology and uses marker-panel feature list naming. | commit with terminology/downstream changes if approved | none |
| `scripts/plot_pan_cancer_two_panel.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff changes lineage labels to cancer_type labels for network figures, outside aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/plot_pan_cancer_umap_deg_set.R` | D | archive_retirement | Not dirty in pre-task status; intentionally retired from active scripts tree during aggregation terminology cleanup. | commit deletion with archive copy if approved | none |
| `scripts/prepare_deseq2_inputs.R` | M | pre_existing_unrelated_change | Dirty before task; upstream DESeq2 input preparation changes, not part of aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/run_gprofiler_from_manifest.R` | M | pre_existing_unrelated_change | Dirty before task; current diff belongs to live g:Profiler runner behaviour, not required for aggregation builder commit. | leave uncommitted pending user review | none |
| `scripts/select_gprofiler_terms.py` | M | pre_existing_unrelated_change | Dirty before task; inspected diff is legacy selector comment cleanup, not required by aggregation schema. | leave uncommitted pending user review | none |
| `scripts/validate_deseq2_staged_inputs.R` | M | pre_existing_unrelated_change | Dirty before task; upstream DESeq2 validation changes, not part of aggregation redesign. | leave uncommitted pending user review | none |
| `scripts/write_deseq2_directional_marker_tables.R` | M | pre_existing_unrelated_change | Dirty before task; inspected diff rewrites upstream directional marker collation. User requested this remain separate from pan-cancer recurrence. | leave uncommitted pending user review | none |
| `.cursor/projects/1769865452490/agent-transcripts/1a994437-e884-4296-ba37-5e96a18f85c4.txt` | ?? | generated_cache | Pre-existing untracked Cursor agent transcript exposed by `.gitignore` change; direct inventory under `.cursor`. | ignore or restore `.gitignore` pattern after approval | none |
| `.cursor/projects/1769865452490/agent-transcripts/3bc656c0-2e16-456b-a993-f049bdb87c87.txt` | ?? | generated_cache | Pre-existing untracked Cursor agent transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/1769865452490/agent-transcripts/66fc9769-aa44-47d5-bd59-603265e20cf1.txt` | ?? | generated_cache | Pre-existing untracked Cursor agent transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/1770100137482/agent-transcripts/9c868c70-5045-4e9a-9d36-b8b66ba80616.txt` | ?? | generated_cache | Pre-existing untracked Cursor agent transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-TCGA-hema-pipeline/terminals/3.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-bcc-analysis-unsupervised-main-scripts-consensusclust/terminals/1.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-bcc-analysis-unsupervised/terminals/1.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-dsmz-rbl-analysis-rbl-processing/terminals/1.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-leu-lym-dsmz-analysis/terminals/1.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-pipeline-pipeline-unsupervised/mcps/cursor-browser-extension/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-pipeline-pipeline-unsupervised/mcps/cursor-ide-browser/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-pipeline/terminals/2.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-pipeline/terminals/6.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/home-chu25-dsmz-pipeline/terminals/7.txt` | ?? | generated_cache | Pre-existing untracked Cursor terminal transcript exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/agent-tools/53317590-b9b2-4477-ac97-b07a16a453bc.txt` | ?? | generated_cache | Pre-existing local Cursor agent-tool output exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/mcps/cursor-browser-extension/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/tmp-379c2f06-fa50-4885-baab-92ee628b99ab/mcps/cursor-ide-browser/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/tmp-c77eb14d-8aa5-48e8-ab49-fb918d645372/mcps/cursor-browser-extension/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/projects/tmp-c77eb14d-8aa5-48e8-ab49-fb918d645372/mcps/cursor-ide-browser/INSTRUCTIONS.md` | ?? | generated_cache | Pre-existing local Cursor MCP instruction file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/skills-cursor/create-rule/SKILL.md` | ?? | generated_cache | Pre-existing local Cursor skill file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/skills-cursor/create-skill/SKILL.md` | ?? | generated_cache | Pre-existing local Cursor skill file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/skills-cursor/create-subagent/SKILL.md` | ?? | generated_cache | Pre-existing local Cursor skill file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/skills-cursor/migrate-to-skills/SKILL.md` | ?? | generated_cache | Pre-existing local Cursor skill file exposed by `.gitignore` change. | ignore after approval | none |
| `.cursor/skills-cursor/update-cursor-settings/SKILL.md` | ?? | generated_cache | Pre-existing local Cursor skill file exposed by `.gitignore` change. | ignore after approval | none |
| `.github/workflows/deploy-pages.yml` | ?? | unknown_requires_user_review | Pre-existing untracked GitHub Pages workflow; inspected header shows deployment config unrelated to aggregation. | user review; do not commit with aggregation | none |
| `AGENT.md` | ?? | unknown_requires_user_review | Pre-existing untracked pipeline execution contract, 834 lines; unrelated to aggregation update. | user review; do not commit with aggregation | none |
| `PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md` | ?? | handoff_file | Pre-existing untracked audit handoff; inspected header records prior audit, not the implementation handoff. | archive or commit separately if desired | none |
| `PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md` | ?? | handoff_file | Created/updated by this task; contains implementation and dirty-state handoff. | commit with documentation/test commit or archive, pending approval | none |
| `SRP409177.runinfo.csv` | ?? | generated_cache | Pre-existing SRA runinfo CSV exposed by `.gitignore` change; inspected as data/scratch input listing. | ignore or archive separately; do not commit with aggregation | none |
| `archive/Snakefile.bak_20260503_062313` | ?? | pre_existing_unrelated_change | Pre-existing archive backup exposed by `.gitignore` change. | ignore or commit archive separately after review | none |
| `archive/Snakefile.bak_20260503_063432` | ?? | pre_existing_unrelated_change | Pre-existing archive backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/Snakefile.bak_20260503_064528` | ?? | pre_existing_unrelated_change | Pre-existing archive backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/Snakefile.bak_20260503_064759` | ?? | pre_existing_unrelated_change | Pre-existing archive backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_REPORT.md` | ?? | pre_existing_unrelated_change | Pre-existing archive report exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_Snakefile.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_multicohort_communities.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_plot_consensus.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_publication.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_scan_ambiguous_rbl.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/NAMING_PASS_diff_visualize_resolved.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_REPORT.md` | ?? | pre_existing_unrelated_change | Pre-existing archive report exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_Snakefile.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_add_component_to_metadata.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_compute_and_plot_multicohort_cancer_communities.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_plot_consensus_graph.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_plot_publication_cell_line_similarity_and_resolved_networks.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_prepare_deseq2_inputs.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_scan_ambiguous_rbl_cellline_ids.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_verify_brca_pipeline_results.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PATIENT_REF_PASS_diff_visualize_resolved_dsmz_graph.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PLOTTING_PASS_REPORT.md` | ?? | pre_existing_unrelated_change | Pre-existing archive report exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PLOTTING_PASS_diff_consensus.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/PLOTTING_PASS_diff_resolved.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/REPRODUCIBILITY_PASS_REPORT.md` | ?? | pre_existing_unrelated_change | Pre-existing archive report exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/REPRODUCIBILITY_PASS_diff_Snakefile.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/REPRODUCIBILITY_PASS_diff_smk_yaml.patch` | ?? | pre_existing_unrelated_change | Pre-existing archive patch exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/agent_passes_202605/semantic_correctness_audit_report.md` | ?? | pre_existing_unrelated_change | Pre-existing archive report exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/dsmz_cellline_graph_node_stats.tsv` | ?? | pre_existing_unrelated_change | Pre-existing archived TSV exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_anchor_audit_20260527_1442` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_anchor_schema_20260527_1513` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_naming_20260527_1629` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_naming_20260528` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_patient_ref_20260528` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_rename_20260528` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/snakefile_backups_202605/Snakefile.bak_reproducibility_20260530` | ?? | pre_existing_unrelated_change | Pre-existing Snakefile backup exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/validation/01_permutation_test_neighbourhood.R` | ?? | pre_existing_unrelated_change | Pre-existing archived validation script exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/validation/02_random_baseline_comparison.R` | ?? | pre_existing_unrelated_change | Pre-existing archived validation script exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/validation/03_silhouette_report.R` | ?? | pre_existing_unrelated_change | Pre-existing archived validation script exposed by `.gitignore` change. | ignore/archive after review | none |
| `archive/validation/04_model_selection_summary.R` | ?? | pre_existing_unrelated_change | Pre-existing archived validation script exposed by `.gitignore` change. | ignore/archive after review | none |
| `docs/version1/DGE_REPORT.md` | ?? | documentation_change | Pre-existing untracked file; inspected header and revised method section. It now contains revised Version 1 wording but also untracked pre-existing report content. | user review; commit only if entire new report should enter repo | none |
| `docs/version1/ecdf.md` | ?? | pre_existing_unrelated_change | Pre-existing untracked thesis/figure text; inspected header, unrelated to aggregation implementation. | user review; do not commit with aggregation | none |
| `docs/version2/version2_implementation_status.md` | ?? | pre_existing_unrelated_change | Pre-existing untracked Version 2 development status; inspected header, explicitly separate from revised Version 1 aggregation. | user review; do not commit with aggregation | none |
| `nbl_21_fastq_fast_check.tsv` | ?? | generated_cache | Pre-existing NBL FASTQ check table exposed by `.gitignore` change. | ignore or archive separately | none |
| `nbl_missing_valid_count_samples.txt` | ?? | generated_cache | Pre-existing NBL sample list exposed by `.gitignore` change. | ignore or archive separately | none |
| `preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt` | ?? | generated_cache | Pre-existing RBL Snakemake output log exposed by `.gitignore` change. | ignore or archive separately | none |
| `scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set` | ?? | archive_retirement | Created by this aggregation update as the preserved copy of the retired legacy DEG-set UMAP script. | commit with retirement deletion if approved | none |
| `scripts/derive_resolved_graph_node_statistics.py` | ?? | unknown_requires_user_review | Pre-existing untracked script; inspected header shows resolved graph node statistics, unrelated to aggregation update. | user review; do not commit with aggregation | none |
| `supplementary/data_description/README.md` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description output exposed by `.gitignore` change. | user review; likely separate supplementary commit | none |
| `supplementary/data_description/analysis_input_files.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/cell_line_sample_ids_by_cohort.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/data_availability_statement.md` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description markdown exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/dsmz_167_cell_line_ids.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/dsmz_167_cell_line_metadata.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/dsmz_rbl_cell_line_subset.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/excluded_or_not_available_sources.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/processing_summary.md` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description markdown exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/reproducibility_checksums.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary checksums exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/scripts/export_cell_line_ids_by_cohort.R` | ?? | pre_existing_unrelated_change | Pre-existing supplementary export script exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/scripts/export_dsmz_cell_line_ids.R` | ?? | pre_existing_unrelated_change | Pre-existing supplementary export script exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/scripts/export_tumour_purity_passed_ids.R` | ?? | pre_existing_unrelated_change | Pre-existing supplementary export script exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `supplementary/data_description/tumour_purity_passed_sample_ids.tsv` | ?? | pre_existing_unrelated_change | Pre-existing supplementary data-description TSV exposed by `.gitignore` change. | user review; separate from aggregation | none |
| `sync_logs/jarvis_hpc_mac_sync_20260627_125015.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260627_224654.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260628_034822.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260628_042608.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260704_170747.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260705_061842.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `sync_logs/jarvis_hpc_mac_sync_20260706_233722.txt` | ?? | generated_cache | Pre-existing sync log exposed by `.gitignore` change. | ignore after approval | none |
| `test_logs/profile_refactor_test_20260623_141241.txt` | ?? | generated_cache | Pre-existing test log exposed by `.gitignore` change. | ignore after approval | none |
| `tests/test_pan_cancer_features.py` | ?? | test_change | Created by this aggregation update; inspected file contains focused recurrence, threshold-boundary, acceptance, direct-retention, and deduplication tests. | commit with revised V1 test/documentation commit if approved | none |

Ignored generated/cache paths seen in `git status --ignored --short` include `.DS_Store`, `.claude/`, `.codex_remote_edits/`, `.local_review/`, `.snakemake/`, `scripts/__pycache__/`, `data/`, `derived_results/`, `results/`, `reports/`, `reference_data/`, and `supplementary_data/`. These are not in normal `git status --short`. No ignored path was removed during this pass.

## Proposed Cleanup Plan Before Execution

Recommended plan: **PLAN C**.

Reason: the worktree was already broadly dirty before the aggregation implementation. Several files that now contain intended revised Version 1 aggregation changes were already modified before the task (`Snakefile`, `config/config.yaml`, `README.md`, `scripts/build_pan_cancer_features.py`, report/enrichment scripts). Whole-file commits would mix the revised aggregation implementation with unrelated pre-existing user changes.

Proposed execution after approval:

1. Preserve current state without deleting anything:
   - leave generated/cache untracked files in place;
   - do not use `git reset --hard`;
   - do not use `git clean -fd`;
   - do not alter `.gitignore` unless explicitly approved.
2. Create a safety WIP branch or preservation commit for the current tracked dirty state only, excluding untracked generated/cache directories and large scratch outputs.
3. Create a clean implementation branch from `13592660ac3eea3b62d5b31b053fbd9eb3decee8`.
4. Cherry-pick or apply only approved revised Version 1 aggregation hunks into focused commits:
   - Commit 1: pan-cancer config and Snakefile rule/output wiring.
   - Commit 2: `scripts/build_pan_cancer_features.py` recurrence plus all-three candidate acceptance.
   - Commit 3: downstream enrichment/report/schema consumers.
   - Commit 4: revised Version 1 terminology, documentation, handoff, and tests.
   - Commit 5: retire `scripts/plot_pan_cancer_umap_deg_set.R` with preserved archive copy.
5. Before each commit, run and record:
   - `git diff --cached --stat`
   - `git diff --cached --name-status`
6. Leave all files classified as `pre_existing_unrelated_change`, `generated_cache`, or `unknown_requires_user_review` uncommitted unless the user explicitly approves including them.

No cleanup plan has been executed yet.

## Safety Branch Snapshot - 2026-07-11

Created branch: `safety/wip-pan-cancer-aggregation-dirty-state-20260711`.

### `git branch --show-current`

```text
safety/wip-pan-cancer-aggregation-dirty-state-20260711
```

### `git rev-parse HEAD`

```text
13592660ac3eea3b62d5b31b053fbd9eb3decee8
```

### `git status --short`

```text
 M .gitignore
 M README.md
 M Snakefile
 M config/config.yaml
 M docs/version1/PIPELINE_DOCUMENTATION.md
 M scripts/build_enrichment_query_sets.R
 M scripts/build_enrichment_summary_top_terms.R
 M scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
 M scripts/build_pan_cancer_cell_line_similarity_graph.R
 M scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
 M scripts/build_pan_cancer_features.py
 M scripts/build_pan_cancer_graph.R
 M scripts/build_ranked_marker_source_panel_enrichment_queries.py
 M scripts/compute_cell_line_louvain_resolution_sweep.R
 M scripts/compute_pan_cancer_cell_line_communities.R
 M scripts/compute_pan_cancer_cell_line_validation.R
 M scripts/compute_pan_cancer_communities.R
 M scripts/deseq2_component_vs_rest.R
 M scripts/deseq2_isolate_degs.R
 M scripts/export_ranked_marker_source_panel_matrices.R
 M scripts/finalize_ranked_marker_source_panel_reports.py
 M scripts/package_ranked_marker_source_panel_outputs.py
 M scripts/parse_gprofiler_export_to_enrichment_summary.py
 M scripts/plot_deg_eval_pan_cancer.R
 M scripts/plot_enrichment_top_terms_heatmap.R
 M scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
 M scripts/plot_pan_cancer_two_panel.R
 D scripts/plot_pan_cancer_umap_deg_set.R
 M scripts/prepare_deseq2_inputs.R
 M scripts/run_gprofiler_from_manifest.R
 M scripts/select_gprofiler_terms.py
 M scripts/validate_deseq2_staged_inputs.R
 M scripts/write_deseq2_directional_marker_tables.R
?? .cursor/
?? .github/
?? AGENT.md
?? PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
?? PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
?? SRP409177.runinfo.csv
?? archive/
?? docs/version1/DGE_REPORT.md
?? docs/version1/ecdf.md
?? docs/version2/version2_implementation_status.md
?? nbl_21_fastq_fast_check.tsv
?? nbl_missing_valid_count_samples.txt
?? preprocessing_and_quality_control/rbl/rbl_snakefile_output.txt
?? scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
?? scripts/derive_resolved_graph_node_statistics.py
?? supplementary/
?? sync_logs/
?? test_logs/
?? tests/
```

### `git diff --stat`

```text
 .gitignore                                         |   19 -
 README.md                                          |   31 +-
 Snakefile                                          | 1318 ++++-----
 config/config.yaml                                 |  203 +-
 docs/version1/PIPELINE_DOCUMENTATION.md            |   73 +-
 scripts/build_enrichment_query_sets.R              |   32 +-
 scripts/build_enrichment_summary_top_terms.R       |   14 +-
 ...profiler_query_sets_from_pan_cancer_features.py |    7 +-
 .../build_pan_cancer_cell_line_similarity_graph.R  |  123 +-
 ...ancer_enrichment_marker_framework_query_sets.py |  117 +-
 scripts/build_pan_cancer_features.py               | 3009 +++++++++-----------
 scripts/build_pan_cancer_graph.R                   |  190 +-
 ...anked_marker_source_panel_enrichment_queries.py |   44 +-
 .../compute_cell_line_louvain_resolution_sweep.R   |   21 +-
 scripts/compute_pan_cancer_cell_line_communities.R |   85 +-
 scripts/compute_pan_cancer_cell_line_validation.R  |   50 +-
 scripts/compute_pan_cancer_communities.R           |   35 +-
 scripts/deseq2_component_vs_rest.R                 |  761 ++---
 scripts/deseq2_isolate_degs.R                      | 1259 ++++----
 .../export_ranked_marker_source_panel_matrices.R   |   10 +-
 .../finalize_ranked_marker_source_panel_reports.py |  690 ++---
 .../package_ranked_marker_source_panel_outputs.py  |   60 +-
 ...parse_gprofiler_export_to_enrichment_summary.py |   15 +-
 scripts/plot_deg_eval_pan_cancer.R                 |   87 +-
 scripts/plot_enrichment_top_terms_heatmap.R        |  140 +-
 ...ot_pan_cancer_tumour_cell_line_alignment_umap.R |   74 +-
 scripts/plot_pan_cancer_two_panel.R                |   58 +-
 scripts/plot_pan_cancer_umap_deg_set.R             |  776 -----
 scripts/prepare_deseq2_inputs.R                    |  483 +++-
 scripts/run_gprofiler_from_manifest.R              |  475 +--
 scripts/select_gprofiler_terms.py                  |   14 +-
 scripts/validate_deseq2_staged_inputs.R            |  511 ++--
 scripts/write_deseq2_directional_marker_tables.R   |  401 +--
 33 files changed, 5091 insertions(+), 6094 deletions(-)
```

### `git diff --name-status`

```text
M	.gitignore
M	README.md
M	Snakefile
M	config/config.yaml
M	docs/version1/PIPELINE_DOCUMENTATION.md
M	scripts/build_enrichment_query_sets.R
M	scripts/build_enrichment_summary_top_terms.R
M	scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
M	scripts/build_pan_cancer_cell_line_similarity_graph.R
M	scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
M	scripts/build_pan_cancer_features.py
M	scripts/build_pan_cancer_graph.R
M	scripts/build_ranked_marker_source_panel_enrichment_queries.py
M	scripts/compute_cell_line_louvain_resolution_sweep.R
M	scripts/compute_pan_cancer_cell_line_communities.R
M	scripts/compute_pan_cancer_cell_line_validation.R
M	scripts/compute_pan_cancer_communities.R
M	scripts/deseq2_component_vs_rest.R
M	scripts/deseq2_isolate_degs.R
M	scripts/export_ranked_marker_source_panel_matrices.R
M	scripts/finalize_ranked_marker_source_panel_reports.py
M	scripts/package_ranked_marker_source_panel_outputs.py
M	scripts/parse_gprofiler_export_to_enrichment_summary.py
M	scripts/plot_deg_eval_pan_cancer.R
M	scripts/plot_enrichment_top_terms_heatmap.R
M	scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
M	scripts/plot_pan_cancer_two_panel.R
D	scripts/plot_pan_cancer_umap_deg_set.R
M	scripts/prepare_deseq2_inputs.R
M	scripts/run_gprofiler_from_manifest.R
M	scripts/select_gprofiler_terms.py
M	scripts/validate_deseq2_staged_inputs.R
M	scripts/write_deseq2_directional_marker_tables.R
```

### `git ls-files --others --exclude-standard`

The untracked set matches the dirty-state snapshot above, including local Cursor files, untracked archives, handoff/audit files, NBL/RBL scratch files, supplementary data-description files, sync/test logs, the retired UMAP archive copy, and `tests/test_pan_cancer_features.py`. Generated/cache/scratch paths remain uncommitted unless explicitly staged below.

## Safety WIP Staged Files - 2026-07-11

Staging commands used:

```bash
git add -u
git add PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md tests/test_pan_cancer_features.py scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
```

No generated/cache/scratch trees were staged.

### `git diff --cached --stat`

```text
 .gitignore                                         |   19 -
 PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md           | 1043 +++++++
 README.md                                          |   31 +-
 Snakefile                                          | 1318 ++++-----
 config/config.yaml                                 |  203 +-
 docs/version1/PIPELINE_DOCUMENTATION.md            |   73 +-
 ...n_cancer_umap_deg_set.R.retired_legacy_deg_set} |    4 +-
 scripts/build_enrichment_query_sets.R              |   32 +-
 scripts/build_enrichment_summary_top_terms.R       |   14 +-
 ...profiler_query_sets_from_pan_cancer_features.py |    7 +-
 .../build_pan_cancer_cell_line_similarity_graph.R  |  123 +-
 ...ancer_enrichment_marker_framework_query_sets.py |  117 +-
 scripts/build_pan_cancer_features.py               | 3009 +++++++++-----------
 scripts/build_pan_cancer_graph.R                   |  190 +-
 ...anked_marker_source_panel_enrichment_queries.py |   44 +-
 .../compute_cell_line_louvain_resolution_sweep.R   |   21 +-
 scripts/compute_pan_cancer_cell_line_communities.R |   85 +-
 scripts/compute_pan_cancer_cell_line_validation.R  |   50 +-
 scripts/compute_pan_cancer_communities.R           |   35 +-
 scripts/deseq2_component_vs_rest.R                 |  761 ++---
 scripts/deseq2_isolate_degs.R                      | 1259 ++++----
 .../export_ranked_marker_source_panel_matrices.R   |   10 +-
 .../finalize_ranked_marker_source_panel_reports.py |  690 ++---
 .../package_ranked_marker_source_panel_outputs.py  |   60 +-
 ...parse_gprofiler_export_to_enrichment_summary.py |   15 +-
 scripts/plot_deg_eval_pan_cancer.R                 |   87 +-
 scripts/plot_enrichment_top_terms_heatmap.R        |  140 +-
 ...ot_pan_cancer_tumour_cell_line_alignment_umap.R |   74 +-
 scripts/plot_pan_cancer_two_panel.R                |   58 +-
 scripts/prepare_deseq2_inputs.R                    |  483 +++-
 scripts/run_gprofiler_from_manifest.R              |  475 +--
 scripts/select_gprofiler_terms.py                  |   14 +-
 scripts/validate_deseq2_staged_inputs.R            |  511 ++--
 scripts/write_deseq2_directional_marker_tables.R   |  401 +--
 tests/test_pan_cancer_features.py                  |  175 ++
 35 files changed, 6311 insertions(+), 5320 deletions(-)
```

### `git diff --cached --name-status`

```text
M	.gitignore
A	PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
M	README.md
M	Snakefile
M	config/config.yaml
M	docs/version1/PIPELINE_DOCUMENTATION.md
R099	scripts/plot_pan_cancer_umap_deg_set.R	scripts/archive/plot_pan_cancer_umap_deg_set.R.retired_legacy_deg_set
M	scripts/build_enrichment_query_sets.R
M	scripts/build_enrichment_summary_top_terms.R
M	scripts/build_gprofiler_query_sets_from_pan_cancer_features.py
M	scripts/build_pan_cancer_cell_line_similarity_graph.R
M	scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
M	scripts/build_pan_cancer_features.py
M	scripts/build_pan_cancer_graph.R
M	scripts/build_ranked_marker_source_panel_enrichment_queries.py
M	scripts/compute_cell_line_louvain_resolution_sweep.R
M	scripts/compute_pan_cancer_cell_line_communities.R
M	scripts/compute_pan_cancer_cell_line_validation.R
M	scripts/compute_pan_cancer_communities.R
M	scripts/deseq2_component_vs_rest.R
M	scripts/deseq2_isolate_degs.R
M	scripts/export_ranked_marker_source_panel_matrices.R
M	scripts/finalize_ranked_marker_source_panel_reports.py
M	scripts/package_ranked_marker_source_panel_outputs.py
M	scripts/parse_gprofiler_export_to_enrichment_summary.py
M	scripts/plot_deg_eval_pan_cancer.R
M	scripts/plot_enrichment_top_terms_heatmap.R
M	scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
M	scripts/plot_pan_cancer_two_panel.R
M	scripts/prepare_deseq2_inputs.R
M	scripts/run_gprofiler_from_manifest.R
M	scripts/select_gprofiler_terms.py
M	scripts/validate_deseq2_staged_inputs.R
M	scripts/write_deseq2_directional_marker_tables.R
A	tests/test_pan_cancer_features.py
```

This safety WIP commit intentionally preserves broad tracked dirt; it is not the reviewed revised V1 implementation commit set.

## Clean Implementation Branch Start - 2026-07-11

Created branch from recorded base:

```bash
git switch -c feature/revised-v1-pan-cancer-aggregation 13592660ac3eea3b62d5b31b053fbd9eb3decee8
```

### Branch and HEAD

```text
branch=feature/revised-v1-pan-cancer-aggregation
HEAD=13592660ac3eea3b62d5b31b053fbd9eb3decee8
```

### `git status --short`

```text
?? .github/
?? AGENT.md
?? PAN_CANCER_AGGREGATION_AUDIT_HANDOFF.md
?? docs/version1/DGE_REPORT.md
?? docs/version1/ecdf.md
?? docs/version2/version2_implementation_status.md
?? scripts/derive_resolved_graph_node_statistics.py
```

The clean branch did not inherit the dirty `.gitignore` change. The listed untracked files are pre-existing unrelated or user-review files and remain uncommitted.

## Clean Branch Commit 1 - Config and Snakefile Wiring

Applied only revised V1 pan-cancer aggregation wiring hunks:

- `config/config.yaml`
- `Snakefile`

Excluded unrelated pre-existing config/Snakefile hunks, including DESeq2 input refactors, profile grid edits, graph/community changes, and `.gitignore` changes.

### `git diff --cached --stat`

```text
 Snakefile          | 182 ++++++++++++++++++++++++++---------------------------
 config/config.yaml |  56 ++++-------------
 2 files changed, 100 insertions(+), 138 deletions(-)
```

### `git diff --cached --name-status`

```text
M	Snakefile
M	config/config.yaml
```

### Validation

```text
ruby -e 'require "yaml"; YAML.load_file("config/config.yaml"); puts "YAML OK"'
YAML OK
```

```text
/Users/eltonugbogu/miniforge3/envs/snakemake/bin/snakemake --runtime-source-cache-path /private/tmp/snakemake_runtime_cache --dry-run --cores 1 construct_pan_cancer_feature_panel --config pipeline_profile=pan_cancer
PASS: dry-run parsed the DAG and planned ensure_profile_deseq2_markers x3 plus construct_pan_cancer_feature_panel x1.
The dry-run reported missing canonical contrast manifests as expected in this checkout.
```

## Clean Branch Commit 2 - Core Builder Refactor

Applied `scripts/build_pan_cancer_features.py` from the inspected safety WIP state. The diff was classified as core aggregation work: recurrence invariant, recurrent/singleton/non-recurrent candidate classification, contextual empirical thresholds, all-three candidate acceptance, direct recurrent retention, selected evidence rows, deterministic deduplication, revised summaries, revised validation, and removal of old sequential/minimum-size/fallback/previous-panel builder interfaces.

### `git diff --cached --stat`

```text
 scripts/build_pan_cancer_features.py | 3009 +++++++++++++++-------------------
 1 file changed, 1324 insertions(+), 1685 deletions(-)
```

### `git diff --cached --name-status`

```text
M	scripts/build_pan_cancer_features.py
```

### Validation

```text
PYTHONPYCACHEPREFIX=/private/tmp/codex_pycache python -m py_compile scripts/build_pan_cancer_features.py
PASS
```

Obsolete active builder term check:

```text
rg -n "selected_rule_key|selected_sequential_rule|minimum_unique_gene_count|predefined_rule_order|fallback_rule|force_exact_size|previous177|previous125|alternative_expression_quantile|all_three_median_expression|all_three_quantile_expression|direction_aware_two_of_three_median_expression" scripts/build_pan_cancer_features.py
PASS: no matches
```

## Clean Branch Commit 4 - Documentation, Terminology, and Tests

Applied only revised V1 documentation/test terminology changes:

- `README.md`: active panel method now labelled `graph_derived_pan_cancer_feature_selection_v1_revised`.
- `scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R`: comments and active CLI terminology now describe a graph-informed DESeq2 marker-derived feature panel, not a DEG set.
- `Snakefile` and `config/config.yaml`: active UMAP rule/config now use `PAN_CANCER_MARKER_PANEL`, `--feature_list`, and `pan_cancer_alignment` consistently with the R script.
- `tests/test_pan_cancer_features.py`: focused recurrence, threshold-boundary, all-three acceptance, direct recurrent retention, deduplication, and threshold-context fixtures.

Excluded broad pre-existing documentation and workflow changes not needed for the revised V1 aggregation extraction.

### Validation

```text
PYTHONPYCACHEPREFIX=/private/tmp/codex_pycache python -m py_compile tests/test_pan_cancer_features.py
PASS
```

```text
PYTHONPYCACHEPREFIX=/private/tmp/codex_pycache python -c "<load tests/test_pan_cancer_features.py and run test_* functions>"
tests ok
```

```text
Rscript -e "invisible(parse('scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R'))"
PASS
```

```text
ruby -e 'require "yaml"; YAML.load_file("config/config.yaml"); puts "YAML OK"'
YAML OK
```

```text
/Users/eltonugbogu/miniforge3/envs/snakemake/bin/snakemake --runtime-source-cache-path /private/tmp/snakemake_runtime_cache --dry-run --cores 1 construct_pan_cancer_feature_panel --config pipeline_profile=pan_cancer
PASS: dry-run parsed the DAG and planned ensure_profile_deseq2_markers x3 plus construct_pan_cancer_feature_panel x1.
The dry-run reported missing canonical contrast manifests as expected in this checkout.
```

### Stale-Term Notes

- Active UMAP/config/Snakefile stale terms `--deg_set`, `deg_set_alignment`, and `DEG_SET` were removed from the current alignment rule.
- Remaining `deg_set` matches are confined to `scripts/plot_pan_cancer_umap_deg_set.R`, which is scheduled for separate legacy retirement.
- Version 2 / probability-graph matches are confined to explicit Version 2 development files or a schema-version comment and are not used as the revised V1 aggregation method label.

### `git diff --cached --stat`

```text
 PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md           | 1337 ++++++++++++++++++++
 README.md                                          |    6 +-
 Snakefile                                          |    9 +-
 config/config.yaml                                 |    3 +-
 ...ot_pan_cancer_tumour_cell_line_alignment_umap.R |   74 +-
 tests/test_pan_cancer_features.py                  |  175 +++
 6 files changed, 1551 insertions(+), 53 deletions(-)
```

### `git diff --cached --name-status`

```text
A	PAN_CANCER_AGGREGATION_UPDATE_HANDOFF.md
M	README.md
M	Snakefile
M	config/config.yaml
M	scripts/plot_pan_cancer_tumour_cell_line_alignment_umap.R
A	tests/test_pan_cancer_features.py
```

## Clean Branch Commit 3 - Downstream Schema and Enrichment Consumers

Applied downstream consumer changes only:

- `scripts/build_enrichment_query_sets.R`
- `scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py`
- `scripts/build_ranked_marker_source_panel_enrichment_queries.py`
- `scripts/export_ranked_marker_source_panel_matrices.R`
- `scripts/finalize_ranked_marker_source_panel_reports.py`
- `scripts/package_ranked_marker_source_panel_outputs.py`
- `scripts/plot_enrichment_top_terms_heatmap.R`

Excluded unrelated pre-existing graph/community/upstream DESeq2/g:Profiler runner/parser changes.

### `git diff --cached --stat`

```text
 scripts/build_enrichment_query_sets.R              |  32 +-
 ...ancer_enrichment_marker_framework_query_sets.py | 117 ++--
 ...anked_marker_source_panel_enrichment_queries.py |  44 +-
 .../export_ranked_marker_source_panel_matrices.R   |  10 +-
 .../finalize_ranked_marker_source_panel_reports.py | 690 ++++++++++-----------
 .../package_ranked_marker_source_panel_outputs.py  |  60 +-
 scripts/plot_enrichment_top_terms_heatmap.R        | 144 +++--
 7 files changed, 544 insertions(+), 553 deletions(-)
```

### `git diff --cached --name-status`

```text
M	scripts/build_enrichment_query_sets.R
M	scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py
M	scripts/build_ranked_marker_source_panel_enrichment_queries.py
M	scripts/export_ranked_marker_source_panel_matrices.R
M	scripts/finalize_ranked_marker_source_panel_reports.py
M	scripts/package_ranked_marker_source_panel_outputs.py
M	scripts/plot_enrichment_top_terms_heatmap.R
```

### Validation

```text
PYTHONPYCACHEPREFIX=/private/tmp/codex_pycache python -m py_compile scripts/build_pan_cancer_enrichment_marker_framework_query_sets.py scripts/build_ranked_marker_source_panel_enrichment_queries.py scripts/package_ranked_marker_source_panel_outputs.py scripts/finalize_ranked_marker_source_panel_reports.py
PASS
```

```text
Rscript -e "invisible(parse('scripts/build_enrichment_query_sets.R')); invisible(parse('scripts/export_ranked_marker_source_panel_matrices.R')); invisible(parse('scripts/plot_enrichment_top_terms_heatmap.R'))"
PASS
```

Stale downstream term check:

```text
rg -n "Legacy|recurrent_or_core|isolate_extension|isolate_rescued|isolate_rescue|ranked_nonrecurrent_marker|Ranked nonrecurrent|Anchor singleton|relaxed_iqr|strict_iqr|selected_rule|379|ranked_marker_source_pan_cancer_panel|DEG_SET|deg_set" <commit-3-files>
PASS: no matches
```
