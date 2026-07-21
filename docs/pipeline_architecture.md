# Pipeline architecture bridge

## Purpose and review scope

This document bridges three different architectural views of the repository:

1. **Source structure** — R/Python/shell files, functions, calls and imports.
   Code Review Graph (CRG) is useful primarily at this layer.
2. **Workflow orchestration** — Snakemake rules, `include:` statements, rule
   inputs/outputs, wildcards, profile configuration, shell invocations and
   cluster settings.
3. **Scientific data flow** — matrices, resolved graph tables, contrast
   manifests, retained markers, feature panels, ranking outputs, community
   assignments and enrichment query/background files.

CRG is a navigation and preliminary change-impact index. It does not replace
inspection of the root `Snakefile`, included `.smk` modules, configuration,
file schemas, scientific methods, statistical assumptions, tests or security
boundaries. Absence of an edge in CRG is not evidence that no dependency exists.

The machine-readable companion is
`reports/agent_report/pipeline_architecture_map.tsv`. Reports are intentionally
local/ignored under the repository report policy; this tracked document carries
the durable architectural interpretation.

## Entry points and profiles

| Scope | Entry point or profile | Role |
| --- | --- | --- |
| Shared/cohort analysis | `Snakefile` | Profile-aware graph, DESeq2, pan-cancer, ranking, community and enrichment workflow. |
| NBL preprocessing | `preprocessing_and_quality_control/nbl/Snakefile` | NBL acquisition, count preparation, QC and cohort-local preprocessing. |
| NBL normalisation | `preprocessing_and_quality_control/nbl/batch_corr_and_normalisation/Snakefile` | NBL batch correction and normalisation subworkflow. |
| RBL preprocessing | `preprocessing_and_quality_control/rbl/Snakefile` | RBL acquisition, count preparation, QC and cohort-local preprocessing. |
| Local execution | `profiles/local/config.yaml` | Local Snakemake execution settings. |
| Cluster execution | `profiles/slurm/config.yaml` | SLURM executor/profile settings; review separately before cluster use. |
| Shared configuration | `config/config.yaml` | Profile paths, thresholds, feature/ranking/community/enrichment parameters. |

BRCA, NBL and RBL are active patient-referenced/DESeq2 profiles and feed the
shared pan-cancer stages. HEME is an auxiliary/comparator analysis, not another
interchangeable member of the three-profile aggregation contract.

## High-level scientific flow

```text
cohort expression and clustering evidence
  -> patient tumour neighbourhoods (HC and, where configured, k-means)
  -> within-representation p-consensus
  -> per-representation patient-referenced cell-line similarity graphs
  -> cross-representation support summaries
  -> resolved cell-line neighbours
  -> node/component/isolate/anchor definitions
  -> isolate and anchor DESeq2 contrasts
  -> retained marker tables plus canonical contrast manifests
  -> BRCA/NBL/RBL pan-cancer feature aggregation
  -> selected feature panel
       -> tumour-to-cell-line bidirectional ranking
       -> selected-feature pan-cancer cell-line graph and communities
       -> graph-derived functional-enrichment query/background manifests
  -> figures, summaries and audit artefacts
```

The selected-feature community analysis is supplemented by a full-expression
DSMZ graph and Louvain-resolution sweep. This sensitivity path must not be
confused with the primary feature-panel graph.

## Rule-to-script bridge

The table below shows principal contracts, not every auxiliary figure or audit
sidecar. Exact path construction remains in the active workflow.

| Workflow/profile | Rule | Script or shell command | Principal inputs | Principal outputs | Downstream consumers | Scientific purpose |
| --- | --- | --- | --- | --- | --- | --- |
| Cohort | `tumour_nh_hc` / `tumour_nh_km` | `scripts/compute_tumour_neighbourhoods.R` | Consensus cluster RDS files and profile config | Direction-local neighbourhood tables and sentinels | `tumour_nh_consensus` | Rank patient tumours around each cell line in each representation. |
| Cohort | `tumour_nh_consensus` | `scripts/tumour_neighbourhood_p_consensus.R` | HC/KM neighbourhood outputs | Final per-direction consensus RDS/TSV | `cell_line_similarity_graph`, `summarize_p_consensus_all` | Aggregate patient-neighbour support within a representation. |
| Cohort | `cell_line_similarity_graph` | `scripts/compute_cell_line_similarity.R` | Per-direction consensus RDS | Similarity matrix/pairs, graph edges, nodes and community sidecars | Graph resolution and support-network rules | Convert patient-referenced evidence into a cell-line graph per representation. |
| Cohort | `summarize_p_consensus_all` | `scripts/summarize_p_consensus_all.R` | All direction consensus RDS files and configured threshold | Direction summaries, winner maps, rankings and weights | `resolve_dsmz_graph_neighbours` | Select globally and locally supported representations. |
| Cohort | `resolve_dsmz_graph_neighbours` | `scripts/resolve_dsmz_graph_neighbours.R` | Winners, direction summary and all direction edge TSVs | `resolved_dsmz_neighbours.tsv` | Node statistics, graph plots and DESeq2 grouping | Intersect global-best neighbours with the union of tied local-winner neighbours. |
| Cohort | `derive_resolved_graph_node_statistics` | **Missing:** `scripts/derive_resolved_graph_node_statistics.py` | Resolved-neighbour TSV | Edge, node-statistics, short-name and anchor-audit TSVs | Plotting and isolate/component/anchor derivation | Intended to materialise graph-node analytical contracts. Currently blocked; see “Verified contract defects.” |
| Cohort | `validate_deseq2_inputs` / `prepare_deseq2_inputs` | `scripts/validate_deseq2_staged_inputs.R`; `scripts/prepare_deseq2_inputs.R` | Raw counts/metadata plus graph definitions | Validation manifest, cell-line-only counts and metadata | Graph-derived contrasts | Validate integer-count/sample contracts and stage the cell-line dataset. |
| Cohort | `derive_isolate_list` / `derive_components_list` / `derive_anchor_list` | Embedded shell/awk | Node-statistics and anchor-audit TSVs | Isolate, component, anchor and anchor-component definitions | DESeq2 contrast rule | Convert graph topology/centrality into focal/reference groups. |
| Cohort | `run_isolate_and_anchor_deseq2_contrasts` | `scripts/deseq2_isolate_degs.R` | Cell-line counts/metadata and graph-derived definitions | Contrast result tables, retained markers, QC and a declared canonical manifest | Directional tables, recurrence, feature aggregation | Run isolate-versus-rest and anchor-versus-outside-component Wald tests. |
| Cohort | `write_deseq2_directional_marker_tables` | `scripts/write_deseq2_directional_marker_tables.R` | Canonical contrast manifest and retained/result tables | Canonical directional marker tables | Profile readiness and shared aggregation | Standardise contrast evidence for shared consumers. |
| Shared | `construct_pan_cancer_feature_panel` | `scripts/build_pan_cancer_features.py` | BRCA/NBL/RBL canonical contrast manifests and retained markers | Final feature panel/lists plus evidence, threshold and provenance audits | Expression matrix, ranking, graph and enrichment | Select recurrent and empirically accepted marker evidence. |
| Shared | `build_pan_cancer_feature_expression_matrix` | `scripts/build_pan_cancer_feature_expression_matrix.R` | Profile VST matrices/metadata and feature list | Joint selected-feature expression matrix and metadata | Alignment and ranking | Create a common tumour/cell-line representation. |
| Shared | `score_tumour_cellline_mapping` | `scripts/score_tumour_cellline_mapping.R` | Selected features, joint expression and metadata | Bidirectional score/rank matrices, top-k/MRR and confidence summaries | Ranking plots and diagnostics | Rank cell-line groups against patient tumours (configured Spearman; cosine sensitivity). |
| Shared | `build_pan_cancer_cell_line_similarity_graph` | `scripts/build_pan_cancer_cell_line_similarity_graph.R` | Cell-line-only selected-feature expression and metadata | Weighted graph, edge table and node metadata | Community rules | Build the primary pan-cancer cell-line graph. |
| Shared | `compute_pan_cancer_cell_line_communities` | `scripts/compute_pan_cancer_cell_line_communities.R` | Selected-feature graph and metadata | Louvain/Leiden assignments, metrics and discordant profiles | Sweeps, validation and plots | Characterise lineage/community concordance. |
| Shared sensitivity | Leiden/Louvain sweep rules | Community script; `scripts/compute_cell_line_louvain_resolution_sweep.R` | Selected-feature graph | Resolution assignments, summaries and diagnostics | Comparative reporting | Test community stability across resolutions. |
| Shared sensitivity | Full-expression object/graph/sweep rules | `scripts/build_dsmz_joint_expression_cell_line_object.R`; graph builder; Louvain sweep | Full VST expression | Full-expression graph and sweep outputs | Comparative reporting | Test whether community conclusions depend on feature selection. |
| Shared enrichment | `construct_graph_derived_enrichment_queries` | `scripts/construct_graph_derived_enrichment_queries.py` | Feature panel, eligible background and profile contrast manifests | Query/background files, query manifest and integrity/recurrence tables | `run_gprofiler_query_batch` | Construct explicit query/universe contracts. |
| Shared enrichment | `run_gprofiler_query_batch` | `scripts/run_gprofiler_from_manifest.R` | Query manifest and query/background files | Primary/IEA term tables, fingerprints and provenance | `aggregate_gprofiler_results` | Call g:Profiler with a custom background. This is a networked scientific stage. |
| Shared enrichment | Aggregate/report/plot rules | Aggregate, summary and plotting R scripts | Per-query g:Profiler corpus and query manifest | Significant-term, sensitivity, support and figure artefacts | Reporting | Aggregate, filter and visualise enrichment evidence. |

## Key data contracts

### Resolved-neighbour table

`scripts/resolve_dsmz_graph_neighbours.R` writes one row per canonical cell
line, including isolates. Required semantic fields are:

- `cell_line` — canonical node identifier;
- `final_neighbors` — semicolon-delimited resolved neighbours, empty for an isolate;
- `n_final` — resolved degree; zero defines an isolate;
- `reps_used` — human-readable global/local evidence description;
- `best_overall_dir` and `winner_dir` — provenance for the intersection rule.

The cell-line universe is the union of winner-map IDs and all edge endpoints,
so a degree-zero node must remain as a row rather than disappearing.

### Node statistics and anchor audit

The intended node-statistics stage should produce component/isolate fields and
an anchor audit containing at least node/component identifiers, degree,
normalised and unnormalised betweenness, degree selection, bridge selection,
`anchor_selected` and `anchor_selection_reason`. The implementation currently
referenced by the rule is absent; the similar logic in
`scripts/visualize_resolved_dsmz_graph.py` is not a verified substitute for the
missing active command.

### Canonical contrast manifest

Shared feature aggregation requires these 20 columns:

```text
cancer_type, contrast_id, contrast_type, marker_evidence_stratum,
focal_profile_id, focal_component_id, reference_definition,
marker_table_path, marker_gene_list_path, result_table_path,
n_result_genes, n_pvalue_nonmissing, n_padj_nonmissing,
n_significant_before_effect_filter, n_markers_before_cap,
n_markers_after_cap, adjusted_p_value_threshold, minimum_base_mean,
minimum_absolute_shrunken_log2fc, maximum_markers_per_contrast
```

`contrast_id` must be unique. Allowed contrast types are
`isolate_focal_vs_other_same_cancer` and
`anchor_focal_vs_outside_focal_component`; anchor rows require a component ID.
Paths are resolved relative to the manifest contract, not inferred by CRG.

### Retained marker tables

The feature builder requires one row per unique cleaned Ensembl gene ID and at
least:

```text
gene_id, baseMean, wald_statistic, p_value, adjusted_p_value,
log2_fold_change_unshrunk, log2_fold_change_shrunken,
log2_fold_change_posterior_sd, absolute_shrunken_log2_fold_change,
effect_direction, contrast_marker_rank
```

Version stripping must not create ID collisions. Marker gene lists and result
tables referenced by the manifest are part of the same contract.

### Feature panel

The final panel uses unique canonical `gene_id` values, deterministic unique
`selection_rank` values, non-empty marker-evidence provenance, recurrence and
candidate-pool classification, effect-direction consistency and attributed
cancer-type/source fields. Recurrent evidence requires at least two retained
marker-list occurrences in the active implementation. Singleton/non-recurrent
candidates are not unconditional “rescues”: all statistical, effect-magnitude
and expression-evidence criteria must pass.

### Ranking outputs

The ranking stage writes both tumour-to-cell-line and cell-line-to-tumour
contracts, including score/rank tables, group summaries, top-k consistency,
MRR/metric tables, group/replicate mapping and low-confidence cases. Consumers
must preserve canonical tumour, cell-line group, lineage, score, rank and metric
identifiers. Euclidean distance is used elsewhere for representation/alignment;
it is not the primary ranking metric in this rule.

### Community assignments

Primary community outputs identify `profile_id`, collapsed `cell_line`,
`cancer_type`, community ID, community size, majority cancer type, purity,
singleton status, degree and weighted degree. Resolution-sweep outputs add a
resolution value and stability/degeneracy summaries. The selected-feature and
full-expression sources must stay distinguishable.

### Enrichment query manifest and backgrounds

The g:Profiler runner validates at least:

```text
query_id, query_family, gene_count, gene_list_path,
background_count, background_path, background_sha256,
ordered, ranking_used_for_enrichment, ranked_genes_path,
query_execution_status
```

Every query must be a subset of its explicit custom background; counts and the
background checksum must match. Ordered queries require a ranked-gene file.
Primary and IEA-sensitivity results use fingerprints/provenance and must not be
silently interchanged.

## Cohort-specific versus shared ownership

- **BRCA:** root-profile inputs and BRCA preprocessing helpers; contributes one
  canonical contrast manifest to shared aggregation.
- **NBL:** its own preprocessing/batch-correction workflows and tumour-purity
  analysis; contributes one canonical manifest.
- **RBL:** its own preprocessing/batch-correction workflow and tumour-purity
  analysis; contributes one canonical manifest.
- **HEME:** auxiliary comparator/extraction code. Do not make it an implicit
  fourth pan-cancer profile without an explicit scientific/workflow decision.
- **Shared:** root pan-cancer feature, ranking, cell-line network/community and
  enrichment stages. Shared rules consume profile contracts; they do not own
  cohort-specific raw-data acquisition.

## Verified contract defects at 2026-07-21

These are source/workflow findings, not CRG inferences:

1. `derive_resolved_graph_node_statistics` declares and executes
   `scripts/derive_resolved_graph_node_statistics.py`, but that tracked file is
   absent. The rule and consumers of its four outputs cannot execute from this
   checkout as written.
2. `run_isolate_and_anchor_deseq2_contrasts` declares and tests
   `markers/contrast_level_marker_manifest.tsv`, but
   `scripts/deseq2_isolate_degs.R` writes only the legacy
   `markers/marker_sets_manifest.tsv` with four columns. The downstream feature
   and enrichment readers require the canonical schema above. No compatible
   active writer was found.

These defects should be resolved through a scientifically reviewed workflow
change plus synthetic contract tests. They were not changed merely to improve
CRG metrics.

## Review procedure and known CRG limitations

For a repository-wide review:

1. Confirm the CRG graph SHA matches `HEAD`.
2. Use minimal context and a focused source/file query to locate likely R,
   Python or shell implementations.
3. Open every returned source before making a claim.
4. Inspect the corresponding rule in `Snakefile` or `rules/*.smk`.
5. Follow concrete input/output paths and inspect the producer/consumer schema.
6. Inspect relevant YAML/profile keys and SLURM settings manually.
7. Use change-impact output as a conservative hint; verify it against the Git
   diff and this bridge.

The current CRG build does not index Snakefiles, `.smk` or generic YAML; it does
not resolve wildcard paths, manifests, file-output-to-input dependencies or
cluster execution. Its detected flows are same-file call paths, not Snakemake
execution flows. A custom “Python-like” Snakemake language mapping is therefore
not retained: syntactic nodes without rule/file semantics would create false
confidence. The reliable operating model is supported source indexing plus
this maintained architecture bridge.
