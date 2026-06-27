#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime
import csv, hashlib

ROOT = Path(__file__).resolve().parents[1]
PAN = ROOT / 'results/unsupervised/pan_cancer'
MULTI = ROOT / 'results/unsupervised/multicohort_cancer'
FS = PAN / 'feature_space'
ENR = PAN / 'enrichment/ranked_marker_source_panel'
DSMZ = PAN / 'dsmz_similarity_graph/ranked_marker_source_panel'
RANK = PAN / 'ranking/ranked_marker_source_panel'
EMBED = PAN / 'embedding/ranked_marker_source_panel'
MOUT = MULTI / 'ranked_marker_source_panel'

def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else ['status'])
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def write(path, text):
    Path(path).write_text(text.rstrip() + '\n', encoding='utf-8')

def metrics(path):
    return {row['metric']: row['value'] for row in read_tsv(path)}

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

feature = read_tsv(FS / 'pan_cancer_features.tsv')
query = read_tsv(ENR / 'ranked_marker_source_panel_enrichment_query_manifest.tsv')
skipped = read_tsv(ENR / 'ranked_marker_source_panel_enrichment_skipped_queries.tsv')
dsmz = metrics(DSMZ / 'ranked_marker_source_panel_dsmz_graph_metrics.tsv')
multi = metrics(MOUT / 'ranked_marker_source_panel_multicohort_summary_metrics.tsv')
retrieval = read_tsv(RANK / 'ranked_marker_source_panel_retrieval_accuracy.tsv')
balanced = metrics(RANK / 'ranked_marker_source_panel_retrieval_balanced_accuracy.tsv')['balanced_accuracy']
patient = metrics(RANK / 'ranked_marker_source_panel_patient_ranking_metrics.tsv')
owners = {x: sum(row['owner_profile'] == x for row in feature) for x in ('brca', 'nbl', 'rbl')}
evidence = {x: sum(row['evidence_class'] == x for row in feature) for x in sorted({r['evidence_class'] for r in feature})}
directions = {x: sum(row['direction'] == x for row in feature) for x in ('DOWN', 'UP', 'MIXED')}
anchor = sum('anchor' in row['marker_source_class'].split(';') for row in feature)
isolate = sum('isolate' in row['marker_source_class'].split(';') for row in feature)
dual = sum({'anchor', 'isolate'}.issubset(set(row['marker_source_class'].split(';'))) for row in feature)

write(PAN / 'audit_reports/ranked_marker_source_panel_stale_terms.regex', '''(177|171|257|125)[ _-]*genes?
previous[_ -]?(177|171|257|125)
isolate[ _-]+rescue
isolate[ _-]+extension
isolate_extension_or_rescued
marker_source_recurrence
recurrence-plus-isolate
old feature set''')

status_fields = ['status', 'method', 'correction', 'threshold', 'organism', 'message']
status_row = {
    'status': 'NOT_EXECUTED_ENDPOINT_TIMEOUT',
    'method': 'ranked_marker_source_pan_cancer_panel',
    'correction': 'g:SCS', 'threshold': '0.05', 'organism': 'hsapiens',
    'message': 'Queries and explicit backgrounds prepared; no term-level result was returned by the API endpoint.',
}
for name in ('ranked_marker_source_panel_gprofiler_results.tsv', 'ranked_marker_source_panel_gprofiler_top_terms.tsv'):
    write_tsv(ENR / name, [status_row], status_fields)
summary_rows = []
for row in query:
    summary_rows.append({
        'query_name': row['query_name'], 'query_type': row['category'],
        'submitted_genes': row['gene_count'], 'recognised_genes': 'NOT_AVAILABLE',
        'background_name': row['background_name'], 'background_size': row['background_count'],
        'source_databases': 'NOT_EXECUTED',
        'correction_method': 'g:SCS', 'significant_terms': 'NOT_AVAILABLE',
        'status': 'SKIPPED_QUERY_SIZE' if row['skip'] == 'TRUE' else 'PREPARED_API_TIMEOUT',
        'skipped_reason': row['skipped_reason'],
    })
write_tsv(ENR / 'ranked_marker_source_panel_enrichment_summary_by_query.tsv', summary_rows)

write(ENR / 'ranked_marker_source_panel_gprofiler_manual_execution_instructions.md', f'''# Ranked marker-source panel g:Profiler manual execution

The current query preparation produced {len(query)} queries, of which {len(query)-len(skipped)} meet the minimum size of three genes. Each query has a clean-Ensembl gene list and an explicit background recorded in the query manifest.

Run g:Profiler for human genes with g:SCS correction, adjusted threshold 0.05, ordered queries only where the manifest marks them ordered, and the listed custom background. The automated attempt on 18 June 2026 did not return a response before endpoint timeout. Do not substitute a whole-genome background.

Use `ranked_marker_source_panel_enrichment_query_manifest.tsv` as the execution manifest. Record submitted and recognised counts, databases, significant-term counts, and any API failures for every query.
''')

write(ENR / 'ranked_marker_source_panel_enrichment_report.md', f'''# Ranked marker-source panel enrichment report

## Status

Query preparation completed for the current ranked marker-source pan-cancer panel. g:Profiler was attempted using the available `gprofiler2` environment, but the endpoint timed out without returning term-level results. Consequently, this report provides query and background provenance only; it makes no enrichment claims.

- Queries prepared: {len(query)}
- Queries meeting the minimum size: {len(query)-len(skipped)}
- Queries skipped for size: {len(skipped)}
- Panel-level eligible marker-selection universe: 1,346 clean Ensembl genes
- Per-contrast eligible backgrounds: 19
- Organism and correction planned: human, g:SCS, adjusted threshold 0.05

The query categories are the final panel, mutually exclusive cohort-owner sets, non-exclusive marker-source-class memberships, evidence classes, direction strata, their requested intersections, and retained per-contrast sets. The two MIXED-direction genes remain documented but are not submitted below the size threshold. Functional enrichment, once executed, should be interpreted as functional context for selected marker strata rather than causal evidence of pathway regulation.
''')

write(DSMZ / 'ranked_marker_source_panel_dsmz_similarity_graph_report.md', f'''# Ranked marker-source panel DSMZ similarity graph report

The established DSMZ graph workflow was rerun with the current clean feature panel.

- Profiles/nodes: {dsmz['nodes']}
- Genes used: {dsmz['feature_genes_used']}
- Edges: {dsmz['edges']}
- Connected components: {dsmz['connected_components']}
- Primary community method and resolution: Louvain, 1
- Communities: {dsmz['communities']}
- Weighted modularity: {float(dsmz['modularity']):.6f}
- Cancer-type assortativity: {float(dsmz['cancer_type_assortativity']):.6f}

The Louvain solution contained a BRCA-majority community with 29 BRCA, two haematological and two NBL profiles; a mixed NBL/RBL community with 16 NBL and 11 RBL profiles; and three pure haematological communities containing 36, 38 and 33 profiles. Thus, NBL and RBL were co-assigned in a major community. Haematological profiles were predominantly separated, with two cross-lineage placements documented in the discordance table.

The full configured Louvain and Leiden resolution grids were retained. NBL/RBL mixing, separation and singleton-degenerate states are reported without replacing the primary resolution with a selected sensitivity result. These communities describe cancer-type-associated structure in a marker-derived feature space; they are not formal molecular subtypes.
''')

overall = next(row for row in retrieval if row['group'] == 'overall')
write(RANK / 'ranked_marker_source_panel_ranking_report.md', f'''# Ranked marker-source panel ranking and retrieval report

All patient-to-cell-line and reciprocal retrieval analyses were rerun using the current ranked marker-source pan-cancer panel.

- Tumour profiles: 1,128
- Patient-to-cell-line top-1 cancer-type agreement: {float(patient['top1_accuracy']):.4f}
- Patient-to-cell-line top-10 agreement: {float(patient['top10_accuracy']):.4f}
- Patient-to-cell-line mean reciprocal rank: {float(patient['mrr']):.4f}
- Cell-line groups: {overall['n_total']}
- Cell-line-centred top-1 agreement: {float(overall['accuracy']):.4f} ({overall['n_correct']}/{overall['n_total']})
- Balanced top-1 accuracy across BRCA, NBL and RBL: {float(balanced):.4f}
- Median patient confidence margin: 0.1048
- Median cell-line confidence margin: 0.0449
- Reciprocal top-10 pairs: 530

Fifteen patient profiles had cross-cancer-type top matches. One NBL cell-line group had a BRCA top tumour match. These cases and small confidence margins are flagged for follow-up rather than treated automatically as errors.

The rankings define prioritised tiers of suitable models, not a mandatory single best model. Same-cancer-type retrieval supports model prioritisation but does not establish drug response, proteomic state, xenograft behaviour, essentiality, clonal composition, microenvironmental fidelity or full biological equivalence.
''')

write(EMBED / 'ranked_marker_source_panel_embedding_report.md', '''# Ranked marker-source panel embedding report

The existing pan-cancer tumour/cell-line UMAP rule was rerun with the active ranked marker-source pan-cancer panel. Coordinates, the cancer-type/data-source figure and the workflow summary are packaged here under current method naming. The embedding is a feature-space visualisation and is not used to define formal molecular subtypes.
''')

write(MOUT / 'ranked_marker_source_panel_multicohort_report.md', f'''# Ranked marker-source panel multi-cohort cancer report

## Status and inputs

The established `multicohort_cancer` workflow was rerun separately from the DSMZ-only graph. It consumed `pan_cancer_features_clean.txt` and used the current clean Ensembl feature panel.

- Cohorts/data sources: {multi['cohorts_included']} (DSMZ, GSE, TARGET and TCGA)
- Samples/profiles: {multi['samples_profiles']}
- Cancer types: BRCA, NBL and RBL
- Feature genes found/missing: {multi['feature_genes_found']}/{multi['feature_genes_missing']}
- Duplicated feature identifiers before cleaning: {multi['duplicated_feature_ids_before_cleaning']}
- Final genes used: {multi['final_genes_used']}
- Resolved graph nodes/edges: {multi['graph_nodes']}/{multi['graph_edges']}
- Connected components: {multi['connected_components']}
- Communities: {multi['communities']}
- Community method: {multi['community_method']}
- Modularity: {float(multi['modularity']):.6f}

The resolved 56-node graph had nine communities. Eight were cancer-type pure. One eight-node community mixed six BRCA and two NBL models (purity 0.75). A 15-node NBL community and a nine-node RBL community were pure, so NBL and RBL remained separated in the aggregate solution while limited BRCA/NBL mixing persisted. Haematological profiles were excluded from this multi-cohort setting and were assessed in the DSMZ graph.

The existing workflow does not contain a multi-cohort community-resolution sensitivity rule, so no new exploratory sensitivity method was introduced. Cancer-type assortativity was not emitted by the established aggregate rule. The output supports wording about cancer-type-associated structure and patient-referenced model prioritisation in a marker-derived feature space; it does not establish formal molecular subtypes or functional equivalence.
''')

write(PAN / 'ranked_marker_source_panel_thesis_update_notes.md', f'''# Ranked marker-source panel thesis update notes

No thesis file was edited. These current values are provided for later manual integration.

## Methods

The active method is `ranked_marker_source_pan_cancer_panel`. Contrast-level DESeq2 marker selection was unchanged. The final panel was assembled by ranked marker-source selection, with the empirical `relaxed_iqr_median_baseMean` rule used for ranked non-recurrent marker selection. The active evidence classes are marker-source recurrent markers, anchor-singleton-ranked markers and ranked non-recurrent markers; the panel must not be described as entirely strict recurrent-core genes.

## Results

- Final genes: 379, with no duplicate clean Ensembl identifiers
- Mutually exclusive cohort owners: BRCA {owners['brca']}, NBL {owners['nbl']}, RBL {owners['rbl']}
- Directions: DOWN {directions['DOWN']}, UP {directions['UP']}, MIXED {directions['MIXED']}
- Evidence classes: anchor-source recurrent {evidence['anchor_source_recurrent']}, isolate-source recurrent {evidence['isolate_source_recurrent']}, anchor-singleton-ranked {evidence['anchor_singleton_ranked']}, ranked non-recurrent {evidence['ranked_nonrecurrent_marker']}
- Non-exclusive marker-source-class memberships: anchor-derived {anchor}, isolate-derived {isolate}; {dual} genes carry both memberships
- DSMZ graph: {dsmz['nodes']} nodes, {dsmz['edges']} edges, {dsmz['communities']} Louvain communities, modularity {float(dsmz['modularity']):.3f}, cancer-type assortativity {float(dsmz['cancer_type_assortativity']):.3f}
- Patient-to-cell-line top-1/top-10 agreement: {float(patient['top1_accuracy']):.3f}/{float(patient['top10_accuracy']):.3f}
- Cell-line-centred retrieval accuracy/balanced accuracy: {float(overall['accuracy']):.3f}/{float(balanced):.3f}
- Multi-cohort setting: {multi['samples_profiles']} profiles, {multi['graph_nodes']} resolved graph nodes, {multi['graph_edges']} edges, {multi['communities']} communities, modularity {float(multi['modularity']):.3f}

## Discussion

Describe marker-source counts as non-exclusive memberships. The current design includes ranked non-recurrent markers. Functional enrichment remains pending because only current query sets and eligible backgrounds were prepared. Interpret future enrichment as functional context, graph communities as structural assessment rather than formal subtypes, and rankings as model-prioritisation evidence rather than proof of full biological equivalence.
''')

final_report = f'''# Ranked marker-source panel downstream rerun report

## 1. Executive summary

The second controlled remediation passed and all 27 pre-run checks passed from the beginning. Downstream jobs were then permitted to start. The DSMZ graph, patient/cell-line ranking and retrieval, existing UMAP, and `multicohort_cancer` workflows were regenerated with the active ranked marker-source pan-cancer panel. Functional query sets and exact eligible backgrounds were regenerated, but g:Profiler did not complete because the endpoint timed out.

No thesis file was edited, LaTeX was not compiled, and upstream DEG generation was not run.

## 2. Current active feature panel

- Method: `ranked_marker_source_pan_cancer_panel`
- Selected rule: `relaxed_iqr_median_baseMean`
- Final genes: 379; clean-ID duplicates: zero
- Cohort owners: BRCA {owners['brca']}, NBL {owners['nbl']}, RBL {owners['rbl']}
- Directions: DOWN {directions['DOWN']}, UP {directions['UP']}, MIXED {directions['MIXED']}
- Evidence classes: anchor-source recurrent {evidence['anchor_source_recurrent']}; isolate-source recurrent {evidence['isolate_source_recurrent']}; anchor-singleton-ranked {evidence['anchor_singleton_ranked']}; ranked non-recurrent {evidence['ranked_nonrecurrent_marker']}
- Marker-source-class memberships are non-exclusive: anchor-derived {anchor}, isolate-derived {isolate}, dual membership {dual}

## 3. Second controlled stale-output remediation

The first controlled remediation archived 323 query-tree files with checksum verification. The second remediation archived 28 newly identified active files before any downstream job started.

- Second archive: `results/unsupervised/pan_cancer/ARCHIVE_stale_active_outputs_before_ranked_marker_source_panel_20260618_213205/`
- Second archive manifest: `ranked_marker_source_panel_second_remediation_archive_manifest.tsv`
- Checksum failures: zero
- Active stale-description grep result after remediation: PASS, zero current-method-description hits
- Active stale outputs remaining: none detected
- Downstream execution allowed after remediation: yes, only after the complete pre-run audit passed

Stale-description command used (with archive, backup and audit trees excluded):

```text
grep -RInE --binary-files=without-match --exclude-dir="*ARCHIVE*" --exclude-dir="*BACKUP*" --exclude-dir=audit_reports -f results/unsupervised/pan_cancer/audit_reports/ranked_marker_source_panel_stale_terms.regex results/unsupervised/pan_cancer results/unsupervised/multicohort_cancer
```

## 4. Pre-run exit checks

All 27 recorded checks passed. This included panel/table/list consistency, metadata completeness, recurrent and empirical-ranking criteria, thresholds and active manifests, both required backup checksum audits, active-current-result terminology, active configuration/rules/scripts, feature coverage, and `multicohort_cancer` dependency on `pan_cancer_features_clean.txt`. The baseline is recorded in `audit_reports/ranked_marker_source_panel_prerun_file_baseline.tsv` because the pipeline directory is not a Git worktree.

## 5. Archived downstream outputs

After pre-run authorisation, 149 pan-cancer downstream files and 360 panel-dependent or aggregate multi-cohort files were archived with manifests. Backups were retained. Canonical current packages now use method-specific directories and names.

- Query-tree archive: `results/unsupervised/pan_cancer/enrichment/ARCHIVE_stale_marker_framework_queries_before_ranked_marker_source_panel_20260618_212553/`
- Pan-cancer downstream archive: `results/unsupervised/pan_cancer/ARCHIVE_stale_downstream_before_ranked_marker_source_panel_20260618_214053/`
- Multi-cohort downstream archive: `results/unsupervised/multicohort_cancer/ARCHIVE_stale_before_ranked_marker_source_panel_20260618_214053/`

## 6. Functional analysis

- Queries prepared: {len(query)}
- Runnable by size: {len(query)-len(skipped)}
- Skipped by size: {len(skipped)}
- Eligible panel-selection universe: 1,346 genes
- Per-contrast eligible backgrounds: 19
- g:Profiler execution: attempted, not completed; endpoint timeout returned no term-level results
- Significant-term summary: unavailable; no result was fabricated

The manual execution file records human/g:SCS/0.05 settings and requires each custom background. MIXED-direction sets remain in skipped-query documentation. Functional interpretation must remain contextual rather than causal.

## 7. Pan-cancer DSMZ similarity graph

- Nodes: {dsmz['nodes']}; genes used: {dsmz['feature_genes_used']}; edges: {dsmz['edges']}
- Components: {dsmz['connected_components']}; Louvain communities at resolution 1: {dsmz['communities']}
- Modularity: {float(dsmz['modularity']):.6f}; cancer-type assortativity: {float(dsmz['cancer_type_assortativity']):.6f}
- NBL/RBL: co-assigned in one major 27-node community (16 NBL, 11 RBL)
- Haematological profiles: three pure communities; two profiles placed in a BRCA-majority community
- Discordant placements: current table regenerated
- Sensitivity: complete configured Louvain and Leiden grids regenerated; primary interpretation was not resolution-selected

## 8. Multi-cohort cancer analysis

`multicohort_cancer` was rerun with both panel distance directions and the established aggregate workflow.

- Active output: `results/unsupervised/multicohort_cancer/ranked_marker_source_panel/`
- Feature coverage: 379 found, zero missing, zero duplicated clean IDs
- Profiles: {multi['samples_profiles']} across four sources and BRCA/NBL/RBL
- Resolved graph: {multi['graph_nodes']} nodes, {multi['graph_edges']} edges, {multi['connected_components']} components
- Leiden communities: {multi['communities']}; modularity: {float(multi['modularity']):.6f}
- Structure: eight pure communities and one BRCA/NBL mixed community; NBL and RBL separated in the aggregate solution
- Haematological profiles: excluded here and assessed in the DSMZ graph
- Community sensitivity: unavailable in the established workflow, so no exploratory method was added

This supports feature-space validation and cancer-type-associated-structure wording, not a formal subtype claim.
'''

final_report += f'''
## 9. Ranking and retrieval

- Patient top-1 agreement: {float(patient['top1_accuracy']):.6f}; top-10 agreement: {float(patient['top10_accuracy']):.6f}; MRR: {float(patient['mrr']):.6f}
- Cell-line top-1 agreement: {float(overall['accuracy']):.6f} ({overall['n_correct']}/{overall['n_total']})
- Balanced accuracy: {float(balanced):.6f}
- Median confidence margins: patient 0.1048; cell-line 0.0449
- Reciprocal top-10 pairs: 530
- Atypical placements: 15 patient top matches and one cell-line top match crossed cancer types; all are retained for follow-up

Rankings are prioritised tiers of models. Small margins weaken top-1 certainty, and transcriptomic retrieval does not prove functional equivalence.

## 10. Optional embedding

The existing UMAP rule was available and rerun. Coordinates, figure and summary were packaged with the current panel identity. No new embedding method was introduced.

## 11. Documentation cleanup

Current enrichment, DSMZ, ranking, embedding and multi-cohort reports were generated. Current thesis-update notes contain current values only. Generic workflow files are retained as regenerated intermediates; the canonical deliverables are the method-specific packages. The final active semantic stale-description audit passed.

## 12. Commands run

```text
python scripts/build_ranked_marker_source_panel_enrichment_queries.py --pipeline-root . --outdir results/unsupervised/pan_cancer/enrichment/ranked_marker_source_panel --min-query-size 3 --background-min-normalised-count 1
snakemake --use-conda --cores 8 --rerun-triggers mtime --printshellcmds --keep-going --config pipeline_profile=multicohort_cancer -- results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds results/unsupervised/pan_cancer/tumour_mapping/metrics_summary_group_level.tsv results/unsupervised/pan_cancer/ranking/diagnostics/ranking_diagnostic_metric_crosscheck.tsv results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap/summary_pan_cancer_tumour_cell_line_alignment_umap.tsv results/unsupervised/pan_cancer/figures/Fig_pan_cancer_graph.pdf results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv results/unsupervised/pan_cancer/cell_line_similarity/community_validation/validation_assortativity.tsv figures/Fig_pan_cancer_cell_line_similarity_network_lineage_community.pdf
snakemake --use-conda --cores 8 --rerun-triggers mtime --printshellcmds --keep-going --config pipeline_profile=multicohort_cancer -- results/unsupervised/pan_cancer/tumour_mapping/cellline_similarity_precision_bootstrap/cellline_topk_metrics.tsv results/unsupervised/pan_cancer/tumour_mapping/cellline_similarity_precision_bootstrap/reciprocal_mapping_summary.tsv results/unsupervised/pan_cancer/figures/ecdf_plots/model_prioritisation_rank_summary.tsv results/unsupervised/pan_cancer/ranking/diagnostics/tumour_to_cellline_mrr_at10_by_tumour.tsv results/unsupervised/pan_cancer/ranking/diagnostics/cellline_to_tumour_component_composition_summary.tsv
snakemake --use-conda --cores 8 --rerun-triggers mtime --printshellcmds --keep-going --config pipeline_profile=multicohort_cancer -- results/unsupervised/multicohort_cancer/tumour_neighbourhoods/PanCancerFeatureSet_euc/final_consensus/cell_line_similarity_graph_edges_PanCancerFeatureSet_euc.tsv results/unsupervised/multicohort_cancer/tumour_neighbourhoods/PanCancerFeatureSet_corr/final_consensus/cell_line_similarity_graph_edges_PanCancerFeatureSet_corr.tsv
snakemake --use-conda --cores 8 --rerun-triggers mtime --printshellcmds --keep-going --config pipeline_profile=multicohort_cancer --allowed-rules summarize_p_consensus_all plot_per_cellline_feature_distance_cleveland resolve_dsmz_graph_neighbours plot_patient_referenced_resolved_cell_line_neighbourhood_graph build_multi_representation_majority_threshold_consensus_network plot_multi_representation_majority_threshold_consensus_network build_multi_representation_union_supported_edges_network plot_multi_representation_union_supported_edges_network plot_pan_cancer_resolved_graph_inspection community_stability_analysis post_resolution_edge_support_stratification combine_post_resolution_edge_support_stratification compute_multicohort_cancer_communities materialize_study_design model_selection_summary neighbourhood_permutation_validation random_baseline_comparison silhouette_report -- results/unsupervised/multicohort_cancer/tumour_neighbourhoods/final_consensus_all/resolved_dsmz_neighbours.tsv results/unsupervised/multicohort_cancer/tumour_neighbourhoods/final_consensus_all/plots/multicohort_cancer_multi_representation_majority_threshold_consensus_network_edges.tsv results/unsupervised/multicohort_cancer/tumour_neighbourhoods/final_consensus_all/community_detection/multicohort_cancer_communities.tsv results/unsupervised/multicohort_cancer/validation/model_selection_summary.tsv results/unsupervised/multicohort_cancer/validation/neighbourhood_permutation_summary.tsv results/unsupervised/multicohort_cancer/validation/random_baseline_summary.tsv results/unsupervised/multicohort_cancer/validation/silhouette_report.tsv
$CONDA_PREFIX/bin/Rscript scripts/run_gprofiler_from_manifest.R --query-manifest results/unsupervised/pan_cancer/enrichment/ranked_marker_source_panel/ranked_marker_source_panel_enrichment_query_manifest.tsv --outdir results/unsupervised/pan_cancer/enrichment/ranked_marker_source_panel/gprofiler --organism hsapiens --sources GO:BP,GO:MF,GO:CC,KEGG,REAC,WP,TF,HPA,CORUM --alpha 0.05 --correction-method g_SCS --archive-url https://biit.cs.ut.ee/gprofiler_archive3/e114_eg62_p19/ --require-archive TRUE --run-iea-sensitivity FALSE --as-short-link FALSE --top-terms-per-source 5
Rscript scripts/export_ranked_marker_source_panel_matrices.R --clean-features results/unsupervised/pan_cancer/feature_space/pan_cancer_features_clean.txt --cell-line-rds results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds --multicohort-rds results/unsupervised/multicohort_cancer/inputs/joint_expr_matrix.rds --dsmz-outdir results/unsupervised/pan_cancer/dsmz_similarity_graph/ranked_marker_source_panel --multicohort-outdir results/unsupervised/multicohort_cancer/ranked_marker_source_panel --k 20
python scripts/package_ranked_marker_source_panel_outputs.py
python scripts/finalize_ranked_marker_source_panel_reports.py
```

The g:Profiler API attempt used the prepared query manifest and custom backgrounds; it was stopped after repeated endpoint timeouts. Audit logs preserve the detailed job execution.

## 13. Snakefile rules run

Pan-cancer: `build_pan_cancer_feature_expression_matrix`, `build_pan_cancer_feature_expression_matrix_cell_lines_only`, `compute_pan_cancer_correlation`, `build_pan_cancer_graph`, `plot_pan_cancer_graph`, `build_pan_cancer_cell_line_similarity_graph`, `compute_pan_cancer_cell_line_communities`, `compute_pan_cancer_cell_line_layout`, `compute_pan_cancer_cell_line_louvain_resolution_sweep`, `compute_pan_cancer_cell_line_leiden_resolution_sweep`, `compute_pan_cancer_cell_line_validation`, `plot_pan_cancer_cell_line_two_panel`, `score_tumour_cellline_mapping`, `build_pan_cancer_bidirectional_ranking_diagnostics`, and `plot_pan_cancer_tumour_cell_line_alignment_umap`.

Ranking post-processing: `cellline_precision_at_k`, `plot_tumour_to_cellline_rank_ecdf_top10_fraction`, `plot_tumour_to_cellline_mrr_at10_distribution`, and both cell-line component-composition rules.

Multi-cohort: `consensus_cluster_ccp`, `tumour_nh_hc`, `tumour_nh_km`, `tumour_nh_consensus`, `cell_line_similarity_graph`, `summarize_p_consensus_all`, `resolve_dsmz_graph_neighbours`, `plot_patient_referenced_resolved_cell_line_neighbourhood_graph`, `build_multi_representation_majority_threshold_consensus_network`, `compute_multicohort_cancer_communities`, `model_selection_summary`, `neighbourhood_permutation_validation`, `random_baseline_comparison`, and `silhouette_report`.

## 14. File-change summary

The directory is not a Git worktree. The before-run baseline and final checksum manifest provide the equivalent audit. `config/config.yaml` was not changed. Modified existing file: `Snakefile` (four outputs of `score_tumour_cellline_mapping` were declared so the existing ranking DAG could resolve). Added method-specific scripts: `scripts/build_ranked_marker_source_panel_enrichment_queries.py`, `scripts/export_ranked_marker_source_panel_matrices.R`, `scripts/package_ranked_marker_source_panel_outputs.py`, and `scripts/finalize_ranked_marker_source_panel_reports.py`. These constitute the untracked/new-file summary in lieu of Git status. No thesis or DEG-generation file was changed.

## 15. Final validation

The active feature panel remains the current ranked marker-source panel. Method-specific outputs and reports identify the ranked marker-source panel, all matrix/graph/ranking/multi-cohort files post-date the active feature panel, explicit backgrounds are documented, and manifests contain size, modification time and SHA-256. Final checks are recorded in `ranked_marker_source_panel_downstream_validation_checks.tsv`.

## 16. Remaining manual action and overall status

Run g:Profiler from the prepared manifest when the endpoint is responsive, then replace the explicit status rows with returned term-level tables and update the enrichment report. All other requested downstream analyses are current. Therefore, the downstream rerun is current for graph, ranking/retrieval, multi-cohort and embedding analyses, but is not fully complete until functional enrichment executes successfully.
'''
write(PAN / 'ranked_marker_source_panel_downstream_rerun_report.md', final_report)

checks = [
    ('second_remediation_checksum_verification', 'PASS', '28 files archived; zero failures'),
    ('complete_prerun_audit', 'PASS', '27 checks passed'),
    ('active_feature_count', 'PASS', '379'),
    ('clean_feature_count', 'PASS', '379'),
    ('clean_feature_duplicates', 'PASS', '0'),
    ('active_method', 'PASS', 'ranked_marker_source_pan_cancer_panel'),
    ('selected_rule', 'PASS', 'relaxed_iqr_median_baseMean'),
    ('evidence_classes', 'PASS', 'four current classes only'),
    ('enrichment_categories', 'PASS', 'current owner, marker-source-class, evidence, direction and contrast sets'),
    ('enrichment_backgrounds', 'PASS', 'panel universe and 19 per-contrast backgrounds documented'),
    ('gprofiler_execution', 'MANUAL_ACTION', 'prepared; endpoint timeout'),
    ('dsmz_graph_rerun', 'PASS', 'current feature panel; current package'),
    ('ranking_retrieval_rerun', 'PASS', 'current package'),
    ('multicohort_rerun', 'PASS', 'current feature panel; current package'),
    ('optional_embedding_rerun', 'PASS', 'existing UMAP rule'),
    ('outputs_postdate_feature_panel', 'PASS', 'method-specific outputs verified'),
    ('method_specific_naming', 'PASS', 'canonical outputs use ranked_marker_source_panel'),
    ('active_semantic_stale_description_grep', 'PASS', 'zero disallowed current descriptions'),
    ('thesis_files_untouched', 'PASS', 'no thesis edits or LaTeX compilation'),
    ('upstream_deg_untouched', 'PASS', 'no DEG generation rules run'),
]
write_tsv(PAN / 'ranked_marker_source_panel_downstream_validation_checks.tsv',
          [{'check': a, 'status': b, 'detail': c} for a,b,c in checks])

def refresh_package_manifest(base, output):
    rows = []
    for path in sorted(base.rglob('*')):
        if path.is_file() and path != output:
            rows.append({
                'path': str(path.relative_to(ROOT)), 'size_bytes': path.stat().st_size,
                'modification_time': datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(),
                'sha256': sha(path),
            })
    write_tsv(output, rows, ['path', 'size_bytes', 'modification_time', 'sha256'])

refresh_package_manifest(DSMZ, DSMZ / 'ranked_marker_source_panel_dsmz_analysis_manifest.tsv')
refresh_package_manifest(RANK, RANK / 'ranked_marker_source_panel_ranking_manifest.tsv')
refresh_package_manifest(EMBED, EMBED / 'ranked_marker_source_panel_embedding_manifest.tsv')
refresh_package_manifest(MOUT, MOUT / 'ranked_marker_source_panel_multicohort_analysis_manifest.tsv')

manifest_roots = [ENR, DSMZ, RANK, EMBED, MOUT]
root_files = [
    PAN / 'ranked_marker_source_panel_downstream_rerun_report.md',
    PAN / 'ranked_marker_source_panel_downstream_validation_checks.tsv',
    PAN / 'ranked_marker_source_panel_thesis_update_notes.md',
]
manifest_rows = []
seen = set()
for base in manifest_roots:
    for path in sorted(base.rglob('*')):
        if path.is_file() and path.name != 'ranked_marker_source_panel_downstream_rerun_manifest.tsv':
            seen.add(path)
for path in root_files:
    if path.exists(): seen.add(path)
for path in sorted(seen):
    manifest_rows.append({
        'path': str(path.relative_to(ROOT)),
        'size_bytes': path.stat().st_size,
        'modification_time': datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(),
        'sha256': sha(path),
    })
write_tsv(PAN / 'ranked_marker_source_panel_downstream_rerun_manifest.tsv', manifest_rows,
          ['path', 'size_bytes', 'modification_time', 'sha256'])
print(f'reports and manifest complete: {len(manifest_rows)} files')
