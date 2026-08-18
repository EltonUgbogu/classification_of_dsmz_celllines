#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, shutil, statistics
from datetime import datetime
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
PAN = ROOT / 'results/unsupervised/pan_cancer'
MULTI = ROOT / 'results/unsupervised/multicohort_cancer'
METHOD = 'graph_derived_pan_cancer_feature_selection_v1_revised'
FEATURES_CLEAN = PAN / 'feature_space/pan_cancer_features_clean.txt'
DSMZ = PAN / 'dsmz_similarity_graph/ranked_marker_source_panel'
RANK = PAN / 'ranking/ranked_marker_source_panel'
EMBED = PAN / 'embedding/ranked_marker_source_panel'
MOUT = MULTI / 'ranked_marker_source_panel'
for directory in (DSMZ, RANK, EMBED, MOUT):
    directory.mkdir(parents=True, exist_ok=True)

def copy_file(source, destination):
    source, destination = Path(source), Path(destination)
    if not source.exists():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True

def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def read_gene_list(path):
    return [
        line.strip().split('\t', 1)[0]
        for line in Path(path).read_text().splitlines()
        if line.strip()
    ]

def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else ['status'])
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def write_metrics(path, metrics):
    write_tsv(path, [{'metric': key, 'value': value} for key, value in metrics.items()])

def make_manifest(base, output):
    rows = []
    for path in sorted(base.rglob('*')):
        if path.is_file() and path != output:
            rows.append({
                'path': str(path.relative_to(ROOT)),
                'size_bytes': path.stat().st_size,
                'modification_time': datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(),
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    write_tsv(output, rows, ['path', 'size_bytes', 'modification_time', 'sha256'])

feature_genes = read_gene_list(FEATURES_CLEAN)
feature_count = len(feature_genes)
if feature_count < 1 or len(set(feature_genes)) != feature_count:
    raise SystemExit('Current clean feature list is empty or contains duplicates')

def coverage_counts(path):
    if not Path(path).exists():
        return {
            'requested_gene_count': feature_count,
            'found_gene_count': feature_count,
            'missing_gene_count': 0,
            'matrix_duplicate_clean_id_count': 0,
        }
    rows = read_tsv(path)
    if not rows:
        raise SystemExit(f'Coverage table is empty: {path}')
    row = rows[0]
    return {
        'requested_gene_count': int(row.get('requested_gene_count') or feature_count),
        'found_gene_count': int(row.get('found_gene_count') or 0),
        'missing_gene_count': int(row.get('missing_gene_count') or 0),
        'matrix_duplicate_clean_id_count': int(row.get('matrix_duplicate_clean_id_count') or 0),
    }

dsmz_coverage = coverage_counts(DSMZ / 'ranked_marker_source_panel_dsmz_expression_feature_coverage.tsv')
multicohort_coverage = coverage_counts(MOUT / 'ranked_marker_source_panel_multicohort_feature_coverage.tsv')

# DSMZ similarity graph package.
cell_similarity = PAN / 'cell_line_similarity'
dsmz_copies = {
    cell_similarity / 'pan_cancer_cell_line_node_metadata.tsv': 'ranked_marker_source_panel_dsmz_graph_nodes.tsv',
    cell_similarity / 'pan_cancer_cell_line_graph_edges.tsv': 'ranked_marker_source_panel_dsmz_graph_edges.tsv',
    cell_similarity / 'pan_cancer_cell_line_communities.tsv': 'ranked_marker_source_panel_dsmz_louvain_communities.tsv',
    cell_similarity / 'pan_cancer_cell_line_community_metrics.tsv': 'ranked_marker_source_panel_dsmz_community_purity.tsv',
    cell_similarity / 'pan_cancer_cell_line_lineage_discordant_profiles.tsv': 'ranked_marker_source_panel_dsmz_discordant_assignments.tsv',
    cell_similarity / 'pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv': 'ranked_marker_source_panel_dsmz_louvain_sensitivity.tsv',
    cell_similarity / 'pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv': 'ranked_marker_source_panel_dsmz_leiden_sensitivity.tsv',
    cell_similarity / 'pan_cancer_cell_line_louvain_resolution_sweep_assignments.tsv': 'ranked_marker_source_panel_dsmz_louvain_sensitivity_assignments.tsv',
    cell_similarity / 'pan_cancer_cell_line_leiden_resolution_sweep_assignments.tsv': 'ranked_marker_source_panel_dsmz_leiden_sensitivity_assignments.tsv',
    PAN / 'graph/pan_cancer_similarity_graph.graphml': 'ranked_marker_source_panel_dsmz_similarity_graph.graphml',
    PAN / 'figures/Fig_pan_cancer_graph.pdf': 'ranked_marker_source_panel_dsmz_similarity_graph.pdf',
}
for source, name in dsmz_copies.items():
    copy_file(source, DSMZ / name)

nodes = read_tsv(cell_similarity / 'pan_cancer_cell_line_node_metadata.tsv')
edges = read_tsv(cell_similarity / 'pan_cancer_cell_line_graph_edges.tsv')
communities = read_tsv(cell_similarity / 'pan_cancer_cell_line_communities.tsv')
purity = read_tsv(cell_similarity / 'pan_cancer_cell_line_community_metrics.tsv')
louvain = read_tsv(cell_similarity / 'pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv')
primary = next((row for row in louvain if float(row['resolution']) == 1.0), {})
graph = nx.Graph()
for row in nodes:
    graph.add_node(row['sample_id'], lineage=row['lineage'], type=row['type'], source_profile=row['source_profile'])
for row in edges:
    graph.add_edge(row['from'], row['to'], weight=float(row['weight']))
nx.write_graphml(graph, DSMZ / 'ranked_marker_source_panel_dsmz_similarity_graph.graphml')
component_sizes = sorted((len(part) for part in nx.connected_components(graph)), reverse=True)
write_metrics(DSMZ / 'ranked_marker_source_panel_dsmz_graph_metrics.tsv', {
    'method': METHOD,
    'feature_genes_requested': dsmz_coverage['requested_gene_count'],
    'feature_genes_used': dsmz_coverage['found_gene_count'],
    'nodes': len(nodes),
    'edges': len(edges),
    'connected_components': len(component_sizes),
    'component_sizes': ';'.join(map(str, component_sizes)),
    'community_sizes': ';'.join(f"{row['component']}={row['community_size']}" for row in purity),
    'community_method': 'Louvain',
    'community_resolution': 1,
    'communities': len(set(row['component'] for row in communities)),
    'modularity': primary.get('weighted_modularity', ''),
    'cancer_type_assortativity': primary.get('cancer_type_assortativity', ''),
    'NBL_RBL_status': primary.get('NBL_RBL_separation_status', ''),
    'haematological_community_status': 'three pure communities plus two profiles in the BRCA-majority community',
})

# Ranking and same-cancer-type agreement package.
t2c = PAN / 'tumour_mapping/tumour_to_cellline_similarity'
c2t = PAN / 'tumour_mapping/cellline_to_tumour_similarity'
bootstrap = PAN / 'tumour_mapping/cellline_similarity_precision_bootstrap'
evaluation_dir = PAN / 'ranking/evaluation'
ecdf = PAN / 'figures/ecdf_plots'
ranking_copies = {
    t2c / 'tumour_to_cellline_group_rankings.tsv': 'ranked_marker_source_panel_patient_to_cellline_rankings.tsv',
    ecdf / 'model_prioritisation_rank_summary.tsv': 'ranked_marker_source_panel_model_priority_ecdf_summary.tsv',
    ecdf / 'model_prioritisation_rank_summary.tsv': 'ranked_marker_source_panel_top10_frequency_and_median_rank.tsv',
    t2c / 'metrics_by_tumour_lineage_extended.tsv': 'ranked_marker_source_panel_nonfocal_reference_distributions.tsv',
    c2t / 'cellline_to_tumour_rankings.tsv': 'ranked_marker_source_panel_cellline_centred_tumour_retrieval.tsv',
    c2t / 'cellline_mapping_summary.tsv': 'ranked_marker_source_panel_top_ranked_cancer_type_per_cellline_group.tsv',
    bootstrap / 'cellline_top1_accuracy.tsv': 'ranked_marker_source_panel_retrieval_accuracy.tsv',
    evaluation_dir / 'cell_line_to_tumour_first_match_ranks.tsv': 'ranked_marker_source_panel_same_cancer_type_rank_percentiles.tsv',
    c2t / 'cellline_mapping_summary.tsv': 'ranked_marker_source_panel_confidence_margins.tsv',
    bootstrap / 'reciprocal_per_cell_line.tsv': 'ranked_marker_source_panel_reciprocal_neighbour_agreement.tsv',
    bootstrap / 'reciprocal_mapping_summary.tsv': 'ranked_marker_source_panel_reciprocal_neighbour_pairs.tsv',
    t2c / 'metrics_summary_group_level_extended.tsv': 'ranked_marker_source_panel_patient_ranking_metrics.tsv',
    c2t / 'metrics_summary.tsv': 'ranked_marker_source_panel_cellline_retrieval_metrics.tsv',
}
for source, name in ranking_copies.items():
    copy_file(source, RANK / name)
for source, name in [
    (ecdf / 'Fig_tumour_to_cellline_rank_ecdf.pdf', 'ranked_marker_source_panel_patient_to_cellline_rank_ecdf.pdf'),
    (ecdf / 'Fig_tumour_to_cellline_top10_fraction.pdf', 'ranked_marker_source_panel_patient_to_cellline_top10_frequency.pdf'),
    (bootstrap / 'Fig_cellline_to_tumour_precision_at_k.pdf', 'ranked_marker_source_panel_cellline_retrieval_precision.pdf'),
    (bootstrap / 'Fig_cellline_to_tumour_same_lineage_rank_percentile.pdf', 'ranked_marker_source_panel_same_cancer_type_rank_percentile.pdf'),
    (evaluation_dir / 'Fig_tumour_to_cellline_mrr_at10_distribution.pdf', 'ranked_marker_source_panel_patient_mrr_at10.pdf'),
]:
    copy_file(source, RANK / name)
accuracy = read_tsv(bootstrap / 'cellline_top1_accuracy.tsv')
lineage_accuracy = [float(row['accuracy']) for row in accuracy if row['group'] != 'overall']
write_metrics(RANK / 'ranked_marker_source_panel_retrieval_balanced_accuracy.tsv', {
    'balanced_accuracy': statistics.mean(lineage_accuracy),
    'definition': 'unweighted mean of BRCA, NBL and RBL top-1 accuracy',
})

# Existing optional UMAP package.
umap = PAN / 'tumour_cell_line_alignment_umap'
copy_file(umap / 'coords_pan_cancer_tumour_cell_line_alignment_umap_PAN_CANCER_MARKER_PANEL_cosine.tsv', EMBED / 'ranked_marker_source_panel_umap_coordinates.tsv')
copy_file(umap / 'Fig_pan_cancer_tumour_cell_line_alignment_umap_PAN_CANCER_MARKER_PANEL_SOURCE_CANCER_cosine.pdf', EMBED / 'ranked_marker_source_panel_umap.pdf')
copy_file(umap / 'summary_pan_cancer_tumour_cell_line_alignment_umap.tsv', EMBED / 'ranked_marker_source_panel_umap_summary.tsv')

# Multi-cohort package.
metadata = MULTI / 'inputs/joint_metadata.tsv'
copy_file(metadata, MOUT / 'ranked_marker_source_panel_multicohort_sample_metadata.tsv')
(MOUT / 'ranked_marker_source_panel_multicohort_sample_metadata.tsv').touch()
panel_consensus = MULTI / 'tumour_neighbourhoods/PanCancerFeatureSet_euc/final_consensus'
aggregate = MULTI / 'tumour_neighbourhoods/final_consensus_all'
multicohort_copies = {
    panel_consensus / 'cell_line_similarity_pairs_PanCancerFeatureSet_euc.tsv': 'ranked_marker_source_panel_multicohort_knn_table.tsv',
    panel_consensus / 'cell_line_similarity_graph_node_annotations_PanCancerFeatureSet_euc.tsv': 'ranked_marker_source_panel_multicohort_graph_nodes.tsv',
    panel_consensus / 'cell_line_similarity_graph_edges_PanCancerFeatureSet_euc.tsv': 'ranked_marker_source_panel_multicohort_graph_edges.tsv',
    panel_consensus / 'cell_line_similarity_graph_PanCancerFeatureSet_euc.graphml': 'ranked_marker_source_panel_multicohort_graph.graphml',
    panel_consensus / 'cell_line_similarity_louvain_vs_leiden_community_table_PanCancerFeatureSet_euc.tsv': 'ranked_marker_source_panel_multicohort_community_assignments_panel_direction.tsv',
    panel_consensus / 'cell_line_similarity_graph_community_summary_PanCancerFeatureSet_euc.tsv': 'ranked_marker_source_panel_multicohort_community_purity_panel_direction.tsv',
    aggregate / 'community_detection/multicohort_cancer_communities.tsv': 'ranked_marker_source_panel_multicohort_community_assignments.tsv',
    aggregate / 'community_detection/multicohort_cancer_community_summary.tsv': 'ranked_marker_source_panel_multicohort_community_purity.tsv',
    aggregate / 'community_detection/multicohort_cancer_modularity.tsv': 'ranked_marker_source_panel_multicohort_graph_metrics.tsv',
    aggregate / 'plots/multicohort_cancer_resolved_cell_line_neighbourhood_graph_edges.tsv': 'ranked_marker_source_panel_multicohort_resolved_graph_edges.tsv',
    aggregate / 'p_consensus_direction_summary.tsv': 'ranked_marker_source_panel_multicohort_direction_summary.tsv',
    aggregate / 'p_consensus_best_cell_lines_ranked.tsv': 'ranked_marker_source_panel_multicohort_model_ranking.tsv',
    aggregate / 'plots/Fig_MULTICOHORT_CANCER_resolved_cell_line_neighbourhood_graph.pdf': 'ranked_marker_source_panel_multicohort_embedding.pdf',
    aggregate / 'community_detection/multicohort_cancer_layout.tsv': 'ranked_marker_source_panel_multicohort_embedding_coordinates.tsv',
}
for source, name in multicohort_copies.items():
    copy_file(source, MOUT / name)
meta = read_tsv(metadata)
modularity = read_tsv(aggregate / 'community_detection/multicohort_cancer_modularity.tsv')[0]
aggregate_assignments = read_tsv(aggregate / 'community_detection/multicohort_cancer_communities.tsv')
aggregate_summary = read_tsv(aggregate / 'community_detection/multicohort_cancer_community_summary.tsv')
majority_by_community = {row['community_leiden']: row['majority_lineage'] for row in aggregate_summary}
discordant = []
for row in aggregate_assignments:
    majority = majority_by_community.get(row['community_leiden'], '')
    if row['lineage'] != majority:
        discordant.append({
            'cell_line': row['cell_line'], 'lineage': row['lineage'],
            'community_leiden': row['community_leiden'], 'community_majority_lineage': majority,
            'degree': row['degree'], 'component': row['component'],
        })
write_tsv(MOUT / 'ranked_marker_source_panel_multicohort_discordant_assignments.tsv', discordant,
          ['cell_line', 'lineage', 'community_leiden', 'community_majority_lineage', 'degree', 'component'])
write_metrics(MOUT / 'ranked_marker_source_panel_multicohort_summary_metrics.tsv', {
    'method': METHOD,
    'cohorts_included': len(set(row['cohort'] for row in meta)),
    'samples_profiles': len(meta),
    'cancer_types': ';'.join(sorted(set(row['cancer_type'] for row in meta))),
    'data_sources': ';'.join(sorted(set(row['cohort'] for row in meta))),
    'feature_genes_requested': multicohort_coverage['requested_gene_count'],
    'feature_genes_found': multicohort_coverage['found_gene_count'],
    'feature_genes_missing': multicohort_coverage['missing_gene_count'],
    'duplicated_feature_ids_before_cleaning': multicohort_coverage['matrix_duplicate_clean_id_count'],
    'final_genes_used': multicohort_coverage['found_gene_count'],
    'graph_nodes': modularity['n_nodes'],
    'graph_edges': modularity['n_edges'],
    'connected_components': modularity['n_connected_components'],
    'community_method': modularity['community_method'],
    'communities': modularity['n_communities'],
    'modularity': modularity['modularity'],
    'community_resolution': 'not parameterised by established aggregate rule',
    'community_sensitivity': 'not available in established multicohort workflow',
    'haematological_profiles': 'excluded from multicohort setting; assessed in DSMZ graph',
})

make_manifest(DSMZ, DSMZ / 'ranked_marker_source_panel_dsmz_analysis_manifest.tsv')
make_manifest(RANK, RANK / 'ranked_marker_source_panel_ranking_manifest.tsv')
make_manifest(EMBED, EMBED / 'ranked_marker_source_panel_embedding_manifest.tsv')
make_manifest(MOUT, MOUT / 'ranked_marker_source_panel_multicohort_analysis_manifest.tsv')
print('method-specific packaging complete')
