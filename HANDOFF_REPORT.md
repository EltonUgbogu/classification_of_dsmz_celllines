Use this as the **final handoff instruction** to the coding agent. I have made it strict and unambiguous.

The key correction is that the **cell line similarity graph must represent cell line transcriptomic similarity derived from clinical patient sample neighbourhood consensus profiles**, not per direction community stability. The Methods ground truth defines each cell line as a vector of consensus probabilities over clinical patient samples, then computes cell line to cell line Pearson similarity from those vectors.  The attached plotting scripts should be treated as implementation references for graph layout, edge support, node statistics, and betweenness based highlighting.  

```text
FINAL CODING AGENT INSTRUCTION

You are updating the thesis network figures. The previous community stability heatmap direction was rejected. The final figures must follow the ground truth methods and the attached reference scripts.

Core conceptual correction
==========================

There are two valid graph levels only:

1. Cell line similarity network construction
2. Graph based consensus resolution of DSMZ cell line neighbourhoods

Do not merge these graph levels.

Do not use the previous “community stability across feature selection and distance representations” framing as the final thesis figure. That framing is no longer valid for the final thesis plots.

The graphs should represent transcriptomic similarity of DSMZ cell lines in relation to clinical patient samples.

Important:
- The graph nodes are cell lines.
- The graph does not need to contain clinical patient samples as nodes.
- Clinical patient samples enter through the cell line representation.
- Each cell line is represented by its clinical patient sample neighbourhood consensus profile.
- Edges between cell lines represent similarity between those clinical patient sample referenced profiles.

For the cell line similarity network construction level:
- Each cell line is represented by a vector of p_consensus values over clinical patient samples.
- Pairwise cell line similarity is Pearson correlation between these consensus profile vectors.
- The similarity matrix is thresholded at the empirical 90th percentile of pairwise similarities.
- Edge weights are the Pearson similarity values.
- Leiden community detection is the primary community detection method.
- Louvain is applied in parallel only as a comparison.
- Leiden and Louvain are not to be treated as per direction or per feature selection community assignments.
- Do not compare raw community IDs across feature selection directions.
- Do not generate or use pairwise community co assignment stability heatmaps as a main thesis output.

For the graph based consensus resolution level:
- Use the resolved neighbour relationships from the graph based consensus resolution.
- The resolved graph structure must match the previously accepted per disease resolved component graphs.
- Do not replace the old resolved graph topology with a different topology produced by the rejected combined script.
- Nodes are cell lines.
- Edges are final stable neighbour relationships.
- Connected components and isolates are reported from this final resolved graph.
- Betweenness centrality is computed within each connected component using unweighted shortest paths.
- Within each non isolate connected component, highlight exactly the node with the highest betweenness centrality.
- This node is the component central node.
- Do not replace this with top degree, top weighted degree, top 1 to 3 nodes per cohort, or top N percent across the whole graph.

Plots to keep
=============

Keep and regenerate the following final thesis figures.

Main Figure 1
-------------

File name:

Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.pdf
Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.svg
Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.png

Title:

Clinical patient sample referenced DSMZ cell line similarity networks

Panels:

A. BRCA cell line similarity network
B. NBL cell line similarity network
C. RBL cell line similarity network
D. Similarity network construction summary

Use the attached per disease consensus graph outputs as the visual and structural reference:
- Fig_dsmz_brca_consensus_graph
- Fig_nbl_consensus_dsmz_graph
- Fig_rbl_dsmz_consensus_graph

Panel A to C visual encoding:

- node = DSMZ cell line
- edge = Pearson similarity between clinical patient sample neighbourhood consensus profiles
- node colour = Leiden community
- node size = graph degree
- node shape = outlier status, where outlier means degree zero in the cell line similarity network
- edge width = Pearson similarity or edge support if available
- edge style = multi supported versus single supported only if this information is available in the input edge table
- optional node border = high betweenness node in the similarity network, but only if clearly defined and documented

Required legend entries:

- Leiden community
- Outlier, degree zero
- Edge weight
- Multi supported edge, if shown
- Single supported edge, if shown

Panel D must summarise per cohort:

- number of cell lines
- number of retained edges
- Pearson similarity threshold value
- threshold rule, fixed text: empirical 90th percentile
- number of Leiden communities
- number of Louvain communities
- Leiden modularity, if available
- Louvain modularity, if available
- number of outliers, degree zero
- median degree

Machine readable source tables for Main Figure 1:

cell_line_similarity_network_node_annotations.tsv
cell_line_similarity_network_edges.tsv
cell_line_similarity_network_summary.tsv
cell_line_similarity_network_layout_coordinates.tsv

Required columns for cell_line_similarity_network_node_annotations.tsv:

cohort
cell_line_raw_id
cell_line_display
degree
weighted_degree
betweenness
is_outlier_degree_zero
leiden_community
louvain_community
display_label
x
y

Required columns for cell_line_similarity_network_edges.tsv:

cohort
node1
node2
node1_display
node2_display
pearson_similarity
edge_weight
similarity_threshold
threshold_rule
support_directions, if available
edge_support_label, if available

Required columns for cell_line_similarity_network_summary.tsv:

cohort
n_cell_lines
n_edges
similarity_metric
similarity_threshold_rule
similarity_threshold_value
n_leiden_communities
n_louvain_communities
leiden_modularity
louvain_modularity
n_outliers_degree_zero
median_degree
mean_degree
median_edge_weight

Main Figure 2
-------------

File name:

Fig_consensus_resolved_DSMZ_neighbourhoods_combined.pdf
Fig_consensus_resolved_DSMZ_neighbourhoods_combined.svg
Fig_consensus_resolved_DSMZ_neighbourhoods_combined.png

Title:

Graph based consensus resolved DSMZ cell line neighbourhoods

Panels:

A. BRCA resolved DSMZ cell line neighbourhoods
B. NBL resolved DSMZ cell line neighbourhoods
C. RBL resolved DSMZ cell line neighbourhoods
D. Resolved component summary

Use the attached previously accepted resolved component plots as the structural ground truth:
- breast_cancer_resolved_components
- nbl_resolved_components
- rbl_resolved_components

The regenerated combined figure must match these graph structures. If the current combined resolved graph has a different topology from these accepted per disease resolved graphs, discard the current combined graph and regenerate it from the correct resolved neighbour inputs.

Panel A to C visual encoding:

- node = DSMZ cell line
- edge = final stable neighbour relationship
- node colour = final connected component
- node shape = isolate status
- node size = degree
- thick black outline or clear halo = component central node
- component central node = highest betweenness centrality node within that non isolate connected component
- edge width = support directions, if available
- edge style = multi supported versus single supported, if available

Do not use node colour for centrality.
Do not use Leiden or Louvain community colours in the resolved graph.
Do not use the same visual encoding for anchor status and centrality.
If anchor status is shown, use a separate marker that does not conflict with the component central node marker.
If this becomes visually crowded, omit anchor status from the main figure and keep it in the table.

Centrality rule:

For each cohort:
1. Build the resolved graph from final neighbour relationships.
2. Decompose the graph into connected components.
3. Ignore isolate components for central node highlighting.
4. For each non isolate component, compute betweenness centrality using unweighted shortest paths.
5. Select the node with the highest betweenness centrality.
6. If there is a tie, break ties by higher degree, then alphabetical cell line display label.
7. Mark only that node as is_component_central_node = TRUE.

Panel D must summarise per cohort:

- number of cell lines
- number of edges
- number of connected components
- number of isolates
- largest component size
- median degree
- number of component central nodes
- median edge support, if available

Machine readable source tables for Main Figure 2:

resolved_DSMZ_neighbourhood_node_summary.tsv
resolved_DSMZ_component_annotations.tsv
resolved_DSMZ_graph_edges.tsv
resolved_DSMZ_centrality_rankings.tsv
resolved_DSMZ_graph_layout_coordinates.tsv

Required columns for resolved_DSMZ_neighbourhood_node_summary.tsv:

cohort
cell_line_raw_id
cell_line_display
connected_component_id
connected_component_size
degree
weighted_degree
betweenness
closeness, if valid
is_isolate
is_component_central_node
centrality_metric_used
centrality_rank_within_component
final_neighbours
n_final_neighbours

Required columns for resolved_DSMZ_component_annotations.tsv:

cohort
connected_component_id
component_label
n_cell_lines
member_cell_lines
n_edges
density
median_degree
max_degree
median_betweenness
max_betweenness
component_central_node
component_central_node_betweenness
component_central_node_degree
n_isolates
annotation_status

Required columns for resolved_DSMZ_graph_edges.tsv:

cohort
node1
node2
node1_display
node2_display
support_directions, if available
support_weight_mean, if available
edge_support_label
edge_source

Required columns for resolved_DSMZ_centrality_rankings.tsv:

cohort
cell_line_display
connected_component_id
connected_component_size
degree
weighted_degree
betweenness
closeness, if valid
is_isolate
is_component_central_node
centrality_metric_used
centrality_rank_within_component
centrality_label

Allowed values:

centrality_metric_used = betweenness

centrality_label:
- component central node
- not highlighted
- isolate

Supplementary Figure 1
----------------------

File name:

Fig_supp_pan_cancer_cell_line_transcriptomic_similarity_graph.pdf
Fig_supp_pan_cancer_cell_line_transcriptomic_similarity_graph.svg
Fig_supp_pan_cancer_cell_line_transcriptomic_similarity_graph.png

This figure is optional.

It is valid only if labelled correctly.

Correct title:

Pan cancer DSMZ cell line transcriptomic similarity graph

Do not call this a pan cancer cell line to clinical patient sample graph.
Do not imply that this graph contains clinical patient samples as nodes.
Do not imply that this graph is the same as the per disease clinical patient sample referenced similarity graph.

Valid interpretation:

- nodes = cell lines
- edges = transcriptomic similarity between cell lines
- communities may show pan cancer cell line transcriptomic structure
- this is a broader cell line similarity graph, not a clinical patient sample referenced graph

Supplementary Figure 2
----------------------

File name:

Fig_supp_resolved_DSMZ_component_graphs.pdf
Fig_supp_resolved_DSMZ_component_graphs.svg
Fig_supp_resolved_DSMZ_component_graphs.png

This figure should show larger per disease resolved graphs if the combined main figure is too compact.

Panels:

A. BRCA detailed resolved component graph
B. NBL detailed resolved component graph
C. RBL detailed resolved component graph

Use the older accepted resolved component plots as the visual reference.
These supplementary graphs should make node labels and component central nodes easier to read.

Plots to discard
================

Discard the following from the final thesis output:

1. Fig_cell_line_community_stability.pdf
2. Fig_cell_line_community_stability.svg
3. Fig_cell_line_community_stability.png
4. Fig_cell_line_community_stability_leiden.pdf
5. Fig_cell_line_community_stability_leiden.svg
6. Fig_cell_line_community_stability_leiden.png

Reason:
These figures frame the result as community stability across feature selection and distance representations. That is not the final ground truth interpretation.

Also discard or stop producing these tables if they are only used for the rejected community stability framing:

- cell_line_method_community_assignments.tsv
- cell_line_louvain_coassignment_matrix.tsv
- cell_line_leiden_coassignment_matrix.tsv
- component_community_overlap_summary.tsv
- supp_cell_line_community_stability.tsv

Do not include pairwise community co assignment heatmaps in the main thesis.

Do not include a 3 by 2 Louvain and Leiden heatmap figure.

Do not include per direction Leiden or Louvain community assignments as a main result.

Do not present raw Leiden or Louvain community IDs as comparable across feature selection methods or directions.

Relabelling rules
=================

Apply these relabellings consistently.

Replace:

“DSMZ Cell line Consensus Graph”

with:

“Clinical patient sample referenced DSMZ cell line similarity network”

Replace:

“Consensus graph”

with:

“Cell line similarity network constructed from clinical patient sample neighbourhood profiles”

Replace:

“Resolved cell line similarity network”

with:

“Graph based consensus resolved DSMZ cell line neighbourhoods”

Replace:

“High betweenness node”

with:

“Component central node”

Replace:

“Co assignment frequency”

with nothing in the final main figures, because the co assignment heatmap figure is discarded.

Replace:

“Pan cancer cell line to clinical patient sample graph”

with:

“Pan cancer DSMZ cell line transcriptomic similarity graph”

Use British English:

- tumour
- neighbourhood
- colour
- visualisation

Avoid unnecessary hyphens:

Use:
- cell line
- pan cancer
- graph based
- clinical patient sample referenced
- multi supported
- single supported

Validation rules
================

The script must stop with an informative error if:

1. BRCA, NBL, or RBL disease specific inputs are missing for either main figure.
2. The resolved graph topology differs from the accepted resolved component graph inputs.
3. The similarity network figure is built from per direction community co assignment matrices instead of clinical patient sample consensus profile based similarity.
4. Leiden or Louvain communities are loaded as per direction or per feature selection outputs for the final main similarity network figure.
5. The resolved graph is coloured by Leiden or Louvain communities.
6. Component central nodes are not computed from within component betweenness centrality.
7. Isolates are marked as component central nodes.
8. RBL replicate labels are ambiguous.
9. Mixed short and long identifiers are used in the same final figure.
10. Any plotted value cannot be traced to a machine readable TSV.

Required provenance
===================

Create:

figure_provenance.tsv

Required columns:

figure_file
figure_role
graph_level
script
input_files
date_time
git_commit, if available
layout_seed
n_nodes
n_edges
visual_encodings
output_width
output_height
source_tables
notes

Allowed graph_level values:

- clinical_patient_sample_referenced_similarity_network
- resolved_neighbour_graph
- pan_cancer_cell_line_similarity_graph

Acceptance criteria
===================

The task is complete only when:

1. The final main figure set contains exactly two main figures:
   - clinical patient sample referenced DSMZ cell line similarity networks
   - graph based consensus resolved DSMZ cell line neighbourhoods

2. The disease specific cell line similarity network figure has:
   - A. BRCA
   - B. NBL
   - C. RBL
   - D. similarity network construction summary

3. The disease specific resolved graph figure has:
   - A. BRCA
   - B. NBL
   - C. RBL
   - D. resolved component summary

4. The rejected community stability heatmap figures are not included as final thesis figures.

5. Leiden is used as the primary community detection display for the cell line similarity networks.

6. Louvain is retained only as a comparison table or optional supplementary annotation, not as a per direction stability plot.

7. The resolved graph topology matches the previously accepted resolved component graphs.

8. Each non isolate connected component in the resolved graph has exactly one component central node.

9. Component central nodes are selected by highest within component betweenness centrality using unweighted shortest paths.

10. The pan cancer graph, if produced, is labelled only as a pan cancer cell line transcriptomic similarity graph.

11. Every plotted value is traceable to a TSV.

12. All outputs are written as PDF, SVG, and 300 dpi PNG.

13. All figure labels use British English and avoid unnecessary hyphens.

Final report to return
======================

After implementation, report:

1. Exact output directory.
2. Exact paths of all final figures.
3. Exact paths of all source TSVs.
4. Number of BRCA, NBL, and RBL nodes in the clinical patient sample referenced similarity networks.
5. Number of BRCA, NBL, and RBL edges in the clinical patient sample referenced similarity networks.
6. Number of Leiden and Louvain communities per cohort.
7. Number of BRCA, NBL, and RBL nodes in the resolved neighbour graphs.
8. Number of BRCA, NBL, and RBL edges in the resolved neighbour graphs.
9. Number of connected components and isolates per cohort.
10. Component central nodes per cohort and component.
11. Confirmation that the rejected community stability heatmap figures are no longer generated as final thesis outputs.
12. Confirmation that the pan cancer graph, if generated, is not labelled as a cell line to clinical patient sample graph.
```


relavant paths:
/Users/eltonugbogu/Downloads/thesis_latex/stab_test
/Users/eltonugbogu/Downloads/thesis_latex/methods.tex
/Users/eltonugbogu/Downloads/thesis_latex/figures
'/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all'
'/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/nbl/tumour_neighbourhoods/final_consensus_all'
'/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all'
'/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/pan_cancer/inputs'