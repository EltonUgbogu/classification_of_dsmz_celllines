# Version 1 inputs used by Version 2

Version 2 reads the following protected Version 1 outputs. `Snakefile.v2` does not regenerate or overwrite these files.

| Version 1 input | Role in Version 2 | Required Version 1 target if missing |
|---|---|---|
| `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_graph_edges.tsv` | Cell-line graph edge support used for weighted consensus evidence. | `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_graph_edges.tsv` |
| `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_node_metadata.tsv` | Cell-line lineage/source metadata used for node annotation. | `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_node_metadata.tsv` |
| `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_communities.tsv` | Louvain community assignments used for boundary and bridge-like node classification. | `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_communities.tsv` |
| `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_leiden_communities.tsv` | Leiden assignments used as a comparison layer for mixed membership checks. | `results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_leiden_communities.tsv` |
| `results/unsupervised/pan_cancer/tumour_mapping/tumour_to_cellline_similarity/tumour_to_cellline_group_rankings.tsv` | Patient tumour to cell-line ranking evidence used for agreement diagnostics. | `results/unsupervised/pan_cancer/tumour_mapping/tumour_to_cellline_similarity/tumour_to_cellline_group_rankings.tsv` |
| `results/unsupervised/pan_cancer/tumour_mapping/cellline_to_tumour_similarity/cellline_to_tumour_rankings.tsv` | Cell-line to tumour ranking evidence used for reciprocal agreement diagnostics. | `results/unsupervised/pan_cancer/tumour_mapping/cellline_to_tumour_similarity/cellline_to_tumour_rankings.tsv` |
| `results/unsupervised/pan_cancer/tumour_mapping/cellline_similarity_precision_bootstrap/cellline_centred_rank_summary.tsv` | Cell-line centred decision margins used for ranking boundary classification. | `results/unsupervised/pan_cancer/tumour_mapping/cellline_similarity_precision_bootstrap/cellline_centred_rank_summary.tsv` |

If one of these files is missing, `Snakefile.v2` fails at workflow load time with a message naming the missing Version 1 file and the Version 1 target to run first.

