# HANDOFF REPORT — Cell line similarity and consensus resolved network main figures

**Status:** completed under the corrected final thesis instruction.
**Generated:** 2026-05-08
**Output directory:** `/Users/eltonugbogu/Downloads/thesis_latex/stab_test/`
**Script:** `scripts/plot_publication_cell_line_similarity_and_resolved_networks.py`
**Git commit:** `141f7cc2`

---

## 1. Scope correction

The previous implementation was rejected. It produced a community-stability
output set framed around per-direction Louvain/Leiden co-assignment. That
direction has been **discarded**. The previous script
`scripts/plot_publication_network_components_and_community_stability.py` has
been **deleted**, and the following obsolete outputs are no longer present in
`stab_test/`:

* `Fig_cell_line_community_stability.{pdf,svg,png}`
* `Fig_cell_line_community_stability_leiden.{pdf,svg,png}`
* `cell_line_louvain_coassignment_matrix.tsv`
* `cell_line_leiden_coassignment_matrix.tsv`
* `component_community_overlap_summary.tsv`
* `supp_cell_line_community_stability.tsv`
* `cell_line_method_community_assignments.tsv`
* `community_stability_results_summary.txt`
* `community_input_validation.tsv`

The `stab_test/` directory has been emptied and rebuilt with the corrected
output set only.

---

## 2. Final thesis output set

Two main figures plus their supporting tables and provenance.

### Main Figure 1 — clinical patient sample referenced cell line similarity networks

```
Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.pdf
Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.svg
Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.png
```

Panels:

* **A.** BRCA clinical patient sample referenced cell line similarity network
* **B.** NBL  clinical patient sample referenced cell line similarity network
* **C.** RBL  clinical patient sample referenced cell line similarity network
* **D.** Similarity network construction summary

Encodings:

* nodes  = DSMZ cell lines
* edges  = Pearson similarity between clinical patient sample
  neighbourhood consensus profiles
* node colour = Leiden community (per-cohort palette; extended categorical
  scale so each Leiden community is uniquely coloured)
* Louvain is reported only as a comparison count in Panel D
* edge width and alpha are scaled by Pearson similarity
* isolates (degree 0) are drawn as diamond markers
* no co-assignment heatmaps and no per-direction community-stability content

Per-cohort winning direction (used for edges and node colour):

| Cohort | Winning direction | Source                                                                             |
|:-------|:------------------|:-----------------------------------------------------------------------------------|
| BRCA   | `PCA_corr`        | `winning_direction.txt`                                                            |
| NBL    | `kTotal_euc`      | `winning_direction.txt`                                                            |
| RBL    | `hvg_euc`         | `p_consensus_direction_summary.tsv` top by `frac_ge_thr` covering all resolved cells |

### Main Figure 2 — graph based consensus resolved DSMZ neighbourhoods

```
Fig_consensus_resolved_DSMZ_neighbourhoods_combined.pdf
Fig_consensus_resolved_DSMZ_neighbourhoods_combined.svg
Fig_consensus_resolved_DSMZ_neighbourhoods_combined.png
```

Panels:

* **A.** BRCA graph based consensus resolved DSMZ cell line neighbourhoods
* **B.** NBL  graph based consensus resolved DSMZ cell line neighbourhoods
* **C.** RBL  graph based consensus resolved DSMZ cell line neighbourhoods
* **D.** Resolved component summary

Encodings:

* nodes = DSMZ cell lines
* edges = final stable neighbour relationships from `resolved_dsmz_neighbours.tsv`
  (BRCA, NBL) and `resolved_dsmz_neighbors.tsv` (RBL)
* node colour = final connected component (Wong 8-colour CB-safe palette;
  isolates and singletons in grey)
* component central node = highest betweenness centrality node within each
  non-isolate connected component, highlighted with a thick red border
  (4.5 pt, `#B22222`)
* edge style: solid for `support_directions ≥ 2`, dashed for
  `support_directions == 1`; width scales with `support_directions`
* isolates drawn as diamonds
* topology matches the accepted per-disease resolved component graphs

---

## 3. Per-cohort numbers

### Main Figure 1 — similarity networks

| Cohort | Direction      | Cell lines | Edges | Isolates | Leiden communities | Louvain communities | Median similarity | Min similarity |
|:------:|:---------------|-----------:|------:|---------:|-------------------:|--------------------:|------------------:|---------------:|
| BRCA   | `PCA_corr`     |       29   |   42  |     3    |               29   |                7    |             0.680 |          0.595 |
| NBL    | `kTotal_euc`   |       18   |   17  |     4    |               18   |                8    |             0.557 |          0.468 |
| RBL    | `hvg_euc`      |       11   |    6  |     5    |               11   |                7    |             0.830 |          0.828 |

> Note. In all three cohorts the upstream Leiden detector returns
> singleton-dominated communities at the pipeline resolution (n_Leiden
> equals n_cell_lines in each cohort). Each cell line is therefore its
> own Leiden community and is given a unique colour in panels A–C. The
> Louvain count column in Panel D shows the alternative (multi-member)
> grouping for comparison only; raw community IDs are not compared
> across cohorts.

### Main Figure 2 — resolved networks

| Cohort | Nodes | Edges | Components | Multi-node | Isolates | Largest component | Median degree | Central nodes (component → cell line)              |
|:------:|------:|------:|-----------:|-----------:|---------:|------------------:|--------------:|:----------------------------------------------------|
| BRCA   |  29   |  39   |     7      |     4      |    3     |        12         |     2.0       | 0:HCC_1143; 1:BT_474; 2:EFM_192B; 3:MDA_MB_468      |
| NBL    |  18   |  17   |     6      |     2      |    4     |         7         |     2.0       | 0:CHP_134; 1:LAN1                                   |
| RBL    |  11   |   6   |     6      |     2      |    4     |         5         |     1.0       | 0:RBL_18_r1; 1:RBL_15_r4                            |

Central nodes are the highest-betweenness-centrality node within each
non-isolate component. They are exported in
`resolved_component_central_nodes.tsv` for cross-referencing.

---

## 4. Files in `stab_test/`

### Figures

| File | Format | Purpose |
|:--|:--|:--|
| `Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.pdf/.svg/.png` | vector + raster | Main Figure 1 |
| `Fig_consensus_resolved_DSMZ_neighbourhoods_combined.pdf/.svg/.png`                   | vector + raster | Main Figure 2 |

### Supporting tables

| File | Rows | Description |
|:--|--:|:--|
| `cell_line_similarity_network_edges.tsv`            |  65 | Per-cohort edge list with Pearson similarity, Leiden/Louvain ID per endpoint |
| `cell_line_similarity_network_nodes.tsv`            |  58 | Per-cohort node table with Leiden/Louvain ID, degree, mean/max edge similarity, isolate flag |
| `similarity_network_construction_summary.tsv`       |   3 | Cohort × winning-direction summary used in Panel D of Main Figure 1 |
| `resolved_DSMZ_neighbourhood_node_summary.tsv`      |  58 | Per-cohort resolved node table with component, central-node flag, betweenness, anchor flag, median incident edge support |
| `resolved_DSMZ_component_annotations.tsv`           |  19 | Per-component annotations (members, central node, anchor list, density, edge support) |
| `supp_resolved_DSMZ_component_annotations.tsv`      |  19 | Slimmer supplementary version of the component annotations |
| `resolved_component_central_nodes.tsv`              |   8 | Central node per non-isolate component with betweenness value |
| `cell_line_network_layout_coordinates.tsv`          | 116 | Reproducible 2-D layout coordinates for both figures |
| `network_input_validation.tsv`                      |   6 | One row per cohort × scope; all rows `PASS` |
| `network_results_summary.txt`                       |   – | Plain-text Results-section summary |
| `figure_provenance.tsv`                             |   1 | Inputs, outputs, encodings, palette, layout seed, git commit |
| `plot_publication_cell_line_similarity_and_resolved_networks.log` | – | Run log |

---

## 5. Reproducing the run

```bash
cd /Users/eltonugbogu/classification_of_dsmz_celllines
RBL_BASE="/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/rbl/tumour_neighbourhoods"

python scripts/plot_publication_cell_line_similarity_and_resolved_networks.py \
  --brca-tn-dir results/unsupervised/brca/tumour_neighbourhoods \
  --brca-resolved results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/resolved_dsmz_neighbours.tsv \
  --nbl-tn-dir  results/unsupervised/nbl/tumour_neighbourhoods \
  --nbl-resolved results/unsupervised/nbl/tumour_neighbourhoods/final_consensus_all/resolved_dsmz_neighbours.tsv \
  --rbl-tn-dir  "$RBL_BASE" \
  --rbl-resolved "$RBL_BASE/final_consensus_all/resolved_dsmz_neighbors.tsv" \
  --rbl-anchors  "$RBL_BASE/final_consensus_all/anchor_components.tsv" \
  --outdir /Users/eltonugbogu/Downloads/thesis_latex/stab_test \
  --layout-seed 42
```

The script auto-detects each cohort's winning direction from
`final_consensus_all/winning_direction.txt`. RBL has no such file; the
script falls back to the highest-`frac_ge_thr` direction in
`p_consensus_direction_summary.tsv` that covers all 11 resolved cells
(`hvg_euc`). Override with `--<cohort>-direction NAME` if required.

---

## 6. What is intentionally absent

* No co-assignment heatmaps (Louvain or Leiden).
* No per-direction Louvain or Leiden stability tables.
* No within- vs between-component co-assignment summary.
* No `Fig_cell_line_community_stability*` outputs.
* No script with "community stability" in its name.

---

## 7. Acceptance check against the corrected instruction

| Requirement                                                                                                            | Status   | Evidence |
|:-----------------------------------------------------------------------------------------------------------------------|:--------:|:--|
| Discard `Fig_cell_line_community_stability.{pdf,svg,png}`                                                              | done     | files removed; no longer produced |
| Discard `Fig_cell_line_community_stability_leiden.{pdf,svg,png}`                                                       | done     | files removed; no longer produced |
| Discard `cell_line_louvain_coassignment_matrix.tsv`                                                                    | done     | file removed; no longer produced |
| Discard `cell_line_leiden_coassignment_matrix.tsv`                                                                     | done     | file removed; no longer produced |
| Discard `component_community_overlap_summary.tsv`                                                                      | done     | file removed; no longer produced |
| Discard `supp_cell_line_community_stability.tsv`                                                                       | done     | file removed; no longer produced |
| Discard `cell_line_method_community_assignments.tsv`                                                                   | done     | file removed; no longer produced |
| Generate `Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.{pdf,svg,png}` with panels A/B/C/D    | done     | three formats present in `stab_test/` |
| Main Figure 1 nodes = DSMZ cell lines                                                                                  | done     | one cell line per node |
| Main Figure 1 edges = Pearson similarity between clinical patient sample neighbourhood consensus profiles              | done     | values from `cell_line_similarity_graph_edges_<dir>.tsv` (BRCA, NBL) and `DSMZ_DSMZ_graph_edges_<dir>.tsv` (RBL); winning direction per cohort |
| Main Figure 1 node colour = Leiden community                                                                           | done     | extended categorical palette; one colour per Leiden ID |
| Main Figure 1: Louvain comparison only                                                                                 | done     | reported only as a count in Panel D and in the supporting nodes TSV |
| Main Figure 1: no co-assignment heatmaps; no per-direction community-stability interpretation                          | done     | no such panels and no such outputs are produced |
| Generate `Fig_consensus_resolved_DSMZ_neighbourhoods_combined.{pdf,svg,png}` with panels A/B/C/D                       | done     | three formats present in `stab_test/` |
| Main Figure 2 nodes = DSMZ cell lines                                                                                  | done     | one cell line per node |
| Main Figure 2 edges = final stable neighbour relationships                                                             | done     | parsed from `final_neighbours` / `final_neighbors` columns in the resolved files |
| Main Figure 2 node colour = final connected component                                                                  | done     | per-cohort components from `resolved_dsmz_neighbours.tsv` |
| Main Figure 2 component central node = highest betweenness centrality within each non-isolate connected component     | done     | thick red border in panels; values in `resolved_component_central_nodes.tsv` and `resolved_DSMZ_component_annotations.tsv` |
| Main Figure 2 topology matches the accepted per-disease resolved component graphs                                      | done     | edges and components computed directly from the accepted resolved files |
| Script no longer named after "community stability"                                                                     | done     | `plot_publication_cell_line_similarity_and_resolved_networks.py`; previous script deleted |

All items pass.
