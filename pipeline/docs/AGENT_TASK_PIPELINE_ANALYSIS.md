# Agent Task: Unsupervised Clustering and Tumour Neighbourhood Pipeline — Step-by-Step Analysis

This document is a **complete task file for an agent**, detailing the analysis progress and which scripts to use. It is derived from `/work/ugbogu/pipeline/Snakefile`.

---

## Prerequisites

- **Config:** `config/config.yaml` must define a profile (e.g. `heme`, `brca`, `nbl`, `rbl`) with:
  - `paths.vst_joint_rds`, `paths.cell_vst_rds`, `paths.tumour_vst_rds`
  - `paths.unsup_root`, `paths.dsmz_meta_csv` (or `meta_tsv` for joint metadata)
  - `feature_sets.methods` and `feature_sets.distances` (or `tumour_neighbourhoods.directions` for an explicit list)
  - `features.hvg_final_gene_list` and optionally `features.hvg_joint_scores_tsv` if using HVG
- **Conda env:** Rules use `envs/tcga-r-env.yaml` (R + optparse, yaml, dplyr, etc.).
- **Run from pipeline root:** Use `--snakefile Snakefile` (or full path to Snakefile) and `--configfile config/config.yaml --config pipeline_profile=<name>`.

---

## Pipeline Overview (per-profile, e.g. HEME)

| Stage | Snakemake target / rule | Purpose |
|-------|-------------------------|--------|
| **PRESTEP** | Feature selection + optional featureset matrices | Gene lists and per-feature expression submatrices |
| **STEP 1** | Agnostic clustering | Per-direction filtered matrices and HC/kmeans clusters |
| **STEP 2** | Consensus clustering | CCP consensus RDS per direction × kind |
| **STEP 3** | Tumour neighbourhoods | p-consensus and final neighbourhood RDS/TSV per direction |
| **STEP 4** | DSMZ–DSMZ similarity | Graphs and plots per direction |
| **STEP 5** | Cross-direction summary | Winner direction, summaries, figures |
| **STEP 6** | QC / UMAP | QC TSV and UMAP PDF per direction |
| **FINAL** | `pipeline_targets` | Touches `results/unsupervised/<profile>/pipeline_targets.done` |

---

## PRESTEP: Feature selection and featureset matrices

### 1. Feature selection (gene lists only)

- **Script:** `scripts/feature_selection_unsupervised.R`
- **Snakemake rule:** `feature_selection_unsupervised`
- **Inputs:** `paths.vst_joint_rds`, `config/config.yaml`
- **Outputs:**
  - `results/unsupervised/<profile>/feature_selection_unsupervised/feature_sets_top500/genes_top500_{Variance,MAD,Entropy,PCA,Spearman,MX,kTotal}.txt`
  - Optional side outputs: `hvg_joint_scores.tsv`, `hvg_final_gene_list.txt`, UpSet plots, etc.
- **What it does:** Ranks genes by multiple methods, writes top-N gene lists per method. Required for all downstream steps that use “direction” (e.g. Variance_euc, hvg_corr).

**Run:**
```bash
snakemake feature_selection_unsupervised \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 1 --use-conda
```

### 2. Featureset expression matrices (per-method submatrices)

- **Script:** `scripts/build_tumour_neighbourhood_input.R`
- **Snakemake rule:** `build_tumour_neighbourhood_input_featureset` (one job per feature, e.g. Entropy, MAD, Variance, MX, kTotal). Optional convenience rule: `build_featureset_matrices` (if present in Snakefile).
- **Inputs:** `paths.vst_joint_rds`, `feature_sets_top500/genes_top500_{feature}.txt`, `paths.dsmz_meta_csv` (or joint metadata TSV with `sample_id` or `sample_name` and a cell-line-like column, e.g. `lineage` or `Cell_line`).
- **Outputs:**
  - `results/unsupervised/<profile>/feature_selection_unsupervised/featuresets/<Feature>/expr_submatrix.rds`
  - `results/unsupervised/<profile>/tumour_neighbourhoods_input/cell_line_to_original_sample_id_{feature}.rds`
- **What it does:** Subsets the joint VST to the gene list for that feature, writes samples × genes RDS and a DSMZ sample → cell-line mapping. Required for consensus and tumour neighbourhood steps that use non-HVG/non-PAM50 directions.

**Run (all featureset matrices):**
```bash
snakemake build_featureset_matrices \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 8 --use-conda --rerun-incomplete
```
If `build_featureset_matrices` does not exist, request the individual `expr_submatrix.rds` outputs for each feature instead.

---

## STEP 1: Agnostic clustering

- **Scripts:**  
  - `scripts/build_agnostic_direction_mats.R`  
  - `scripts/pca_hc_cell.R`, `pca_hc_tumour.R`, `pca_hc_cell_tumour.R`  
  - `scripts/pca_kmeans_cell.R`, `pca_kmeans_tumour.R`, `pca_kmeans_cell_tumour.R`  
  - `scripts/hc_cell.R`, `hc_tumour.R`, `hc_cell_tumour.R`  
  - `scripts/kmeans_cell.R`, `kmeans_tumour.R`, `kmeans_cell_tumour.R`
- **Snakemake rule (aggregate):** `agnostic_all_directions`
- **Inputs:**  
  - Cell/tumour VST RDS; per-direction gene list from feature selection (or HVG list for hvg_*).  
  - Directions come from `tumour_neighbourhoods.directions` or from `feature_sets.methods` × `feature_sets.distances` (e.g. Variance_euc, hvg_corr).
- **Outputs:**  
  - `results/unsupervised/<profile>/agnostic_clustering/<direction>/inputs/cell_expr.rds`, `tumour_expr.rds`  
  - For each direction: `pca_hc_cell`, `pca_hc_tumour`, `pca_hc_cell_tumour`, `hc_cell`, `hc_tumour`, `hc_cell_tumour`, and for _euc only: `pca_kmeans_*`, `kmeans_*` → `*_clusters_optimal.rds`

**Run:**
```bash
snakemake agnostic_all_directions \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 8 --use-conda
```

---

## STEP 2: Consensus clustering (CCP)

- **Script:** `scripts/consensus_ccp_cell_tumour.R`
- **Snakemake rule (aggregate):** `consensus_all_directions`
- **Inputs:** Config, profile, per-direction agnostic cluster RDS; HVG list path (e.g. `feature_sets_top500/genes_top500_MX.txt`).
- **Outputs:**  
  - `results/unsupervised/<profile>/consensus/<direction>/<kind>/<kind>_clusters_optimal.rds`  
  - Kinds: e.g. `ccp_hc_expr_cell_tumour`, `ccp_kmeans_expr_cell_tumour`, and other `ccp_*` kinds defined in the Snakefile.

**Run:**
```bash
snakemake consensus_all_directions \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 8 --use-conda
```

---

## STEP 3: Tumour neighbourhood consensus

### 3a. Build tumour neighbourhood inputs (HVG and PAM50, if used)

- **Script:** `scripts/build_tumour_neighbourhood_input.R`
- **Rules:** `build_tumour_neighbourhood_input_hvg`, and optionally `build_tumour_neighbourhood_input_pam50` if PAM50 is enabled.
- **Outputs:**  
  - HVG: `tumour_neighbourhoods_input/expr_hvg.rds`, `cell_line_to_original_sample_id_hvg.rds`  
  - PAM50 (if used): analogous PAM50 RDS files.

### 3b. Compute tumour neighbourhoods (HC then kmeans for _euc)

- **Script:** `scripts/comp_tumour_neighbourhoods.R`
- **Rules:** `tumour_neighbourhoods_hc_hvg`, `tumour_neighbourhoods_hc_featureset`, `tumour_neighbourhoods_km_any`, etc.
- **Inputs:** Per-direction expr RDS and mapping RDS (from featuresets or HVG/PAM50).
- **Outputs:**  
  - `tumour_neighbourhoods/<direction>/.tumour_neighbourhoods_done`  
  - Various `Top_m_long_*.csv` and related files under each direction.

### 3c. p-consensus (final neighbourhood consensus per direction)

- **Script:** `scripts/tumour_neighbourhood_p_consensus.R`
- **Rules:** `tumour_neighbourhood_p_consensus_hvg`, `tumour_neighbourhood_p_consensus_featureset` (and PAM50 if enabled).
- **Outputs:**  
  - `tumour_neighbourhoods/<direction>/final_consensus/Final_consensus_tumour_neighbourhoods_<direction>.rds`  
  - `Final_consensus_tumour_neighbourhoods_<direction>.tsv`

**Run (all of STEP 3 is demanded by `pipeline_targets`):** run `pipeline_targets` or the individual rules above.

---

## STEP 4: DSMZ–DSMZ similarity graph

- **Script:** `scripts/compute_dsmz_dsmz_similarity.R`
- **Snakemake rule:** `dsmz_dsmz_similarity_graph` (one per direction).
- **Inputs:** `Final_consensus_tumour_neighbourhoods_<direction>.rds` (and config/profile).
- **Outputs:**  
  - Under `tumour_neighbourhoods/<direction>/final_consensus/`:  
    - `DSMZ_DSMZ_similarity_matrix_<direction>.rds`, `DSMZ_DSMZ_similarity_pairs_<direction>.tsv`  
    - `DSMZ_DSMZ_graph_edges_<direction>.tsv`, node summaries, community tables  
    - `Fig_DSMZ_DSMZ_*.pdf` (histogram, scatter, heatmap, Leiden/Louvain graphs)

**Run:**
```bash
snakemake dsmz_dsmz_similarity_graph \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 8 --use-conda
```
(Or rely on `pipeline_targets` to request these outputs.)

---

## STEP 5: Cross-direction summary

- **Script:** `scripts/summarize_p_consensus_all.R`
- **Snakemake rule:** `summarize_p_consensus_all`
- **Inputs:** All `Final_consensus_tumour_neighbourhoods_<direction>.rds`, config (and `tumour_neighbourhoods.p_consensus_threshold`).
- **Outputs:**  
  - Under `tumour_neighbourhoods/final_consensus_all/`:  
    - `p_consensus_direction_summary.tsv`, `p_consensus_cellline_direction_summary.tsv`  
    - `p_consensus_winner_map_per_cell_line.tsv`, `p_consensus_winners_by_max_p.tsv`, `p_consensus_winners_by_frac_ge_thr.tsv`  
    - `p_consensus_best_cell_lines_ranked.tsv`, composite weights, top fraction/score TSVs  
    - `Fig_p_consensus_direction_comparison.pdf`, `Fig_p_consensus_composite_PCA_*.pdf`  
  - And e.g. `winning_direction.txt` (used by later characterisation).

**Run:**
```bash
snakemake summarize_p_consensus_all \
  --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=heme -j 1 --use-conda
```

---

## STEP 6: QC and UMAP (per direction)

- **Script:** `scripts/tumour_neighbourhood_qc_umap.R` (or equivalent; invoked by rules such as `tumour_neighbourhood_qc_umap_hvg`, `tumour_neighbourhood_qc_umap_featureset`).
- **Outputs:**  
  - `tumour_neighbourhoods/<direction>/qc/nh_qc_summary.tsv`  
  - `tumour_neighbourhoods/<direction>/qc/nh_umap.tsv`  
  - `tumour_neighbourhoods/<direction>/qc/nh_umap.pdf`

These are requested by `pipeline_targets` as the final QC/UMAP deliverables.

---

## Characterisation (Layer 3, optional)

- **Scripts:**  
  - `scripts/characterize_clusters.R` (cluster characterisation and cell-line characterisation for the selected direction).  
  - `scripts/compute_representativeness_score.R` (representativeness score from cell-line characterisation).
- **Snakemake rules:** `lock_selected_direction`, `representativeness_score` (and any rule that runs `characterize_clusters.R`).
- **Inputs:** Winner direction from `final_consensus_all/winning_direction.txt`, and the corresponding tumour neighbourhood and cluster outputs.
- **Outputs:**  
  - `characterization/tumour/SELECTED/`, `cell_line_representativeness_ranked.tsv`, `cell_line_representativeness_top10_per_cluster.tsv`, etc.

---

## Full pipeline run (single profile)

To run the **entire** per-profile pipeline up to the main target:

```bash
cd /work/ugbogu/pipeline
SNAKEMAKE_PROFILE= snakemake pipeline_targets \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --config pipeline_profile=heme \
  -j 8 --use-conda --rerun-incomplete
```

This will:

1. Run feature selection (gene lists).
2. Build featureset matrices only if some downstream rule requests them (e.g. `build_tumour_neighbourhood_input_featureset`). If your Snakefile has `build_featureset_matrices`, run that first if you want all featureset expr/map RDS files before the rest.
3. Run agnostic clustering for all directions.
4. Run consensus clustering for all directions × kinds.
5. Build tumour neighbourhood inputs (HVG and per-feature), run tumour neighbourhoods and p-consensus.
6. Run DSMZ–DSMZ similarity graphs.
7. Run summarize_p_consensus_all.
8. Run QC/UMAP steps.
9. Touch `results/unsupervised/<profile>/pipeline_targets.done`.

---

## Script → step quick reference

| Script | Stage |
|--------|--------|
| `feature_selection_unsupervised.R` | PRESTEP: gene lists |
| `build_tumour_neighbourhood_input.R` | PRESTEP: featureset matrices; STEP 3a: HVG (and optionally PAM50) inputs |
| `build_agnostic_direction_mats.R` | STEP 1: filtered cell/tumour matrices per direction |
| `pca_hc_*.R`, `hc_*.R`, `pca_kmeans_*.R`, `kmeans_*.R` | STEP 1: agnostic clustering |
| `consensus_ccp_cell_tumour.R` | STEP 2: CCP consensus |
| `comp_tumour_neighbourhoods.R` | STEP 3b: tumour neighbourhoods |
| `tumour_neighbourhood_p_consensus.R` | STEP 3c: p-consensus |
| `compute_dsmz_dsmz_similarity.R` | STEP 4: DSMZ–DSMZ graphs |
| `summarize_p_consensus_all.R` | STEP 5: cross-direction summary |
| `tumour_neighbourhood_qc_umap*.R` | STEP 6: QC/UMAP |
| `characterize_clusters.R`, `compute_representativeness_score.R` | Layer 3: characterisation |

---

## Notes for agents

- Always pass `--config pipeline_profile=<profile>` (e.g. `heme`) when running Snakemake; the Snakefile expects this.
- If you see “No rule to produce” for a path under `feature_selection_unsupervised/featuresets/`, ensure either `build_featureset_matrices` or the individual `build_tumour_neighbourhood_input_featureset` outputs are requested; the rule uses a `{feature}` wildcard.
- Metadata for `build_tumour_neighbourhood_input.R` must contain a sample-ID column (`sample_name` or `sample_id`) and a cell-line–like column (`Cell_line`, `lineage`, or `source`); the script supports TSV and CSV.
- Logs are written under `logs/<profile>/` (e.g. `logs/heme/`). Check the corresponding `.log` file for any rule that fails.
