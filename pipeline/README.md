Unsupervised Clustering and Tumour Neighbourhood Analysis Pipeline

This repository provides a Snakemake workflow for unsupervised analysis of cancer cell line and tumour **bulk RNA-seq** gene expression data. The pipeline selects informative gene sets, applies multiple clustering strategies, and constructs tumour neighbourhood profiles to classify cell lines and quantify cell line–tumour similarity.

## Contents

* Overview
* Requirements
* Installation
* Configuration
* Running the pipeline
* Pipeline outputs
* Output layout
* Troubleshooting
* Citation
* License
* Disclaimer

---

## Overview

The pipeline is designed to compare cancer cell lines (for example, DSMZ collections) with primary tumour samples (for example, TCGA cohorts). The primary objective is to identify cell lines whose transcriptional profiles best match tumour samples and to assess the robustness of these matches across multiple analysis settings.

### Main functions

* Generates gene sets using eight feature-selection methods
* Performs hierarchical clustering and k-means in expression space and PCA space
* Computes stable cluster assignments using consensus clustering
* Derives tumour neighbourhood profiles for each cell line
* Builds cell line similarity graphs based on tumour neighbourhood profiles
* Supports multiple analysis directions (HVG vs PAM50; Euclidean vs correlation distance)

### Supported directions

| Direction    | Gene set        | Distance             |
| ------------ | --------------- | -------------------- |
| `hvg_euc`    | HVG set (Top N) | Euclidean            |
| `hvg_corr`   | HVG set (Top N) | Correlation distance |
| `pam50_euc`  | PAM50 genes     | Euclidean            |
| `pam50_corr` | PAM50 genes     | Correlation distance |

Note: k-means clustering is applied only to Euclidean-distance directions.

---

## Requirements

### Software

* Snakemake (≥ 7)
* R (≥ 4.3)
* Conda or Mamba (recommended for reproducible environments)

### R packages

Required R packages are installed via the Conda environments provided in `envs/`. Commonly used packages include:

* `yaml`, `dplyr`, `readr`, `tibble`, `ggplot2`
* `matrixStats`, `stats`, `cluster`
* `ConsensusClusterPlus`
* `WGCNA`, `igraph`
* `UpSetR`, `pheatmap`, `uwot`, `ggrepel`

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-org/unsupervised-clustering-pipeline.git
cd unsupervised-clustering-pipeline
```

### 2. Create Conda environments

```bash
conda env create -f envs/tcga-r-env.yaml
conda activate tcga-r-env

conda env create -f envs/tumour_nh_qc.yaml
```

### 3. Verify installation

```bash
snakemake --version
Rscript -e "library(WGCNA); library(ConsensusClusterPlus); cat('OK\n')"
```

---

## Configuration

All pipeline settings are defined in `config/config.yaml`.

At minimum, the pipeline requires:

* VST-normalised expression matrices (`genes × samples`)

  * joint matrix (cell lines + tumours)
  * cell line-only matrix
  * tumour-only matrix
* Cell line metadata table
* Optional PAM50 inputs for PAM50-based analyses

Example (shortened):

```yaml
paths:
  unsup_root: "/path/to/results/unsupervised"
  vst_joint_rds: "/path/to/data/vst_joint_matrix.rds"
  cell_vst_rds: "/path/to/data/cell_lines_vst.rds"
  tumour_vst_rds: "/path/to/data/tumours_vst.rds"
  dsmz_meta_csv: "/path/to/metadata/dsmz_cell_lines.csv"

features:
  pam50_ensembl_gene_list: "/path/to/resources/pam50_ensembl_ids.txt"

feature_selection:
  top_n_method: 3000
  final_top: 500

tumour_neighbourhoods:
  directions: ["hvg_euc","hvg_corr","pam50_euc","pam50_corr"]
  p_consensus_threshold: 0.7
```

---

## Running the pipeline

### Local execution

```bash
conda activate tcga-r-env

# Dry run
snakemake -n all

# Run with 8 cores
snakemake -j 8 all

# Run with per-rule Conda environments
snakemake -j 8 --use-conda all
```

### Run selected stages

```bash
snakemake -j 4 feature_selection_unsupervised
snakemake -j 8 agnostic_all_directions
snakemake -j 8 consensus_all_directions
```

### SLURM execution

```bash
snakemake -j 100 --profile slurm all
```

---

## Pipeline outputs

For each configured direction, the pipeline generates:

* Feature-selection outputs (rank tables, gene lists, overlap plots)
* Clustering results for all methods and parameter settings
* Consensus clustering outputs
* Tumour neighbourhood profiles (cell line × tumour)
* Cell line similarity graphs and summary tables
* QC tables and visualisations (including UMAP where configured)

---

## Output layout

A typical output structure is:

* `results/unsupervised/feature_selection_unsupervised/`
* `results/unsupervised/agnostic_clustering/<direction>/`
* `results/unsupervised/consensus/<direction>/`
* `results/unsupervised/tumour_neighbourhoods/<direction>/`
* `results/unsupervised/logs/`

Exact file names and subdirectories may vary depending on configuration and enabled rules.

---

## Troubleshooting

Common issues and checks:

* **Missing inputs**: verify file paths and permissions in `config/config.yaml`
* **Empty outputs**: inspect log files in `results/unsupervised/logs/`
* **Memory errors**: increase available RAM or reduce block sizes in feature selection
* **Environment issues**: recreate Conda environments and ensure `conda init` is configured

Useful commands:

```bash
snakemake -n -p <target>
snakemake --unlock
snakemake -j 4 --forcerun <rule_name> all
```

---

## Citation

If this pipeline is used in academic work, please cite the relevant software:

* Snakemake: Mölder et al. (2021), *F1000Research*
* WGCNA: Langfelder & Horvath (2008), *BMC Bioinformatics*
* ConsensusClusterPlus: Wilkerson & Hayes (2010), *Bioinformatics*

---

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.

---

## Disclaimer
Parts of the codebase were initially drafted, reviewed, or refactored with assistance from AI-based tools. 