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
* Supports multiple analysis directions built from your configured feature-set methods/distances; PAM50 is optional and only used when enabled for BRCA

### Supported directions

Directions are generated automatically from the feature-set grid in `config.yaml`:

* `feature_sets.methods` × `feature_sets.distances` (defaults: 8 methods × 2 distances → 16 directions)
* Optional PAM50 directions (`pam50_euc`, `pam50_corr`) are appended only when `analysis.use_pam50: true` (BRCA runs).

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

### Profiles

The configuration file now supports multiple **profiles** so you can keep BRCA (default), RBL, and NBL settings in one place. Set the `profile` key (or override on the command line) to select which block is used.

* `brca` (default): BRCA/DSMZ run with PAM50 (optional) plus the full feature-set grid.
* `rbl`: RBL analysis using the full feature-set grid (PAM50 disabled).
* `nbl`: NBL analysis using the full feature-set grid (PAM50 disabled).

To run with a different profile, pass `--config profile=<name>`:

```bash
snakemake -j 8 --config profile=rbl all
```

### Minimum inputs

At minimum, the pipeline requires:

* VST-normalised expression matrices (`genes × samples`)
  * joint matrix (cell lines + tumours)
  * cell line-only matrix
  * tumour-only matrix
* Cell line metadata table
* Optional PAM50 inputs for PAM50-based analyses (BRCA profile only)

### Testing with or without data

* **Dry-run / rule graph**: You can check configuration and rule wiring without any data by running `snakemake -n all` or `snakemake --rulegraph`. This validates the config structure and environment activation but does not touch the 5.6 GB inputs.
* **Smoke test**: Running any rule for real requires the data files defined in your selected profile (for example, the 5.6 GB VST matrices). If you want a lighter run, point the config to a downsampled subset with the same filenames/columns; the pipeline will treat it the same way.
* **Full run**: Uses the full-size VST matrices and metadata. Ensure the configured paths are accessible and that you have enough disk and memory for intermediate outputs.

The repository does **not** include sample input data. To run anything beyond a dry-run, supply your own matrices/metadata (or lightweight substitutes) and update `config/config.yaml` accordingly.

Example (shortened) for a single profile:

```yaml
profile: brca

profiles:
  brca:
    analysis:
      cancer_type: BRCA
      use_pam50: true

    paths:
      unsup_root: "/path/to/results/unsupervised"
      vst_joint_rds: "/path/to/data/vst_joint_matrix.rds"
      cell_vst_rds: "/path/to/data/cell_lines_vst.rds"
      tumour_vst_rds: "/path/to/data/tumours_vst.rds"
      dsmz_meta_csv: "/path/to/metadata/dsmz_cell_lines.csv"

    pam50:
      ensembl_gene_list: "/path/to/resources/pam50_ensembl_ids.txt"

    feature_sets:
      methods: ["Variance", "MAD", "MeanAbsDev", "Entropy", "PCA", "Spearman", "MX", "kTotal"]
      distances: ["euc", "corr"]

    tumour_neighbourhoods:
      # Directions will be constructed from feature_sets + use_pam50
      p_consensus_threshold: 0.7
```

---

## Optional utilities

These helpers are **not** wired into the Snakefile but can be useful for one-off tasks or exploratory analyses:

* **`rebuild_joint_vst.R`** — builds the joint VST matrix from the cohort-specific VST RDS inputs if it does not exist.
* **`make_cancer_types_table.R`** — builds the DSMZ cell type table used in the HTML report.
* **`compute_similarity_graph.R`** — builds a sample–sample similarity matrix and kNN edge list from any expression matrix (genes × samples). It is **optional** and not called by the Snakefile. Point `--expr_rds` at any pipeline-produced expression matrix (e.g., a feature-set submatrix under `results/feature_selection_unsupervised` or a tumour-neighbourhood input RDS), keep the default `--similarity spearman` unless you want Pearson on log-scale data, and set `--k` to something smaller than your sample count (e.g., `min(30, n/3)`) to avoid kNN failures.

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

### Running with profiles

```bash
# Use default profile (brca)
snakemake -j 8 all

# Override profile on command line
snakemake -j 8 --config profile=rbl all
snakemake -j 8 --config profile=nbl all
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



---

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.

---

## Disclaimer
This codebase was drafted, reviewed, or refactored with assistance from AI-based tools.


