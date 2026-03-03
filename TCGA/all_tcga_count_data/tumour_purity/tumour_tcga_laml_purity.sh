#!/bin/bash
#SBATCH --job-name=tumour_laml_purity
#SBATCH --chdir=/home/chu25/TCGA
#SBATCH --output=/home/chu25/TCGA/logs/tumour_laml_purity_analysis_%x_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/tumour_laml_purity_analysis_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=168:00:00

set -euo pipefail

# Ensure log dir exists
mkdir -p /home/chu25/TCGA/logs

# Activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tidyestimate || { echo "Failed to activate Conda environment tidyestimate"; exit 1; }

# Run R script
R_SCRIPT="/home/chu25/TCGA/all_tcga_count_data/tumour_purity/TCGA_LAML.R"
Rscript "$R_SCRIPT"
