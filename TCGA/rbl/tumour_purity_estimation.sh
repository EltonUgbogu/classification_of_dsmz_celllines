#!/bin/bash
#SBATCH --job-name=rbl_purity_analysis
#SBATCH --chdir=/home/chu25/TCGA/rbl
#SBATCH --output=/home/chu25/TCGA/rbl/logs/%x_%j.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=168:00:00

# Activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tidyestimate || { echo "Failed to activate Conda environment tidyestimate"; exit 1; }

# Ensure input R script exists
R_SCRIPT="/home/chu25/TCGA/rbl/tumour_purity_estimation.R"
if [[ ! -f "$R_SCRIPT" ]]; then
    echo "Error: R script $R_SCRIPT not found"
    exit 1
fi

echo "Running R script $R_SCRIPT"

Rscript "$R_SCRIPT" || { echo 'R script execution failed'; exit 1; }

echo "Job completed successfully."