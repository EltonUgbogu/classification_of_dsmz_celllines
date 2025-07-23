#!/bin/bash

# #SBATCH --job-name=download_all_tcga_star_count_data
# #SBATCH --output=logs/all_tcga_se%j.log
# #SBATCH --ntasks=1
# #SBATCH --time=24:00:00
# #SBATCH --mem=32G

echo "Running R package installation script on HPC"

echo "Start time: $(date)"

# Setup Conda and activate environment
source ~/miniconda3/etc/profile.d/conda.sh
# Activate my conda enviroment
conda activate tcga-r-env

# Verify Rscript is in PATH
which Rscript || { echo "Rscript not found in PATH"; exit 1; }

# Make sure logs folder exists
mkdir -p logs

# Run the R script
Rscript all_tcga_se.r > logs/all_tcga_se.log 2>&1

echo "Finished at: $(date)"
echo "Check logs/all_tcga_se.log for details"
