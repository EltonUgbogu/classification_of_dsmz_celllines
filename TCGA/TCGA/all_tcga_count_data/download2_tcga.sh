#!/bin/bash
#SBATCH --job-name=tcga_count_data
#SBATCH --output=/home/chu25/TCGA/logs/tcga_count_data_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/tcga_count_data_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00

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
Rscript tcga_all_data.r > /home/chu25/TCGA/logs/tcga_all_data.log 2>&1

echo "Finished at: $(date)"
echo "Check logs/tcga_all_data.log for details"
