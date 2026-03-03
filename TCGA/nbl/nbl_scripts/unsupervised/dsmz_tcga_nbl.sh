#!/bin/bash

# SLURM directives: Define job parameters for resource allocation and logging.
#SBATCH --job-name=dsmz_tcga_nbl
#SBATCH --output=/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/logs/dsmz_tcga_nbl_%j.out
#SBATCH --error=/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/logs/dsmz_tcga_nbl_%j.err
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    # Memory limit: 64GB (sufficient for R analysis with large datasets).

# Initialize Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env || { echo "Failed to activate 'tcga-r-env' env"; exit 1; }


# Change to the scripts directory
cd /home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised || { echo "Failed to change to scripts directory"; exit 1; }

# Run R script
echo "Starting DSMZ-TCGA Breast Cancer analysis..."
Rscript nbl_main.R
STATUS=$?
exit $STATUS

