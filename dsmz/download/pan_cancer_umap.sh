#!/bin/bash

# SLURM directives: Define job parameters for resource allocation and logging.
#SBATCH --job-name=pan_cancer_umap
#SBATCH --output=/home/chu25/colDataAnnData/logs/pan_cancer_umap_%j.out
#SBATCH --error=/home/chu25/colDataAnnData/logs/pan_cancer_umap_%j.err
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    # Memory limit: 64GB (sufficient for R analysis with large datasets).

# Initialize Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env || { echo "Failed to activate 'tcga-r-env' env"; exit 1; }


# Change to the scripts directory
cd /home/chu25/colDataAnnData/scripts || { echo "Failed to change to scripts directory"; exit 1; }

# Run R script
echo "Starting pan_cancer_umap processing..."
Rscript pan_cancer_umap.R
STATUS=$?
exit $STATUS

