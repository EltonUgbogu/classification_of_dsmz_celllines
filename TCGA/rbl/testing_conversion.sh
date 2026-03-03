#!/bin/bash

# SLURM directives: Define job parameters for resource allocation and logging.
#SBATCH --job-name=testing_conversion
#SBATCH --output=/home/chu25/TCGA/rbl/logs/testing_conversion_%j.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/testing_conversion_%j.err
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    # Memory limit: 64GB (sufficient for R analysis with large datasets).

# Initialize Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env || { echo "Failed to activate 'tcga-r-env' env"; exit 1; }


# Change to the scripts directory
cd /home/chu25/TCGA/rbl || { echo "Failed to change to scripts directory"; exit 1; }

# Run R script
echo "Starting testing conversion..."
Rscript testing_conversion.R
STATUS=$?
exit $STATUS

