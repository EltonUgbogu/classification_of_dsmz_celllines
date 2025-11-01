#!/bin/bash
#
#SBATCH --job-name=dsmz_tcga_brca_subtypes
#SBATCH --output=/home/chu25/dsmz/scripts/brca_analysis/logs/dsmz_tcga_brca_subtypes3_%j.out
#SBATCH --error=/home/chu25/dsmz/scripts/brca_analysis/logs/dsmz_tcga_brca_subtypes3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=120:00:00

# --- Conda activation (no nounset here) ---
# don't use `set -u` yet; conda's activate.d hooks expect some vars to be unset
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env || { echo "[FATAL] Could not activate 'tcga-r-env'"; exit 1; }

# Change to the scripts directory
cd /home/chu25/dsmz/scripts/brca_analysis/ || { echo "Failed to change to scripts directory"; exit 1; }

# Run R script
echo "Starting DSMZ-TCGA breast cancer subtypes analysis via R script..."
Rscript brca_analysis3.R
STATUS=$?
exit $STATUS