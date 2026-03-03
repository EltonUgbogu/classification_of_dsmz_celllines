#!/bin/bash
# SLURM job script: Check for strand-specificity batch effects in NBL cohort
#SBATCH --job-name=nbl_strand_check
#SBATCH --output=/home/chu25/TCGA/nbl/logs/strand_batch_check_%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/strand_batch_check_%j.err
#SBATCH --ntasks=1
#SBATCH --time=02:00:00           # Only 2 hours needed (fast heatmap)
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G                 # 32GB is more than enough for VST matrix + pheatmap

# ==============================
# Activate conda environment
# ==============================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env || { echo "Failed to activate 'tcga-r-env'"; exit 1; }   

# ==============================
# Set up paths
# ==============================
SCRIPT_DIR="/home/chu25/TCGA/nbl/preprocessing"
LOG_DIR="${SCRIPT_DIR}/logs"
RESULTS_DIR="${SCRIPT_DIR}/results"
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

# ==============================
# Run strand batch effect check
# ==============================
echo "Starting strand batch effect check..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

Rscript strand_batch_effect.R

STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "Strand batch effect check completed successfully!"
else
    echo "ERROR: strand_batch_effect.R failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS