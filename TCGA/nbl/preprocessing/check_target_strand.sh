#!/bin/bash
# SLURM job script: Check for strandedness in TARGET NBL cohort
#SBATCH --job-name=nbl_target_strand_check
#SBATCH --output=/home/chu25/TCGA/nbl/logs/nbl_target_strand_check_%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/nbl_target_strand_check_%j.err
#SBATCH --ntasks=1
#SBATCH --time=01:00:00           # Only 1 hour needed
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
# Run NBL target strand check
# ==============================
echo "Starting NBL target strand check..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

Rscript check_target_nbl_strand.R

STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "NBL target strand check completed successfully!"
else
    echo "ERROR: check_target_nbl_strand.R failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS