#!/bin/bash
# SLURM job script: Merge GEO and TARGET NBL data
#SBATCH --job-name=nbl_merge_geo_target
#SBATCH --output=/home/chu25/TCGA/nbl/logs/merge_geo_target_%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/merge_geo_target_%j.err
#SBATCH --ntasks=1
#SBATCH --time=04:00:00           # 4 hours needed for merging
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G                 # 64GB is more than enough for merging

# ==============================
# Activate conda environment
# ==============================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate parquet-r || { echo "Failed to activate 'parquet-r'"; exit 1; }   

# ==============================
# Set up paths
# ==============================
SCRIPT_DIR="/home/chu25/TCGA/nbl/preprocessing"
LOG_DIR="${SCRIPT_DIR}/logs"
RESULTS_DIR="${SCRIPT_DIR}/results"
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

# ==============================
# Run merge GEO and TARGET NBL data
# ==============================
echo "Starting merge GEO and TARGET NBL data..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

Rscript merge_geo_target.R

STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "Merge GEO and TARGET NBL data completed successfully!"
else
    echo "ERROR: merge_geo_target.R failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS