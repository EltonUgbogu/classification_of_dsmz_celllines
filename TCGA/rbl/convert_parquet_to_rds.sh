#!/bin/bash
# SLURM job script: Convert parquet to rds  
#SBATCH --job-name=convert_rbl_parquet_to_rds
#SBATCH --output=/home/chu25/TCGA/rbl/logs/convert_rbl_parquet_to_rds_%j.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/convert_rbl_parquet_to_rds_%j.err
#SBATCH --ntasks=1
#SBATCH --time=02:00:00           # Only 2 hours needed (fast conversion)
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G                 # 32GB is more than enough for conversion

# ==============================
# Activate conda environment
# ==============================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate parquet-r || { echo "Failed to activate 'parquet-r'"; exit 1; }	
# ==============================
# Set up paths
# ==============================
SCRIPT_DIR="/home/chu25/TCGA/rbl"
LOG_DIR="${SCRIPT_DIR}/logs"
RESULTS_DIR="${SCRIPT_DIR}/results"
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

# ==============================
# Run convert RBL parquet to rds
# ==============================
echo "Starting convert RBL parquet to rds..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

Rscript convert_parquet_to_rds.R
STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "Convert RBL parquet to rds completed successfully!"
else
    echo "ERROR: convert_rbl_parquet_to_rds.R failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS