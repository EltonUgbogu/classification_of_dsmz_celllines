#!/bin/bash
# SLURM job script to run NBL cohort count matrix aggregation
#SBATCH --job-name=nbl_cohort_agg
#SBATCH --output=/home/chu25/TCGA/nbl/logs/nbl_cohort_agg_%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/nbl_cohort_agg_%j.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# ==============================
# Load modules & activate conda
# ==============================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv || { echo "Failed to activate 'myenv'"; exit 1; }

# ==============================
# Set up paths
# ==============================
SCRIPT_DIR="/home/chu25/TCGA/nbl/preprocessing"
SCRIPT="nbl_cohort_aggregator.py"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"

cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

# ==============================
# Run the aggregator
# ==============================
echo "Starting NBL cohort aggregation (nbl_cohort_aggregator.py)..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

python "$SCRIPT"

STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "NBL cohort aggregation completed successfully!"
else
    echo "ERROR: nbl_cohort_aggregator.py failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS