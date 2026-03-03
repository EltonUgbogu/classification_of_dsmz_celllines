#!/bin/bash
# SLURM job script to run RBL cohort count matrix aggregation
#SBATCH --job-name=rbl_cohort_agg
#SBATCH --output=/home/chu25/TCGA/rbl/logs/rbl_cohort_agg_%j.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/rbl_cohort_agg_%j.err
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
SCRIPT_DIR="/home/chu25/TCGA/rbl"
SCRIPT="rbl_cohort_aggregator.py"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"

cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

# ==============================
# Run the aggregator
# ==============================
echo "Starting RBL cohort aggregation (rbl_cohort_aggregator.py)..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

python "$SCRIPT"

STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "RBL cohort aggregation completed successfully!"
else
    echo "ERROR: rbl_cohort_aggregator.py failed with exit code $STATUS"
fi

echo "End time: $(date)"
exit $STATUS