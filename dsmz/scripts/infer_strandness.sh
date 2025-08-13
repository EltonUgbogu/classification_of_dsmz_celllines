#!/bin/bash
#SBATCH --job-name=infer_strandness                            # Job name
#SBATCH --output=/home/chu25/dsmz/logs/strandness/infer_strandness_%j.out  # STDOUT log file (%j = job ID)
#SBATCH --error=/home/chu25/dsmz/logs/strandness/infer_strandness_%j.err   # STDERR log file
#SBATCH --time=24:00:00                                        # Max runtime
#SBATCH --cpus-per-task=4                                      # Number of CPUs
#SBATCH --mem=8G                                               # Memory allocation

# -------------------------------------------------------------------------------
# SLURM wrapper script to run strandness inference using RSeQC's infer_experiment.py
# through your batch Python script: run_infer_experiment.py
# -------------------------------------------------------------------------------

# Exit immediately on error (-e), undefined variable (-u), or pipe error (-o pipefail)
set -euo pipefail

# === Load Conda environment ===
echo "[INFO] Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh || {
    echo "[ERROR] Cannot source Conda profile script"
    exit 1
}

conda activate rseqc || {
    echo "[ERROR] Failed to activate 'rseqc' environment"
    exit 1
}

# === Define script path ===
SCRIPT_PATH="/home/chu25/dsmz/scripts/run_infer_experiment.py"

# Ensure the Python script exists before running
if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "[ERROR] Python script not found: $SCRIPT_PATH"
    exit 1
fi

# === Define and create log directory ===
LOG_DIR="/home/chu25/dsmz/logs/strandness"
mkdir -p "$LOG_DIR"

# === Run the script ===
echo "[INFO] Starting strandness inference using: $SCRIPT_PATH"
python3 "$SCRIPT_PATH"

# === Completion message ===
echo "[SUCCESS] Strandness inference completed for all BAM files in batch."
