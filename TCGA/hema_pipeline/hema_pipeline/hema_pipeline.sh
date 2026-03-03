#!/bin/bash
#SBATCH --job-name=hema_pipeline
#SBATCH --output=/home/chu25/TCGA/hema_pipeline/hema_pipeline/logs/hema_pipeline_%j.out
#SBATCH --error=/home/chu25/TCGA/hema_pipeline/hema_pipeline/logs/hema_pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# ------------------------------------------------------------------
# 1. Setup: conda + project directory
# ------------------------------------------------------------------
# Ensure logs directory exists
mkdir -p /home/chu25/TCGA/hema_pipeline/hema_pipeline/logs

# Load conda
if [ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    echo "Error: Conda profile script not found at $HOME/miniconda3/etc/profile.d/conda.sh"
    exit 1
fi
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate environment with Snakemake
conda activate snakemake || { echo "Failed to activate 'snakemake' env"; exit 1; }

# Go to pipeline directory
cd /home/chu25/TCGA/hema_pipeline/hema_pipeline || { echo "Failed to change to project directory"; exit 1; }

# Use classic conda solver for --use-conda
export CONDA_FRONTEND=conda

echo "Running from directory: $PWD"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK}"

# ------------------------------------------------------------------
# 2. Unlock any previous Snakemake run
# ------------------------------------------------------------------
echo "Unlocking Snakemake working directory (if locked)..."
snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --unlock || true

# ------------------------------------------------------------------
# 3. Run Snakemake pipeline
# ------------------------------------------------------------------
echo "Starting Snakemake pipeline..."

snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --use-conda --conda-frontend "$CONDA_FRONTEND" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going

STATUS=$?
echo "Snakemake pipeline completed with status: $STATUS"
exit $STATUS
