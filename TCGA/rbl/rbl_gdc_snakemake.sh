#!/bin/bash

# Purpose: SLURM job submission script for executing a Snakemake-based RNA-seq pipeline on HPC.
# This script configures resources, initializes the environment, runs the analysis, and dispatches a GitHub event
# on completion (success or failure) to trigger CI/CD workflows (e.g., result sync).

# SLURM directives: Define job parameters for resource allocation and logging.
#SBATCH --job-name=rbl_gdc
#SBATCH --output=/home/chu25/TCGA/rbl/logs/rbl_gdc_%j.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/rbl_gdc_%j.err
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# --- Initialization and Setup ---

# Initialize Conda environment (vanilla Conda; no Mamba). Adjust path if your installation differs.
if [ ! -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    echo "Error: Conda profile script not found at ~/miniconda3/etc/profile.d/conda.sh"
    exit 1
fi
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake || { echo "Failed to activate 'snakemake' env"; exit 1; }

# Source GitHub environment setup for dispatch variables.
if [ ! -f /home/chu25/actions.sh ]; then
    echo "Error: GitHub actions script not found at /home/chu25/actions.sh"
    exit 1
fi
source /home/chu25/actions.sh

# Change to the project directory
cd /home/chu25/TCGA/rbl/ || { echo "Failed to change to project directory"; exit 1; }

# Optional: Export additional context for enriched event payloads.
export JOB_LOG="$PWD/slurm-$SLURM_JOB_ID.out"
export HPC_LABEL="rbl_gdc"
# Force classic conda solver (disable libmamba for this job)
export CONDA_FRONTEND=conda

# --- Snakemake Execution ---

# 1. Ensure the directory is unlocked from any previous failed jobs
# This addresses the LockException seen in previous logs.
echo "Unlocking Snakemake directory..."
snakemake --snakefile rbl_gdc_rnaseq_pipeline.smk --unlock

# 2. Run Snakemake with robust options for HPC environments
# Note: Rerunning incomplete jobs is included in the --keep-going logic on failure,
# but we prioritize --latency-wait to avoid false "missing file" errors.
echo "Starting Snakemake pipeline..."
snakemake \
  --snakefile rbl_gdc_rnaseq_pipeline.smk \
  --configfile config/config.yaml \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --use-conda --conda-frontend "$CONDA_FRONTEND" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going

# Capture exit status and trigger GitHub repository_dispatch event.
STATUS=$?
JOB_ID=$SLURM_JOB_ID
echo "Snakemake pipeline completed with status: $STATUS"

# --- GitHub Dispatch ---
# Uses conditional bash for event_type and status based on $STATUS.
if [ "$STATUS" -eq 0 ]; then
    EVENT_TYPE='slurm-job-complete'
    JOB_STATUS='success'
else
    EVENT_TYPE='slurm-job-failed'
    JOB_STATUS='failure'
fi

echo "Sending GitHub dispatch event..."
curl -X POST \
  -H "Authorization: token $GH_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/repos/$GH_OWNER/$GH_REPO/dispatches \
  -d "{\"event_type\":\"$EVENT_TYPE\",\"client_payload\":{\"status\":\"$JOB_STATUS\",\"job_id\":\"$JOB_ID\",\"scheduler\":\"slurm\",\"host\":\"$(hostname)\",\"user\":\"$USER\",\"note\":\"Alignment and QC complete\",\"label\":\"$HPC_LABEL\"}}" \
  --fail --silent --show-error || echo "Warning: Failed to send GitHub dispatch event"

# Exit with the original status to propagate success/failure in SLURM accounting.
echo "Job completed. Final status: $STATUS"
exit $STATUS
