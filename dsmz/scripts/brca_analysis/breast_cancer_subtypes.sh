#!/bin/bash

# Purpose: SLURM job submission script for executing breast cancer subtypes analysis on HPC.
# This script configures resources, initializes the environment, runs the R analysis, and dispatches a GitHub event
# on completion (success or failure) to trigger CI/CD workflows (e.g., result sync).
# Usage: sbatch script.sh # Submits the job; monitor with squeue -u $USER.
# Prerequisites:
#   - SLURM scheduler on HPC.
#   - Conda installed (e.g., Miniconda3 at ~/miniconda3);
#   - tcga-r-env conda environment with required R packages (tidyverse, SummarizedExperiment, etc.)
#   - ~/.secrets/github_token exists with a valid GitHub PAT (Personal Access Token) for repo dispatches.
#   - actions.sh sourced or GitHub vars exported for dispatch (GH_OWNER, GH_REPO, GH_TOKEN).
# Security Note: Token is loaded securely; ensure file permissions are 600. Never echo or log the token.

# SLURM directives: Define job parameters for resource allocation and logging.
#SBATCH --job-name=brca_subtypes
#SBATCH --output=/home/chu25/dsmz/logs/brca_subtypes_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/brca_subtypes_%j.err
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    # Memory limit: 64GB (sufficient for R analysis with large datasets).

# Initialize Conda environment for R
source ~/miniconda3/etc/profile.d/conda.sh  # Initializes Conda (creates 'conda' command).
conda activate snakemake || { echo "Failed to activate 'snakemake' env"; exit 1; }

# Source GitHub environment setup for dispatch variables.
# actions.sh loads GH_OWNER, GH_REPO, GH_TOKEN securely.
source /home/chu25/actions.sh

# Optional: Export additional context for enriched event payloads.
export JOB_LOG="$PWD/slurm-$SLURM_JOB_ID.out"  # Path to job log for reference in dispatches.
export HPC_LABEL="brca_subtypes_analysis"


# Change to the project directory and run the Snakemake pipeline.
cd /home/chu25/dsmz/scripts/brca_analysis/ # Update to actual project path

# Default conda frontend if not provided
: "${CONDA_FRONTEND:=conda}"

# Run Snakemake with robust options for HPC environments
  snakemake \
  --snakefile breast_cancer_subtypes.smk \
  --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --use-conda --conda-frontend "$CONDA_FRONTEND" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going

# Create output directories
OUT_DIR="/home/chu25/dsmz/results/brca"
LOG_DIR="/home/chu25/dsmz/logs"
mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

# Debug: Show current environment and paths
echo "Current directory: $(pwd)"
echo "Conda environment: $CONDA_DEFAULT_ENV"

echo "Snakemake analysis completed with exit status: $?"

# Capture exit status and trigger GitHub repository_dispatch event.
STATUS=$?  # Stores 0 (success) or non-zero (failure) from Snakemake analysis.
JOB_ID=$SLURM_JOB_ID  # SLURM job ID for payload.

# Dispatch event to GitHub: Triggers workflows like HPC Job Results Sync.
# Uses conditional bash for event_type and status based on $STATUS.
curl -X POST \
  -H "Authorization: token $GH_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/repos/$GH_OWNER/$GH_REPO/dispatches \
  -d "{\"event_type\":\"$([ $STATUS -eq 0 ] && echo 'slurm-job-complete' || echo 'slurm-job-failed')\",\"client_payload\":{\"status\":\"$([ $STATUS -eq 0 ] && echo 'success' || echo 'failure')\",\"job_id\":\"$JOB_ID\",\"scheduler\":\"slurm\",\"host\":\"$(hostname)\",\"user\":\"$USER\",\"note\":\"Analysis complete\",\"label\":\"$HPC_LABEL\"}}"

# Exit with the original status to propagate success/failure in SLURM accounting.
exit $STATUS