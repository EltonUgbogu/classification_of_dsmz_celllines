#!/bin/bash
#SBATCH --job-name=nbl_gdc
#SBATCH --output=/home/chu25/TCGA/nbl/logs/%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=120:00:00

set -euo pipefail

# --- paths ---
WORKDIR="/home/chu25/TCGA/nbl"
LOGDIR="$WORKDIR/logs"
SNAKEFILE="$WORKDIR/nbl_gdc_rnaseq_pipeline.smk"
CFG="$WORKDIR/config/config.yaml"

mkdir -p "$LOGDIR"
cd "$WORKDIR"

# --- conda/mamba ---
# Adjust this if your base is elsewhere
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake || { echo "Failed to activate 'snakemake' env"; exit 1; }

# --- run snakemake ---
snakemake \
  --snakefile "$SNAKEFILE" \
  --configfile "$CFG" \
  --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --use-conda --conda-frontend "$CONDA_FRONTEND" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going
