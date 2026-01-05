#!/bin/bash
#SBATCH --job-name=unsupervised_pipeline
#SBATCH --output=/work/ugbogu/pipeline/logs/unsupervised_pipeline_%j.out
#SBATCH --error=/work/ugbogu/pipeline/logs/unsupervised_pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/work/ugbogu/pipeline
#SBATCH --requeue

# Optional: set these if your cluster needs them (leave commented otherwise)
# Override walltime at submission: sbatch --time=24:00:00 unsupervised_pipeline.sh
##SBATCH --partition=cpu
##SBATCH --qos=normal
##SBATCH --account=YOUR_ACCOUNT

set -euo pipefail

# ------------------------------------------------------------------
# 0. Paths (single source of truth)
# ------------------------------------------------------------------
# Prefer the directory you submitted from (correct on clusters),
# otherwise fall back to current working directory (set by --chdir above)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"

SNAKEFILE="${SNAKEFILE:-$PROJECT_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$PROJECT_DIR/config/config.yaml}"
LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"

mkdir -p "$LOG_DIR"

echo "Running from directory: $PROJECT_DIR"
echo "Snakefile: $SNAKEFILE"
echo "Config:    $CONFIGFILE"
echo "Logs:      $LOG_DIR"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-8}"
echo "SLURM_JOB_ID = ${SLURM_JOB_ID:-NA}"

cd "$PROJECT_DIR" || { echo "Failed to cd into: $PROJECT_DIR"; exit 1; }

# ------------------------------------------------------------------
# 1. Snakemake executable (no conda activate; no --use-conda)
# ------------------------------------------------------------------
SMK_ENV_PREFIX="${SMK_ENV_PREFIX:-/work/ugbogu/.conda/envs/smk}"
SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-$SMK_ENV_PREFIX/bin/snakemake}"

if [ ! -x "$SNAKEMAKE_BIN" ]; then
  echo "[ERROR] snakemake not found/executable at: $SNAKEMAKE_BIN"
  echo "        Set SMK_ENV_PREFIX or SNAKEMAKE_BIN to the correct path."
  exit 1
fi
echo "[INFO] Using snakemake: $SNAKEMAKE_BIN"

# ------------------------------------------------------------------
# 2. Unlock any previous Snakemake run (safe if not locked)
# ------------------------------------------------------------------
echo "Unlocking Snakemake working directory (if locked)..."
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --unlock || true

# ------------------------------------------------------------------
# 3. Run Snakemake pipeline
# ------------------------------------------------------------------
echo "Starting Snakemake pipeline..."

N_CORES="${SLURM_CPUS_PER_TASK:-8}"

set +e
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --cores "$N_CORES" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  all
STATUS=$?
set -e

echo "Snakemake pipeline completed with status: $STATUS"
exit $STATUS

