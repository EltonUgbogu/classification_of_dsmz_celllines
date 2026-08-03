#!/bin/bash
#SBATCH --job-name=rbl_counts_light
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -euo pipefail

# ------------------------------------------------------------------
# 0. Paths (single source of truth - auto-detected)
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"
PIPELINE_DIR="$(cd "$PROJECT_DIR/../.." && pwd)"

SNAKEFILE="${SNAKEFILE:-$PROJECT_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-}"
LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"
DATA_ROOT="${DATA_ROOT:-$REPO_ROOT/data/rbl}"
export DATA_ROOT

if [ -z "$CONFIGFILE" ]; then
  for candidate in \
    "$PROJECT_DIR/config/config.yaml" \
    "$PIPELINE_DIR/preprocessing_and_quality_control/rbl/config/config.yaml" \
    "$PIPELINE_DIR/config/config.yaml"
  do
    if [ -f "$candidate" ]; then
      CONFIGFILE="$candidate"
      break
    fi
  done
fi
export CONFIGFILE

# Auto-detect envs directory and tolerate stale ENVS_DIR exports.
for candidate in "${ENVS_DIR:-}" \
                 "$PROJECT_DIR/envs" \
                 "$PIPELINE_DIR/envs" \
                 "$REPO_ROOT/envs"
do
  if [ -n "$candidate" ] && [ -f "$candidate/smk.yaml" ]; then
    ENVS_DIR="$candidate"
    break
  fi
done
ENVS_DIR="${ENVS_DIR:-$REPO_ROOT/envs}"
SMK_ENV_YAML="${SMK_ENV_YAML:-$ENVS_DIR/smk.yaml}"

# Always ensure logs directory exists in the project
mkdir -p "$LOG_DIR"


echo "[INFO] Running from directory: $PROJECT_DIR"
echo "[INFO] Snakefile: $SNAKEFILE"
echo "[INFO] Config:    $CONFIGFILE"
echo "[INFO] Logs:      $LOG_DIR"
echo "[INFO] ENVS_DIR:  $ENVS_DIR"
echo "[INFO] DATA_ROOT: $DATA_ROOT"
echo "[INFO] SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-8}"
echo "[INFO] SLURM_JOB_ID = ${SLURM_JOB_ID:-NA}"

cd "$PROJECT_DIR" || { echo "[ERROR] Failed to cd into: $PROJECT_DIR"; exit 1; }

if [ ! -f "$SNAKEFILE" ]; then
  echo "[ERROR] Snakefile not found: $SNAKEFILE"
  exit 1
fi

if [ ! -f "$CONFIGFILE" ]; then
  echo "[ERROR] Config file not found: $CONFIGFILE"
  exit 1
fi

# ------------------------------------------------------------------
# 1. Conda-only setup + create/activate Snakemake env from ENVS_DIR
# ------------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] 'conda' not found in PATH."
  echo "        Load/activate conda first (module load / source init)."
  exit 1
fi

if [ ! -f "$SMK_ENV_YAML" ]; then
  echo "[ERROR] Snakemake environment YAML not found: $SMK_ENV_YAML"
  echo "        Expected location: $ENVS_DIR/smk.yaml"
  exit 1
fi

# Create a local runtime env inside preprocessing project (portable)
SMK_ENV_PATH="${SMK_ENV_PATH:-$PROJECT_DIR/envs/.conda/smk}"
mkdir -p "$(dirname "$SMK_ENV_PATH")"
SKIP_ENV_CREATE="${SKIP_ENV_CREATE:-0}"

# Load conda shell function
# shellcheck disable=SC1091
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  echo "[ERROR] Could not find conda.sh at: $(conda info --base)/etc/profile.d/conda.sh"
  exit 1
fi

# Create or reuse env at envs/.conda/smk.
# Useful on limited compute nodes: set SKIP_ENV_CREATE=1 to avoid network calls.
if [ "$SKIP_ENV_CREATE" = "1" ]; then
  if [ -d "$SMK_ENV_PATH" ]; then
    echo "[INFO] SKIP_ENV_CREATE=1; reusing existing env: $SMK_ENV_PATH"
  else
    echo "[ERROR] SKIP_ENV_CREATE=1 but env does not exist: $SMK_ENV_PATH"
    echo "        Create it first on a node with network access, then resubmit."
    exit 1
  fi
else
  if [ ! -d "$SMK_ENV_PATH" ]; then
    echo "[INFO] Creating Snakemake env at: $SMK_ENV_PATH"
    conda env create -f "$SMK_ENV_YAML" -p "$SMK_ENV_PATH"
  else
    echo "[INFO] Snakemake env exists: $SMK_ENV_PATH"
  fi
fi

echo "[INFO] Activating Snakemake env: $SMK_ENV_PATH"
conda activate "$SMK_ENV_PATH" || {
  echo "[ERROR] Failed to activate conda env: $SMK_ENV_PATH"
  exit 1
}

SNAKEMAKE_BIN="$(command -v snakemake || true)"
if [ -z "$SNAKEMAKE_BIN" ]; then
  echo "[ERROR] snakemake not found in activated env ($SMK_ENV_PATH)"
  exit 1
fi
echo "[INFO] Using snakemake: $SNAKEMAKE_BIN ($(snakemake --version))"

# ------------------------------------------------------------------
# 2. Unlock any previous Snakemake run (safe if not locked)
# ------------------------------------------------------------------
echo "[INFO] Unlocking Snakemake working directory (if locked)..."
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --config data_root="$DATA_ROOT" \
  --use-conda \
  --conda-frontend conda \
  --unlock || true

# ------------------------------------------------------------------
# 3. Run RBL tumour count generation only
# ------------------------------------------------------------------
echo "[INFO] Starting RBL tumour count generation..."
N_CORES="${SLURM_CPUS_PER_TASK:-8}"

set +e
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --config data_root="$DATA_ROOT" \
  --cores "$N_CORES" \
  --use-conda \
  --conda-frontend conda \
  --printshellcmds -p \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --latency-wait 300 \
  "$DATA_ROOT/count_data/rbl_tumour_count.rds"
STATUS=$?
set -e

echo "[INFO] RBL tumour count generation completed with status: $STATUS"
exit $STATUS
