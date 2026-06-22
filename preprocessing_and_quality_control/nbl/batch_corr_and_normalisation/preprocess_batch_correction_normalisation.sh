#!/bin/bash
#SBATCH --job-name=nbl_batch_corr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --requeue
# NOTE: SLURM log paths are set at submit time.  Submit from the project
# directory so that logs/ resolves correctly, or pass absolute paths:
#   sbatch --output=/abs/path/logs/slurm-%j.out --error=... this_script.sh
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -euo pipefail

# ------------------------------------------------------------------
# 0. Paths (single source of truth - auto-detected)
# ------------------------------------------------------------------
# SCRIPT_DIR is always the directory containing this file, regardless of
# where sbatch was invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"
# nbl/ is one level up, preprocessing_and_quality_control/ is two levels up
NBL_DIR="$(cd "$PROJECT_DIR/.." && pwd)"
PREPROC_ROOT="$(cd "$PROJECT_DIR/../.." && pwd)"

SNAKEFILE="${SNAKEFILE:-$PROJECT_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$PROJECT_DIR/config/config.yaml}"
LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"

# Auto-detect envs directory containing smk.yaml (Snakemake runner env).
# Search order: env var override → sibling preprocessing envs → pipeline root.
for candidate in "${ENVS_DIR:-}" \
                 "$PREPROC_ROOT/envs" \
                 "$NBL_DIR/envs" \
                 "$PROJECT_DIR/envs" \
                 "$REPO_ROOT/envs"
do
  if [ -n "$candidate" ] && [ -f "$candidate/smk.yaml" ]; then
    ENVS_DIR="$candidate"
    break
  fi
done
ENVS_DIR="${ENVS_DIR:-$PREPROC_ROOT/envs}"
SMK_ENV_YAML="${SMK_ENV_YAML:-$ENVS_DIR/smk.yaml}"

# Ensure log directories required by the Snakefile rules exist.
mkdir -p "$LOG_DIR" "$PROJECT_DIR/logs/nbl"

echo "[INFO] Running from directory: $PROJECT_DIR"
echo "[INFO] Snakefile: $SNAKEFILE"
echo "[INFO] Config:    $CONFIGFILE"
echo "[INFO] Logs:      $LOG_DIR"
echo "[INFO] ENVS_DIR:  $ENVS_DIR"
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
# Useful on restricted compute nodes: set SKIP_ENV_CREATE=1 to avoid network calls.
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
  --use-conda \
  --unlock || true

# ------------------------------------------------------------------
# 3. Run Snakemake pipeline
# ------------------------------------------------------------------
echo "[INFO] Starting Snakemake pipeline..."
N_CORES="${SLURM_CPUS_PER_TASK:-8}"

set +e
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --cores "$N_CORES" \
  --use-conda \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  all
STATUS=$?
set -e

echo "[INFO] Snakemake pipeline completed with status: $STATUS"
exit $STATUS
