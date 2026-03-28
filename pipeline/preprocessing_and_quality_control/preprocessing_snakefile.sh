#!/bin/bash
#SBATCH --job-name=preprocessing_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --requeue
# Note: --output and --error can be overridden at submission:
#   sbatch --output=/path/to/logs/%j.out --error=/path/to/logs/%j.err preprocessing_snakefile.sh

set -euo pipefail

# ------------------------------------------------------------------
# 0. Paths (single source of truth - auto-detected)
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"
PIPELINE_DIR="$(cd "$PROJECT_DIR/.." && pwd)"

SNAKEFILE="${SNAKEFILE:-$PROJECT_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$PROJECT_DIR/config/config.yaml}"
LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"

# Auto-detect envs directory and tolerate stale ENVS_DIR exports.
for candidate in "${ENVS_DIR:-}" \
                 "$PROJECT_DIR/envs" \
                 "$PIPELINE_DIR/envs" \
                 "/work/ugbogu/pipeline/envs"
do
  if [ -n "$candidate" ] && [ -f "$candidate/smk.yaml" ]; then
    ENVS_DIR="$candidate"
    break
  fi
done
ENVS_DIR="${ENVS_DIR:-$PROJECT_DIR/envs}"
SMK_ENV_YAML="${SMK_ENV_YAML:-$ENVS_DIR/smk.yaml}"

# Always ensure logs directory exists in the project
mkdir -p "$LOG_DIR"

echo "[INFO] Running from directory: $PROJECT_DIR"
echo "[INFO] Snakefile: $SNAKEFILE"
echo "[INFO] Config:    $CONFIGFILE"
echo "[INFO] Logs:      $LOG_DIR"
echo "[INFO] ENVS_DIR:  $ENVS_DIR"
echo "[INFO] SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-8}"
echo "[INFO] SLURM_JOB_ID = ${SLURM_JOB_ID:-NA}"

cd "$PROJECT_DIR" || { echo "[ERROR] Failed to cd into: $PROJECT_DIR"; exit 1; }

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

# Load conda shell function
# shellcheck disable=SC1091
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  echo "[ERROR] Could not find conda.sh at: $(conda info --base)/etc/profile.d/conda.sh"
  exit 1
fi

# Create or update the env at envs/.conda/smk
if [ ! -d "$SMK_ENV_PATH" ]; then
  echo "[INFO] Creating Snakemake env at: $SMK_ENV_PATH"
  conda env create -f "$SMK_ENV_YAML" -p "$SMK_ENV_PATH"
else
  echo "[INFO] Snakemake env exists: $SMK_ENV_PATH"
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
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  all
STATUS=$?
set -e

echo "[INFO] Snakemake pipeline completed with status: $STATUS"
exit $STATUS
