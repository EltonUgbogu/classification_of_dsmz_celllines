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
# 0. Repository-relative path resolution
# ------------------------------------------------------------------
# Resolve the repository root from the launcher location. Execution is routed
# through the parent NBL workflow so nbl/config/config.yaml remains the single
# active configuration.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NBL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PREPROC_ROOT="$(cd "$NBL_DIR/.." && pwd)"
REPO_ROOT="$(cd "$PREPROC_ROOT/.." && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$NBL_DIR}"

SNAKEFILE="${SNAKEFILE:-$NBL_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$NBL_DIR/config/config.yaml}"
TARGET="${TARGET:-nbl_count_based_preprocessing}"
LOG_DIR="${LOG_DIR:-$NBL_DIR/logs}"

# Resolve the Snakemake runner environment from an explicit override,
# preprocessing-level definitions, or the repository-level environment.
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
echo "[INFO] Target:    $TARGET"
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
# 1. Configurable Snakemake execution environment
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

# Keep the Snakemake runner environment inside the preprocessing workspace
# unless the caller supplies SMK_ENV_PATH.
SMK_ENV_PATH="${SMK_ENV_PATH:-$PROJECT_DIR/envs/.conda/smk}"
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$PROJECT_DIR/.snakemake/conda}"
mkdir -p "$(dirname "$SMK_ENV_PATH")"
SKIP_ENV_CREATE="${SKIP_ENV_CREATE:-0}"

# Initialise conda shell integration before activating the runner environment.
# shellcheck disable=SC1091
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  echo "[ERROR] Could not find conda.sh at: $(conda info --base)/etc/profile.d/conda.sh"
  exit 1
fi

# Create or reuse the configured runner environment. On compute nodes without
# network access, set SKIP_ENV_CREATE=1 to require an existing environment.
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

SNAKEMAKE="${SNAKEMAKE:-${SNAKEMAKE_BIN:-$(command -v snakemake || true)}}"
if [ -z "$SNAKEMAKE" ]; then
  echo "[ERROR] snakemake not found in activated env ($SMK_ENV_PATH)"
  exit 1
fi
echo "[INFO] Using snakemake: $SNAKEMAKE ($("$SNAKEMAKE" --version))"
echo "[INFO] Snakemake conda prefix: $SNAKEMAKE_CONDA_PREFIX"

# ------------------------------------------------------------------
# 2. Clear any stale Snakemake working-directory lock
# ------------------------------------------------------------------
echo "[INFO] Unlocking Snakemake working directory (if locked)..."
"$SNAKEMAKE" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --use-conda \
  --conda-frontend conda \
  --unlock || true

# ------------------------------------------------------------------
# 3. Execute the configured workflow target
# ------------------------------------------------------------------
echo "[INFO] Starting Snakemake pipeline..."
CORES="${CORES:-${N_CORES:-${SLURM_CPUS_PER_TASK:-8}}}"

set +e
"$SNAKEMAKE" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --cores "$CORES" \
  --use-conda \
  --conda-frontend conda \
  --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  "$TARGET"
STATUS=$?
set -e

echo "[INFO] Snakemake pipeline completed with status: $STATUS"
exit $STATUS
