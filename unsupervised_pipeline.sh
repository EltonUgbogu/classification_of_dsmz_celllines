#!/bin/bash
#SBATCH --job-name=unsupervised_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
# Note: --output and --error can be overridden at submission:
#   sbatch --output=/path/to/logs/%j.out --error=/path/to/logs/%j.err unsupervised_pipeline.sh

set -euo pipefail

# ------------------------------------------------------------------
# 0. Paths (single source of truth - auto-detected)
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"

SNAKEFILE="${SNAKEFILE:-$PROJECT_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$PROJECT_DIR/config/config.yaml}"
LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"
PIPELINE_PROFILE="${PIPELINE_PROFILE:-}"
PIPELINE_UNLOCK="${PIPELINE_UNLOCK:-0}"
PIPELINE_NOLOCK="${PIPELINE_NOLOCK:-1}"

# Prefer the explicit positional profile supplied to this wrapper over any
# exported PIPELINE_PROFILE inherited by sbatch from the submit shell.
if [[ $# -gt 0 ]]; then
  PIPELINE_PROFILE="$1"
  shift
fi

if [[ -z "$PIPELINE_PROFILE" ]]; then
  echo "[ERROR] Missing pipeline profile."
  echo "        Usage: bash unsupervised_pipeline.sh <multicohort_cancer|brca|nbl|rbl|heme> [snakemake targets...]"
  echo "        Or set PIPELINE_PROFILE in the environment."
  exit 1
fi

TARGETS=("$@")
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=(all)
fi

ENVS_DIR="${ENVS_DIR:-$PROJECT_DIR/envs}"
SMK_ENV_YAML="${SMK_ENV_YAML:-$ENVS_DIR/smk.yaml}"

# Always ensure logs directory exists in the project
mkdir -p "$LOG_DIR"
RUN_LOCK_DIR="${RUN_LOCK_DIR:-$LOG_DIR/run_locks}"
mkdir -p "$RUN_LOCK_DIR"

echo "[INFO] Running from directory: $PROJECT_DIR"
echo "[INFO] Snakefile: $SNAKEFILE"
echo "[INFO] Config:    $CONFIGFILE"
echo "[INFO] Logs:      $LOG_DIR"
echo "[INFO] ENVS_DIR:  $ENVS_DIR"
echo "[INFO] PIPELINE_PROFILE = ${PIPELINE_PROFILE:-NA}"
echo "[INFO] TARGETS = ${TARGETS[*]}"
echo "[INFO] PIPELINE_UNLOCK = $PIPELINE_UNLOCK"
echo "[INFO] PIPELINE_NOLOCK = $PIPELINE_NOLOCK"
echo "[INFO] SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-8}"
echo "[INFO] SLURM_JOB_ID = ${SLURM_JOB_ID:-NA}"

cd "$PROJECT_DIR" || { echo "[ERROR] Failed to cd into: $PROJECT_DIR"; exit 1; }

# Profile-level lock. Different profiles may run concurrently, but a second
# run of the same profile would target the same output tree and is blocked.
LOCK_SAFE_PROFILE="$(printf '%s' "$PIPELINE_PROFILE" | tr -c 'A-Za-z0-9_.-' '_')"
RUN_LOCK_FILE="$RUN_LOCK_DIR/${LOCK_SAFE_PROFILE}.lock"
exec 9>"$RUN_LOCK_FILE"
if ! flock -n 9; then
  echo "[ERROR] Another run for profile '$PIPELINE_PROFILE' is already active."
  echo "        Lock file: $RUN_LOCK_FILE"
  echo "        Different profiles can run in parallel; duplicate profile runs cannot."
  exit 1
fi
echo "[INFO] Acquired profile run lock: $RUN_LOCK_FILE"

# ------------------------------------------------------------------
# 1. Conda-only setup + create/activate Snakemake env from ./envs/
# ------------------------------------------------------------------
# Bootstrap conda for non-interactive HPC batch shells.
if ! command -v conda >/dev/null 2>&1; then
  for conda_sh in \
    "${CONDA_SH_PATH:-}" \
    "${CONDA_BASE:-}/etc/profile.d/conda.sh" \
    "${HOME}/miniforge3/etc/profile.d/conda.sh" \
    "${HOME}/miniconda3/etc/profile.d/conda.sh" \
    "${HOME}/anaconda3/etc/profile.d/conda.sh"
  do
    if [ -n "${conda_sh}" ] && [ -f "${conda_sh}" ]; then
      # shellcheck disable=SC1090
      source "${conda_sh}"
      break
    fi
  done
fi

if ! command -v conda >/dev/null 2>&1 && [ -n "${CONDA_EXE:-}" ]; then
  CONDA_BASE_FROM_EXE="$(cd "$(dirname "${CONDA_EXE}")/.." && pwd)"
  if [ -f "${CONDA_BASE_FROM_EXE}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_BASE_FROM_EXE}/etc/profile.d/conda.sh"
  fi
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] 'conda' not found in PATH."
  echo "        Set CONDA_SH_PATH/CONDA_BASE, or load your site conda module before sbatch."
  exit 1
fi

if [ ! -f "$SMK_ENV_YAML" ]; then
  echo "[ERROR] Snakemake environment YAML not found: $SMK_ENV_YAML"
  echo "        Expected location: $ENVS_DIR/smk.yaml"
  exit 1
fi

# Create a local env *inside the project* (portable)
SMK_ENV_PATH="${SMK_ENV_PATH:-$ENVS_DIR/.conda/smk}"
mkdir -p "$(dirname "$SMK_ENV_PATH")"

# Load conda shell function after bootstrap.
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
  echo "[INFO] Updating env from: $SMK_ENV_YAML"
  conda env update -f "$SMK_ENV_YAML" -p "$SMK_ENV_PATH" || {
    echo "[WARNING] Env update failed; continuing with existing env."
  }
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
# 2. Optional unlock
# ------------------------------------------------------------------
if [[ "$PIPELINE_UNLOCK" == "1" || "$PIPELINE_UNLOCK" == "true" ]]; then
  echo "[INFO] Unlocking Snakemake working directory by explicit request..."
  "$SNAKEMAKE_BIN" \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    ${PIPELINE_PROFILE:+--config pipeline_profile="$PIPELINE_PROFILE"} \
    --use-conda \
    --unlock || true
else
  echo "[INFO] Skipping automatic Snakemake unlock. Set PIPELINE_UNLOCK=1 only after confirming no runs are active."
fi

# ------------------------------------------------------------------
# 3. Run Snakemake pipeline
# ------------------------------------------------------------------
echo "[INFO] Starting Snakemake pipeline..."
N_CORES="${SLURM_CPUS_PER_TASK:-8}"
SNAKEMAKE_LOCK_ARGS=()
if [[ "$PIPELINE_NOLOCK" == "1" || "$PIPELINE_NOLOCK" == "true" ]]; then
  SNAKEMAKE_LOCK_ARGS+=(--nolock)
fi

set +e
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  ${PIPELINE_PROFILE:+--config pipeline_profile="$PIPELINE_PROFILE"} \
  --cores "$N_CORES" \
  --use-conda \
  "${SNAKEMAKE_LOCK_ARGS[@]}" \
  --printshellcmds -p \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  "${TARGETS[@]}"
STATUS=$?
set -e

echo "[INFO] Snakemake pipeline completed with status: $STATUS"
exit $STATUS
