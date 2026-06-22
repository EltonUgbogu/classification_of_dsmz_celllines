#!/bin/bash
#SBATCH --job-name=rbl_batch_corr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --output=../logs/slurm-%j.out
#SBATCH --error=../logs/slurm-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
RBL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PREPROC_ROOT="$(cd "$RBL_DIR/.." && pwd)"

SNAKEFILE="${SNAKEFILE:-$RBL_DIR/Snakefile}"
CONFIGFILE="${CONFIGFILE:-$RBL_DIR/config/config.yaml}"
LOG_DIR="${LOG_DIR:-$RBL_DIR/logs}"
TARGET="${TARGET:-preprocess_integrate_rbl}"

for candidate in "${ENVS_DIR:-}" \
                 "$PREPROC_ROOT/envs" \
                 "$RBL_DIR/envs" \
                 "$SCRIPT_DIR/envs" \
                 "$REPO_ROOT/envs"
do
  if [ -n "$candidate" ] && [ -f "$candidate/smk.yaml" ]; then
    ENVS_DIR="$candidate"
    break
  fi
done
ENVS_DIR="${ENVS_DIR:-$PREPROC_ROOT/envs}"
SMK_ENV_YAML="${SMK_ENV_YAML:-$ENVS_DIR/smk.yaml}"
SMK_ENV_PATH="${SMK_ENV_PATH:-$RBL_DIR/envs/.conda/smk}"
SKIP_ENV_CREATE="${SKIP_ENV_CREATE:-0}"

mkdir -p "$LOG_DIR" "$RBL_DIR/logs/rbl" "$(dirname "$SMK_ENV_PATH")"

echo "[INFO] RBL directory: $RBL_DIR"
echo "[INFO] Snakefile: $SNAKEFILE"
echo "[INFO] Config:    $CONFIGFILE"
echo "[INFO] Target:    $TARGET"
echo "[INFO] Logs:      $LOG_DIR"
echo "[INFO] ENVS_DIR:  $ENVS_DIR"
echo "[INFO] SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-8}"
echo "[INFO] SLURM_JOB_ID = ${SLURM_JOB_ID:-NA}"

cd "$RBL_DIR"

if [ ! -f "$SNAKEFILE" ]; then
  echo "[ERROR] Snakefile not found: $SNAKEFILE"
  false
fi

if [ ! -f "$CONFIGFILE" ]; then
  echo "[ERROR] Config file not found: $CONFIGFILE"
  false
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] 'conda' not found in PATH."
  echo "        Load or initialise conda before submitting this launcher."
  false
fi

if [ ! -f "$SMK_ENV_YAML" ]; then
  echo "[ERROR] Snakemake environment YAML not found: $SMK_ENV_YAML"
  false
fi

# shellcheck disable=SC1091
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  echo "[ERROR] Could not find conda.sh at: $(conda info --base)/etc/profile.d/conda.sh"
  false
fi

if [ "$SKIP_ENV_CREATE" = "1" ]; then
  if [ -d "$SMK_ENV_PATH" ]; then
    echo "[INFO] SKIP_ENV_CREATE=1; reusing env: $SMK_ENV_PATH"
  else
    echo "[ERROR] SKIP_ENV_CREATE=1 but env does not exist: $SMK_ENV_PATH"
    false
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
conda activate "$SMK_ENV_PATH"

SNAKEMAKE_BIN="$(command -v snakemake || true)"
if [ -z "$SNAKEMAKE_BIN" ]; then
  echo "[ERROR] snakemake not found in activated env ($SMK_ENV_PATH)"
  false
fi
echo "[INFO] Using snakemake: $SNAKEMAKE_BIN ($("$SNAKEMAKE_BIN" --version))"

echo "[INFO] Unlocking Snakemake working directory if needed..."
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --use-conda \
  --unlock || true

echo "[INFO] Starting Snakemake target: $TARGET"
N_CORES="${SLURM_CPUS_PER_TASK:-8}"
"$SNAKEMAKE_BIN" \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --cores "$N_CORES" \
  --use-conda \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 300 \
  --keep-going \
  "$TARGET"

echo "[INFO] Snakemake target completed: $TARGET"
