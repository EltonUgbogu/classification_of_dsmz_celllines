#!/bin/bash
#SBATCH --job-name=download_gse100148_missing
#SBATCH --chdir=.
#SBATCH --output=download_gse100148_missing_%j.out
#SBATCH --error=download_gse100148_missing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

echo "== $(date) :: START =="
echo "Host: $(hostname)"

ROOT="${NBL_DATA_ROOT:-$REPO_ROOT/data/nbl/new_data}"
export NBL_DATA_ROOT="$ROOT"
LOG_DIR="${ROOT}/logs"
SCRIPT="$SCRIPT_DIR/download_nbl_tumour_sample_srr_ids.py"

mkdir -p "$LOG_DIR"

echo "[INFO] ROOT=$ROOT"
echo "[INFO] LOG_DIR=$LOG_DIR"
echo "[INFO] SCRIPT=$SCRIPT"

CONDA_SH_PATH="${CONDA_SH_PATH:-${HOME}/miniforge3/etc/profile.d/conda.sh}"
if [[ ! -f "$CONDA_SH_PATH" ]] && command -v conda >/dev/null 2>&1; then
    CONDA_SH_PATH="$(conda info --base)/etc/profile.d/conda.sh"
fi
source "$CONDA_SH_PATH"
conda activate sra3 || true

export PATH="$HOME/edirect:$PATH"

echo "[INFO] Tool checks"
which python3
python3 --version

which prefetch
prefetch --version || true

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Script not found: $SCRIPT" >&2
    exit 2
fi

echo "[INFO] Running targeted GSE100148 missing raw-data download"
python3 "$SCRIPT" download-gse100148-missing

echo "[INFO] Checking downloaded SRA files"
find "$ROOT/patient_tumour_recovery" -name "*.sra" -ls || true

echo "[INFO] Manifest:"
cat "$ROOT/patient_tumour_recovery/GSE100148_missing_download_manifest.tsv" || true

echo "== $(date) :: DONE =="
