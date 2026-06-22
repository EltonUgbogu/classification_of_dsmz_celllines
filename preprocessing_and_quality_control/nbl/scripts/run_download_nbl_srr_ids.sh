#!/bin/bash
#SBATCH --job-name=run_download_nbl_srr_ids
#SBATCH --chdir=.
#SBATCH --output=run_download_nbl_srr_ids_%j.out
#SBATCH --error=run_download_nbl_srr_ids_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

echo "== $(date) :: START =="
echo "Host: $(hostname)"

ROOT="${NBL_DATA_ROOT:-$REPO_ROOT/data/nbl}"
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
which esearch
which efetch

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Script not found: $SCRIPT" >&2
    exit 2
fi

echo "[INFO] Running SRR extraction for all NBL datasets"
python3 "$SCRIPT" all

echo "== $(date) :: DONE =="
