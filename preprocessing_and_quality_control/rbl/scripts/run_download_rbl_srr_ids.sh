#!/bin/bash
#SBATCH --job-name=extract_srr_ids
#SBATCH --output=extract_srr_%j.out
#SBATCH --error=extract_srr_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

echo "== $(date) :: START =="
echo "Host: $(hostname)"

SCRIPT="$SCRIPT_DIR/download_rbl_tumour_sample_srr_ids.py"
export RBL_DATA_ROOT="${RBL_DATA_ROOT:-$REPO_ROOT/data/rbl}"

CONDA_SH_PATH="${CONDA_SH_PATH:-${HOME}/miniforge3/etc/profile.d/conda.sh}"
if [[ ! -f "$CONDA_SH_PATH" ]] && command -v conda >/dev/null 2>&1; then
    CONDA_SH_PATH="$(conda info --base)/etc/profile.d/conda.sh"
fi
source "$CONDA_SH_PATH"
conda activate sra3 || true
export PATH="$HOME/edirect:$PATH"

which python3
python3 --version
which esearch
which efetch

echo "[INFO] Running SRR extraction for all datasets"
python3 "$SCRIPT" all

echo "== $(date) :: DONE =="
