#!/bin/bash
#SBATCH --job-name=run_download_nbl_srr_ids
#SBATCH --chdir=/work/ugbogu/pipeline/data/nbl
#SBATCH --output=/work/ugbogu/pipeline/data/nbl/logs/run_download_nbl_srr_ids_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/nbl/logs/run_download_nbl_srr_ids_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

ROOT=/work/ugbogu/pipeline/data/nbl
LOG_DIR="${ROOT}/logs"
SCRIPT=/work/ugbogu/pipeline/preprocessing_and_quality_control/scripts/download_nbl_tumour_sample_srr_ids.py

mkdir -p "$LOG_DIR"

echo "[INFO] ROOT=$ROOT"
echo "[INFO] LOG_DIR=$LOG_DIR"
echo "[INFO] SCRIPT=$SCRIPT"

source /home/ugbogu/miniforge3/etc/profile.d/conda.sh
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