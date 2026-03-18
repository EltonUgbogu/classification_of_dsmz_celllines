#!/bin/bash
#SBATCH --job-name=extract_srr_ids
#SBATCH --output=/work/ugbogu/pipeline/data/rbl/logs/extract_srr_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/rbl/logs/extract_srr_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

SCRIPT=/work/ugbogu/pipeline/preprocessing_and_quality_control/scripts/download_rbl_tumour_sample_srr_ids.py

source /home/ugbogu/miniforge3/etc/profile.d/conda.sh
conda activate sra3 || true
export PATH="$HOME/edirect:$PATH"

which python3
python3 --version
which esearch
which efetch

echo "[INFO] Running SRR extraction for all datasets"
python3 "$SCRIPT" all

echo "== $(date) :: DONE =="