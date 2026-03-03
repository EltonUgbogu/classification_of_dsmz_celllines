#!/bin/bash
#SBATCH --job-name=download_mmml_rnaseq
#SBATCH --output=/home/chu25/TCGA/icgc_mmml_bl_dlbcl_fl/logs/download_mmml_rnaseq_%j.out
#SBATCH --error=/home/chu25/TCGA/icgc_mmml_bl_dlbcl_fl/logs/download_mmml_rnaseq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00

set -euo pipefail

echo "=== [$HOSTNAME] Starting EGAD00001001672 download at $(date) ==="

OUT_DIR="/home/chu25/data/icgc_mmml/EGAD00001001672/fastq"
CRED="/home/chu25/.ega/credentials.json"

mkdir -p "$OUT_DIR"
mkdir -p /home/chu25/TCGA/icgc_mmml_bl_dlbcl_fl/logs

# Activate env that has pyega3
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv || { echo "[ERROR] Failed to activate conda env"; exit 1; }

echo "[INFO] Using conda env: $CONDA_DEFAULT_ENV"
echo "[INFO] Output dir: $OUT_DIR"

# Optional: check pyega3 version & help
echo "[INFO] pyega3 in use: $(which pyega3)"
pyega3 -h | head -n 5 || echo "[WARN] Could not show pyega3 help"

# Move into output directory because older pyega3 writes to CWD
cd "$OUT_DIR" || { echo "[ERROR] Cannot cd to $OUT_DIR"; exit 1; }

# Optional: list files in dataset first
# pyega3 -cf "$CRED" files EGAD00001001672

# Download all files in the dataset (into OUT_DIR/CWD)
pyega3 -cf "$CRED" fetch EGAD00001001672

echo "=== Finished EGAD00001001672 download at $(date) ==="
