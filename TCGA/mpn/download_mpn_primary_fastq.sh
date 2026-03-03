#!/bin/bash
#SBATCH --job-name=download_mpn_primary_fastq
#SBATCH --output=/home/chu25/TCGA/mpn/logs/download_mpn_primary_fastq_%j.out
#SBATCH --error=/home/chu25/TCGA/mpn/logs/download_mpn_primary_fastq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=120:00:00

set -euo pipefail

echo "=== [$HOSTNAME] Starting MPN monocyte bulk RNA-seq download at $(date) ==="

LOG_DIR="/home/chu25/TCGA/mpn/logs"
OUT_DIR="/home/chu25/data/mpn/GSE228758/fastq"

mkdir -p "$LOG_DIR" "$OUT_DIR"

# Activate conda env with SRA Toolkit 3.x
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3 || { echo "[ERROR] Failed to activate conda env 'sra3'"; exit 1; }

echo "[INFO] Using conda env: $CONDA_DEFAULT_ENV"
echo "[INFO] Output directory: $OUT_DIR"

echo "[INFO] which prefetch: $(which prefetch || echo 'prefetch not found')"
echo "[INFO] which fasterq-dump: $(which fasterq-dump || echo 'fasterq-dump not found')"

prefetch --version || echo "[WARN] Could not get prefetch version"
fasterq-dump --version || echo "[WARN] Could not get fasterq-dump version"

cd "$OUT_DIR" || { echo "[ERROR] Could not cd into $OUT_DIR"; exit 1; }

# Use pigz if available for faster compression
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -p 8"
  echo "[INFO] Using pigz for compression"
else
  COMPRESS="gzip"
  echo "[INFO] pigz not found, using gzip"
fi

# MPN monocyte samples from PRJNA951430 / SRP430505
SRRS=(
  SRR24043724
  SRR24043725
  SRR24043726
  SRR24043727
  SRR24043728
  SRR24043729
)

for SRRID in "${SRRS[@]}"; do
  echo "============================================================"
  echo "[INFO] Downloading and converting $SRRID ..."
  echo "============================================================"

  # Download .sra via prefetch
  prefetch "$SRRID"
  if [ $? -ne 0 ]; then
    echo "[ERROR] prefetch failed for $SRRID, skipping."
    continue
  fi
  echo "[INFO] prefetch complete for $SRRID"

  # Convert to FASTQ (paired-end, uncompressed)
  fasterq-dump "$SRRID" --split-files --threads 8 -O "$OUT_DIR"
  if [ $? -ne 0 ]; then
    echo "[ERROR] fasterq-dump failed for $SRRID"
    continue
  fi
  echo "[INFO] FASTQ conversion complete for $SRRID"

  # Compress the FASTQ files
  if [ -f "${SRRID}_1.fastq" ]; then
    echo "[INFO] Compressing ${SRRID}_1.fastq and ${SRRID}_2.fastq ..."
    $COMPRESS "${SRRID}_1.fastq" "${SRRID}_2.fastq"
  else
    echo "[WARN] Expected FASTQs for ${SRRID} not found; skipping compression."
  fi

  # Optional: remove cached .sra to save space
  # CACHE="$HOME/ncbi/public/sra/${SRRID}.sra"
  # if [ -f "$CACHE" ]; then
  #   echo "[INFO] Removing cached SRA file $CACHE"
  #   rm -f "$CACHE"
  # fi

done

echo "=== Finished MPN monocyte bulk RNA-seq downloads at $(date) ==="
exit 0
