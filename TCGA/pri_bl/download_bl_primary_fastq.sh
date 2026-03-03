#!/bin/bash
#SBATCH --job-name=download_bl_primary_fastq
#SBATCH --output=/home/chu25/TCGA/pri_bl/logs/download_bl_primary_fastq_%j.out
#SBATCH --error=/home/chu25/TCGA/pri_bl/logs/download_bl_primary_fastq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=120:00:00

set -euo pipefail

echo "=== [$HOSTNAME] Starting download BL primary FASTQ job at $(date) ==="

# Make sure log and output dirs exist
mkdir -p /home/chu25/TCGA/pri_bl/logs

OUT_DIR="/home/chu25/data/pri_bl/GSE199413/fastq"
mkdir -p "$OUT_DIR"

# Activate conda env with NEW SRA toolkit
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3 || { echo "[ERROR] Failed to activate conda env 'sra3'"; exit 1; }

echo "[INFO] Using conda env: $CONDA_DEFAULT_ENV"
echo "[INFO] Output directory: $OUT_DIR"

# Show which binaries we are really using
echo "[INFO] which prefetch: $(which prefetch || echo 'prefetch not found')"
echo "[INFO] which fasterq-dump: $(which fasterq-dump || echo 'fasterq-dump not found')"

# Optional: print SRA Toolkit version
prefetch --version || echo "[WARN] Could not get prefetch version"
fasterq-dump --version || echo "[WARN] Could not get fasterq-dump version"

# List of primary BL runs (patients AA1912 & AA2000)
SRRS=(
  SRR18482594
  SRR18482595
  SRR18482596
  SRR18482597
  SRR18482598
  SRR18482599
  SRR18482600
  SRR18482601
)

cd "$OUT_DIR" || { echo "[ERROR] Could not cd into $OUT_DIR"; exit 1; }

# Use pigz if available, otherwise fallback to gzip
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -p 8"
  echo "[INFO] Using pigz for compression"
else
  COMPRESS="gzip"
  echo "[INFO] pigz not found, using gzip"
fi

for SRRID in "${SRRS[@]}"; do
  echo "============================================================"
  echo "[INFO] Downloading and converting $SRRID ..."
  echo "============================================================"

  # Download .sra file (SRA Toolkit: prefetch)
  prefetch "$SRRID"
  if [ $? -ne 0 ]; then
    echo "[ERROR] prefetch failed for $SRRID, skipping."
    continue
  fi
  echo "[INFO] prefetch complete for $SRRID"

  # Convert to FASTQ (paired-end, split files, UNCOMPRESSED)
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

  # Optional: remove cached .sra file to save space
  # DEFAULT CACHE: "$HOME/ncbi/public/sra"
  # rm -f "$HOME/ncbi/public/sra/${SRRID}.sra" 2>/dev/null || true

done

echo "=== Finished BL primary FASTQ downloads at $(date) ==="
exit 0
