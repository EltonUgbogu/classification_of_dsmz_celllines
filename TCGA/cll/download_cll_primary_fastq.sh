#!/bin/bash
#SBATCH --job-name=download_cll_primary_fastq
#SBATCH --output=/home/chu25/TCGA/cll/logs/download_cll_primary_fastq_%j.out
#SBATCH --error=/home/chu25/TCGA/cll/logs/download_cll_primary_fastq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=120:00:00

set -euo pipefail

echo "=== [$HOSTNAME] Starting CLL bulk RNA-seq download at $(date) ==="

LOG_DIR="/home/chu25/TCGA/cll/logs"
OUT_DIR="/home/chu25/data/cll/GSE237907/fastq"

mkdir -p "$LOG_DIR" "$OUT_DIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3 || { echo "[ERROR] Failed to activate conda env 'sra3'"; exit 1; }

echo "[INFO] Using conda env: $CONDA_DEFAULT_ENV"
echo "[INFO] Output directory: $OUT_DIR"
echo "[INFO] which prefetch: $(which prefetch || echo 'not found')"
echo "[INFO] which fasterq-dump: $(which fasterq-dump || echo 'not found')"

prefetch --version || true
fasterq-dump --version || true

SRRS=(
  SRR25381845
  SRR25381846
  SRR25381847
  SRR25381848
  SRR25381849
  SRR25381850
  SRR25381851
  SRR25381852
  SRR25381853
  SRR25381854
)

cd "$OUT_DIR" || { echo "[ERROR] Could not cd into $OUT_DIR"; exit 1; }

if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -p 8"
  echo "[INFO] Using pigz for compression"
else
  COMPRESS="gzip"
  echo "[INFO] Using gzip for compression"
fi

for SRRID in "${SRRS[@]}"; do
  echo "============================================================"
  echo "[INFO] Downloading and converting $SRRID ..."
  echo "============================================================"

  prefetch "$SRRID" || { echo "[ERROR] prefetch failed for $SRRID, skipping"; continue; }
  echo "[INFO] prefetch complete for $SRRID"

  fasterq-dump "$SRRID" --split-files --threads 8 -O "$OUT_DIR" || {
    echo "[ERROR] fasterq-dump failed for $SRRID"; continue;
  }
  echo "[INFO] FASTQ conversion complete for $SRRID"

  if [ -f "${SRRID}_1.fastq" ]; then
    echo "[INFO] Compressing ${SRRID}_1.fastq and ${SRRID}_2.fastq ..."
    $COMPRESS "${SRRID}_1.fastq" "${SRRID}_2.fastq"
  else
    echo "[WARN] ${SRRID}_1.fastq not found; skipping compression."
  fi
done

echo "=== Finished CLL bulk RNA-seq downloads at $(date) ==="
exit 0
