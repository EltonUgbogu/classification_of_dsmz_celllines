#!/bin/bash
#SBATCH --job-name=download_retinoblastoma_nbln
#SBATCH --chdir=/home/chu25/TCGA
#SBATCH --output=/home/chu25/TCGA/logs/download_retinoblastoma_nbl/%x_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/download_retinoblastoma_nbl/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --array=1-3%3
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=168:00:00

set -euo pipefail

# --- env ---
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3

# TLS bundle (pick what exists)
if [ -f /etc/ssl/certs/ca-certificates.crt ]; then
  export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
elif [ -f /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ]; then
  export SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem
fi

BASE=/home/chu25/data/rbl
SRA_DIR="$BASE/sra"
FASTQ_DIR="$BASE/fastq"
LOG_DIR="$BASE/logs"
TMPDIR="${SLURM_TMPDIR:-$BASE/tmp}"
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMPDIR"

# The three big SRRs:
SRRS=( SRR17960480 SRR17960482 SRR17960483 )
SRR="${SRRS[$((SLURM_ARRAY_TASK_ID-1))]}"

THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[$(date)] Start ${SRR} on $(hostname); threads=${THREADS}; tmp=${TMPDIR}"

# Prefetch (no size cap)
if [[ ! -d "${SRA_DIR}/${SRR}" && ! -f "${SRA_DIR}/${SRR}.sra" ]]; then
  prefetch --max-size 0 -O "$SRA_DIR" "$SRR"
fi

# Source for fasterq-dump
if   [[ -d "${SRA_DIR}/${SRR}" ]]; then SRC="${SRA_DIR}/${SRR}"
elif [[ -f "${SRA_DIR}/${SRR}.sra" ]]; then SRC="${SRA_DIR}/${SRR}.sra"
else SRC="${SRR}"
fi

# Convert; skip technical to avoid giant *_3.fastq
fasterq-dump "$SRC" \
  --split-files \
  --skip-technical \
  --threads "$THREADS" \
  --temp "$TMPDIR" \
  --outdir "$FASTQ_DIR"

# Compress only this SRR’s fastqs
if compgen -G "${FASTQ_DIR}/${SRR}_*.fastq" > /dev/null; then
  pigz -p "$THREADS" "${FASTQ_DIR}/${SRR}"_*.fastq
fi

echo "[$(date)] Done ${SRR