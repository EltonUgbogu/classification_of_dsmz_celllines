#!/bin/bash
#SBATCH --job-name=download_retinoblastoma_nbln2
#SBATCH --chdir=/home/chu25/TCGA
#SBATCH --output=/home/chu25/TCGA/logs/download_retinoblastoma_nbln2/%x_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/download_retinoblastoma_nbln2/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --array=1-3%3
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=168:00:00

set -euo pipefail

# --- env ---
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3

# TLS bundle
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
RUNINFO="$BASE/SRP133555_RunInfo.csv"
LIST="$BASE/rb_GSE111168_tumor.txt"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMPDIR"

# Build tumor list if missing
if [ ! -s "$LIST" ]; then
  wget -q -O "$RUNINFO" \
    'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP133555'
  awk -F',' '
    BEGIN{IGNORECASE=1}
    NR==1{for(i=1;i<=NF;i++){if($i=="Run") r=i; if($i=="LibraryStrategy") ls=i; if($i=="SampleName") s=i; if($i=="LibraryName") l=i;} next}
    $ls ~ /RNA-Seq/ { txt=$s" "$l; if (txt ~ /tumou?r|retinoblastoma|(^|[^a-z0-9])rb([^a-z0-9]|$)/ && txt !~ /normal|control|retina/) print $r }
  ' "$RUNINFO" > "$LIST"
  # fallback: keep all if filter empty
  if [ ! -s "$LIST" ]; then awk -F',' 'NR>1{print $1}' "$RUNINFO" > "$LIST"; fi
fi

# Pick SRR for this array task
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
if [[ -z "${SRR:-}" ]]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds list length ($(wc -l < "$LIST"))."
  exit 0
fi

THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[$(date)] Start ${SRR} on $(hostname); threads=${THREADS}; tmp=${TMPDIR}"

# Prefetch (no size cap)
if [[ ! -d "${SRA_DIR}/${SRR}" && ! -f "${SRA_DIR}/${SRR}.sra" ]]; then
  prefetch --max-size 0 -O "$SRA_DIR" "$SRR"
fi

# Source
if   [[ -d "${SRA_DIR}/${SRR}" ]]; then SRC="${SRA_DIR}/${SRR}"
elif [[ -f "${SRA_DIR}/${SRR}.sra" ]]; then SRC="${SRA_DIR}/${SRR}.sra"
else SRC="${SRR}"
fi

# Convert
fasterq-dump "$SRC" \
  --split-files \
  --skip-technical \
  --threads "$THREADS" \
  --temp "$TMPDIR" \
  --outdir "$FASTQ_DIR"

# Compress
if compgen -G "${FASTQ_DIR}/${SRR}_*.fastq" > /dev/null; then
  pigz -p "$THREADS" "${FASTQ_DIR}/${SRR}"_*.fastq
fi

echo "[$(date)] Done ${SRR}"