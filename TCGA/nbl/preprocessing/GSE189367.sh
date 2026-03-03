#!/bin/bash
#SBATCH --job-name=GSE189367	
#SBATCH --chdir=/home/chu25
#SBATCH --output=/home/chu25/TCGA/logs/GSE189367/%j.out
#SBATCH --error=/home/chu25/TCGA/logs/GSE189367/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

# ---------- config ----------
LIST=/home/chu25/data/nbl/GSE189367/SRP347311.txt  # one SRR per line
BASE=/home/chu25/data/nbl/GSE189367
SRA_DIR="${BASE}/sra"
FASTQ_DIR="${BASE}/fastq"
LOG_ROOT=/home/chu25/TCGA/logs
JOB_LOG_DIR="${LOG_ROOT}/GSE189367"
TMPDIR="${SLURM_TMPDIR:-${BASE}/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$JOB_LOG_DIR" "$TMPDIR"

# ---------- env: conda sra3 ----------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3

# TLS CA for NCBI HTTPS (avoid cert errors)
if [ -f /etc/ssl/certs/ca-certificates.crt ]; then
  export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
elif [ -f /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ]; then
  export SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem
fi

# ---------- sanity checks ----------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }
need prefetch
need fasterq-dump
need pigz

# ---------- helpers ----------
prefetch_if_missing() {
  local srr="$1"
  if [[ -d "${SRA_DIR}/${srr}" || -f "${SRA_DIR}/${srr}.sra" ]]; then
    echo "[prefetch] skip: ${srr} already present"
    return 0
  fi
  echo "[prefetch] ${srr} → ${SRA_DIR}"
  prefetch --max-size 0 -O "$SRA_DIR" "$srr"
}

convert_one() {
  local srr="$1"
  echo "---- $(date) :: START ${srr} ----"

  # skip if gzipped FASTQs already exist
  if compgen -G "${FASTQ_DIR}/${srr}_*.fastq.gz" >/dev/null; then
    echo "[skip] ${srr}: *.fastq.gz found"
    echo "---- $(date) :: DONE  ${srr} ----"
    return 0
  fi

  prefetch_if_missing "$srr"

  # select source for fasterq-dump
  local src
  if   [[ -d "${SRA_DIR}/${srr}" ]]; then src="${SRA_DIR}/${srr}"
  elif [[ -f "${SRA_DIR}/${srr}.sra" ]]; then src="${SRA_DIR}/${srr}.sra"
  else src="${srr}"
  fi

  echo "[fasterq-dump] SRC=${src} threads=${THREADS} tmp=${TMPDIR}"
  fasterq-dump "$src" \
    --split-files \
    --skip-technical \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    --outdir "$FASTQ_DIR"

  # compress any produced FASTQs for this SRR
  if compgen -G "${FASTQ_DIR}/${srr}_*.fastq" >/dev/null; then
    echo "[pigz] compressing ${srr}_*.fastq"
    pigz -p "$THREADS" "${FASTQ_DIR}/${srr}"_*.fastq
  else
    echo "[WARN] no *.fastq produced for ${srr}"
  fi

  echo "---- $(date) :: DONE  ${srr} ----"
}

# ---------- run sequentially ----------
if [[ ! -s "$LIST" ]]; then
  echo "ERROR: List file not found or empty: $LIST"
  exit 2
fi

# count non-empty, non-comment lines
TOTAL=$(grep -Evc '^\s*(#|$)' "$LIST" || true)
echo "[INFO] LIST=$LIST  TOTAL=$TOTAL  OUT=$FASTQ_DIR  SRA=$SRA_DIR  TMP=$TMPDIR"

i=0
while IFS= read -r line; do
  # strip whitespace and skip blanks/comments
  srr="${line//[[:space:]]/}"
  [[ -z "$srr" || "${srr:0:1}" == "#" ]] && continue
  i=$((i+1))
  echo "[SEQ] ($i/$TOTAL) $srr"
  convert_one "$srr"
done < "$LIST"

echo "== $(date) :: ALL DONE =="
