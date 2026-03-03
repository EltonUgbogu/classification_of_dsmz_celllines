#!/bin/bash
#SBATCH --job-name=SRP409177_dl
#SBATCH --chdir=/home/chu25
#SBATCH --output=/home/chu25/TCGA/logs/SRP409177_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/SRP409177_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

# ---------- config ----------
LIST=/home/chu25/data/nbl/SRP409177_SRRs.txt   # one SRR per line
BASE=/home/chu25/data/nbl/SRP409177
SRA_DIR=${BASE}/sra
FASTQ_DIR=${BASE}/fastq
LOG_DIR=/home/chu25/TCGA/logs
TMPDIR="${SLURM_TMPDIR:-${BASE}/tmp}"
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMPDIR"

# ---------- env: conda sra3 ----------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3

# TLS CA for NCBI HTTPS (avoid cert errors)
if [ -f /etc/ssl/certs/ca-certificates.crt ]; then
  export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
elif [ -f /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ]; then
  export SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem
fi

# ---------- helpers ----------
prefetch_if_missing() {
  local SRR="$1"
  if [[ -d "${SRA_DIR}/${SRR}" || -f "${SRA_DIR}/${SRR}.sra" ]]; then
    echo "[prefetch] skip: ${SRR} already present"
    return 0
  fi
  echo "[prefetch] ${SRR} → ${SRA_DIR}"
  prefetch --max-size 0 -O "$SRA_DIR" "$SRR"
}

convert_one() {
  local SRR="$1"
  echo "---- $(date) :: START ${SRR} ----"

  # skip if gzipped FASTQs already exist
  if compgen -G "${FASTQ_DIR}/${SRR}_*.fastq.gz" > /dev/null; then
    echo "[skip] ${SRR}: *.fastq.gz found"
    echo "---- $(date) :: DONE  ${SRR} ----"
    return 0
  fi

  prefetch_if_missing "$SRR"

  # choose source for fasterq-dump
  local SRC
  if   [[ -d "${SRA_DIR}/${SRR}" ]]; then SRC="${SRA_DIR}/${SRR}"
  elif [[ -f "${SRA_DIR}/${SRR}.sra" ]]; then SRC="${SRA_DIR}/${SRR}.sra"
  else SRC="${SRR}"
  fi

  echo "[fasterq-dump] SRC=${SRC} threads=${THREADS} tmp=${TMPDIR}"
  fasterq-dump "$SRC" \
    --split-files \
    --skip-technical \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    --outdir "$FASTQ_DIR"

  if compgen -G "${FASTQ_DIR}/${SRR}_*.fastq" > /dev/null; then
    echo "[pigz] compressing ${SRR}_*.fastq"
    pigz -p "$THREADS" "${FASTQ_DIR}/${SRR}"_*.fastq
  else
    echo "[WARN] no *.fastq produced for ${SRR}"
  fi

  echo "---- $(date) :: DONE  ${SRR} ----"
}

# ---------- run sequentially ----------
TOTAL=$(grep -cve '^\s*$' "$LIST" || true)
echo "[INFO] LIST=$LIST  TOTAL=$TOTAL"
i=0
while IFS= read -r SRR; do
  SRR="${SRR//[[:space:]]/}"
  [[ -z "$SRR" ]] && continue
  i=$((i+1))
  echo "[SEQ] ($i/$TOTAL) $SRR"
  convert_one "$SRR"
done < "$LIST"

echo "== $(date) :: ALL DONE =="
