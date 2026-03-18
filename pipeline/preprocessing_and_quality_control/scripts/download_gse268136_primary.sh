#!/bin/bash
#SBATCH --job-name=GSE268136_primary
#SBATCH --chdir=/work/ugbogu/pipeline/data/rbl
#SBATCH --output=/work/ugbogu/pipeline/data/rbl/GSE268136_primary_tumours/logs/GSE268136_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/rbl/GSE268136_primary_tumours/logs/GSE268136_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

# ---------- config ----------
BASE=/work/ugbogu/pipeline/data/rbl/GSE268136_primary_tumours
LIST="${BASE}/all_primary_tumour_pe.txt"
SRA_DIR="${BASE}/sra"
FASTQ_DIR="${BASE}/fastq"
LOG_ROOT="${BASE}/logs"
JOB_LOG_DIR="${LOG_ROOT}"
TMPDIR="${SLURM_TMPDIR:-${BASE}/tmp/${SLURM_JOB_ID:-manual}}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$JOB_LOG_DIR" "$TMPDIR"

# ---------- env ----------
source /home/ugbogu/miniforge3/etc/profile.d/conda.sh
conda activate sra3

# ---------- TLS CA for NCBI HTTPS ----------
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
  echo "[prefetch] ${srr} -> ${SRA_DIR}"
  prefetch --max-size 0 -O "$SRA_DIR" "$srr"
}

convert_one() {
  local srr="$1"
  echo "---- $(date) :: START ${srr} ----"

  if [[ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" && -s "${FASTQ_DIR}/${srr}_2.fastq.gz" ]]; then
    echo "[skip] ${srr}: paired FASTQs already exist"
    echo "---- $(date) :: DONE ${srr} ----"
    return 0
  fi

  rm -f \
    "${FASTQ_DIR}/${srr}.fastq" \
    "${FASTQ_DIR}/${srr}_1.fastq" \
    "${FASTQ_DIR}/${srr}_2.fastq" \
    "${FASTQ_DIR}/${srr}.fastq.gz" \
    "${FASTQ_DIR}/${srr}_1.fastq.gz" \
    "${FASTQ_DIR}/${srr}_2.fastq.gz"

  prefetch_if_missing "$srr"

  local src
  if [[ -d "${SRA_DIR}/${srr}" ]]; then
    src="${SRA_DIR}/${srr}"
  elif [[ -f "${SRA_DIR}/${srr}.sra" ]]; then
    src="${SRA_DIR}/${srr}.sra"
  else
    src="${srr}"
  fi

  echo "[fasterq-dump] SRC=${src} threads=${THREADS} tmp=${TMPDIR}"
  fasterq-dump "$src" \
    --split-files \
    --skip-technical \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    --outdir "$FASTQ_DIR"

  if [[ -s "${FASTQ_DIR}/${srr}_1.fastq" && -s "${FASTQ_DIR}/${srr}_2.fastq" ]]; then
    echo "[pigz] compressing ${srr}_1.fastq and ${srr}_2.fastq"
    pigz -p "$THREADS" "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${srr}_2.fastq"
  else
    echo "[ERROR] expected paired FASTQs not produced for ${srr}" >&2
    return 1
  fi

  echo "---- $(date) :: DONE ${srr} ----"
}

# ---------- run sequentially ----------
if [[ ! -s "$LIST" ]]; then
  echo "ERROR: List file not found or empty: $LIST"
  exit 2
fi

TOTAL=$(grep -Evc '^\s*(#|$)' "$LIST")
echo "[INFO] LIST=$LIST TOTAL=$TOTAL OUT=$FASTQ_DIR SRA=$SRA_DIR TMP=$TMPDIR"

i=0
while IFS= read -r line; do
  srr="${line//[[:space:]]/}"
  [[ -z "$srr" || "${srr:0:1}" == "#" ]] && continue
  i=$((i+1))
  echo "[SEQ] ($i/$TOTAL) $srr"
  convert_one "$srr"
done < "$LIST"

echo "== $(date) :: ALL DONE =="