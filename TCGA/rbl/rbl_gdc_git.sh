#!/bin/bash
#SBATCH --job-name=download_rbl_all
#SBATCH --chdir=/home/chu25/data/rbl
#SBATCH --output=/home/chu25/TCGA/rbl/logs/rbl_all_%A_%a.out
#SBATCH --error=/home/chu25/TCGA/rbl/logs/rbl_all_%A_%a.err
#SBATCH --array=1-100%10  # Adjust max based on total SRRs; %10 for concurrency limit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

# ---------- config ----------
BASE=/home/chu25/data/rbl
LOG_DIR="$BASE/logs"
TMP_BASE="$BASE/tmp"
MASTER_LIST="$BASE/all_rbl_srrs.txt"
JOB_LOG_DIR="$LOG_DIR"

# Source GitHub environment setup for dispatch variables.
# actions.sh loads GH_OWNER, GH_REPO, GH_TOKEN securely.
source /home/chu25/actions.sh  

# Optional: Export additional context for enriched event payloads.
export JOB_LOG="$PWD/rbl_all_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"  # Path to job log for reference in dispatches.
export HPC_LABEL="rbl_gse"             

JOB_ID=$SLURM_JOB_ID  # SLURM job ID for payload.

mkdir -p "$LOG_DIR"

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
  local sra_dir="$2"
  if [[ -d "${sra_dir}/${srr}" || -f "${sra_dir}/${srr}.sra" ]]; then
    echo "[prefetch] skip: ${srr} already present"
    return 0
  fi
  echo "[prefetch] ${srr} → ${sra_dir}"
  prefetch --max-size 0 -O "$sra_dir" "$srr"
}

convert_one() {
  local srr="$1"
  local gse_dir="$2"
  local sra_dir="$3"
  local fastq_dir="$4"
  local tmpdir="$5"
  local threads="$6"
  echo "---- $(date) :: START ${srr} ----"

  # skip if gzipped FASTQs already exist
  if compgen -G "${fastq_dir}/${srr}_*.fastq.gz" >/dev/null; then
    echo "[skip] ${srr}: *.fastq.gz found"
    echo "---- $(date) :: DONE  ${srr} ----"
    return 0
  fi

  prefetch_if_missing "$srr" "$sra_dir"

  # select source for fasterq-dump
  local src
  if   [[ -d "${sra_dir}/${srr}" ]]; then src="${sra_dir}/${srr}"
  elif [[ -f "${sra_dir}/${srr}.sra" ]]; then src="${sra_dir}/${srr}.sra"
  else src="${srr}"
  fi

  echo "[fasterq-dump] SRC=${src} threads=${threads} tmp=${tmpdir}"
  fasterq-dump "$src" \
    --split-files \
    --skip-technical \
    --threads "$threads" \
    --temp "$tmpdir" \
    --outdir "$fastq_dir"

  # compress any produced FASTQs for this SRR
  if compgen -G "${fastq_dir}/${srr}_*.fastq" >/dev/null; then
    echo "[pigz] compressing ${srr}_*.fastq"
    pigz -p "$threads" "${fastq_dir}/${srr}"_*.fastq
  else
    echo "[WARN] no *.fastq produced for ${srr}"
  fi

  echo "---- $(date) :: DONE  ${srr} ----"
}

# Build master list of unique SRR with GSE dir if missing
if [ ! -s "$MASTER_LIST" ]; then
  echo "Building master SRR list from GSE folders..."
  > "$MASTER_LIST"  # Clear if exists
  for gse_dir in "$BASE"/GSE*; do
    if [ -d "$gse_dir" ]; then
      gse=$(basename "$gse_dir")
      echo "Processing $gse_dir..."
      txt_file=$(find "$gse_dir" -maxdepth 1 -name "*.txt" -print -quit 2>/dev/null)
      if [ -n "$txt_file" ] && [ -s "$txt_file" ]; then
        echo "Found $txt_file in $gse"
        grep -v '^$' "$txt_file" | while read -r srr; do
          echo "$gse_dir $srr" >> "$MASTER_LIST"
        done
      else
        echo "Warning: No .txt file found in $gse; skipping this directory"
      fi
    fi
  done
  # Remove duplicates (based on SRR only, but keep first GSE dir occurrence) and empty lines
  awk '!seen[$2]++' "$MASTER_LIST" > "${MASTER_LIST}.tmp" && mv "${MASTER_LIST}.tmp" "$MASTER_LIST"
  echo "Master list built with $(wc -l < "$MASTER_LIST") unique SRRs"
fi

if [[ ! -s "$MASTER_LIST" ]]; then
  echo "ERROR: Master list file not found or empty: $MASTER_LIST"
  exit 2
fi

# count non-empty, non-comment lines
TOTAL=$(grep -Evc '^\s*(#|$)' "$MASTER_LIST" || true)
echo "[INFO] LIST=$MASTER_LIST  TOTAL=$TOTAL"

# Pick line for this array task
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL" ]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds list length ($TOTAL)."
  exit 0
fi
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MASTER_LIST")
GSE_DIR=$(echo "$line" | awk '{print $1}')
SRR=$(echo "$line" | awk '{print $2}')

SRA_DIR="$GSE_DIR/sra"
FASTQ_DIR="$GSE_DIR/fastq"
TMPDIR="${SLURM_TMPDIR:-$GSE_DIR/tmp}"
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMPDIR"

THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[ARRAY] ($SLURM_ARRAY_TASK_ID/$TOTAL) $SRR in $GSE_DIR"

convert_one "$SRR" "$GSE_DIR" "$SRA_DIR" "$FASTQ_DIR" "$TMPDIR" "$THREADS"

echo "== $(date) :: ALL DONE =="

# Capture exit status (will be 0 if all succeeded, else non-zero due to set -e)
STATUS=$?

# Dispatch event to GitHub: Triggers workflows like HPC Job Results Sync.
# Uses conditional bash for event_type and status based on $STATUS.
curl -X POST \
  -H "Authorization: token $GH_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/repos/$GH_OWNER/$GH_REPO/dispatches \
  -d "{\"event_type\":\"$([ $STATUS -eq 0 ] && echo 'slurm-job-complete' || echo 'slurm-job-failed')\",\"client_payload\":{\"status\":\"$([ $STATUS -eq 0 ] && echo 'success' || echo 'failure')\",\"job_id\":\"$JOB_ID\",\"scheduler\":\"slurm\",\"host\":\"$(hostname)\",\"user\":\"$USER\",\"note\":\"SRR download complete\",\"label\":\"$HPC_LABEL\"}}"

# Exit with the original status to propagate success/failure in SLURM accounting.
exit $STATUS