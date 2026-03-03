#!/bin/bash
#SBATCH --job-name=download_retinoblastoma_nbl
#SBATCH --chdir=/home/chu25/TCGA
#SBATCH --output=/home/chu25/TCGA/logs/download_retinoblastoma_nbl/%x_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/download_retinoblastoma_nbl/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=168:00:00

set -euo pipefail

# --- environment setup ---
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3

# ====== CONFIG ======
BASE="/home/chu25/data"
TARGET_DIR="${BASE}/nbl"     # TARGET-NBL 
RBL_DIR="${BASE}/rbl"        # Retinoblastoma root
SRA_DIR="${RBL_DIR}/sra"     # SRA cache/output
FASTQ_DIR="${RBL_DIR}/fastq" # FASTQs here
LOG_DIR="${RBL_DIR}/logs"    # logs

# Use all CPUs Slurm gives us
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Temp dir: prefer node-local SSD if Slurm provides it
if [[ -n "${SLURM_TMPDIR:-}" && -d "$SLURM_TMPDIR" ]]; then
  TMPDIR="$SLURM_TMPDIR"
else
  TMPDIR="${RBL_DIR}/tmp"
fi

mkdir -p "$TARGET_DIR" "$RBL_DIR" "$SRA_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMPDIR"

echo "== $(date) =="
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID:-NA}  CPUs: $THREADS  MEM: ${SLURM_MEM_PER_NODE:-64G}"
echo "BASE: $BASE"
echo "TMPDIR: $TMPDIR"

# Point SRA toolkit under BASE (avoid cluttering $HOME)
export NCBI_SETTINGS="${RBL_DIR}/ncbi_settings"
export NCBI_VDB_HTTPS=YES

retry_prefetch() {
  # retry prefetch up to 3x on transient failures
  local listfile="$1"; local outdir="$2"; local attempt=1
  until prefetch --option-file "$listfile" --output-directory "$outdir"; do
    if (( attempt >= 3 )); then
      echo "prefetch failed after $attempt attempts for $listfile" >&2
      return 1
    fi
    echo "prefetch failed (attempt $attempt). Sleeping 60s and retrying…"
    sleep 60
    ((attempt++))
  done
}

gz_all_fastq() {
  # use pigz if available; else gzip
  if command -v pigz >/dev/null 2>&1; then
    pigz -p "${THREADS}" "$FASTQ_DIR"/*.fastq
  else
    parallel --will-cite -j "${THREADS}" gzip ::: "$FASTQ_DIR"/*.fastq || \
    xargs -n1 -P"${THREADS}" gzip <<<"$(printf "%s\n" "$FASTQ_DIR"/*.fastq)"
  fi
}

# ====== PART 1 : TARGET-NBL via TCGAbiolinks ======
#echo "--- Part 1 : TARGET-NBL with TCGAbiolinks ---"
#cat > "${TARGET_DIR}/download_target_nbl.R" <<'EOF'
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
#suppressPackageStartupMessages(library(TCGAbiolinks))
#qry <- GDCquery(project = "TARGET-NBL",
#                data.category = "Transcriptome Profiling",
#                data.type = "Gene Expression Quantification",
#                workflow.type = "STAR - Counts")
#GDCdownload(qry)
#nbl <- GDCprepare(qry)
#saveRDS(nbl, file = "TARGET_NBL_STAR_Counts_data.rds")
#cat("Saved TARGET_NBL_STAR_Counts_data.rds\n")
#EOF

#(
#  cd "$TARGET_DIR"
#  Rscript download_target_nbl.R | tee "${LOG_DIR}/target_nbl.log"
#)

# ====== PART 2: GSE196420 (RB tumors) ======
echo
echo "--- Part 2: GSE196420 (primary RB tumors) ---"

cat > "${RBL_DIR}/rb_GSE196420.txt" <<'EOF'
SRR17960480
SRR17960481
SRR17960482
SRR17960483
SRR17960497
SRR17960498
SRR17960499
SRR17960500
SRR17960501
SRR17960502
SRR17960503
SRR17960504
SRR17960505
SRR17960506
SRR17960507
SRR17960508
SRR17960509
SRR17960510
SRR17960511
SRR17960512
EOF

echo "[GSE196420] prefetch SRA → ${SRA_DIR}"
retry_prefetch "${RBL_DIR}/rb_GSE196420.txt" "${SRA_DIR}" | tee "${LOG_DIR}/prefetch_GSE196420.log"

echo "[GSE196420] fasterq-dump → FASTQ"
while read -r SRR; do
  SRA="${SRA_DIR}/${SRR}.sra"
  if [[ -f "$SRA" ]]; then
    fasterq-dump "$SRA" \
      --split-files \
      --threads "${THREADS}" \
      --temp "${TMPDIR}" \
      --outdir "${FASTQ_DIR}" | tee -a "${LOG_DIR}/fasterq_GSE196420.log"
  else
    echo "WARNING: $SRA not found" | tee -a "${LOG_DIR}/fasterq_GSE196420.log"
  fi
done < "${RBL_DIR}/rb_GSE196420.txt"

echo "[GSE196420] compress FASTQs"
gz_all_fastq

# ====== PART 3: GSE111168 (RB tumors only) ======
echo
echo "--- Part 3: GSE111168 (tumor-only) ---"
RUNINFO="${RBL_DIR}/SRP133555_RunInfo.csv"
wget -q -O "${RUNINFO}" 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP133555'

# Filter: include tumor/RB; exclude obvious normals/retina controls
awk -F',' '
  BEGIN{IGNORECASE=1}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="Run") r=i
      if($i=="LibraryStrategy") ls=i
      if($i=="SampleName") s=i
      if($i=="LibraryName") l=i
    } next
  }
  $ls ~ /RNA-Seq/ {
    sn = ($s=="" ? "" : $s)
    ln = ($l=="" ? "" : $l)
    txt = sn " " ln
    # include tumor/RB patterns, exclude normals
    if (txt ~ /tumou?r|retinoblastoma|(^|[^A-Za-z0-9])rb([^A-Za-z0-9]|$)/ && txt !~ /normal|control|retina/) print $r
  }
' "${RUNINFO}" > "${RBL_DIR}/rb_GSE111168_tumor.txt"

echo "[GSE111168] Selected tumor SRRs:"
cat "${RBL_DIR}/rb_GSE111168_tumor.txt" | tee "${LOG_DIR}/GSE111168_tumor_SRRs.log"

if [[ ! -s "${RBL_DIR}/rb_GSE111168_tumor.txt" ]]; then
  echo "Filter yielded no SRRs — falling back to all SRP133555 runs (you can drop normals later)." | tee -a "${LOG_DIR}/GSE111168_tumor_SRRs.log"
  awk -F',' 'NR>1{print $1}' "${RUNINFO}" > "${RBL_DIR}/rb_GSE111168_tumor.txt"
fi

echo "[GSE111168] prefetch tumor SRRs → ${SRA_DIR}"
retry_prefetch "${RBL_DIR}/rb_GSE111168_tumor.txt" "${SRA_DIR}" | tee "${LOG_DIR}/prefetch_GSE111168.log"

echo "[GSE111168] fasterq-dump → FASTQ"
while read -r SRR; do
  SRA="${SRA_DIR}/${SRR}.sra"
  if [[ -f "$SRA" ]]; then
    fasterq-dump "$SRA" \
      --split-files \
      --threads "${THREADS}" \
      --temp "${TMPDIR}" \
      --outdir "${FASTQ_DIR}" | tee -a "${LOG_DIR}/fasterq_GSE111168.log"
  fi
done < "${RBL_DIR}/rb_GSE111168_tumor.txt"

echo "[GSE111168] compress FASTQs"
gz_all_fastq

echo
echo "=== DONE ==="
echo "TARGET-NBL: ${TARGET_DIR}/TARGET_NBL_STAR_Counts_data.rds"
echo "RB FASTQs: ${FASTQ_DIR} (gzipped)"
echo "SRA objs:  ${SRA_DIR}"
echo "Logs:      ${LOG_DIR}"
