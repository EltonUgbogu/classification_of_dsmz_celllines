#!/bin/bash
# =============================================================================
# download_nbl_primary_pe.sh
#
# Downloads, extracts, and compresses paired-end SRA data for neuroblastoma
# primary tumour cohorts. For each SRR accession in a plain-text list, the
# script runs:
#   1. prefetch      — downloads the SRA archive from NCBI.
#   2. fasterq-dump  — extracts paired FASTQ files (_1.fastq, _2.fastq).
#   3. pigz          — compresses both files in parallel.
#
# Supported datasets:
#   GSE100148
#   GSE189367
#   SRP409177
#
# Expected per-dataset input:
#   /work/ugbogu/pipeline/data/nbl/<DATASET>/all_primary_tumour_pe.txt
#
# Output directories created automatically:
#   /work/ugbogu/pipeline/data/nbl/<DATASET>/sra
#   /work/ugbogu/pipeline/data/nbl/<DATASET>/fastq
#   /work/ugbogu/pipeline/data/nbl/<DATASET>/logs
#   /work/ugbogu/pipeline/data/nbl/<DATASET>/tmp/<jobid>
#
# The script is idempotent: completed samples are skipped on resubmission.
#
# Usage:
#   sbatch download_nbl_primary_pe.sh [DATASET]
#
# Example:
#   sbatch download_nbl_primary_pe.sh GSE100148
#   sbatch download_nbl_primary_pe.sh GSE189367
#   sbatch download_nbl_primary_pe.sh SRP409177
# =============================================================================

# SLURM directives — command-line sbatch flags override these defaults.
#SBATCH --job-name=download_nbl_primary_pe
#SBATCH --chdir=/work/ugbogu/pipeline/data/nbl
#SBATCH --output=/work/ugbogu/pipeline/data/nbl/%x_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/nbl/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail   # Exit on error, treat unset variables as errors, propagate pipe failures.


# =============================================================================
# Section 1 — Configuration
# =============================================================================
echo "== $(date) :: START =="
echo "Host: $(hostname)"

ROOT=/work/ugbogu/pipeline/data/nbl
DATASET="${1:-GSE100148}"
BASE="${ROOT}/${DATASET}"
LIST="${BASE}/all_primary_tumour_pe.txt"    # One SRR accession per line.
SRA_DIR="${BASE}/sra"                       # prefetch output directory.
FASTQ_DIR="${BASE}/fastq"                   # fasterq-dump and pigz output directory.
LOG_DIR="${BASE}/logs"
TMPDIR="${SLURM_TMPDIR:-${BASE}/tmp/${SLURM_JOB_ID:-manual}}"  # Job-specific scratch space.
THREADS="${SLURM_CPUS_PER_TASK:-8}"         # Passed to fasterq-dump and pigz.

case "$DATASET" in
    GSE100148|GSE189367|SRP409177)
        ;;
    *)
        echo "ERROR: Unsupported dataset '$DATASET'" >&2
        echo "Supported datasets: GSE100148, GSE189367, SRP409177" >&2
        exit 2
        ;;
esac

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMPDIR"

echo "[INFO] ROOT=$ROOT"
echo "[INFO] DATASET=$DATASET"
echo "[INFO] BASE=$BASE"
echo "[INFO] LIST=$LIST"
echo "[INFO] SRA_DIR=$SRA_DIR"
echo "[INFO] FASTQ_DIR=$FASTQ_DIR"
echo "[INFO] LOG_DIR=$LOG_DIR"
echo "[INFO] TMPDIR=$TMPDIR"
echo "[INFO] THREADS=$THREADS"


# =============================================================================
# Section 2 — Conda Environment Activation
#
# The sra3 environment provides prefetch, fasterq-dump, and pigz. The profile
# script must be sourced explicitly because SLURM shells are non-interactive.
# =============================================================================
source /home/ugbogu/miniforge3/etc/profile.d/conda.sh
conda activate sra3


# =============================================================================
# Section 3 — TLS Certificate Configuration
#
# Sets SSL_CERT_FILE to the system CA bundle so prefetch can validate NCBI
# HTTPS endpoints. Covers both Debian (/etc/ssl/) and RHEL (/etc/pki/) systems.
# =============================================================================
if [[ -f /etc/ssl/certs/ca-certificates.crt ]]; then
    export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
elif [[ -f /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ]]; then
    export SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem
fi


# =============================================================================
# Section 4 — Dependency Validation
#
# Checks that required tools are on PATH before processing begins. Exit 127
# is the POSIX convention for 'command not found'.
# =============================================================================
need() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 127; }
}

need prefetch
need fasterq-dump
need pigz

if [[ ! -d "$BASE" ]]; then
    echo "ERROR: Dataset directory does not exist: $BASE"
    exit 2
fi

if [[ ! -s "$LIST" ]]; then
    echo "ERROR: List file not found or empty: $LIST"
    exit 2
fi


# =============================================================================
# Section 5 — Helper Functions
# =============================================================================

# prefetch_if_missing()
# Downloads the SRA archive for $srr unless it already exists locally.
# Recognises both the directory layout (SRA Toolkit ≥ 2.10) and the flat
# .sra file layout (older versions). --max-size 0 disables the size cap.
prefetch_if_missing() {
    local srr="$1"

    if [[ -d "${SRA_DIR}/${srr}" || -f "${SRA_DIR}/${srr}.sra" ]]; then
        echo "[prefetch] skip: ${srr} already present"
        return 0
    fi

    echo "[prefetch] downloading ${srr} -> ${SRA_DIR}"
    prefetch --max-size 0 -O "$SRA_DIR" "$srr"
}

# resolve_sra_source()
# Returns the correct source argument for fasterq-dump based on which archive
# layout is present: directory form, flat .sra file, or bare accession.
resolve_sra_source() {
    local srr="$1"

    if [[ -d "${SRA_DIR}/${srr}" && -f "${SRA_DIR}/${srr}/${srr}.sra" ]]; then
        echo "${SRA_DIR}/${srr}"
    elif [[ -f "${SRA_DIR}/${srr}.sra" ]]; then
        echo "${SRA_DIR}/${srr}.sra"
    else
        echo "${srr}"
    fi
}

# convert_one()
# Runs the full prefetch → fasterq-dump → pigz pipeline for a single accession.
#
# Guard 1: skips the accession if both compressed paired FASTQs already exist.
# Guard 2: removes incomplete intermediates before re-running.
# Guard 3: verifies paired FASTQs were produced before compression.
convert_one() {
    local srr="$1"
    echo "---- $(date) :: START ${srr} ----"

    # Guard 1: skip if compressed paired FASTQs already exist.
    if [[ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" && -s "${FASTQ_DIR}/${srr}_2.fastq.gz" ]]; then
        echo "[skip] ${srr}: paired FASTQs already exist"
        echo "---- $(date) :: DONE ${srr} ----"
        return 0
    fi

    # Guard 2: remove incomplete files from any prior failed attempt.
    rm -f \
        "${FASTQ_DIR}/${srr}.fastq" \
        "${FASTQ_DIR}/${srr}_1.fastq" \
        "${FASTQ_DIR}/${srr}_2.fastq" \
        "${FASTQ_DIR}/${srr}.fastq.gz" \
        "${FASTQ_DIR}/${srr}_1.fastq.gz" \
        "${FASTQ_DIR}/${srr}_2.fastq.gz"

    prefetch_if_missing "$srr"

    local src
    src="$(resolve_sra_source "$srr")"

    echo "[fasterq-dump] SRC=${src}"
    echo "[fasterq-dump] threads=${THREADS}"
    echo "[fasterq-dump] tmp=${TMPDIR}"
    echo "[fasterq-dump] outdir=${FASTQ_DIR}"

    fasterq-dump "$src" \
        --split-files \
        --skip-technical \
        --threads "$THREADS" \
        --temp "$TMPDIR" \
        --outdir "$FASTQ_DIR"

    # Guard 3: verify paired output before compression.
    if [[ -s "${FASTQ_DIR}/${srr}_1.fastq" && -s "${FASTQ_DIR}/${srr}_2.fastq" ]]; then
        echo "[pigz] compressing ${srr}_1.fastq and ${srr}_2.fastq"
        pigz -p "$THREADS" "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${srr}_2.fastq"
    else
        echo "[ERROR] expected paired FASTQs not produced for ${srr}" >&2
        return 1
    fi

    echo "---- $(date) :: DONE ${srr} ----"
}


# =============================================================================
# Section 6 — Sequential Processing Loop
#
# Accessions are processed one at a time to avoid saturating network bandwidth.
# Blank lines and comment lines (#) in LIST are skipped.
# =============================================================================
TOTAL=$(grep -Evc '^[[:space:]]*(#|$)' "$LIST")
echo "[INFO] TOTAL SRRs TO PROCESS: $TOTAL"

i=0
while IFS= read -r line; do
    srr="$(echo "$line" | tr -d '[:space:]')"
    [[ -z "$srr" || "${srr:0:1}" == "#" ]] && continue

    i=$((i + 1))
    echo "[SEQ] ($i/$TOTAL) $srr"
    convert_one "$srr"
done < "$LIST"

echo "== $(date) :: ALL DONE =="