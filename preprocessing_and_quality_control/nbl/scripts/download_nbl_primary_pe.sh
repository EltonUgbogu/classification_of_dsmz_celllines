#!/bin/bash
# =============================================================================
# download_nbl_primary_pe.sh
#
# Downloads, extracts, and compresses paired-end SRA data for neuroblastoma
# primary tumour cohorts. For each SRR accession in a plain-text list, the
# script executes a three-stage pipeline:
#   1. prefetch      — downloads the SRA archive from NCBI.
#   2. fasterq-dump  — extracts paired FASTQ files (_1.fastq, _2.fastq).
#   3. pigz          — compresses both files in parallel.
#
# The script is idempotent: accessions whose compressed paired FASTQs already
# exist are skipped on resubmission without re-running any pipeline stage.
#
# Failed accessions are recorded with a structured failure reason in:
#   <DATASET>/removed_nbl_srr_id_primary_tumour_pe.txt
#
# Four distinct failure reasons are distinguished:
#   prefetch_failed_unresolvable  — SRA archive could not be downloaded.
#   fasterq_dump_failed           — FASTQ extraction exited with a non-zero code.
#   paired_fastq_not_produced     — fasterq-dump succeeded but _1/_2 files are absent or empty.
#   pigz_failed                   — compression exited with a non-zero code or produced no output.
#
# Supported datasets:
#   GSE100148
#   GSE189367
#   SRP409177
#
# Usage:
#   sbatch download_nbl_primary_pe.sh [DATASET|all]
#   DATASET defaults to GSE100148 if not supplied.
#   Use "all" to process every supported dataset sequentially within one Slurm job.
#
# Slurm #SBATCH --output/--error are set to /dev/null because paths cannot depend
# on DATASET; after startup the script redirects into each dataset's logs/ dir.
# =============================================================================

#SBATCH --job-name=download_nbl_primary_pe
#SBATCH --chdir=.
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
ROOT="${NBL_DATA_ROOT:-$REPO_ROOT/data/nbl/new_data}"
ARG="${1:-GSE100148}"

need() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: '$1' not found in PATH" >&2
        return 127
    }
}

run_one_dataset() {
    local DATASET="$1"
    local BASE="${ROOT}/${DATASET}"
    local LIST="${BASE}/all_primary_tumour_pe.txt"
    local SRA_DIR="${BASE}/sra"
    local FASTQ_DIR="${BASE}/fastq"
    local LOG_DIR="${BASE}/logs"
    local TMPDIR="${SLURM_TMPDIR:-${BASE}/stage/${SLURM_JOB_ID:-manual}}"
    local THREADS="${SLURM_CPUS_PER_TASK:-8}"

    local REMOVED_FILE="${BASE}/removed_nbl_srr_id_primary_tumour_pe.txt"
    local PREFETCH_LOG_DIR="${LOG_DIR}/prefetch"
    local FQDUMP_LOG_DIR="${LOG_DIR}/fasterq_dump"
    local PIGZ_LOG_DIR="${LOG_DIR}/pigz"

    case "$DATASET" in
        GSE100148|GSE189367|SRP409177)
            ;;
        *)
            mkdir -p "${ROOT}/logs"
            echo "ERROR: Unsupported dataset '$DATASET'" >> "${ROOT}/logs/download_nbl_primary_pe_${SLURM_JOB_ID:-manual}.err"
            echo "Supported datasets: GSE100148, GSE189367, SRP409177" >> "${ROOT}/logs/download_nbl_primary_pe_${SLURM_JOB_ID:-manual}.err"
            return 2
            ;;
    esac

    if [[ ! -d "$BASE" ]]; then
        mkdir -p "${ROOT}/logs"
        echo "ERROR: Dataset directory does not exist: $BASE" >> "${ROOT}/logs/download_nbl_primary_pe_${SLURM_JOB_ID:-manual}.err"
        return 2
    fi

    if [[ ! -s "$LIST" ]]; then
        mkdir -p "$LOG_DIR"
        echo "ERROR: List file not found or empty: $LIST" >> "${LOG_DIR}/${DATASET}_${SLURM_JOB_ID:-manual}.err"
        return 2
    fi

    mkdir -p \
        "$SRA_DIR" \
        "$FASTQ_DIR" \
        "$LOG_DIR" \
        "$TMPDIR" \
        "$PREFETCH_LOG_DIR" \
        "$FQDUMP_LOG_DIR" \
        "$PIGZ_LOG_DIR"

    touch "$REMOVED_FILE"

    local JOB_STDOUT="${LOG_DIR}/${DATASET}_${SLURM_JOB_ID:-manual}.out"
    local JOB_STDERR="${LOG_DIR}/${DATASET}_${SLURM_JOB_ID:-manual}.err"

    (
        exec >>"$JOB_STDOUT" 2>>"$JOB_STDERR"

        echo "== $(date) :: START =="
        echo "Host: $(hostname)"
        echo "[INFO] ROOT=$ROOT"
        echo "[INFO] DATASET=$DATASET"
        echo "[INFO] BASE=$BASE"
        echo "[INFO] LIST=$LIST"
        echo "[INFO] SRA_DIR=$SRA_DIR"
        echo "[INFO] FASTQ_DIR=$FASTQ_DIR"
        echo "[INFO] LOG_DIR=$LOG_DIR"
        echo "[INFO] TMPDIR=$TMPDIR"
        echo "[INFO] THREADS=$THREADS"
        echo "[INFO] REMOVED_FILE=$REMOVED_FILE"
        echo "[INFO] JOB_STDOUT=$JOB_STDOUT"
        echo "[INFO] JOB_STDERR=$JOB_STDERR"

        CONDA_SH_PATH="${CONDA_SH_PATH:-${HOME}/miniforge3/etc/profile.d/conda.sh}"
        if [[ ! -f "$CONDA_SH_PATH" ]] && command -v conda >/dev/null 2>&1; then
            CONDA_SH_PATH="$(conda info --base)/etc/profile.d/conda.sh"
        fi
        source "$CONDA_SH_PATH"
        conda activate sra3
        hash -r

        if [[ -f /etc/ssl/certs/ca-certificates.crt ]]; then
            export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
        elif [[ -f /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ]]; then
            export SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem
        fi

        need prefetch
        need fasterq-dump
        need pigz

        record_failure() {
            local srr="$1"
            local reason="$2"
            if [[ -f "$REMOVED_FILE" ]]; then
                grep -v -E "^${srr}[[:space:]]" "$REMOVED_FILE" > "${REMOVED_FILE}.tmp" || true
                mv "${REMOVED_FILE}.tmp" "$REMOVED_FILE"
            fi
            printf "%s\t%s\n" "$srr" "$reason" >> "$REMOVED_FILE"
            echo "[REMOVED] ${srr} -> ${reason}"
        }

        clear_failure_record() {
            local srr="$1"
            if [[ -f "$REMOVED_FILE" ]]; then
                grep -v -E "^${srr}[[:space:]]" "$REMOVED_FILE" > "${REMOVED_FILE}.tmp" || true
                mv "${REMOVED_FILE}.tmp" "$REMOVED_FILE"
            fi
        }

        prefetch_if_missing() {
            local srr="$1"
            local prefetch_log="${PREFETCH_LOG_DIR}/${srr}.prefetch.log"

            if [[ -d "${SRA_DIR}/${srr}" || -f "${SRA_DIR}/${srr}.sra" ]]; then
                echo "[prefetch] skip: ${srr} already present"
                return 0
            fi

            echo "[prefetch] downloading ${srr} -> ${SRA_DIR}"
            if prefetch --max-size 0 -O "$SRA_DIR" "$srr" >"$prefetch_log" 2>&1; then
                return 0
            fi

            echo "[ERROR] prefetch failed for ${srr}; see ${prefetch_log}" >&2
            return 1
        }

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

        cleanup_partial_outputs() {
            local srr="$1"
            rm -f \
                "${FASTQ_DIR}/${srr}.fastq" \
                "${FASTQ_DIR}/${srr}_1.fastq" \
                "${FASTQ_DIR}/${srr}_2.fastq" \
                "${FASTQ_DIR}/${srr}.fastq.gz" \
                "${FASTQ_DIR}/${srr}_1.fastq.gz" \
                "${FASTQ_DIR}/${srr}_2.fastq.gz"
        }

        convert_one() {
            local srr="$1"
            local src
            local fasterq_log="${FQDUMP_LOG_DIR}/${srr}.fasterq_dump.log"
            local pigz_log="${PIGZ_LOG_DIR}/${srr}.pigz.log"

            echo "---- $(date) :: START ${srr} ----"

            if [[ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" && -s "${FASTQ_DIR}/${srr}_2.fastq.gz" ]]; then
                echo "[skip] ${srr}: paired FASTQs already exist"
                clear_failure_record "$srr"
                echo "---- $(date) :: DONE ${srr} ----"
                return 0
            fi

            cleanup_partial_outputs "$srr"

            if ! prefetch_if_missing "$srr"; then
                cleanup_partial_outputs "$srr"
                record_failure "$srr" "prefetch_failed_unresolvable"
                echo "---- $(date) :: FAIL ${srr} ----"
                return 1
            fi

            src="$(resolve_sra_source "$srr")"

            echo "[fasterq-dump] SRC=${src}"
            echo "[fasterq-dump] threads=${THREADS}"
            echo "[fasterq-dump] tmp=${TMPDIR}"
            echo "[fasterq-dump] outdir=${FASTQ_DIR}"

            if ! fasterq-dump "$src" \
                --split-files \
                --skip-technical \
                --threads "$THREADS" \
                --temp "$TMPDIR" \
                --outdir "$FASTQ_DIR" >"$fasterq_log" 2>&1; then
                echo "[ERROR] fasterq-dump failed for ${srr}; see ${fasterq_log}" >&2
                cleanup_partial_outputs "$srr"
                record_failure "$srr" "fasterq_dump_failed"
                echo "---- $(date) :: FAIL ${srr} ----"
                return 1
            fi

            if [[ ! -s "${FASTQ_DIR}/${srr}_1.fastq" || ! -s "${FASTQ_DIR}/${srr}_2.fastq" ]]; then
                echo "[ERROR] expected paired FASTQs not produced for ${srr}" >&2
                cleanup_partial_outputs "$srr"
                record_failure "$srr" "paired_fastq_not_produced"
                echo "---- $(date) :: FAIL ${srr} ----"
                return 1
            fi

            echo "[pigz] compressing ${srr}_1.fastq and ${srr}_2.fastq"
            if ! pigz -p "$THREADS" "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${srr}_2.fastq" >"$pigz_log" 2>&1; then
                echo "[ERROR] pigz failed for ${srr}; see ${pigz_log}" >&2
                cleanup_partial_outputs "$srr"
                record_failure "$srr" "pigz_failed"
                echo "---- $(date) :: FAIL ${srr} ----"
                return 1
            fi

            if [[ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" && -s "${FASTQ_DIR}/${srr}_2.fastq.gz" ]]; then
                clear_failure_record "$srr"
                echo "---- $(date) :: DONE ${srr} ----"
                return 0
            else
                echo "[ERROR] pigz completed but compressed paired FASTQs missing for ${srr}" >&2
                cleanup_partial_outputs "$srr"
                record_failure "$srr" "pigz_failed"
                echo "---- $(date) :: FAIL ${srr} ----"
                return 1
            fi
        }

        TOTAL=$(grep -Evc '^[[:space:]]*(#|$)' "$LIST")
        echo "[INFO] TOTAL SRRs TO PROCESS: $TOTAL"

        success_count=0
        fail_count=0
        skip_count=0
        i=0

        while IFS= read -r line; do
            srr="$(echo "$line" | tr -d '[:space:]')"
            [[ -z "$srr" || "${srr:0:1}" == "#" ]] && continue

            i=$((i + 1))
            echo "[SEQ] ($i/$TOTAL) $srr"

            if [[ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" && -s "${FASTQ_DIR}/${srr}_2.fastq.gz" ]]; then
                echo "[skip] ${srr}: paired FASTQs already exist"
                clear_failure_record "$srr"
                skip_count=$((skip_count + 1))
                continue
            fi

            if convert_one "$srr"; then
                success_count=$((success_count + 1))
            else
                fail_count=$((fail_count + 1))
            fi
        done < "$LIST"

        echo "== $(date) :: ALL DONE =="
        echo "[SUMMARY] total=${TOTAL}"
        echo "[SUMMARY] success=${success_count}"
        echo "[SUMMARY] skipped=${skip_count}"
        echo "[SUMMARY] failed=${fail_count}"
        echo "[SUMMARY] removed_file=${REMOVED_FILE}"
    )
}

if [[ "$ARG" == "all" ]]; then
    overall_status=0
    for ds in GSE100148 GSE189367 SRP409177; do
        if ! run_one_dataset "$ds"; then
            overall_status=1
        fi
    done
    exit "$overall_status"
else
    run_one_dataset "$ARG"
    exit "$?"
fi
