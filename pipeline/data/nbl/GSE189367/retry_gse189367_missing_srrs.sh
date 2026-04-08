#!/bin/bash
#SBATCH --job-name=retry_gse189367_missing
#SBATCH --chdir=/work/ugbogu/pipeline/data/nbl/GSE189367
#SBATCH --output=/work/ugbogu/pipeline/data/nbl/GSE189367/logs/retry_missing_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/nbl/GSE189367/logs/retry_missing_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -u -o pipefail

cd /work/ugbogu/pipeline/data/nbl/GSE189367
mkdir -p logs sra fastq tmp

THREADS="${SLURM_CPUS_PER_TASK:-8}"

PREFETCH="/work/ugbogu/.conda/envs/sra3/bin/prefetch"
FASTERQ_DUMP="/work/ugbogu/.conda/envs/sra3/bin/fasterq-dump"
PIGZ="/work/ugbogu/.conda/envs/sra3/bin/pigz"

MISSING_FILE="removed_nbl_srr_id_primary_tumour_pe.txt"
SUCCESS_LOG="logs/retry_success_$(date +%Y%m%d_%H%M%S).txt"
FAIL_LOG="logs/retry_fail_$(date +%Y%m%d_%H%M%S).txt"

echo "== $(date) :: START =="
echo "Host: $(hostname)"
echo "[INFO] THREADS=${THREADS}"
echo "[INFO] PREFETCH=${PREFETCH}"
echo "[INFO] FASTERQ_DUMP=${FASTERQ_DUMP}"
echo "[INFO] PIGZ=${PIGZ}"
echo "[INFO] MISSING_FILE=${MISSING_FILE}"

"${PREFETCH}" --version
"${FASTERQ_DUMP}" --version
"${PIGZ}" --version

if [[ ! -s "${MISSING_FILE}" ]]; then
  echo "[ERROR] Missing retry list: ${MISSING_FILE}"
  exit 1
fi

while IFS=$'\t' read -r SRR REASON; do
  [[ -z "${SRR}" ]] && continue

  echo "== Retrying ${SRR} (${REASON}) =="

  rm -rf "sra/${SRR}"
  rm -rf "tmp/manual_retry_${SRR}"
  rm -f "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"
  rm -f "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"

  mkdir -p "tmp/manual_retry_${SRR}" sra

  if ! "${PREFETCH}" --max-size 0 --output-directory sra "${SRR}"; then
    echo -e "${SRR}\tprefetch_failed" | tee -a "${FAIL_LOG}"
    continue
  fi

  SRC="sra/${SRR}"
  if [[ -f "sra/${SRR}.sra" ]]; then
    SRC="sra/${SRR}.sra"
  elif [[ -f "sra/${SRR}/${SRR}.sra" ]]; then
    SRC="sra/${SRR}"
  fi

  if ! "${FASTERQ_DUMP}" "${SRC}" \
    --split-files \
    --skip-technical \
    --threads "${THREADS}" \
    --temp "tmp/manual_retry_${SRR}" \
    --outdir fastq
  then
    echo -e "${SRR}\tfasterq_dump_failed" | tee -a "${FAIL_LOG}"
    continue
  fi

  if [[ ! -s "fastq/${SRR}_1.fastq" || ! -s "fastq/${SRR}_2.fastq" ]]; then
    echo -e "${SRR}\tpaired_fastq_not_produced" | tee -a "${FAIL_LOG}"
    continue
  fi

  if ! "${PIGZ}" -p "${THREADS}" "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"; then
    echo -e "${SRR}\tpigz_failed" | tee -a "${FAIL_LOG}"
    continue
  fi

  if [[ -s "fastq/${SRR}_1.fastq.gz" && -s "fastq/${SRR}_2.fastq.gz" ]]; then
    echo -e "${SRR}\tsuccess" | tee -a "${SUCCESS_LOG}"
    ls -lh "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"
  else
    echo -e "${SRR}\tcompressed_fastq_missing" | tee -a "${FAIL_LOG}"
  fi

done < "${MISSING_FILE}"

echo "== $(date) :: DONE =="
echo "[INFO] SUCCESS_LOG=${SUCCESS_LOG}"
echo "[INFO] FAIL_LOG=${FAIL_LOG}"