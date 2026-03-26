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

echo "== $(date) :: START =="
echo "Host: $(hostname)"
echo "[INFO] THREADS=${THREADS}"
echo "[INFO] PREFETCH=${PREFETCH}"
echo "[INFO] FASTERQ_DUMP=${FASTERQ_DUMP}"
echo "[INFO] PIGZ=${PIGZ}"

"${PREFETCH}" --version
"${FASTERQ_DUMP}" --version
"${PIGZ}" --version

for SRR in \
  SRR17010969 \
  SRR17010970 \
  SRR17010973 \
  SRR17010975 \
  SRR17010976 \
  SRR17010978 \
  SRR17010982 \
  SRR17010991
do
  echo "== Retrying ${SRR} =="

  rm -rf "sra/${SRR}"
  rm -rf "tmp/manual_retry_${SRR}"
  rm -f "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"
  rm -f "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"

  mkdir -p "tmp/manual_retry_${SRR}" sra

  if ! "${PREFETCH}" "${SRR}" --output-directory sra; then
    echo "[FAIL] prefetch ${SRR}"
    continue
  fi

  if ! "${FASTERQ_DUMP}" "sra/${SRR}" \
    --threads "${THREADS}" \
    --temp "tmp/manual_retry_${SRR}" \
    --outdir fastq
  then
    echo "[FAIL] fasterq-dump ${SRR}"
    continue
  fi

  if ! "${PIGZ}" -p "${THREADS}" "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"; then
    echo "[FAIL] pigz ${SRR}"
    continue
  fi

  ls -lh "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"
done

echo "== $(date) :: DONE =="