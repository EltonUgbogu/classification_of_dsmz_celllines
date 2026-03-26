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

echo "== $(date) :: START =="
echo "Host: $(hostname)"
echo "[INFO] THREADS=${THREADS}"

for SRR in \
  SRR17010972 \
  SRR17010974 \
  SRR17010977 \
  SRR17010979 \
  SRR17010980 \
  SRR17010981 \
  SRR17010983 \
  SRR17010984 \
  SRR17010985 \
  SRR17010986 \
  SRR17010990
do
  echo "== Retrying ${SRR} =="

  rm -rf "sra/${SRR}"
  rm -rf "tmp/manual_retry_${SRR}"
  rm -f "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"
  rm -f "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"

  mkdir -p "tmp/manual_retry_${SRR}" sra

  if ! prefetch "${SRR}" --output-directory sra; then
    echo "[FAIL] prefetch ${SRR}"
    continue
  fi

  if ! fasterq-dump "sra/${SRR}" \
    --threads "${THREADS}" \
    --temp "tmp/manual_retry_${SRR}" \
    --outdir fastq
  then
    echo "[FAIL] fasterq-dump ${SRR}"
    continue
  fi

  if ! pigz -p "${THREADS}" "fastq/${SRR}_1.fastq" "fastq/${SRR}_2.fastq"; then
    echo "[FAIL] pigz ${SRR}"
    continue
  fi

  ls -lh "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"
done

echo "== $(date) :: DONE =="
