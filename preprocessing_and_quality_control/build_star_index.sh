#!/bin/bash
#SBATCH --job-name=star_index_gencode_v44
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/reference_build/slurm-star-index-%j.out
#SBATCH --error=logs/reference_build/slurm-star-index-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REFERENCE_ROOT="${REFERENCE_ROOT:-${PIPELINE_ROOT}/data/reference}"
REF_DIR="${REFERENCE_ROOT}/gencode_v44"
INDEX_DIR="${STAR_INDEX_DIR:-${REFERENCE_ROOT}/star_gencode_v44}"
GTF_PATH="${REF_DIR}/gencode.v44.basic.annotation.gtf"
FASTA_PATH="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
LOG_DIR="${PIPELINE_ROOT}/logs/reference_build"

mkdir -p "${LOG_DIR}" "${REF_DIR}" "${INDEX_DIR}"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq-star

if [[ ! -s "${FASTA_PATH}" ]]; then
    wget -c -O "${FASTA_PATH}.gz" \
      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
    gunzip -f "${FASTA_PATH}.gz"
fi

if [[ ! -s "${GTF_PATH}" ]]; then
    wget -c -O "${GTF_PATH}.gz" \
      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
    gunzip -f "${GTF_PATH}.gz"
fi

STAR \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --genomeFastaFiles "${FASTA_PATH}" \
  --sjdbGTFfile "${GTF_PATH}" \
  --sjdbOverhang 100 \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}" \
  --limitGenomeGenerateRAM 60000000000
