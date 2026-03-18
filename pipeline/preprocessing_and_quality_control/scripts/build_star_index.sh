#!/bin/bash
#SBATCH --job-name=star_index_gencode_v44
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-star-index-%j.out
#SBATCH --error=slurm-star-index-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REF_DIR="${PIPELINE_ROOT}/data/reference/gencode_v44"
INDEX_DIR="${PIPELINE_ROOT}/data/reference/star_gencode_v44"
GTF_PATH="${PIPELINE_ROOT}/data/gtf/gencode.v44.basic.annotation.gtf"
LOG_DIR="${PIPELINE_ROOT}/logs"

mkdir -p "${LOG_DIR}" "${INDEX_DIR}"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq-star

STAR --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --genomeFastaFiles "${REF_DIR}/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "${GTF_PATH}" \
  --sjdbOverhang 100 \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}"