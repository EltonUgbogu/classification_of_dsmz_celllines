#!/bin/bash
#SBATCH --job-name=dsmz_tcga_correlation
#SBATCH --output=/home/chu25/dsmz/logs/dsmz_tcga_correlation/%x_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/dsmz_tcga_correlation/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

echo "=== Preliminary analysis R job started ==="
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "CWD:  $(pwd)"

# Ensure output dirs exist
mkdir -p /home/chu25/dsmz/logs/dsmz_tcga_correlation
mkdir -p /home/chu25/dsmz/results/dsmz_tcga_correlation

# Conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env

# BLAS threads = CPUs requested (prevents oversubscription)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Add TinyTeX / user bin, if you’ve installed it
export PATH="$HOME/bin:$PATH"

# Decide output format automatically:
# Use PDF if both pandoc and a LaTeX engine are available, else HTML.
FMT=pdf_document
echo "Using rmarkdown output_format: ${FMT}"

# Render
Rscript -e "rmarkdown::render(
  'dsmz_tcga_correlation.R',
  output_format = 'pdf_document',
  output_dir    = '/home/chu25/dsmz/results/dsmz_tcga_correlation',
  clean         = FALSE,
  quiet         = FALSE
)"

STATUS=$?
echo "=== Job finished at $(date) with status ${STATUS} ==="
exit ${STATUS}
