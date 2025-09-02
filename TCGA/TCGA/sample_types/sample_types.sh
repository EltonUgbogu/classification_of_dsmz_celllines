#!/bin/bash
#SBATCH --job-name=sample_types
#SBATCH --output=logs/sample_types_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

echo "=== Sample Types R job started ==="
echo "Start time: $(date)"
echo "Running in: $(pwd)"

mkdir -p logs

# Setup Conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env

# Add TinyTeX bin dir to PATH
export PATH="$HOME/bin:$PATH"

# Render the R Markdown document
Rscript -e "rmarkdown::render('sample_types.Rmd', output_format = 'pdf_document')"

echo "=== Job finished at $(date) ==="

