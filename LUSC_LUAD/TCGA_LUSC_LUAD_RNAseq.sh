#!/bin/bash
#SBATCH --job-name=TCGA_data_analysis
#SBATCH --output=logs/TCGA_data_analysis_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

echo "=== TCGA R job started ==="
echo "Start time: $(date)"
echo "Running in: $(pwd)"
# Setup Conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env

# Add TinyTeX bin dir to PATH
export PATH="$HOME/bin:$PATH"

# Render the R Markdown document
Rscript -e "rmarkdown::render('TCGA_bulk_RNAseq_data.Rmd', output_format = 'pdf_document')"

echo "=== Job finished at $(date) ==="

