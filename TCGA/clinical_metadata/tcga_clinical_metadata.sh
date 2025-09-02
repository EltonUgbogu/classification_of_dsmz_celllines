#!/bin/bash
#SBATCH --job-name=tcga_clinical_metadata
#SBATCH --output=/home/chu25/TCGA/logs/tcga_clinical_metadata_%j.log
#SBATCH --error=/home/chu25/TCGA/logs/tcga_clinical_metadata_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00

echo "=== TCGA Metadata R job started ==="
echo "Start time: $(date)"
echo "Running in: $(pwd)"

mkdir -p logs

# Setup Conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env

# Add TinyTeX bin dir to PATH
export PATH="$HOME/bin:$PATH"

# Render the R Markdown document
Rscript -e "rmarkdown::render('tcga_metadata.Rmd', output_format = 'pdf_document')"

echo "=== Job finished at $(date) ==="

