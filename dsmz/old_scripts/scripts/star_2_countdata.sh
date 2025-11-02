#!/bin/bash

#SBATCH --job-name=star_2_countdata
#SBATCH --output=/home/chu25/dsmz/logs/star_2_countdata/star_2_countdata_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/star_2_countdata/star_2_countdata_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00

echo "=== Star_2_countdata R job started ==="
echo "Start time: $(date)"
echo "Running in: $(pwd)"
# Setup Conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env

# Add TinyTeX bin dir to PATH
export PATH="$HOME/bin:$PATH"

# Render the R Markdown document
Rscript -e "rmarkdown::render('star_2_countdata.Rmd', output_format = 'pdf_document')"

echo "=== Job finished at $(date) ==="
