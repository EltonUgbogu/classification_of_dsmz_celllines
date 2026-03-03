#!/bin/bash

#SBATCH --job-name=nbl_fastq_2_countdata
#SBATCH --output=/home/chu25/TCGA/nbl/logs/nbl_fastq_2_countdata_%j.out
#SBATCH --error=/home/chu25/TCGA/nbl/logs/nbl_fastq_2_countdata_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00

echo "=== NBL fastq_2_countdata R job started ==="
echo "Start time: $(date)"
echo "Running in: $(pwd)"
# Setup Conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcga-r-env

# Render the R Markdown document
Rscript fastq_2_countdata.R

echo "=== Job finished at $(date) ==="