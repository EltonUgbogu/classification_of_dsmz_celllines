#!/bin/bash
#SBATCH --job-name=multiqc_summary
#SBATCH --output=/home/chu25/dsmz/logs/multiqc/multiqc_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/multiqc/multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10:30:00

# Load conda and activate Snakemake
echo "[INFO] Activating Snakemake conda environment..."
source /home/chu25/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# Run the standalone Snakefile
snakemake \
  --snakefile /home/chu25/dsmz/scripts/multiqc/Snakefile \
  --cores 1 \
  --use-conda 
