#!/bin/bash
#SBATCH --job-name=check_cov
#SBATCH --output=/home/chu25/dsmz/logs/check_coverage/check_cov_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/check_coverage/check_cov_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=300:30:00

# Load conda
echo "[INFO] Activating conda environment..."
source /home/chu25/miniconda3/etc/profile.d/conda.sh
conda activate rseqc


# === Run the Python script ===
echo "[INFO] Starting BAM coverage check..."
python3	bam_coverage.py
echo "[INFO] Finished."
