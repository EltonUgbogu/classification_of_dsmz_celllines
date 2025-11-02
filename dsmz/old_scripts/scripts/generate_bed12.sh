#!/bin/bash
#SBATCH --job-name=gtf2bed12
#SBATCH --output=/home/chu25/dsmz/logs/gtf2bed/gtf2bed12_%j.out
#SBATCH --error=/home/chu25/dsmz/logs/gtf2bed/gtf2bed12_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Load environment
set -euo pipefail

# Activate environment
echo "[INFO] Activating ucsc-tools environment..."
source /home/chu25/miniconda3/etc/profile.d/conda.sh

conda activate ucsc-tools

# Paths
GTF="/home/chu25/data/GTF/Homo_sapiens.GRCh38.107.gtf"
BED_OUT="/home/chu25/data/hg38_bed12/hg38.refseq.bed12"

# Run conversion
gtfToGenePred "$GTF" stdout | genePredToBed stdin "$BED_OUT"

echo "[INFO] BED12 file written to $BED_OUT"

