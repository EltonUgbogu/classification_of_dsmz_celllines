#!/bin/bash
#SBATCH --job-name=fix_SRR25381852
#SBATCH --output=/home/chu25/TCGA/hema_pipeline/hema_pipeline/logs/fix_SRR25381852_%j.out
#SBATCH --error=/home/chu25/TCGA/hema_pipeline/hema_pipeline/logs/fix_SRR25381852_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

cd /home/chu25/TCGA/hema_pipeline/hema_pipeline

snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --use-conda --conda-frontend conda \
  --cores 8 \
  /home/chu25/data/cll/GSE237907/results/star/SRR25381852.Aligned.sortedByCoord.out.bam \
  /home/chu25/data/cll/GSE237907/results/star/valid_count/SRR25381852.tab
