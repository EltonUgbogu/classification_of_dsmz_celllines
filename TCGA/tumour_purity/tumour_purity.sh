#!/bin/bash
#SBATCH --job-name=tcga_estimate_purity
#SBATCH --chdir=/home/chu25/TCGA/tumour_purity
#SBATCH --output=/home/chu25/TCGA/logs/tcga_estimate_purity_%x_%j.out
#SBATCH --error=/home/chu25/TCGA/logs/tcga_estimate_purity_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=168:00:00

# Activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tcga-r-env || { echo "Failed to activate Conda environment tcga-r-env"; exit 1; }

# Output folder
OUT_DIR="/home/chu25/TCGA/results/tumour_purity_results"
mkdir -p "$OUT_DIR" || { echo "Failed to create output directory $OUT_DIR"; exit 1; }

# Ensure input R script exists
R_SCRIPT="/home/chu25/TCGA/tumour_purity/tcga_estimate_purity.Rmd"
if [[ ! -f "$R_SCRIPT" ]]; then
    echo "Error: R script $R_SCRIPT not found"
    exit 1
fi

echo "Running R script $R_SCRIPT"
# Run R script to render Rmd to PDF
Rscript -e "library(rmarkdown); rmarkdown::render(
  input = '$R_SCRIPT',
  output_format = rmarkdown::pdf_document(latex_engine = 'xelatex'),
  output_dir = '$OUT_DIR',
  clean = FALSE,
  quiet = FALSE
)" || { echo 'R script execution failed'; exit 1; }

echo "Job completed successfully. Output saved in $OUT_DIR"
