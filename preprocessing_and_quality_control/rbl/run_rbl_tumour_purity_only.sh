#!/usr/bin/env bash
#SBATCH --job-name=rbl_purity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Avoid accidental inherited data-root variables.
unset DATA_ROOT

SMK_ENV="${SMK_ENV:-envs/.conda/smk}"

if [ ! -d "$SMK_ENV" ]; then
  echo "[ERROR] Snakemake env not found: $SMK_ENV"
  exit 1
fi

if [ -n "${CONDA_SH_PATH:-}" ]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH_PATH"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
else
  echo "[ERROR] conda not found. Set CONDA_SH_PATH or add conda to PATH." >&2
  exit 1
fi

conda activate "$SMK_ENV"

echo "[INFO] Job ID: ${SLURM_JOB_ID:-manual}"
echo "[INFO] Running RBL tumour purity only"
echo "[INFO] Snakemake: $(command -v snakemake) ($(snakemake --version))"

snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --use-conda \
  --conda-frontend conda \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --latency-wait 300 \
  -p \
  rbl_tumour_purity_analysis
