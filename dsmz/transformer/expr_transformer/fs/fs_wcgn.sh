#!/bin/bash
#SBATCH --job-name=WGCNA_Analysis
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/work/ugbogu/expr_transformer/fs/logs/wgcna_%j.out
#SBATCH --error=/work/ugbogu/expr_transformer/fs/logs/wgcna_%j.err

set -euo pipefail

cd /work/ugbogu/expr_transformer/fs
mkdir -p /work/ugbogu/expr_transformer/fs/logs

CONDA_BASE="/srv/software/spack/packages/linux-rocky9-zen4/gcc-14.2.0/miniforge3-4.8.3-4-Linux-x86_64-wceylqrywx6ca5fhcwpemr7ishol4nrv"
export PS1="${PS1-}"   # Avoids "PS1: unbound variable" error under set -u in non-interactive jobs
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate r_env

start=$(date +%s)
echo "[START] $(date -Is)"
echo "[NODE] ${SLURM_NODELIST:-NA}  [CPUS] ${SLURM_CPUS_PER_TASK:-NA}  [JOBID] ${SLURM_JOB_ID:-NA}"

Rscript fs_wcgna.R

end=$(date +%s)
elapsed=$((end - start))
printf "[END]   %s\n" "$(date -Is)"
printf "[TIME]  %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
