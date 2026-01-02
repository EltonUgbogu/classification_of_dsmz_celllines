#!/bin/bash
#SBATCH --job-name=exprtf_brca
#SBATCH --partition=gpu
#SBATCH --gres=gpu:h100_94gb:1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/work/ugbogu/expr_transformer
#SBATCH --output=/work/ugbogu/expr_transformer/logs/exprtf_brca_%j.out
#SBATCH --error=/work/ugbogu/expr_transformer/logs/exprtf_brca_%j.err

set -euo pipefail
mkdir -p logs ckpt embeddings data

module load arch/r1/zen4
module load miniforge3/4.8.3-4

# conda init
export PS1=${PS1:-}
export PS1=${PS1:-}
source /srv/software/spack/packages/linux-rocky9-zen4/gcc-14.2.0/miniforge3-4.8.3-4-Linux-x86_64-wceylqrywx6ca5fhcwpemr7ishol4nrv/etc/profile.d/conda.sh

conda activate exprtf

python src/prep_rds_to_npz.py
python src/train_mgm.py
python src/export_embeddings.py
