#!/bin/bash

#train_exprtf.sh

#SBATCH --job-name=exprtf_brca
#SBATCH --output=exprtf_brca_%j.out
#SBATCH --error=exprtf_brca_%j.err
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${EXPRTF_CONDA_ENV:-exprtf}"

python "${PROJECT_ROOT}/src/train_mgm.py"
python "${PROJECT_ROOT}/src/export_embeddings.py"
