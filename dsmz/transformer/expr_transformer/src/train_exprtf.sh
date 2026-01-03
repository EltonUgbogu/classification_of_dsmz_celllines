#!/bin/bash
#SBATCH --job-name=exprtf_brca
#SBATCH --partition=gpu
#SBATCH --qos=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/work/%u/expr_transformer/logs/%x_%j.out
#SBATCH --error=/work/%u/expr_transformer/logs/%x_%j.err

set -euo pipefail

# --- paths ---
PROJECT_ROOT="/work/$USER/expr_transformer"
LOGDIR="$PROJECT_ROOT/logs"
mkdir -p "$LOGDIR"

# --- conda (match what worked on hgc02) ---
source "/work/$USER/software/miniforge3/etc/profile.d/conda.sh"
conda activate "${EXPRTF_CONDA_ENV:-exprtf}"

# --- GPU sanity checks (fail fast) ---
nvidia-smi
python - <<'PY'
import sys
import torch

print("cuda?", torch.cuda.is_available())
print("torch_cuda", torch.version.cuda)
if not torch.cuda.is_available():
    sys.exit("ERROR: CUDA is not available. This job must run on a GPU node.")
print("gpu:", torch.cuda.get_device_name(0))
PY

cd "$PROJECT_ROOT"

python src/train_mgm.py
python src/export_embeddings.py
