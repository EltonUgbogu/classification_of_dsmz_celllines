#!/bin/bash

#train_exprtf.sh

#SBATCH --job-name=exprtf_brca
#SBATCH --output=/home/ugbogu/projects/expr_transformer/logs/exprtf_brca_%j.out
#SBATCH --error=/home/ugbogu/projects/expr_transformer/logs/exprtf_brca_%j.err
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

mkdir -p /home/ugbogu/projects/expr_transformer/logs
source ~/miniconda3/etc/profile.d/conda.sh
conda activate exprtf

python /home/ugbogu/projects/expr_transformer/src/train_mgm.py
python /home/ugbogu/projects/expr_transformer/src/export_embeddings.py
