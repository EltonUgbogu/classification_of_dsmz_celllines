#!/bin/bash
#SBATCH --job-name=rb_srr6785x
#SBATCH --chdir=/home/chu25/data/rbl
#SBATCH --output=/home/chu25/data/rbl/logs/srr6785x_%A_%a.out
#SBATCH --error=/home/chu25/data/rbl/logs/srr6785x_%A_%a.err
#SBATCH --array=1-6%3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sra3
[ -f /etc/ssl/certs/ca-certificates.crt ] && export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

BASE=/home/chu25/data/rbl
SRA_DIR=$BASE/sra
FASTQ_DIR=$BASE/fastq
TMPDIR=${SLURM_TMPDIR:-$BASE/tmp}
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMPDIR"

# list
SRRS=( SRR6785775 SRR6785776 SRR6785777 SRR6785778 SRR6785779 SRR6785780 )
SRR="${SRRS[$((SLURM_ARRAY_TASK_ID-1))]}"

# prefetch (remove size cap)
if [[ ! -d "$SRA_DIR/$SRR" && ! -f "$SRA_DIR/$SRR.sra" ]]; then
  prefetch --max-size 0 -O "$SRA_DIR" "$SRR"
fi

SRC=$SRR
[[ -d "$SRA_DIR/$SRR" ]] && SRC="$SRA_DIR/$SRR"
[[ -f "$SRA_DIR/$SRR.sra" ]] && SRC="$SRA_DIR/$SRR.sra"

fasterq-dump "$SRC" --split-files --skip-technical --threads "$SLURM_CPUS_PER_TASK" \
  --temp "$TMPDIR" --outdir "$FASTQ_DIR"

pigz -p "$SLURM_CPUS_PER_TASK" "$FASTQ_DIR/${SRR}"_*.fastq
