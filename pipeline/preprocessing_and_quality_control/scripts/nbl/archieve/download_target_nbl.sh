#!/bin/bash
#SBATCH --job-name=download_target_nbl
#SBATCH --chdir=/work/ugbogu/pipeline/data/nbl
#SBATCH --output=/work/ugbogu/pipeline/data/nbl/target_nbl/logs/download_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/nbl/target_nbl/logs/download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"

ROOT=/work/ugbogu/pipeline/data/nbl/target_nbl
MANIFEST="${ROOT}/manifest.tsv"
OUTDIR="${ROOT}/count_data"
LOGDIR="${ROOT}/logs"
GDC_CLIENT=/work/ugbogu/tools/gdc-client

mkdir -p "$ROOT" "$OUTDIR" "$LOGDIR"

echo "[INFO] ROOT=$ROOT"
echo "[INFO] MANIFEST=$MANIFEST"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] LOGDIR=$LOGDIR"
echo "[INFO] GDC_CLIENT=$GDC_CLIENT"

if [[ ! -x "$GDC_CLIENT" ]]; then
    echo "ERROR: gdc-client not found or not executable: $GDC_CLIENT" >&2
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: manifest.tsv not found: $MANIFEST" >&2
    exit 1
fi

cd "$ROOT"

echo "[INFO] Starting GDC download"
"$GDC_CLIENT" download \
    -m "$MANIFEST" \
    -d "$OUTDIR" \
    --n-processes 4

echo "== $(date) :: DONE =="