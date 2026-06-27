#!/usr/bin/env bash
# Run tumour neighbourhoods for all feature-set directions and metrics.
# Fail fast if inputs missing; log each run for debugging.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$SCRIPT_DIR}"
cd "$REPO_ROOT"
set -euo pipefail

mkdir -p logs/heme/manual_tumour_neighbourhoods

for d in Entropy MAD MeanAbsDev MX PCA Spearman Variance kTotal; do
  expr="results/unsupervised/heme/feature_selection_unsupervised/featuresets/${d}/expr_submatrix.rds"
  map="results/unsupervised/heme/tumour_neighbourhoods_input/cell_line_to_original_sample_id_${d}.rds"

  # Hard checks (stop immediately if something is missing)
  [[ -s "$expr" ]] || { echo "[ERROR] Missing expr: $expr" >&2; exit 1; }
  [[ -s "$map"  ]] || { echo "[ERROR] Missing map : $map"  >&2; exit 1; }

  for metric in euc corr; do
    direction="${d}_${metric}"
    log="logs/heme/manual_tumour_neighbourhoods/comp_tumour_neighbourhoods_${direction}.log"
    echo "[RUN] $direction  (log: $log)"

    Rscript scripts/comp_tumour_neighbourhoods.R \
      --config config/config.yaml \
      --profile heme \
      --direction "$direction" \
      --expr-rds "$expr" \
      --mapping-rds "$map" \
      > "$log" 2>&1
  done
done

echo "[OK] All requested directions completed."
