#!/usr/bin/env bash
# Rebuild all RBL per-direction DSMZ graph edge TSVs (canonical IDs).
# Rule: dsmz_dsmz_similarity_graph
# Output pattern (from Snakefile): {TUMOUR_NH_ROOT}/{direction}/final_consensus/DSMZ_DSMZ_graph_edges_{direction}.tsv
# Run from pipeline root: bash scripts/rebuild_rbl_direction_edges.sh

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Exact path prefix matching Snakefile (TUMOUR_NH_ROOT for RBL = results/unsupervised/rbl/tumour_neighbourhoods)
NH="results/unsupervised/rbl/tumour_neighbourhoods"

# Directions that exist under results/unsupervised/rbl/tumour_neighbourhoods/ (from config + ls)
DIRECTIONS=(
  Entropy_corr Entropy_euc MAD_corr MAD_euc MX_corr MX_euc
  MeanAbsDev_corr MeanAbsDev_euc PCA_corr PCA_euc
  Spearman_corr Spearman_euc Variance_corr Variance_euc
  HVG_euc HVG_corr kTotal_corr kTotal_euc
)

# Build target list: one edge TSV per direction (same pattern as rule output)
TARGETS=()
for d in "${DIRECTIONS[@]}"; do
  TARGETS+=( "$NH/$d/final_consensus/DSMZ_DSMZ_graph_edges_$d.tsv" )
done

echo "[INFO] Requesting ${#TARGETS[@]} targets from rule dsmz_dsmz_similarity_graph"
snakemake \
  --configfile config/config.yaml \
  --config pipeline_profile=rbl \
  --cores 1 \
  "${TARGETS[@]}"

echo "[OK] Direction edge TSVs rebuilt. Run: python3 scan_ambiguous_ids.py"
