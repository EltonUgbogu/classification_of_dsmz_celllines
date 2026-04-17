#!/usr/bin/env bash
# =============================================================================
# export_direction_edges_to_final_consensus_all.sh
# =============================================================================
#
# Exports per-direction DSMZ–DSMZ edge TSVs into final_consensus_all/ so that:
#   - All downstream artifacts live in one stable folder (same style as
#     resolve_dsmz_graph_neighbors.R → resolved_dsmz_neighbors.tsv).
#   - Optional: merge into cell_line_similarity_graph_edges_ALL.tsv.
#
# Usage:
#   ./scripts/export_direction_edges_to_final_consensus_all.sh PROFILE
#   e.g. ./scripts/export_direction_edges_to_final_consensus_all.sh rbl
#
# Writes:
#   results/unsupervised/<PROFILE>/tumour_neighbourhoods/final_consensus_all/
#     cell_line_similarity_graph_edges_<direction>.tsv   (one per direction)
#     cell_line_similarity_graph_edges_ALL.tsv           (optional merged file)
#
# Use COPY=1 for real copies instead of symlinks (safer across filesystems).
# =============================================================================

set -euo pipefail

PROFILE="${1:-}"
if [[ -z "$PROFILE" ]]; then
  echo "Usage: $0 PROFILE [COPY=0]" >&2
  echo "  PROFILE  e.g. rbl, brca" >&2
  echo "  COPY=1   copy files instead of symlinks" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="${BASE}/results/unsupervised/${PROFILE}/tumour_neighbourhoods/final_consensus_all"
NH_ROOT="${BASE}/results/unsupervised/${PROFILE}/tumour_neighbourhoods"

mkdir -p "$OUTDIR"

# Adjust glob to where per-direction edge files actually live
count=0
for src in "${NH_ROOT}"/*/final_consensus/cell_line_similarity_graph_edges_*.tsv; do
  [[ -e "$src" ]] || continue
  bn=$(basename "$src")
  if [[ "${COPY:-0}" == "1" ]]; then
    cp -f "$src" "${OUTDIR}/${bn}"
  else
    rel_src=$(python3 -c "import os; print(os.path.relpath(os.path.abspath('$src'), '$OUTDIR'))")
    ln -sf "$rel_src" "${OUTDIR}/${bn}"
  fi
  (( count++ )) || true
done

echo "[OK] Exported $count per-direction edge file(s) to ${OUTDIR}"

# Optional: merged file (requires R or Python; here a simple cat header + body)
# Uncomment and adapt if you want cell_line_similarity_graph_edges_ALL.tsv
# ALL="${OUTDIR}/cell_line_similarity_graph_edges_ALL.tsv"
# head -1 "${NH_ROOT}"/*/final_consensus/cell_line_similarity_graph_edges_*.tsv 2>/dev/null | head -1 > "$ALL"
# tail -n +2 -q "${NH_ROOT}"/*/final_consensus/cell_line_similarity_graph_edges_*.tsv 2>/dev/null >> "$ALL" || true
# echo "[OK] Merged edges: $ALL"
