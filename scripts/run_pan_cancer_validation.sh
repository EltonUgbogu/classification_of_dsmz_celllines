#!/bin/bash
# =============================================================================
# run_pan_cancer_validation.sh
# =============================================================================
#
# Master script to run Steps 2 and 3 of pan-cancer validation pipeline
#
# This script:
#   1. Builds combined expression matrix (Step 2)
#   2. Computes correlation matrix (Step 3.1)
#   3. Computes NN purity test (Step 3.2)
#   4. Builds resolved neighbourhood graph (Step 3.3)
#   5. Quantifies cross-lineage edges (Step 3.4)
#
# USAGE:
#   bash scripts/run_pan_cancer_validation.sh
#
# =============================================================================

set -e  # Exit on error

# Base paths (derive from script location for portability)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
RESULTS_DIR="${PIPELINE_DIR}/results/unsupervised/pan_cancer"

# Input files (adjust paths as needed)
BRCA_VST="${PIPELINE_DIR}/results/brca/inputs/vst_joint.rds"
NBL_VST="${PIPELINE_DIR}/results/nbl/inputs/vst_joint.rds"
RBL_VST="${PIPELINE_DIR}/results/rbl/inputs/vst_joint.rds"
HEME_VST="${PIPELINE_DIR}/results/heme/inputs/vst_joint.rds"  # Optional
METADATA="${PIPELINE_DIR}/data/dsmz/DSMZ_metadata.csv"
GENES="${RESULTS_DIR}/pan_cancer_features_clean.txt"

# Output files
EXPR_OUT="${RESULTS_DIR}/pan_cancer_expr.rds"
COR_OUT="${RESULTS_DIR}/pan_cancer_cor.rds"
NN_PURITY_OUT="${RESULTS_DIR}/nn_purity.tsv"
GRAPH_DIR="${RESULTS_DIR}/graph"
CROSS_EDGES_OUT="${RESULTS_DIR}/cross_lineage_edges.tsv"

echo "=========================================="
echo "PAN-CANCER VALIDATION PIPELINE"
echo "=========================================="
echo ""

# Step 2: Build combined expression matrix
echo "[STEP 2] Building combined expression matrix..."
if [ -f "${HEME_VST}" ]; then
  Rscript ${SCRIPTS_DIR}/build_pan_cancer_expression_matrix.R \
    --brca-vst "${BRCA_VST}" \
    --nbl-vst "${NBL_VST}" \
    --rbl-vst "${RBL_VST}" \
    --heme-vst "${HEME_VST}" \
    --metadata "${METADATA}" \
    --genes "${GENES}" \
    --output "${EXPR_OUT}"
else
  echo "  Warning: HEME VST not found, building without hematologic samples"
  Rscript ${SCRIPTS_DIR}/build_pan_cancer_expression_matrix.R \
    --brca-vst "${BRCA_VST}" \
    --nbl-vst "${NBL_VST}" \
    --rbl-vst "${RBL_VST}" \
    --metadata "${METADATA}" \
    --genes "${GENES}" \
    --output "${EXPR_OUT}"
fi

if [ ! -f "${EXPR_OUT}" ]; then
  echo "ERROR: Expression matrix not created"
  exit 1
fi
echo "  ✓ Expression matrix created: ${EXPR_OUT}"
echo ""

# Step 3.1: Compute correlation matrix
echo "[STEP 3.1] Computing Spearman correlation matrix..."
Rscript ${SCRIPTS_DIR}/compute_pan_cancer_correlation.R \
  --input "${EXPR_OUT}" \
  --output "${COR_OUT}" \
  --method spearman

if [ ! -f "${COR_OUT}" ]; then
  echo "ERROR: Correlation matrix not created"
  exit 1
fi
echo "  ✓ Correlation matrix created: ${COR_OUT}"
echo ""

# Step 3.2: Compute NN purity test
echo "[STEP 3.2] Computing nearest-neighbour purity test..."
Rscript ${SCRIPTS_DIR}/compute_nn_purity.R \
  --cor-matrix "${COR_OUT}" \
  --expr-data "${EXPR_OUT}" \
  --k 10 \
  --output "${NN_PURITY_OUT}"

if [ ! -f "${NN_PURITY_OUT}" ]; then
  echo "ERROR: NN purity table not created"
  exit 1
fi
echo "  ✓ NN purity table created: ${NN_PURITY_OUT}"
echo ""

# Step 3.3: Build resolved neighbourhood graph
echo "[STEP 3.3] Building resolved neighbourhood graph..."
mkdir -p "${GRAPH_DIR}"
Rscript ${SCRIPTS_DIR}/build_pan_cancer_graph.R \
  --cor-matrix "${COR_OUT}" \
  --k 20 \
  --output-dir "${GRAPH_DIR}"

if [ ! -f "${GRAPH_DIR}/pan_cancer_graph_edges.tsv" ]; then
  echo "ERROR: Graph edges not created"
  exit 1
fi
echo "  ✓ Graph created in: ${GRAPH_DIR}"
echo ""

# Step 3.4: Quantify cross-lineage edges
echo "[STEP 3.4] Quantifying cross-lineage edges..."
Rscript ${SCRIPTS_DIR}/quantify_cross_lineage_edges.R \
  --edges "${GRAPH_DIR}/pan_cancer_graph_edges.tsv" \
  --expr-data "${EXPR_OUT}" \
  --output "${CROSS_EDGES_OUT}"

if [ ! -f "${CROSS_EDGES_OUT}" ]; then
  echo "ERROR: Cross-lineage edges not analyzed"
  exit 1
fi
echo "  ✓ Cross-lineage analysis complete: ${CROSS_EDGES_OUT}"
echo ""

echo "=========================================="
echo "VALIDATION PIPELINE COMPLETE"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Check NN purity: cat ${NN_PURITY_OUT} | grep HEME"
echo "  2. Check cross-lineage edges: wc -l ${CROSS_EDGES_OUT}"
echo "  3. Review graph components: cat ${GRAPH_DIR}/pan_cancer_components.tsv"
echo ""
