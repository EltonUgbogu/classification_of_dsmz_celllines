#!/usr/bin/env Rscript

# =============================================================================
# score_tumour_cellline_mapping.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script quantifies transcriptomic similarity between tumour samples and
# cell lines to identify which cell lines best represent the molecular profiles
# of individual tumours. The analysis employs correlation-based scoring in a
# curated pan-cancer feature space to rank cell lines by their expression
# similarity to each tumour sample.
#
# BIOLOGICAL RATIONALE
# --------------------
# Cancer cell lines serve as essential in vitro model systems for studying
# tumour biology and evaluating therapeutic responses. However, the fidelity
# with which cell lines recapitulate the molecular characteristics of primary
# tumours varies considerably. This script addresses a fundamental question in
# translational cancer research: given a patient tumour sample, which available
# cell lines most faithfully represent its transcriptomic profile?
#
# The matching problem is non-trivial because:
#
#   1. Cell line derivation and prolonged culture introduce transcriptomic
#      drift, potentially altering expression patterns relative to the
#      originating tumour tissue.
#
#   2. Tumours exhibit substantial intra- and inter-patient heterogeneity,
#      meaning no single cell line can perfectly represent all tumours of
#      a given cancer type.
#
#   3. Molecular subtypes within cancer types (e.g., PAM50 subtypes in breast
#      cancer) may be differentially represented in cell line panels.
#
# By systematically quantifying tumour-cell line similarity, this analysis
# enables evidence-based selection of model systems for downstream functional
# studies, drug screening, and mechanistic investigations.
#
#
# FEATURE SPACE DEFINITION
# ------------------------
# Rather than using the full transcriptome (~20,000 protein-coding genes), the
# analysis operates within a curated pan-cancer feature space. This approach
# offers several advantages:
#
#   Noise Reduction:
#     The transcriptome contains many genes with low or highly variable
#     expression that contribute noise rather than signal. Excluding these
#     uninformative features improves the signal-to-noise ratio in similarity
#     calculations.
#
#   Biological Focus:
#     Curated feature sets typically capture cancer-relevant pathways and
#     biological processes (e.g., proliferation, immune signalling, metabolic
#     reprogramming). This ensures that similarity scores reflect biologically
#     meaningful variation rather than technical or stochastic differences.
#
#   Cross-Cohort Comparability:
#     A common feature set enables consistent comparisons across cancer types
#     (pan-cancer analysis), facilitating identification of cell lines that
#     may serve as models for multiple tumour lineages.
#
#   Computational Efficiency:
#     Reduced dimensionality accelerates pairwise similarity calculations,
#     which scale quadratically with the number of features.
#
#
# SPEARMAN RANK CORRELATION
# -------------------------
# Similarity between tumour and cell line expression profiles is quantified
# using Spearman rank correlation (rho). This non-parametric measure offers
# key advantages for transcriptomic comparisons:
#
#   Robustness to Outliers:
#     Gene expression data frequently contain extreme values due to biological
#     variation or technical artefacts. Rank transformation mitigates the
#     influence of outliers on similarity scores.
#
#   Monotonic Relationship Detection:
#     Spearman correlation captures any monotonic association between profiles,
#     not assuming linearity. This is appropriate for expression data, where
#     relationships between genes may be non-linear.
#
#   Scale Invariance:
#     The rank-based approach is invariant to differences in expression
#     magnitude between tumours and cell lines, which may arise from
#     platform differences, normalisation procedures, or biological factors.
#
#   Distribution-Free:
#     No assumptions about normality are required, accommodating the
#     heterogeneous distributional properties of gene expression data.
#
# For tumour t and cell line c with expression vectors x_t and x_c across n
# features, Spearman correlation is computed on rank-transformed values,
# yielding rho in [-1, +1] where +1 indicates perfect positive monotonic
# association.
#
#
# RANKING AND EVALUATION METRICS
# ------------------------------
# For each tumour, cell lines are ranked by descending correlation score.
# The quality of mappings is evaluated using complementary metrics:
#
#   Top-1 Accuracy:
#     Proportion of tumours whose single best-matching cell line originates
#     from the same cancer type (lineage). This stringent metric requires
#     the highest-ranked cell line to be a "correct" match.
#
#     Interpretation: High top-1 accuracy indicates that same-lineage cell
#     lines consistently outperform cross-lineage alternatives, suggesting
#     preserved lineage-specific transcriptomic programmes.
#
#   Top-k Accuracy:
#     Proportion of tumours with at least one correct-lineage cell line
#     among the top k rankings. This relaxed metric acknowledges that
#     multiple cell lines may be reasonable models.
#
#     Interpretation: Top-k >> Top-1 suggests multiple competitive candidates,
#     while Top-k ~ Top-1 indicates a single dominant match per tumour.
#
#   Mean Reciprocal Rank (MRR):
#     For each tumour, compute 1/rank for the first correct-lineage cell line,
#     then average across all tumours. MRR rewards correct matches appearing
#     earlier in rankings:
#       - MRR = 1.0: All tumours have correct lineage at rank 1
#       - MRR = 0.5: Correct lineage consistently appears at rank 2
#       - MRR -> 0: Correct lineage appears deep in rankings
#
#     Interpretation: MRR provides a single summary statistic that captures
#     both accuracy and ranking quality.
#
# These metrics assume that tumour-cell line pairs from the same cancer type
# represent "correct" matches. This is a reasonable approximation, though
# molecular subtypes within a cancer type may be better matched by specific
# cell lines (e.g., HER2+ tumours matching HER2-amplified cell lines).
#
#
# INPUT REQUIREMENTS
# ------------------
# 1. Pan-cancer expression RDS (--pan-cancer-expr):
#      - Contains cell line expression matrix (genes x samples) and metadata
#      - Expression values should be normalised (e.g., VST, log-TPM)
#      - Metadata must include sample_id and lineage columns
#
# 2. Joint VST files (--joint-vst-brca, --joint-vst-nbl, --joint-vst-rbl):
#      - Batch-corrected expression matrices containing both cell lines and
#        tumours from the same cancer type
#      - Joint normalisation ensures tumour and cell line expression scales
#        are directly comparable
#      - Tumour samples identified by exclusion from cell line metadata
#
# 3. Pan-cancer feature genes (--genes):
#      - Text file with one gene identifier per line (Ensembl IDs)
#      - Typically derived from variance-based or biological feature selection
#      - Version suffixes (e.g., .16) are stripped for compatibility
#
#
# OUTPUT FILES
# ------------
# The script produces several output files for downstream analysis:
#
#   tumour_to_cellline_rankings.tsv
#       Full ranking table with columns: tumour, cell_line, rank, score,
#       cell_lineage, tum_lineage, is_correct. Contains top-k cell lines
#       per tumour for detailed inspection of matches.
#
#   tumour_mapping_summary.tsv
#       Per-tumour summary including best match, confidence metrics
#       (score delta, fraction same lineage), and top-k cell line lists.
#       Facilitates identification of tumours with clear vs. ambiguous matches.
#
#   metrics_summary.tsv
#       Overall evaluation metrics (top-1 accuracy, top-k accuracy, MRR).
#       Provides aggregate assessment of mapping quality.
#
#   metrics_by_lineage.tsv
#       Per-cancer-type breakdown of mapping accuracy. Reveals whether
#       certain lineages are better represented in the cell line panel.
#
#   tumour_cellline_scores.rds
#       Full tumour-by-cell-line score matrix (tumours x cell lines).
#       Enables downstream analyses such as hierarchical clustering,
#       visualisation, or alternative ranking strategies.
#
#   Extended metrics files (from Step 5.5):
#       - metrics_summary_extended.tsv: Comprehensive tumour-driven metrics
#       - metrics_by_tumour_lineage_extended.tsv: Per-lineage detailed metrics
#       - cellline_metrics_by_tumour_cohort.tsv: Cell-line-driven metrics
#       - cellline_metrics_overall.tsv: Pooled cell line performance
#         (includes topk_lineage_entropy and top1_margin_z)
#       - cellline_topk_lineage_enrichment.tsv: Lineage fractions in top-k tumours
#       - tumour_cellline_scores_long.tsv.gz: Long-format score table
#
#
# USAGE EXAMPLE
# -------------
#   Rscript score_tumour_cellline_mapping.R \
#     --pan-cancer-expr results/unsupervised/pan_cancer/pan_cancer_expr.rds \
#     --joint-vst-brca data/brca/bcc_analysis_results/vst_dsmz_bcc_tcga_brca_batch_corrected.rds \
#     --joint-vst-nbl data/nbl/preprocessing_results/vst_NBL_joint_batch_corrected.rds \
#     --joint-vst-rbl data/rbl/preprocessing_results/vst_RBL_joint_batch_corrected.rds \
#     --genes results/unsupervised/pan_cancer/pan_cancer_features_clean.txt \
#     --output-dir results/unsupervised/pan_cancer/tumour_mapping \
#     --top-k 10
#
#
# DEPENDENCIES
# ------------
# R packages:
#   - data.table: High-performance data manipulation for large matrices
#   - optparse:   Command-line argument parsing
#
# =============================================================================


# -----------------------------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------------------------
# Package loading is wrapped in suppressPackageStartupMessages() to provide
# clean console output, particularly important when the script is executed
# as part of automated pipelines (e.g., Snakemake workflows).

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})


# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# The script uses the optparse package for robust command-line argument
# handling. This approach enables:
#
#   - Integration with workflow managers (Snakemake, Nextflow)
#   - Self-documenting parameter requirements
#   - Default values for optional parameters
#   - Automatic help message generation (--help)
#
# Arguments are organised into logical groups:
#
#   Input Data:
#     --pan-cancer-expr: Cell line expression and metadata (RDS format)
#     --joint-vst-*:     Batch-corrected joint expression per cancer type
#     --genes:           Pan-cancer feature gene list (text file)
#
#   Output:
#     --output-dir:      Directory for all output files (created if absent)
#
#   Parameters:
#     --top-k:            Number of top cell lines to rank per tumour
#     --similarity:       spearman or cosine similarity metric
#     --zscore-cohort:    Gene-wise z-scoring within cell lines and tumours
#     --drop-global-genes:Remove ribosomal/mitochondrial/histone genes (if symbols)
#     --global-gene-regex:Regex for global-program gene symbols
#     --gene-symbol-map:  Optional Ensembl -> symbol map for filtering
#     --use-signatures:   Placeholder for optional signature scoring

option_list <- list(
  make_option(c("--pan-cancer-expr"), type = "character", default = NULL,
              help = "Path to pan-cancer expression RDS (cell lines only)"),
  make_option(c("--joint-vst-brca"), type = "character", default = NULL,
              help = "Path to BRCA joint VST (cell lines + tumours)"),
  make_option(c("--joint-vst-nbl"), type = "character", default = NULL,
              help = "Path to NBL joint VST (cell lines + tumours)"),
  make_option(c("--joint-vst-rbl"), type = "character", default = NULL,
              help = "Path to RBL joint VST (cell lines + tumours)"),
  make_option(c("--genes"), type = "character", default = NULL,
              help = "Path to pan-cancer feature gene list"),
  make_option(c("--output-dir"), type = "character", default = NULL,
              help = "Output directory for results"),
  make_option(c("--top-k"), type = "integer", default = 10,
              help = "Number of top cell lines to rank (default: 10)"),
  make_option(c("--similarity"), type = "character", default = "spearman",
              help = "Similarity metric: spearman or cosine (default: spearman)"),
  make_option(c("--zscore-cohort"), type = "logical", default = TRUE,
              help = "Gene-wise z-score within cohorts (cell lines and tumours)"),
  make_option(c("--drop-global-genes"), type = "logical", default = TRUE,
              help = "Drop ribosomal/mitochondrial/histone genes when symbols available"),
  make_option(c("--global-gene-regex"), type = "character",
              default = "^RPL|^RPS|^MT-|^HIST",
              help = "Regex for global-program gene symbols to drop"),
  make_option(c("--gene-symbol-map"), type = "character", default = NULL,
              help = "Optional TSV/CSV mapping with columns ensembl_id and symbol"),
  make_option(c("--use-signatures"), action = "store_true", default = FALSE,
              help = "Also compute signature scores (anchor/isolate programmes)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

valid_similarity <- c("spearman", "cosine")
opt$similarity <- tolower(opt$similarity)
if (!opt$similarity %in% valid_similarity) {
  stop("--similarity must be one of: ", paste(valid_similarity, collapse = ", "))
}


# -----------------------------------------------------------------------------
# Input Validation
# -----------------------------------------------------------------------------
# Critical inputs are validated before proceeding to prevent cryptic errors
# during execution. The script requires:
#
#   1. Feature-restricted pan-cancer expression (--pan-cancer-expr)
#   2. Output directory (--output-dir)
#
# Joint VST inputs are retained for backwards compatibility but ignored in the
# default config-driven pipeline.

if (is.null(opt$`pan-cancer-expr`) || is.null(opt$`output-dir`)) {
  stop("--pan-cancer-expr and --output-dir are required")
}

# Create output directory with recursive = TRUE to handle nested paths
# showWarnings = FALSE suppresses messages if directory already exists
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

# Create subfolders for direction-specific similarity outputs
t2c_dir <- file.path(opt$`output-dir`, "tumour_to_cellline_similarity")
c2t_dir <- file.path(opt$`output-dir`, "cellline_to_tumour_similarity")
dir.create(t2c_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(c2t_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("TUMOUR - CELL-LINE MAPPING\n")
cat("========================================\n\n")


# =============================================================================
# STEP 1: Load Feature-Restricted Pan-Cancer Expression Object
# =============================================================================
cat("[1] Loading pan-cancer expression object...\n")
pan_obj <- readRDS(opt$`pan-cancer-expr`)
if (!is.list(pan_obj) || is.null(pan_obj$expr) || is.null(pan_obj$meta)) {
  stop("--pan-cancer-expr must be an RDS containing list(expr=..., meta=...)")
}
expr_all <- pan_obj$expr
meta_all <- as.data.table(pan_obj$meta)

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

required_meta_cols <- c("sample_id", "lineage", "type")
missing <- setdiff(required_meta_cols, colnames(meta_all))
if (length(missing) > 0) {
  stop("Pan-cancer metadata missing columns: ", paste(missing, collapse = ", "))
}

if (!all(meta_all$sample_id %in% colnames(expr_all))) {
  bad <- setdiff(meta_all$sample_id, colnames(expr_all))
  stop("Metadata sample IDs not found in expression matrix: ",
       paste(head(bad, 20), collapse = ", "))
}
expr_all <- expr_all[, meta_all$sample_id, drop = FALSE]

if (!is.null(pan_obj$genes)) {
  feature_genes <- strip_version(pan_obj$genes)
  cat("  Feature genes from RDS: ", length(feature_genes), "\n", sep = "")
} else if (!is.null(opt$genes) && nzchar(opt$genes)) {
  fallback_genes <- readLines(opt$genes)
  fallback_genes <- fallback_genes[fallback_genes != ""]
  feature_genes <- strip_version(fallback_genes)
  cat("  Feature genes from --genes: ", length(feature_genes), "\n", sep = "")
} else {
  feature_genes <- strip_version(rownames(expr_all))
  cat("  Feature genes inferred from expression matrix: ", length(feature_genes), "\n", sep = "")
}


# =============================================================================
# STEP 2: Split Cell Lines and Tumours
# =============================================================================
cat("\n[2] Splitting cell lines and tumours...\n")
type_norm <- tolower(meta_all$type)
type_norm[type_norm %in% c("cellline", "cell line", "cells", "cell_lines")] <- "cell_line"
type_norm[type_norm %in% c("tumor", "tumours", "tumor_sample", "tumour")] <- "tumour"
meta_all[, type := type_norm]

meta_cell <- meta_all[type == "cell_line"]
meta_tum <- meta_all[type == "tumour"]

if (nrow(meta_cell) == 0) {
  stop("Pan-cancer expression object contains no cell line samples")
}

expr_cell <- expr_all[, meta_cell$sample_id, drop = FALSE]
expr_tum <- expr_all[, meta_tum$sample_id, drop = FALSE]

cat("  Cell lines: ", ncol(expr_cell), " samples\n", sep = "")
cat("  Cell line lineages:\n")
print(table(meta_cell$lineage))


if (nrow(meta_tum) == 0) {
  cat("  No tumour samples present in pan-cancer expression object. Writing empty outputs...\n")

  empty_rank_t2c <- data.table(
    tumour = character(),
    rank = integer(),
    cell_line = character(),
    score = numeric(),
    cell_lineage = character(),
    tum_lineage = character(),
    is_correct = logical()
  )
  empty_summary_t2c <- data.table(
    tumour_id = character(),
    tumour_lineage = character(),
    best_cell_line = character(),
    best_score = numeric(),
    best_cell_line_lineage = character(),
    is_correct = logical(),
    topk_cell_lines = character(),
    topk_scores = character(),
    topk_lineages = character(),
    n_correct_in_topk = integer(),
    confidence_delta = numeric(),
    confidence_frac_same_lineage = numeric()
  )
  empty_rank_c2t <- data.table(
    cell_line = character(),
    rank = integer(),
    tumour = character(),
    score = numeric(),
    tumour_lineage = character(),
    cell_lineage = character(),
    is_correct = logical()
  )
  empty_summary_c2t <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    best_tumour = character(),
    best_score = numeric(),
    best_tumour_lineage = character(),
    is_correct = logical(),
    topk_tumours = character(),
    topk_scores = character(),
    topk_lineages = character(),
    n_correct_in_topk = integer(),
    confidence_delta = numeric(),
    confidence_frac_same_lineage = numeric()
  )

  empty_metrics <- data.table(metric = character(), value = numeric())
  empty_metrics_by_lineage <- data.table(
    lineage = character(),
    n_tumours = integer(),
    top1_accuracy = numeric()
  )
  empty_metrics_by_lineage_cl <- data.table(
    lineage = character(),
    n_cell_lines = integer(),
    top1_accuracy = numeric()
  )

  empty_scores_long_t2c <- data.table(
    tumour = character(),
    tum_lineage = character(),
    cell_line = character(),
    cell_lineage = character(),
    score = numeric(),
    rank_in_tumour = integer(),
    in_topk = logical()
  )
  empty_scores_long_c2t <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    tumour = character(),
    tum_lineage = character(),
    score = numeric(),
    rank_in_tumour = integer(),
    in_topk = logical()
  )

  empty_metrics_extended <- data.table(metric = character(), value = numeric())
  empty_metrics_by_tumour_lineage <- data.table(
    tum_lineage = character(),
    n_tumours = integer(),
    top1_accuracy = numeric(),
    topk_accuracy = numeric(),
    mrr = numeric(),
    mean_best_score = numeric(),
    median_best_score = numeric(),
    mean_delta_s1_sk = numeric(),
    median_delta_s1_sk = numeric(),
    mean_frac_same_lineage_topk = numeric(),
    median_frac_same_lineage_topk = numeric()
  )
  empty_cellline_by_cohort <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    tum_lineage = character(),
    n_tumours = integer(),
    mean_score = numeric(),
    median_score = numeric(),
    mean_rank = numeric(),
    median_rank = numeric(),
    pct_in_top1 = numeric(),
    pct_in_topk = numeric(),
    best_tumour = character(),
    best_score = numeric(),
    best_rank = integer()
  )
  empty_cellline_overall <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    n_tumours = integer(),
    mean_score = numeric(),
    median_score = numeric(),
    mean_rank = numeric(),
    median_rank = numeric(),
    pct_in_top1 = numeric(),
    pct_in_topk = numeric(),
    best_tumour = character(),
    best_score = numeric(),
    best_rank = integer(),
    best_cohort = character(),
    topk_lineage_entropy = numeric(),
    top1_margin_z = numeric()
  )
  empty_topk_enrich <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    topk_lineage_entropy = numeric()
  )
  empty_topk_consistency <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    n_total = integer(),
    n_same_lineage = integer(),
    frac_same = numeric(),
    all_same = logical()
  )
  empty_low_confidence <- data.table(
    cell_line = character(),
    cell_lineage = character(),
    best_tumour_lineage = character(),
    is_correct = logical(),
    confidence_delta = numeric(),
    confidence_frac_same_lineage = numeric()
  )

  fwrite(empty_rank_t2c, file = file.path(t2c_dir, "tumour_to_cellline_rankings.tsv"), sep = "\t")
  fwrite(empty_summary_t2c, file = file.path(t2c_dir, "tumour_mapping_summary.tsv"), sep = "\t")
  fwrite(empty_metrics, file = file.path(t2c_dir, "metrics_summary.tsv"), sep = "\t")
  fwrite(empty_metrics_by_lineage, file = file.path(t2c_dir, "metrics_by_lineage.tsv"), sep = "\t")
  fwrite(empty_metrics_extended, file = file.path(t2c_dir, "metrics_summary_extended.tsv"), sep = "\t")
  fwrite(empty_metrics_by_tumour_lineage, file = file.path(t2c_dir, "metrics_by_tumour_lineage_extended.tsv"), sep = "\t")
  fwrite(empty_scores_long_t2c, file = file.path(t2c_dir, "tumour_cellline_scores_long.tsv.gz"), sep = "\t")
  saveRDS(matrix(numeric(0), nrow = 0, ncol = 0),
          file = file.path(t2c_dir, "tumour_cellline_scores.rds"))

  fwrite(empty_rank_c2t, file = file.path(c2t_dir, "cellline_to_tumour_rankings.tsv"), sep = "\t")
  fwrite(empty_summary_c2t, file = file.path(c2t_dir, "cellline_mapping_summary.tsv"), sep = "\t")
  fwrite(empty_metrics, file = file.path(c2t_dir, "metrics_summary.tsv"), sep = "\t")
  fwrite(empty_metrics_by_lineage_cl, file = file.path(c2t_dir, "metrics_by_lineage.tsv"), sep = "\t")
  fwrite(empty_scores_long_c2t, file = file.path(c2t_dir, "cellline_tumour_scores_long.tsv.gz"), sep = "\t")
  saveRDS(matrix(numeric(0), nrow = 0, ncol = 0),
          file = file.path(c2t_dir, "cellline_tumour_scores.rds"))

  fwrite(empty_cellline_by_cohort, file = file.path(c2t_dir, "cellline_metrics_by_tumour_cohort.tsv"), sep = "\t")
  fwrite(empty_cellline_overall, file = file.path(c2t_dir, "cellline_metrics_overall.tsv"), sep = "\t")
  fwrite(empty_topk_enrich, file = file.path(c2t_dir, "cellline_topk_lineage_enrichment.tsv"), sep = "\t")
  fwrite(empty_topk_consistency, file = file.path(c2t_dir, "topk_lineage_consistency.tsv"), sep = "\t")
  fwrite(empty_low_confidence, file = file.path(c2t_dir, "low_confidence_cases.tsv"), sep = "\t")

  quit(status = 0)
}

# Validate tumour metadata schema
req_cols_t <- c("sample_id", "lineage")
missing_t <- setdiff(req_cols_t, colnames(meta_tum))
if (length(missing_t) > 0) {
  stop("meta_tum is missing required columns: ", paste(missing_t, collapse = ", "))
}

# Validate expression-metadata alignment for tumours
if (!all(meta_tum$sample_id %in% colnames(expr_tum))) {
  bad <- setdiff(meta_tum$sample_id, colnames(expr_tum))
  stop("These meta_tum$sample_id are missing in expr_tum colnames: ",
       paste(head(bad, 20), collapse = ", "),
       if (length(bad) > 20) paste0(" ... (+", length(bad) - 20, " more)") else "")
}

# Ensure column order matches metadata
expr_tum <- expr_tum[, meta_tum$sample_id, drop = FALSE]
stopifnot(identical(colnames(expr_tum), meta_tum$sample_id))

cat("\n  Total tumours: ", ncol(expr_tum), "\n", sep = "")
cat("  Tumour lineages:\n")
print(table(meta_tum$lineage))


# =============================================================================
# STEP 4: Filter to Pan-Cancer Features and Align Gene Sets
# =============================================================================
# The expression matrices are filtered to retain only genes present in all
# three sources:
#
#   1. The pan-cancer feature list (curated discriminative genes)
#   2. The cell line expression matrix (reference database)
#   3. The tumour expression matrix (query samples)
#
# Gene Set Intersection Rationale:
# This three-way intersection ensures that similarity calculations are
# performed on identical gene sets, avoiding several potential issues:
#
#   Missing Data Handling:
#     Genes absent from one matrix would require imputation or exclusion.
#     Imputation introduces assumptions that may bias similarity scores.
#
#   Annotation Consistency:
#     Different RNA-seq processing pipelines may use different gene
#     annotations, leading to non-overlapping identifiers.
#
#   Comparability:
#     Using the same genes for all comparisons ensures scores are directly
#     comparable across tumours and cell lines.
#
# The intersection operation may reduce the feature set substantially if
# annotations differ between datasets. The script reports the number of
# common genes to enable assessment of data compatibility.

cat("\n[4] Filtering to pan-cancer features...\n")

common_genes <- intersect(rownames(expr_cell), rownames(expr_tum))
common_genes <- intersect(common_genes, feature_genes)

if (length(common_genes) == 0) {
  stop("No common genes shared by cell lines, tumours, and the feature list")
}
# Subset matrices to common genes
# drop = FALSE ensures matrix structure is preserved even with one sample
expr_cell_filt <- expr_cell[common_genes, , drop = FALSE]
expr_tum_filt <- expr_tum[common_genes, , drop = FALSE]

cat("  Common genes: ", length(common_genes), "\n", sep = "")
cat("  Cell line matrix: ", nrow(expr_cell_filt), " genes x ",
    ncol(expr_cell_filt), " samples\n", sep = "")
cat("  Tumour matrix: ", nrow(expr_tum_filt), " genes x ",
    ncol(expr_tum_filt), " samples\n", sep = "")


# =============================================================================
# STEP 4.5: Optional Global-Gene Filtering + Cohort Z-Scoring
# =============================================================================
# Global transcriptional programs (ribosomal, mitochondrial, histone) can
# dominate similarity scores and pull unrelated lineages together. This step
# optionally removes those genes when gene symbols are available, then performs
# gene-wise z-scoring separately within cell lines and tumours to reduce
# cohort-level shifts.

cat("\n[4.5] Optional global-gene filtering and cohort z-scoring...\n")

# Helper: read gene symbol map if provided
read_symbol_map <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  map <- fread(path)
  cols <- tolower(names(map))
  ens_col <- which(cols %in% c("ensembl_id", "ensembl", "gene_id", "geneid", "ensg", "ensembl_gene_id"))
  sym_col <- which(cols %in% c("symbol", "gene_symbol", "genesymbol", "hgnc_symbol"))
  if (length(ens_col) == 0 || length(sym_col) == 0) {
    stop("--gene-symbol-map must contain columns for Ensembl ID and gene symbol")
  }
  ens_col <- ens_col[1]
  sym_col <- sym_col[1]
  setnames(map, c(names(map)[ens_col], names(map)[sym_col]), c("ensembl_id", "symbol"))
  map[, ensembl_id := sub("\\.[0-9]+$", "", ensembl_id)]
  unique(map[, .(ensembl_id, symbol)])
}

# Helper: gene-wise z-score
zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sigma <- apply(mat, 1, sd, na.rm = TRUE)
  sigma[is.na(sigma) | sigma == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sigma, "/")
}

# Optional global-gene filtering
if (isTRUE(opt$`drop-global-genes`)) {
  gene_symbols <- rownames(expr_cell_filt)
  looks_ensembl <- any(grepl("^ENSG", gene_symbols))
  symbol_map <- read_symbol_map(opt$`gene-symbol-map`)
  if (!is.null(symbol_map)) {
    gene_symbols <- symbol_map$symbol[match(sub("\\.[0-9]+$", "", rownames(expr_cell_filt)),
                                           symbol_map$ensembl_id)]
  } else if (looks_ensembl) {
    cat("  [WARN] Ensembl IDs detected but no --gene-symbol-map provided; skipping global-gene filter.\n")
  }

  if (!looks_ensembl || !is.null(symbol_map)) {
    drop_regex <- opt$`global-gene-regex`
    drop_idx <- !is.na(gene_symbols) & grepl(drop_regex, gene_symbols)
    n_drop <- sum(drop_idx)
    if (n_drop > 0) {
      expr_cell_filt <- expr_cell_filt[!drop_idx, , drop = FALSE]
      expr_tum_filt  <- expr_tum_filt[!drop_idx, , drop = FALSE]
      cat("  Dropped ", n_drop, " global-program genes using regex: ", drop_regex, "\n", sep = "")
    } else {
      cat("  No global-program genes matched regex: ", drop_regex, "\n", sep = "")
    }
  }
}

# Optional cohort-wise z-scoring
if (isTRUE(opt$`zscore-cohort`)) {
  expr_cell_filt <- zscore_by_gene(expr_cell_filt)
  expr_tum_filt  <- zscore_by_gene(expr_tum_filt)
  cat("  Applied gene-wise z-scoring within cell lines and tumours\n")
}

# =============================================================================
# STEP 5: Compute Similarity Scores
# =============================================================================
# Similarity is computed between each tumour sample and each cell line.
# The method is controlled by --similarity:
#   - spearman: rank-based correlation (robust to outliers)
#   - cosine:   normalized dot product (stable after z-scoring)
#
# Matrix Orientation:
# The cor() / cosine implementations operate on columns of two matrices.
# With expr_tum_filt (genes x tumours) and expr_cell_filt (genes x cell lines),
# the resulting score matrix has:
#   - Rows: Tumour samples (from first matrix columns)
#   - Columns: Cell line samples (from second matrix columns)
#
# Interpretation:
#   - Higher values indicate greater transcriptomic similarity.
#   - Scores are comparable within a run using a fixed preprocessing + method.
#
# Computational Complexity:
# Computing n_tumour x n_cellline similarities each over n_genes features
cat("\n[5] Computing similarity scores...\n")

if (opt$similarity == "spearman") {
  # cor() with method = "spearman" performs rank transformation internally
  # Result: tumours (rows) x cell lines (columns)
  score_cor <- cor(expr_tum_filt, expr_cell_filt, method = "spearman")
} else {
  # Cosine similarity on column vectors (genes x samples)
  tum_norm <- sqrt(colSums(expr_tum_filt^2, na.rm = TRUE))
  cell_norm <- sqrt(colSums(expr_cell_filt^2, na.rm = TRUE))
  tum_norm[tum_norm == 0 | is.na(tum_norm)] <- NA_real_
  cell_norm[cell_norm == 0 | is.na(cell_norm)] <- NA_real_
  score_cor <- crossprod(expr_tum_filt, expr_cell_filt)
  denom <- outer(tum_norm, cell_norm, "*")
  score_cor <- score_cor / denom
}

cat("  Score matrix: ", nrow(score_cor), " tumours x ",
    ncol(score_cor), " cell lines\n", sep = "")
cat("  Score range: [", round(min(score_cor, na.rm = TRUE), 3), ", ",
    round(max(score_cor, na.rm = TRUE), 3), "]\n", sep = "")


# =============================================================================
# STEP 5.5: Extended Metrics (Tumour-driven + Cell-line-driven)
# =============================================================================
# This section computes comprehensive metrics from the score matrix, analysing
# the tumour-cell line relationships from two complementary perspectives:
#
# A) Tumour-Driven Metrics:
#    Answer the question "For a given tumour, how well do cell lines match?"
#    - Top-1/Top-k accuracy: Does the best match come from the correct lineage?
#    - MRR: How highly ranked is the first correct-lineage cell line?
#    - Confidence metrics: How much better is the top match than alternatives?
#
# B) Cell-Line-Driven Metrics:
#    Answer the question "For a given cell line, how well does it represent tumours?"
#    - Mean/median rank: On average, where does this cell line rank?
#    - Selection frequency: How often is this cell line in the top-k?
#    - Best tumour match: Which tumour is most similar to this cell line?
#
# These dual perspectives enable both:
#   - Selection of appropriate cell lines for a specific tumour (personalised medicine)
#   - Assessment of which cell lines are broadly representative (general models)

cat("\n[5.5] Computing extended tumour- and cell-line-driven metrics...\n")

top_k_val <- opt$`top-k`


# -----------------------------------------------------------------------------
# Convert Score Matrix to Long Format
# -----------------------------------------------------------------------------
# The wide score matrix (tumours x cell lines) is converted to long format
# (one row per tumour-cell line pair) using data.table's as.table() conversion.
# Long format enables efficient grouped operations using data.table syntax.

dt <- as.data.table(as.table(score_cor))
setnames(dt, c("tumour", "cell_line", "score"))

# Create lookup maps for lineage annotation
# Named vectors enable O(1) lookup by sample ID
meta_cell_map <- setNames(meta_cell$lineage, meta_cell$sample_id)
meta_tum_map  <- setNames(meta_tum$lineage,  meta_tum$sample_id)

# Annotate each tumour-cell line pair with lineage information
dt[, tum_lineage  := meta_tum_map[tumour]]
dt[, cell_lineage := meta_cell_map[cell_line]]

# Compute rank of each cell line within each tumour
# frank() with ties.method = "average" handles tied scores appropriately
# Negative score ensures higher scores get lower (better) ranks
dt[, rank_in_tumour := frank(-score, ties.method = "average"), by = tumour]


# -----------------------------------------------------------------------------
# (A) Tumour-Driven Metrics
# -----------------------------------------------------------------------------
# These metrics evaluate mapping quality from the tumour's perspective:
# "Given this tumour, how well can we identify an appropriate cell line model?"


# --- Top-1 Analysis ---
# For each tumour, identify the single best-matching cell line (highest score)
# Using .I[which.max(score)] ensures robust tie handling

top1_dt <- dt[dt[, .I[which.max(score)], by = tumour]$V1]
top1_dt[, is_correct := (cell_lineage == tum_lineage)]


# --- Top-k Analysis ---
# Flag cell lines within the top-k rankings for each tumour
dt[, in_topk := rank_in_tumour <= top_k_val]


# =============================================================================
# PLOT: ECDF of tumour-level ranks for top BRCA model candidates
# =============================================================================
# This plot supports the statistical story: "choose the model with the most
# consistently low ranks across tumours" (median_rank + top-k frequency).

clean_cl_label <- function(x) {
  x <- gsub("^NG-[0-9]+_", "", x)  # drop NG-xxxxx_ prefix
  x <- gsub("_lib.*$", "", x)     # drop _lib... suffix
  x <- gsub("_", "-", x)          # nicer separators
  x
}

draw_ecdf_polyline <- function(x, col, lwd = 3, lty = 1, max_pts = 250) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(invisible(NULL))

  e <- ecdf(x)

  # Always anchor at rank 1 so all curves start from (1, 0),
  # even if the cell line never achieves rank 1 (e.g. negative controls)
  xs <- sort(unique(c(1, x)))

  # Downsample if too many unique x values (prevents "broken" PDF rendering)
  if (length(xs) > max_pts) {
    xs <- as.numeric(quantile(x, probs = seq(0, 1, length.out = max_pts), names = FALSE, type = 8))
    xs <- sort(unique(c(1, xs)))
  }

  ys <- e(xs)
  lines(xs, ys, type = "l", col = col, lwd = lwd, lty = lty)
}

plot_ecdf_rank_brca <- function(dt, outdir, top_k_val, n_show = 4, add_neg_control = TRUE) {

  # Keep only BRCA tumours
  dt_brca <- dt[tum_lineage == "BRCA" & !is.na(rank_in_tumour)]

  # Candidate pool: BRCA cell lines only
  pool <- dt_brca[cell_lineage == "BRCA",
                  .(
                    n_tumours = .N,
                    median_rank = median(rank_in_tumour, na.rm = TRUE),
                    median_score = median(score, na.rm = TRUE),
                    pct_in_topk = mean(rank_in_tumour <= top_k_val, na.rm = TRUE)
                  ),
                  by = cell_line][order(median_rank, -median_score)]

  if (nrow(pool) == 0) {
    cat("  [WARN] No BRCA candidates found for ECDF plot.\n")
    return(invisible(NULL))
  }

  # Pick top N by median_rank (tie-break by median_score already in ordering)
  candidates <- head(pool$cell_line, n_show)

  # Optional: add one negative control (best non-BRCA by median_rank)
  neg_added <- FALSE
  if (isTRUE(add_neg_control)) {
    neg <- dt_brca[cell_lineage != "BRCA",
                   .(median_rank = median(rank_in_tumour, na.rm = TRUE)),
                   by = cell_line][order(median_rank)]
    if (nrow(neg) > 0) {
      candidates <- c(candidates, neg$cell_line[1])
      neg_added <- TRUE
    }
  }

  # --- Styling ---
  cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#666666")  # 4 candidates + neg ctrl
  if (length(candidates) > length(cols)) {
    cols <- rep(cols, length.out = length(candidates))
  }
  lwd_main <- 3
  lwd_neg  <- 2
  lty_neg  <- 1

  # Prepare output path
  figdir <- file.path(outdir, "figures")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(figdir, sprintf("Fig_BRCA_rank_ECDF_top%d_k%d.pdf", n_show, top_k_val))

  # Compute max rank for x-axis limit
  max_rank <- ceiling(max(dt_brca[cell_line %in% candidates, rank_in_tumour], na.rm = TRUE))

  pdf(outfile, width = 7.5, height = 5.5)
  on.exit(dev.off(), add = TRUE)

  # Plot first ECDF (continuous polyline)
  first <- candidates[1]
  x1 <- dt_brca[cell_line == first, rank_in_tumour]

  plot(NULL,
       xlim = c(1, max_rank),
       ylim = c(0, 1),
       main = sprintf("BRCA tumours (n=%d): ECDF of tumour-level ranks",
                      length(unique(dt_brca$tumour))),
       xlab = "Cell line rank per breast cancer clinical patient sample (1 = best)",
       ylab = "Fraction of tumours with rank \u2264 x")

  draw_ecdf_polyline(x1, col = cols[1], lwd = lwd_main, lty = 1)

  # Optional: log unique ranks (BRCA often has many -> downsampling helps PDF)
  cat("  BRCA unique ranks for ", first, ": ", length(unique(x1)), "\n", sep = "")

  # Add others
  if (length(candidates) > 1) {
    for (i in 2:length(candidates)) {
      cl <- candidates[i]
      x <- dt_brca[cell_line == cl, rank_in_tumour]
      is_neg <- neg_added && (i == length(candidates))  # last is neg control if added
      draw_ecdf_polyline(x,
                         col = cols[i],
                         lwd = if (is_neg) lwd_neg else lwd_main,
                         lty = 1)
    }
  }

  # Vertical line at top-k threshold
  abline(v = top_k_val, lty = 2, col = "grey50")

  # Horizontal line at y = 0.5 (median rank reference)
  # Its intersection with each curve shows the median rank visually
  abline(h = 0.5, lty = 2, col = "grey50")

  # Legend: cell line name only (no stats); sorted by median_rank (best first).
  # Negative control gets an explicit label instead of its raw cell line ID.
  leg_candidates <- if (neg_added) candidates[-length(candidates)] else candidates
  leg_dt <- pool[cell_line %in% leg_candidates,
                 .(cell_line, median_rank)][order(median_rank)]
  leg_labels <- clean_cl_label(leg_dt$cell_line)
  leg_cols   <- cols[match(leg_dt$cell_line, candidates)]
  leg_lwd    <- rep(lwd_main, nrow(leg_dt))
  leg_lty    <- rep(1, nrow(leg_dt))

  # Append negative control entry last
  if (neg_added) {
    leg_labels <- c(leg_labels, "Negative control (non-BRCA)")
    leg_cols   <- c(leg_cols, cols[length(candidates)])
    leg_lwd    <- c(leg_lwd, lwd_neg)
    leg_lty    <- c(leg_lty, lty_neg)
  }

  legend("bottomright",
         legend = leg_labels, bty = "n",
         col = leg_cols,
         lwd = leg_lwd,
         lty = leg_lty)
  cat("  Saved ECDF rank plot:", outfile, "\n")
}

# Call ECDF plot after dt and top-k flags exist
plot_ecdf_rank_brca(dt = dt, outdir = c2t_dir, top_k_val = top_k_val,
                    n_show = 4, add_neg_control = TRUE)

plot_ecdf_rank_lineage <- function(dt, lineage, outdir, top_k_val,
                                   n_show = 4, add_neg_control = TRUE) {

  dt_lin <- dt[tum_lineage == lineage & !is.na(rank_in_tumour)]

  pool <- dt_lin[cell_lineage == lineage,
                 .(
                   n_tumours    = .N,
                   median_rank  = median(rank_in_tumour, na.rm = TRUE),
                   median_score = median(score, na.rm = TRUE),
                   pct_in_topk  = mean(rank_in_tumour <= top_k_val, na.rm = TRUE)
                 ),
                 by = cell_line][order(median_rank, -median_score)]

  if (nrow(pool) == 0) {
    cat("  [WARN] No", lineage, "candidates found for ECDF plot.\n")
    return(invisible(NULL))
  }

  candidates <- head(pool$cell_line, n_show)

  neg_added <- FALSE
  if (isTRUE(add_neg_control)) {
    neg <- dt_lin[cell_lineage != lineage,
                  .(median_rank = median(rank_in_tumour, na.rm = TRUE)),
                  by = cell_line][order(median_rank)]
    if (nrow(neg) > 0) {
      candidates <- c(candidates, neg$cell_line[1])
      neg_added <- TRUE
    }
  }

  # --- Styling ---
  cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#666666")  # 4 candidates + neg ctrl
  if (length(candidates) > length(cols)) {
    cols <- rep(cols, length.out = length(candidates))
  }
  lwd_main <- 3
  lwd_neg  <- 2
  lty_neg  <- 1

  figdir <- file.path(outdir, "figures")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(figdir, sprintf("Fig_%s_rank_ECDF_top%d_k%d.pdf",
                                       lineage, n_show, top_k_val))

  max_rank <- ceiling(max(dt_lin[cell_line %in% candidates, rank_in_tumour], na.rm = TRUE))

  pdf(outfile, width = 7.5, height = 5.5)
  on.exit(dev.off(), add = TRUE)

  first <- candidates[1]
  x1 <- dt_lin[cell_line == first, rank_in_tumour]

  plot(NULL,
       xlim = c(1, max_rank),
       ylim = c(0, 1),
       main = sprintf("%s tumours (n=%d): ECDF of tumour-level ranks",
                      lineage, length(unique(dt_lin$tumour))),
       xlab = sprintf("Cell line rank per %s clinical patient sample (1 = best)",
                      switch(lineage,
                             "NBL" = "neuroblastoma",
                             "RBL" = "retinoblastoma",
                             lineage)),
       ylab = "Fraction of tumours with rank \u2264 x")

  draw_ecdf_polyline(x1, col = cols[1], lwd = lwd_main, lty = 1)

  if (length(candidates) > 1) {
    for (i in 2:length(candidates)) {
      cl <- candidates[i]
      x <- dt_lin[cell_line == cl, rank_in_tumour]
      is_neg <- neg_added && (i == length(candidates))  # last is neg control if added
      draw_ecdf_polyline(x,
                         col = cols[i],
                         lwd = if (is_neg) lwd_neg else lwd_main,
                         lty = 1)
    }
  }

  abline(v = top_k_val, lty = 2, col = "grey50")

  # Horizontal line at y = 0.5 (median rank reference)
  abline(h = 0.5, lty = 2, col = "grey50")

  # Legend: cell line name only (no stats); sorted by median_rank (best first).
  # Negative control gets an explicit label instead of its raw cell line ID.
  leg_candidates <- if (neg_added) candidates[-length(candidates)] else candidates
  leg_dt <- pool[cell_line %in% leg_candidates,
                 .(cell_line, median_rank)][order(median_rank)]
  leg_labels <- clean_cl_label(leg_dt$cell_line)
  leg_cols   <- cols[match(leg_dt$cell_line, candidates)]
  leg_lwd    <- rep(lwd_main, nrow(leg_dt))
  leg_lty    <- rep(1, nrow(leg_dt))

  if (neg_added) {
    leg_labels <- c(leg_labels, sprintf("Negative control (non-%s)", lineage))
    leg_cols   <- c(leg_cols, cols[length(candidates)])
    leg_lwd    <- c(leg_lwd, lwd_neg)
    leg_lty    <- c(leg_lty, lty_neg)
  }

  legend("bottomright",
         legend = leg_labels, bty = "n",
         col = leg_cols,
         lwd = leg_lwd,
         lty = leg_lty)
  cat("  Saved ECDF rank plot:", outfile, "\n")
}

plot_ecdf_rank_lineage(dt, "NBL",  outdir = c2t_dir, top_k_val = top_k_val,
                       n_show = 4, add_neg_control = TRUE)
plot_ecdf_rank_lineage(dt, "RBL",  outdir = c2t_dir, top_k_val = top_k_val,
                       n_show = 4, add_neg_control = TRUE)

# Per-tumour top-k summary:
#   any_correct_topk: Is there at least one correct-lineage cell line in top-k?
#   frac_same_lineage_topk: What fraction of top-k are from correct lineage?
tumour_hit_dt <- dt[in_topk == TRUE, .(
  any_correct_topk = any(cell_lineage == tum_lineage, na.rm = TRUE),
  frac_same_lineage_topk = mean(cell_lineage == tum_lineage, na.rm = TRUE)
), by = tumour]


# --- Mean Reciprocal Rank (MRR) ---
# For each tumour, find the rank of the first correct-lineage cell line
# If no correct match exists, reciprocal rank is 0

mrr_dt <- dt[, .(
  first_correct_rank = {
    # Get ranks of all correct-lineage cell lines
    rr <- rank_in_tumour[cell_lineage == tum_lineage]
    # Return minimum rank (first/best correct match) or NA if none
    if (length(rr) == 0) NA_real_ else min(rr, na.rm = TRUE)
  }
), by = tumour]
# Compute reciprocal rank (0 if no correct match found)
mrr_dt[, rr := ifelse(is.na(first_correct_rank), 0, 1 / first_correct_rank)]


# --- Confidence Metrics ---
# Confidence delta: Difference between best score (s1) and k-th best score (sk)
# Large delta indicates a clear best match; small delta suggests ambiguity

conf_dt <- dt[, {
  ord <- order(-score)  # Sort by descending score
  s1 <- score[ord][1]   # Best (highest) score
  # k-th score, or last available if fewer than k cell lines
  sk <- score[ord][min(top_k_val, .N)]
  list(s1 = s1, sk = sk, delta_s1_sk = s1 - sk)
}, by = tumour]


# --- Merge Tumour-Level Summaries ---
# Combine all per-tumour metrics into a single comprehensive table

tumour_level <- merge(
  top1_dt[, .(tumour, tum_lineage = tum_lineage, top1_cell_line = cell_line,
              top1_score = score, top1_cell_lineage = cell_lineage, 
              top1_correct = is_correct)],
  tumour_hit_dt, by = "tumour", all.x = TRUE
)
tumour_level <- merge(tumour_level, mrr_dt[, .(tumour, rr)], 
                      by = "tumour", all.x = TRUE)
tumour_level <- merge(tumour_level, conf_dt[, .(tumour, delta_s1_sk)], 
                      by = "tumour", all.x = TRUE)


# --- Overall Tumour-Driven Metrics ---
# Aggregate metrics across all tumours for global assessment

metrics_overall <- data.table(
  metric = c(
    "top1_accuracy",
    paste0("top", top_k_val, "_accuracy"),
    "mrr",
    "mean_best_score",
    "median_best_score",
    "mean_delta_s1_sk",
    "median_delta_s1_sk",
    paste0("mean_frac_same_lineage_top", top_k_val),
    paste0("median_frac_same_lineage_top", top_k_val)
  ),
  value = c(
    mean(tumour_level$top1_correct, na.rm = TRUE),
    mean(tumour_level$any_correct_topk, na.rm = TRUE),
    mean(tumour_level$rr, na.rm = TRUE),
    mean(tumour_level$top1_score, na.rm = TRUE),
    median(tumour_level$top1_score, na.rm = TRUE),
    mean(tumour_level$delta_s1_sk, na.rm = TRUE),
    median(tumour_level$delta_s1_sk, na.rm = TRUE),
    mean(tumour_level$frac_same_lineage_topk, na.rm = TRUE),
    median(tumour_level$frac_same_lineage_topk, na.rm = TRUE)
  )
)

metrics_overall_file <- file.path(t2c_dir, "metrics_summary_extended.tsv")
fwrite(metrics_overall, metrics_overall_file, sep = "\t")
cat("  Saved:", metrics_overall_file, "\n")


# --- Metrics by Tumour Lineage ---
# Break down metrics by cancer type to identify lineages with better/worse
# cell line representation

metrics_by_tumour_lineage <- tumour_level[, .(
  n_tumours = .N,
  top1_accuracy = mean(top1_correct, na.rm = TRUE),
  topk_accuracy = mean(any_correct_topk, na.rm = TRUE),
  mrr = mean(rr, na.rm = TRUE),
  mean_best_score = mean(top1_score, na.rm = TRUE),
  median_best_score = median(top1_score, na.rm = TRUE),
  mean_delta_s1_sk = mean(delta_s1_sk, na.rm = TRUE),
  median_delta_s1_sk = median(delta_s1_sk, na.rm = TRUE),
  mean_frac_same_lineage_topk = mean(frac_same_lineage_topk, na.rm = TRUE),
  median_frac_same_lineage_topk = median(frac_same_lineage_topk, na.rm = TRUE)
), by = tum_lineage][order(tum_lineage)]

metrics_by_lineage_file <- file.path(t2c_dir, 
                                      "metrics_by_tumour_lineage_extended.tsv")
fwrite(metrics_by_tumour_lineage, metrics_by_lineage_file, sep = "\t")
cat("  Saved:", metrics_by_lineage_file, "\n")


# -----------------------------------------------------------------------------
# (B) Cell-Line-Driven Metrics
# -----------------------------------------------------------------------------
# These metrics evaluate each cell line's utility as a tumour model:
# "How well does this cell line represent tumours across the cohort?"
#
# This perspective is valuable for:
#   - Identifying broadly representative cell lines for general studies
#   - Detecting cell lines that may have drifted from tumour-like profiles
#   - Prioritising cell lines for expansion of model panels


# --- Flag Top-1 Matches ---
# For downstream analysis, mark which cell line is the top match for each tumour

top1_pairs <- top1_dt[, .(tumour, cell_line)]
top1_pairs[, is_top1 := TRUE]
setkey(top1_pairs, tumour, cell_line)

dt[, is_top1_for_tumour := FALSE]
setkey(dt, tumour, cell_line)
dt[top1_pairs, is_top1_for_tumour := TRUE]


# --- Cell Line Metrics by Tumour Cohort ---
# For each cell line, compute performance metrics separately for each tumour
# cohort (lineage). This reveals whether cell lines preferentially match
# their own lineage or cross-lineage tumours.

cellline_by_cohort <- dt[, .(
  n_tumours = .N,
  mean_score = mean(score, na.rm = TRUE),
  median_score = median(score, na.rm = TRUE),
  mean_rank = mean(rank_in_tumour, na.rm = TRUE),
  median_rank = median(rank_in_tumour, na.rm = TRUE),
  pct_in_top1 = mean(is_top1_for_tumour, na.rm = TRUE),  # Fraction of times selected as top-1
  pct_in_topk = mean(in_topk, na.rm = TRUE),              # Fraction of times in top-k
  best_tumour = tumour[which.max(score)],                 # Most similar tumour
  best_score = max(score, na.rm = TRUE),                  # Highest similarity score
  best_rank = rank_in_tumour[which.max(score)]            # Rank for best tumour
), by = .(cell_line, cell_lineage, tum_lineage)][order(tum_lineage, -mean_score)]

cellline_by_cohort_file <- file.path(c2t_dir,
                                      "cellline_metrics_by_tumour_cohort.tsv")
fwrite(cellline_by_cohort, cellline_by_cohort_file, sep = "\t")
cat("  Saved:", cellline_by_cohort_file, "\n")


# --- Cell Line Metrics Overall ---
# Aggregate cell line performance across all tumours, regardless of lineage

cellline_overall <- dt[, .(
  n_tumours = .N,
  mean_score = mean(score, na.rm = TRUE),
  median_score = median(score, na.rm = TRUE),
  mean_rank = mean(rank_in_tumour, na.rm = TRUE),
  median_rank = median(rank_in_tumour, na.rm = TRUE),
  pct_in_top1 = mean(is_top1_for_tumour, na.rm = TRUE),
  pct_in_topk = mean(in_topk, na.rm = TRUE),
  best_tumour = tumour[which.max(score)],
  best_score = max(score, na.rm = TRUE),
  best_rank = rank_in_tumour[which.max(score)],
  best_cohort = tum_lineage[which.max(score)]  # Lineage of best-matching tumour
), by = .(cell_line, cell_lineage)][order(-mean_score)]

# --- Cell Line Top-k Lineage Enrichment + Entropy ---
# For each cell line, compute the lineage composition of its top-k tumour matches.

cellline_topk_list <- lapply(colnames(score_cor), function(cid) {
  scores <- score_cor[, cid]
  ord <- order(scores, decreasing = TRUE, na.last = NA)
  top_n <- min(top_k_val, length(ord))
  if (top_n == 0) return(NULL)
  top_idx <- ord[seq_len(top_n)]
  top_tums <- rownames(score_cor)[top_idx]
  data.table(
    cell_line = cid,
    cell_lineage = meta_cell_map[cid],
    tumour = top_tums,
    tum_lineage = meta_tum_map[top_tums],
    score = scores[top_tums],
    rank = seq_len(top_n)
  )
})
cellline_topk <- rbindlist(cellline_topk_list, fill = TRUE)

if (nrow(cellline_topk) > 0) {
  n_topk <- cellline_topk[, .(n_topk = .N), by = .(cell_line, cell_lineage)]
  topk_counts <- cellline_topk[, .(n_in_topk = .N), by = .(cell_line, cell_lineage, tum_lineage)]
  topk_counts <- merge(topk_counts, n_topk, by = c("cell_line", "cell_lineage"))
  topk_counts[, frac_in_topk := n_in_topk / n_topk]

  entropy_dt <- topk_counts[, .(
    topk_lineage_entropy = -sum(frac_in_topk * log(frac_in_topk))
  ), by = .(cell_line, cell_lineage)]

  topk_enrich_file <- file.path(c2t_dir, "cellline_topk_lineage_enrichment.tsv")
  fwrite(topk_counts, topk_enrich_file, sep = "\t")
  cat("  Saved:", topk_enrich_file, "\n")
} else {
  entropy_dt <- data.table(cell_line = character(), cell_lineage = character(),
                           topk_lineage_entropy = numeric())
}

# --- Cell Line Margin Z-Score (top1 vs distribution) ---
cellline_margin <- data.table(
  cell_line = colnames(score_cor),
  top1_margin_z = vapply(colnames(score_cor), function(cid) {
    scores <- score_cor[, cid]
    if (all(is.na(scores))) return(NA_real_)
    mu <- mean(scores, na.rm = TRUE)
    sig <- sd(scores, na.rm = TRUE)
    if (is.na(sig) || sig == 0) return(NA_real_)
    (max(scores, na.rm = TRUE) - mu) / sig
  }, numeric(1))
)

# Merge entropy + margin into overall cell line metrics
cellline_overall <- merge(cellline_overall, entropy_dt, by = c("cell_line", "cell_lineage"), all.x = TRUE)
cellline_overall <- merge(cellline_overall, cellline_margin, by = "cell_line", all.x = TRUE)

cellline_overall_file <- file.path(c2t_dir, "cellline_metrics_overall.tsv")
fwrite(cellline_overall, cellline_overall_file, sep = "\t")
cat("  Saved:", cellline_overall_file, "\n")


# --- Save Long-Format Score Tables ---
# Compressed output for detailed downstream analysis or debugging

# Tumour -> cell line (primary direction)
dt_file <- file.path(t2c_dir, "tumour_cellline_scores_long.tsv.gz")
fwrite(dt[, .(tumour, tum_lineage, cell_line, cell_lineage, score,
              rank_in_tumour, in_topk)],
       dt_file, sep = "\t")
cat("  Saved:", dt_file, "\n")

# Cell line -> tumour (flipped direction)
dt_rev <- dt[, .(
  cell_line,
  cell_lineage,
  tumour,
  tum_lineage,
  score,
  rank_in_tumour,
  in_topk
)]
dt_rev_file <- file.path(c2t_dir, "cellline_tumour_scores_long.tsv.gz")
fwrite(dt_rev, dt_rev_file, sep = "\t")
cat("  Saved:", dt_rev_file, "\n")

cat("  [5.5] Extended metrics computation complete.\n")


# =============================================================================
# STEP 6: Rank Top-k Cell Lines per Tumour
# =============================================================================
# For each tumour, cell lines are ranked by descending correlation score.
# The top-k cell lines represent the most transcriptomically similar models
# and serve as candidate model systems for functional studies.
#
# Ranking vs. Thresholding:
# Ranking is preferred over absolute score thresholding because:
#
#   Score Dependency on Data Characteristics:
#     Absolute correlation values depend on the feature set size, expression
#     normalisation, and inherent biological variability. A correlation of
#     0.7 may indicate excellent matching in one context but mediocre
#     matching in another.
#
#   Cross-Tumour Comparability:
#     Relative rankings are comparable across tumours regardless of their
#     absolute similarity levels. A tumour with a challenging profile
#     (low correlations overall) can still have a meaningful "best" match.
#
#   Actionable Recommendations:
#     Top-k lists provide directly usable recommendations for cell line
#     selection, facilitating translation to experimental validation.

cat("\n[6] Ranking top-k cell lines per tumour...\n")

# Helper function to extract top-k cell line names by score
top_k <- function(score_row, k) {
  ord <- order(score_row, decreasing = TRUE)
  names(score_row)[ord][seq_len(min(k, length(ord)))]
}

# Apply to each tumour (row of score matrix)
rank_list <- lapply(rownames(score_cor), function(tid) {
  scores <- score_cor[tid, ]
  tops <- top_k(scores, k = opt$`top-k`)
  data.table(
    tumour = tid,
    rank = seq_along(tops),
    cell_line = tops,
    score = scores[tops]
  )
})

# Combine into single data.table
rank_df <- rbindlist(rank_list)


# -----------------------------------------------------------------------------
# Add Lineage Information for Evaluation
# -----------------------------------------------------------------------------
# Lineage (cancer type) information enables evaluation of mapping accuracy
# under the assumption that tumour-cell line pairs from the same cancer type
# represent "correct" matches.
#
# Limitations of Lineage-Based Evaluation:
# This is a coarse-grained assessment. More refined evaluation would consider:
#   - Molecular subtypes (e.g., PAM50 for breast cancer)
#   - Driver mutation status (e.g., MYCN amplification in neuroblastoma)
#   - Clinical features (grade, stage, treatment response)
#
# However, comprehensive molecular annotations are often unavailable for both
# tumours and cell lines, making lineage matching a practical compromise.

meta_cell_map <- setNames(meta_cell$lineage, meta_cell$sample_id)
meta_tum_map <- setNames(meta_tum$lineage, meta_tum$sample_id)

rank_df[, cell_lineage := meta_cell_map[cell_line]]
rank_df[, tum_lineage := meta_tum_map[tumour]]
rank_df[, is_correct := cell_lineage == tum_lineage]

cat("  Generated rankings for", length(unique(rank_df$tumour)), "tumours\n")

# -----------------------------------------------------------------------------
# Cell line -> tumour rankings (reverse direction)
# -----------------------------------------------------------------------------
rank_list_cl <- lapply(colnames(score_cor), function(cid) {
  scores <- score_cor[, cid]
  tops <- top_k(scores, k = opt$`top-k`)
  data.table(
    cell_line = cid,
    rank = seq_along(tops),
    tumour = tops,
    score = scores[tops]
  )
})

rank_df_cl <- rbindlist(rank_list_cl)
rank_df_cl[, tumour_lineage := meta_tum_map[tumour]]
rank_df_cl[, cell_lineage := meta_cell_map[cell_line]]
rank_df_cl[, is_correct := tumour_lineage == cell_lineage]

cat("  Generated rankings for", length(unique(rank_df_cl$cell_line)), "cell lines\n")

# -----------------------------------------------------------------------------
# Evaluation-only rankings (exclude HEME/OTHER from BRCA/NBL/RBL summaries)
# -----------------------------------------------------------------------------
eval_lineages <- c("BRCA", "NBL", "RBL")
eval_cell_ids <- meta_cell[lineage %in% eval_lineages, sample_id]
eval_tumour_ids <- meta_tum[lineage %in% eval_lineages, sample_id]

if (length(eval_cell_ids) == 0 || length(eval_tumour_ids) == 0) {
  stop("No BRCA/NBL/RBL samples available for evaluation metrics.")
}

score_cor_eval <- score_cor[eval_tumour_ids, eval_cell_ids, drop = FALSE]

rank_list_eval <- lapply(rownames(score_cor_eval), function(tid) {
  scores <- score_cor_eval[tid, ]
  tops <- top_k(scores, k = opt$`top-k`)
  data.table(
    tumour = tid,
    rank = seq_along(tops),
    cell_line = tops,
    score = scores[tops]
  )
})

rank_df_eval <- rbindlist(rank_list_eval)
rank_df_eval[, cell_lineage := meta_cell_map[cell_line]]
rank_df_eval[, tum_lineage := meta_tum_map[tumour]]
rank_df_eval[, is_correct := cell_lineage == tum_lineage]

rank_list_cl_eval <- lapply(colnames(score_cor_eval), function(cid) {
  scores <- score_cor_eval[, cid]
  tops <- top_k(scores, k = opt$`top-k`)
  data.table(
    cell_line = cid,
    rank = seq_along(tops),
    tumour = tops,
    score = scores[tops]
  )
})

rank_df_cl_eval <- rbindlist(rank_list_cl_eval)
rank_df_cl_eval[, tumour_lineage := meta_tum_map[tumour]]
rank_df_cl_eval[, cell_lineage := meta_cell_map[cell_line]]
rank_df_cl_eval[, is_correct := tumour_lineage == cell_lineage]

# Filter NBL to only DSMZ cell lines (starting with NG-) for consistency with figure validation
rank_df_cl_eval <- rank_df_cl_eval[
  !(cell_lineage == "NBL" & !grepl("^NG-", cell_line))
]


# =============================================================================
# STEP 7: Compute Evaluation Metrics
# =============================================================================
# Three complementary metrics quantify overall mapping quality across the
# full tumour cohort. These metrics enable comparison between different
# matching strategies, feature sets, or cell line panels.
#
# Metric Interpretation Guidelines:
#
#   Top-1 Accuracy (Stringent):
#     The proportion of tumours whose single best-matching cell line is from
#     the correct lineage. This is the most stringent metric.
#
#     Interpretation:
#       > 80%: Excellent - same-lineage cell lines consistently dominate
#       60-80%: Good - majority of tumours match correctly
#       < 60%: Poor - cross-lineage matching is common
#
#   Top-k Accuracy (Relaxed):
#     The proportion of tumours with at least one correct-lineage cell line
#     in the top-k rankings.
#
#     Interpretation:
#       Top-k >> Top-1: Multiple competitive candidates exist
#       Top-k ~ Top-1: Single dominant match per tumour
#
#   Mean Reciprocal Rank (MRR):
#     Average of 1/rank for the first correct match. MRR integrates both
#     accuracy (whether correct match exists) and ranking quality (how
#     highly it ranks).
#
#     Interpretation:
#       MRR > 0.8: Correct matches typically in top 1-2 positions
#       MRR 0.5-0.8: Correct matches typically in top 2-5 positions
#       MRR < 0.5: Correct matches often ranked lower than position 5

cat("\n[7] Computing evaluation metrics...\n")

# Top-1 accuracy: fraction of tumours with correct lineage at rank 1
top1 <- rank_df_eval[rank == 1]
top1_acc <- mean(top1$is_correct, na.rm = TRUE)

# Top-k accuracy: fraction of tumours with any correct lineage in top-k
topk_acc <- mean(tapply(rank_df_eval$is_correct, rank_df_eval$tumour, any), na.rm = TRUE)

# Mean Reciprocal Rank: average of 1/rank for first correct match
mrr <- mean(tapply(rank_df_eval$is_correct, rank_df_eval$tumour, function(v) {
  hit <- which(v)[1]  # Position of first TRUE
  if (is.na(hit)) 0 else 1 / hit  # 0 if no correct match in top-k
}), na.rm = TRUE)

# -----------------------------------------------------------------------------
# Per-Lineage Metrics
# -----------------------------------------------------------------------------
# Breaking down metrics by cancer type reveals whether certain lineages
# are better or worse represented in the cell line panel. This can guide:
#   - Prioritisation of cell line development for underrepresented lineages
#   - Identification of lineages where cross-lineage models may be acceptable
#   - Understanding of lineage-specific transcriptomic distinctiveness

lineage_metrics <- lapply(eval_lineages, function(lin) {
  tum_subset <- top1[tum_lineage == lin]
  if (nrow(tum_subset) == 0) return(NULL)
  data.table(
    lineage = lin,
    n_tumours = nrow(tum_subset),
    top1_accuracy = mean(tum_subset$is_correct, na.rm = TRUE)
  )
})
lineage_metrics <- rbindlist(lineage_metrics)

cat("  Top-1 accuracy: ", round(top1_acc * 100, 2), "%\n", sep = "")
cat("  Top-", opt$`top-k`, " accuracy: ", round(topk_acc * 100, 2), "%\n", sep = "")
cat("  MRR: ", round(mrr, 3), "\n", sep = "")


# =============================================================================
# STEP 8: Create Per-Tumour Summary Table
# =============================================================================
# The summary table provides a concise overview of each tumour's mapping
# results, facilitating:
#
#   - Identification of tumours with clear best matches (high confidence)
#   - Detection of tumours with ambiguous matches (low confidence)
#   - Selection of appropriate cell line models for specific tumours
#
# Confidence Metrics Explanation:
#
#   Score Delta (s1 - sk):
#     The difference between the best score and the k-th best score.
#     Large deltas indicate a clear best match that substantially outperforms
#     alternatives. Small deltas suggest multiple similarly-matched cell lines.
#
#     Interpretation:
#       delta > 0.1: Strong confidence in top match
#       delta 0.05-0.1: Moderate confidence
#       delta < 0.05: Multiple competitive candidates
#
#   Fraction Same Lineage:
#     The proportion of top-k cell lines from the correct (same) cancer type.
#     High values indicate robust lineage-specific matching; low values
#     suggest the tumour profile is better matched by cross-lineage cell lines.
#
#     Interpretation:
#       > 0.8: Strong lineage-specific signal
#       0.5-0.8: Mixed lineage matching
#       < 0.5: Cross-lineage cell lines may be more representative

cat("\n[8] Creating tumour summary table...\n")

# Base summary: best match information for each tumour
tumour_summary <- rank_df[rank == 1, .(
  tumour_id = tumour,
  tumour_lineage = tum_lineage,
  best_cell_line = cell_line,
  best_score = score,
  best_cell_line_lineage = cell_lineage,
  is_correct = is_correct
)]

# Add top-k cell lines as semicolon-separated lists for compact representation
topk_cells <- rank_df[rank <= opt$`top-k`, .(
  topk_cell_lines = paste(cell_line, collapse = ";"),
  topk_scores = paste(round(score, 3), collapse = ";"),
  topk_lineages = paste(cell_lineage, collapse = ";"),
  n_correct_in_topk = sum(is_correct)
), by = tumour]

tumour_summary <- merge(tumour_summary, topk_cells,
                        by.x = "tumour_id", by.y = "tumour", all.x = TRUE)

# Compute confidence metrics
# Improved delta calculation handles cases where top_k > number of cell lines
confidence <- rank_df[, {
  ord <- order(-score)
  s1 <- score[ord][1]
  # Use min() to handle cases with fewer cell lines than top_k
  sk <- score[ord][min(opt$`top-k`, .N)]
  list(s1 = s1, sk = sk, confidence_delta = s1 - sk)
}, by = tumour]

# Safer confidence fraction calculation (avoid row order dependency)
frac_same <- rank_df[rank <= opt$`top-k`,
  .(confidence_frac_same_lineage = mean(is_correct, na.rm = TRUE)),
  by = tumour]
confidence <- merge(confidence, frac_same, by = "tumour", all.x = TRUE)

tumour_summary <- merge(
  tumour_summary,
  confidence[, .(tumour, confidence_delta, confidence_frac_same_lineage)],
  by.x = "tumour_id", by.y = "tumour", all.x = TRUE
)

# -----------------------------------------------------------------------------
# Cell line summary table (reverse direction)
# -----------------------------------------------------------------------------
cat("\n[8.1] Creating cell line summary table...\n")

cellline_summary <- rank_df_cl[rank == 1, .(
  cell_line = cell_line,
  cell_lineage = cell_lineage,
  best_tumour = tumour,
  best_score = score,
  best_tumour_lineage = tumour_lineage,
  is_correct = is_correct
)]

topk_tumours <- rank_df_cl[rank <= opt$`top-k`, .(
  topk_tumours = paste(tumour, collapse = ";"),
  topk_scores = paste(round(score, 3), collapse = ";"),
  topk_lineages = paste(tumour_lineage, collapse = ";"),
  n_correct_in_topk = sum(is_correct)
), by = cell_line]

cellline_summary <- merge(cellline_summary, topk_tumours,
                          by = "cell_line", all.x = TRUE)

confidence_cl <- rank_df_cl[, {
  ord <- order(-score)
  s1 <- score[ord][1]
  sk <- score[ord][min(opt$`top-k`, .N)]
  list(s1 = s1, sk = sk, confidence_delta = s1 - sk)
}, by = cell_line]

frac_same_cl <- rank_df_cl[rank <= opt$`top-k`,
  .(confidence_frac_same_lineage = mean(is_correct, na.rm = TRUE)),
  by = cell_line]
confidence_cl <- merge(confidence_cl, frac_same_cl, by = "cell_line", all.x = TRUE)

cellline_summary <- merge(
  cellline_summary,
  confidence_cl[, .(cell_line, confidence_delta, confidence_frac_same_lineage)],
  by = "cell_line", all.x = TRUE
)

# -----------------------------------------------------------------------------
# Cell line -> tumour evaluation metrics (reverse direction)
# -----------------------------------------------------------------------------
rank_df_cl_eval[, cell_line_group := ifelse(
  cell_lineage == "RBL",
  cell_line,
  sub("_lib.*$", "", cell_line)
)]

# Collapse library-level replicates before C2T evaluation metrics
rank_df_cl_collapsed <- rank_df_cl_eval[
  order(cell_line_group, -score),
  .SD[1],
  by = .(cell_line_group, tumour)
]
rank_df_cl_collapsed[, rank_eval := frank(-score, ties.method = "first"),
                     by = cell_line_group]
rank_df_cl_collapsed[, is_correct := (tumour_lineage == cell_lineage)]

top1_cl <- rank_df_cl_collapsed[rank_eval == 1]
top1_acc_cl <- mean(top1_cl$is_correct, na.rm = TRUE)

topk_acc_cl <- mean(tapply(rank_df_cl_collapsed$is_correct,
                           rank_df_cl_collapsed$cell_line_group, any),
                    na.rm = TRUE)

mrr_cl <- mean(tapply(rank_df_cl_collapsed$is_correct,
                      rank_df_cl_collapsed$cell_line_group, function(v) {
  hit <- which(v)[1]
  if (is.na(hit)) 0 else 1 / hit
}), na.rm = TRUE)

lineage_metrics_cl <- lapply(eval_lineages, function(lin) {
  cl_subset <- top1_cl[cell_lineage == lin]
  if (nrow(cl_subset) == 0) return(NULL)
  data.table(
    lineage = lin,
    n_cell_lines = nrow(cl_subset),
    top1_accuracy = mean(cl_subset$is_correct, na.rm = TRUE)
  )
})
lineage_metrics_cl <- rbindlist(lineage_metrics_cl)

# Top-k lineage consistency (cell-line -> tumour), collapsed to biological units
topk_consistency_cl <- rank_df_cl_collapsed[rank_eval <= opt$`top-k`, .(
  n_total = .N,
  n_same_lineage = sum(tumour_lineage == cell_lineage, na.rm = TRUE),
  frac_same = mean(tumour_lineage == cell_lineage, na.rm = TRUE)
), by = .(cell_line_group, cell_lineage)]
topk_consistency_cl[, all_same := n_same_lineage == n_total]

# Low confidence cases (cell-line -> tumour) on collapsed units
confidence_cl_group <- rank_df_cl_collapsed[, {
  ord <- order(-score)
  s1 <- score[ord][1]
  sk <- score[ord][min(opt$`top-k`, .N)]
  list(s1 = s1, sk = sk, confidence_delta = s1 - sk)
}, by = cell_line_group]
frac_same_cl_group <- rank_df_cl_collapsed[rank_eval <= opt$`top-k`,
  .(confidence_frac_same_lineage = mean(is_correct, na.rm = TRUE)),
  by = cell_line_group]
cellline_summary_group <- merge(confidence_cl_group, frac_same_cl_group,
                                by = "cell_line_group", all.x = TRUE)
cellline_summary_group <- merge(cellline_summary_group, top1_cl[, .(
  cell_line_group,
  cell_lineage,
  best_tumour = tumour,
  best_tumour_lineage = tumour_lineage,
  is_correct
)], by = "cell_line_group", all.x = TRUE)

low_confidence_cl <- cellline_summary_group[confidence_delta < 0.03]


# =============================================================================
# STEP 9: Save Outputs
# =============================================================================
# All results are saved in tab-separated format (TSV) for broad compatibility
# with downstream analysis tools (R, Python, Excel, command-line utilities).
# The full score matrix is saved as RDS for efficient storage and fast
# reloading in R.
#
# Output File Descriptions:
#
#   tumour_to_cellline_rankings.tsv:
#     Complete ranking table enabling detailed analysis of all matches.
#     Useful for: Custom filtering, alternative ranking strategies, QC.
#
#   tumour_mapping_summary.tsv:
#     One row per tumour with key metrics. Suitable for integration with
#     clinical metadata and selection of experimental models.
#
#   metrics_summary.tsv:
#     Aggregate metrics for reporting and comparison between analyses.
#
#   metrics_by_lineage.tsv:
#     Per-lineage breakdown identifying cancer types with good/poor
#     cell line representation.
#
#   tumour_cellline_scores.rds:
#     Raw tumour→cell-line score matrix for advanced analyses.
#
#   cellline_tumour_scores.rds:
#     Transposed cell-line→tumour score matrix for symmetric analyses.

cat("\n[9] Saving outputs...\n")

# Full rankings (all top-k cell lines per tumour)
rank_file <- file.path(t2c_dir, "tumour_to_cellline_rankings.tsv")
fwrite(rank_df, file = rank_file, sep = "\t")
cat("  Rankings saved to:", rank_file, "\n")

# Per-tumour summary
summary_file <- file.path(t2c_dir, "tumour_mapping_summary.tsv")
fwrite(tumour_summary, file = summary_file, sep = "\t")
cat("  Summary saved to:", summary_file, "\n")

# Cell line -> tumour rankings
# Filter NBL to only DSMZ cell lines (starting with NG-) for consistency with figure validation
rank_df_cl_filtered <- rank_df_cl[
  !(cell_lineage == "NBL" & !grepl("^NG-", cell_line))
]
rank_cl_file <- file.path(c2t_dir, "cellline_to_tumour_rankings.tsv")
fwrite(rank_df_cl_filtered, file = rank_cl_file, sep = "\t")
cat("  Cell line rankings saved to:", rank_cl_file, "\n")

# Per-cell-line summary
# Filter NBL to only DSMZ cell lines (starting with NG-) for consistency with figure validation
cellline_summary_filtered <- cellline_summary[
  !(cell_lineage == "NBL" & !grepl("^NG-", cell_line))
]
summary_cl_file <- file.path(c2t_dir, "cellline_mapping_summary.tsv")
fwrite(cellline_summary_filtered, file = summary_cl_file, sep = "\t")
cat("  Cell line summary saved to:", summary_cl_file, "\n")

# Cell line -> tumour metrics (summary + by lineage)
metrics_cl <- data.table(
  metric = c("top1_accuracy", paste0("top", opt$`top-k`, "_accuracy"), "mrr"),
  value = c(top1_acc_cl, topk_acc_cl, mrr_cl)
)
metrics_cl_file <- file.path(c2t_dir, "metrics_summary.tsv")
fwrite(metrics_cl, file = metrics_cl_file, sep = "\t")
cat("  Cell line metrics saved to:", metrics_cl_file, "\n")

if (nrow(lineage_metrics_cl) > 0) {
  lineage_cl_file <- file.path(c2t_dir, "metrics_by_lineage.tsv")
  fwrite(lineage_metrics_cl, file = lineage_cl_file, sep = "\t")
  cat("  Cell line per-lineage metrics saved to:", lineage_cl_file, "\n")
}

# Cell line -> tumour top-k consistency and low-confidence cases
topk_cl_file <- file.path(c2t_dir, "topk_lineage_consistency.tsv")
topk_consistency_cl_out <- copy(topk_consistency_cl)
setnames(topk_consistency_cl_out, "cell_line_group", "cell_line")
fwrite(topk_consistency_cl_out, file = topk_cl_file, sep = "\t")
cat("  Cell line top-k consistency saved to:", topk_cl_file, "\n")

low_conf_cl_file <- file.path(c2t_dir, "low_confidence_cases.tsv")
fwrite(low_confidence_cl[, .(
  cell_line = cell_line_group,
  cell_lineage,
  best_tumour_lineage,
  is_correct,
  confidence_delta,
  confidence_frac_same_lineage
)], file = low_conf_cl_file, sep = "\t")
cat("  Cell line low-confidence cases saved to:", low_conf_cl_file, "\n")

# Overall metrics
metrics <- data.table(
  metric = c("top1_accuracy", paste0("top", opt$`top-k`, "_accuracy"), "mrr"),
  value = c(top1_acc, topk_acc, mrr)
)
metrics_file <- file.path(t2c_dir, "metrics_summary.tsv")
fwrite(metrics, file = metrics_file, sep = "\t")
cat("  Metrics saved to:", metrics_file, "\n")

# Per-lineage metrics
if (nrow(lineage_metrics) > 0) {
  lineage_file <- file.path(t2c_dir, "metrics_by_lineage.tsv")
  fwrite(lineage_metrics, file = lineage_file, sep = "\t")
  cat("  Per-lineage metrics saved to:", lineage_file, "\n")
}

# Full score matrices (RDS for efficient R-native storage)
# Tumour -> cell line
score_file <- file.path(t2c_dir, "tumour_cellline_scores.rds")
saveRDS(score_cor, file = score_file)
cat("  Score matrix saved to:", score_file, "\n")

# Cell line -> tumour (transpose)
score_rev_file <- file.path(c2t_dir, "cellline_tumour_scores.rds")
saveRDS(t(score_cor), file = score_rev_file)
cat("  Score matrix saved to:", score_rev_file, "\n")

# -----------------------------------------------------------------------------
# Move tumour→cell-line analysis docs into direction-specific folder (if present)
# -----------------------------------------------------------------------------
# These files are authored alongside tumour-to-cell-line analyses and should
# live under tumour_to_cellline_similarity for clarity and downstream use.
doc_files <- c(
  "FIGURE_CAPTIONS.md",
  "FINAL_ASSESSMENT.md",
  "INTERPRETATION_GUIDE.md",
  "MAPPING_RESULTS.md",
  "METHODS_AND_RESULTS.md",
  "PUBLICATION_READY.md"
)

for (doc in doc_files) {
  src <- file.path(opt$`output-dir`, doc)
  dest <- file.path(t2c_dir, doc)
  if (file.exists(src)) {
    if (file.exists(dest)) {
      cat("  Skipped doc move (exists):", dest, "\n")
    } else {
      ok <- file.rename(src, dest)
      if (isTRUE(ok)) {
        cat("  Moved doc to:", dest, "\n")
      } else {
        cat("  Warning: failed to move doc:", src, "\n")
      }
    }
  }
}

# Move inspection folder if present (diagnostics for tumour→cell-line mapping)
insp_src <- file.path(opt$`output-dir`, "inspection")
insp_dest <- file.path(t2c_dir, "inspection")
if (dir.exists(insp_src)) {
  if (dir.exists(insp_dest)) {
    cat("  Skipped inspection move (exists):", insp_dest, "\n")
  } else {
    ok <- file.rename(insp_src, insp_dest)
    if (isTRUE(ok)) {
      cat("  Moved inspection to:", insp_dest, "\n")
    } else {
      cat("  Warning: failed to move inspection:", insp_src, "\n")
    }
  }
}


# =============================================================================
# Execution Summary
# =============================================================================
# Final summary output provides a quick overview of the analysis results,
# enabling rapid assessment of mapping quality without examining output files.

cat("\n========================================\n")
cat("MAPPING COMPLETE\n")
cat("========================================\n")
cat("\nSummary:\n")
cat("  Tumours mapped: ", nrow(tumour_summary), "\n", sep = "")
cat("  Top-1 accuracy: ", round(top1_acc * 100, 2), "%\n", sep = "")
cat("  Top-", opt$`top-k`, " accuracy: ", round(topk_acc * 100, 2), "%\n", sep = "")
cat("  MRR: ", round(mrr, 3), "\n", sep = "")


# =============================================================================
# END OF SCRIPT
# =============================================================================
# 
# METHODOLOGICAL NOTES FOR THESIS DOCUMENTATION
# ---------------------------------------------
# This script implements a correlation-based tumour-cell line matching approach
# that addresses several key challenges in cancer model selection:
#
# 1. Batch Effect Mitigation:
#    By using jointly normalised expression matrices (tumours + cell lines
#    from the same cancer type), systematic technical differences between
#    datasets are minimised, ensuring similarity scores reflect biology
#    rather than batch effects.
#
# 2. Feature Selection:
#    Operating within a curated pan-cancer feature space focuses the analysis
#    on genes that capture biologically meaningful variation, improving
#    signal-to-noise ratio and enabling cross-lineage comparisons.
#
# 3. Non-Parametric Similarity:
#    Spearman rank correlation provides a robust, distribution-free measure
#    of transcriptomic similarity that is appropriate for heterogeneous
#    expression data.
#
# 4. Comprehensive Evaluation:
#    Multiple complementary metrics (top-1, top-k, MRR) provide nuanced
#    assessment of mapping quality, while confidence metrics enable
#    identification of tumours with clear vs. ambiguous matches.
#
# 5. Dual Perspective Analysis:
#    Both tumour-driven and cell-line-driven metrics are computed, enabling
#    both personalised model selection and general cell line characterisation.
#
# Limitations and Future Directions:
#
#   - Lineage-based evaluation is a coarse approximation; molecular subtype
#     annotation would enable finer-grained assessment.
#
#   - The analysis uses transcriptomic similarity alone; integration of
#     genomic (mutation, copy number) or epigenomic data could improve
#     matching accuracy.
#
#   - Correlation-based approaches may miss important differences in
#     specific gene sets; pathway-aware or signature-based methods could
#     complement this analysis.
#
# =============================================================================
