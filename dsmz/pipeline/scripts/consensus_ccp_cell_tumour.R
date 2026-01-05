#!/usr/bin/env Rscript

# ==============================================================================
# consensus_ccp_cell_tumour.R
# Consensus Clustering Pipeline for Tumour and Cell Line Integration
# ==============================================================================
#
# PURPOSE:
# This script performs consensus clustering using the ConsensusClusterPlus 
# algorithm to identify robust sample groupings across cancer 
# tumour samples and DSMZ cancer cell lines. The pipeline supports 
# multiple feature sets, distance metrics, and clustering algorithms.
#
# USAGE:
# Rscript consensus_ccp_cell_tumour.R \
#   --config config/config.yaml \
#   --direction <direction_name> \
#   --kind <cell|tumour|cell_tumour> \
#   --mode <expr|pca> \
#   --alg <km|hc> \
#   --outdir <output_directory> \
#   --cluster_rds <output_rds_path>
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The suppressPackageStartupMessages() function prevents verbose startup 
# messages from cluttering the console output, which is particularly useful 
# when running the script in automated pipelines.

suppressPackageStartupMessages({
  library(yaml)                 # YAML configuration file parsing
  library(optparse)             # Command-line argument parsing
  library(dplyr)                # Data manipulation utilities
  library(readr)                # Efficient CSV reading and writing
  library(tibble)               # Modern data frame implementation
  library(ConsensusClusterPlus) # Core consensus clustering algorithm
  library(matrixStats)          # Efficient row/column statistics (e.g., rowVars)
})

# ------------------------------------------------------------------------------
# SECTION 2: UTILITY FUNCTIONS
# ------------------------------------------------------------------------------
# Two helper functions are defined for standardised logging and error handling.
# These ensure consistent output formatting throughout the pipeline.

# info(): Prints informational messages with a consistent prefix.
# The sprintf() function enables formatted string construction with variable
# substitution, similar to C-style printf formatting.
info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")

# stop_with(): Terminates execution with an error message.
# Setting call. = FALSE suppresses the call stack in the error output,
# producing cleaner error messages for end users.
stop_with <- function(...) stop(paste0(...), call. = FALSE)

# ------------------------------------------------------------------------------
# SECTION 3: COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
# The optparse package provides a systematic approach to defining and parsing
# command-line arguments. Each make_option() call defines a single argument
# with its type and default value.

option_list <- list(
  make_option("--config",      type = "character", default = "config/config.yaml"),
  make_option("--direction",   type = "character", default = NULL),
  make_option("--kind",        type = "character", default = NULL),
  make_option("--mode",        type = "character", default = NULL),  # expr | pca
  make_option("--alg",         type = "character", default = NULL),  # km | hc
  make_option("--outdir",      type = "character", default = NULL),
  make_option("--cluster_rds", type = "character", default = NULL)
)

# parse_args() processes the command-line arguments according to the defined
# option_list and returns a named list containing the parsed values.
opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------------------------
# SECTION 4: INPUT VALIDATION
# ------------------------------------------------------------------------------
# Rigorous input validation prevents cryptic downstream errors by catching
# invalid or missing parameters early in execution. Each required argument
# is checked for NULL values, and enumerated arguments are validated against
# their allowed values.

if (is.null(opt$direction))   stop_with("Missing --direction")
if (is.null(opt$kind))        stop_with("Missing --kind")
if (is.null(opt$mode))        stop_with("Missing --mode (expr|pca)")
if (is.null(opt$alg))         stop_with("Missing --alg (km|hc)")
if (is.null(opt$outdir))      stop_with("Missing --outdir")
if (is.null(opt$cluster_rds)) stop_with("Missing --cluster_rds")

# Validate enumerated parameters against allowed values.
# The %in% operator checks membership in a character vector.
if (!opt$mode %in% c("expr","pca")) stop_with("--mode must be expr or pca")
if (!opt$alg  %in% c("km","hc"))    stop_with("--alg must be km or hc")

# Create the output directory if it does not exist.
# Setting recursive = TRUE allows creation of nested directory structures.
# Setting showWarnings = FALSE suppresses warnings if the directory exists.
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# SECTION 5: CONFIGURATION LOADING
# ------------------------------------------------------------------------------
# The YAML configuration file centralises parameters that may change between
# runs, separating configuration from code logic.

cfg <- yaml::read_yaml(opt$config)

# Validate that the specified direction exists in the configuration.
# This prevents cryptic NULL-related errors during subsequent operations.
if (is.null(cfg$methods[[opt$direction]])) {
  stop_with("Unknown direction in config methods: ", opt$direction)
}

# Extract feature type and distance metric from the configuration.
# These determine which gene set is used and how sample similarity is computed.
feature  <- cfg$methods[[opt$direction]]$feature   # PAM50 | HVG
distance <- cfg$methods[[opt$direction]]$distance  # euclidean | correlation

# ------------------------------------------------------------------------------
# SECTION 6: CLUSTERING PARAMETER EXTRACTION
# ------------------------------------------------------------------------------
# The k_grid defines the range of cluster numbers to evaluate.
# ConsensusClusterPlus requires at least two k values to compute delta AUC.

k_grid <- cfg$clustering$k_grid
if (is.null(k_grid) || length(k_grid) < 2) stop_with("config: clustering.k_grid must have >=2 values")

# Convert to sorted integer vector and extract maximum k.
# The unlist() function flattens any nested list structure from YAML parsing.
k_grid <- sort(as.integer(unlist(k_grid)))
maxK   <- max(k_grid)

# Set random seed for reproducibility.
# This ensures that bootstrap sampling produces identical results across runs.
seed <- 42
if (!is.null(cfg$agnostic_clustering$seed)) seed <- as.integer(cfg$agnostic_clustering$seed)

# Set number of principal components for PCA mode.
# This parameter controls dimensionality reduction when mode = "pca".
n_pcs <- 20
if (!is.null(cfg$agnostic_clustering$n_pcs)) n_pcs <- as.integer(cfg$agnostic_clustering$n_pcs)

# ------------------------------------------------------------------------------
# SECTION 7: DISTANCE METRIC CONSIDERATIONS
# ------------------------------------------------------------------------------
# IMPORTANT METHODOLOGICAL NOTE:
# ConsensusClusterPlus is configured to use Euclidean distance with Ward.D2
# linkage for all analyses. When correlation-based clustering is desired,
# the input data undergoes a geometric transformation (see Section 11) that
# converts Euclidean distance calculations to correlation-equivalent geometry.
# This approach ensures compatibility with Ward's minimum variance criterion,
# which requires Euclidean distance for mathematical validity.
# ------------------------------------------------------------------------------

ccp_distance <- "euclidean"  # Always Euclidean for CCP internal calculations
innerLinkage <- "ward.D2"    # Linkage for consensus matrix construction
finalLinkage <- "ward.D2"    # Linkage for final cluster assignment

# Log the analysis configuration for reproducibility documentation.
info("direction=%s | feature=%s | distance=%s | mode=%s | alg=%s | kind=%s",
     opt$direction, feature, distance, opt$mode, opt$alg, opt$kind)

# ------------------------------------------------------------------------------
# SECTION 8: EXPRESSION MATRIX CONSTRUCTION FUNCTION
# ------------------------------------------------------------------------------
# The build_expr_mat() function handles the complexity of loading and 
# integrating expression data from multiple sources. It returns a list
# containing the expression matrix and dataset origin annotations.
#
# PARAMETERS:
#   feature: The feature set to use ("PAM50" or "HVG")
#   cfg:     The parsed configuration object
#   kind:    The sample scope ("cell", "tumour", or "cell_tumour")
#
# RETURNS:
#   A list with two elements:
#     - expr_mat: A genes × samples matrix of expression values
#     - dataset:  A named vector indicating sample origin (TCGA or DSMZ)

build_expr_mat <- function(feature, cfg, kind) {
  
  # Determine the sample scope from the kind argument using pattern matching.
  # The grepl() function performs regular expression matching.
  is_cell_tumour <- grepl("_cell_tumour", kind)
  is_cell_only   <- grepl("_cell$", kind) && !is_cell_tumour
  is_tumour_only <- grepl("_tumour$", kind) && !is_cell_tumour

  # --------------------------------------------------------------------------
  # PAM50 FEATURE SET HANDLING
  # --------------------------------------------------------------------------
  # PAM50 is a clinically validated 50-gene signature for breast cancer

  # subtype classification. Pre-computed expression matrices for these genes
  # are loaded from RDS files specified in the configuration.
  
  if (feature == "PAM50") {
    tcga_path <- cfg$paths$tcga_brca_pam50_expr
    dsmz_path <- cfg$paths$dsmz_bcc_pam50_expr
    
    if (is.null(tcga_path) || is.null(dsmz_path)) {
      stop_with("Missing pam50 paths: paths.tcga_brca_pam50_expr or paths.dsmz_bcc_pam50_expr")
    }
    
    # Load pre-computed PAM50 expression matrices.
    # These are expected to be in genes × samples orientation.
    tcga_mat <- readRDS(tcga_path)
    dsmz_mat <- readRDS(dsmz_path)

    # Identify genes present in both matrices.
    # The intersect() function returns the common elements of two vectors.
    common_genes <- intersect(rownames(tcga_mat), rownames(dsmz_mat))
    
    if (length(common_genes) < 20) {
      stop_with("Too few common genes between TCGA and DSMZ PAM50 matrices.")
    }

    # Construct the expression matrix and dataset annotations based on scope.
    if (is_cell_tumour) {
      # Combine both datasets using cbind() for column-wise binding.
      expr_mat <- cbind(tcga_mat[common_genes, , drop = FALSE],
                        dsmz_mat[common_genes, , drop = FALSE])
      
      # Create dataset origin vector using conditional assignment.
      # The ifelse() function performs element-wise conditional evaluation.
      dataset <- ifelse(colnames(expr_mat) %in% colnames(dsmz_mat), "DSMZ", "TCGA")
      names(dataset) <- colnames(expr_mat)
      
    } else if (is_cell_only) {
      expr_mat <- dsmz_mat[common_genes, , drop = FALSE]
      dataset <- setNames(rep("DSMZ", ncol(expr_mat)), colnames(expr_mat))
      
    } else if (is_tumour_only) {
      expr_mat <- tcga_mat[common_genes, , drop = FALSE]
      dataset <- setNames(rep("TCGA", ncol(expr_mat)), colnames(expr_mat))
      
    } else {
      stop_with("Cannot determine scope from kind: ", kind)
    }

    return(list(expr_mat = expr_mat, dataset = dataset))
  }

  # --------------------------------------------------------------------------
  # HVG (HIGHLY VARIABLE GENES) FEATURE SET HANDLING
  # --------------------------------------------------------------------------
  # Highly variable genes are identified through variance-based feature
  # selection, capturing genes with the greatest expression variability
  # across samples. This data-driven approach complements the biologically-
  # informed PAM50 signature.
  
  if (feature == "HVG") {
    cell_path   <- cfg$paths$cell_vst_rds
    tumour_path <- cfg$paths$tumour_vst_rds
    hvg_list    <- cfg$features$hvg_final_gene_list

    if (is.null(cell_path) || is.null(tumour_path)) {
      stop_with("Missing: paths.cell_vst_rds and/or paths.tumour_vst_rds")
    }
    
    if (is.null(hvg_list) || !file.exists(hvg_list)) {
      stop_with("Missing HVG list: features.hvg_final_gene_list")
    }

    # Load variance-stabilised transformed expression matrices.
    # VST transformation normalises variance across the expression range,
    # making genes with different expression levels more comparable.
    cell_mat   <- readRDS(cell_path)
    tumour_mat <- readRDS(tumour_path)

    # Load the HVG list, filtering out empty lines.
    # The nzchar() function tests for non-zero-length character strings.
    genes <- readLines(hvg_list)
    genes <- genes[nzchar(genes)]

    if (is_cell_tumour) {
      # Retain only HVGs present in both cell line and tumour matrices.
      keep <- intersect(genes, intersect(rownames(cell_mat), rownames(tumour_mat)))
      
      if (length(keep) < 50) {
        stop_with("Too few HVG genes found in both cell and tumour matrices.")
      }
      
      expr_mat <- cbind(tumour_mat[keep, , drop = FALSE],
                        cell_mat[keep, , drop = FALSE])
      
      # Construct dataset origin vector for the combined matrix.
      dataset <- c(
        setNames(rep("TCGA", ncol(tumour_mat[keep, , drop = FALSE])), colnames(tumour_mat)),
        setNames(rep("DSMZ", ncol(cell_mat[keep, , drop = FALSE])), colnames(cell_mat))
      )
      dataset <- dataset[colnames(expr_mat)]
      names(dataset) <- colnames(expr_mat)
      
    } else if (is_cell_only) {
      keep <- intersect(genes, rownames(cell_mat))
      
      if (length(keep) < 50) {
        stop_with("Too few HVG genes found in cell matrix.")
      }
      
      expr_mat <- cell_mat[keep, , drop = FALSE]
      dataset <- setNames(rep("DSMZ", ncol(expr_mat)), colnames(expr_mat))
      
    } else if (is_tumour_only) {
      keep <- intersect(genes, rownames(tumour_mat))
      
      if (length(keep) < 50) {
        stop_with("Too few HVG genes found in tumour matrix.")
      }
      
      expr_mat <- tumour_mat[keep, , drop = FALSE]
      dataset <- setNames(rep("TCGA", ncol(expr_mat)), colnames(expr_mat))
      
    } else {
      stop_with("Cannot determine scope from kind: ", kind)
    }

    return(list(expr_mat = expr_mat, dataset = dataset))
  }

  stop_with("Unsupported feature: ", feature)
}

# ------------------------------------------------------------------------------
# SECTION 9: CORRELATION-GEOMETRY TRANSFORMATION FUNCTION
# ------------------------------------------------------------------------------
# This function transforms the input matrix such that Euclidean distance
# between columns (samples) becomes equivalent to correlation-based distance.
#
# MATHEMATICAL BASIS:
# For two vectors x and y that are centred (mean = 0) and L2-normalised
# (||x|| = ||y|| = 1), the Euclidean distance relates to Pearson correlation:
#   ||x - y||² = 2(1 - cor(x, y))
#
# This transformation enables the use of Ward's minimum variance criterion,
# which requires Euclidean distance, while preserving correlation-based
# sample relationships.
#
# PARAMETERS:
#   X: A features × samples matrix
#
# RETURNS:
#   The transformed matrix with centred, unit-normalised columns

corr_geometry_transform <- function(X) {
  # Centre each column (sample) by subtracting its mean.
  # The sweep() function applies an operation across matrix margins.
  # MARGIN = 2 indicates column-wise operation.
  X <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  
  # Compute the L2 norm (Euclidean length) of each column.
  norms <- sqrt(colSums(X^2, na.rm = TRUE))
  
  # Handle edge cases where the norm is zero or non-finite.
  # This prevents division by zero for constant columns.
  norms[!is.finite(norms) | norms == 0] <- 1
  
  # Normalise each column to unit length.
  X <- sweep(X, 2, norms, "/")
  
  X
}

# ------------------------------------------------------------------------------
# SECTION 10: EXPRESSION MATRIX CONSTRUCTION AND PREPROCESSING
# ------------------------------------------------------------------------------
# The expression matrix is constructed using the build_expr_mat() function,
# which handles the complexity of loading and integrating data from multiple
# sources based on the specified feature set and sample scope.

built    <- build_expr_mat(feature, cfg, opt$kind)
expr_mat <- built$expr_mat
dataset  <- built$dataset

# Remove genes with zero variance (constant expression).
# Such genes provide no discriminatory information for clustering and can
# cause numerical issues in downstream calculations.
# The rowVars() function from matrixStats efficiently computes row variances.
v <- rowVars(expr_mat, na.rm = TRUE)
expr_mat <- expr_mat[v > 0 & is.finite(v), , drop = FALSE]

# Verify sample identifier uniqueness.
# Duplicate identifiers would corrupt clustering results and downstream
# sample-to-cluster mappings.
if (anyDuplicated(colnames(expr_mat))) {
  stop_with("Duplicate sample IDs detected in integrated matrix.")
}

# ------------------------------------------------------------------------------
# SECTION 11: INPUT DATA PREPARATION FOR CONSENSUS CLUSTERING
# ------------------------------------------------------------------------------
# ConsensusClusterPlus expects a matrix with features as rows and samples
# (items to be clustered) as columns. Two input modes are supported:
#
# 1. Expression mode (expr): Uses gene expression values directly as features.
#    This preserves the original high-dimensional gene space.
#
# 2. PCA mode (pca): Reduces dimensionality by projecting samples onto the
#    top principal components. This can reduce noise and computational cost.

if (opt$mode == "expr") {
  # Expression mode: use the gene expression matrix directly.
  # Matrix orientation is already genes × samples.
  data_for_ccp <- expr_mat
  info("CCP input = expr (genes × samples): %d × %d", nrow(data_for_ccp), ncol(data_for_ccp))
  
} else {
  # PCA mode: perform principal component analysis on samples.
  # First transpose to samples × genes orientation for prcomp().
  samples_mat <- t(expr_mat)
  
  # Perform PCA with centring but without scaling.
  # Centring removes the mean of each gene across samples.
  # Scaling is omitted to preserve relative variance differences between genes.
  pca <- prcomp(samples_mat, center = TRUE, scale. = FALSE)
  
  # Extract the top n_pcs principal component scores.
  # The $x component contains the sample scores (projections).
  pcs_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]  # samples × PCs
  
  # Transpose to PCs × samples orientation for ConsensusClusterPlus.
  data_for_ccp <- t(pcs_mat)
  info("CCP input = pca (PCs × samples): %d × %d", nrow(data_for_ccp), ncol(data_for_ccp))
}

# Apply correlation-geometry transformation if correlation distance is specified.
# This transforms the data such that Euclidean distance becomes equivalent to
# correlation-based distance, enabling the use of Ward.D2 linkage.
if (distance == "correlation") {
  info("Applying correlation-geometry transform (center + unit norm) so Ward.D2 uses euclidean safely.")
  data_for_ccp <- corr_geometry_transform(data_for_ccp)
}

# ------------------------------------------------------------------------------
# SECTION 12: CONSENSUS CLUSTERING EXECUTION
# ------------------------------------------------------------------------------
# ConsensusClusterPlus outputs multiple diagnostic files (plots, CSVs) to
# the working directory. A dedicated subdirectory is created to isolate
# these outputs, and the working directory is temporarily changed.

ccp_wd <- file.path(opt$outdir, "ccp")
dir.create(ccp_wd, showWarnings = FALSE, recursive = TRUE)

# Save the current working directory and switch to the CCP output directory.
# The on.exit() call ensures the original directory is restored even if
# an error occurs during execution.
old_wd <- getwd()
setwd(ccp_wd)
on.exit(setwd(old_wd), add = TRUE)

# Construct the argument list for ConsensusClusterPlus.
# Key parameters:
#   d:            Input data matrix (features × samples)
#   maxK:         Maximum number of clusters to evaluate
#   reps:         Number of bootstrap resampling iterations
#   pItem:        Proportion of samples to include in each bootstrap
#   pFeature:     Proportion of features to include (1.0 = all features)
#   clusterAlg:   Clustering algorithm ("km" = k-means, "hc" = hierarchical)
#   distance:     Distance metric for consensus matrix calculation
#   innerLinkage: Linkage method for consensus matrix construction
#   finalLinkage: Linkage method for final cluster assignment
#   seed:         Random seed for reproducibility
#   title:        Prefix for output file names
#   plot:         Output format for diagnostic plots
#   writeTable:   Whether to write cluster assignments to CSV files
ccp_args <- list(
  d            = data_for_ccp,
  maxK         = maxK,
  reps         = 1000,
  pItem        = 0.8,
  pFeature     = 1.0,
  clusterAlg   = opt$alg,
  distance     = ccp_distance,   # Always Euclidean (see Section 7)
  innerLinkage = innerLinkage,   # Ward.D2 for robust clustering
  finalLinkage = finalLinkage,   # Ward.D2 for final assignment
  seed         = seed,
  title        = paste0("CCP_", opt$direction, "_", opt$kind),
  plot         = "png",
  writeTable   = TRUE,
  verbose      = FALSE
)

# Execute ConsensusClusterPlus using do.call() for argument unpacking.
# The result is a list indexed by k, where each element contains:
#   - consensusMatrix: The n × n sample co-clustering frequency matrix
#   - consensusClass:  The cluster assignment vector
#   - consensusTree:   The hierarchical clustering dendrogram
ccp_res <- do.call(ConsensusClusterPlus, ccp_args)

# ------------------------------------------------------------------------------
# SECTION 13: OPTIMAL CLUSTER NUMBER SELECTION
# ------------------------------------------------------------------------------
# The optimal number of clusters is determined using the area under the
# cumulative distribution function (CDF) of consensus values. This metric
# quantifies clustering stability: higher AUC indicates more consistent
# sample co-clustering across bootstrap iterations.
#
# The delta AUC (ΔAUC) measures the improvement in clustering stability
# when increasing from k-1 to k clusters. The k with maximum ΔAUC
# represents the point of greatest stability improvement.

# Define a function to compute AUC for a given k.
# The function extracts consensus values from the upper triangle of the
# consensus matrix and integrates the empirical CDF using the trapezoidal rule.
auc_at_k <- function(k) {
  # Extract the consensus matrix for the specified k.
  cm <- ccp_res[[k]]$consensusMatrix
  
  # Extract the upper triangle values (excluding the diagonal).
  # The consensus matrix is symmetric, so only the upper triangle is needed.
  vals <- cm[upper.tri(cm)]
  
  # Compute the empirical cumulative distribution function.
  ec <- ecdf(vals)
  
  # Evaluate the CDF at 500 evenly-spaced points and integrate.
  # Trapezoidal integration: average of consecutive points times interval width.
  x <- seq(0, 1, length.out = 500)
  sum(diff(x) * (ec(head(x, -1)) + ec(tail(x, -1))) / 2)
}

# Compute AUC values across the k grid.
# vapply() is used for type-safe iteration with a numeric return value.
auc_values <- vapply(k_grid, auc_at_k, numeric(1))

# Compute delta AUC (change in AUC between successive k values).
delta_auc <- diff(auc_values)

# Select the k corresponding to maximum delta AUC.
# The +1 offset accounts for diff() reducing the vector length by one.
best_k <- k_grid[which.max(delta_auc) + 1]

# Log the selection metrics for transparency and debugging.
info("k_grid: %s", paste(k_grid, collapse = ","))
info("AUC:    %s", paste(round(auc_values, 4), collapse = " "))
info("ΔAUC:   %s", paste(round(delta_auc, 4), collapse = " "))
info("best_k: %d", best_k)

# ------------------------------------------------------------------------------
# SECTION 14: FINAL CLUSTER ASSIGNMENT EXTRACTION
# ------------------------------------------------------------------------------
# Extract the cluster assignments from the optimal k solution and format
# them with a descriptive prefix for clarity in downstream analyses.

final_clusters <- ccp_res[[best_k]]$consensusClass

# Add "C" prefix to cluster numbers for human readability.
# The paste0() function concatenates without a separator.
final_clusters <- setNames(paste0("C", final_clusters), names(final_clusters))

# ------------------------------------------------------------------------------
# SECTION 15: OUTPUT FILE GENERATION
# ------------------------------------------------------------------------------
# Two output files are generated:
# 1. A CSV file containing sample-level cluster assignments for easy inspection
# 2. An RDS file containing the complete analysis results for programmatic access

# Write the cluster assignments to a CSV file.
# The tibble structure provides clean formatting with sample ID, cluster,
# and dataset origin columns.
clusters_csv <- file.path(opt$outdir, sprintf("tcga_dsmz_ccp_clusters_K%d.csv", best_k))

write_csv(
  tibble(
    sample_id = names(final_clusters),
    cluster   = unname(final_clusters),
    dataset   = dataset[names(final_clusters)]
  ),
  clusters_csv
)

# Save the comprehensive results object as an RDS file.
# This file contains all information needed to reproduce or extend the analysis:
#   - Input parameters (direction, kind, feature, distance, mode, alg)
#   - Clustering configuration (ccp_distance, linkage methods, k_grid)
#   - Selection metrics (auc_values, delta_auc, best_k)
#   - Results (clusters, dataset annotations, full CCP output)
#   - Metadata (timestamp for provenance tracking)
saveRDS(
  list(
    direction      = opt$direction,
    kind           = opt$kind,
    feature        = feature,
    distance       = distance,
    mode           = opt$mode,
    alg            = opt$alg,
    ccp_distance   = ccp_distance,
    innerLinkage   = innerLinkage,
    finalLinkage   = finalLinkage,
    k_grid         = k_grid,
    auc_values     = setNames(auc_values, k_grid),
    delta_auc      = delta_auc,
    best_k         = best_k,
    clusters       = final_clusters,
    dataset        = dataset,
    ccp_results    = ccp_res,
    timestamp      = Sys.time()
  ),
  opt$cluster_rds
)

# Log the output file paths for user reference.
info("Wrote: %s", clusters_csv)
info("Wrote: %s", opt$cluster_rds)
info("Done.")