#!/usr/bin/env Rscript

# ==============================================================================
# CONSENSUS CLUSTERING FOR TUMOUR AND CELL LINE INTEGRATION
# ==============================================================================
#
# SCRIPT OVERVIEW
# ---------------
# This script implements consensus clustering using the ConsensusClusterPlus
# algorithm to identify robust molecular groupings within and across cancer
# tumour samples and DSMZ (Deutsche Sammlung von Mikroorganismen und
# Zellkulturen) cancer cell lines. The methodology enables systematic
# evaluation of cluster stability through bootstrap resampling, providing
# confidence measures for the identified sample groupings.
#
# BIOLOGICAL CONTEXT
# ------------------
# Cancer is characterised by substantial molecular heterogeneity, with tumours
# of the same histological type often exhibiting distinct transcriptomic
# profiles that correlate with clinical outcomes and treatment responses.
# Identifying these molecular subtypes through unsupervised clustering is
# fundamental to precision oncology.
#
# Cell lines derived from primary tumours serve as essential models for
# cancer research and drug development. However, the fidelity with which
# these cell lines represent the molecular diversity of primary tumours
# remains an open question. By clustering cell lines together with primary
# tumour samples, this analysis assesses whether cell lines integrate into
# tumour-defined molecular subtypes or form distinct, potentially divergent,
# groupings.
#
# CONSENSUS CLUSTERING RATIONALE
# ------------------------------
# Standard clustering algorithms produce a single partition of samples, but
# this partition may be sensitive to:
#   - Random initialisation (e.g., k-means centroid selection)
#   - Specific samples included in the analysis
#   - Choice of distance metric and linkage method
#
# Consensus clustering addresses these concerns by repeatedly clustering
# random subsamples of the data and recording how frequently each pair of
# samples is assigned to the same cluster. The resulting consensus matrix
# provides a robust measure of sample relationships that is less sensitive
# to individual algorithmic choices.
#
# The ConsensusClusterPlus algorithm (Wilkerson & Hayes, 2010) implements
# this approach with the following procedure:
#   1. For each of B bootstrap iterations (typically B=1000):
#      a. Randomly subsample a fraction (typically 80%) of samples
#      b. Perform clustering on the subsample
#      c. Record which sample pairs were assigned to the same cluster
#   2. Compute the consensus matrix: M[i,j] = (co-clustering count) / (co-sampled count)
#   3. Cluster the consensus matrix to obtain final assignments
#   4. Evaluate cluster stability using the cumulative distribution function (CDF)
#
# FEATURE SETS
# ------------
# The script supports multiple feature sets for clustering:
#
#   PAM50: The Prediction Analysis of Microarray 50 (PAM50) gene signature
#   is a clinically validated classifier for breast cancer molecular subtypes.
#   Developed by Parker et al. (2009), PAM50 identifies five intrinsic subtypes:
#   Luminal A, Luminal B, HER2-enriched, Basal-like, and Normal-like. These
#   subtypes have distinct prognoses and treatment implications.
#
#   HVG (Highly Variable Genes): A data-driven feature set identified through
#   variance-based selection. Genes with high expression variability across
#   samples are more likely to capture biologically meaningful differences
#   between sample groups.
#
#   Method-specific feature sets (Variance, MAD, Entropy, PCA, Spearman, etc.):
#   These represent alternative feature selection strategies, each capturing
#   different aspects of expression variability. Evaluating multiple feature
#   sets enables assessment of result robustness to methodological choices.
#
# DISTANCE METRICS
# ----------------
# Two distance metrics are supported:
#
#   Euclidean distance: The straight-line distance in gene expression space.
#   Appropriate for variance-stabilised transformed (VST) data where the
#   variance is approximately constant across the expression range.
#
#   Correlation-based distance: Captures similarity in relative expression
#   patterns, making samples with proportional expression profiles appear
#   similar regardless of absolute expression levels. Implemented through
#   a geometric transformation (see Section 9) that converts Euclidean
#   calculations to correlation-equivalent geometry.
#
# SAMPLE VIEWS
# ------------
# Three clustering configurations are supported:
#
#   cell: Clusters cell lines independently to reveal intrinsic groupings
#   based on transcriptomic similarity. This may reflect tissue of origin,
#   driver mutation status, or culture-induced adaptations.
#
#   tumour: Clusters primary tumour samples independently to identify
#   molecular subtypes. These groupings often correspond to clinically
#   relevant categories with distinct prognoses.
#
#   cell_tumour: Clusters cell lines and tumours jointly, enabling direct
#   assessment of whether cell lines integrate into tumour-defined subtypes.
#   This is the primary configuration for evaluating cell line representativeness.
#
# DIMENSIONALITY REDUCTION
# ------------------------
# Two input modes are supported:
#
#   expr: Uses the full gene expression matrix as clustering features.
#   This preserves all gene-level information but may be susceptible to
#   noise from uninformative genes.
#
#   pca: Projects samples onto the top principal components before clustering.
#   This reduces noise by discarding low-variance components and improves
#   computational efficiency, but may lose gene-level interpretability.
#
# OUTPUT
# ------
# The script produces:
#   - Cluster assignments for each sample at the optimal k
#   - Consensus matrices visualising sample co-clustering frequencies
#   - CDF plots for cluster stability assessment
#   - A comprehensive RDS object containing all analysis parameters and results
# ==============================================================================
# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The following R packages provide essential functionality for this analysis.
# The suppressPackageStartupMessages() wrapper prevents verbose startup
# messages from cluttering log files during automated pipeline execution.

suppressPackageStartupMessages({
  
  # yaml: Provides functions for reading YAML configuration files.
  # YAML (YAML Ain't Markup Language) is a human-readable data serialisation
  # format commonly used for configuration files in bioinformatics pipelines.
  library(yaml)
  
  # optparse: Enables systematic command-line argument parsing.
  # This package follows the GNU-style option conventions, supporting both
  # short (-x) and long (--option) argument formats with type validation.
  library(optparse)
  
  # dplyr: Provides a grammar of data manipulation with functions for
  # filtering, selecting, mutating, and summarising data frames.
  # Part of the tidyverse ecosystem for data science in R.
  library(dplyr)
  
  # readr: Implements fast and consistent functions for reading and writing
  # rectangular data (CSV, TSV). Provides better defaults than base R
  # functions for handling character encodings and column type inference.
  library(readr)
  
  # tibble: A modern reimplementation of data frames with stricter subsetting
  # rules and better printing methods. Tibbles never convert strings to
  # factors and never create row names automatically.
  library(tibble)
  
  # ConsensusClusterPlus: The core algorithm implementation for consensus
  # clustering. This Bioconductor package provides robust cluster discovery
  # with confidence assessment through bootstrap resampling.
  library(ConsensusClusterPlus)
  
  # matrixStats: Provides highly optimised functions for row and column
  # operations on matrices. The rowVars() function used in this script
  # is substantially faster than apply(mat, 1, var) for large matrices.
  library(matrixStats)
})

# ------------------------------------------------------------------------------
# SECTION 2: UTILITY FUNCTIONS
# ------------------------------------------------------------------------------
# Two helper functions standardise logging and error handling throughout
# the script. Consistent formatting improves log readability and facilitates
# automated parsing of pipeline outputs.

info <- function(...) {
  # info(): Outputs informational messages with a standardised prefix.
  #
  # The sprintf() function enables C-style formatted string construction,

  # supporting format specifiers such as %s (string), %d (integer), and
  # %f (floating point). This is more readable than paste() for complex
  # messages with multiple variable substitutions.
  #

  # Parameters:
  #   ...: Format string followed by values to substitute
  #
  # Example usage:
  #   info("Processing %d samples from %s dataset", n_samples, dataset_name)
  
  cat("[INFO] ", sprintf(...), "\n", sep = "")
}

stop_with <- function(...) {
  # stop_with(): Terminates execution with a formatted error message.
  #
  # The call. = FALSE argument suppresses the function call stack in the
  # error output. This produces cleaner error messages for end users who
  # do not need to see internal R function calls.
  #
  # Parameters:
  #   ...: Components of the error message (concatenated with paste0)
  #
  # Example usage:
  #   stop_with("Missing required file: ", filepath)
  
  stop(paste0(...), call. = FALSE)
}

# ------------------------------------------------------------------------------
# SECTION 3: COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
# The optparse package provides a systematic framework for defining and
# parsing command-line arguments. Each argument is specified with its
# name, type, default value, and optional help text.
#
# This approach enables the script to be called from automated pipelines
# with different parameter combinations without modifying the source code.

option_list <- list(
  
  # --config: Path to the YAML configuration file containing analysis parameters.
  # The configuration file centralises settings that may vary between runs,
  # including file paths, algorithm parameters, and feature set definitions.
  make_option("--config", type = "character", default = "config/config.yaml"),
  

  # --profile: The configuration profile to use (e.g., 'brca', 'nbl', 'rbl').
  # Profiles enable disease-specific parameter sets within a single config file.
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  
  # --direction: Specifies the feature set and distance metric combination.
  # Format: {feature_set}_{distance_metric} (e.g., 'HVG_euc', 'MAD_corr')
  make_option("--direction", type = "character", default = NULL),

  # --feature_list: Optional explicit gene list file for the requested feature.
  # Required for HVG directions so the script uses the actual HVG-derived set
  # rather than falling back to legacy aliases.
  make_option("--feature_list", type = "character", default = NULL,
              help = "Path to the feature gene list (required for HVG directions)"),
  
  # --kind: Specifies the sample scope for clustering.
  # Valid values: 'cell' (cell lines only), 'tumour' (tumours only),
  # 'cell_tumour' (combined analysis)
  make_option("--kind", type = "character", default = NULL),
  
  # --mode: Specifies the input data representation.
  # 'expr': Use gene expression values directly
  # 'pca': Use principal component scores (dimensionality reduction)
  make_option("--mode", type = "character", default = NULL),
  
  # --alg: Specifies the clustering algorithm.
  # 'km': k-means clustering (partitional, Euclidean distance only)
  # 'hc': Hierarchical clustering (agglomerative, supports multiple distances)
  make_option("--alg", type = "character", default = NULL),
  
  # --outdir: Directory for output files (created if it does not exist).
  make_option("--outdir", type = "character", default = NULL),
  
  # --cluster_rds: Path for the output RDS file containing all results.
  make_option("--cluster_rds", type = "character", default = NULL)
)

# Parse command-line arguments according to the defined option specifications.
# The parse_args() function returns a named list where each element
# corresponds to an option with its parsed (or default) value.
opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------------------------
# SECTION 4: INPUT VALIDATION
# ------------------------------------------------------------------------------
# Rigorous input validation prevents cryptic downstream errors by detecting
# invalid or missing parameters early in execution. Each required argument
# is checked, and enumerated arguments are validated against their allowed
# values.
#
# This defensive programming approach is essential for scripts that run
# in automated pipelines where immediate user intervention is not possible.

# Validate that all required arguments have been provided.
# The is.null() check detects arguments that were not specified on the
# command line and have no default value.
if (is.null(opt$direction))   stop_with("Missing --direction")
if (is.null(opt$kind))        stop_with("Missing --kind")
if (is.null(opt$mode))        stop_with("Missing --mode (expr|pca)")
if (is.null(opt$alg))         stop_with("Missing --alg (km|hc)")
if (is.null(opt$outdir))      stop_with("Missing --outdir")
if (is.null(opt$cluster_rds)) stop_with("Missing --cluster_rds")

# Validate enumerated parameters against their allowed values.
# The %in% operator tests whether the left operand is an element of the
# right operand vector, returning TRUE or FALSE.
if (!opt$mode %in% c("expr","pca")) stop_with("--mode must be expr or pca")
if (!opt$alg  %in% c("km","hc"))    stop_with("--alg must be km or hc")

# Create the output directory structure.
# The recursive = TRUE argument enables creation of nested directories
# (e.g., creating both 'results' and 'results/clustering' in one call).
# The showWarnings = FALSE argument suppresses warnings if the directory
# already exists, which is common in pipeline re-runs.
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# SECTION 5: CONFIGURATION LOADING
# ------------------------------------------------------------------------------
# The YAML configuration file centralises parameters that may change between
# pipeline runs, enabling separation of configuration from code logic.
# This design pattern improves reproducibility by documenting all analysis
# parameters in a single, version-controlled file.
#
# The configuration system supports hierarchical profiles, allowing
# disease-specific parameter overrides while inheriting common defaults.

# Determine the script directory relative to the configuration file.
# This enables loading of helper libraries that reside alongside the
# main analysis scripts.
config_dir <- dirname(normalizePath(opt$config))
script_dir <- file.path(config_dir, "..", "scripts")
lib_config_path <- file.path(script_dir, "lib_config.R")

# Attempt to load the profiled configuration helper library.
# This library provides functions for merging default parameters with
# profile-specific overrides, mimicking the Snakemake configuration system.
if (file.exists(lib_config_path)) {
  source(lib_config_path)
  
  # Determine the profile to use with the following priority:
  # 1. Command-line --profile argument (highest priority)
  # 2. SNAKEMAKE_PROFILE environment variable (set by Snakemake)
  # 3. Default profile specified in the configuration file
  profile_override <- if (!is.null(opt$profile)) opt$profile else Sys.getenv("SNAKEMAKE_PROFILE", unset = "")
  if (identical(profile_override, "")) profile_override <- NULL
  
  # Load and merge the configuration using the helper function.
  cfg <- read_profiled_config(opt$config, profile_override = profile_override)
  
} else {
  # Fallback: Direct YAML parsing without profile merging.
  # This maintains backward compatibility with older pipeline versions.
  cfg <- yaml::read_yaml(opt$config)
  
  # If the configuration has a defaults section, perform manual merging.
  if (!is.null(cfg$defaults)) {
    defaults <- cfg$defaults
    profile <- cfg$profile %||% names(cfg$profiles)[[1]] %||% "default"
    
    if (!is.null(cfg$profiles) && !is.null(cfg$profiles[[profile]])) {
      profcfg <- cfg$profiles[[profile]]
      
      # Simple merge: profile-specific values override defaults.
      # The c() function combines named lists, with later values taking
      # precedence for duplicate names.
      cfg <- list(
        paths = c(defaults$paths %||% list(), profcfg$paths %||% list()),
        methods = defaults$methods %||% list(),
        analysis = c(defaults$analysis %||% list(), profcfg$analysis %||% list()),
        features = c(defaults$features %||% list(), profcfg$features %||% list()),
        agnostic_clustering = defaults$agnostic_clustering %||% list(),
        clustering = defaults$clustering %||% list(),
        tumour_neighbourhoods = c(defaults$tumour_neighbourhoods %||% list(), 
                                  profcfg$tumour_neighbourhoods %||% list())
      )
    } else {
      cfg <- defaults
    }
  }
}

# ------------------------------------------------------------------------------
# PATH RESOLUTION
# ------------------------------------------------------------------------------
# File paths in the configuration may be specified relative to the pipeline
# root directory. This section resolves relative paths to absolute paths,
# ensuring correct file access regardless of the working directory from
# which the script is invoked.

pipe_root <- normalizePath(file.path(dirname(opt$config), ".."))

abs_from_root <- function(p) {
  # abs_from_root(): Converts a potentially relative path to an absolute path.
  #
  # If the path is already absolute (starts with '/' on Unix or a drive letter
  # on Windows), it is returned unchanged. Otherwise, it is resolved relative
  # to the pipeline root directory.
  #
  # Parameters:
  #   p: A file system path (character string)
  #
  # Returns:
  #   The absolute path (character string)
  
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(pipe_root, p)
}

# Resolve all paths in the configuration to absolute paths.
# This is performed recursively for all character values in the paths section.
if (!is.null(cfg$paths)) {
  for (key in names(cfg$paths)) {
    if (!is.null(cfg$paths[[key]]) && is.character(cfg$paths[[key]])) {
      cfg$paths[[key]] <- abs_from_root(cfg$paths[[key]])
    }
  }
}

# ------------------------------------------------------------------------------
# DIRECTION PARSING
# ------------------------------------------------------------------------------
# The direction string encodes both the feature selection method and the
# distance metric in a single identifier. This function parses the direction
# to extract these components for downstream processing.

parse_direction <- function(direction) {
  # parse_direction(): Extracts feature set and distance metric from direction string.
  #
  # Direction strings follow the pattern: {feature}_{distance}
  # Examples: 'HVG_euc', 'pam50_corr', 'Variance_euc', 'MAD_corr'
  #
  # Parameters:
  #   direction: The direction string to parse
  #
  # Returns:
  #   A list with two elements:
  #     - feature: The feature set name (e.g., 'HVG', 'PAM50', 'Variance')
  #     - distance: The distance metric (e.g., 'euclidean', 'correlation')
  
  # Parse general pattern: split on underscore and extract components.
  parts <- strsplit(direction, "_")[[1]]
  if (length(parts) < 2) {
    stop_with("Invalid direction: ", direction, " (expected pattern: {feature}_{distance})")
  }

  # The distance metric is always the last component.
  dist_tag <- tail(parts, 1)

  # The feature name may contain underscores (e.g., 'MeanAbsDev').
  feat <- paste(head(parts, -1), collapse = "_")

  # Normalise feature names for consistency.
  if (feat == "pam50") {
    feat <- "PAM50"
  }
  
  # Normalise distance metric names.
  dist <- switch(dist_tag,
                 euc = "euclidean",
                 corr = "correlation",
                 cos = "cosine",
                 dist_tag)
  
  list(feature = feat, distance = dist)
}

# Parse the direction to extract feature set and distance metric.
pd <- parse_direction(opt$direction)
feature <- pd$feature
distance <- pd$distance

# ------------------------------------------------------------------------------
# SECTION 6: CLUSTERING PARAMETER EXTRACTION
# ------------------------------------------------------------------------------
# The k_grid specifies the range of cluster numbers (k values) to evaluate.
# ConsensusClusterPlus computes consensus matrices for each k value and
# provides metrics for selecting the optimal number of clusters.

# Extract the k_grid from configuration and validate.
k_grid <- cfg$clustering$k_grid
if (is.null(k_grid) || length(k_grid) < 2) {
  stop_with("config: clustering.k_grid must have >=2 values")
}

# Convert to a sorted integer vector.
# The unlist() function handles any nested list structure from YAML parsing.
k_grid <- sort(as.integer(unlist(k_grid)))

# Set the random seed for reproducibility.
# This seed controls bootstrap sample selection in ConsensusClusterPlus,
# ensuring identical results across pipeline re-runs.
seed <- 42
if (!is.null(cfg$agnostic_clustering$seed)) {
  seed <- as.integer(cfg$agnostic_clustering$seed)
}

# Set the number of principal components for PCA mode.
# This parameter controls the dimensionality of the reduced feature space.
# Typical values range from 10 to 50, with 20 being a common default.
n_pcs <- 20
if (!is.null(cfg$agnostic_clustering$n_pcs)) {
  n_pcs <- as.integer(cfg$agnostic_clustering$n_pcs)
}

# ------------------------------------------------------------------------------
# SECTION 7: DISTANCE METRIC AND LINKAGE CONFIGURATION
# ------------------------------------------------------------------------------
# IMPORTANT METHODOLOGICAL NOTE ON DISTANCE METRICS
#
# ConsensusClusterPlus internally uses Euclidean distance for all calculations.
# When correlation-based clustering is desired, the input data undergoes a
# geometric transformation (see Section 9) that converts Euclidean distance
# to correlation-equivalent geometry.
#
# This approach is mathematically necessary because Ward's minimum variance
# criterion, used for linkage, is only valid with Euclidean distance. The
# transformation preserves the relative sample relationships while enabling
# the use of Ward's robust agglomerative method.
#
# MATHEMATICAL BASIS
# For two vectors x and y that are centred (mean = 0) and L2-normalised
# (||x|| = ||y|| = 1), the Euclidean distance relates to Pearson correlation:
#
#   ||x - y||^2 = 2(1 - cor(x, y))
#
# Therefore, minimising Euclidean distance on transformed data is equivalent
# to maximising Pearson correlation on the original data.

ccp_distance <- "euclidean"

# Ward.D2 linkage is used for both consensus matrix construction and final
# cluster assignment. Ward.D2 implements the correct Ward algorithm that
# minimises the total within-cluster variance at each merge step.
#
# Note: The original Ward algorithm as implemented in some software packages
# (ward.D) uses an incorrect formula. Ward.D2 provides the mathematically
# correct implementation.
innerLinkage <- "ward.D2"
finalLinkage <- "ward.D2"

# Log the analysis configuration for reproducibility documentation.
info("direction=%s | feature=%s | distance=%s | mode=%s | alg=%s | kind=%s",
     opt$direction, feature, distance, opt$mode, opt$alg, opt$kind)

# ------------------------------------------------------------------------------
# SECTION 8: EXPRESSION MATRIX CONSTRUCTION FUNCTION
# ------------------------------------------------------------------------------
# The build_expr_mat() function handles the complexity of loading and
# integrating expression data from multiple sources. It supports different
# feature sets (PAM50, HVG, method-specific) and sample scopes (cell-only,
# tumour-only, combined).
#
# The function returns a standardised output format regardless of input
# source, enabling consistent downstream processing.

# Helper function for detecting gene identifiers.
is_gene_id <- function(x) grepl("^ENSG", x)

# Helper function for extracting sample identifiers from a matrix.
get_sample_ids <- function(mat) {
  rn <- rownames(mat)
  cn <- colnames(mat)
  
  # Detect matrix orientation by checking which dimension contains gene IDs.
  if (!is.null(cn) && length(cn) > 0 && mean(is_gene_id(cn), na.rm = TRUE) > 0.5) {
    return(rn)  # Genes in columns, samples in rows
  }
  return(cn)  # Samples in columns, genes in rows
}

# Helper function for standardising matrix orientation.
ensure_samples_in_rows <- function(mat) {
  # ensure_samples_in_rows(): Transposes matrix if necessary to place samples in rows.
  #
  # The function detects the current orientation by examining whether column
  # names look like Ensembl gene identifiers (ENSG prefix). If genes are in
  # columns, samples must be in rows, and no transposition is needed.
  #
  # Parameters:
  #   mat: An expression matrix
  #
  # Returns:
  #   The matrix with samples as rows and genes as columns
  
  rn <- rownames(mat)
  cn <- colnames(mat)
  
  if (!is.null(cn) && length(cn) > 0 && mean(is_gene_id(cn), na.rm = TRUE) > 0.5) {
    return(as.matrix(mat))  # Already correct orientation
  }
  return(t(as.matrix(mat)))  # Transpose to put samples in rows
}

build_expr_mat <- function(feature, cfg, kind) {
  # build_expr_mat(): Constructs the expression matrix for clustering.
  #
  # This function handles the complexity of loading expression data from
  # multiple sources (TCGA tumours, DSMZ cell lines) and integrating them
  # according to the specified feature set and sample scope.
  #
  # Parameters:
  #   feature: The feature set to use ('PAM50', 'HVG', or method-specific)
  #   cfg:     The parsed configuration object
  #   kind:    The sample scope ('cell', 'tumour', or 'cell_tumour')
  #
  # Returns:
  #   A list with two elements:
  #     - expr_mat: A genes x samples matrix of expression values
  #     - dataset:  A named vector indicating sample origin ('TCGA' or 'DSMZ')
  
  # Determine the sample scope from the kind argument.
  is_cell_tumour <- grepl("_cell_tumour", kind)
  is_cell_only   <- grepl("_cell$", kind) && !is_cell_tumour
  is_tumour_only <- grepl("_tumour$", kind) && !is_cell_tumour

  # --------------------------------------------------------------------------
  # PAM50 FEATURE SET HANDLING
  # --------------------------------------------------------------------------
  # Pre-computed expression matrices containing only PAM50 genes are loaded
  # from RDS files. These matrices are generated upstream in the pipeline
  # to ensure consistent gene selection across analyses.
  
  if (feature == "PAM50" && (is.null(opt$feature_list) || !nzchar(opt$feature_list))) {
    tcga_path <- cfg$paths$tcga_brca_pam50_expr
    dsmz_path <- cfg$paths$dsmz_bcc_pam50_expr
    
    if (is.null(tcga_path) || is.null(dsmz_path)) {
      stop_with("Missing pam50 paths: paths.tcga_brca_pam50_expr or paths.dsmz_bcc_pam50_expr")
    }
    
    # Load pre-computed PAM50 expression matrices.
    tcga_mat <- readRDS(tcga_path)
    dsmz_mat <- readRDS(dsmz_path)

    # Identify genes present in both matrices.
    common_genes <- intersect(rownames(tcga_mat), rownames(dsmz_mat))
    
    if (length(common_genes) < 20) {
      stop_with("Too few common genes between TCGA and DSMZ PAM50 matrices.")
    }

    # Construct the output based on the specified sample scope.
    if (is_cell_tumour) {
      # Combine both datasets column-wise (samples are columns).
      expr_mat <- cbind(tcga_mat[common_genes, , drop = FALSE],
                        dsmz_mat[common_genes, , drop = FALSE])
      
      # Create dataset origin annotation vector.
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
  
  if (feature == "HVG") {
    cell_path   <- cfg$paths$cell_vst_rds
    tumour_path <- cfg$paths$tumour_vst_rds
    
    # Determine the gene list source.
    # Command-line argument takes precedence over configuration.
    hvg_list <- opt$feature_list
    if (!is.null(hvg_list) && nzchar(hvg_list)) {
      hvg_list <- abs_from_root(hvg_list)
    }

    if (is.null(cell_path) || is.null(tumour_path)) {
      stop_with("Missing: paths.cell_vst_rds and/or paths.tumour_vst_rds")
    }
    
    if (is.null(hvg_list) || !file.exists(hvg_list)) {
      stop_with("Missing HVG list file. Pass --feature_list pointing to the HVG gene list.")
    }

    # Load variance-stabilised transformed expression matrices.
    # VST transformation (from DESeq2) normalises variance across the
    # expression range, making genes with different expression levels
    # more comparable for distance calculations.
    cell_mat   <- readRDS(cell_path)
    tumour_mat <- readRDS(tumour_path)

    # Load the HVG list, filtering out empty lines.
    genes <- readLines(hvg_list)
    genes <- genes[nzchar(genes)]

    # Construct the output based on the specified sample scope.
    if (is_cell_tumour) {
      keep <- intersect(genes, intersect(rownames(cell_mat), rownames(tumour_mat)))
      
      if (length(keep) < 50) {
        stop_with("Too few HVG genes found in both cell and tumour matrices.")
      }
      
      expr_mat <- cbind(tumour_mat[keep, , drop = FALSE],
                        cell_mat[keep, , drop = FALSE])
      
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

  # --------------------------------------------------------------------------
  # METHOD-SPECIFIC FEATURE SETS
  # --------------------------------------------------------------------------
  cell_path <- cfg$paths$cell_vst_rds
  tumour_path <- cfg$paths$tumour_vst_rds
  feature_list <- opt$feature_list
  if (!is.null(feature_list) && nzchar(feature_list)) {
    feature_list <- abs_from_root(feature_list)
  }

  if (is.null(cell_path) || is.null(tumour_path)) {
    stop_with("Missing: paths.cell_vst_rds and/or paths.tumour_vst_rds")
  }

  if (is.null(feature_list) || !file.exists(feature_list)) {
    stop_with("Missing feature list file. Pass --feature_list pointing to the gene list.")
  }

  cell_mat <- readRDS(cell_path)
  tumour_mat <- readRDS(tumour_path)

  genes <- readLines(feature_list)
  genes <- genes[nzchar(genes)]

  if (is_cell_tumour) {
    keep <- intersect(genes, intersect(rownames(cell_mat), rownames(tumour_mat)))
    if (length(keep) < 50) {
      stop_with("Too few feature genes found in both cell and tumour matrices.")
    }
    expr_mat <- cbind(tumour_mat[keep, , drop = FALSE],
                      cell_mat[keep, , drop = FALSE])
    dataset <- c(
      setNames(rep("TCGA", ncol(tumour_mat[keep, , drop = FALSE])), colnames(tumour_mat)),
      setNames(rep("DSMZ", ncol(cell_mat[keep, , drop = FALSE])), colnames(cell_mat))
    )
    dataset <- dataset[colnames(expr_mat)]
    names(dataset) <- colnames(expr_mat)
  } else if (is_cell_only) {
    keep <- intersect(genes, rownames(cell_mat))
    if (length(keep) < 50) {
      stop_with("Too few feature genes found in cell matrix.")
    }
    expr_mat <- cell_mat[keep, , drop = FALSE]
    dataset <- setNames(rep("DSMZ", ncol(expr_mat)), colnames(expr_mat))
  } else if (is_tumour_only) {
    keep <- intersect(genes, rownames(tumour_mat))
    if (length(keep) < 50) {
      stop_with("Too few feature genes found in tumour matrix.")
    }
    expr_mat <- tumour_mat[keep, , drop = FALSE]
    dataset <- setNames(rep("TCGA", ncol(expr_mat)), colnames(expr_mat))
  } else {
    stop_with("Cannot determine scope from kind: ", kind)
  }

  return(list(expr_mat = expr_mat, dataset = dataset))
}

# ------------------------------------------------------------------------------
# SECTION 9: CORRELATION-GEOMETRY TRANSFORMATION
# ------------------------------------------------------------------------------
# This function implements the mathematical transformation that enables
# correlation-based clustering using Euclidean distance calculations.
#
# MATHEMATICAL BACKGROUND
# -----------------------
# Pearson correlation between two vectors x and y is defined as:
#
#   cor(x, y) = sum((x - mean(x)) * (y - mean(y))) / (||x - mean(x)|| * ||y - mean(y)||)
#
# If we centre each vector (subtract its mean) and normalise to unit length,
# the correlation simplifies to the dot product:
#
#   cor(x_transformed, y_transformed) = x_transformed . y_transformed
#
# Furthermore, for unit-length vectors, Euclidean distance relates to the
# dot product (and thus correlation) as:
#
#   ||x - y||^2 = ||x||^2 + ||y||^2 - 2(x . y) = 2 - 2(x . y) = 2(1 - cor(x, y))
#
# Therefore, minimising Euclidean distance on transformed data is equivalent
# to maximising Pearson correlation on the original data.

corr_geometry_transform <- function(X) {
  # corr_geometry_transform(): Transforms data for correlation-based clustering.
  #
  # The transformation consists of two steps:
  # 1. Centre each column (sample) by subtracting its mean
  # 2. Normalise each column to unit L2 length
  #
  # After transformation, Euclidean distance between columns is equivalent
  # to correlation-based distance on the original data.
  #
  # Parameters:
  #   X: A features x samples matrix (genes as rows, samples as columns)
  #
  # Returns:
  #   The transformed matrix with centred, unit-normalised columns
  
  # Centre each column by subtracting its mean.
  # The sweep() function applies an operation across matrix margins.
  # MARGIN = 2 specifies column-wise operation.
  X <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  
  # Compute the L2 norm (Euclidean length) of each column.
  norms <- sqrt(colSums(X^2, na.rm = TRUE))
  
  # Handle edge cases where the norm is zero or non-finite.
  # This prevents division by zero for constant columns (zero variance).
  norms[!is.finite(norms) | norms == 0] <- 1
  
  # Normalise each column to unit length.
  X <- sweep(X, 2, norms, "/")
  
  X
}

# ------------------------------------------------------------------------------
# SECTION 10: EXPRESSION MATRIX CONSTRUCTION AND PREPROCESSING
# ------------------------------------------------------------------------------
# The expression matrix is constructed using the build_expr_mat() function,
# which handles the complexity of loading and integrating data from
# multiple sources based on the specified feature set and sample scope.

built    <- build_expr_mat(feature, cfg, opt$kind)
expr_mat <- built$expr_mat
dataset  <- built$dataset

# Remove genes with zero variance (constant expression across all samples).
# Such genes provide no discriminatory information for clustering and can
# cause numerical issues in distance calculations and PCA.
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
# as columns. Two input modes are supported:
#
# Expression mode (expr): Uses gene expression values directly as features.
# This preserves the full gene-level information but may include noise from
# uninformative genes.
#
# PCA mode (pca): Reduces dimensionality by projecting samples onto the
# top principal components. This can improve signal-to-noise ratio by
# concentrating variance into a smaller number of features, but loses
# direct gene-level interpretability.

if (opt$mode == "expr") {
  # Expression mode: use the gene expression matrix directly.
  data_for_ccp <- expr_mat
  info("CCP input = expr (genes x samples): %d x %d", nrow(data_for_ccp), ncol(data_for_ccp))
  
} else {
  # PCA mode: perform principal component analysis on samples.
  
  # Transpose to samples x genes orientation for prcomp().
  samples_mat <- t(expr_mat)
  
  # Perform PCA with centring but without scaling.
  # Centring removes the mean of each gene across samples.
  # Scaling is omitted to preserve relative variance differences between genes,
  # which may carry biological significance.
  pca <- prcomp(samples_mat, center = TRUE, scale. = FALSE)
  
  # Extract the top n_pcs principal component scores.
  # The $x component contains the sample scores (projections onto PC axes).
  pcs_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
  
  # Transpose to PCs x samples orientation for ConsensusClusterPlus.
  data_for_ccp <- t(pcs_mat)
  info("CCP input = pca (PCs x samples): %d x %d", nrow(data_for_ccp), ncol(data_for_ccp))
}

# Apply correlation-geometry transformation if correlation distance is specified.
if (distance == "correlation") {
  info("Applying correlation-geometry transform (center + unit norm) for Ward.D2 compatibility.")
  data_for_ccp <- corr_geometry_transform(data_for_ccp)
}

# ------------------------------------------------------------------------------
# K-GRID VALIDATION
# ------------------------------------------------------------------------------
# The k_grid must be capped based on sample size to prevent algorithm failures.
# K-means clustering requires k < n_samples, and with bootstrap subsampling
# (pItem < 1.0), the effective sample size is further reduced.

n_samples <- ncol(data_for_ccp)
n_features <- nrow(data_for_ccp)

# Determine maximum allowed k based on algorithm constraints.
k_max_allowed <- if (opt$alg == "km") {
  min(n_features, n_samples) - 1
} else {
  n_samples - 1
}

# Filter k_grid to valid values.
k_grid_effective <- k_grid[k_grid <= k_max_allowed]

if (length(k_grid_effective) < 2) {
  if (k_max_allowed >= 2) {
    k_grid_effective <- c(2, min(3, k_max_allowed))
    info("Sample size too small for full k_grid. Using minimal grid: %s", 
         paste(k_grid_effective, collapse = ","))
  } else {
    stop_with("Too few samples for consensus clustering. Need at least k_max=2.")
  }
}

# Account for bootstrap subsampling in maximum k calculation.
pItem <- 0.8
n_samples_effective <- floor(pItem * n_samples)
maxK_safe <- min(max(k_grid_effective), n_samples_effective - 1, n_features - 1)

if (maxK_safe < 2) {
  maxK_safe <- 2
}

# Final k_grid filtering.
k_grid_effective <- k_grid_effective[k_grid_effective <= maxK_safe]

if (length(k_grid_effective) < 2) {
  if (maxK_safe >= 2) {
    k_grid_effective <- c(2, min(3, maxK_safe))
  } else {
    stop_with("Too few samples for consensus clustering.")
  }
}

info("Final k_grid: %s (maxK=%d)", paste(k_grid_effective, collapse = ","), maxK_safe)

# ------------------------------------------------------------------------------
# SECTION 12: CONSENSUS CLUSTERING EXECUTION
# ------------------------------------------------------------------------------
# ConsensusClusterPlus outputs multiple diagnostic files to the working
# directory. A dedicated subdirectory is created to isolate these outputs.

ccp_wd <- file.path(opt$outdir, "ccp")
dir.create(ccp_wd, showWarnings = FALSE, recursive = TRUE)

# Save the current working directory and switch to the CCP output directory.
old_wd <- getwd()
original_wd <- old_wd
setwd(ccp_wd)
on.exit(setwd(old_wd), add = TRUE)

# Construct the argument list for ConsensusClusterPlus.
ccp_args <- list(
  d            = data_for_ccp,     # Input data matrix (features x samples)
  maxK         = maxK_safe,        # Maximum number of clusters to evaluate
  reps         = 1000,             # Number of bootstrap iterations
  pItem        = pItem,            # Proportion of samples per bootstrap
  pFeature     = 1.0,              # Proportion of features (1.0 = all)
  clusterAlg   = opt$alg,          # Clustering algorithm (km or hc)
  distance     = ccp_distance,     # Distance metric (always Euclidean)
  innerLinkage = innerLinkage,     # Linkage for consensus matrix
  finalLinkage = finalLinkage,     # Linkage for final assignment
  seed         = seed,             # Random seed for reproducibility
  title        = paste0("CCP_", opt$direction, "_", opt$kind),
  plot         = "png",            # Output format for diagnostic plots
  writeTable   = TRUE,             # Write cluster assignments to CSV
  verbose      = FALSE
)

# Execute ConsensusClusterPlus.
# The result is a list indexed by k, containing:
#   - consensusMatrix: The n x n sample co-clustering frequency matrix
#   - consensusClass:  The cluster assignment vector
#   - consensusTree:   The hierarchical clustering dendrogram
ccp_res <- do.call(ConsensusClusterPlus, ccp_args)

# ------------------------------------------------------------------------------
# SECTION 13: OPTIMAL CLUSTER NUMBER SELECTION
# ------------------------------------------------------------------------------
# The optimal number of clusters is determined using the area under the
# cumulative distribution function (CDF) of consensus values.
#
# INTERPRETATION
# --------------
# The consensus matrix M[i,j] contains values between 0 and 1 representing
# how frequently samples i and j were assigned to the same cluster across
# bootstrap iterations.
#
# Ideal clustering produces a bimodal distribution:
#   - M[i,j] near 0 for samples in different clusters (never co-cluster)
#   - M[i,j] near 1 for samples in the same cluster (always co-cluster)
#
# The CDF of consensus values approaches a step function for ideal clustering.
# The area under the CDF (AUC) quantifies this: higher AUC indicates sharper
# separation between co-clustered and non-co-clustered sample pairs.
#
# Delta AUC (the change in AUC between successive k values) identifies the
# k at which cluster stability improves most dramatically.

auc_at_k <- function(k) {
  # auc_at_k(): Computes AUC of the consensus CDF for a given k.
  #
  # Parameters:
  #   k: The number of clusters
  #
  # Returns:
  #   The area under the empirical CDF of consensus values
  
  cm <- ccp_res[[k]]$consensusMatrix
  vals <- cm[upper.tri(cm)]
  ec <- ecdf(vals)
  x <- seq(0, 1, length.out = 500)
  sum(diff(x) * (ec(head(x, -1)) + ec(tail(x, -1))) / 2)
}

# Compute AUC for each k value.
auc_values <- vapply(k_grid_effective, auc_at_k, numeric(1))

# Compute delta AUC (improvement in AUC between successive k values).
delta_auc <- diff(auc_values)

# Select k with maximum delta AUC.
best_k <- k_grid_effective[which.max(delta_auc) + 1]

info("k_grid: %s", paste(k_grid_effective, collapse = ","))
info("AUC:    %s", paste(round(auc_values, 4), collapse = " "))
info("delta_AUC: %s", paste(round(delta_auc, 4), collapse = " "))
info("best_k: %d", best_k)

# ------------------------------------------------------------------------------
# SECTION 14: FINAL CLUSTER ASSIGNMENT EXTRACTION
# ------------------------------------------------------------------------------
# Extract cluster assignments from the optimal k solution and format them
# with descriptive labels for downstream analyses.

final_clusters <- ccp_res[[best_k]]$consensusClass

# Add "C" prefix to cluster numbers for human readability.
final_clusters <- setNames(paste0("C", final_clusters), names(final_clusters))

# ------------------------------------------------------------------------------
# SECTION 15: OUTPUT FILE GENERATION
# ------------------------------------------------------------------------------
# Two output files are generated:
# 1. A CSV file containing sample-level cluster assignments
# 2. An RDS file containing the complete analysis results

# Ensure output directory exists with absolute path.
outdir_abs <- normalizePath(opt$outdir, mustWork = FALSE)
if (!dir.exists(outdir_abs)) {
  dir.create(outdir_abs, showWarnings = FALSE, recursive = TRUE)
}

# Write cluster assignments to CSV.
clusters_csv <- file.path(outdir_abs, sprintf("tcga_dsmz_ccp_clusters_K%d.csv", best_k))

write_csv(
  tibble(
    sample_id = names(final_clusters),
    cluster   = unname(final_clusters),
    dataset   = dataset[names(final_clusters)]
  ),
  clusters_csv
)

# Save comprehensive results object as RDS.
# This file enables programmatic access to all analysis parameters and results.
cluster_rds_abs <- if (startsWith(opt$cluster_rds, "/")) {
  opt$cluster_rds
} else {
  file.path(original_wd, opt$cluster_rds)
}

cluster_rds_dir <- dirname(cluster_rds_abs)
if (!dir.exists(cluster_rds_dir)) {
  dir.create(cluster_rds_dir, showWarnings = FALSE, recursive = TRUE)
}

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
    k_grid         = k_grid_effective,
    auc_values     = setNames(auc_values, k_grid_effective),
    delta_auc      = delta_auc,
    best_k         = best_k,
    clusters       = final_clusters,
    dataset        = dataset,
    ccp_results    = ccp_res,
    timestamp      = Sys.time()
  ),
  cluster_rds_abs
)

info("Wrote: %s", clusters_csv)
info("Wrote: %s", opt$cluster_rds)
info("Done.")
