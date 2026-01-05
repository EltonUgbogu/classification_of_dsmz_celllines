#!/usr/bin/env Rscript

# ==============================================================================
# comp_tumour_neighbourhoods.R
# Tumour Neighbourhood Computation Pipeline
# ==============================================================================
#
# PURPOSE:
# This script computes tumour neighbourhoods for each cell line in an integrated
# expression matrix. For each DSMZ breast cancer cell line, the pipeline
# identifies the most similar TCGA tumour samples based on expression distance,
# enabling biological interpretation of cell line-tumour relationships within
# the context of consensus clustering results.
#
# USAGE:
# Rscript comp_tumour_neighbourhoods.R \
#   --config config/config.yaml \
#   --direction <pam50_euc|pam50_corr|hvg_euc|hvg_corr> \
#   [--expr-rds <path>] \
#   [--mapping-rds <path>]
#
# ==============================================================================

# Disable automatic string-to-factor conversion for consistent data handling.
# This is a legacy R behaviour that can cause unexpected issues with character data.
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The suppressPackageStartupMessages() function prevents verbose startup
# messages from cluttering the console output, which is particularly useful
# when running the script in automated pipelines.

suppressPackageStartupMessages({
  library(dplyr)    # Data manipulation with tidyverse grammar
  library(readr)    # Efficient CSV reading and writing
  library(yaml)     # YAML configuration file parsing
  library(optparse) # Command-line argument parsing
})

cat("=== Starting computation of tumour neighbourhoods ===\n")

# ------------------------------------------------------------------------------
# SECTION 2: COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
# The optparse package provides a systematic approach to defining and parsing
# command-line arguments. Each make_option() call defines a single argument
# with its type, default value, and help text for user documentation.

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to YAML config file [default: %default]"),
  make_option("--expr-rds", type = "character", default = NULL,
              help = "Optional override: integrated DSMZ+TCGA expression RDS"),
  make_option("--mapping-rds", type = "character", default = NULL,
              help = "Optional override: cell_line → original_sample_id mapping RDS"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr)")
)

# parse_args() processes the command-line arguments according to the defined
# option_list and returns a named list containing the parsed values.
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Using config file:", opt$config, "\n")
cfg <- yaml::read_yaml(opt$config)

# ------------------------------------------------------------------------------
# SECTION 3: NULL-COALESCING OPERATOR DEFINITION
# ------------------------------------------------------------------------------
# The %||% operator provides a concise syntax for supplying default values
# when a variable may be NULL. This pattern is commonly used in R for
# handling optional configuration parameters.
#
# USAGE: value %||% default
#   Returns 'value' if it is not NULL, otherwise returns 'default'.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# SECTION 4: HELPER FUNCTION LOADING
# ------------------------------------------------------------------------------
# The pipeline relies on modular helper functions stored in a separate
# directory. This approach promotes code reuse and maintainability by
# separating core algorithms from pipeline orchestration logic.

if (is.null(cfg$paths$tumour_nh_base_functions)) {
  stop("Config must define paths$tumour_nh_base_functions")
}

base_fun_dir <- cfg$paths$tumour_nh_base_functions

if (!dir.exists(base_fun_dir)) {
  stop("Base-functions directory does not exist: ", base_fun_dir)
}

# Discover all R script files in the helper directory.
# The list.files() function with pattern matching identifies files by extension.
helper_files <- list.files(
  base_fun_dir,
  pattern = "\\.R$",
  full.names = TRUE
)

if (length(helper_files) == 0L) {
  stop("No .R helper files found in base_functions dir: ", base_fun_dir)
}

# Source all helper files to load their functions into the global environment.
# The invisible() wrapper suppresses the return value output from source().
cat("[INFO] Loading", length(helper_files), "helper files from:", base_fun_dir, "\n")
invisible(lapply(helper_files, source))

# ------------------------------------------------------------------------------
# SECTION 5: DIRECTION PARAMETER VALIDATION
# ------------------------------------------------------------------------------
# The "direction" parameter specifies the feature set and distance metric
# combination used for neighbourhood computation. Valid directions are
# defined in the configuration file, allowing project-specific customisation.
#
# Direction format: <feature_set>_<distance_metric>
#   - Feature sets: pam50 (PAM50 gene signature), hvg (highly variable genes)
#   - Distance metrics: euc (Euclidean), corr (correlation)

if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Retrieve allowed directions from configuration, with fallback defaults.
# This design allows different cancer types to support different feature sets
# (e.g., PAM50 is specific to breast cancer).
directions <- cfg$tumour_neighbourhoods$directions
if (is.null(directions) || length(directions) == 0) {
  directions <- c("hvg_euc", "hvg_corr")  # Fallback for non-BRCA analyses
}

if (!direction %in% directions) {
  stop("Invalid --direction: ", direction,
       "\nAllowed (from config tumour_neighbourhoods$directions): ",
       paste(directions, collapse = ", "))
}

# ------------------------------------------------------------------------------
# SECTION 6: DIRECTION PARSING AND PATH CONFIGURATION
# ------------------------------------------------------------------------------
# The direction string is parsed to extract the gene set and distance metric,
# which determine input file paths and computational parameters.

unsup_root <- cfg$paths$unsup_root

# Extract gene set from direction prefix.
# The startsWith() function performs efficient prefix matching.
gene_set <- if (startsWith(direction, "pam50")) "PAM50" else "HVG500"

# Extract distance type from direction suffix.
# The grepl() function performs regular expression matching.
dist_type <- if (grepl("_corr$", direction)) "correlation" else "euclidean"

cat("[INFO] Direction: ", direction, "\n", sep = "")
cat("[INFO] Gene set: ", gene_set, "\n", sep = "")
cat("[INFO] Distance metric: ", dist_type, "\n", sep = "")

# Define output directory structure within the unsupervised analysis tree.
base_functions_dir <- cfg$paths$tumour_nh_base_functions
tn_results_root    <- file.path(unsup_root, "tumour_neighbourhoods", direction)

# ------------------------------------------------------------------------------
# SECTION 7: INPUT FILE PATH RESOLUTION
# ------------------------------------------------------------------------------
# Input file paths are resolved using a priority system:
#   1. Command-line overrides (highest priority)
#   2. Configuration file paths (default)
#
# This design allows flexible execution in both interactive and pipeline contexts.

# Resolve expression matrix path based on gene set.
expr_mat_path <- opt$expr_rds %||%
  if (gene_set == "PAM50") {
    cfg$paths$tumour_nh_expr_pam50
  } else if (gene_set == "HVG500") {
    cfg$paths$tumour_nh_expr_hvg
  } else {
    stop("Unsupported gene_set: ", gene_set, ". Must be PAM50 or HVG500.")
  }

# Resolve sample ID mapping path.
# The mapping file links technical sample identifiers to canonical cell line names.
mapping_path <- opt$mapping_rds %||% {
  # Infer mapping file location from expression matrix directory.
  map_suffix <- if (startsWith(direction, "pam50")) {
    "cell_line_to_original_sample_id_pam50.rds"
  } else {
    "cell_line_to_original_sample_id_hvg.rds"
  }
  file.path(dirname(expr_mat_path), map_suffix)
}

cat("[INFO] Expression matrix path:\n  ", expr_mat_path, "\n", sep = "")
cat("[INFO] Cell line mapping path:\n  ", mapping_path, "\n", sep = "")

# Validate that required input files exist before proceeding.
# Early validation prevents cryptic errors during downstream processing.
if (!file.exists(expr_mat_path)) {
  stop("Expression matrix RDS not found at: ", expr_mat_path)
}
if (!file.exists(mapping_path)) {
  stop("cell_line_to_original_sample_id RDS not found at: ", mapping_path)
}

# ------------------------------------------------------------------------------
# SECTION 8: CORE HELPER FUNCTION LOADING
# ------------------------------------------------------------------------------
# Two essential helper modules are loaded:
#   - tumour_nh_io.R: Input/output utilities and path management
#   - tumour_neighbourhood.R: Core neighbourhood computation algorithm

cat("[INFO] Sourcing base functions from:\n  ", base_functions_dir, "\n", sep = "")

for (f in c("tumour_nh_io.R", "tumour_neighbourhood.R")) {
  src_path <- file.path(base_functions_dir, f)
  if (!file.exists(src_path)) {
    stop("Required helper not found: ", src_path)
  }
  # Source without local=TRUE to make functions available in global scope.
  source(src_path)
}

# Verify that expected functions were successfully loaded.
# These safety checks catch configuration or file structure errors early.
if (!exists("make_nh_paths") && !exists("nh_paths")) {
  stop("Neither nh_paths nor make_nh_paths() is defined after sourcing tumour_nh_io.R")
}
if (!exists("compute_tumour_neighbourhoods")) {
  stop("compute_tumour_neighbourhoods() not found after sourcing tumour_neighbourhood.R")
}

# Initialise path configuration object if factory function is available.
if (exists("make_nh_paths")) {
  nh_paths <- make_nh_paths(unsup_root)
  cat("[INFO] Built nh_paths using unsup_root:", unsup_root, "\n")
}

# ------------------------------------------------------------------------------
# SECTION 9: SAMPLE IDENTIFIER NORMALISATION FUNCTION
# ------------------------------------------------------------------------------
# Sample identifiers from different sources (TCGA, DSMZ) use inconsistent
# formatting conventions. This function standardises identifiers to enable
# reliable matching across datasets.
#
# TRANSFORMATIONS:
#   - Removes dataset prefixes (CELL:, TUMOUR:)
#   - Converts hyphens to underscores (CAL-120 → CAL_120)
#   - Converts spaces to underscores
#
# PARAMETERS:
#   x: A character vector of sample identifiers
#
# RETURNS:
#   A character vector with normalised identifiers

normalize_id <- function(x) {
  # Remove AGN (Anthropic Genomics Notation) prefixes if present.
  # These prefixes may be added during earlier pipeline stages.
  x <- sub("^(CELL:|TUMOUR:)", "", x)
  
  # Standardise delimiter characters to underscores.
  x <- gsub("-", "_", x)
  x <- gsub(" ", "_", x)
  x
}

# ------------------------------------------------------------------------------
# SECTION 10: TECHNICAL IDENTIFIER PARSING FUNCTION
# ------------------------------------------------------------------------------
# DSMZ samples have technical identifiers that include library preparation
# information. This function extracts the canonical cell line name from
# these verbose identifiers.
#
# EXAMPLE:
#   Input:  "NG_29643_CAL_120_lib581301_8005_3"
#   Output: "CAL_120"
#
# PARAMETERS:
#   x: A character vector of technical DSMZ identifiers (already normalised)
#
# RETURNS:
#   A character vector of extracted cell line names

extract_cell_line_from_technical_id <- function(x) {
  # Remove the NG_XXXXX_ prefix (sample accession number).
  x <- gsub("^NG_[0-9]+_", "", x)
  
  # Remove the _libXXXXXX_... suffix (library preparation metadata).
  x <- gsub("_lib[0-9_]+$", "", x)
  x
}

# ------------------------------------------------------------------------------
# SECTION 11: EXPRESSION MATRIX LOADING
# ------------------------------------------------------------------------------
# The integrated expression matrix contains both TCGA tumour samples and
# DSMZ cell line samples in a unified feature space. The matrix is stored
# in RDS format for efficient R-native serialisation.

expr_mat <- readRDS(expr_mat_path)
cat("[INFO] Expression matrix:", nrow(expr_mat), "samples x", ncol(expr_mat), "genes\n")

# ------------------------------------------------------------------------------
# SECTION 12: SAMPLE IDENTIFIER MAPPING LOADING
# ------------------------------------------------------------------------------
# The mapping file provides a lookup table from technical sample identifiers
# (used in raw data files) to canonical cell line names (used for biological
# interpretation).
#
# STRUCTURE:
#   names(mapping)  = Original technical identifiers (NG-29643_CAL_120_lib...)
#   values(mapping) = Canonical cell line names (CAL-120, CAL-51, ...)

orig_to_cellline <- readRDS(mapping_path)
cat("[INFO] Loaded original_sample_id → cell_line mapping (",
    length(orig_to_cellline), " entries)\n", sep = "")

# Extract mapping components for processing.
original_ids_raw    <- names(orig_to_cellline)
cell_line_names_raw <- unname(orig_to_cellline)

# Log example mappings for debugging and verification.
cat("[DEBUG] Example mapping keys (original IDs):\n")
print(head(original_ids_raw, 5))
cat("[DEBUG] Example mapping values (cell lines):\n")
print(head(cell_line_names_raw, 5))

current_rownames <- rownames(expr_mat)
cat("[DEBUG] Example expr_mat rownames:\n")
print(head(current_rownames, 10))

# ------------------------------------------------------------------------------
# SECTION 13: DATASET ORIGIN ANNOTATION
# ------------------------------------------------------------------------------
# Each sample is annotated with its source dataset (TCGA or DSMZ) BEFORE
# any identifier renaming occurs. This order of operations is critical
# because the original identifiers are needed to correctly identify DSMZ
# samples via the mapping file.
#
# IMPORTANT: Dataset annotation must precede identifier normalisation to
# ensure accurate provenance tracking.

# Create a logical mask identifying DSMZ samples by matching against mapping keys.
dsmz_mask <- current_rownames %in% original_ids_raw
n_dsmz    <- sum(dsmz_mask)

if (n_dsmz == 0) {
  stop("FATAL: No DSMZ samples found in expr_mat rownames using mapping values.\n",
       "  Example rownames: ", paste(head(current_rownames, 5), collapse = ", "),
       "\n  Example mapping values: ", paste(head(original_ids_raw, 5), collapse = ", "))
}

cat("[INFO] Detected ", n_dsmz, " DSMZ samples in expr_mat via mapping\n", sep = "")

# Assign dataset labels using vectorised conditional assignment.
# The ifelse() function evaluates element-wise, returning "DSMZ" for TRUE
# positions and "TCGA" for FALSE positions.
dataset_vec <- ifelse(dsmz_mask, "DSMZ", "TCGA")
names(dataset_vec) <- current_rownames

# ------------------------------------------------------------------------------
# SECTION 14: IDENTIFIER MAPPING CONSTRUCTION
# ------------------------------------------------------------------------------
# Two mapping objects are created:
#   1. orig_to_cell: Maps original IDs to normalised cell line names
#   2. orig_to_cell_norm: Maps normalised original IDs to cell line names
#
# The second mapping handles edge cases where input identifiers have already
# undergone partial normalisation.

# Normalise cell line names for consistent downstream matching.
cell_line_norm <- normalize_id(cell_line_names_raw)
orig_to_cell   <- setNames(cell_line_norm, original_ids_raw)

# Create a parallel mapping using normalised original IDs as keys.
orig_ids_norm <- normalize_id(original_ids_raw)
orig_to_cell_norm <- setNames(cell_line_norm, orig_ids_norm)

# ------------------------------------------------------------------------------
# SECTION 15: EXPRESSION MATRIX IDENTIFIER RENAMING
# ------------------------------------------------------------------------------
# DSMZ sample identifiers in the expression matrix are replaced with their
# canonical cell line names. This transformation simplifies downstream
# interpretation and enables consistent cross-dataset comparisons.

# Identify the subset of rownames that require renaming (DSMZ samples only).
ids_to_replace <- current_rownames[dsmz_mask]

# Apply the mapping to replace technical IDs with cell line names.
current_rownames[dsmz_mask] <- orig_to_cell[ids_to_replace]

# Update the expression matrix rownames.
rownames(expr_mat) <- current_rownames

# Apply final normalisation to all rownames and synchronise dataset vector.
rownames(expr_mat) <- normalize_id(rownames(expr_mat))
names(dataset_vec) <- rownames(expr_mat)

cat("[INFO] Successfully renamed ", n_dsmz, " DSMZ samples in expr_mat.\n", sep = "")
cat("Verification of DSMZ sample names (first 5):\n")
print(head(rownames(expr_mat)[dataset_vec == "DSMZ"], 5))

cat("[INFO] Total samples: ", nrow(expr_mat),
    " ( ", sum(dataset_vec == "DSMZ"), " DSMZ )\n", sep = "")

# ------------------------------------------------------------------------------
# SECTION 16: CLUSTERING METHOD DISCOVERY
# ------------------------------------------------------------------------------
# The pipeline dynamically discovers available clustering results rather than
# relying on hardcoded paths. The get_nh_methods() function enumerates actual
# pipeline outputs, ensuring robustness to changes in upstream processing.

methods <- get_nh_methods(unsup_root = unsup_root, direction = direction)

cat("\n[INFO] Methods discovered:\n")
print(methods)

# Filter to only methods with existing cluster files.
methods_exist <- methods %>% dplyr::filter(exists)

if (nrow(methods_exist) == 0) {
  stop("FATAL: No cluster RDS found for direction=", direction,
       "\nChecked:\n", paste(methods$path, collapse = "\n"))
}

cat(sprintf("\n[INFO] %d/%d methods available (existing files)\n",
            nrow(methods_exist), nrow(methods)))

# ------------------------------------------------------------------------------
# SECTION 17: SINGLE-METHOD NEIGHBOURHOOD COMPUTATION FUNCTION
# ------------------------------------------------------------------------------
# This function encapsulates the complete workflow for computing tumour
# neighbourhoods using a single clustering method. It handles cluster loading,
# identifier alignment, neighbourhood computation, and result serialisation.
#
# PARAMETERS:
#   path:      File path to the cluster assignment RDS file
#   method_id: Unique identifier for the clustering method
#   outdir:    Output directory for results
#
# RETURNS:
#   A list containing neighbourhood results, or NULL if processing fails

run_single_neighbourhood <- function(path, method_id, outdir) {
  cat("\n============================================================\n")
  cat("Method:", method_id, "\n")
  cat("Loading clusters from:", path, "\n")

  # --------------------------------------------------------------------------
  # STEP 1: CLUSTER FILE VALIDATION AND LOADING
  # --------------------------------------------------------------------------
  
  # Guard clause: skip if cluster file does not exist.
  if (!file.exists(path)) {
    message("[WARN] Cluster file not found for method ", method_id, ": ", path,
            " — skipping.")
    return(NULL)
  }

  # Load the cluster object and validate its structure.
  clust_obj <- readRDS(path)
  if (is.null(clust_obj$clusters)) {
    stop("Object at ", path, " does not have a $clusters element.")
  }

  cluster_vec <- clust_obj$clusters
  if (is.null(names(cluster_vec))) {
    stop("Cluster vector in ", path, " has no names.")
  }

  # --------------------------------------------------------------------------
  # STEP 2: CLUSTER IDENTIFIER NORMALISATION
  # --------------------------------------------------------------------------
  # Cluster vector identifiers must match expression matrix rownames exactly.
  # Multiple normalisation strategies are applied to handle various input formats.

  # Apply standard normalisation (hyphens to underscores, etc.).
  names(cluster_vec) <- normalize_id(names(cluster_vec))

  # Attempt to map any remaining original NG IDs to cell line names.
  nm  <- names(cluster_vec)
  hit <- nm %in% names(orig_to_cell_norm)
  if (any(hit)) {
    nm[hit] <- orig_to_cell_norm[nm[hit]]
    names(cluster_vec) <- nm
    cat("[INFO] Mapped", sum(hit),
        "cluster IDs from NG_* to cell-line names for", method_id, "\n")
  }
  
  # Fallback: extract cell line names from any remaining technical IDs.
  # This handles cases where the mapping may be incomplete.
  is_ng <- grepl("^NG_[0-9]+_", nm)
  if (any(is_ng)) {
    nm[is_ng] <- extract_cell_line_from_technical_id(nm[is_ng])
    names(cluster_vec) <- nm
    cat("[INFO] Fallback: extracted cell-line names from", sum(is_ng),
        "NG_* technical IDs for", method_id, "\n")
  }

  # Log DSMZ sample names for verification.
  dsmz_names_in_clusters <- names(cluster_vec)[!grepl("^TCGA-", names(cluster_vec))]
  cat("[INFO] First DSMZ names in cluster_vec: ",
      paste(head(dsmz_names_in_clusters, 5), collapse = ", "), "\n")

  # --------------------------------------------------------------------------
  # STEP 3: SAMPLE ALIGNMENT BETWEEN EXPRESSION AND CLUSTER DATA
  # --------------------------------------------------------------------------
  # The expression matrix and cluster vector may contain different sample sets.
  # Only samples present in both are used for neighbourhood computation.

  common_ids <- intersect(rownames(expr_mat), names(cluster_vec))

  if (length(common_ids) == 0) {
    stop("No overlapping sample IDs between expr_mat and cluster_vec for ", method_id,
         "\n  expr_mat sample IDs (first 5): ", paste(head(rownames(expr_mat), 5), collapse = ", "),
         "\n  cluster_vec sample IDs (first 5): ", paste(head(names(cluster_vec), 5), collapse = ", "))
  }

  # Report samples present in expression matrix but missing from clusters.
  missing_in_clusters <- setdiff(rownames(expr_mat), names(cluster_vec))
  if (length(missing_in_clusters) > 0) {
    cat("[WARNING] ", length(missing_in_clusters),
        " samples in expr_mat not found in cluster_vec for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_clusters, 5), collapse = ", "), "\n", sep = "")
  }

  # Report samples present in clusters but missing from expression matrix.
  missing_in_expr <- setdiff(names(cluster_vec), rownames(expr_mat))
  if (length(missing_in_expr) > 0) {
    dsmz_only <- missing_in_expr[!grepl("^TCGA-", missing_in_expr)]
    cat("[WARNING] ", length(missing_in_expr),
        " samples in cluster_vec not found in expr_mat for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_expr, 5), collapse = ", "), "\n", sep = "")
    if (length(dsmz_only) > 0) {
      cat("  DSMZ-only in cluster_vec (not in expr_mat): ",
          paste(head(dsmz_only, 10), collapse = ", "), "\n", sep = "")
    }
  }

  # Subset both data structures to common samples only.
  if (length(common_ids) < nrow(expr_mat)) {
    cat("  Using only ", length(common_ids), " common samples\n", sep = "")
    expr_mat_subset    <- expr_mat[common_ids, , drop = FALSE]
    dataset_vec_subset <- dataset_vec[common_ids]
  } else {
    expr_mat_subset    <- expr_mat
    dataset_vec_subset <- dataset_vec
  }

  # Align cluster vector to expression matrix row order.
  cluster_vec <- cluster_vec[rownames(expr_mat_subset)]
  
  # Verify complete alignment (no NA values after subsetting).
  if (any(is.na(cluster_vec))) {
    missing_ids <- rownames(expr_mat_subset)[is.na(cluster_vec)]
    stop("Internal error: Some samples still missing after subsetting for ", method_id,
         ". Examples: ", paste(head(missing_ids, 5), collapse = ", "))
  }

  # --------------------------------------------------------------------------
  # STEP 4: NEIGHBOURHOOD COMPUTATION
  # --------------------------------------------------------------------------
  # The core neighbourhood algorithm computes pairwise distances between all
  # cell lines and tumour samples, then identifies the nearest neighbours
  # for each cell line.

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Execute the neighbourhood computation with configured parameters.
  # PARAMETERS:
  #   emb_mat:   Expression matrix (samples × genes)
  #   cluster_m: Named vector of cluster assignments
  #   dataset:   Named vector of dataset origin labels
  #   top_frac:  Fraction of tumours to include as neighbours (10%)
  #   top_n_min: Minimum number of neighbours (ensures robustness)
  #   top_n_max: Maximum number of neighbours (prevents oversized neighbourhoods)
  #   method_id: Identifier for output labelling
  #   distance:  Distance metric (euclidean or correlation)
  res <- compute_tumour_neighbourhoods(
    emb_mat   = expr_mat_subset,
    cluster_m = cluster_vec,
    dataset   = dataset_vec_subset,
    top_frac  = 0.10,
    top_n_min = 30,
    top_n_max = 200,
    method_id = method_id,
    distance  = dist_type
  )

  # --------------------------------------------------------------------------
  # STEP 5: RESULT SERIALISATION
  # --------------------------------------------------------------------------
  # Results are saved in multiple formats:
  #   - RDS: Native R format for efficient programmatic access
  #   - CSV: Human-readable format for inspection and external tools

  nh_rds   <- file.path(outdir, paste0("Top_m_neighbourhoods_", method_id, ".rds"))
  long_rds <- file.path(outdir, paste0("Top_m_long_", method_id, ".rds"))
  long_csv <- file.path(outdir, paste0("Top_m_long_", method_id, ".csv"))

  saveRDS(res$neighbourhoods, nh_rds)
  saveRDS(res$long_df,       long_rds)
  write_csv(res$long_df,     long_csv)

  # --------------------------------------------------------------------------
  # STEP 6: RESULT SUMMARY LOGGING
  # --------------------------------------------------------------------------
  
  cat("\nTumour neighbourhoods computed successfully!\n")
  cat("Method      :", res$method_id, "\n")
  cat("Cell lines  :", length(res$neighbourhoods), "\n")
  cat("Total pairs :", nrow(res$long_df), "\n\n")

  # Display neighbourhood membership summary by cell line.
  print(
    res$long_df %>%
      count(cell_line, in_top) %>%
      arrange(desc(in_top), cell_line)
  )

  cat("\nAll results saved to:", outdir, "\n")
  
  # Return results invisibly (suppresses console output but allows assignment).
  invisible(res)
}

# ------------------------------------------------------------------------------
# SECTION 18: BATCH EXECUTION ACROSS ALL CLUSTERING METHODS
# ------------------------------------------------------------------------------
# The pipeline iterates over all discovered clustering methods, executing
# neighbourhood computation for each. Results are aggregated into a list
# for potential downstream analysis.

# Pre-allocate result storage for efficiency.
all_results <- vector("list", nrow(methods_exist))
names(all_results) <- methods_exist$method_id

# Initialise counters for execution summary.
n_success <- 0
n_skipped <- 0

# Execute neighbourhood computation for each method.
for (i in seq_len(nrow(methods_exist))) {
  m <- methods_exist[i, ]
  result <- run_single_neighbourhood(
    path      = m$path,
    method_id = m$method_id,
    outdir    = m$outdir
  )
  
  # Track execution outcomes.
  if (is.null(result)) {
    n_skipped <- n_skipped + 1
  } else {
    all_results[[m$method_id]] <- result
    n_success <- n_success + 1
  }
}

# Remove NULL entries from results list (failed methods).
# The Filter() function with Negate(is.null) retains only non-NULL elements.
all_results <- Filter(Negate(is.null), all_results)

# ------------------------------------------------------------------------------
# SECTION 19: EXECUTION SUMMARY
# ------------------------------------------------------------------------------
# A final summary reports the number of methods attempted, successful, and
# skipped, providing a quick overview of pipeline execution status.

cat("\n============================================================\n")
cat("SUMMARY:\n")
cat(sprintf("  Methods attempted: %d\n", nrow(methods_exist)))
cat(sprintf("  Successful: %d\n", n_success))
cat(sprintf("  Skipped (missing files): %d\n", n_skipped))

# Fail if no methods completed successfully.
if (n_success == 0) {
  stop("FATAL: No clustering methods completed successfully. Check cluster file paths.")
}

cat("\nAll available methods completed.\n")