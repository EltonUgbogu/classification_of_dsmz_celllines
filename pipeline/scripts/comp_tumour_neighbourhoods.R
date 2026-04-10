#!/usr/bin/env Rscript

# ==============================================================================
# comp_tumour_neighbourhoods.R
# Tumour Neighbourhood Computation Pipeline
# ==============================================================================
#
# PURPOSE:
# This script computes tumour neighbourhoods for each cell line in an integrated
# expression matrix. For each DSMZ cell line, the pipeline identifies the most
# similar tumour samples (non-DSMZ cohort) based on expression distance, enabling
# biological interpretation of cell line-tumour relationships within the context
# of consensus clustering results.
#
# USAGE:
# Rscript comp_tumour_neighbourhoods.R \
#   --config config/config.yaml \
#   --direction <pam50_euc|pam50_corr|HVG_euc|HVG_corr|MX_euc|...> \
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
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  make_option("--expr-rds", type = "character", default = NULL,
              help = "Optional override: integrated DSMZ+tumour expression RDS"),
  make_option("--mapping-rds", type = "character", default = NULL,
              help = "Optional override: cell_line → original_sample_id mapping RDS"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (e.g. HVG_euc, HVG_corr, MX_euc, pam50_euc)")
)

# parse_args() processes the command-line arguments according to the defined
# option_list and returns a named list containing the parsed values.
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Using config file:", opt$config, "\n")

# Load profiled config helper (merges defaults + profile)
# Determine script directory relative to config file
config_dir <- dirname(normalizePath(opt$config))
script_dir <- file.path(config_dir, "..", "scripts")
lib_config_path <- file.path(script_dir, "lib_config.R")
if (file.exists(lib_config_path)) {
  source(lib_config_path)
  # Use profiled config to merge defaults + profile (same as Snakemake)
  # Priority: opt$profile -> SNAKEMAKE_PROFILE env var -> config default
  # Note: %||% is defined in lib_config.R, so we can use it after sourcing
  profile_override <- if (!is.null(opt$profile)) opt$profile else Sys.getenv("SNAKEMAKE_PROFILE", unset = "")
  if (identical(profile_override, "")) profile_override <- NULL
  cfg <- read_profiled_config(opt$config, profile_override = profile_override)
} else {
  # Fallback: use direct yaml read (for backward compatibility)
  cfg <- yaml::read_yaml(opt$config)
  # If config has defaults, merge them
  if (!is.null(cfg$defaults)) {
    defaults <- cfg$defaults
    profile <- cfg$profile %||% names(cfg$profiles)[[1]] %||% "default"
    if (!is.null(cfg$profiles) && !is.null(cfg$profiles[[profile]])) {
      profcfg <- cfg$profiles[[profile]]
      # Simple merge: profile overrides defaults
      cfg <- list(
        paths = c(defaults$paths %||% list(), profcfg$paths %||% list()),
        analysis = c(defaults$analysis %||% list(), profcfg$analysis %||% list()),
        features = c(defaults$features %||% list(), profcfg$features %||% list()),
        tumour_neighbourhoods = c(defaults$tumour_neighbourhoods %||% list(), profcfg$tumour_neighbourhoods %||% list())
      )
    } else {
      cfg <- defaults
    }
  }
}

# --- Resolve relative paths against pipeline root (derived from config path) ---
# This ensures paths work correctly when Snakemake runs with --directory
pipe_root <- normalizePath(file.path(dirname(opt$config), ".."))

abs_from_root <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)  # already absolute
  file.path(pipe_root, p)
}

# Patch any config paths that might be relative
if (!is.null(cfg$paths)) {
  # Recursively absolutize all paths
  for (key in names(cfg$paths)) {
    if (!is.null(cfg$paths[[key]]) && is.character(cfg$paths[[key]])) {
      cfg$paths[[key]] <- abs_from_root(cfg$paths[[key]])
    }
  }
}

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
#   - Feature sets: pam50 (PAM50 gene signature), HVG (trend-corrected HVG), MX, etc.
#   - Distance metrics: euc (Euclidean), corr (correlation)

if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Validate direction (suffix-based, consistent with get_nh_methods)
if (!grepl("_(euc|corr)$", direction)) {
  stop("Invalid --direction: ", direction,
       "\nExpected a direction ending in _euc or _corr (e.g. HVG_euc, pam50_corr, MAD_euc, Spearman_corr).")
}

# ------------------------------------------------------------------------------
# SECTION 6: DIRECTION PARSING AND PATH CONFIGURATION
# ------------------------------------------------------------------------------
# The direction string is parsed to extract the gene set and distance metric,
# which determine input file paths and computational parameters.

# --- Derive profile-scoped unsup_root from the provided inputs ---
# The config may have a generic root, but we need the profile-specific one
# (e.g., results/unsupervised/brca, not just results/unsupervised)
guess_unsup_root_from_expr <- function(expr_rds) {
  # expr_rds like: 
  #   - results/unsupervised/brca/feature_selection_unsupervised/featuresets/Variance/expr_submatrix.rds
  #   - results/unsupervised/brca/feature_selection_unsupervised/featuresets/HVG/expr_submatrix.rds
  # Return: results/unsupervised/brca
  if (is.null(expr_rds) || !nzchar(expr_rds)) return(NULL)
  # Remove everything after /feature_selection_unsupervised/ or /tumour_neighbourhoods_input/
  result <- sub("/feature_selection_unsupervised/.*$", "", expr_rds)
  result <- sub("/tumour_neighbourhoods_input/.*$", "", result)
  result
}

guess_unsup_root_from_mapping <- function(mapping_rds) {
  # mapping_rds like: results/unsupervised/brca/tumour_neighbourhoods_input/cell_line_to_original_sample_id_Variance.rds
  # Return: results/unsupervised/brca
  if (is.null(mapping_rds) || !nzchar(mapping_rds)) return(NULL)
  sub("/tumour_neighbourhoods_input/.*$", "", mapping_rds)
}

# Try to derive from inputs first (most reliable)
# optparse keeps hyphenated option names, so access with backticks
unsup_root_expr <- guess_unsup_root_from_expr(opt$`expr-rds`)
unsup_root_map  <- guess_unsup_root_from_mapping(opt$`mapping-rds`)

# Use expr-derived if available, else mapping-derived, else fall back to config
unsup_root <- if (!is.null(unsup_root_expr) && nzchar(unsup_root_expr)) {
  unsup_root_expr
} else if (!is.null(unsup_root_map) && nzchar(unsup_root_map)) {
  unsup_root_map
} else {
  cfg$paths$unsup_root %||% "results/unsupervised"
}

cat("[INFO] Using unsup_root:", unsup_root, "\n")

# Extract feature method from direction (everything before last _euc/_corr).
direction_feature <- sub("_(euc|corr)$", "", direction)
# gene_set label: PAM50 for pam50, otherwise use the feature method name.
gene_set <- if (startsWith(direction, "pam50")) "PAM50" else toupper(direction_feature)

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
# optparse keeps hyphenated option names, so access with backticks
expr_mat_path <- opt$`expr-rds` %||%
  if (gene_set == "PAM50") {
    cfg$paths$tumour_nh_expr_pam50
  } else {
    stop("--expr-rds must be supplied for non-PAM50 directions (direction: ", direction, ")")
  }

# Resolve sample ID mapping path.
# The mapping file links technical sample identifiers to canonical cell line names.
mapping_path <- opt$`mapping-rds` %||% {
  # Infer mapping file location from expression matrix directory.
  map_suffix <- if (startsWith(direction, "pam50")) {
    "cell_line_to_original_sample_id_pam50.rds"
  } else {
    paste0("cell_line_to_original_sample_id_", direction_feature, ".rds")
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

# --- Robustness: accept either nh_paths or make_nh_paths conventions ---
# Handle cases where the file defines one but the script expects the other
if (!exists("make_nh_paths", inherits = TRUE) && exists("nh_paths", inherits = TRUE)) {
  # If nh_paths exists but make_nh_paths doesn't, create make_nh_paths as an alias
  make_nh_paths <- function(unsup_root) nh_paths
}
if (!exists("nh_paths", inherits = TRUE) && exists("make_nh_paths", inherits = TRUE)) {
  # If make_nh_paths exists but nh_paths doesn't, create nh_paths when make_nh_paths is called
  # This will be set below when make_nh_paths is called
}

# Verify that expected functions were successfully loaded.
# These safety checks catch configuration or file structure errors early.
if (!exists("make_nh_paths", inherits = TRUE) && !exists("nh_paths", inherits = TRUE)) {
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
# Sample identifiers from different sources (tumour cohorts, DSMZ) use inconsistent
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
  # These prefixes may be added during earlier pipeline stages.
  x <- sub("^(CELL:|TUMOUR:)", "", x)
  
  # Standardise delimiter characters to underscores.
  x <- gsub("-", "_", x)
  x <- gsub(" ", "_", x)
  x
}

# ------------------------------------------------------------------------------
# SECTION 9B: DSMZ ID PARSING FUNCTION
# ------------------------------------------------------------------------------
# Parses DSMZ technical IDs to extract canonical cell line and signature components.
# Expected format: NG-30919_RBL_15_lib626622_10102_4_1
#
# PARAMETERS:
#   x: A character vector of DSMZ technical identifiers
#
# RETURNS:
#   A list with components: ok, ng, cell, lib, lane, rep, sig

parse_dsmz_id <- function(x) {
  # Expected format: NG-30919_RBL_15_lib626622_10102_4_1
  # Or more generally: NG-XXXXX_CELL_LINE_libXXXXXX_XXXXX_XX_XX
  m <- regexec("^(NG-[0-9]+)_(.+)_(lib[0-9]+)_([0-9]+)_([0-9]+)_[0-9]+$", x)
  hit <- regmatches(x, m)[[1]]
  
  if (length(hit) == 0) {
    return(list(ok = FALSE, ng = NA_character_, cell = NA_character_, lib = NA_character_, 
                lane = NA_character_, rep = NA_character_,
                sig = substr(x, 1, 30)))
  }
  
  ng   <- hit[2]
  cell <- hit[3]  # Cell line (disease-agnostic: RBL_15, CAL_120, etc.)
  lib  <- hit[4]
  lane <- hit[5]
  rep  <- hit[6]
  
  list(
    ok = TRUE,
    ng = ng,
    cell = cell,
    lib = lib,
    lane = lane,
    rep = rep,
    sig = paste0(ng, "|", lib, "|", lane, "_", rep)
  )
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
# The integrated expression matrix contains both tumour samples (non-DSMZ cohort)
# and DSMZ cell line samples in a unified feature space. The matrix is stored
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
# First 10 lines of readRDS(mapping_path) as keys + values.
n_show <- min(10L, length(orig_to_cellline))
cat("[DEBUG] First ", n_show, " mapping entries (key → value):\n", sep = "")
for (j in seq_len(n_show)) {
  cat("  ", original_ids_raw[j], " → ", cell_line_names_raw[j], "\n", sep = "")
}

current_ids <- rownames(expr_mat)
cat("[DEBUG] Example expr_mat rownames:\n")
print(head(current_ids, 10))

# ------------------------------------------------------------------------------
# SECTION 13: DATASET ORIGIN ANNOTATION
# ------------------------------------------------------------------------------
# Each sample is annotated with its source dataset (DSMZ or TUMOUR) BEFORE
# any identifier renaming occurs. This order of operations is critical
# because the original identifiers are needed to correctly identify DSMZ
# samples via the mapping file.
#
# IMPORTANT: Dataset annotation must precede identifier normalisation to
# ensure accurate provenance tracking.
#
# Bind dsmz_mask to current_ids (technical IDs) for consistency.

# Create a logical mask identifying DSMZ samples by matching against mapping keys.
dsmz_mask <- current_ids %in% original_ids_raw
n_dsmz    <- sum(dsmz_mask)

if (n_dsmz == 0) {
  stop("FATAL: No DSMZ samples found in expr_mat rownames using mapping values.\n",
       "  Example rownames: ", paste(head(current_ids, 5), collapse = ", "),
       "\n  Example mapping values: ", paste(head(original_ids_raw, 5), collapse = ", "))
}

cat("[INFO] Detected ", n_dsmz, " DSMZ samples in expr_mat via mapping\n", sep = "")

# Assign dataset labels using vectorised conditional assignment.
# The ifelse() function evaluates element-wise, returning "DSMZ" for TRUE
# positions and "TUMOUR" for FALSE positions (cohort-agnostic).
dataset_vec <- ifelse(dsmz_mask, "DSMZ", "TUMOUR")
names(dataset_vec) <- current_ids

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
# SECTION 15: OPTION A - COLLAPSE DSMZ REPLICATES TO CANONICAL IDs
# ------------------------------------------------------------------------------
# Collapse technical DSMZ replicate IDs to canonical cell line IDs, and
# aggregate expression rows by mean. Tumour IDs remain unchanged.
#
# This produces outputs where DSMZ nodes are labeled as RBL_XX (canonical).

collapse_matrix_by_id <- function(mat, ids) {
  if (length(ids) != nrow(mat)) {
    stop("collapse_matrix_by_id: ids length does not match nrow(mat).")
  }
  sums <- rowsum(mat, group = ids, reorder = FALSE)
  counts <- as.integer(table(ids)[rownames(sums)])
  sums / counts
}

collapse_dataset_by_id <- function(vec, ids) {
  if (length(ids) != length(vec)) {
    stop("collapse_dataset_by_id: ids length does not match length(vec).")
  }
  groups <- split(vec, ids)
  vals <- vapply(names(groups), function(k) {
    u <- unique(groups[[k]])
    if (length(u) != 1L) {
      stop("Mixed dataset labels for collapsed ID '", k, "': ",
           paste(u, collapse = ", "))
    }
    u
  }, character(1))
  vals
}

collapse_cluster_vec_by_id <- function(vec) {
  if (is.null(names(vec))) {
    stop("collapse_cluster_vec_by_id: cluster_vec must be named.")
  }
  split_vec <- split(as.character(vec), names(vec))
  collapsed <- vapply(names(split_vec), function(k) {
    tab <- sort(table(split_vec[[k]]), decreasing = TRUE)
    if (length(tab) > 1 && tab[1] == tab[2]) {
      warning("Tied cluster assignments for ", k,
              " (", paste(names(tab)[tab == tab[1]], collapse = ", "),
              "); using first.")
    }
    names(tab)[1]
  }, character(1))
  collapsed
}

current_ids_raw <- rownames(expr_mat)

# Parse DSMZ IDs to extract components
parsed <- lapply(current_ids_raw, parse_dsmz_id)

# Extract canonical cell line names
# Prefer mapping if available; otherwise parse from technical ID
canonical_from_map <- orig_to_cell[current_ids_raw]  # NA for tumour samples
canonical_from_parse <- vapply(parsed, function(z) {
  if (isTRUE(z$ok)) z$cell else NA_character_
}, character(1))

# Choose canonical for DSMZ: mapping takes precedence, fallback to parse
cell_line_canonical_raw <- ifelse(
  dsmz_mask,
  ifelse(is.na(canonical_from_map), canonical_from_parse, canonical_from_map),
  NA_character_
)
names(cell_line_canonical_raw) <- current_ids_raw

# Validate that all DSMZ samples have canonical labels
if (any(dsmz_mask & is.na(cell_line_canonical_raw))) {
  bad <- current_ids_raw[dsmz_mask & is.na(cell_line_canonical_raw)]
  stop("FATAL: Cannot determine canonical cell line for DSMZ IDs (first 10): ",
       paste(head(bad, 10), collapse = ", "))
}

# Build collapsed sample IDs (DSMZ -> CELL:canonical to avoid tumour ID collisions, tumour -> original)
sample_id_collapsed <- ifelse(
  dsmz_mask,
  paste0("CELL:", cell_line_canonical_raw),
  current_ids_raw
)
names(sample_id_collapsed) <- current_ids_raw

# Guard: Tumour IDs should remain unique after collapse
tumour_dups <- sample_id_collapsed[!dsmz_mask & duplicated(sample_id_collapsed)]
if (length(tumour_dups) > 0) {
  stop("Unexpected duplicate tumour IDs after collapse (first 5): ",
       paste(head(unique(tumour_dups), 5), collapse = ", "))
}

# Collapse expression matrix and dataset labels
expr_mat <- collapse_matrix_by_id(expr_mat, sample_id_collapsed)
dataset_vec <- collapse_dataset_by_id(dataset_vec, sample_id_collapsed)
dataset_vec <- dataset_vec[rownames(expr_mat)]

current_ids <- rownames(expr_mat)
dsmz_mask <- dataset_vec == "DSMZ"

# Build mapping from technical IDs to collapsed IDs (raw + normalized)
id_to_collapsed <- setNames(sample_id_collapsed, current_ids_raw)
id_to_collapsed_norm <- setNames(sample_id_collapsed, normalize_id(current_ids_raw))

# Canonical and display maps for outputs (DSMZ: strip CELL: prefix for clean names)
cell_line_canonical <- ifelse(dsmz_mask, sub("^CELL:", "", current_ids), NA_character_)
names(cell_line_canonical) <- current_ids
cell_line_display <- ifelse(dsmz_mask, sub("^CELL:", "", current_ids), current_ids)
names(cell_line_display) <- current_ids

cat("[INFO] Option A: DSMZ replicates collapsed to canonical IDs.\n")
cat("[INFO] Total samples after collapse: ", nrow(expr_mat),
    " ( ", sum(dataset_vec == "DSMZ"), " DSMZ )\n", sep = "")
cat("[DEBUG] Example DSMZ labels (collapsed):\n")
print(head(data.frame(
  canon  = current_ids[dsmz_mask],
  stringsAsFactors = FALSE
), 6))

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
  # STEP 2: CLUSTER IDENTIFIER NORMALISATION (COLLAPSED IDs)
  # --------------------------------------------------------------------------
  # Cluster vector identifiers are technical IDs; map them to collapsed IDs
  # (DSMZ -> canonical, tumour unchanged). Use raw and normalized mapping and
  # keep whichever yields better overlap.

  nm_raw <- names(cluster_vec)

  mapped_raw <- id_to_collapsed[nm_raw]
  n_raw <- sum(!is.na(mapped_raw) & mapped_raw %in% current_ids)

  nm_norm <- normalize_id(nm_raw)
  mapped_norm <- id_to_collapsed_norm[nm_norm]
  n_norm <- sum(!is.na(mapped_norm) & mapped_norm %in% current_ids)

  if (n_norm > n_raw) {
    names(cluster_vec) <- mapped_norm
    mode <- "normalized"
    n_best <- n_norm
  } else {
    names(cluster_vec) <- mapped_raw
    mode <- "raw"
    n_best <- n_raw
  }

  # Drop entries that could not be mapped
  valid_name <- !is.na(names(cluster_vec)) & nzchar(names(cluster_vec))
  if (any(!valid_name)) {
    cluster_vec <- cluster_vec[valid_name]
  }

  # Collapse duplicate IDs by majority vote (replicates)
  if (anyDuplicated(names(cluster_vec)) > 0) {
    cluster_vec <- collapse_cluster_vec_by_id(cluster_vec)
  }

  cat("[INFO] Cluster ID matching mode for ", method_id, ": ", mode, "\n", sep = "")
  cat("[INFO] Matched ", n_best, " / ", length(nm_raw),
      " cluster IDs to collapsed expr_mat rownames\n", sep = "")

  if (n_best == 0) {
    stop("No cluster IDs match expr_mat rownames for ", method_id,
         "\nExamples cluster IDs: ", paste(head(nm_raw, 5), collapse = ", "),
         "\nExamples expr IDs: ", paste(head(current_ids, 5), collapse = ", "))
  }

  # Quick DSMZ-specific sanity: do we match DSMZ at all?
  matched_ids <- intersect(current_ids, names(cluster_vec))
  matched_dsmz <- sum(dsmz_mask[current_ids %in% matched_ids])
  cat("[INFO] Matched DSMZ samples: ", matched_dsmz, " / ", sum(dsmz_mask), "\n", sep = "")

  # Log DSMZ sample names for verification (cohort-agnostic: use dataset_vec, not pattern matching).
  # Note: cluster_vec names are collapsed IDs, so we need to check against current_ids dataset_vec
  matched_ids_in_clusters <- intersect(names(cluster_vec), current_ids)
  dsmz_names_in_clusters <- matched_ids_in_clusters[dataset_vec[matched_ids_in_clusters] == "DSMZ"]
  cat("[INFO] First DSMZ names in cluster_vec: ",
      paste(head(dsmz_names_in_clusters, 5), collapse = ", "), "\n")

  # --------------------------------------------------------------------------
  # STEP 3: SAMPLE ALIGNMENT BETWEEN EXPRESSION AND CLUSTER DATA
  # --------------------------------------------------------------------------
  # After collapsing, expr_mat and cluster_vec should share canonical IDs.
  # Only samples present in both are used for neighbourhood computation.

  common_ids <- intersect(current_ids, names(cluster_vec))

  if (length(common_ids) == 0) {
    stop("No overlapping sample IDs between expr_mat and cluster_vec for ", method_id,
         "\n  expr_mat sample IDs (first 5): ", paste(head(current_ids, 5), collapse = ", "),
         "\n  cluster_vec sample IDs (first 5): ", paste(head(names(cluster_vec), 5), collapse = ", "))
  }

  # Report samples present in expression matrix but missing from clusters.
  missing_in_clusters <- setdiff(current_ids, names(cluster_vec))
  if (length(missing_in_clusters) > 0) {
    cat("[WARNING] ", length(missing_in_clusters),
        " samples in expr_mat not found in cluster_vec for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_clusters, 5), collapse = ", "), "\n", sep = "")
  }

  # Report samples present in clusters but missing from expression matrix.
  missing_in_expr <- setdiff(names(cluster_vec), current_ids)
  if (length(missing_in_expr) > 0) {
    # Identify potential DSMZ samples by checking if they match canonical cell line names
    # (cohort-agnostic: uses mapping-based detection, not pattern matching)
    missing_norm <- normalize_id(missing_in_expr)
    potential_dsmz <- missing_in_expr[missing_norm %in% normalize_id(cell_line_names_raw)]
    cat("[WARNING] ", length(missing_in_expr),
        " samples in cluster_vec not found in expr_mat for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_expr, 5), collapse = ", "), "\n", sep = "")
    if (length(potential_dsmz) > 0) {
      cat("  Potential DSMZ samples in cluster_vec (not in expr_mat): ",
          paste(head(potential_dsmz, 10), collapse = ", "), "\n", sep = "")
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

  # Align cluster vector to expression matrix row order (collapsed IDs).
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
  # TARGET-style output: minimal schema (method, cell_line, tumor_id, rank, distance, in_top, cluster).

  # Ensure long_df has a cell identifier column
  if (!"cell_tech_id" %in% colnames(res$long_df)) {
    if ("cell_line" %in% colnames(res$long_df)) {
      res$long_df <- dplyr::rename(res$long_df, cell_tech_id = cell_line)
    } else {
      stop("res$long_df missing cell identifier column (expected cell_tech_id or cell_line).")
    }
  }

  # Sanity check: if orig_to_cell keys don't match cell_tech_id, canonical mapping will be NA.
  cat("[DEBUG] Example cell_tech_id values:\n")
  print(head(res$long_df$cell_tech_id, 5))
  cat("[DEBUG] Example mapped canonical values (orig_to_cell):\n")
  print(head(orig_to_cell[head(res$long_df$cell_tech_id, 5)], 5))

  # TARGET-style final table: single cell_line column (use orig_to_cell when keys match, else strip CELL: for DSMZ, else tech ID).
  out_long <- res$long_df %>%
    dplyr::mutate(
      cell_line = dplyr::coalesce(orig_to_cell[cell_tech_id], sub("^CELL:", "", cell_tech_id), cell_tech_id)
    ) %>%
    dplyr::select(method, cell_line, tumor_id, rank, distance, in_top, cluster)

  nh_rds   <- file.path(outdir, paste0("Top_m_neighbourhoods_", method_id, ".rds"))
  long_rds <- file.path(outdir, paste0("Top_m_long_", method_id, ".rds"))
  long_csv <- file.path(outdir, paste0("Top_m_long_", method_id, ".csv"))

  saveRDS(res$neighbourhoods, nh_rds)
  saveRDS(out_long,          long_rds)
  write_csv(out_long,        long_csv)

  # --------------------------------------------------------------------------
  # STEP 6: RESULT SUMMARY LOGGING
  # --------------------------------------------------------------------------
  
  cat("\nTumour neighbourhoods computed successfully!\n")
  cat("Method      :", res$method_id, "\n")
  cat("Cell lines  :", length(res$neighbourhoods), "\n")
  cat("Total pairs :", nrow(out_long), "\n\n")

  # Display neighbourhood membership summary by cell line.
  cat("\nSummary by cell line ID:\n")
  print(
    out_long %>%
      dplyr::count(cell_line, in_top) %>%
      dplyr::arrange(dplyr::desc(in_top), cell_line)
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