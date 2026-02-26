#!/usr/bin/env Rscript

# tumour_neighbourhood_p_consensus.R

# ==============================================================================
# TUMOUR NEIGHBOURHOOD CONSENSUS COMPUTATION
# ==============================================================================
#
# SCRIPT OVERVIEW
# ---------------
# This script aggregates tumour neighbourhood assignments across multiple
# clustering methods to compute the consensus probability p_consensus(c, t)
# for each cell line-tumour pair. The p_consensus metric quantifies the
# robustness of the relationship between a cancer cell line (c) and a
# primary tumour sample (t) by measuring the proportion of clustering
# analyses that place them in the same molecular cluster.
#
# BIOLOGICAL CONTEXT
# ------------------
# Cancer cell lines serve as essential laboratory models for understanding
# tumour biology and testing therapeutic interventions. However, the degree
# to which any given cell line faithfully represents the transcriptomic
# landscape of primary tumours remains a critical question for translational
# research.
#
# The tumour neighbourhood concept addresses this by identifying, for each
# cell line, the set of primary tumours that exhibit similar gene expression
# profiles. A cell line's "neighbourhood" consists of tumours that cluster
# together with it, suggesting shared molecular characteristics.
#
# THE CONSENSUS PROBLEM
# ---------------------
# Any single clustering analysis is subject to methodological choices:
#   - Feature selection method (variance, MAD, entropy, PCA loadings, etc.)
#   - Distance metric (Euclidean, correlation-based)
#   - Clustering algorithm (hierarchical, k-means)
#   - Number of clusters (k)
#
# Different choices may yield different cluster assignments, raising the
# question: which tumours are "truly" in a cell line's neighbourhood?
#
# CONSENSUS SOLUTION
# ------------------
# The p_consensus metric resolves this ambiguity by aggregating results
# across multiple clustering configurations. For a cell line c and tumour t:
#
#   p_consensus(c, t) = (number of methods placing t in c's neighbourhood) /
#                       (total number of clustering methods evaluated)
#
# This value ranges from 0 (no methods agree) to 1 (all methods agree).
# High p_consensus indicates a robust relationship that persists regardless
# of analytical choices, while low values suggest method-dependent associations.
#
# INTERPRETATION GUIDELINES
# -------------------------
#   p_consensus >= 0.7: Strong consensus; the cell line-tumour relationship
#                       is robust across the majority of analytical methods.
#
#   p_consensus = 1.0:  Perfect consensus; all evaluated methods agree that
#                       the tumour belongs to the cell line's neighbourhood.
#
#   p_consensus < 0.5:  Weak consensus; the relationship is observed in fewer
#                       than half of the methods and may be method-specific.
#
# OPERATIONAL MODES
# -----------------
# The script supports two operational modes:
#
#   Per-disease mode: Reads pre-computed neighbourhood tables from individual
#   clustering analyses and aggregates them. This is the standard mode for
#   single-cancer-type analyses (e.g., breast cancer only).
#
#   Pan-cancer mode: Directly processes clustering RDS files from pan-cancer
#   analyses where multiple cancer types are analysed jointly. This mode
#   is activated when --ccp_hc_rds and --ccp_kmeans_rds arguments are provided.
#
# OUTPUT FILES
# ------------
# The script generates:
#   1. Consensus table (TSV and RDS): p_consensus values for all cell line-tumour pairs
#   2. Distribution histogram: Visualises the distribution of p_consensus values
#   3. Per-cell-line summary plot: Shows maximum p_consensus for each cell line
#   4. Scatter plot: Displays relationship between max p_consensus and neighbourhood breadth
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The following R packages provide essential functionality for data manipulation,
# file I/O, and visualisation. The suppressPackageStartupMessages() wrapper
# prevents verbose startup messages from cluttering log files.

suppressPackageStartupMessages({
  
  # dplyr: Provides a grammar of data manipulation with intuitive verbs
  # (filter, select, mutate, summarise, arrange) for transforming data frames.
  # Essential for the aggregation and summarisation operations in this script.
  library(dplyr)
  
  # readr: Implements fast and consistent functions for reading and writing
  # delimited files. Provides better defaults than base R for handling

  # character encodings and column type inference.
  library(readr)
  
  # purrr: Provides functional programming tools for working with functions
  # and vectors. The map() family of functions enables elegant iteration
  # over lists and data frame rows.
  library(purrr)
  
  # tidyr: Provides tools for reshaping data between wide and long formats.
  # The unnest() function is used to expand nested data frames.
  library(tidyr)
  
  # ggplot2: The grammar of graphics implementation for R. Provides a
  # consistent framework for building layered statistical visualisations.
  library(ggplot2)
  
  # ggrepel: Extends ggplot2 with text labels that automatically avoid
  # overlapping. Essential for creating readable scatter plots with
  # many labelled points.
  library(ggrepel)
  
  # yaml: Provides functions for reading YAML configuration files.
  # YAML is used throughout the pipeline for parameter specification.
  library(yaml)
})

# ------------------------------------------------------------------------------
# NULL-COALESCING OPERATOR
# ------------------------------------------------------------------------------
# The %||% operator returns the left operand if it is not NULL, otherwise
# returns the right operand. This pattern simplifies handling of optional
# values and default assignments throughout the script.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# SECTION 2: COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
# The optparse package enables systematic parsing of command-line arguments,
# allowing the script to be called from automated pipelines with different
# parameter combinations without modifying the source code.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  
  # --config: Path to the YAML configuration file containing analysis parameters.
  # This centralises settings including file paths and algorithm parameters.
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  
  # --profile: The configuration profile to use (e.g., 'brca', 'nbl', 'rbl').
  # Profiles enable disease-specific parameter sets within a single config file.
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  
  # --direction: Specifies the feature set and distance metric combination.
  # Format: {feature_set}_{distance_metric}
  # Examples: 'hvg_euc', 'pam50_corr', 'Entropy_euclidean', 'Variance_cosine'
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr, etc.)"),
  
  # --workdir: Pipeline working directory for anchoring relative paths.
  # This ensures correct file resolution regardless of invocation context.
  make_option("--workdir", type = "character", default = ".",
              help = "Pipeline working directory [default: current directory]"),
  
  # Pan-cancer mode arguments (activated when --ccp_hc_rds is provided):
  
  # --ccp_hc_rds: Path to hierarchical clustering results from ConsensusClusterPlus.
  # When provided, enables pan-cancer mode for cross-disease analysis.
  make_option("--ccp_hc_rds", type = "character", default = NULL,
              help = "Pan-cancer: CCP HC clustering RDS (optional)"),
  
  # --ccp_kmeans_rds: Path to k-means clustering results from ConsensusClusterPlus.
  make_option("--ccp_kmeans_rds", type = "character", default = NULL,
              help = "Pan-cancer: CCP KMeans clustering RDS (optional)"),
  
  # --meta_tsv: Path to joint metadata file containing sample annotations.
  # Required in pan-cancer mode to distinguish cell lines from tumours.
  make_option("--meta_tsv", type = "character", default = NULL,
              help = "Pan-cancer: Joint metadata TSV (required if --ccp_hc_rds provided)"),
  
  # --out_rds: Output path for the RDS file containing consensus results.
  make_option("--out_rds", type = "character", default = NULL,
              help = "Pan-cancer: Output RDS path (required if --ccp_hc_rds provided)"),
  
  # --out_tsv: Output path for the TSV file containing consensus results.
  make_option("--out_tsv", type = "character", default = NULL,
              help = "Pan-cancer: Output TSV path (required if --ccp_hc_rds provided)")
)

# Parse command-line arguments according to the defined option specifications.
opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------------------------
# SECTION 3: PATH RESOLUTION AND CONFIGURATION LOADING
# ------------------------------------------------------------------------------
# File paths in command-line arguments or configuration files may be specified
# as relative paths. This section establishes the working directory context
# and provides utilities for resolving paths to absolute form.

# Normalise the working directory path.
# The normalizePath() function resolves symbolic links and relative components
# to produce a canonical absolute path.
workdir <- normalizePath(opt$workdir, winslash = "/", mustWork = TRUE)
cat("[INFO] getwd(): ", getwd(), "\n", sep = "")
cat("[INFO] workdir: ", workdir, "\n", sep = "")

abs_path <- function(p) {
  # abs_path(): Converts a potentially relative path to an absolute path.
  #
  # If the path is already absolute (starts with '/' on Unix or a drive letter
  # on Windows), it is returned unchanged. Otherwise, it is resolved relative
  # to the pipeline working directory.
  #
  # Parameters:
  #   p: A file system path (character string)
  #
  # Returns:
  #   The absolute path (character string), or the input unchanged if NULL/empty
  
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(workdir, p)
}

# Load shared plotting helper functions.
# These functions provide consistent visual styling across pipeline outputs.
# Note: unsup_plot_helpers.R contains plotting utilities (theme_dimred, etc.)
# and shared helpers used across the pipeline. It's sourced here for the
# plotting theme functions used throughout this script.
pipe_root <- normalizePath(file.path(dirname(opt$config), ".."), mustWork = TRUE)
source(file.path(pipe_root, "R", "unsup_plot_helpers.R"))

# ------------------------------------------------------------------------------
# DIRECTION VALIDATION
# ------------------------------------------------------------------------------
# The direction parameter encodes the feature selection method and distance
# metric. This validation ensures the direction follows expected naming
# conventions before proceeding with analysis.

if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Validate direction suffix format.
# Accepted suffixes indicate the distance metric:
#   _euc / _euclidean: Euclidean distance
#   _corr / _correlation: Correlation-based distance
#   _cos / _cosine: Cosine distance
if (!grepl("_(euc|corr|euclidean|cosine)$", direction)) {
  stop("Invalid --direction: ", direction,
       "\nExpected a direction ending in one of: _euc, _corr, _euclidean, _cosine.",
       "\nExamples: hvg_euc, pam50_corr, Entropy_euclidean, Variance_cosine.")
}

# ------------------------------------------------------------------------------
# CONFIGURATION LOADING
# ------------------------------------------------------------------------------
# The configuration file centralises analysis parameters. The profiled
# configuration system enables disease-specific parameter overrides while
# inheriting common default values.

config_dir <- dirname(normalizePath(opt$config))
script_dir <- file.path(config_dir, "..", "scripts")
lib_config_path <- file.path(script_dir, "lib_config.R")

if (file.exists(lib_config_path)) {
  source(lib_config_path)
  
  # Determine the profile to use with priority:
  # 1. Command-line --profile argument
  # 2. SNAKEMAKE_PROFILE environment variable
  # 3. Default profile from configuration
  profile_override <- if (!is.null(opt$profile)) opt$profile else Sys.getenv("SNAKEMAKE_PROFILE", unset = "")
  if (identical(profile_override, "")) profile_override <- NULL
  
  cfg <- read_profiled_config(opt$config, profile_override = profile_override)
} else {
  # Fallback: direct YAML parsing for backward compatibility
  cfg <- yaml::read_yaml(opt$config)
}

cat("[INFO] Using config file: ", opt$config, "\n", sep = "")
cat("[INFO] Direction: ", direction, "\n", sep = "")

# Extract gene set label for plot annotations.
# The gene set is inferred from the direction prefix.
gene_set <- if (startsWith(direction, "pam50")) "PAM50" else "HVG500"
cat("[INFO] Gene set: ", gene_set, "\n", sep = "")

# Validate that the unsupervised analysis root directory is configured.
if (is.null(cfg$paths$unsup_root) || !nzchar(cfg$paths$unsup_root)) {
  profile_used <- cfg$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "not set")
  stop("cfg$paths$unsup_root is empty — check profile merge and SNAKEMAKE_PROFILE.\n",
       "  Profile used: ", profile_used, "\n",
       "  Config file: ", opt$config)
}

# Resolve paths to absolute form.
unsup_root <- cfg$paths$unsup_root
unsup_root <- abs_path(unsup_root)
base_dir   <- file.path(unsup_root, "tumour_neighbourhoods", direction)

cat("[INFO] unsup_root: ", unsup_root, "\n", sep = "")
cat("[INFO] tumour_neighbourhoods base_dir: ", base_dir, "\n", sep = "")

# ==============================================================================
# SECTION 4: PAN-CANCER MODE
# ==============================================================================
# Pan-cancer mode is activated when clustering RDS files are provided directly.
# This mode supports cross-disease analyses where tumours and cell lines from
# multiple cancer types are analysed jointly.
#
# In pan-cancer mode, the script:
# 1. Loads pre-computed clustering results from ConsensusClusterPlus
# 2. Identifies cell lines and tumours using metadata annotations
# 3. Computes neighbourhoods by finding tumours in the same cluster as each cell line
# 4. Aggregates across clustering methods to compute p_consensus

if (!is.null(opt$ccp_hc_rds) && !is.null(opt$ccp_kmeans_rds)) {
  cat("[INFO] Pan-cancer mode: computing neighbourhoods from clustering RDS\n")
  
  # Validate required arguments for pan-cancer mode.
  if (is.null(opt$meta_tsv) || is.null(opt$out_rds) || is.null(opt$out_tsv)) {
    stop("Pan-cancer mode requires --meta_tsv, --out_rds, and --out_tsv")
  }
  
  # ---------------------------------------------------------------------------
  # SAMPLE ID NORMALISATION
  # ---------------------------------------------------------------------------
  # Sample identifiers may carry prefixes indicating their type (CELL:, TUMOUR:).
  # These prefixes are stripped to enable matching across different data sources.
  
  strip_prefix <- function(x) {
    # strip_prefix(): Removes type prefixes from sample identifiers.
    #
    # Prefixes such as "CELL:", "TUMOUR:", or "TUMOR:" are used in some
    # data sources to distinguish sample types. This function removes them
    # to enable consistent identifier matching.
    #
    # Parameters:
    #   x: A character vector of sample identifiers
    #
    # Returns:
    #   The identifiers with prefixes removed
    
    sub("^(CELL:|TUMOUR:|TUMOR:)", "", x)
  }
  
  # Load and normalise metadata.
  meta <- readr::read_tsv(opt$meta_tsv, show_col_types = FALSE) %>%
    mutate(sample_id = strip_prefix(as.character(sample_id)))
  
  # ---------------------------------------------------------------------------
  # CLUSTER EXTRACTION UTILITY
  # ---------------------------------------------------------------------------
  # Clustering results may be stored in various formats depending on the
  # upstream analysis tool. This function handles multiple common formats.
  
  extract_clusters <- function(obj, label) {
    # extract_clusters(): Extracts cluster assignments from various object formats.
    #
    # Different clustering tools store results in different structures. This
    # function attempts to locate and extract the cluster assignment vector
    # from several common formats.
    #
    # Parameters:
    #   obj:   The clustering result object (atomic vector or list)
    #   label: A descriptive label for error messages
    #
    # Returns:
    #   A named vector where names are sample IDs and values are cluster labels
    
    cand <- NULL
    
    if (is.atomic(obj) && !is.null(names(obj))) {
      # Direct named vector format
      cand <- obj
    } else if (is.list(obj)) {
      # Search common slot names in list structures
      for (nm in c("clusters", "clusters_labeled", "final_clusters", 
                   "consensusClass", "cluster", "labels")) {
        if (!is.null(obj[[nm]])) {
          cand <- obj[[nm]]
          break
        }
      }
    }
    
    if (is.null(cand)) stop("Could not extract clusters from ", label)
    if (is.null(names(cand))) stop("Extracted clusters from ", label, " have no names")
    cand
  }
  
  # ---------------------------------------------------------------------------
  # LOAD CLUSTERING RESULTS
  # ---------------------------------------------------------------------------
  # Two clustering methods are loaded: hierarchical clustering (HC) and
  # k-means (KM). Each provides an independent view of sample groupings.
  
  hc_obj <- readRDS(opt$ccp_hc_rds)
  km_obj <- readRDS(opt$ccp_kmeans_rds)
  
  # Extract cluster assignments as named character vectors.
  hc_clusters <- hc_obj$clusters
  km_clusters <- km_obj$clusters
  
  if (is.null(hc_clusters) || is.null(names(hc_clusters))) {
    stop("HC clusters missing or unnamed. Expected named character vector.")
  }
  if (is.null(km_clusters) || is.null(names(km_clusters))) {
    stop("KMeans clusters missing or unnamed. Expected named character vector.")
  }
  
  # ---------------------------------------------------------------------------
  # SAMPLE TYPE INFERENCE
  # ---------------------------------------------------------------------------
  # Distinguish cell lines from tumours using metadata annotations.
  # This is essential for computing neighbourhoods (which map cell lines to tumours).
  
  pick_col <- function(df, candidates) {
    # pick_col(): Selects the first matching column name from candidates.
    #
    # Parameters:
    #   df:         A data frame
    #   candidates: A character vector of potential column names
    #
    # Returns:
    #   The first matching column name, or NULL if none match
    
    hit <- candidates[candidates %in% colnames(df)][1]
    if (is.na(hit)) return(NULL)
    hit
  }
  
  # Identify the sample ID column.
  id_col <- pick_col(meta, c("sample_id", "sample", "id", "barcode", "sample_name", "cell_line"))
  if (is.null(id_col)) stop("Could not find sample ID column in metadata")
  
  infer_is_tumour <- function(meta) {
    # infer_is_tumour(): Determines which samples are tumours vs cell lines.
    #
    # Examines metadata columns for sample type annotations and returns
    # a logical vector indicating tumour status.
    #
    # Parameters:
    #   meta: A data frame containing sample metadata
    #
    # Returns:
    #   A logical vector (TRUE for tumours, FALSE for cell lines),
    #   or NULL if type cannot be inferred
    
    type_col <- pick_col(meta, c("sample_type", "type", "dataset_type"))
    if (!is.null(type_col)) {
      x <- tolower(trimws(as.character(meta[[type_col]])))
      return(x %in% c("tumour", "tumor", "tcga", "primary"))
    }
    NULL
  }
  
  is_tumour <- infer_is_tumour(meta)
  if (is.null(is_tumour)) stop("Cannot infer tumour vs cell from metadata")
  
  # Partition sample IDs by type.
  meta_ids <- as.character(meta[[id_col]])
  cell_ids <- meta_ids[!is_tumour]
  tumour_ids <- meta_ids[is_tumour]
  
  # ---------------------------------------------------------------------------
  # NEIGHBOURHOOD COMPUTATION
  # ---------------------------------------------------------------------------
  # For each clustering method, a tumour is in a cell line's neighbourhood
  # if they are assigned to the same cluster. This function identifies all
  # such cell line-tumour pairs.
  
  compute_neighbourhoods <- function(clusters, method_name) {
    # compute_neighbourhoods(): Identifies tumours in each cell line's neighbourhood.
    #
    # A tumour belongs to a cell line's neighbourhood if they share the same
    # cluster assignment. This function performs a join operation to find
    # all such matching pairs.
    #
    # Parameters:
    #   clusters:    A named vector of cluster assignments (names = sample IDs)
    #   method_name: A descriptive label for the clustering method
    #
    # Returns:
    #   A data frame with columns: cell_line, tumor_id, method, in_top
    
    # Create data frame with normalised sample IDs.
    clust_df <- data.frame(
      sample_id = strip_prefix(names(clusters)),
      cluster = as.character(clusters),
      stringsAsFactors = FALSE
    )
    
    # Separate cell lines and tumours.
    cell_clust <- clust_df %>% 
      filter(sample_id %in% cell_ids) %>%
      select(sample_id, cluster)
    
    tumour_clust <- clust_df %>% 
      filter(sample_id %in% tumour_ids) %>%
      select(sample_id, cluster)
    
    # Join on cluster: tumours in the same cluster as a cell line are neighbours.
    # The many-to-many relationship is expected since multiple cell lines and
    # tumours may share a cluster.
    cell_clust %>%
      inner_join(tumour_clust, by = "cluster", 
                 suffix = c("_cell", "_tumour"), 
                 relationship = "many-to-many") %>%
      select(cell_line = sample_id_cell, tumor_id = sample_id_tumour, cluster) %>%
      mutate(method = method_name, in_top = TRUE) %>%
      select(cell_line, tumor_id, method, in_top)
  }
  
  # Compute neighbourhoods for each clustering method.
  nh_hc <- compute_neighbourhoods(hc_clusters, "ccp_hc_expr_cell_tumour")
  nh_km <- compute_neighbourhoods(km_clusters, "ccp_kmeans_expr_cell_tumour")
  
  # Combine results from all methods.
  all_long <- dplyr::bind_rows(nh_hc, nh_km)
  
  # Add technical ID column for aggregation.
  all_long <- all_long %>%
    mutate(cell_tech_id = cell_line)
  
  n_methods_total <- length(unique(all_long$method))
  
  # ---------------------------------------------------------------------------
  # P_CONSENSUS COMPUTATION
  # ---------------------------------------------------------------------------
  # Aggregate across methods to compute the consensus probability.
  
  consensus_pairs <- all_long %>%
    filter(in_top) %>%
    count(cell_tech_id, tumor_id, name = "n_methods") %>%
    mutate(p_consensus = n_methods / n_methods_total) %>%
    arrange(cell_tech_id, desc(p_consensus))
  
  # Add display labels if available.
  if ("cell_line_canonical" %in% colnames(all_long)) {
    labels <- all_long %>%
      select(cell_tech_id, cell_line_canonical, cell_line_display) %>%
      distinct()
    consensus_pairs <- consensus_pairs %>% 
      left_join(labels, by = "cell_tech_id")
  }
  
  # Write outputs.
  dir.create(dirname(opt$out_rds), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(consensus_pairs, opt$out_tsv)
  saveRDS(consensus_pairs, opt$out_rds)
  
  cat("[SUCCESS] Pan-cancer p_consensus computed and written\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# SECTION 5: PER-DISEASE MODE - METHOD DISCOVERY
# ==============================================================================
# In per-disease mode, the script discovers available clustering methods by
# scanning the directory structure for pre-computed neighbourhood tables.
# This dynamic discovery approach ensures that all available methods contribute
# to the consensus calculation without requiring explicit configuration.
#
# The expected directory structure is:
#   {unsup_root}/tumour_neighbourhoods/{direction}/{method_id}/Top_m_long_*.csv

# List all method subdirectories.
method_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
method_dirs <- method_dirs[method_dirs != "final_consensus"]  # Exclude output directory

# Build a table of discovered methods and their associated files.
methods_tbl <- tibble::tibble(
  method_id = character(0),
  subdir    = character(0),
  file      = character(0),
  path      = character(0)
)

for (subdir in method_dirs) {
  method_dir <- file.path(base_dir, subdir)
  csv_files <- list.files(method_dir, pattern = "^Top_m_long_.*\\.csv$", full.names = FALSE)
  
  for (csv_file in csv_files) {
    method_id <- subdir
    full_path <- file.path(method_dir, csv_file)
    
    methods_tbl <- methods_tbl %>%
      dplyr::add_row(
        method_id = method_id,
        subdir    = subdir,
        file      = csv_file,
        path      = full_path
      )
  }
}

# Filter to only existing files.
methods_tbl <- methods_tbl %>%
  dplyr::mutate(exists = file.exists(path)) %>%
  dplyr::filter(exists)

n_methods_total <- nrow(methods_tbl)
cat("[INFO] Number of clustering methods contributing to p_consensus: ",
    n_methods_total, "\n", sep = "")

if (n_methods_total == 0) {
  stop("FATAL: No Top_m_long_*.csv files found in ", base_dir,
       "\nExpected files under: ", base_dir, "/*/Top_m_long_*.csv")
}

cat("[INFO] Methods discovered:\n")
print(methods_tbl %>% dplyr::select(method_id, subdir, file))

# ==============================================================================
# SECTION 6: NEIGHBOURHOOD DATA LOADING
# ==============================================================================
# Load neighbourhood tables from all discovered methods. Each table contains
# the tumour neighbourhood assignments for one clustering configuration.
#
# The tables are expected to contain:
#   - cell_line or cell_tech_id: Cell line identifier
#   - tumor_id: Tumour sample identifier
#   - in_top: Boolean indicating whether the tumour is in the neighbourhood
#   - Optional: rank, distance, display labels

neigh_list <- methods_tbl %>%
  mutate(
    data = map2(path, method_id, ~ {
      
      if (!file.exists(.x)) {
        stop("Neighbourhood file not found: ", .x)
      }
      
      read_csv(.x, show_col_types = FALSE) %>%
        mutate(method = .y) %>%
        # Ensure cell_tech_id exists as the primary key.
        # This handles variation in column naming across input files.
        {
          df <- .
          if ("cell_tech_id" %in% names(df)) {
            df %>% mutate(cell_tech_id = as.character(cell_tech_id))
          } else if ("cell_line" %in% names(df)) {
            df %>% mutate(cell_tech_id = as.character(cell_line))
          } else {
            stop("CSV file missing both cell_tech_id and cell_line columns: ", .x)
          }
        } %>%
        # Maintain backward compatibility with legacy column name.
        mutate(cell_line = cell_tech_id) %>%
        # Standardise column ordering.
        select(
          method, cell_tech_id, cell_line, tumor_id, in_top,
          any_of(c("rank", "distance", "cell_line_canonical", "cell_line_display")),
          everything()
        ) %>%
        # Remove potentially inconsistent cluster columns.
        select(-any_of(c("cluster", "cluster_id", "cluster_label", "subtype")))
    })
  )

# Combine all method results into a single long-format table.
all_long <- neigh_list %>%
  select(method_id, data) %>%
  unnest(data)

# Diagnostic output: verify data structure.
cat("=== Sanity check: counts per method / cell line / in_top ===\n")
all_long %>%
  count(method, cell_tech_id, in_top) %>%
  print(n = 12)

# ==============================================================================
# SECTION 7: P_CONSENSUS COMPUTATION
# ==============================================================================
# The p_consensus metric is computed as the proportion of clustering methods
# that identify a given tumour as being in a cell line's neighbourhood.
#
# FORMULA:
#   p_consensus(c, t) = |{methods where t in Top_m(c)}| / |{all methods}|
#
# where Top_m(c) denotes the set of top-m tumour neighbours for cell line c
# in a given clustering analysis.

consensus_pairs <- all_long %>%
  filter(in_top) %>%                           # Consider only neighbourhood members
  count(cell_tech_id, tumor_id, name = "n_methods") %>%  # Count contributing methods
  mutate(
    p_consensus = n_methods / n_methods_total  # Compute proportion
  ) %>%
  arrange(cell_tech_id, desc(p_consensus))     # Sort for readability

# Add display labels for visualisation if available.
if ("cell_line_canonical" %in% colnames(all_long)) {
  labels <- all_long %>%
    select(cell_tech_id, cell_line_canonical, cell_line_display) %>%
    distinct()
  consensus_pairs <- consensus_pairs %>% 
    left_join(labels, by = "cell_tech_id")
}

# ---------------------------------------------------------------------------
# Cohort-aware labels for plots (prevents "TCGA" wording when tumours are not TCGA)
# ---------------------------------------------------------------------------
# Prefer config overrides if you have them; otherwise infer from the data.
tumour_label <- cfg$labels$tumour %||% NA_character_
cell_label   <- cfg$labels$cell_line %||% "DSMZ cell lines"

if (is.na(tumour_label)) {
  # Infer from tumour IDs seen in results (works for TCGA and non-TCGA cohorts)
  tumour_ids_seen <- unique(consensus_pairs$tumor_id)
  
  if (length(tumour_ids_seen) == 0) {
    tumour_label <- "tumours"
  } else if (all(grepl("^TCGA-", tumour_ids_seen))) {
    tumour_label <- "TCGA tumours"
  } else if (any(grepl("^TCGA-", tumour_ids_seen))) {
    tumour_label <- "tumours (mixed cohorts)"
  } else {
    tumour_label <- "tumours (non-TCGA cohort)"
  }
}

# Optional: disease label (only for text, not logic)
disease_label <- cfg$labels$disease %||% "cell lines"

# Create output directory and write results.
cons_dir <- file.path(base_dir, "final_consensus")
dir.create(cons_dir, showWarnings = FALSE, recursive = TRUE)

cons_basename <- sprintf("Final_consensus_tumour_neighbourhoods_%s", direction)

cons_rds_path <- file.path(cons_dir, paste0(cons_basename, ".rds"))
cons_tsv_path <- file.path(cons_dir, paste0(cons_basename, ".tsv"))

write_tsv(consensus_pairs, cons_tsv_path)
saveRDS(consensus_pairs, cons_rds_path)

# Display summary statistics.
cat("\n=== Consensus summary ===\n")
consensus_pairs %>% glimpse()

cat("\n=== Cell lines with perfect consensus (p_consensus == 1) ===\n")
consensus_pairs %>%
  filter(p_consensus == 1) %>%
  count(cell_tech_id, sort = TRUE) %>%
  print()

# ==============================================================================
# SECTION 8: CONSENSUS DISTRIBUTION VISUALISATION
# ==============================================================================
# The histogram of p_consensus values provides insight into the overall
# robustness of cell line-tumour relationships. A distribution skewed toward
# high values indicates that most relationships are robust to methodological
# choices, while a uniform or low-skewed distribution suggests high sensitivity
# to analytical parameters.

# Filter invalid values before visualisation.
consensus_pairs <- consensus_pairs %>%
  filter(!is.na(p_consensus), p_consensus >= 0, p_consensus <= 1)

# Compute summary statistics for annotation.
n_pairs      <- nrow(consensus_pairs)
frac_ge_0_7  <- mean(consensus_pairs$p_consensus >= 0.7)
frac_eq_1    <- mean(consensus_pairs$p_consensus == 1)

cat("\n=== Summary of consensus strength ===\n")
cat(sprintf("Total pairs:        %d\n", n_pairs))
cat(sprintf("p_consensus >= 0.7 : %.1f%% (%d pairs)\n",
            100 * frac_ge_0_7,
            sum(consensus_pairs$p_consensus >= 0.7)))
cat(sprintf("p_consensus = 1.0 : %.1f%% (%d pairs)\n",
            100 * frac_eq_1,
            sum(consensus_pairs$p_consensus == 1)))

# Prepare annotation text for the plot.
annot_text <- sprintf(
  "Total pairs: %d\n>= 0.7: %.1f%%\n= 1.0: %.1f%%",
  n_pairs,
  100 * frac_ge_0_7,
  100 * frac_eq_1
)

# Create the histogram.
p_hist <- ggplot(consensus_pairs, aes(x = p_consensus)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),  # Convert to proportions
    binwidth = 0.05,
    colour   = "white"
  ) +
  geom_vline(xintercept = 0.7, linetype = "dashed") +  # Threshold reference line
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_dimred(base_size = 14) +
  labs(
    title    = sprintf(
      "Distribution of consensus strength across %d clustering methods (%s)",
      n_methods_total, gene_set
    ),
    subtitle = sprintf(
      "%d (cell line, tumour) pairs | %s gene set | %s",
      n_pairs, gene_set, tumour_label
    ),
    x        = expression(p[consensus]),
    y        = "Proportion of pairs"
  ) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  annotate(
    "text",
    x      = 0.98,
    y      = 0.95,
    label  = annot_text,
    hjust  = 1,
    vjust  = 1,
    size   = 3.5
  )

# Save the histogram.
plot_path_hist <- file.path(
  cons_dir,
  sprintf("Fig_p_consensus_distribution_%s.pdf", direction)
)
ggsave(plot_path_hist, p_hist, width = 8, height = 5)

cat("\nHistogram saved to:\n  ", plot_path_hist, "\n")

# ==============================================================================
# SECTION 9: PER-CELL-LINE SUMMARY
# ==============================================================================
# This summary provides a cell line-centric view of neighbourhood robustness.
# For each cell line, the following metrics are computed:
#   - n_pairs: Number of tumours in the cell line's neighbourhood (across any method)
#   - max_p: Maximum p_consensus value achieved with any tumour
#   - median_p: Median p_consensus across all tumour neighbours
#   - frac_ge_0_7: Fraction of neighbours with strong consensus (p >= 0.7)
#
# Cell lines with high max_p values have at least one tumour with which they
# show robust transcriptomic similarity. Those with high frac_ge_0_7 have
# broadly consistent neighbourhoods.

cell_summary <- consensus_pairs %>%
  group_by(cell_tech_id) %>%
  summarise(
    n_pairs      = n(),
    max_p        = max(p_consensus),
    median_p     = median(p_consensus),
    frac_ge_0_7  = mean(p_consensus >= 0.7),
    .groups      = "drop"
  ) %>%
  arrange(desc(max_p))

# Add display labels for plotting.
if ("cell_line_display" %in% colnames(consensus_pairs)) {
  labels <- consensus_pairs %>%
    select(cell_tech_id, cell_line_canonical, cell_line_display) %>%
    distinct()
  cell_summary <- cell_summary %>%
    left_join(labels, by = "cell_tech_id") %>%
    mutate(cell_line = ifelse(!is.na(cell_line_display), cell_line_display, cell_tech_id))
} else {
  cell_summary <- cell_summary %>%
    mutate(cell_line = cell_tech_id)
}

cat("\n=== Per-cell-line summary (top 10 by max p_consensus) ===\n")
cell_summary %>% head(10) %>% print()

# Create bar chart of per-cell-line maximum p_consensus.
p_cell <- ggplot(
  cell_summary,
  aes(x = reorder(cell_line, max_p), y = max_p,
      fill = max_p >= 0.7)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#bdbdbd"),
    labels = c("FALSE" = "< 0.7", "TRUE" = ">= 0.7"),
    name   = "Max p_consensus"
  ) +
  theme_dimred(base_size = 12) +
  labs(
    title    = sprintf(
      "Per-cell-line consensus with %s (%s)",
      tumour_label, direction
    ),
    subtitle = "Bars show max p_consensus across all tumour neighbours",
    x        = cell_label,
    y        = expression(max~p[consensus])
  ) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_path_cell <- file.path(
  cons_dir,
  sprintf("Fig_p_consensus_per_cell_line_%s.pdf", direction)
)
ggsave(plot_path_cell, p_cell, width = 8, height = 10)

cat("\nPer-cell-line summary plot saved to:\n  ", plot_path_cell, "\n")

# ==============================================================================
# SECTION 10: SCATTER PLOT - ANCHORING STRENGTH ANALYSIS
# ==============================================================================
# This visualisation examines the relationship between two aspects of cell
# line representativeness:
#
#   X-axis (max_p): The strongest consensus achieved with any single tumour.
#   This measures whether the cell line has at least one "anchor" tumour
#   with which it consistently clusters.
#
#   Y-axis (frac_ge_0_7): The proportion of neighbourhood tumours showing
#   strong consensus. This measures the breadth of robust representation.
#
# INTERPRETATION:
#   - High max_p, high frac_ge_0_7: Well-anchored cell line with broad tumour
#     representation. These are ideal models for the represented tumour subtype.
#
#   - High max_p, low frac_ge_0_7: Cell line has one or few strong tumour
#     matches but many weak associations. May represent a specific tumour niche.
#
#   - Low max_p, any frac_ge_0_7: Cell line lacks robust tumour anchors.
#     May have diverged from primary tumour biology during culture.

cell_summary2 <- cell_summary %>%
  mutate(
    frac_ge_0_7 = frac_ge_0_7,
    highlight   = (max_p >= 0.9 | frac_ge_0_7 >= 0.8)  # Label outstanding cell lines
  )

# Print highlighted cell lines for reference.
cell_summary2 %>%
  filter(highlight) %>%
  arrange(desc(max_p)) %>%
  select(cell_line, max_p, frac_ge_0_7) %>%
  print(n = 30)

# Create the scatter plot.
p_scatter <- ggplot(
  cell_summary2,
  aes(x = max_p, y = frac_ge_0_7)
) +
  geom_point(
    aes(size = n_pairs, fill = max_p >= 0.7),
    shape = 21,
    colour = "black",
    alpha  = 0.8
  ) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +   # Reference: 50% strong neighbours
  geom_vline(xintercept = 0.7, linetype = "dashed") +   # Reference: strong consensus threshold
  geom_text_repel(
    data = subset(cell_summary2, highlight),
    aes(label = cell_line),
    size              = 3.5,
    max.overlaps      = Inf,
    box.padding       = 0.4,
    point.padding     = 0.3,
    segment.size      = 0.3,
    min.segment.length = 0
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::percent_format(accuracy = 0)
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#bdbdbd"),
    labels = c("FALSE" = "< 0.7", "TRUE" = ">= 0.7"),
    name   = "Max p_consensus"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name  = "Number of tumour neighbours"
  ) +
  theme_dimred(base_size = 14) +
  labs(
    title    = sprintf(
      "Anchoring strength of %s (%s)",
      cell_label, direction
    ),
    subtitle = expression(
      x == max~p[consensus]~" per line; " ~
      y == "fraction of neighbours with " ~ p[consensus] >= 0.7
    ),
    x = expression(max~p[consensus]),
    y = "Tumour neighbours with strong consensus (%)"
  ) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

plot_path_scatter <- file.path(
  cons_dir,
  sprintf("Fig_p_consensus_cell_scatter_%s.pdf", direction)
)
ggsave(plot_path_scatter, p_scatter, width = 7, height = 6)

cat("\nScatter plot (max_p vs frac_ge_0_7) saved to:\n  ", plot_path_scatter, "\n")
cat("\n[SUCCESS] p_consensus computation finished for gene set: ", gene_set, "\n", sep = "")