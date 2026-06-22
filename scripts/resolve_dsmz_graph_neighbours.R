#!/usr/bin/env Rscript

# =============================================================================
# resolve_dsmz_graph_neighbors.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script resolves stable neighbours from DSMZ–DSMZ similarity graphs by
# computing the intersection of two complementary selection criteria:
#
#   (1) The best overall-performing feature–distance direction (global best),
#       determined by ranking all directions on aggregate consensus metrics.
#
#   (2) The best-performing feature–distance direction for each individual
#       cell line (cell-line winner), capturing per-sample optimal behaviour.
#
# The rationale for this intersection approach is to identify neighbours that
# are robust across both global and local optimality criteria, thereby reducing
# spurious connections that may arise from a single ranking scheme.
#
# OUTPUT
# ------
# The script produces a TSV file containing one row per cell line with:
#   - cell_line:        Cell line identifier
#   - final_neighbors:  Semicolon-separated list of stable neighbours
#   - n_final:          Count of stable neighbours (0 for isolates)
#   - best_overall_dir: The global best direction used
#   - reps_used:        Description of the intersection rule applied
#
# Isolates (cell lines with no stable neighbours under the intersection rule)
# are retained in the output with n_final = 0 for completeness.
#
#
# USAGE
# -----
# Rscript resolve_dsmz_graph_neighbors.R \
#     --config config/config.yaml \
#     --profile brca \
#     --winners_tsv results/unsupervised/tumour_neighbourhoods/final_consensus_all/p_consensus_winners_by_frac_ge_thr.tsv
#
# DEPENDENCIES
# ------------
# R packages: dplyr, readr, purrr, yaml, optparse
# External:   lib_config.R (shared configuration utilities)
#
# AUTHOR
# ------
# Developed as part of the unsupervised clustering pipeline for multi-cohort
# tumour neighbourhood analysis.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------
# The script loads required packages with startup messages suppressed to
# maintain clean log output in batch processing environments.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(yaml)
  library(optparse)
})

# -----------------------------------------------------------------------------
# Shared Configuration Utilities
# -----------------------------------------------------------------------------
# The script sources lib_config.R, which provides:
#   - read_profiled_config(): Profile-aware YAML configuration parsing
#   - %||%: Null-coalescing operator for default value assignment
#
# The script directory is resolved dynamically to support execution from
# different working directories.

# Resolve script directory more robustly
# Try multiple methods to find the script location
script_dir <- NULL

# Method 1: Check if script was called with --file argument
script_args <- commandArgs()
file_arg_idx <- grep("^--file=", script_args)
if (length(file_arg_idx) > 0) {
  script_path <- sub("^--file=", "", script_args[file_arg_idx[1]])
  if (file.exists(script_path)) {
    script_dir <- dirname(normalizePath(script_path))
  }
}

# Method 2: Try to get from Rscript's file argument (without =)
if (is.null(script_dir)) {
  file_arg_idx <- grep("--file", script_args)
  if (length(file_arg_idx) > 0 && file_arg_idx[1] < length(script_args)) {
    script_path <- script_args[file_arg_idx[1] + 1]
    if (file.exists(script_path)) {
      script_dir <- dirname(normalizePath(script_path))
    }
  }
}

# Method 3: Fallback to current working directory
if (is.null(script_dir) || length(script_dir) == 0) {
  script_dir <- getwd()
}

lib_config_path <- file.path(script_dir, "lib_config.R")
if (file.exists(lib_config_path)) {
  source(lib_config_path)
} else {
  stop("Cannot find lib_config.R at: ", lib_config_path)
}

# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# The script accepts the following arguments:
#
# Required:
#   --winners_tsv        Path to per-cell-line best direction results
#
# Configuration:
#   --config             Path to config.yaml (default: config/config.yaml)
#   --profile            Profile name (brca, nbl, rbl, or multicohort_cancer)
#
# Optional overrides:
#   --direction_summary_tsv   Path to direction summary for global best selection
#   --best_overall_dir        Manual override for global best direction
#   --graph_root              Custom graph root directory
#   --output_tsv              Custom output path

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml"),
  make_option("--profile", type = "character", default = NULL,
              help = "Profile name (brca, nbl, rbl, or multicohort_cancer)"),

  make_option("--winners_tsv", type = "character", default = NULL,
              help = "Path to p_consensus_winners_by_frac_ge_thr.tsv (per-cell-line best direction)"),

  make_option("--direction_summary_tsv", type = "character", default = NULL,
              help = "Optional: Path to p_consensus_direction_summary.tsv (for selecting best overall direction). If not provided, it is auto-resolved from output directory."),

  make_option("--best_overall_dir", type = "character", default = NULL,
              help = "Optional: Manually specify best overall direction (overrides direction_summary_tsv)."),

  make_option("--graph_root", type = "character", default = NULL,
              help = "Custom graph root directory. For pan-cancer typically: results/.../graphs/cell_line_similarity/within"),

  make_option("--output_tsv", type = "character", default = NULL,
              help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -----------------------------------------------------------------------------
# Configuration Loading
# -----------------------------------------------------------------------------
# The configuration is loaded with profile awareness. The profile may be
# specified via:
#   1. Command-line argument (--profile)
#   2. Environment variable (SNAKEMAKE_PROFILE)
#   3. Default profile in config.yaml
#
# This hierarchy supports both interactive use and Snakemake integration.

profile_override <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", unset = NA)
if (is.na(profile_override)) profile_override <- NULL

cfg <- read_profiled_config(opt$config, profile_override = profile_override)
profile <- cfg$profile

# -----------------------------------------------------------------------------
# Path Resolution Utility
# -----------------------------------------------------------------------------
# This helper function resolves paths that may be specified as either:
#   - Absolute paths (returned as-is after normalisation)
#   - Relative paths (resolved relative to the config file directory)
#
# This approach ensures consistent behaviour regardless of the working
# directory from which the script is invoked.

resolve_path_if_needed <- function(p, config_path) {
  if (is.null(p) || !nzchar(p)) return(p)
  if (file.exists(p) || dir.exists(p)) return(normalizePath(p, mustWork = FALSE))
  config_dir <- dirname(normalizePath(config_path))
  normalizePath(file.path(config_dir, "..", p), mustWork = FALSE)
}

# -----------------------------------------------------------------------------
# Tumour Neighbourhood Root Resolution
# -----------------------------------------------------------------------------
# The graph root directory contains the DSMZ–DSMZ similarity graph edge files.
# Resolution follows this precedence:
#
#   1. Explicit --graph_root argument (highest priority)
#   2. Profile-specific default paths:
#      - multicohort_cancer: {outdir}/graphs/cell_line_similarity/within
#      - Other profiles: {unsup_root}/tumour_neighbourhoods
#
# The script validates that the resolved directory exists before proceeding.

if (!is.null(opt$graph_root) && nzchar(opt$graph_root)) {
  tumour_nh_root <- resolve_path_if_needed(opt$graph_root, opt$config)
  if (!dir.exists(tumour_nh_root)) {
    stop("Cannot find graph_root directory: ", opt$graph_root,
         " (resolved to: ", tumour_nh_root, ")")
  }
} else {
  if (profile == "multicohort_cancer") {
    # Pan-cancer profile: read from dedicated multicohort_cancer config section
    cfg_full <- yaml::read_yaml(opt$config)
    mc_cfg <- cfg_full$multicohort_cancer %||%
      stop("multicohort_cancer section not found in config")

    outdir <- mc_cfg$outdir %||% "results/multicohort_cancer_benchmark"
    outdir <- resolve_path_if_needed(outdir, opt$config)

    tumour_nh_root <- file.path(outdir, "graphs", "cell_line_similarity", "within")
  } else {
    # Single-cohort profiles: use unsup_root from paths section
    unsup_root <- cfg$paths$unsup_root %||%
      stop("Missing paths.unsup_root in config")
    unsup_root <- resolve_path_if_needed(unsup_root, opt$config)
    if (!dir.exists(unsup_root)) {
      stop("Cannot find unsup_root directory (resolved to: ", unsup_root, ")")
    }
    tumour_nh_root <- file.path(unsup_root, "tumour_neighbourhoods")
  }
}

if (!dir.exists(tumour_nh_root)) {
  stop("Tumour neighbourhood root not found: ", tumour_nh_root)
}

# -----------------------------------------------------------------------------
# Direction Definitions
# -----------------------------------------------------------------------------
# Directions represent combinations of feature selection methods and distance
# metrics used to construct similarity graphs. The available directions depend
# on the analysis profile:
#
# Pan-cancer:
#   - Feature methods: Variance, MAD, MeanAbsDev, Entropy, PCA, Spearman, MX, kTotal, HVG
#   - Distance metrics: euc, corr
#   - Resulting directions: e.g., "Variance_euc", "MAD_corr", "HVG_euc"
#
# Single-cohort (brca, nbl, rbl, heme):
#   - Directions: the 9-method × 2-distance grid declared in
#     tumour_neighbourhoods.directions in config.yaml, e.g.:
#     Variance_euc, Variance_corr, MAD_euc, ..., kTotal_euc, kTotal_corr,
#     HVG_euc, HVG_corr
#   - HVG uses top-3000 genes (LOESS mean-variance trend correction).
#   - PAM50 directions (pam50_euc, pam50_corr) are added only when
#     use_pam50 is enabled in the profile config.

if (profile == "multicohort_cancer") {
  cfg_full <- yaml::read_yaml(opt$config)
  mc_cfg <- cfg_full$multicohort_cancer %||%
    stop("multicohort_cancer section not found in config")
  feature_methods <- mc_cfg$feature_methods %||%
    c("Variance", "MAD", "MeanAbsDev", "Entropy", "PCA", "Spearman", "MX", "kTotal", "HVG")
  dist_metrics <- mc_cfg$dist_metrics %||%
    c("euc", "corr")
  directions <- paste0(rep(feature_methods, each = length(dist_metrics)), "_",
                       rep(dist_metrics, times = length(feature_methods)))
} else {
  # Fall back to the tumour_neighbourhoods.directions from config; if not set,
  # use the 9-method × 2-distance grid that all active single-cohort profiles
  # declare (including HVG_euc / HVG_corr as of the 9-method activation).
  directions <- cfg$tumour_neighbourhoods$directions %||%
    c("Variance_euc", "Variance_corr",
      "MAD_euc",      "MAD_corr",
      "MeanAbsDev_euc", "MeanAbsDev_corr",
      "Entropy_euc",  "Entropy_corr",
      "PCA_euc",      "PCA_corr",
      "Spearman_euc", "Spearman_corr",
      "MX_euc",       "MX_corr",
      "kTotal_euc",   "kTotal_corr",
      "HVG_euc",      "HVG_corr")
}

# -----------------------------------------------------------------------------
# Winners File Loading (Per-Cell-Line Best Directions)
# -----------------------------------------------------------------------------
# The winners file contains the best-performing direction for each cell line,
# as determined by the fraction of consensus clustering solutions meeting a
# quality threshold (frac_ge_thr).
#
# Required columns:
#   - cell_line:            Cell line identifier
#   - best_dir_frac_ge_thr: Best direction for this cell line
#
# The script constructs a named vector (winners_map) for efficient lookup
# during neighbour resolution.

if (is.null(opt$winners_tsv) || !file.exists(opt$winners_tsv)) {
  stop("winners_tsv must be provided and exist. Got: ",
       opt$winners_tsv %||% "<NULL>")
}
winners <- read_tsv(opt$winners_tsv, show_col_types = FALSE)

required_cols <- c("cell_line", "best_dir_frac_ge_thr")
missing_cols <- setdiff(required_cols, colnames(winners))
if (length(missing_cols) > 0) {
  stop("winners_tsv missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

canonical_cell_line_id <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y)] <- ""

  has_lib <- grepl("^NG[-_][^_]+_.+_lib", y)
  y[has_lib] <- sub("^NG[-_][^_]+_(.+?)_lib.*$", "\\1", y[has_lib])

  # Retinoblastoma profiles must remain biological group-level graph nodes.
  y <- sub("^(RBL_[0-9]+)_[0-9]+$", "\\1", y)
  y
}

split_direction_values <- function(x) {
  vals <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  vals[nzchar(vals) & !is.na(vals)]
}

collapse_direction_values <- function(x) {
  vals <- unique(unlist(lapply(x, split_direction_values), use.names = FALSE))
  paste(sort(vals), collapse = ";")
}

if (profile == "multicohort_cancer") {
  before_winner_nodes <- length(unique(winners$cell_line))
  winners <- winners %>%
    mutate(cell_line = canonical_cell_line_id(cell_line)) %>%
    filter(nzchar(cell_line)) %>%
    group_by(cell_line) %>%
    summarise(
      best_dir_frac_ge_thr = collapse_direction_values(best_dir_frac_ge_thr),
      .groups = "drop"
    )
  cat("Canonicalised multicohort winner nodes:",
      before_winner_nodes, "->", length(unique(winners$cell_line)), "\n")
}

winners_map <- setNames(winners$best_dir_frac_ge_thr, winners$cell_line)

# -----------------------------------------------------------------------------
# Neighbour Extraction Utility
# -----------------------------------------------------------------------------
# This function extracts all neighbours of a given node from an edge dataframe.
# The function handles two common edge file formats:
#
#   Format 1: cell_line1, cell_line2 (undirected edges)
#   Format 2: from, to (directed or undirected edges)
#
# For both formats, the function returns all unique nodes connected to the
# query node, treating edges as undirected.
#
# Parameters:
#   edges_df: Dataframe containing edge definitions
#   node:     Node identifier for which to retrieve neighbours
#
# Returns:
#   Character vector of unique neighbour identifiers (empty if none found)

get_neighbors <- function(edges_df, node) {
  if (is.null(edges_df) || nrow(edges_df) == 0) return(character())

  if (all(c("cell_line1", "cell_line2") %in% colnames(edges_df))) {
    neighbors <- c(
      edges_df$cell_line2[edges_df$cell_line1 == node],
      edges_df$cell_line1[edges_df$cell_line2 == node]
    )
  } else if (all(c("from", "to") %in% colnames(edges_df))) {
    neighbors <- c(
      edges_df$to[edges_df$from == node],
      edges_df$from[edges_df$to == node]
    )
  } else {
    stop("Edge file missing 'cell_line1'/'cell_line2' or 'from'/'to' columns")
  }

  unique(neighbors[!is.na(neighbors)])
}

canonicalise_edge_df <- function(edges_df, direction = "") {
  if (is.null(edges_df) || nrow(edges_df) == 0) return(edges_df)

  if (all(c("cell_line1", "cell_line2") %in% colnames(edges_df))) {
    c1 <- "cell_line1"
    c2 <- "cell_line2"
  } else if (all(c("from", "to") %in% colnames(edges_df))) {
    c1 <- "from"
    c2 <- "to"
  } else {
    return(edges_df)
  }

  before_n <- nrow(edges_df)
  out <- edges_df %>%
    mutate(
      .node1 = canonical_cell_line_id(.data[[c1]]),
      .node2 = canonical_cell_line_id(.data[[c2]]),
      .a = pmin(.node1, .node2),
      .b = pmax(.node1, .node2)
    ) %>%
    filter(nzchar(.a), nzchar(.b), .a != .b)

  if ("similarity" %in% colnames(out)) {
    out <- out %>%
      group_by(.a, .b) %>%
      summarise(similarity = mean(suppressWarnings(as.numeric(similarity)), na.rm = TRUE),
                .groups = "drop")
  } else {
    out <- out %>% distinct(.a, .b)
  }

  out[[c1]] <- out$.a
  out[[c2]] <- out$.b
  out <- out %>% select(-.a, -.b)

  if (nrow(out) != before_n) {
    cat("Canonicalised", direction, "edges:", before_n, "->", nrow(out), "\n")
  }
  out
}

# -----------------------------------------------------------------------------
# Edge File Loading
# -----------------------------------------------------------------------------
# The script loads edge files for each direction. Edge files contain pairwise
# similarity relationships between DSMZ cell lines.
#
# File naming conventions differ by profile:
#   - Pan-cancer: {direction}/cell_line_similarity_graph_edges_{direction}.tsv
#   - Single-cohort: {direction}/final_consensus/cell_line_similarity_graph_edges_{direction}.tsv
#
# Missing edge files generate warnings but do not halt execution, allowing
# partial results when some directions are unavailable.

edge_file_for_direction <- function(root, direction) {
  candidates <- c(
    file.path(root, direction, "final_consensus",
              paste0("cell_line_similarity_graph_edges_", direction, ".tsv")),
    file.path(root, direction,
              paste0("cell_line_similarity_graph_edges_", direction, ".tsv"))
  )
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) == 0 || is.na(hit)) {
    candidates[1]
  } else {
    hit
  }
}

# Discover available directions from actual edge files (more robust than config)
# This ensures we use all available directions, not just those in config
discovered_directions <- character()
if (dir.exists(tumour_nh_root)) {
  subdirs <- list.dirs(tumour_nh_root, full.names = FALSE, recursive = FALSE)
  for (subdir in subdirs) {
    edge_file <- edge_file_for_direction(tumour_nh_root, subdir)
    if (file.exists(edge_file)) {
      discovered_directions <- c(discovered_directions, subdir)
    }
  }
}

# Use discovered directions if available, otherwise fall back to config directions
if (length(discovered_directions) > 0) {
  directions <- discovered_directions
  cat("Discovered", length(directions), "directions from edge files\n")
} else {
  cat("No directions discovered from edge files, using config directions\n")
}

edges_list <- list()
for (d in directions) {
  edge_file <- edge_file_for_direction(tumour_nh_root, d)

  if (file.exists(edge_file)) {
    edges_list[[d]] <- read_tsv(edge_file, show_col_types = FALSE)
  } else {
    warning("Graph edge file not found: ", edge_file)
  }
}

if (length(edges_list) == 0) {
  stop("No graph edge files found under: ", tumour_nh_root)
}

if (profile == "multicohort_cancer") {
  edges_list <- imap(edges_list, canonicalise_edge_df)
}

# -----------------------------------------------------------------------------
# Best Overall Direction Selection
# -----------------------------------------------------------------------------
# The best overall direction is the direction that performs best across all
# cell lines according to aggregate consensus metrics. Selection criteria
# (in order of priority):
#
#   1. frac_ge_thr:       Fraction of solutions meeting quality threshold
#   2. median_p_consensus: Median consensus probability across solutions
#   3. mean_p_consensus:   Mean consensus probability (tiebreaker)
#
# The direction may be:
#   - Manually specified via --best_overall_dir (highest priority)
#   - Auto-selected from direction summary TSV (default behaviour)
#
# If neither is available, the script terminates with an error.

best_overall_dir <- NULL

if (!is.null(opt$best_overall_dir) && nzchar(opt$best_overall_dir)) {
  # Manual override takes precedence
  best_overall_dir <- opt$best_overall_dir
} else {
  # Auto-resolve direction summary path
  dir_sum_path <- opt$direction_summary_tsv

  if (is.null(dir_sum_path) || !nzchar(dir_sum_path)) {
    # Default: same directory as winners_tsv
    dir_sum_path <- file.path(dirname(opt$winners_tsv),
                              "p_consensus_direction_summary.tsv")
  }

  if (!file.exists(dir_sum_path)) {
    stop("Direction summary TSV not found. ",
         "Provide --direction_summary_tsv or --best_overall_dir. ",
         "Tried: ", dir_sum_path)
  }

  dir_sum <- read_tsv(dir_sum_path, show_col_types = FALSE)

  # Handle column name variations (backwards compatibility)
  if ("median_p_consensus" %in% colnames(dir_sum)) {
    dir_sum <- dir_sum %>% mutate(median_p_consensus = median_p_consensus)
  } else if ("median_p" %in% colnames(dir_sum)) {
    dir_sum <- dir_sum %>% mutate(median_p_consensus = median_p)
  } else {
    stop("Direction summary TSV missing median_p_consensus or median_p column")
  }
  
  if ("mean_p_consensus" %in% colnames(dir_sum)) {
    dir_sum <- dir_sum %>% mutate(mean_p_consensus = mean_p_consensus)
  } else if ("median_p" %in% colnames(dir_sum)) {
    # Use median_p as fallback for mean_p_consensus if not available
    dir_sum <- dir_sum %>% mutate(mean_p_consensus = median_p_consensus)
  } else {
    stop("Direction summary TSV missing mean_p_consensus or median_p column")
  }

  needed <- c("direction", "frac_ge_thr", "median_p_consensus", "mean_p_consensus")
  miss <- setdiff(needed, colnames(dir_sum))
  if (length(miss) > 0) {
    stop("Direction summary TSV missing required columns: ",
         paste(miss, collapse = ", "))
  }

  # Select best direction by hierarchical sorting
  best_overall_dir <- dir_sum %>%
    arrange(desc(frac_ge_thr),
            desc(median_p_consensus),
            desc(mean_p_consensus)) %>%
    slice(1) %>%
    pull(direction)
}

# Validate that the selected direction has available edge data
if (!best_overall_dir %in% names(edges_list)) {
  warning("Best overall direction '", best_overall_dir,
          "' not available in edges_list. Available directions: ",
          paste(names(edges_list), collapse = ", "))
}

cat("Best overall direction:", best_overall_dir, "\n")

# -----------------------------------------------------------------------------
# Main Neighbour Resolution Function
# -----------------------------------------------------------------------------
# This function implements the core intersection logic for neighbour resolution.
#
# For a given cell line A:
#   - Let B = best_overall_dir (global best direction)
#   - Let F(A) = winners_map[A] (cell-line-specific best direction)
#   - Final neighbours = neighbours_B(A) ∩ neighbours_F(A)(A)
#
# The intersection ensures that only neighbours appearing in BOTH the global
# best and cell-line-specific best graphs are retained, providing a robust
# estimate of true similarity relationships.
#
# Edge cases:
#   - If B is missing from edges_list, an empty set is returned
#   - If F(A) is missing or unavailable, an empty set is returned
#   - If either neighbour set is empty, the intersection is empty
#
# Parameters:
#   cell_line:        Cell line identifier to resolve
#   best_overall_dir: Global best direction identifier
#   winners_map:      Named vector mapping cell lines to their best directions
#   edges_list:       Named list of edge dataframes by direction
#
# Returns:
#   Character vector of stable neighbour identifiers

resolve_neighbors_main <- function(cell_line, best_overall_dir,
                                   winners_map, edges_list) {
  winner_dir <- winners_map[[cell_line]]

  # If winner is missing, the intersection cannot be computed
  if (is.null(winner_dir) || is.na(winner_dir) || !nzchar(winner_dir)) {
    return(list(
      neighbors = character(),
      best_overall = best_overall_dir,
      winner = NA_character_
    ))
  }

  # Handle ties: winner_dir may be semicolon-separated list of tied directions
  winner_dirs <- strsplit(winner_dir, ";", fixed = TRUE)[[1]]
  winner_dirs <- trimws(winner_dirs)
  winner_dirs <- winner_dirs[nzchar(winner_dirs)]  # Remove empty strings

  # If best_overall_dir is not available, intersection is empty
  if (!best_overall_dir %in% names(edges_list)) {
    return(list(
      neighbors = character(),
      best_overall = best_overall_dir,
      winner = winner_dir
    ))
  }

  # Retrieve neighbours from best_overall_dir
  nb_best <- get_neighbors(edges_list[[best_overall_dir]], cell_line)

  # For tied winner directions: take UNION of neighbors from all tied directions
  # Then intersect with best_overall_dir neighbors
  nb_win_all <- character()
  for (wd in winner_dirs) {
    if (wd %in% names(edges_list)) {
      nb_win_all <- unique(c(nb_win_all, get_neighbors(edges_list[[wd]], cell_line)))
    }
  }

  # Return intersection of best_overall and union of winner directions
  list(
    neighbors = intersect(nb_best, nb_win_all),
    best_overall = best_overall_dir,
    winner = winner_dir
  )
}

# -----------------------------------------------------------------------------
# Cell Line Universe Definition
# -----------------------------------------------------------------------------
# The universe of cell lines for neighbour resolution is defined as the union
# of two sources:
#
#   1. All cell lines appearing in the winners TSV (cell_line column)
#   2. All nodes appearing in any edge file (as either endpoint)
#
# This ensures comprehensive coverage even when some cell lines appear only
# in edge files (e.g., as neighbours) but lack direct winner assignments.

get_nodes_from_edges <- function(edges_df) {
  if (is.null(edges_df) || nrow(edges_df) == 0) return(character())
  if (all(c("cell_line1", "cell_line2") %in% colnames(edges_df))) {
    c(edges_df$cell_line1, edges_df$cell_line2)
  } else if (all(c("from", "to") %in% colnames(edges_df))) {
    c(edges_df$from, edges_df$to)
  } else {
    character()
  }
}

all_nodes_edges <- unique(unlist(map(edges_list, get_nodes_from_edges)))
all_nodes_edges <- unique(all_nodes_edges[!is.na(all_nodes_edges)])

all_cell_lines <- sort(unique(c(winners$cell_line, canonical_cell_line_id(all_nodes_edges))))
if (profile == "multicohort_cancer") {
  profile_level_nodes <- grep("^NG[-_]", all_cell_lines, value = TRUE)
  if (length(profile_level_nodes) > 0) {
    stop("Profile-level nodes remain after multicohort canonicalisation: ",
         paste(head(profile_level_nodes, 20), collapse = ", "))
  }
}
cat("Total cell lines (winners ∪ edges):", length(all_cell_lines), "\n")

# -----------------------------------------------------------------------------
# Neighbour Resolution Across All Cell Lines
# -----------------------------------------------------------------------------
# The script applies the intersection-based resolution rule to each cell line
# in the universe. Results are collected into a long-format dataframe with
# one row per (cell_line, neighbour) pair.
#
# Cell lines with no stable neighbours produce a single row with an empty
# neighbour field, ensuring all cell lines appear in the output.

resolved_neighbors <- map_dfr(all_cell_lines, function(cl) {
  res <- resolve_neighbors_main(cl, best_overall_dir, winners_map, edges_list)

  neigh <- res$neighbors
  reps_used <- paste0("best_overall=", res$best_overall, ";winner=", res$winner)

  if (length(neigh) == 0) {
    data.frame(
      cell_line = cl,
      neighbor = "",
      reps_used = reps_used,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      cell_line = cl,
      neighbor = neigh,
      reps_used = reps_used,
      stringsAsFactors = FALSE
    )
  }
})

# -----------------------------------------------------------------------------
# Output Formatting
# -----------------------------------------------------------------------------
# The long-format results are aggregated to one row per cell line with:
#   - final_neighbors: Semicolon-separated neighbour list (empty string if none)
#   - n_final:         Count of stable neighbours
#   - best_overall_dir: The global best direction used (for reproducibility)
#   - reps_used:       Human-readable description of the intersection rule
#
# Results are sorted by descending n_final to facilitate review of highly
# connected cell lines, with alphabetical sorting as a secondary criterion.

resolved_neighbors_formatted <- resolved_neighbors %>%
  group_by(cell_line) %>%
  summarise(
    final_neighbors = ifelse(
      all(neighbor == "" | is.na(neighbor)),
      "",
      paste(neighbor[neighbor != "" & !is.na(neighbor)], collapse = ";")
    ),
    n_final = sum(neighbor != "" & !is.na(neighbor)),
    reps_used = dplyr::first(reps_used),
    .groups = "drop"
  ) %>%
  mutate(
    best_overall_dir = sub("^best_overall=([^;]+);.*$", "\\1", reps_used),
    winner_dir = sub("^.*;winner=(.*)$", "\\1", reps_used)
  ) %>%
  arrange(desc(n_final), cell_line)

# -----------------------------------------------------------------------------
# Output Path Resolution
# -----------------------------------------------------------------------------
# The output path follows the same resolution logic as input paths:
#
#   1. Explicit --output_tsv argument (highest priority)
#   2. Profile-specific defaults:
#      - multicohort_cancer: {outdir}/unsupervised/tumour_neighbourhoods/final_consensus_all/
#      - Other: {unsup_root}/tumour_neighbourhoods/final_consensus_all/
#
# Parent directories are created automatically if they do not exist.

if (is.null(opt$output_tsv) || !nzchar(opt$output_tsv)) {
  if (profile == "multicohort_cancer") {
    cfg_full <- yaml::read_yaml(opt$config)
    mc_cfg <- cfg_full$multicohort_cancer %||%
      stop("multicohort_cancer section not found in config")
    outdir <- mc_cfg$outdir %||% "results/multicohort_cancer_benchmark"
    outdir <- resolve_path_if_needed(outdir, opt$config)
    output_tsv <- file.path(outdir, "unsupervised", "tumour_neighbourhoods",
                            "final_consensus_all", "resolved_dsmz_neighbors.tsv")
  } else {
    unsup_root <- cfg$paths$unsup_root %||%
      stop("Missing paths.unsup_root in config")
    unsup_root <- resolve_path_if_needed(unsup_root, opt$config)
    output_tsv <- file.path(unsup_root, "tumour_neighbourhoods",
                            "final_consensus_all", "resolved_dsmz_neighbors.tsv")
  }
} else {
  output_tsv <- opt$output_tsv
}

dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Output Writing
# -----------------------------------------------------------------------------
# The formatted results are written as a tab-separated file for downstream
# consumption by analysis scripts and visualisation tools.

write_tsv(resolved_neighbors_formatted, output_tsv)

# -----------------------------------------------------------------------------
# Execution Summary
# -----------------------------------------------------------------------------
# The script outputs key statistics to facilitate validation and logging:
#   - Output file path
#   - Total cell lines processed
#   - Count of isolates (no stable neighbours)
#   - Count of non-isolates (at least one stable neighbour)

cat("Resolved neighbours written to:", output_tsv, "\n")
cat("Total cell lines:", nrow(resolved_neighbors_formatted), "\n")
isolates_n <- sum(resolved_neighbors_formatted$n_final == 0)
cat("Isolates (n_final = 0):", isolates_n, "\n")
cat("Non-isolates (n_final > 0):",
    sum(resolved_neighbors_formatted$n_final > 0), "\n")
