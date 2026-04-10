#!/usr/bin/env Rscript

# plot_winners_by_frac_ge_thr.R
# Creates a dot plot matrix showing the best feature representation per cell line
# Supports both per-disease and multicohort_cancer profiles

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(yaml)
  library(optparse)
})

# Load shared config utilities
script_dir <- NULL
script_args <- commandArgs()
file_arg <- script_args[grepl("^--file=", script_args)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  if (file.exists(script_path)) {
    script_dir <- dirname(normalizePath(script_path))
  }
}
if (is.null(script_dir) || length(script_dir) == 0) {
  script_dir <- getwd()
}
lib_config_path <- file.path(script_dir, "lib_config.R")
if (file.exists(lib_config_path)) {
  source(lib_config_path)
} else {
  stop("Cannot find lib_config.R at: ", lib_config_path)
}

# Parse command line arguments
option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml"),
  make_option("--profile", type = "character", default = NULL,
              help = "Profile name (brca, nbl, rbl, or multicohort_cancer)"),
  make_option("--winners_tsv", type = "character", default = NULL,
              help = "Path to p_consensus_winners_by_frac_ge_thr.tsv"),
  make_option("--output_prefix", type = "character", default = NULL,
              help = "Output file prefix (without extension)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Read config
profile_override <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", unset = NA)
if (is.na(profile_override)) profile_override <- NULL

cfg <- read_profiled_config(opt$config, profile_override = profile_override)
profile <- cfg$profile

# Determine output directory
if (!is.null(opt$output_prefix)) {
  output_prefix <- opt$output_prefix
} else {
  if (profile == "multicohort_cancer") {
    cfg_full <- yaml::read_yaml(opt$config)
    pan_cfg <- cfg_full$multicohort_cancer %||% stop("multicohort_cancer section not found in config")
    outdir <- pan_cfg$outdir %||% "results/multicohort_cancer_benchmark"
    if (!dir.exists(outdir)) {
      config_dir <- dirname(normalizePath(opt$config))
      outdir <- file.path(config_dir, "..", outdir)
      outdir <- normalizePath(outdir, mustWork = FALSE)
    }
    final_consensus_all_dir <- file.path(outdir, "unsupervised", "tumour_neighbourhoods", "final_consensus_all")
  } else {
    unsup_root <- cfg$paths$unsup_root
    if (is.null(unsup_root)) {
      stop("Missing paths.unsup_root in config")
    }
    if (!dir.exists(unsup_root)) {
      config_dir <- dirname(normalizePath(opt$config))
      unsup_root <- file.path(config_dir, "..", unsup_root)
      unsup_root <- normalizePath(unsup_root, mustWork = FALSE)
    }
    final_consensus_all_dir <- file.path(unsup_root, "tumour_neighbourhoods", "final_consensus_all")
  }
  output_prefix <- file.path(final_consensus_all_dir, "Fig_p_consensus_winners_by_frac_ge_thr")
}

# Read winners file
if (is.null(opt$winners_tsv)) {
  winners_tsv <- file.path(dirname(output_prefix), "p_consensus_winners_by_frac_ge_thr.tsv")
} else {
  winners_tsv <- opt$winners_tsv
}

if (!file.exists(winners_tsv)) {
  stop("Winners TSV not found: ", winners_tsv)
}

winners <- read_tsv(winners_tsv, show_col_types = FALSE)

# Parse direction to extract feature and distance
winners <- winners %>%
  mutate(
    feature = case_when(
      grepl("^pam50_", best_dir_frac_ge_thr) ~ "PAM50",
      grepl("^HVG_", best_dir_frac_ge_thr) ~ "HVG",
      grepl("_euc$|_euclidean$", best_dir_frac_ge_thr) ~ sub("_(euc|euclidean)$", "", best_dir_frac_ge_thr),
      grepl("_corr$|_cosine$", best_dir_frac_ge_thr) ~ sub("_(corr|cosine)$", "", best_dir_frac_ge_thr),
      TRUE ~ best_dir_frac_ge_thr
    ),
    distance = case_when(
      grepl("_euc$|_euclidean$", best_dir_frac_ge_thr) ~ "Euclidean",
      grepl("_corr$|_cosine$", best_dir_frac_ge_thr) ~ "Correlation",
      TRUE ~ "Unknown"
    )
  )

# Create feature family mapping (one color per feature family)
# Use a named character vector so indexing returns characters, not lists.
feature_families <- c(
  "Variance" = "Variance",
  "MAD" = "MAD",
  "MeanAbsDev" = "MAD",
  "Entropy" = "Entropy",
  "Spearman" = "Spearman",
  "MX" = "MX",
  "kTotal" = "kTotal",
  "PAM50" = "PAM50",
  "HVG" = "HVG"
)

winners <- winners %>%
  mutate(
    feature_family = if_else(
      feature %in% names(feature_families),
      unname(feature_families[feature]),
      feature
    )
  )

# Create dot plot matrix
p <- ggplot(winners, aes(x = reorder(cell_line, best_frac_ge_thr), 
                         y = best_frac_ge_thr,
                         color = feature_family,
                         shape = distance)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", "25", "50", "75", "100"),
                     limits = c(0, 1)) +
  labs(
    x = "Cell Line",
    y = "Fraction of (cell line, tumour) pairs with strong consensus (%)",
    title = paste("Best Feature Representation per Cell Line -", toupper(profile)),
    color = "Feature Family",
    shape = "Distance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"  # Remove legend as requested
  )

# Save plots
ggsave(paste0(output_prefix, ".pdf"), p, width = 12, height = 8)
ggsave(paste0(output_prefix, ".png"), p, width = 12, height = 8, dpi = 300)

cat("Plot saved to:", output_prefix, "\n")
