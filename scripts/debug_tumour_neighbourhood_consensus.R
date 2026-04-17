#!/usr/bin/env Rscript

# debug_tumour_neighbourhood_consensus.R
#
# Sanity-check a single direction:
#   - pam50_euc
#   - pam50_corr
#   - HVG_euc
#   - HVG_corr
#
# Uses:
#   config.yaml → paths$unsup_root
#   Final_consensus_tumour_neighbourhoods_<direction>.rds
#
# Prints:
#   - basic dimensions and a glimpse of consensus_pairs
#   - per-cell-line anchoring summary (max p_consensus, #tumours)
#   - top/bottom cell lines by max_p_consensus
#   - optional peek at one Top_m_long_* CSV if it exists

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/lib_config.R")

option_list <- list(
  make_option(
    c("-c", "--config"),
    type    = "character",
    default = "config/config.yaml",
    help    = "Path to config.yaml [default %default]"
  ),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  make_option(
    "--direction",
    type    = "character",
    default = "pam50_euc",
    help    = "Direction: pam50_euc, pam50_corr, HVG_euc, HVG_corr"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== debug_tumour_neighbourhood_consensus.R ===\n")
cat("Config:    ", opt$config,    "\n")
cat("Direction: ", opt$direction, "\n\n")

# ----------------------------------------------------------------------
# 1) Load config, locate consensus file
# ----------------------------------------------------------------------
# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
unsup_root <- cfg$paths$unsup_root

nh_root   <- file.path(unsup_root, "tumour_neighbourhoods", opt$direction)
cons_dir  <- file.path(nh_root, "final_consensus")
cons_rds  <- file.path(
  cons_dir,
  sprintf("Final_consensus_tumour_neighbourhoods_%s.rds", opt$direction)
)

cat("Tumour neighbourhood root:\n  ", nh_root,  "\n")
cat("Final consensus RDS:\n  ", cons_rds, "\n\n")

if (!file.exists(cons_rds)) {
  stop("[ERROR] Consensus RDS not found. Expected at:\n  ", cons_rds,
       "\nMake sure tumour_neighbourhood_p_consensus has run for this direction.")
}

consensus_pairs <- readRDS(cons_rds)

if (!"tumour_id" %in% colnames(consensus_pairs)) {
  if ("tumor_id" %in% colnames(consensus_pairs)) {
    consensus_pairs <- consensus_pairs %>% dplyr::rename(tumour_id = tumor_id)
  } else {
    stop("[ERROR] consensus_pairs missing tumour_id column")
  }
}

cat("=== 1) Basic structure of consensus_pairs ===\n\n")
cat("Class: ", class(consensus_pairs), "\n")
cat("Rows:  ", nrow(consensus_pairs), "\n\n")

cat("Columns:\n")
print(colnames(consensus_pairs))

cat("\nHead:\n")
print(utils::head(consensus_pairs, 10))

# ----------------------------------------------------------------------
# 2) Basic stats: unique counts
# ----------------------------------------------------------------------
cat("\n=== 2) Unique DSMZ cell lines and TCGA tumours ===\n\n")

if (!all(c("cell_line", "tumour_id", "p_consensus") %in% colnames(consensus_pairs))) {
  stop("[ERROR] consensus_pairs must contain columns: cell_line, tumour_id, p_consensus.")
}

n_cell <- length(unique(consensus_pairs$cell_line))
n_tum  <- length(unique(consensus_pairs$tumour_id))

cat("Unique DSMZ cell lines: ", n_cell, "\n")
cat("Unique TCGA tumours:    ", n_tum,  "\n\n")

# ----------------------------------------------------------------------
# 3) Global p_consensus distribution
# ----------------------------------------------------------------------
cat("=== 3) Global p_consensus distribution ===\n\n")

pc <- consensus_pairs$p_consensus

cat("Summary of p_consensus:\n")
print(summary(pc))

cat("\nQuantiles of p_consensus:\n")
print(quantile(pc, probs = seq(0, 1, by = 0.1), na.rm = TRUE))

# ----------------------------------------------------------------------
# 4) Per-cell-line anchoring summary
# ----------------------------------------------------------------------
cat("\n=== 4) Per-cell-line anchoring statistics ===\n\n")

per_cell <- consensus_pairs %>%
  group_by(cell_line) %>%
  summarise(
    n_tumours       = n_distinct(tumour_id),
    max_p_consensus = max(p_consensus, na.rm = TRUE),
    mean_p_consensus = mean(p_consensus, na.rm = TRUE),
    frac_ge_0_7     = mean(p_consensus >= 0.7, na.rm = TRUE),
    frac_ge_0_5     = mean(p_consensus >= 0.5, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  arrange(desc(max_p_consensus))

cat("Top 10 cell lines by max_p_consensus:\n")
print(per_cell %>% slice(1:10))

cat("\nBottom 10 cell lines by max_p_consensus:\n")
print(per_cell %>% slice(max(1, n() - 9):n()))

# Save per-cell summary (optional but handy)
per_cell_out <- file.path(
  cons_dir,
  sprintf("debug_per_cell_anchoring_%s.tsv", opt$direction)
)
readr::write_tsv(per_cell, per_cell_out)
cat("\nPer-cell anchoring summary written to:\n  ", per_cell_out, "\n")

# ----------------------------------------------------------------------
# 5) Optional: peek at one Top_m_long_* CSV if present
# ----------------------------------------------------------------------
cat("\n=== 5) Optional: check one Top_m_long_* file ===\n\n")

example_csv <- file.path(
  nh_root,
  "HC_PAM50_dynamic",
  "Top_m_long_HC_PAM50_dynamic.csv"
)

if (file.exists(example_csv)) {
  cat("Found example Top_m_long CSV:\n  ", example_csv, "\n\n")
  top_long <- readr::read_csv(example_csv, show_col_types = FALSE)

  if (!"tumour_id" %in% colnames(top_long) && "tumor_id" %in% colnames(top_long)) {
    top_long <- top_long %>% dplyr::rename(tumour_id = tumor_id)
  }

  cat("Top_m_long_HC_PAM50_dynamic: rows =", nrow(top_long), "\n")
  cat("Columns:\n")
  print(colnames(top_long))

  if (all(c("cell_line", "tumour_id", "rank") %in% colnames(top_long))) {
    per_cell_top <- top_long %>%
      group_by(cell_line) %>%
      summarise(
        n_top_tumours = n_distinct(tumour_id),
        max_rank      = max(rank, na.rm = TRUE),
        .groups       = "drop"
      ) %>%
      arrange(desc(n_top_tumours))

    cat("\nPer-cell stats in this method (HC_PAM50_dynamic):\n")
    print(utils::head(per_cell_top, 10))
  } else {
    cat("Note: expected columns 'cell_line', 'tumour_id', 'rank' not all found.\n")
  }
} else {
  cat("No example Top_m_long file found at:\n  ", example_csv, "\n")
  cat("This is OK if tumour_neighbourhoods hasn't run yet for this direction.\n")
}

cat("\n=== Done: debug_tumour_neighbourhood_consensus.R ===\n")

