#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(yaml)
  library(optparse)
  library(scales)
  library(rlang)
  library(stringr)
})

# Optional: shared theme if available
# source("R/brca_unsup_methods.R")

# Helper operator for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  make_option("--threshold", type = "double", default = 0.7,
              help = "Strong-consensus threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Source lib_config.R using absolute path (works with --directory)
# Derive pipeline root from config path
if (!is.null(opt$config) && file.exists(opt$config)) {
  config_dir <- dirname(normalizePath(opt$config))
  pipe_root <- normalizePath(file.path(config_dir, ".."))
  lib_config_path <- file.path(pipe_root, "scripts", "lib_config.R")
} else {
  # Fallback: try relative path (for manual runs)
  lib_config_path <- "scripts/lib_config.R"
}
if (file.exists(lib_config_path)) {
  source(lib_config_path)
} else {
  stop("Cannot find lib_config.R. Tried: ", lib_config_path)
}

# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
if (is.null(cfg$paths$unsup_root)) stop("paths$unsup_root is not defined in config.yaml")

# --- Resolve relative paths against pipeline root (derived from config path) ---
# This ensures paths work correctly when Snakemake runs with --directory
pipe_root <- normalizePath(file.path(dirname(opt$config), ".."))

abs_from_root <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)  # already absolute
  file.path(pipe_root, p)
}

unsup_root <- abs_from_root(cfg$paths$unsup_root)
base_dir   <- file.path(unsup_root, "tumour_neighbourhoods")
thr        <- opt$threshold

# Directions present in this project (BRCA may have PAM50; NBL/RBL won't)
default_dirs <- c("hvg_euc", "hvg_corr")

directions <- cfg$tumour_neighbourhoods$directions
if (is.null(directions) || length(directions) == 0) {
  # Auto-discover available directions from filesystem
  if (dir.exists(base_dir)) {
    available_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
    available_dirs <- available_dirs[available_dirs != "final_consensus_all"]
    # Check which ones have final_consensus RDS files
    cons_files <- file.path(base_dir, available_dirs, "final_consensus", 
                           paste0("Final_consensus_tumour_neighbourhoods_", available_dirs, ".rds"))
    existing <- file.exists(cons_files)
    directions <- available_dirs[existing]
    if (length(directions) > 0) {
      cat("[INFO] Auto-discovered", length(directions), "directions from filesystem\n")
    } else {
      directions <- default_dirs
    }
  } else {
    directions <- default_dirs
  }
}
directions <- as.character(directions)

# Filter out PAM50 directions if use_pam50 is false (e.g., for RBL/NBL profiles)
use_pam50 <- cfg$tumour_neighbourhoods$use_pam50 %||% FALSE
if (!use_pam50) {
  directions <- directions[!grepl("^pam50_", directions)]
  if (length(directions) == 0) {
    stop("No valid directions found after filtering PAM50. Check config tumour_neighbourhoods.directions and use_pam50.")
  }
}

# Helper to locate the expected RDS file
cons_rds_path <- function(direction) {
  file.path(
    base_dir, direction, "final_consensus",
    sprintf("Final_consensus_tumour_neighbourhoods_%s.rds", direction)
  )
}

paths <- tibble(
  direction = directions,
  rds = map_chr(directions, cons_rds_path),
  exists = map_lgl(rds, file.exists)
)

if (!all(paths$exists)) {
  missing <- paths %>% filter(!exists)
  stop("Missing final_consensus RDS for:\n", paste(missing$direction, collapse = ", "),
       "\nExpected paths:\n", paste(missing$rds, collapse = "\n"))
}

# Read all consensus tables
dat <- paths %>%
  mutate(df = map(rds, readRDS)) %>%
  select(direction, df) %>%
  unnest(df) %>%
  mutate(direction = factor(direction, levels = directions))

# Ensure cell_line column exists (use cell_tech_id if cell_line is missing)
if (!"cell_line" %in% colnames(dat) && "cell_tech_id" %in% colnames(dat)) {
  dat <- dat %>% mutate(cell_line = cell_tech_id)
}

# Clean p_consensus just in case (prevents those histogram warnings)
dat <- dat %>%
  filter(!is.na(p_consensus), p_consensus >= 0, p_consensus <= 1)

# ---- 1) Direction-level summary ----
# First compute per-cell-line metrics, then aggregate to direction level
# This ensures consistency with winners determination (which uses per-cell-line metrics)
cell_dir_summary_temp <- dat %>%
  group_by(direction, cell_line) %>%
  summarise(
    n_pairs_cell = n(),
    frac_ge_thr_cell = mean(p_consensus >= thr),
    frac_eq_1_cell   = mean(p_consensus == 1),
    median_p_cell    = median(p_consensus),
    p90_cell         = quantile(p_consensus, 0.90),
    .groups = "drop"
  )

# Aggregate direction summary from per-cell-line metrics (weighted by cell line, not pairs)
dir_summary <- cell_dir_summary_temp %>%
  group_by(direction) %>%
  summarise(
    n_pairs = sum(n_pairs_cell),  # Total pairs across all cell lines
    n_cell_lines = n(),  # Number of cell lines
    frac_ge_thr = mean(frac_ge_thr_cell),  # Mean of per-cell-line frac_ge_thr (consistent with winners)
    frac_eq_1   = mean(frac_eq_1_cell),    # Mean of per-cell-line frac_eq_1
    median_p    = median(median_p_cell),   # Median of per-cell-line median_p
    p90         = quantile(p90_cell, 0.90), # 90th percentile of per-cell-line p90
    .groups = "drop"
  ) %>%
  arrange(direction)

# ---- 2) Cell-line x direction summary ----
cell_dir_summary <- dat %>%
  group_by(direction, cell_line) %>%
  summarise(
    n_pairs     = n(),
    max_p       = max(p_consensus),
    median_p    = median(p_consensus),
    frac_ge_thr = mean(p_consensus >= thr),
    frac_eq_1   = mean(p_consensus == 1),
    .groups = "drop"
  ) %>%
  arrange(direction, desc(max_p), desc(frac_ge_thr))

# ---- Standardise column names (backwards compatible) ----
if (!("frac_ge_thr" %in% names(cell_dir_summary)) && ("frac_ge_0_7" %in% names(cell_dir_summary))) {
  cell_dir_summary <- dplyr::rename(cell_dir_summary, frac_ge_thr = frac_ge_0_7)
}

if (!("frac_eq_1" %in% names(cell_dir_summary)) && ("frac_eq_1.0" %in% names(cell_dir_summary))) {
  cell_dir_summary <- dplyr::rename(cell_dir_summary, frac_eq_1 = frac_eq_1.0)
}

cat("[DEBUG] cell_dir_summary columns:\n")
print(names(cell_dir_summary))

# ==============================================================================
# PCA-based composite performance score (direction-agnostic)
# ==============================================================================
# Metrics to use (keep generic — no gene-set assumptions)
metric_cols <- c("frac_ge_thr", "median_p", "max_p", "frac_eq_1")

# Safety check
missing_metrics <- setdiff(metric_cols, colnames(cell_dir_summary))
if (length(missing_metrics) > 0) {
  stop("Missing required metrics in cell_dir_summary: ",
       paste(missing_metrics, collapse = ", "))
}

# ---- Prepare matrix for PCA ----
pca_input <- cell_dir_summary %>%
  select(cell_line, direction, all_of(metric_cols)) %>%
  drop_na()

X <- pca_input %>%
  select(all_of(metric_cols)) %>%
  as.matrix()

# Drop near-constant metrics to avoid NaNs in scaling
sds <- apply(X, 2, sd, na.rm = TRUE)
keep <- which(sds > 1e-12)
if (length(keep) < 2) stop("Not enough variable metrics for PCA composite score.")
X <- X[, keep, drop = FALSE]
metric_cols_used <- colnames(X)

# Standardise metrics (important)
X_scaled <- scale(X)

# ---- PCA ----
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

# Use PC1 loadings as data-driven weights
weights <- abs(pca$rotation[, 1])
weights <- weights / sum(weights)

weights_tbl <- tibble(
  metric = names(weights),
  weight = as.numeric(weights)
)

# ---- PCA diagnostics ----
eig <- (pca$sdev)^2
var_expl <- eig / sum(eig)

p_scree <- tibble(PC = paste0("PC", seq_along(var_expl)), var = var_expl) %>%
  ggplot(aes(x = PC, y = var)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "PCA scree plot", y = "Variance explained", x = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p_load <- weights_tbl %>%
  mutate(metric = factor(metric, levels = metric[order(weight, decreasing = TRUE)])) %>%
  ggplot(aes(x = metric, y = weight)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "PC1-derived metric weights", y = "Weight", x = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1))

cat("\n=== PCA-based composite weights ===\n")
cat("Weights derived from 1st principal component:\n")
print(weights_tbl)
cat(sprintf("PC1 variance explained: %.1f%%\n", 100 * var_expl[1]))

# ---- Composite score per (cell_line, direction) ----
pca_scores <- as.numeric(X_scaled %*% weights)

rank_tbl <- pca_input %>%
  mutate(composite_score = pca_scores)

# Merge composite scores back into cell_dir_summary
cell_dir_summary <- cell_dir_summary %>%
  left_join(
    rank_tbl %>% select(cell_line, direction, composite_score),
    by = c("cell_line", "direction")
  )

# Update rankable table to include composite_score
cell_dir_summary_rankable <- cell_dir_summary %>%
  arrange(cell_line, desc(composite_score), desc(frac_ge_thr), desc(median_p), desc(max_p), direction)

# ---- Best direction per cell line ----
best_per_cell <- rank_tbl %>%
  group_by(cell_line) %>%
  slice_max(order_by = composite_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composite_score)) %>%
  mutate(rank = row_number())

# ---- Final ranked table (clean, thesis-ready) ----
# best_per_cell already contains the metrics (because it comes from rank_tbl = pca_input + composite_score)
ranked_cell_lines <- best_per_cell %>%
  mutate(CPI = composite_score) %>%  # Consensus Performance Index (alias for composite_score)
  transmute(
    rank,
    cell_line,
    best_direction = direction,
    CPI,
    composite_score,  # Keep for backwards compatibility
    frac_ge_thr,
    median_p,
    max_p,
    frac_eq_1
  ) %>%
  arrange(rank)

# ---- Top-N selection helpers ----
top_frac <- cfg$tumour_neighbourhoods$top_fraction %||% 0.10
min_score <- cfg$tumour_neighbourhoods$min_composite_score %||% NA_real_

n_top <- max(1, ceiling(nrow(ranked_cell_lines) * top_frac))

top_by_fraction <- ranked_cell_lines %>% slice_head(n = n_top)

top_by_score <- if (is.na(min_score)) {
  ranked_cell_lines[0, ]  # empty if not set
} else {
  ranked_cell_lines %>% filter(composite_score >= min_score)
}

# Keep best_cell_lines for backward compatibility (same as ranked_cell_lines but without rank)
best_cell_lines <- ranked_cell_lines %>%
  select(-rank) %>%
  arrange(desc(composite_score), desc(frac_ge_thr), desc(median_p), desc(max_p), cell_line)

# Optional: keep a long table (all directions per cell line) but sorted nicely
# (will be updated after composite_score is added)

# ---- 2b) Winner maps per cell line ----
# Tie handling: if ties exist, join multiple directions with ";"

winner_by <- function(df, metric_col) {
  metric_col <- rlang::ensym(metric_col)
  
  df %>%
    group_by(cell_line) %>%
    mutate(metric = !!metric_col) %>%
    summarise(
      best_value = max(metric, na.rm = TRUE),
      best_dirs  = paste(direction[metric == best_value], collapse = ";"),
      .groups = "drop"
    )
}

# Winners
win_max_p <- winner_by(cell_dir_summary, max_p) %>%
  rename(best_max_p = best_value, best_dir_max_p = best_dirs)

win_frac <- winner_by(cell_dir_summary, frac_ge_thr) %>%
  rename(best_frac_ge_thr = best_value, best_dir_frac_ge_thr = best_dirs)

# Wide table for per-direction values (nice for thesis tables)
cell_dir_wide <- cell_dir_summary %>%
  select(direction, cell_line, max_p, frac_ge_thr, median_p, n_pairs) %>%
  pivot_wider(
    names_from = direction,
    values_from = c(max_p, frac_ge_thr, median_p, n_pairs),
    names_glue = "{.value}__{direction}"
  )

# Final "scoreboard"
winner_map <- cell_dir_wide %>%
  left_join(win_max_p, by = "cell_line") %>%
  left_join(win_frac,  by = "cell_line") %>%
  arrange(best_dir_frac_ge_thr, desc(best_frac_ge_thr), best_dir_max_p, desc(best_max_p), cell_line)

# ---- 4) Write outputs ----
out_dir <- file.path(base_dir, "final_consensus_all")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

winner_direction <- dir_summary %>%
  arrange(desc(frac_ge_thr), desc(median_p), desc(p90), direction) %>%
  slice_head(n = 1) %>%
  pull(direction)

if (length(winner_direction) == 0) {
  stop("Unable to determine winning direction from dir_summary.")
}

winner_direction <- as.character(winner_direction)
writeLines(winner_direction, file.path(out_dir, "winning_direction.txt"))

write_tsv(dir_summary, file.path(out_dir, "p_consensus_direction_summary.tsv"))
write_tsv(cell_dir_summary, file.path(out_dir, "p_consensus_cellline_direction_summary.tsv"))
write_tsv(winner_map, file.path(out_dir, "p_consensus_winner_map_per_cell_line.tsv"))
write_tsv(
  win_max_p %>% arrange(desc(best_max_p), cell_line),
  file.path(out_dir, "p_consensus_winners_by_max_p.tsv")
)
write_tsv(
  win_frac %>% arrange(desc(best_frac_ge_thr), cell_line),
  file.path(out_dir, "p_consensus_winners_by_frac_ge_thr.tsv")
)
write_tsv(ranked_cell_lines, file.path(out_dir, "p_consensus_best_cell_lines_ranked.tsv"))
write_tsv(cell_dir_summary_rankable, file.path(out_dir, "p_consensus_cellline_direction_summary.long.tsv"))
write_tsv(weights_tbl, file.path(out_dir, "p_consensus_composite_weights.tsv"))
write_tsv(top_by_fraction, file.path(out_dir, "p_consensus_best_cell_lines_top_fraction.tsv"))
# Always write top_by_score (even if empty) for Snakemake compatibility
write_tsv(top_by_score, file.path(out_dir, "p_consensus_best_cell_lines_top_score.tsv"))

# ---- Write PCA diagnostic plots ----
ggsave(file.path(out_dir, "Fig_p_consensus_composite_PCA_scree.pdf"),
       p_scree, width = 7.5, height = 4.5)

ggsave(file.path(out_dir, "Fig_p_consensus_composite_PCA_weights.pdf"),
       p_load, width = 7.5, height = 4.5)

# ---- 5) Simple comparison plot (thesis-quality) ----
# Canonical x-axis order for direction comparison (same across profiles, e.g. NBL and RBL)
# Feature order then _euc before _corr within each feature; hvg last.
direction_plot_order <- c(
  "Variance_euc", "Variance_corr", "MAD_euc", "MAD_corr",
  "MeanAbsDev_euc", "MeanAbsDev_corr", "Entropy_euc", "Entropy_corr",
  "PCA_euc", "PCA_corr", "Spearman_euc", "Spearman_corr",
  "MX_euc", "MX_corr", "kTotal_euc", "kTotal_corr",
  "hvg_euc", "hvg_corr",
  "pam50_euc", "pam50_corr"  # if present
)
# Use only directions that exist in the data, in canonical order; append any others at end
dirs_in_data <- unique(as.character(dir_summary$direction))
order_used <- c(direction_plot_order[direction_plot_order %in% dirs_in_data],
                 setdiff(dirs_in_data, direction_plot_order))

feature_cols <- c(
  Variance    = "#D4AF37",
  MAD         = "#F28E2B",
  MeanAbsDev  = "#00BFC4",
  Entropy     = "#4682B4",
  PCA         = "#9E9E9E",
  Spearman    = "#5E2B97",
  MX          = "#2E8B57",
  kTotal      = "#B00020",
  hvg         = "#800000",
  pam50       = "#C71585"
)

# Prepare data: extract feature names and add percentage labels; fix direction order for plot
dir_summary_plot <- dir_summary %>%
  mutate(
    # feature = everything before the last _euc/_corr
    feature  = str_replace(direction, "_(euc|corr)$", ""),
    # label inside bar: percent
    pct_lab  = percent(frac_ge_thr, accuracy = 1),
    # enforce canonical x-axis order (same as RBL)
    direction = factor(as.character(direction), levels = order_used)
  )
# Force same legend/order so NBL matches RBL (stable fill mapping); pam50 last for BRCA
feature_levels <- c("Variance", "MAD", "MeanAbsDev", "Entropy", "PCA", "Spearman", "MX", "kTotal", "hvg", "pam50")
dir_summary_plot <- dir_summary_plot %>%
  mutate(feature = factor(as.character(feature), levels = feature_levels))

# Safety check: no feature should silently get grey
missing <- setdiff(levels(dir_summary_plot$feature), names(feature_cols))
if (length(missing) > 0) stop("No colour mapping for: ", paste(missing, collapse = ", "))

# Use fill values in exact same order as factor levels so MeanAbsDev and hvg never swap
fill_values <- feature_cols[levels(dir_summary_plot$feature)]

# Create plot with vertical bars, fixed feature colors, and percentage labels inside bars
p <- ggplot(dir_summary_plot, aes(x = direction, y = frac_ge_thr, fill = feature)) +
  geom_col(width = 0.75) +
  # percentage INSIDE bar
  geom_text(
    aes(label = pct_lab),
    vjust = 1.2,              # pushes text slightly down into bar
    size = 3.2,
    color = "white"
  ) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Direction",
    y = "Fraction of (cell line, tumour) pairs\nwith strong consensus",
    fill = "Feature"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 18)
  )

ggsave(file.path(out_dir, "Fig_p_consensus_direction_comparison.pdf"), p, width = 10, height = 6)
# Save PNG version (high-quality, slide-ready)
ggsave(
  file.path(out_dir, "Fig_p_consensus_direction_comparison.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- 6) Winner map summary statistics ----
# Count winners per direction (handling ties by counting all tied directions)
winner_summary_max_p <- win_max_p %>%
  separate_rows(best_dir_max_p, sep = ";") %>%
  count(best_dir_max_p, name = "n_cell_lines_max_p") %>%
  arrange(desc(n_cell_lines_max_p))

winner_summary_frac <- win_frac %>%
  separate_rows(best_dir_frac_ge_thr, sep = ";") %>%
  count(best_dir_frac_ge_thr, name = "n_cell_lines_frac") %>%
  arrange(desc(n_cell_lines_frac))

cat("\n=== Winner map summary ===\n")
cat("Cell lines with best max_p_consensus per direction:\n")
print(winner_summary_max_p)

cat("\nCell lines with best frac_ge_thr per direction:\n")
print(winner_summary_frac)

cat("\n=== Top 10 cell lines by best max_p_consensus ===\n")
print(winner_map %>% head(10) %>% select(cell_line, best_dir_max_p, best_max_p, best_dir_frac_ge_thr, best_frac_ge_thr))

cat("\n=== Top 20 best-performing cell lines (overall, by composite score) ===\n")
print(ranked_cell_lines %>% head(20))

cat(sprintf("\n=== Top %.0f%% cell lines (n=%d) ===\n", 100 * top_frac, n_top))
print(top_by_fraction)

cat("\n[SUCCESS] Wrote:\n")
cat("  ", file.path(out_dir, "winning_direction.txt"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_direction_summary.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_cellline_direction_summary.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_winner_map_per_cell_line.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_winners_by_max_p.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_winners_by_frac_ge_thr.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_best_cell_lines_ranked.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_cellline_direction_summary.long.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_composite_weights.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_best_cell_lines_top_fraction.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "p_consensus_best_cell_lines_top_score.tsv"), "\n", sep = "")
cat("  ", file.path(out_dir, "Fig_p_consensus_direction_comparison.pdf"), "\n", sep = "")
cat("  ", file.path(out_dir, "Fig_p_consensus_composite_PCA_scree.pdf"), "\n", sep = "")
cat("  ", file.path(out_dir, "Fig_p_consensus_composite_PCA_weights.pdf"), "\n", sep = "")
