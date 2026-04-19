#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(optparse)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml"),
  make_option("--profile", type = "character", default = NULL),
  make_option("--out-tsv", type = "character"),
  make_option("--out-plot", type = "character"),
  make_option("--out-notes", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

config_dir <- dirname(normalizePath(opt$config))
pipe_root <- normalizePath(file.path(config_dir, ".."))
source(file.path(pipe_root, "scripts", "lib_config.R"))
raw_cfg <- yaml::read_yaml(opt$config)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
abs_from_root <- function(p) if (grepl("^/", p)) p else file.path(pipe_root, p)
unsup_root <- abs_from_root(cfg$paths$unsup_root)
final_dir <- file.path(unsup_root, "tumour_neighbourhoods", "final_consensus_all")
shortlist_n <- raw_cfg$validation$shortlist_n %||% 10
n_perm <- raw_cfg$validation$permutation_iterations %||% 1000
seed <- raw_cfg$validation$seed %||% 42
set.seed(seed)

ranked <- read_tsv(file.path(final_dir, "p_consensus_best_cell_lines_ranked.tsv"), show_col_types = FALSE)
long_tbl <- read_tsv(file.path(final_dir, "p_consensus_cellline_direction_summary.long.tsv"), show_col_types = FALSE)
selected_n <- min(shortlist_n, nrow(ranked))
obs <- ranked %>% slice_head(n = selected_n)
obs_mean <- mean(obs$composite_score)
obs_median <- median(obs$composite_score)

perm_stats <- bind_rows(lapply(seq_len(n_perm), function(i) {
  permuted <- long_tbl %>%
    group_by(direction) %>%
    mutate(cell_line_perm = sample(cell_line)) %>%
    ungroup()
  perm_best <- permuted %>%
    transmute(
      cell_line = cell_line_perm,
      direction,
      composite_score,
      frac_ge_thr,
      median_p,
      max_p
    ) %>%
    group_by(cell_line) %>%
    slice_max(order_by = composite_score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(composite_score), desc(frac_ge_thr), desc(median_p), desc(max_p), cell_line) %>%
    slice_head(n = selected_n)
  tibble(iteration = i, mean_selection_score = mean(perm_best$composite_score), median_selection_score = median(perm_best$composite_score))
}))

summary_tbl <- tibble(
  cohort = profile,
  shortlist_n = selected_n,
  permutation_iterations = n_perm,
  observed_mean_selection_score = obs_mean,
  null_mean_selection_score = mean(perm_stats$mean_selection_score),
  null_sd_selection_score = sd(perm_stats$mean_selection_score),
  empirical_p_mean = (sum(perm_stats$mean_selection_score >= obs_mean) + 1) / (n_perm + 1),
  observed_median_selection_score = obs_median,
  empirical_p_median = (sum(perm_stats$median_selection_score >= obs_median) + 1) / (n_perm + 1)
)
write_tsv(summary_tbl, opt$`out-tsv`)

p <- ggplot(perm_stats, aes(x = mean_selection_score)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs_mean, colour = "firebrick", linewidth = 1) +
  labs(
    title = sprintf("Permutation null for shortlist mean score: %s", profile),
    x = "Permuted shortlist mean selection score",
    y = "Count"
  ) +
  theme_bw(base_size = 11)
ggsave(opt$`out-plot`, p, width = 7.5, height = 4.8)

notes <- c(
  sprintf("Profile: %s", profile),
  sprintf("Null size: %d permutations.", n_perm),
  "This is a post-hoc aggregated null built from p_consensus_cellline_direction_summary.long.tsv, not a rerun of the raw-expression clustering workflow.",
  "Within each direction, cell-line labels are permuted while per-direction score distributions are preserved.",
  "Use this as a stability sanity check for shortlist concentration, not as a replacement for the primary ranking rule."
)
writeLines(notes, opt$`out-notes`)
