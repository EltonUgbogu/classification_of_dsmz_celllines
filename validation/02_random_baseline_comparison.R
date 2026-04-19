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
n_iter <- raw_cfg$validation$random_baseline_iterations %||% 1000
seed <- raw_cfg$validation$seed %||% 42
set.seed(seed)

ranked <- read_tsv(file.path(final_dir, "p_consensus_best_cell_lines_ranked.tsv"), show_col_types = FALSE)
selected_n <- min(shortlist_n, nrow(ranked))
obs <- ranked %>% slice_head(n = selected_n)
obs_summary <- obs %>% summarise(
  mean_selection_score = mean(composite_score),
  mean_frac_ge_thr = mean(frac_ge_thr),
  mean_median_p = mean(median_p),
  mean_max_p = mean(max_p)
)

rand_stats <- bind_rows(lapply(seq_len(n_iter), function(i) {
  samp <- ranked %>% slice_sample(n = selected_n)
  tibble(
    iteration = i,
    mean_selection_score = mean(samp$composite_score),
    mean_frac_ge_thr = mean(samp$frac_ge_thr),
    mean_median_p = mean(samp$median_p),
    mean_max_p = mean(samp$max_p)
  )
}))

summary_tbl <- tibble(
  cohort = profile,
  shortlist_n = selected_n,
  random_iterations = n_iter,
  observed_mean_selection_score = obs_summary$mean_selection_score,
  random_mean_selection_score = mean(rand_stats$mean_selection_score),
  z_mean_selection_score = (obs_summary$mean_selection_score - mean(rand_stats$mean_selection_score)) / sd(rand_stats$mean_selection_score),
  empirical_p_mean_selection_score = (sum(rand_stats$mean_selection_score >= obs_summary$mean_selection_score) + 1) / (n_iter + 1),
  observed_mean_frac_ge_thr = obs_summary$mean_frac_ge_thr,
  observed_mean_median_p = obs_summary$mean_median_p,
  observed_mean_max_p = obs_summary$mean_max_p
)
write_tsv(summary_tbl, opt$`out-tsv`)

p <- ggplot(rand_stats, aes(x = mean_selection_score)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs_summary$mean_selection_score, colour = "steelblue4", linewidth = 1) +
  labs(
    title = sprintf("Random shortlist baseline: %s", profile),
    x = "Random-set mean selection score",
    y = "Count"
  ) +
  theme_bw(base_size = 11)
ggsave(opt$`out-plot`, p, width = 7.5, height = 4.8)

notes <- c(
  sprintf("Profile: %s", profile),
  sprintf("Random baseline iterations: %d.", n_iter),
  sprintf("Observed shortlist size: top %d cell lines by the existing composite score.", selected_n),
  "The comparator is a random cell-line set of the same size sampled from the ranked output table.",
  "This baseline is descriptive only; it does not replace the consensus ranking procedure."
)
writeLines(notes, opt$`out-notes`)
