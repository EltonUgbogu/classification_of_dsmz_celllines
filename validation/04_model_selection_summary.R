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

ranked <- read_tsv(file.path(final_dir, "p_consensus_best_cell_lines_ranked.tsv"), show_col_types = FALSE)
long_tbl <- read_tsv(file.path(final_dir, "p_consensus_cellline_direction_summary.long.tsv"), show_col_types = FALSE)
winning_direction <- readLines(file.path(final_dir, "winning_direction.txt"), warn = FALSE)[1] %||% NA_character_
selected <- ranked %>% slice_head(n = min(shortlist_n, n()))

metric_tbl <- long_tbl %>%
  transmute(
    cell_line,
    best_direction = direction,
    n_supported_tumours = n_pairs,
    frac_ge_thr_long = frac_ge_thr,
    median_p_long = median_p,
    max_p_long = max_p,
    composite_score_long = composite_score
  )

load_graph_annotations <- function(direction) {
  p <- file.path(unsup_root, "tumour_neighbourhoods", direction, "final_consensus",
                 sprintf("cell_line_similarity_graph_node_annotations_%s.tsv", direction))
  if (!file.exists(p)) {
    return(tibble(best_direction = character(), cell_line = character(), community_id = character(), graph_component = character()))
  }
  read_tsv(p, show_col_types = FALSE) %>%
    transmute(
      best_direction = direction,
      cell_line,
      community_id = coalesce(as.character(community_leid), as.character(community_louv)),
      graph_component = as.character(component)
    )
}

annot_tbl <- bind_rows(lapply(unique(selected$best_direction), load_graph_annotations))

summary_tbl <- selected %>%
  left_join(metric_tbl, by = c("cell_line", "best_direction")) %>%
  left_join(annot_tbl, by = c("cell_line", "best_direction")) %>%
  transmute(
    cohort = profile,
    cell_line_id = cell_line,
    rank,
    selection_score = coalesce(composite_score, composite_score_long, CPI),
    median_p_consensus = coalesce(median_p, median_p_long),
    max_p_consensus = coalesce(max_p, max_p_long),
    frac_ge_thr = coalesce(frac_ge_thr, frac_ge_thr_long),
    n_supported_tumours,
    best_direction,
    community_id,
    graph_component,
    global_winning_direction = winning_direction
  ) %>%
  arrange(rank)

write_tsv(summary_tbl, opt$`out-tsv`)

plot_tbl <- summary_tbl %>% mutate(cell_line_id = factor(cell_line_id, levels = rev(cell_line_id)))
p <- ggplot(plot_tbl, aes(x = cell_line_id, y = selection_score, fill = best_direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(x = NULL, y = "Selection score", title = sprintf("Top %d candidate models: %s", nrow(plot_tbl), profile)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(opt$`out-plot`, p, width = 8, height = max(4.5, 0.35 * nrow(plot_tbl) + 2))

notes <- c(
  sprintf("Profile: %s", profile),
  sprintf("Shortlist size: top %d ranked cell lines from p_consensus_best_cell_lines_ranked.tsv.", nrow(summary_tbl)),
  "Selection score is the existing composite_score/CPI emitted by summarize_p_consensus_all.R.",
  "n_supported_tumours is taken from the selected cell line's best_direction row in p_consensus_cellline_direction_summary.long.tsv.",
  "community_id is contextual graph metadata from the same best_direction when that annotation exists; it does not redefine ranking."
)
writeLines(notes, opt$`out-notes`)
