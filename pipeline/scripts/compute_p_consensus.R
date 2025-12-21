# ==============================================================================
# compute_p_consensus.R
# Consensus Aggregation Across Clustering Methods
# ==============================================================================
#
# This script aggregates tumour neighbourhoods across multiple clustering
# methods to compute p_consensus(c, t): the proportion of methods in which
# tumour t appears in the neighbourhood of cell line c.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
  library(optparse)
})

source("R/brca_unsup_methods.R")

# ------------------------------------------------------------------------------
# 0) Command-line argument parsing
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$direction)) {
  stop("Please supply --direction (pam50_euc|pam50_corr|hvg_euc|hvg_corr).")
}

direction <- opt$direction
if (!direction %in% c("pam50_euc", "pam50_corr", "hvg_euc", "hvg_corr")) {
  stop("Invalid --direction: ", direction,
       " (allowed: pam50_euc, pam50_corr, hvg_euc, hvg_corr)")
}

cat("[INFO] Using config file: ", opt$config, "\n", sep = "")
cat("[INFO] Direction: ", direction, "\n", sep = "")

gene_set <- if (startsWith(direction, "pam50")) "PAM50" else "HVG500"
cat("[INFO] Gene set: ", gene_set, "\n", sep = "")

cfg <- yaml::read_yaml(opt$config)

if (is.null(cfg$paths$unsup_root)) {
  stop("paths$unsup_root is not defined in config.yaml")
}

unsup_root <- cfg$paths$unsup_root
base_dir   <- file.path(unsup_root, "tumour_neighbourhoods", direction)

cat("[INFO] unsup_root: ", unsup_root, "\n", sep = "")
cat("[INFO] tumour_neighbourhoods base_dir: ", base_dir, "\n", sep = "")

# ------------------------------------------------------------------------------
# 1) Define the 10 clustering methods contributing to consensus
# ------------------------------------------------------------------------------
# Each method produces a neighbourhood file containing (cell_line, tumor_id,
# in_top) tuples indicating which tumours are in the top neighbourhood.

methods_tbl <- tibble::tribble(
  ~method_id,               ~subdir,                    ~file,
  "HC_PAM50_dynamic",       "HC_PAM50_dynamic",         "Top_m_long_HC_PAM50_dynamic.csv",
  "HC_PAM50_optimal",       "HC_PAM50_optimal",         "Top_m_long_HC_PAM50_optimal.csv",
  "HC_PCA_dynamic",         "HC_PCA_dynamic",           "Top_m_long_HC_PCA_dynamic.csv",
  "HC_PCA_optimal",         "HC_PCA_optimal",           "Top_m_long_HC_PCA_optimal.csv",
  "km_PAM50",               "km_PAM50",                 "Top_m_long_km_PAM50.csv",
  "km_PCA",                 "km_PCA",                   "Top_m_long_km_PCA.csv",
  "consensus_HC_PAM50",     "consensus_HC_PAM50",       "Top_m_long_consensus_HC_PAM50.csv",
  "consensus_kmeans_PAM50", "consensus_kmeans_PAM50",   "Top_m_long_consensus_kmeans_PAM50.csv",
  "consensus_HC_PCA",       "consensus_HC_PCA",         "Top_m_long_consensus_HC_PCA.csv",
  "consensus_kmeans_PCA",   "consensus_kmeans_PCA",     "Top_m_long_consensus_kmeans_PCA.csv"
)

n_methods_total <- nrow(methods_tbl)
cat("[INFO] Number of clustering methods contributing to p_consensus: ",
    n_methods_total, "\n", sep = "")

# ------------------------------------------------------------------------------
# 2) Load neighbourhood results from all methods
# ------------------------------------------------------------------------------
# The map2() function iterates over paths and method IDs simultaneously,
# reading each CSV and standardising column order.

neigh_list <- methods_tbl %>%
  mutate(
    path = file.path(base_dir, subdir, file),
    data = map2(path, method_id, ~ {
      if (!file.exists(.x)) {
        stop("Neighbourhood file not found: ", .x)
      }
      read_csv(.x, show_col_types = FALSE) %>%
        mutate(method = .y) %>%
        select(method, cell_line, tumor_id, in_top, rank, distance, everything()) %>%
        select(-any_of(c("cluster", "cluster_id", "cluster_label", "subtype")))
    })
  )

# Unnest into a single long-format table
all_long <- neigh_list %>%
  select(method_id, data) %>%
  unnest(data)

cat("=== Sanity check: counts per method / cell line / in_top ===\n")
all_long %>%
  count(method, cell_line, in_top) %>%
  print(n = 12)

# ------------------------------------------------------------------------------
# 3) Compute p_consensus(c, t)
# ------------------------------------------------------------------------------
# p_consensus is the fraction of methods where tumour t appears in the
# neighbourhood of cell line c. A value of 1.0 indicates perfect agreement
# across all clustering methods.

consensus_pairs <- all_long %>%
  filter(in_top) %>%
  count(cell_line, tumor_id, name = "n_methods") %>%
  mutate(p_consensus = n_methods / n_methods_total) %>%
  arrange(cell_line, desc(p_consensus))

# Write outputs
cons_dir <- file.path(base_dir, "final_consensus")
dir.create(cons_dir, showWarnings = FALSE, recursive = TRUE)

cons_basename <- sprintf("Final_consensus_tumour_neighbourhoods_%s", direction)
cons_rds_path <- file.path(cons_dir, paste0(cons_basename, ".rds"))
cons_csv_path <- file.path(cons_dir, paste0(cons_basename, ".csv"))

write_tsv(consensus_pairs, cons_csv_path)
saveRDS(consensus_pairs, cons_rds_path)

cat("\n=== Consensus summary ===\n")
consensus_pairs %>% glimpse()

cat("\n=== Cell lines with perfect consensus (p_consensus == 1) ===\n")
consensus_pairs %>%
  filter(p_consensus == 1) %>%
  count(cell_line, sort = TRUE) %>%
  print()

# ------------------------------------------------------------------------------
# 4) Histogram of p_consensus distribution
# ------------------------------------------------------------------------------
# Visualises the overall distribution of consensus strength across all
# (cell line, tumour) pairs.

n_pairs      <- nrow(consensus_pairs)
frac_ge_0_7  <- mean(consensus_pairs$p_consensus >= 0.7)
frac_eq_1    <- mean(consensus_pairs$p_consensus == 1)

cat("\n=== Summary of consensus strength ===\n")
cat(sprintf("Total pairs:        %d\n", n_pairs))
cat(sprintf("p_consensus >= 0.7 : %.1f%% (%d pairs)\n",
            100 * frac_ge_0_7, sum(consensus_pairs$p_consensus >= 0.7)))
cat(sprintf("p_consensus = 1.0 : %.1f%% (%d pairs)\n",
            100 * frac_eq_1, sum(consensus_pairs$p_consensus == 1)))

annot_text <- sprintf(
  "Total pairs: %d\n>= 0.7: %.1f%%\n= 1.0: %.1f%%",
  n_pairs, 100 * frac_ge_0_7, 100 * frac_eq_1
)

p_hist <- ggplot(consensus_pairs, aes(x = p_consensus)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 0.05,
    colour   = "white"
  ) +
  geom_vline(xintercept = 0.7, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_dimred(base_size = 14) +
  labs(
    title = sprintf("Distribution of consensus strength across %d clustering methods (%s)",
                    n_methods_total, gene_set),
    subtitle = sprintf("%d (cell line, TCGA tumour) pairs | %s gene set", n_pairs, gene_set),
    x = expression(p[consensus]),
    y = "Proportion of pairs"
  ) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()) +
  annotate("text", x = 0.98, y = 0.95, label = annot_text,
           hjust = 1, vjust = 1, size = 3.5)

plot_path_hist <- file.path(cons_dir, sprintf("Fig_p_consensus_distribution_%s.pdf", direction))
ggsave(plot_path_hist, p_hist, width = 8, height = 5)

cat("\nHistogram saved to:\n  ", plot_path_hist, "\n")

# ------------------------------------------------------------------------------
# 5) Per-cell-line summary statistics
# ------------------------------------------------------------------------------
# Summarises consensus strength for each cell line: maximum p_consensus,
# median p_consensus, and fraction of neighbours with strong consensus.

cell_summary <- consensus_pairs %>%
  group_by(cell_line) %>%
  summarise(
    n_pairs     = n(),
    max_p       = max(p_consensus),
    median_p    = median(p_consensus),
    frac_ge_0_7 = mean(p_consensus >= 0.7),
    .groups     = "drop"
  ) %>%
  arrange(desc(max_p))

cat("\n=== Per-cell-line summary (top 10 by max p_consensus) ===\n")
cell_summary %>% head(10) %>% print()

# Bar plot of max p_consensus per cell line
p_cell <- ggplot(
  cell_summary,
  aes(x = reorder(cell_line, max_p), y = max_p, fill = max_p >= 0.7)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#bdbdbd"),
    labels = c("FALSE" = "< 0.7", "TRUE" = ">= 0.7"),
    name   = "Max p_consensus"
  ) +
  theme_dimred(base_size = 12) +
  labs(
    title = sprintf("Per-cell-line consensus with TCGA tumours (%s)", direction),
    subtitle = "Bars show max p_consensus across all tumour neighbours",
    x = "DSMZ breast cancer cell lines",
    y = expression(max~p[consensus])
  ) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

plot_path_cell <- file.path(cons_dir, sprintf("Fig_p_consensus_per_cell_line_%s.pdf", direction))
ggsave(plot_path_cell, p_cell, width = 8, height = 10)

cat("\nPer-cell-line summary plot saved to:\n  ", plot_path_cell, "\n")

# ------------------------------------------------------------------------------
# 6) Scatter plot: max p_consensus vs fraction of strong neighbours
# ------------------------------------------------------------------------------
# This visualisation identifies cell lines with both high maximum consensus
# and a large proportion of strongly-anchored tumour neighbours.

cell_summary2 <- cell_summary %>%
  mutate(
    frac_ge_0_7 = frac_ge_0_7,
    highlight   = (max_p >= 0.9 | frac_ge_0_7 >= 0.8)
  )

cell_summary2 %>%
  filter(highlight) %>%
  arrange(desc(max_p)) %>%
  select(cell_line, max_p, frac_ge_0_7) %>%
  print(n = 30)

p_scatter <- ggplot(cell_summary2, aes(x = max_p, y = frac_ge_0_7)) +
  geom_point(
    aes(size = n_pairs, fill = max_p >= 0.7),
    shape = 21, colour = "black", alpha = 0.8
  ) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_vline(xintercept = 0.7, linetype = "dashed") +
  geom_text_repel(
    data = subset(cell_summary2, highlight),
    aes(label = cell_line),
    size = 3.5, max.overlaps = Inf, box.padding = 0.4,
    point.padding = 0.3, segment.size = 0.3, min.segment.length = 0
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                     labels = scales::percent_format(accuracy = 0)) +
  scale_fill_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#bdbdbd"),
    labels = c("FALSE" = "< 0.7", "TRUE" = ">= 0.7"),
    name   = "Max p_consensus"
  ) +
  scale_size_continuous(range = c(2, 6), name = "Number of tumour neighbours") +
  theme_dimred(base_size = 14) +
  labs(
    title = sprintf("Anchoring strength of DSMZ breast cancer cell lines (%s)", direction),
    subtitle = expression(
      x == max~p[consensus]~" per line; " ~
      y == "fraction of neighbours with " ~ p[consensus] >= 0.7
    ),
    x = expression(max~p[consensus]),
    y = "Tumour neighbours with strong consensus (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

plot_path_scatter <- file.path(cons_dir, sprintf("Fig_p_consensus_cell_scatter_%s.pdf", direction))
ggsave(plot_path_scatter, p_scatter, width = 7, height = 6)

cat("\nScatter plot (max_p vs frac_ge_0_7) saved to:\n  ", plot_path_scatter, "\n")
cat("\n[SUCCESS] p_consensus computation finished for gene set: ", gene_set, "\n", sep = "")