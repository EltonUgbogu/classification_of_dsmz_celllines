#!/usr/bin/env Rscript

# =============================================================================
# plot_cellline_metrics.R
# =============================================================================
#
# PURPOSE
# -------
# Generate cell-line–driven summary plots from:
#   - cellline_metrics_by_tumour_cohort.tsv
#   - cellline_metrics_overall.tsv
#
# OUTPUT
# ------
# Writes PDFs into --outdir (default: ./figures_cellline):
#   1) Fig_cellline_mean_rank_by_cohort.pdf
#   2) Fig_cellline_pct_in_top1_by_cohort.pdf
#   3) Fig_best_cellline_per_cohort_by_mean_score.pdf
#   4) Fig_best_cellline_per_cohort_by_pct_in_top1.pdf
#   5) Fig_cellline_score_vs_rank_scatter.pdf
#   6) Fig_cellline_representativeness_heatmap_mean_score.pdf
#
# DEPENDENCIES
# ------------
# data.table, ggplot2
#
# USAGE
# -----
# Rscript plot_cellline_metrics.R \
#   --by-cohort results/unsupervised/pan_cancer/tumour_mapping/cellline_metrics_by_tumour_cohort.tsv \
#   --overall  results/unsupervised/pan_cancer/tumour_mapping/cellline_metrics_overall.tsv \
#   --outdir   results/unsupervised/pan_cancer/tumour_mapping/figures
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("--by-cohort"), type="character", default=NULL,
              help="Path to cellline_metrics_by_tumour_cohort.tsv"),
  make_option(c("--overall"), type="character", default=NULL,
              help="Path to cellline_metrics_overall.tsv"),
  make_option(c("--outdir"), type="character", default="figures_cellline",
              help="Output directory for plots (default: figures_cellline)"),
  make_option(c("--top-n"), type="integer", default=25,
              help="Number of top cell lines to show per cohort in bar plots (default: 25)"),
  make_option(c("--cohorts"), type="character", default="BRCA,NBL,RBL",
              help="Comma-separated tumour cohorts to include (default: BRCA,NBL,RBL)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$`by-cohort`) || !file.exists(opt$`by-cohort`)) {
  stop("--by-cohort is required and must exist")
}
if (is.null(opt$overall) || !file.exists(opt$overall)) {
  stop("--overall is required and must exist")
}

outdir <- opt$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cohorts_keep <- trimws(strsplit(opt$cohorts, ",")[[1]])
top_n <- opt$`top-n`

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
save_pdf <- function(p, filename, w=10, h=6) {
  ggsave(filename = filename, plot = p, device = "pdf", width = w, height = h, units = "in")
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
byc <- fread(opt$`by-cohort`)
ov  <- fread(opt$overall)

# Basic schema checks
req_byc <- c("cell_line","cell_lineage","tum_lineage","mean_score","median_score",
             "mean_rank","median_rank","pct_in_top1","pct_in_topk")
miss_byc <- setdiff(req_byc, names(byc))
if (length(miss_byc) > 0) stop("by-cohort file missing columns: ", paste(miss_byc, collapse=", "))

req_ov <- c("cell_line","cell_lineage","mean_score","median_score","mean_rank","median_rank",
            "pct_in_top1","pct_in_topk")
miss_ov <- setdiff(req_ov, names(ov))
if (length(miss_ov) > 0) stop("overall file missing columns: ", paste(miss_ov, collapse=", "))

if (!"cell_line_display" %in% names(byc)) byc[, cell_line_display := cell_line]
if (!"cell_line_display" %in% names(ov)) ov[, cell_line_display := cell_line]

# Filter cohorts
byc <- byc[tum_lineage %in% cohorts_keep]

# Harmonize factor order
byc[, tum_lineage := factor(tum_lineage, levels = cohorts_keep)]
ov[, cell_lineage := as.factor(cell_lineage)]

# =============================================================================
# 1) Mean rank per cell line (per cohort)
# =============================================================================
# Show top-N most representative cell lines per cohort by mean_rank (lower is better)
rank_top <- byc[, .SD[order(mean_rank)][1:min(.N, top_n)], by = tum_lineage]
rank_top[, cell_line_display := factor(cell_line_display, levels = unique(cell_line_display[order(mean_rank)])), by=tum_lineage]

p1 <- ggplot(rank_top, aes(x = reorder(cell_line_display, -mean_rank), y = mean_rank, fill = cell_lineage)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ tum_lineage, scales = "free_y") +
  labs(
    title = sprintf("Cell-line mean rank within tumour cohorts (top %d by lowest mean rank)", top_n),
    x = "Cell line",
    y = "Mean rank (lower = better)",
    fill = "Cell line lineage"
  ) +
  theme_pub()

save_pdf(p1, file.path(outdir, "Fig_cellline_mean_rank_by_cohort.pdf"), w=12, h=8)

# =============================================================================
# 2) % tumours selecting cell line as Top-1 (per cohort)
# =============================================================================
top1_top <- byc[, .SD[order(-pct_in_top1, -mean_score)][1:min(.N, top_n)], by = tum_lineage]

p2 <- ggplot(top1_top, aes(x = reorder(cell_line_display, pct_in_top1), y = pct_in_top1, fill = cell_lineage)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ tum_lineage, scales = "free_y") +
  labs(
    title = sprintf("Cell-line selection frequency as Top-1 (top %d per cohort)", top_n),
    x = "Cell line",
    y = "Fraction of tumours selecting as Top-1",
    fill = "Cell line lineage"
  ) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100*x)) +
  theme_pub()

save_pdf(p2, file.path(outdir, "Fig_cellline_pct_in_top1_by_cohort.pdf"), w=12, h=8)

# =============================================================================
# 3) Best cell line per cohort (bar plot) — by highest mean_score
# =============================================================================
best_mean <- byc[order(tum_lineage, -mean_score, -pct_in_top1, median_rank)][
  , .SD[1], by = tum_lineage
]

p3 <- ggplot(best_mean, aes(x = tum_lineage, y = mean_score, fill = cell_line_display)) +
  geom_col() +
  geom_text(aes(label = cell_line_display), vjust = -0.4, size = 3) +
  labs(
    title = "Best cell line per tumour cohort (highest mean Spearman score)",
    x = "Tumour cohort",
    y = "Mean Spearman score",
    fill = "Cell line"
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_pdf(p3, file.path(outdir, "Fig_best_cellline_per_cohort_by_mean_score.pdf"), w=8, h=5)

# =============================================================================
# 4) Best cell line per cohort — by highest pct_in_top1
# =============================================================================
best_top1 <- byc[order(tum_lineage, -pct_in_top1, -mean_score, median_rank)][
  , .SD[1], by = tum_lineage
]

p4 <- ggplot(best_top1, aes(x = tum_lineage, y = pct_in_top1, fill = cell_line_display)) +
  geom_col() +
  geom_text(aes(label = cell_line_display), vjust = -0.4, size = 3) +
  labs(
    title = "Best cell line per tumour cohort (highest Top-1 selection frequency)",
    x = "Tumour cohort",
    y = "Fraction selected as Top-1",
    fill = "Cell line"
  ) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100*x)) +
  theme_pub() +
  theme(legend.position = "none")

save_pdf(p4, file.path(outdir, "Fig_best_cellline_per_cohort_by_pct_in_top1.pdf"), w=8, h=5)

# =============================================================================
# 5) Cell line score vs rank scatter (overall)
# =============================================================================
p5 <- ggplot(ov, aes(x = mean_rank, y = mean_score, shape = cell_lineage)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "Cell-line representativeness (overall): mean Spearman score vs mean rank",
    x = "Mean rank across tumours (lower = better)",
    y = "Mean Spearman score",
    shape = "Cell line lineage"
  ) +
  theme_pub()

save_pdf(p5, file.path(outdir, "Fig_cellline_score_vs_rank_scatter.pdf"), w=9, h=6)

# =============================================================================
# 6) Heatmap: cell line × cohort representativeness (mean_score)
# =============================================================================
# Build matrix-like long format for ggplot heatmap.
# To keep this readable, plot top-N cell lines per cohort by mean_score,
# then take union across cohorts.
top_union <- unique(byc[order(-mean_score), head(cell_line, top_n)])
hm <- byc[cell_line %in% top_union]

# Order cell lines by overall mean_score (within these)
cell_order <- hm[, .(mean_score_all = mean(mean_score, na.rm=TRUE)), by=cell_line_display][order(-mean_score_all)]$cell_line_display
hm[, cell_line_display := factor(cell_line_display, levels = cell_order)]

p6 <- ggplot(hm, aes(x = tum_lineage, y = cell_line_display, fill = mean_score)) +
  geom_tile() +
  labs(
    title = sprintf("Cell line × tumour cohort representativeness (mean score; union of top %d per cohort)", top_n),
    x = "Tumour cohort",
    y = "Cell line",
    fill = "Mean score"
  ) +
  theme_pub() +
  theme(panel.grid = element_blank())

save_pdf(p6, file.path(outdir, "Fig_cellline_representativeness_heatmap_mean_score.pdf"), w=9, h=10)

cat("\nDone. Plots saved to: ", outdir, "\n", sep="")
