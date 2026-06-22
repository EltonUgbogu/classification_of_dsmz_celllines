#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      opts[[name]] <- TRUE
      i <- i + 1L
    } else {
      opts[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  opts
}

required_arg <- function(opts, name) {
  value <- opts[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required argument --", name)
  value
}

as_logical_column <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- tolower(as.character(x))
  ifelse(x_chr %in% c("true", "t", "1", "yes"), TRUE,
         ifelse(x_chr %in% c("false", "f", "0", "no"), FALSE, NA))
}

bootstrap_mean_ci <- function(x, n_resamples, conf_level = 0.95) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) == 1L || n_resamples <= 0L) return(rep(mean(x), 2L))
  boot <- replicate(n_resamples, mean(sample(x, size = length(x), replace = TRUE)))
  alpha <- (1 - conf_level) / 2
  as.numeric(stats::quantile(boot, probs = c(alpha, 1 - alpha), names = FALSE, na.rm = TRUE))
}

tumour_to_cellline_mrr_plot_config <- list(
  top_k = 10L,
  expected_biological_cell_line_groups = 56L,
  plot_width_in = 7.8,
  plot_height_in = 5.6,
  png_dpi = 300L,
  base_font_size = 15,
  axis_title_font_size = 17,
  axis_text_font_size = 14,
  legend_text_font_size = 14,
  title_font_size = 18,
  subtitle_font_size = 15,
  point_size = 2.6,
  point_alpha = 0.72,
  jitter_width = 0.16,
  mean_marker_shape = 23,
  mean_marker_size = 5.2,
  mean_marker_stroke = 0.9,
  show_confidence_intervals = TRUE,
  bootstrap_resamples = 2000L,
  bootstrap_seed = 20260618L,
  cohort_levels = c("BRCA", "NBL", "RBL"),
  cohort_colours = c(BRCA = "#4E79A7", NBL = "#F28E2B", RBL = "#59A14F")
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
rankings_path <- required_arg(opts, "rankings")
metrics_path <- required_arg(opts, "metrics")
by_tumour_path <- required_arg(opts, "by-tumour")
by_cohort_path <- required_arg(opts, "by-cohort")
pdf_path <- required_arg(opts, "pdf")
png_path <- required_arg(opts, "png")

if (!is.null(opts[["top-k"]])) tumour_to_cellline_mrr_plot_config$top_k <- as.integer(opts[["top-k"]])
if (!is.null(opts[["expected-groups"]])) {
  tumour_to_cellline_mrr_plot_config$expected_biological_cell_line_groups <- as.integer(opts[["expected-groups"]])
}

top_k <- tumour_to_cellline_mrr_plot_config$top_k
expected_groups <- tumour_to_cellline_mrr_plot_config$expected_biological_cell_line_groups

for (path in c(rankings_path, metrics_path)) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}
for (path in c(by_tumour_path, by_cohort_path, pdf_path, png_path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

# Load patient tumour-to-cell-line group rankings
rankings <- fread(rankings_path)
required_cols <- c("tumour", "rank", "cell_line", "cell_lineage", "tum_lineage")
missing_cols <- setdiff(required_cols, names(rankings))
if (length(missing_cols)) stop("Ranking table is missing columns: ", paste(missing_cols, collapse = ", "))
rankings[, rank := as.integer(rank)]
rankings[, tumour := as.character(tumour)]
rankings[, cell_line := as.character(cell_line)]
rankings[, cell_lineage := as.character(cell_lineage)]
rankings[, tum_lineage := as.character(tum_lineage)]
if ("is_correct" %in% names(rankings)) {
  rankings[, is_correct_logical := as_logical_column(is_correct)]
} else {
  rankings[, is_correct_logical := cell_lineage == tum_lineage]
}

candidate_universe_count <- uniqueN(rankings$cell_line)
if (candidate_universe_count != expected_groups) {
  stop(
    "Expected a ", expected_groups, " biological cell-line group candidate universe, but observed ",
    candidate_universe_count, " unique groups in the tumour-to-cell-line ranking table."
  )
}

ranked_group_counts <- rankings[, .(
  n_ranked_biological_cell_line_groups = expected_groups,
  n_displayed_top10_biological_cell_line_groups = uniqueN(cell_line),
  max_rank = max(rank, na.rm = TRUE),
  duplicate_rank_rows = .N - uniqueN(rank)
), by = tumour]
if (any(ranked_group_counts$max_rank != top_k | ranked_group_counts$n_displayed_top10_biological_cell_line_groups != top_k)) {
  bad <- ranked_group_counts[max_rank != top_k | n_displayed_top10_biological_cell_line_groups != top_k]
  stop(
    "Expected exactly ", top_k, " displayed group-level ranks per tumour for RR@10. Example tumours: ",
    paste(head(bad$tumour, 10L), collapse = ", ")
  )
}
if (any(ranked_group_counts$duplicate_rank_rows != 0L)) {
  bad <- ranked_group_counts[duplicate_rank_rows != 0L]
  stop("Ranking table has duplicate rank rows for tumours: ", paste(head(bad$tumour, 10L), collapse = ", "))
}

# Identify the first same-cancer-type cell-line group within the top 10
top_hits <- rankings[rank <= top_k & cell_lineage == tum_lineage, .(
  first_matching_cell_line_group_rank_at10 = min(rank, na.rm = TRUE)
), by = tumour]
top1 <- rankings[rank == 1L, .(
  tumour,
  tumour_cancer_type = tum_lineage,
  tumour_cohort = tum_lineage,
  top_ranked_cell_line_group = cell_line,
  top_ranked_cell_line_group_cancer_type = cell_lineage,
  rank1_correct = cell_lineage == tum_lineage
)]
if (nrow(top1) != uniqueN(rankings$tumour)) stop("Could not identify exactly one rank-1 candidate per tumour.")

by_tumour <- merge(top1, top_hits, by = "tumour", all.x = TRUE, sort = FALSE)
by_tumour <- merge(
  by_tumour,
  ranked_group_counts[, .(tumour, n_ranked_biological_cell_line_groups, n_displayed_top10_biological_cell_line_groups)],
  by = "tumour",
  all.x = TRUE,
  sort = FALSE
)

# Compute reciprocal rank at 10 for each patient tumour
by_tumour[, reciprocal_rank_at10 := fifelse(
  is.na(first_matching_cell_line_group_rank_at10),
  0,
  1 / first_matching_cell_line_group_rank_at10
)]
by_tumour[, first_matching_cell_line_group_rank_at10 := as.integer(first_matching_cell_line_group_rank_at10)]
setcolorder(by_tumour, c(
  "tumour",
  "tumour_cancer_type",
  "tumour_cohort",
  "first_matching_cell_line_group_rank_at10",
  "reciprocal_rank_at10",
  "top_ranked_cell_line_group",
  "top_ranked_cell_line_group_cancer_type",
  "rank1_correct",
  "n_ranked_biological_cell_line_groups",
  "n_displayed_top10_biological_cell_line_groups"
))
by_tumour[, tumour_cohort_sort := match(tumour_cohort, tumour_to_cellline_mrr_plot_config$cohort_levels)]
setorder(by_tumour, tumour_cohort_sort, tumour)
by_tumour[, tumour_cohort_sort := NULL]
fwrite(by_tumour, by_tumour_path, sep = "\t")

# Summarise MRR@10 by tumour cohort
set.seed(tumour_to_cellline_mrr_plot_config$bootstrap_seed)
summarise_rr <- function(dt, cohort_label) {
  rr <- dt$reciprocal_rank_at10
  ci <- if (tumour_to_cellline_mrr_plot_config$show_confidence_intervals) {
    bootstrap_mean_ci(rr, tumour_to_cellline_mrr_plot_config$bootstrap_resamples)
  } else {
    c(NA_real_, NA_real_)
  }
  data.table(
    cohort = cohort_label,
    n_tumours = length(rr),
    mean_mrr_at10 = mean(rr, na.rm = TRUE),
    median_rr_at10 = median(rr, na.rm = TRUE),
    fraction_rr_at10_gt0 = mean(rr > 0, na.rm = TRUE),
    mean_mrr_at10_ci_low = ci[[1]],
    mean_mrr_at10_ci_high = ci[[2]],
    bootstrap_resamples = tumour_to_cellline_mrr_plot_config$bootstrap_resamples
  )
}
cohort_summary <- rbindlist(lapply(
  tumour_to_cellline_mrr_plot_config$cohort_levels,
  function(cohort) summarise_rr(by_tumour[tumour_cohort == cohort], cohort)
))
cohort_summary <- rbind(cohort_summary, summarise_rr(by_tumour, "Overall"), fill = TRUE)

metrics <- fread(metrics_path)
reported_mrr <- if (all(c("metric", "value") %in% names(metrics)) && any(metrics$metric == "mrr")) {
  as.numeric(metrics[metric == "mrr", value][[1]])
} else {
  NA_real_
}
cohort_summary[, reported_existing_mrr := NA_real_]
cohort_summary[cohort == "Overall", reported_existing_mrr := reported_mrr]
cohort_summary[, delta_vs_existing_mrr := mean_mrr_at10 - reported_existing_mrr]
cohort_summary[, matches_existing_reported_mrr := fifelse(
  is.na(delta_vs_existing_mrr),
  NA,
  abs(delta_vs_existing_mrr) < 1e-12
)]
fwrite(cohort_summary, by_cohort_path, sep = "\t")

plot_dt <- copy(by_tumour)
plot_dt[, tumour_cohort := factor(tumour_cohort, levels = tumour_to_cellline_mrr_plot_config$cohort_levels)]
summary_dt <- cohort_summary[cohort %in% tumour_to_cellline_mrr_plot_config$cohort_levels]
summary_dt[, cohort := factor(cohort, levels = tumour_to_cellline_mrr_plot_config$cohort_levels)]

# Plot tumour-level reciprocal rank at 10 by cohort
p <- ggplot(plot_dt, aes(x = tumour_cohort, y = reciprocal_rank_at10, colour = tumour_cohort)) +
  geom_jitter(
    width = tumour_to_cellline_mrr_plot_config$jitter_width,
    height = 0,
    size = tumour_to_cellline_mrr_plot_config$point_size,
    alpha = tumour_to_cellline_mrr_plot_config$point_alpha,
    show.legend = FALSE
  ) +
  {
    if (tumour_to_cellline_mrr_plot_config$show_confidence_intervals) {
      geom_errorbar(
        data = summary_dt,
        aes(x = cohort, ymin = mean_mrr_at10_ci_low, ymax = mean_mrr_at10_ci_high),
        inherit.aes = FALSE,
        width = 0.18,
        linewidth = 0.85,
        colour = "#222222"
      )
    }
  } +
  geom_point(
    data = summary_dt,
    aes(x = cohort, y = mean_mrr_at10),
    inherit.aes = FALSE,
    shape = tumour_to_cellline_mrr_plot_config$mean_marker_shape,
    size = tumour_to_cellline_mrr_plot_config$mean_marker_size,
    stroke = tumour_to_cellline_mrr_plot_config$mean_marker_stroke,
    colour = "#222222",
    fill = "white"
  ) +
  scale_colour_manual(values = tumour_to_cellline_mrr_plot_config$cohort_colours, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Tumour-to-cell-line MRR@10 by cancer type",
    subtitle = "Each point is one patient tumour; 0 indicates no same-cancer-type cell-line group within the top 10",
    x = "Tumour cohort",
    y = "Reciprocal rank at 10"
  ) +
  theme_bw(base_size = tumour_to_cellline_mrr_plot_config$base_font_size) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = tumour_to_cellline_mrr_plot_config$title_font_size,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = tumour_to_cellline_mrr_plot_config$subtitle_font_size,
      hjust = 0,
      margin = margin(b = 8)
    ),
    axis.title = element_text(
      face = "bold",
      size = tumour_to_cellline_mrr_plot_config$axis_title_font_size
    ),
    axis.text = element_text(size = tumour_to_cellline_mrr_plot_config$axis_text_font_size),
    legend.text = element_text(size = tumour_to_cellline_mrr_plot_config$legend_text_font_size),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 12, r = 18, b = 12, l = 16)
  )

ggsave(
  pdf_path,
  plot = p,
  width = tumour_to_cellline_mrr_plot_config$plot_width_in,
  height = tumour_to_cellline_mrr_plot_config$plot_height_in,
  units = "in"
)
ggsave(
  png_path,
  plot = p,
  width = tumour_to_cellline_mrr_plot_config$plot_width_in,
  height = tumour_to_cellline_mrr_plot_config$plot_height_in,
  units = "in",
  dpi = tumour_to_cellline_mrr_plot_config$png_dpi
)

cat("Wrote tumour-to-cell-line MRR@10 per-tumour table: ", by_tumour_path, "\n", sep = "")
cat("Wrote tumour-to-cell-line MRR@10 cohort table: ", by_cohort_path, "\n", sep = "")
cat("Wrote tumour-to-cell-line MRR@10 plot PDF: ", pdf_path, "\n", sep = "")
cat("Wrote tumour-to-cell-line MRR@10 plot PNG: ", png_path, "\n", sep = "")
if (is.finite(reported_mrr)) {
  cat(sprintf(
    "Computed Overall MRR@10 %.15f; existing reported MRR %.15f; delta %.15f\n",
    cohort_summary[cohort == "Overall", mean_mrr_at10],
    reported_mrr,
    cohort_summary[cohort == "Overall", delta_vs_existing_mrr]
  ))
}
