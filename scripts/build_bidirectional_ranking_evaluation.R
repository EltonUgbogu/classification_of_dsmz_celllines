#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

required_arg <- function(opts, name) {
  value <- opts[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

metric_value <- function(dt, metric_name) {
  value <- dt[metric == metric_name, value]
  if (length(value) == 0) return(NA_real_)
  as.numeric(value[[1]])
}

write_count_matrix <- function(tab, path) {
  dt <- data.table(true_cancer_type = rownames(tab), as.data.frame.matrix(tab, stringsAsFactors = FALSE))
  fwrite(dt, path, sep = "\t")
  dt
}

write_fraction_matrix <- function(count_dt, labels, path) {
  frac <- copy(count_dt)
  row_totals <- rowSums(as.matrix(frac[, ..labels]))
  for (label in labels) {
    frac[[label]] <- ifelse(row_totals > 0, frac[[label]] / row_totals, NA_real_)
  }
  fwrite(frac, path, sep = "\t")
  frac
}

plot_ecdf <- function(ecdf_dt, title, path_pdf, path_png) {
  p <- ggplot(ecdf_dt, aes(
    x = rank,
    y = fraction_with_first_same_lineage_at_or_below_rank
  )) +
    geom_step(linewidth = 0.75, colour = "#1f4e79") +
    geom_point(
      data = ecdf_dt[rank %in% unique(c(1L, 10L, max(rank)))],
      size = 1.6,
      colour = "#1f4e79"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = function(x) sprintf("%.0f%%", 100 * x),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(
      breaks = pretty(ecdf_dt$rank, n = 8),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = title,
      x = "Rank threshold",
      y = "Samples with first same-lineage candidate"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  ggsave(path_pdf, p, width = 5.2, height = 3.8, units = "in")
  ggsave(path_png, p, width = 5.2, height = 3.8, units = "in", dpi = 300)
}

plot_confusion <- function(count_dt, frac_dt, labels, title, path_pdf, path_png) {
  count_long <- melt(
    count_dt,
    id.vars = "true_cancer_type",
    variable.name = "top1_candidate_cancer_type",
    value.name = "n"
  )
  frac_long <- melt(
    frac_dt,
    id.vars = "true_cancer_type",
    variable.name = "top1_candidate_cancer_type",
    value.name = "row_fraction"
  )
  plot_dt <- merge(count_long, frac_long,
                   by = c("true_cancer_type", "top1_candidate_cancer_type"))
  plot_dt[, label := sprintf("%d\n%.1f%%", n, 100 * row_fraction)]
  plot_dt[, true_cancer_type := factor(true_cancer_type, levels = labels)]
  plot_dt[, top1_candidate_cancer_type := factor(top1_candidate_cancer_type, levels = labels)]

  p <- ggplot(plot_dt, aes(
    x = top1_candidate_cancer_type,
    y = true_cancer_type,
    fill = row_fraction
  )) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 3.2, lineheight = 0.95) +
    scale_fill_gradient(
      low = "#f7fbff",
      high = "#08519c",
      limits = c(0, 1),
      labels = function(x) sprintf("%.0f%%", 100 * x)
    ) +
    coord_equal() +
    labs(
      title = title,
      x = "Top-1 candidate cancer type",
      y = "Sample cancer type",
      fill = "Row fraction"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave(path_pdf, p, width = 4.6, height = 4.0, units = "in")
  ggsave(path_png, p, width = 4.6, height = 4.0, units = "in", dpi = 300)
}

make_first_match_table <- function(score_mat, sample_labels, candidate_labels,
                                   direction, sample_display = NULL) {
  sample_ids <- rownames(score_mat)
  candidate_ids <- colnames(score_mat)
  rows <- vector("list", length(sample_ids))
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[[i]]
    scores <- as.numeric(score_mat[i, ])
    ord <- order(scores, decreasing = TRUE, na.last = NA)
    ranked_candidates <- candidate_ids[ord]
    ranked_candidate_types <- unname(candidate_labels[ranked_candidates])
    sample_type <- unname(sample_labels[sample_id])
    hit <- which(ranked_candidate_types == sample_type)
    first_rank <- if (length(hit) == 0) NA_integer_ else as.integer(hit[[1]])
    top1_candidate_id <- ranked_candidates[[1]]
    top1_candidate_type <- ranked_candidate_types[[1]]
    display <- sample_id
    if (!is.null(sample_display) && sample_id %in% names(sample_display)) {
      display <- unname(sample_display[sample_id])
    }
    rows[[i]] <- data.table(
      direction = direction,
      sample_id = sample_id,
      sample_label = display,
      sample_cancer_type = sample_type,
      first_same_lineage_rank = first_rank,
      top1_candidate_id = top1_candidate_id,
      top1_candidate_cancer_type = top1_candidate_type,
      top1_correct_cancer_type = isTRUE(top1_candidate_type == sample_type),
      n_ranked_candidates = length(candidate_ids),
      notes = "Cancer-type same-lineage evaluation; no exact one-to-one label is used."
    )
  }
  rbindlist(rows)
}

make_ecdf_table <- function(first_dt, max_rank, direction) {
  ranks <- seq_len(max_rank)
  first_ranks <- first_dt$first_same_lineage_rank
  data.table(
    direction = direction,
    rank = ranks,
    n_total = nrow(first_dt),
    n_with_first_same_lineage_at_or_below_rank =
      vapply(ranks, function(r) sum(!is.na(first_ranks) & first_ranks <= r), integer(1))
  )[, fraction_with_first_same_lineage_at_or_below_rank :=
       n_with_first_same_lineage_at_or_below_rank / n_total][]
}

make_confusion_outputs <- function(first_dt, prefix, title, outdir) {
  labels <- sort(unique(c(first_dt$sample_cancer_type, first_dt$top1_candidate_cancer_type)))
  tab <- table(
    factor(first_dt$sample_cancer_type, levels = labels),
    factor(first_dt$top1_candidate_cancer_type, levels = labels)
  )
  count_path <- file.path(outdir, paste0(prefix, "_top1_confusion_matrix_counts.tsv"))
  frac_path <- file.path(outdir, paste0(prefix, "_top1_confusion_matrix_row_fraction.tsv"))
  pdf_path <- file.path(outdir, paste0(prefix, "_top1_confusion_matrix.pdf"))
  png_path <- file.path(outdir, paste0(prefix, "_top1_confusion_matrix.png"))
  count_dt <- write_count_matrix(tab, count_path)
  frac_dt <- write_fraction_matrix(count_dt, labels, frac_path)
  plot_confusion(count_dt, frac_dt, labels, title, pdf_path, png_path)
  list(
    labels = labels,
    accuracy = sum(diag(tab)) / sum(tab),
    count_path = count_path,
    frac_path = frac_path,
    pdf_path = pdf_path,
    png_path = png_path
  )
}

register_manifest <- function(rows, name, direction, paths, file_types, claim, notes) {
  rbind(
    rows,
    data.table(
      evaluation_name = name,
      direction = direction,
      file_path = paths,
      file_type = file_types,
      workflow_managed = "yes",
      supports_result_claim = claim,
      notes = notes
    )
  )
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
tumour_score_rds <- required_arg(opts, "tumour-score-rds")
cellline_score_rds <- required_arg(opts, "cellline-score-rds")
tumour_rankings <- required_arg(opts, "tumour-rankings")
cellline_rankings <- required_arg(opts, "cellline-rankings")
tumour_summary_path <- required_arg(opts, "tumour-summary")
cellline_summary_path <- required_arg(opts, "cellline-summary")
tumour_metrics_path <- required_arg(opts, "tumour-metrics")
cellline_metrics_path <- required_arg(opts, "cellline-metrics")
outdir <- required_arg(opts, "outdir")

for (path in c(
  tumour_score_rds, cellline_score_rds, tumour_rankings, cellline_rankings,
  tumour_summary_path, cellline_summary_path, tumour_metrics_path,
  cellline_metrics_path
)) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tumour_scores <- readRDS(tumour_score_rds)
cellline_scores <- readRDS(cellline_score_rds)
tumour_summary <- fread(tumour_summary_path)
cellline_summary <- fread(cellline_summary_path)
tumour_metrics <- fread(tumour_metrics_path)
cellline_metrics <- fread(cellline_metrics_path)

tumour_type <- setNames(tumour_summary$tumour_lineage, tumour_summary$tumour_id)
cellline_type <- setNames(cellline_summary$cell_lineage, cellline_summary$cell_line)
cellline_display <- setNames(cellline_summary$cell_line_display, cellline_summary$cell_line)

if (!all(rownames(tumour_scores) %in% names(tumour_type))) {
  stop("Tumour score matrix contains samples without tumour cancer-type labels.")
}
if (!all(colnames(tumour_scores) %in% names(cellline_type))) {
  stop("Tumour score matrix contains candidates without cell-line cancer-type labels.")
}
if (!all(rownames(cellline_scores) %in% names(cellline_type))) {
  stop("Cell-line score matrix contains samples without cell-line cancer-type labels.")
}
if (!all(colnames(cellline_scores) %in% names(tumour_type))) {
  stop("Cell-line score matrix contains candidates without tumour cancer-type labels.")
}

tumour_first <- make_first_match_table(
  score_mat = tumour_scores,
  sample_labels = tumour_type,
  candidate_labels = cellline_type,
  direction = "tumour_to_cell_line"
)
cellline_first <- make_first_match_table(
  score_mat = cellline_scores,
  sample_labels = cellline_type,
  candidate_labels = tumour_type,
  direction = "cell_line_to_tumour",
  sample_display = cellline_display
)

tumour_first_path <- file.path(outdir, "tumour_to_cell_line_first_match_ranks.tsv")
cellline_first_path <- file.path(outdir, "cell_line_to_tumour_first_match_ranks.tsv")
fwrite(tumour_first, tumour_first_path, sep = "\t")
fwrite(cellline_first, cellline_first_path, sep = "\t")

tumour_ecdf <- make_ecdf_table(tumour_first, ncol(tumour_scores), "tumour_to_cell_line")
cellline_ecdf <- make_ecdf_table(cellline_first, ncol(cellline_scores), "cell_line_to_tumour")
tumour_ecdf_path <- file.path(outdir, "tumour_to_cell_line_first_same_lineage_ecdf.tsv")
cellline_ecdf_path <- file.path(outdir, "cell_line_to_tumour_first_same_lineage_ecdf.tsv")
fwrite(tumour_ecdf, tumour_ecdf_path, sep = "\t")
fwrite(cellline_ecdf, cellline_ecdf_path, sep = "\t")

tumour_ecdf_pdf <- file.path(outdir, "tumour_to_cell_line_first_same_lineage_ecdf.pdf")
tumour_ecdf_png <- file.path(outdir, "tumour_to_cell_line_first_same_lineage_ecdf.png")
cellline_ecdf_pdf <- file.path(outdir, "cell_line_to_tumour_first_same_lineage_ecdf.pdf")
cellline_ecdf_png <- file.path(outdir, "cell_line_to_tumour_first_same_lineage_ecdf.png")
plot_ecdf(tumour_ecdf, "Tumour-to-cell-line first same-lineage rank ECDF",
          tumour_ecdf_pdf, tumour_ecdf_png)
plot_ecdf(cellline_ecdf, "Cell-line-to-tumour first same-lineage rank ECDF",
          cellline_ecdf_pdf, cellline_ecdf_png)

tumour_conf <- make_confusion_outputs(
  tumour_first,
  "tumour_to_cell_line",
  "Tumour-to-cell-line top-1 cancer-type agreement",
  outdir
)
cellline_conf <- make_confusion_outputs(
  cellline_first,
  "cell_line_to_tumour",
  "Cell-line-to-tumour top-1 cancer-type agreement",
  outdir
)

crosscheck <- data.table(
  direction = c("tumour_to_cell_line", "cell_line_to_tumour"),
  reported_top1 = c(metric_value(tumour_metrics, "top1_accuracy"),
                    metric_value(cellline_metrics, "top1_accuracy")),
  confusion_matrix_diagonal_accuracy = c(tumour_conf$accuracy, cellline_conf$accuracy),
  reported_top10 = c(metric_value(tumour_metrics, "top10_accuracy"),
                     metric_value(cellline_metrics, "top10_accuracy")),
  reported_mrr = c(metric_value(tumour_metrics, "mrr"),
                   metric_value(cellline_metrics, "mrr")),
  notes = c(
    "Cancer-type top-1 agreement from current tumour-to-cell-line score matrix.",
    "Cancer-type top-1 agreement from current cell-line-to-tumour score matrix."
  )
)
crosscheck[, matches_reported_top1 :=
             abs(reported_top1 - confusion_matrix_diagonal_accuracy) < 1e-12]
setcolorder(crosscheck, c(
  "direction", "reported_top1", "confusion_matrix_diagonal_accuracy",
  "matches_reported_top1", "reported_top10", "reported_mrr", "notes"
))
crosscheck_path <- file.path(outdir, "ranking_metric_crosscheck.tsv")
fwrite(crosscheck, crosscheck_path, sep = "\t")

manifest <- data.table()
manifest <- register_manifest(
  manifest,
  "first_same_lineage_ecdf",
  "tumour_to_cell_line",
  c(tumour_ecdf_path, tumour_ecdf_pdf, tumour_ecdf_png, tumour_first_path),
  c("table", "figure_pdf", "figure_png", "table"),
  "Rank position of first same-lineage cell-line candidate for each tumour sample.",
  "Cancer-type same-lineage evaluation based on current tumour-to-cell-line scores."
)
manifest <- register_manifest(
  manifest,
  "first_same_lineage_ecdf",
  "cell_line_to_tumour",
  c(cellline_ecdf_path, cellline_ecdf_pdf, cellline_ecdf_png, cellline_first_path),
  c("table", "figure_pdf", "figure_png", "table"),
  "Rank position of first same-lineage tumour candidate for each cell-line group.",
  "Cancer-type same-lineage evaluation based on current cell-line-to-tumour scores."
)
manifest <- register_manifest(
  manifest,
  "top1_confusion_matrix",
  "tumour_to_cell_line",
  c(tumour_conf$count_path, tumour_conf$frac_path, tumour_conf$pdf_path, tumour_conf$png_path),
  c("table", "table", "figure_pdf", "figure_png"),
  "Top-1 cancer-type agreement for tumour-to-cell-line ranking.",
  "Rows are tumour cancer type; columns are top-1 cell-line candidate cancer type."
)
manifest <- register_manifest(
  manifest,
  "top1_confusion_matrix",
  "cell_line_to_tumour",
  c(cellline_conf$count_path, cellline_conf$frac_path, cellline_conf$pdf_path, cellline_conf$png_path),
  c("table", "table", "figure_pdf", "figure_png"),
  "Top-1 cancer-type agreement for cell-line-to-tumour ranking.",
  "Rows are cell-line cancer type; columns are top-1 tumour candidate cancer type."
)
manifest <- register_manifest(
  manifest,
  "metric_crosscheck",
  "bidirectional",
  crosscheck_path,
  "table",
  "Checks top-1 cancer-type agreement against reported ranking metrics.",
  "Uses current workflow metric tables as reported-value source."
)
manifest_path <- file.path(outdir, "ranking_evaluation_manifest.tsv")
fwrite(manifest, manifest_path, sep = "\t")

cat("Wrote bidirectional ranking evaluations to ", outdir, "\n", sep = "")
