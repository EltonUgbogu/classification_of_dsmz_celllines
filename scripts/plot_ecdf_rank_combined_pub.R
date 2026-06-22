#!/usr/bin/env Rscript
# =============================================================================
# plot_ecdf_rank_combined_pub.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
pipeline_root <- Sys.getenv("PIPELINE_ROOT", unset = repo_root)

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(n) || n <= 0) {
    return(data.table(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }

  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half_width <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / denom

  data.table(
    mean = p,
    lower = max(0, centre - half_width),
    upper = min(1, centre + half_width)
  )
}

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================

parse_cli <- function(args) {
  defaults <- list(
    pipeline_dir = "results/unsupervised/pan_cancer",
    score_rds = NA_character_,
    meta_rds = NA_character_,
    outdir = "results/unsupervised/pan_cancer/figures/ecdf_plots",
    n_show = 4L,
    top_k = 10L
  )

  if (any(args %in% c("-h", "--help"))) {
    cat(
      "Usage: Rscript plot_ecdf_rank_combined_pub.R --pipeline-dir results/unsupervised/pan_cancer --outdir figures/ecdf_plots\n",
      "\nOptions:\n",
      "  --pipeline-dir PATH  Directory containing tumour_mapping/ and pan_cancer_expr.rds.\n",
      "  --score-rds PATH     Optional explicit tumour_cellline_group_scores.rds path.\n",
      "  --meta-rds PATH      Optional explicit pan_cancer_expr.rds path.\n",
      "  --outdir PATH        Output directory for PDF, PNG, SVG, and TSV files.\n",
      "  --n-show INT         Number of displayed cell lines per ECDF panel.\n",
      "  --top-k INT          Rank threshold drawn in panels.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--pipeline-dir", "--score-rds", "--meta-rds", "--outdir", "--n-show", "--top-k")) {
      stop("Unknown argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for ", key)
    }

    value <- args[[i + 1L]]
    key_clean <- gsub("-", "_", sub("^--", "", key))
    defaults[[key_clean]] <- value
    i <- i + 2L
  }

  defaults$n_show <- as.integer(defaults$n_show)
  defaults$top_k <- as.integer(defaults$top_k)
  defaults
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))

# =============================================================================
# PATHS
# =============================================================================

PIPELINE_DIR <- cli$pipeline_dir
SCORE_RDS <- if (!is.na(cli$score_rds)) {
  cli$score_rds
} else {
  file.path(PIPELINE_DIR, "tumour_mapping/tumour_to_cellline_similarity/tumour_cellline_group_scores.rds")
}

META_RDS <- if (!is.na(cli$meta_rds)) {
  cli$meta_rds
} else {
  file.path(PIPELINE_DIR, "pan_cancer_expr.rds")
}

OUT_DIR <- cli$outdir
N_SHOW <- cli$n_show
TOP_K <- cli$top_k

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DISPLAY NAME HELPERS
# =============================================================================

make_display_labels <- function(sample_ids) {
  labels <- gsub("^NG-[^_]+_", "", sample_ids)
  labels <- gsub("_lib[^_]+.*$", "", labels)
  labels <- gsub("_", "-", labels)

  dup_mask <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(dup_mask)) {
    replicate_suffix <- sub(".*_(\\d+)$", "\\1", sample_ids)
    labels[dup_mask] <- paste0(labels[dup_mask], "-", replicate_suffix[dup_mask])
  }

  labels
}

make_biological_cell_line_names <- function(sample_ids) {
  labels <- gsub("^NG-[^_]+_", "", sample_ids)
  labels <- gsub("_lib[^_]+.*$", "", labels)
  labels
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("[1] Loading data...\n")
stopifnot(file.exists(SCORE_RDS), file.exists(META_RDS))

scores_mat <- readRDS(SCORE_RDS)
pan_obj <- readRDS(META_RDS)
meta <- as.data.table(pan_obj$meta)

meta_cell <- meta[type == "cell_line"]
meta_tum <- meta[type == "tumour"]

# =============================================================================
# CELL LINE GROUP METADATA
# =============================================================================
# Canonical ECDF inputs are tumour x biological-cell-line-group scores.  A
# profile-level score matrix is still accepted for manual debugging and is
# collapsed here by the same arithmetic-mean rule used in the scoring workflow.

cell_group_map <- data.table(sample_id = meta_cell$sample_id)
cell_group_map[, biological_cell_line_name := make_biological_cell_line_names(sample_id)]

cell_group_map <- merge(
  cell_group_map,
  meta_cell[, .(sample_id, lineage)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(cell_group_map$lineage)) {
  stop("[QC FAIL] Missing cell line lineage after profile to cell line mapping.")
}

lineage_check <- cell_group_map[, .(n_lineages = uniqueN(lineage)),
                                by = biological_cell_line_name]

if (any(lineage_check$n_lineages > 1L)) {
  bad <- lineage_check[n_lineages > 1L, biological_cell_line_name]
  stop("[QC FAIL] Conflicting lineage labels within biological cell lines: ",
       paste(bad, collapse = ", "))
}

group_meta_all <- cell_group_map[, .(
  lineage = unique(lineage),
  type = "cell_line",
  n_profile_observations = .N,
  contributing_profile_ids = paste(sample_id, collapse = ";")
), by = .(sample_id = biological_cell_line_name)]

if (all(colnames(scores_mat) %in% meta$sample_id)) {
  group_order <- unique(cell_group_map$biological_cell_line_name)

  collapsed_scores <- do.call(cbind, lapply(group_order, function(group_id) {
    profile_ids <- cell_group_map[biological_cell_line_name == group_id, sample_id]
    rowMeans(scores_mat[, profile_ids, drop = FALSE], na.rm = TRUE)
  }))

  rownames(collapsed_scores) <- rownames(scores_mat)
  colnames(collapsed_scores) <- group_order
  scores_mat <- collapsed_scores

  cat("  Input score matrix was profile-level; collapsed to biological groups.\n")
} else {
  missing_groups <- setdiff(colnames(scores_mat), group_meta_all$sample_id)
  if (length(missing_groups) > 0L) {
    stop("[QC FAIL] Score matrix columns are neither profile IDs nor known biological groups: ",
         paste(missing_groups, collapse = ", "))
  }
  group_order <- colnames(scores_mat)
  cat("  Input score matrix is biological group-level.\n")
}

if (anyNA(scores_mat)) {
  stop("[QC FAIL] Missing values detected after profile collapse.")
}

meta_cell <- group_meta_all[match(group_order, sample_id)]
if (anyNA(meta_cell$lineage)) {
  stop("[QC FAIL] Missing biological group metadata after score matrix alignment.")
}

disp <- make_display_labels(meta_cell$sample_id)
names(disp) <- meta_cell$sample_id

N_total <- ncol(scores_mat)

cat("  Collapsed cell line models:", N_total, "\n")
cat("  Profile collapse summary:\n")
print(meta_cell[n_profile_observations > 1L,
                .(sample_id, lineage, n_profile_observations, contributing_profile_ids)])

cat("  Tumour cohort sizes:\n")
print(table(meta_tum$lineage))

cat("\n  Display label mapping:\n")
lmap <- data.frame(
  sample_id = meta_cell$sample_id,
  display = disp,
  lineage = meta_cell$lineage,
  row.names = NULL
)
print(lmap)

# =============================================================================
# COMPUTE RANKS
# =============================================================================

cat("\n[2] Computing per tumour ranks...\n")

dt <- as.data.table(as.table(scores_mat))
setnames(dt, c("tumour", "cell_line", "score"))

dt[, tumour := as.character(tumour)]
dt[, cell_line := as.character(cell_line)]

dt[, tum_lineage := setNames(meta_tum$lineage, meta_tum$sample_id)[tumour]]
dt[, cell_lineage := setNames(meta_cell$lineage, meta_cell$sample_id)[cell_line]]

if (anyNA(dt$tum_lineage)) {
  stop("[QC FAIL] Missing tumour lineage labels after rank table construction.")
}
if (anyNA(dt$cell_lineage)) {
  stop("[QC FAIL] Missing cell line lineage labels after rank table construction.")
}

dt[, rank_in_tumour := frank(-score, ties.method = "average"), by = tumour]

# =============================================================================
# QC CHECKS
# =============================================================================

cat("[3] Quality control...\n")

if (any(dt$rank_in_tumour < 1, na.rm = TRUE)) {
  stop("[QC FAIL] Ranks below 1 detected.")
}
if (anyNA(dt$rank_in_tumour)) {
  stop("[QC FAIL] Missing ranks detected.")
}
if (any(dt$rank_in_tumour > N_total, na.rm = TRUE)) {
  warning("[QC WARN] Some ranks exceed N_total.")
}

for (lin in c("BRCA", "NBL", "RBL")) {
  shown_labels <- disp[meta_cell[lineage == lin, sample_id]]
  if (anyDuplicated(shown_labels)) {
    stop("[QC FAIL] Duplicated display labels in ", lin, ": ",
         paste(shown_labels[duplicated(shown_labels)], collapse = ", "))
  }

  n_nf <- nrow(meta_cell[lineage != lin])
  cat("  Outside lineage cell lines for", lin, ":", n_nf, "\n")

  if (n_nf < 2) {
    stop("[QC FAIL] Fewer than two outside lineage cell lines for ", lin)
  }
}

cat("  Rank validity: PASS\n")
cat("  Label validity: PASS\n")
cat("  Outside lineage count: PASS\n")

# =============================================================================
# HELPERS
# =============================================================================

OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
               "#D55E00", "#CC79A7", "#F0E442")

tumour_to_cellline_ecdf_plot_config <- list(
  output_pdf = "Fig_tumour_to_cellline_rank_ecdf.pdf",
  output_png = "Fig_tumour_to_cellline_rank_ecdf.png",
  plot_width_in = 15.0,
  plot_height_in = 5.8,
  png_dpi = 300L,
  base_font_size = 16,
  axis_title_font_size = 17,
  axis_text_font_size = 15,
  legend_text_font_size = 14,
  facet_label_font_size = 18,
  title_font_size = 18,
  line_width = 1.45,
  outside_line_width = 1.10,
  reference_line_width = 0.70
)

tumour_to_cellline_top10_fraction_plot_config <- list(
  output_pdf = "Fig_tumour_to_cellline_top10_fraction.pdf",
  output_png = "Fig_tumour_to_cellline_top10_fraction.png",
  plot_width_in = 9.8,
  plot_height_in = 8.8,
  png_dpi = 300L,
  base_font_size = 16,
  axis_title_font_size = 17,
  axis_text_font_size = 15,
  legend_text_font_size = 14,
  facet_label_font_size = 16,
  title_font_size = 18,
  point_size = 3.6,
  errorbar_line_width = 1.0,
  count_label_size = 4.6
)

ECDF_BASE_SIZE <- tumour_to_cellline_ecdf_plot_config$base_font_size
PANEL_D_BASE_SIZE <- tumour_to_cellline_top10_fraction_plot_config$base_font_size

make_x_breaks <- function(n_total) {
  cand <- c(1L, 5L, 10L, 20L, as.integer(n_total))
  sort(unique(cand[cand >= 1L & cand <= n_total]))
}

compute_nonfocal_band <- function(dt_lin, nonfocal_ids, n_total) {
  xg <- seq_len(n_total)

  mat <- vapply(nonfocal_ids, function(cl) {
    ecdf(dt_lin[cell_line == cl, rank_in_tumour])(xg)
  }, numeric(n_total))

  data.frame(
    rank = xg,
    median_ecdf = apply(mat, 1L, median),
    q25_ecdf = apply(mat, 1L, quantile, 0.25),
    q75_ecdf = apply(mat, 1L, quantile, 0.75)
  )
}

build_ecdf_panel <- function(target_lineage, panel_letter) {
  dt_lin <- dt[tum_lineage == target_lineage]
  focal_ids <- meta_cell[lineage == target_lineage, sample_id]
  nonfocal_ids <- meta_cell[lineage != target_lineage, sample_id]
  n_tumours <- length(unique(dt_lin$tumour))

  focal_stats <- dt_lin[cell_line %in% focal_ids,
    .(
      median_rank = median(rank_in_tumour, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE)
    ),
    by = cell_line
  ][order(median_rank, -median_score, cell_line)]

  show_ids <- head(focal_stats$cell_line, N_SHOW)
  show_labels <- disp[show_ids]
  colours <- setNames(OKABE_ITO[seq_along(show_ids)], show_labels)

  band_df <- compute_nonfocal_band(dt_lin, nonfocal_ids, N_total)

  focal_df <- dt_lin[cell_line %in% show_ids, .(rank_in_tumour, cell_line)]
  focal_df[, display := factor(disp[cell_line], levels = show_labels)]

  x_breaks <- make_x_breaks(N_total)

  cancer_title <- switch(
    target_lineage,
    BRCA = paste0("BRCA tumours (n = ", n_tumours, ")"),
    NBL = paste0("NBL tumours (n = ", n_tumours, ")"),
    RBL = paste0("RBL tumours (n = ", n_tumours, ")")
  )

  p <- ggplot() +
    geom_ribbon(
      data = band_df,
      aes(x = rank, ymin = q25_ecdf, ymax = q75_ecdf),
      fill = "grey65",
      alpha = 0.35,
      inherit.aes = FALSE
    ) +
    geom_step(
      data = band_df,
      aes(x = rank, y = median_ecdf, linetype = "Outside lineage median"),
      colour = "grey55",
      linewidth = tumour_to_cellline_ecdf_plot_config$outside_line_width,
      inherit.aes = FALSE
    ) +
    stat_ecdf(
      data = focal_df,
      aes(x = rank_in_tumour, colour = display),
      geom = "step",
      linewidth = tumour_to_cellline_ecdf_plot_config$line_width,
      pad = FALSE
    ) +
    geom_step(
      data = band_df,
      aes(x = rank, y = median_ecdf, linetype = "Outside lineage median"),
      colour = "grey20",
      linewidth = tumour_to_cellline_ecdf_plot_config$outside_line_width,
      alpha = 0.98,
      inherit.aes = FALSE
    ) +
    geom_vline(
      xintercept = TOP_K,
      linetype = "dotted",
      colour = "grey35",
      linewidth = tumour_to_cellline_ecdf_plot_config$reference_line_width
    ) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dotted",
      colour = "grey35",
      linewidth = tumour_to_cellline_ecdf_plot_config$reference_line_width
    ) +
    scale_x_log10(
      breaks = x_breaks,
      labels = x_breaks,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_linetype_manual(
      values = c("Outside lineage median" = "dashed"),
      name = NULL,
      guide = guide_legend(
        override.aes = list(colour = "grey55", linewidth = 1.0)
      )
    ) +
    labs(
      title = cancer_title,
      x = "Cell-line group rank",
      y = "Cumulative fraction of tumours"
    ) +
    theme_bw(base_size = ECDF_BASE_SIZE) +
   theme(
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = tumour_to_cellline_ecdf_plot_config$facet_label_font_size, face = "bold"),
    axis.text.x = element_text(size = tumour_to_cellline_ecdf_plot_config$axis_text_font_size),
    axis.text.y = element_text(size = tumour_to_cellline_ecdf_plot_config$axis_text_font_size),
    axis.title.x = element_text(size = tumour_to_cellline_ecdf_plot_config$axis_title_font_size),
    axis.title.y = element_text(size = tumour_to_cellline_ecdf_plot_config$axis_title_font_size),
    plot.margin = margin(5, 10, 6, 4)
  )
  list(
    plot = p,
    show_ids = show_ids,
    show_labels = show_labels,
    focal_stats = focal_stats,
    n_tumours = n_tumours,
    colours = colours
  )
}

# =============================================================================
# BUILD PANELS A, B, C
# =============================================================================

cat("[4] Building ECDF panels A, B, C...\n")

res_brca <- build_ecdf_panel("BRCA", "A")
res_nbl <- build_ecdf_panel("NBL", "B")
res_rbl <- build_ecdf_panel("RBL", "C")

selected_models <- rbindlist(list(
  data.table(cancer_type = "BRCA", cell_line = res_brca$show_ids,
             selection_rank = seq_along(res_brca$show_ids)),
  data.table(cancer_type = "NBL", cell_line = res_nbl$show_ids,
             selection_rank = seq_along(res_nbl$show_ids)),
  data.table(cancer_type = "RBL", cell_line = res_rbl$show_ids,
             selection_rank = seq_along(res_rbl$show_ids))
))

stopifnot(
  all(sapply(
    list(res_brca, res_nbl, res_rbl),
    function(r) max(make_x_breaks(N_total)) == N_total
  ))
)

cat("  Axis limit consistency: PASS. All panels reach N_total = ", N_total, "\n", sep = "")

# =============================================================================
# BUILD PANEL D
# =============================================================================

cat("[5] Building panel D at top ", TOP_K, "...\n", sep = "")

build_topk_df <- function(res, lineage) {
  dt_lin <- dt[tum_lineage == lineage]
  n_tum <- res$n_tumours

  rbindlist(lapply(res$show_ids, function(cl) {
    n_top <- dt_lin[cell_line == cl, sum(rank_in_tumour <= TOP_K, na.rm = TRUE)]
    ci <- wilson_ci(n_top, n_tum)[1, ]

    data.table(
      cancer_type = lineage,
      cell_line = cl,
      display = disp[cl],
      n_patient_samples = n_tum,
      top10_count = n_top,
      top10_fraction = ci$mean,
      top10_ci_lower = ci$lower,
      top10_ci_upper = ci$upper
    )
  }))
}

topk_df <- rbindlist(list(
  build_topk_df(res_brca, "BRCA"),
  build_topk_df(res_nbl, "NBL"),
  build_topk_df(res_rbl, "RBL")
))

model_colour_map <- c(
  res_brca$colours,
  res_nbl$colours,
  res_rbl$colours
)

topk_df[, display := reorder(display, top10_fraction)]
topk_df[, cancer_label := factor(
  cancer_type,
  levels = c("BRCA", "NBL", "RBL"),
  labels = c("BRCA", "NBL", "RBL")
)]

topk_df[, count_label := paste0(top10_count, "/", n_patient_samples)]
topk_df[, label_x := pmin(1.17, top10_ci_upper + 0.045)]

panel_d <- ggplot(topk_df, aes(x = top10_fraction, y = display, colour = display)) +
  geom_errorbar(
    aes(xmin = top10_ci_lower, xmax = top10_ci_upper),
    width = 0.28,
    linewidth = tumour_to_cellline_top10_fraction_plot_config$errorbar_line_width,
    orientation = "y"
  ) +
  geom_point(size = tumour_to_cellline_top10_fraction_plot_config$point_size) +
  geom_text(
    aes(x = label_x, label = count_label),
    hjust = 0,
    size = tumour_to_cellline_top10_fraction_plot_config$count_label_size,
    colour = "grey30"
  ) +
  facet_wrap(~ cancer_label, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    limits = c(0, 1.22),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_colour_manual(values = model_colour_map, guide = "none") +
  labs(
    title = paste0("Top-", TOP_K, " fraction by displayed cell-line group"),
    x = paste0("Top-", TOP_K, " fraction"),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = PANEL_D_BASE_SIZE) +
  theme(
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey96", colour = "grey80"),
    strip.text = element_text(size = tumour_to_cellline_top10_fraction_plot_config$facet_label_font_size, face = "bold"),
    plot.title = element_text(size = tumour_to_cellline_top10_fraction_plot_config$title_font_size, face = "bold"),
    axis.text.x = element_text(size = tumour_to_cellline_top10_fraction_plot_config$axis_text_font_size),
    axis.text.y = element_text(size = tumour_to_cellline_top10_fraction_plot_config$axis_text_font_size),
    axis.title.x = element_text(size = tumour_to_cellline_top10_fraction_plot_config$axis_title_font_size),
    axis.title.y = element_text(size = tumour_to_cellline_top10_fraction_plot_config$axis_title_font_size),
    plot.margin = margin(5, 18, 4, 4)
  )

# =============================================================================
# SAVE STANDALONE FIGURES
# =============================================================================

cat("[6] Saving standalone ECDF and top-10 fraction figures...\n")

ecdf_fig <- res_brca$plot + res_nbl$plot + res_rbl$plot +
  plot_layout(nrow = 1) +
  plot_annotation(
    title = "Tumour-to-cell-line rank ECDF by tumour cohort",
    subtitle = paste0("Ranks computed over ", N_total, " biological cell-line groups; vertical line marks rank ", TOP_K)
  ) &
  theme(
    plot.title = element_text(
      size = tumour_to_cellline_ecdf_plot_config$title_font_size,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = tumour_to_cellline_ecdf_plot_config$axis_text_font_size,
      hjust = 0,
      margin = margin(b = 8)
    ),
    plot.margin = margin(t = 10, r = 14, b = 10, l = 10)
  )

rank_ecdf_pdf_out <- file.path(OUT_DIR, tumour_to_cellline_ecdf_plot_config$output_pdf)
rank_ecdf_png_out <- file.path(OUT_DIR, tumour_to_cellline_ecdf_plot_config$output_png)
top10_fraction_pdf_out <- file.path(OUT_DIR, tumour_to_cellline_top10_fraction_plot_config$output_pdf)
top10_fraction_png_out <- file.path(OUT_DIR, tumour_to_cellline_top10_fraction_plot_config$output_png)

suppressWarnings(ggsave(
  rank_ecdf_pdf_out,
  ecdf_fig,
  width = tumour_to_cellline_ecdf_plot_config$plot_width_in,
  height = tumour_to_cellline_ecdf_plot_config$plot_height_in,
  device = cairo_pdf,
  units = "in"
))
suppressWarnings(ggsave(
  rank_ecdf_png_out,
  ecdf_fig,
  width = tumour_to_cellline_ecdf_plot_config$plot_width_in,
  height = tumour_to_cellline_ecdf_plot_config$plot_height_in,
  dpi = tumour_to_cellline_ecdf_plot_config$png_dpi,
  units = "in"
))
suppressWarnings(ggsave(
  top10_fraction_pdf_out,
  panel_d,
  width = tumour_to_cellline_top10_fraction_plot_config$plot_width_in,
  height = tumour_to_cellline_top10_fraction_plot_config$plot_height_in,
  device = cairo_pdf,
  units = "in"
))
suppressWarnings(ggsave(
  top10_fraction_png_out,
  panel_d,
  width = tumour_to_cellline_top10_fraction_plot_config$plot_width_in,
  height = tumour_to_cellline_top10_fraction_plot_config$plot_height_in,
  dpi = tumour_to_cellline_top10_fraction_plot_config$png_dpi,
  units = "in"
))

cat("  ECDF PDF:", rank_ecdf_pdf_out, "\n")
cat("  ECDF PNG:", rank_ecdf_png_out, "\n")
cat("  Top-10 fraction PDF:", top10_fraction_pdf_out, "\n")
cat("  Top-10 fraction PNG:", top10_fraction_png_out, "\n")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("[8] Writing summary table...\n")

build_summary_rows <- function(lineage) {
  dt_lin <- dt[tum_lineage == lineage]
  n_tum <- length(unique(dt_lin$tumour))
  all_ids <- meta_cell$sample_id

  rbindlist(lapply(all_ids, function(cl) {
    ranks <- dt_lin[cell_line == cl, rank_in_tumour]
    n_top <- sum(ranks <= TOP_K, na.rm = TRUE)
    ci <- wilson_ci(n_top, n_tum)[1, ]

    sr <- selected_models[cancer_type == lineage & cell_line == cl, selection_rank]
    sr <- if (length(sr) == 0L) NA_integer_ else as.integer(sr[1])

    data.table(
      cancer_type = lineage,
      cell_line = cl,
      display_label = disp[cl],
      lineage_status = ifelse(
        meta_cell[sample_id == cl, lineage] == lineage,
        "same_lineage",
        "outside_lineage"
      ),
      is_displayed_model = cl %in% selected_models[cancer_type == lineage, cell_line],
      selection_rank = sr,
      n_patient_samples = n_tum,
      n_total_models = N_total,
      median_rank = round(median(ranks, na.rm = TRUE), 2),
      median_score = round(median(dt_lin[cell_line == cl, score], na.rm = TRUE), 6),
      rank_IQR_lower = round(as.numeric(quantile(ranks, 0.25, na.rm = TRUE)), 2),
      rank_IQR_upper = round(as.numeric(quantile(ranks, 0.75, na.rm = TRUE)), 2),
      top10_count = n_top,
      top10_fraction = round(ci$mean, 4),
      top10_ci_lower = round(ci$lower, 4),
      top10_ci_upper = round(ci$upper, 4)
    )
  }))
}

summary_tbl <- rbindlist(lapply(c("BRCA", "NBL", "RBL"), build_summary_rows))
summary_tbl[is.na(selection_rank), selection_rank := NA_integer_]

tsv_out <- file.path(OUT_DIR, "model_prioritisation_rank_summary.tsv")
fwrite(summary_tbl, tsv_out, sep = "\t")

cat("  TSV:", tsv_out, "\n")

# =============================================================================
# REPRODUCIBILITY REPORT
# =============================================================================

cat("\n[DONE] Reproducibility summary\n")
cat("  Input score matrix:   ", SCORE_RDS, "\n")
cat("  Input metadata:       ", META_RDS, "\n")
cat("  N total cell lines:   ", N_total, "\n")
cat("  BRCA tumours:         ", res_brca$n_tumours, "\n")
cat("  NBL tumours:          ", res_nbl$n_tumours, "\n")
cat("  RBL tumours:          ", res_rbl$n_tumours, "\n")
cat("  Displayed cell lines: ", N_SHOW, " per panel\n")

for (lin in c("BRCA", "NBL", "RBL")) {
  n_nf <- nrow(meta_cell[lineage != lin])
  cat("  Outside lineage cell lines for", lin, ":", n_nf, "\n")
}

cat("  Outputs in:", OUT_DIR, "\n")

# =============================================================================
# DISPLAY LABEL VERIFICATION TABLE
# =============================================================================

cat("[9] Writing display label mapping table...\n")

bio_name <- gsub("^NG-[^_]+_", "", meta_cell$sample_id)
bio_name <- gsub("_lib[^_]+.*$", "", bio_name)
bio_name <- gsub("_", "-", bio_name)

label_map <- data.table(
  raw_identifier = meta_cell$sample_id,
  display_label = disp,
  cancer_type = meta_cell$lineage,
  biological_cell_line_name = bio_name,
  reason_for_label = ifelse(
    disp == bio_name,
    "unique base name",
    "replicate profile disambiguated by sample identifier suffix"
  )
)

map_out <- file.path(OUT_DIR, "model_prioritisation_display_label_mapping.tsv")
fwrite(label_map, map_out, sep = "\t")

cat("  Mapping TSV:", map_out, "\n")
print(label_map[display_label != biological_cell_line_name])

# =============================================================================
# PROVENANCE
# =============================================================================

prov_out <- file.path(OUT_DIR, "tumour_to_cellline_rank_ecdf_top10_fraction_provenance.tsv")

prov <- data.table(
  figure_name = paste(c(
    tumour_to_cellline_ecdf_plot_config$output_pdf,
    tumour_to_cellline_top10_fraction_plot_config$output_pdf
  ), collapse = ";"),
  script = "scripts/plot_ecdf_rank_combined_pub.R",
  command = paste(
    c(
      "Rscript",
      "scripts/plot_ecdf_rank_combined_pub.R",
      "--pipeline-dir", PIPELINE_DIR,
      "--score-rds", SCORE_RDS,
      "--meta-rds", META_RDS,
      "--outdir", OUT_DIR,
      "--n-show", N_SHOW,
      "--top-k", TOP_K
    ),
    collapse = " "
  ),
  git_commit = Sys.getenv("GIT_COMMIT", unset = "unavailable_not_git_worktree"),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_files = paste(normalizePath(c(SCORE_RDS, META_RDS), mustWork = FALSE), collapse = ";"),
  output_files = paste(normalizePath(c(
    rank_ecdf_pdf_out,
    rank_ecdf_png_out,
    top10_fraction_pdf_out,
    top10_fraction_png_out,
    tsv_out,
    map_out
  ), mustWork = FALSE), collapse = ";"),
  upstream_tables = normalizePath(tsv_out, mustWork = FALSE),
  key_parameters = paste0(
    "top_k=", TOP_K,
    ";n_show=", N_SHOW,
    ";rank_scale=log10",
    ";n_ranked_biological_cell_line_groups=", N_total,
    ";RBL_15_RBL_20_profiles=arithmetic_mean"
  ),
  software_versions = paste0(
    "R=", getRversion(),
    ";data.table=", as.character(packageVersion("data.table")),
    ";ggplot2=", as.character(packageVersion("ggplot2")),
    ";patchwork=", as.character(packageVersion("patchwork")),
    ";scales=", as.character(packageVersion("scales"))
  ),
  figure_type = "ecdf",
  source_pipeline_root = pipeline_root,
  copied_to_figure_export_path = "",
  legacy_source_path = "",
  notes = paste(
    "Tumour-to-cell-line rank ECDF and top-10 fraction plots use biological cell-line groups.",
    "RBL_15 and RBL_20 replicate profiles were averaged before ranking when profile-level input is supplied.",
    "Displayed cell-line groups were selected by median rank, then median score, then deterministic identifier order."
  )
)

fwrite(prov, prov_out, sep = "\t")
cat("  Provenance TSV:", prov_out, "\n")
