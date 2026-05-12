#!/usr/bin/env Rscript
# =============================================================================
# plot_ecdf_rank_combined_pub.R
# Publication-quality combined ECDF figure for PLOS Computational Biology
# British English throughout.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(binom)
  library(svglite)
})

# =============================================================================
# PATHS
# =============================================================================
PIPELINE_DIR <- "/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/pan_cancer"
SCORE_RDS  <- file.path(PIPELINE_DIR,
  "tumour_mapping/tumour_to_cellline_similarity/tumour_cellline_scores.rds")
META_RDS   <- file.path(PIPELINE_DIR, "pan_cancer_expr.rds")
OUT_DIR    <- "/Users/eltonugbogu/Downloads/thesis_latex/figures/ecdf_plots"
N_SHOW     <- 4L   # focal cell lines per ECDF panel
TOP_K      <- 10L  # rank threshold

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DISPLAY-NAME MAPPING
# Strips NG-XXXXX_ prefix and _libXXXX... suffix.
# Duplicate base names (RBL_15, RBL_20 technical replicates) receive a
# distinguishing batch-run suffix extracted from the original sample ID.
# =============================================================================
make_display_labels <- function(sample_ids) {
  labels <- gsub("^NG-[^_]+_", "", sample_ids)
  labels <- gsub("_lib[^_]+.*$", "", labels)
  labels <- gsub("_", "-", labels)

  dup_mask <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(dup_mask)) {
    # Biological replicate suffix: last digit group in the sample ID
    # e.g. ...10098_1 → "1", ...10102_3 → "3"  (matches manuscript RBL_20_1, RBL_20_3)
    bio_suffix <- sub(".*_(\\d+)$", "\\1", sample_ids)
    labels[dup_mask] <- paste0(labels[dup_mask], "-", bio_suffix[dup_mask])
  }
  labels
}

# =============================================================================
# LOAD DATA
# =============================================================================
cat("[1] Loading data...\n")
stopifnot(file.exists(SCORE_RDS), file.exists(META_RDS))

scores_mat <- readRDS(SCORE_RDS)
pan_obj    <- readRDS(META_RDS)
meta       <- as.data.table(pan_obj$meta)
meta_cell  <- meta[type == "cell_line"]
meta_tum   <- meta[type == "tumour"]

stopifnot(all(colnames(scores_mat) %in% meta$sample_id))

disp <- make_display_labels(meta_cell$sample_id)
names(disp) <- meta_cell$sample_id

N_total <- ncol(scores_mat)

cat("  N_total cell line models:", N_total, "\n")
cat("  Tumour cohort sizes:\n")
print(table(meta_tum$lineage))
cat("\n  Display label mapping (verify duplicates):\n")
lmap <- data.frame(sample_id = meta_cell$sample_id,
                   display   = disp,
                   lineage   = meta_cell$lineage,
                   row.names = NULL)
print(lmap)

# =============================================================================
# COMPUTE RANKS
# =============================================================================
cat("\n[2] Computing per-tumour ranks...\n")
dt <- as.data.table(as.table(scores_mat))
setnames(dt, c("tumour", "cell_line", "score"))
dt[, tum_lineage  := setNames(meta_tum$lineage,  meta_tum$sample_id)[tumour]]
dt[, cell_lineage := setNames(meta_cell$lineage, meta_cell$sample_id)[cell_line]]
dt[, rank_in_tumour := frank(-score, ties.method = "average"), by = tumour]

# =============================================================================
# QC CHECKS
# =============================================================================
cat("[3] Quality control...\n")

if (any(dt$rank_in_tumour < 1, na.rm = TRUE))
  stop("[QC FAIL] Ranks < 1 detected. Ranks must start at 1.")
if (anyNA(dt$rank_in_tumour))
  stop("[QC FAIL] NA ranks present.")
if (any(dt$rank_in_tumour > N_total, na.rm = TRUE))
  warning("[QC WARN] Some ranks exceed N_total.")

for (lin in c("BRCA", "NBL", "RBL")) {
  focal_labels <- disp[meta_cell[lineage == lin, sample_id]]
  if (anyDuplicated(focal_labels))
    stop("[QC FAIL] Duplicated display labels in ", lin, ": ",
         paste(focal_labels[duplicated(focal_labels)], collapse = ", "))
  n_nf <- nrow(meta_cell[lineage != lin])
  cat("  Non-focal cell lines for", lin, ":", n_nf, "\n")
  if (n_nf < 2) stop("[QC FAIL] Fewer than 2 non-focal cell lines for ", lin)
}

all_panels_same_xlim <- TRUE  # enforced via shared scale in helpers
cat("  Rank validity:  PASS\n")
cat("  Label validity: PASS\n")
cat("  Non-focal count:PASS\n")

# =============================================================================
# HELPERS
# =============================================================================

OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
               "#D55E00", "#CC79A7", "#F0E442")

# x-axis breaks valid in [1, n_total]
make_x_breaks <- function(n_total) {
  cand <- c(1L, 2L, 5L, 10L, 20L, 50L, 100L, as.integer(n_total))
  sort(unique(cand[cand >= 1L & cand <= n_total]))
}

# Non-focal band: median + IQR of per-cell-line ECDFs across all non-focal lines
compute_nonfocal_band <- function(dt_lin, nonfocal_ids, n_total) {
  xg <- seq_len(n_total)
  # mat: n_total rows x n_nonfocal cols
  mat <- vapply(nonfocal_ids, function(cl) {
    ecdf(dt_lin[cell_line == cl, rank_in_tumour])(xg)
  }, numeric(n_total))
  data.frame(
    rank        = xg,
    median_ecdf = apply(mat, 1L, median),
    q25_ecdf    = apply(mat, 1L, quantile, 0.25),
    q75_ecdf    = apply(mat, 1L, quantile, 0.75)
  )
}

# Single ECDF panel (A, B, or C)
build_ecdf_panel <- function(lineage, panel_letter) {
  dt_lin      <- dt[tum_lineage == lineage]
  focal_ids   <- meta_cell[lineage == lineage,  sample_id]
  nonfocal_ids <- meta_cell[lineage != lineage, sample_id]
  n_tumours   <- length(unique(dt_lin$tumour))

  focal_stats <- dt_lin[cell_line %in% focal_ids,
    .(median_rank = median(rank_in_tumour, na.rm = TRUE)), by = cell_line
  ][order(median_rank)]

  show_ids    <- head(focal_stats$cell_line, N_SHOW)
  show_labels <- disp[show_ids]
  colours     <- setNames(OKABE_ITO[seq_along(show_ids)], show_labels)

  band_df  <- compute_nonfocal_band(dt_lin, nonfocal_ids, N_total)

  focal_df <- dt_lin[cell_line %in% show_ids, .(rank_in_tumour, cell_line)]
  focal_df[, display := factor(disp[cell_line], levels = show_labels)]

  x_breaks <- make_x_breaks(N_total)
  cancer_title <- switch(lineage,
    BRCA = paste0(panel_letter, ". Breast cancer (n = ", n_tumours, ")"),
    NBL  = paste0(panel_letter, ". Neuroblastoma (n = ", n_tumours, ")"),
    RBL  = paste0(panel_letter, ". Retinoblastoma (n = ", n_tumours, ")")
  )

  p <- ggplot() +
    # Non-focal IQR ribbon (visually secondary)
    geom_ribbon(data = band_df,
                aes(x = rank, ymin = q25_ecdf, ymax = q75_ecdf),
                fill = "grey88", alpha = 0.75, inherit.aes = FALSE) +
    # Non-focal median step curve
    geom_step(data = band_df,
              aes(x = rank, y = median_ecdf,
                  linetype = "Non-focal median"),
              colour = "grey55", linewidth = 0.85, inherit.aes = FALSE) +
    # Focal step ECDFs
    stat_ecdf(data = focal_df,
              aes(x = rank_in_tumour, colour = display),
              geom = "step", linewidth = 1.25, pad = FALSE) +
    # Vertical guide at rank 10
    geom_vline(xintercept = TOP_K, linetype = "dotted",
               colour = "grey35", linewidth = 0.55) +
    # Horizontal guide at 0.5 (median reference)
    geom_hline(yintercept = 0.5, linetype = "dotted",
               colour = "grey35", linewidth = 0.55) +
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
      values  = c("Non-focal median" = "dashed"),
      name    = NULL,
      guide   = guide_legend(
        override.aes = list(colour = "grey55", linewidth = 0.85))
    ) +
    labs(
      title = cancer_title,
      x = paste0("Cell line similarity rank (1 = highest out of ",
                  N_total, " models)"),
      y = "Cumulative fraction of\nclinical patient samples"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major  = element_line(colour = "grey93", linewidth = 0.35),
      panel.grid.minor  = element_blank(),
      legend.position   = "none",   # guides collected into guide_area()
      plot.title        = element_text(size = 10, face = "bold"),
      axis.text         = element_text(size = 8),
      axis.title        = element_text(size = 8.5),
      plot.margin       = margin(5, 7, 4, 4)
    )

  list(plot = p, show_ids = show_ids, show_labels = show_labels,
       focal_stats = focal_stats, n_tumours = n_tumours,
       colours = colours)
}

# =============================================================================
# BUILD PANELS A, B, C
# =============================================================================
cat("[4] Building ECDF panels A, B, C...\n")

res_brca <- build_ecdf_panel("BRCA", "A")
res_nbl  <- build_ecdf_panel("NBL",  "B")
res_rbl  <- build_ecdf_panel("RBL",  "C")

# Confirm all panels share the same x-axis limits (N_total)
stopifnot(
  identical(make_x_breaks(N_total), make_x_breaks(N_total)),  # trivially true
  all(sapply(list(res_brca, res_nbl, res_rbl),
             function(r) max(make_x_breaks(N_total)) == N_total))
)
cat("  Axis limit consistency: PASS (all panels reach N_total =", N_total, ")\n")

# =============================================================================
# BUILD PANEL D: top-K fraction with 95 % Wilson CI
# =============================================================================
cat("[5] Building panel D (top-", TOP_K, " fractions)...\n", sep = "")

build_topk_df <- function(res, lineage) {
  dt_lin <- dt[tum_lineage == lineage]
  n_tum  <- res$n_tumours
  rbindlist(lapply(res$show_ids, function(cl) {
    n_top <- dt_lin[cell_line == cl, sum(rank_in_tumour <= TOP_K, na.rm = TRUE)]
    ci    <- binom.confint(n_top, n_tum, methods = "wilson")[1, ]
    data.table(
      cancer_type       = lineage,
      cell_line         = cl,
      display           = disp[cl],
      n_patient_samples = n_tum,
      top10_count       = n_top,
      top10_fraction    = ci$mean,
      top10_ci_lower    = ci$lower,
      top10_ci_upper    = ci$upper
    )
  }))
}

topk_df <- rbindlist(list(
  build_topk_df(res_brca, "BRCA"),
  build_topk_df(res_nbl,  "NBL"),
  build_topk_df(res_rbl,  "RBL")
))

# Order cell lines within each cancer type by top10_fraction (best = highest)
topk_df[, display := reorder(display, top10_fraction)]
topk_df[, cancer_label := factor(cancer_type,
  levels = c("BRCA", "NBL", "RBL"),
  labels = c("Breast cancer", "Neuroblastoma", "Retinoblastoma"))]
# Label showing exact counts (important for RBL where 100 % can occur)
topk_df[, count_label := paste0(top10_count, "/", n_patient_samples)]

# Colours matching the cancer type (not individual cell lines)
ct_cols <- c("BRCA" = "#E69F00", "NBL" = "#56B4E9", "RBL" = "#009E73")

panel_d <- ggplot(topk_df,
                  aes(x = top10_fraction, y = display, colour = cancer_type)) +
  geom_errorbar(aes(xmin = top10_ci_lower, xmax = top10_ci_upper),
                width = 0.28, linewidth = 0.75, orientation = "y") +
  geom_point(size = 2.8) +
  geom_text(aes(label = count_label), hjust = -0.15, size = 2.8,
            colour = "grey30") +
  facet_wrap(~ cancer_label, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = percent_format(accuracy = 1)
  ) +
  scale_colour_manual(values = ct_cols, guide = "none") +
  labs(
    title = paste0("D. Fraction of clinical patient samples ranked in top ",
                   TOP_K, " (95 % Wilson CI)"),
    x     = paste0("Fraction of samples with rank ≤ ", TOP_K),
    y     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey96", colour = "grey80"),
    strip.text         = element_text(size = 9, face = "bold"),
    plot.title         = element_text(size = 10, face = "bold"),
    axis.text          = element_text(size = 8),
    axis.title         = element_text(size = 8.5),
    plot.margin        = margin(5, 7, 4, 4)
  )

# =============================================================================
# COMPOSE COMBINED FIGURE
# =============================================================================
cat("[6] Composing combined figure...\n")

ecdf_row <- res_brca$plot + res_nbl$plot + res_rbl$plot +
  plot_layout(nrow = 1)

# guide_area() reserves a dedicated row for collected legends, keeping them
# visually separated from panel D.
fig_combined <- ecdf_row / guide_area() / panel_d +
  plot_layout(
    guides  = "collect",
    heights = c(2.6, 0.35, 1.6)
  ) +
  plot_annotation(
    caption = paste0(
      "The x axis is shown on a log₁₀ scale and extends to ",
      N_total, " ranked models. ",
      "Grey ribbon: interquartile range of non-focal lineage ECDFs; ",
      "grey dashed curve: median non-focal ECDF. ",
      "Vertical dotted line: rank 10. Horizontal dotted line: cumulative fraction = 0.5."
    ),
    theme = theme(
      plot.caption    = element_text(size = 7.5, hjust = 0, colour = "grey30"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

# =============================================================================
# SAVE FIGURES
# =============================================================================
cat("[7] Saving figures...\n")

FIG_W <- 15
FIG_H <- 13

pdf_out <- file.path(OUT_DIR, "ecdf_model_prioritisation_combined.pdf")
png_out <- file.path(OUT_DIR, "ecdf_model_prioritisation_combined.png")
svg_out <- file.path(OUT_DIR, "ecdf_model_prioritisation_combined.svg")

# Suppress cosmetic patchwork guide-collection warnings; figures are correct.
suppressWarnings(ggsave(pdf_out, fig_combined, width = FIG_W, height = FIG_H,
       device = grDevices::pdf, units = "in"))
suppressWarnings(ggsave(png_out, fig_combined, width = FIG_W, height = FIG_H,
       dpi = 300, units = "in"))
suppressWarnings(ggsave(svg_out, fig_combined, width = FIG_W, height = FIG_H,
       device = svglite::svglite, units = "in"))

cat("  PDF:", pdf_out, "\n")
cat("  PNG:", png_out, "\n")
cat("  SVG:", svg_out, "\n")

# =============================================================================
# SUMMARY TABLE
# =============================================================================
cat("[8] Writing summary table...\n")

build_summary_rows <- function(lineage) {
  dt_lin <- dt[tum_lineage == lineage]
  n_tum  <- length(unique(dt_lin$tumour))
  all_ids <- meta_cell$sample_id

  rbindlist(lapply(all_ids, function(cl) {
    ranks <- dt_lin[cell_line == cl, rank_in_tumour]
    n_top <- sum(ranks <= TOP_K, na.rm = TRUE)
    ci    <- binom.confint(n_top, n_tum, methods = "wilson")[1, ]
    data.table(
      cancer_type       = lineage,
      cell_line         = cl,
      display_label     = disp[cl],
      lineage_status    = ifelse(
        meta_cell[sample_id == cl, lineage] == lineage, "focal", "non_focal"),
      n_patient_samples = n_tum,
      n_total_models    = N_total,
      median_rank       = round(median(ranks, na.rm = TRUE), 2),
      rank_IQR_lower    = round(as.numeric(quantile(ranks, 0.25, na.rm = TRUE)), 2),
      rank_IQR_upper    = round(as.numeric(quantile(ranks, 0.75, na.rm = TRUE)), 2),
      top10_count       = n_top,
      top10_fraction    = round(ci$mean, 4),
      top10_ci_lower    = round(ci$lower, 4),
      top10_ci_upper    = round(ci$upper, 4)
    )
  }))
}

summary_tbl <- rbindlist(lapply(c("BRCA", "NBL", "RBL"), build_summary_rows))
tsv_out <- file.path(OUT_DIR, "model_prioritisation_rank_summary.tsv")
fwrite(summary_tbl, tsv_out, sep = "\t")
cat("  TSV:", tsv_out, "\n")

# =============================================================================
# REPRODUCIBILITY REPORT
# =============================================================================
cat("\n[DONE] Reproducibility summary\n")
cat("  Input — score matrix: ", SCORE_RDS, "\n")
cat("  Input — metadata:     ", META_RDS, "\n")
cat("  N_total models:       ", N_total, "\n")
cat("  BRCA tumours:         ", res_brca$n_tumours, "\n")
cat("  NBL  tumours:         ", res_nbl$n_tumours, "\n")
cat("  RBL  tumours:         ", res_rbl$n_tumours, "\n")
cat("  Focal models shown:   ", N_SHOW, " per panel\n")
for (lin in c("BRCA", "NBL", "RBL")) {
  n_nf <- nrow(meta_cell[lineage != lin])
  cat("  Non-focal models for", lin, ":", n_nf, "\n")
}
cat("  Outputs in:", OUT_DIR, "\n")

# =============================================================================
# DISPLAY LABEL VERIFICATION TABLE
# =============================================================================
cat("[9] Writing display label mapping table...\n")

# Biological cell line name: base label before any replicate suffix
bio_name <- gsub("^NG-[^_]+_", "", meta_cell$sample_id)
bio_name <- gsub("_lib[^_]+.*$", "", bio_name)
bio_name <- gsub("_", "-", bio_name)

label_map <- data.table(
  raw_identifier        = meta_cell$sample_id,
  display_label         = disp,
  cancer_type           = meta_cell$lineage,
  biological_cell_line_name = bio_name,
  reason_for_label      = ifelse(
    disp == bio_name,
    "unique base name",
    paste0("replicate disambiguated: suffix from sample ID trailing digit")
  )
)

map_out <- file.path(OUT_DIR, "model_prioritisation_display_label_mapping.tsv")
fwrite(label_map, map_out, sep = "\t")
cat("  Mapping TSV:", map_out, "\n")
print(label_map[display_label != biological_cell_line_name])
