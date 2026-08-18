#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

# ==============================================================================
# plot_per_cellline_feature_distance_cleveland.R
# ------------------------------------------------------------------------------
# Generate per-cell-line Cleveland plots from p_consensus summary tables.
#
# The plot summarises, for each cell line, the fraction of cell-line--tumour
# neighbourhood pairs passing a consensus threshold across feature-selection
# methods and distance metrics. Distance metrics are shown as facets, feature-
# selection methods are encoded by point colour, and the winning representation
# for each cell line is highlighted by point type.
#
# Inputs:
#   --summary-tsv
#       Table with one row per cell_line x feature--distance direction.
#       Required columns: direction, cell_line, frac_ge_thr.
#
#   --winners-tsv
#       Table identifying the best feature--distance representation(s) per
#       cell line.
#       Required columns: cell_line, best_frac_ge_thr, best_dir_frac_ge_thr.
#
# Outputs:
#   --out-pdf
#       Vector PDF plot.
#
#   --out-png
#       Raster PNG plot.
#
# Cohort behaviour:
#   The script is driven by the directions present in the summary table.
#   BRCA may include PAM50 when PAM50_euc/PAM50_corr are present.
#   NBL and RBL usually contain only the unsupervised feature-selection methods.
#
# Validation:
#   --validate-only performs input parsing, column checks, direction parsing,
#   winner matching, and plot-dimension checks without writing outputs.
# ==============================================================================


# ------------------------------------------------------------------------------
# Command-line interface
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--summary-tsv", type = "character",
              help = "p_consensus_cellline_direction_summary.tsv"),
  make_option("--winners-tsv", type = "character",
              help = "p_consensus_winners_by_frac_ge_thr.tsv"),
  make_option("--out-pdf", type = "character", default = NA,
              help = "Output PDF path (ignored with --validate-only)"),
  make_option("--out-png", type = "character", default = NA,
              help = "Output PNG path (ignored with --validate-only)"),
  make_option("--disease-label", type = "character",
              help = "Full disease label for visible plot text"),
  # Configuration-owned: the workflow passes
  # patient_referenced_graph.p_consensus_threshold. No default is declared here,
  # so a figure can never be drawn against a threshold this script chose.
  make_option("--threshold", type = "double", default = NULL,
              help = paste("p_consensus threshold used for the plotted fraction;",
                           "supplied by the workflow from",
                           "patient_referenced_graph.p_consensus_threshold")),
  make_option("--validate-only", action = "store_true", default = FALSE,
              help = "Validate inputs and print checks without writing outputs")
)

opt <- parse_args(OptionParser(option_list = option_list))

# The p-consensus fraction threshold is configuration-owned and supplied by the
# workflow; there is no built-in default here, so a figure can never be produced
# against a different threshold from the one the graph stage applies.
if (is.null(opt$threshold) || !is.finite(opt$threshold) ||
    opt$threshold <= 0 || opt$threshold > 1) {
  stop(
    "A p-consensus fraction threshold in (0, 1] must be supplied via --threshold; ",
    "the workflow passes patient_referenced_graph.p_consensus_threshold."
  )
}


# ------------------------------------------------------------------------------
# Required argument and file validation
# ------------------------------------------------------------------------------

required <- c("summary-tsv", "winners-tsv", "disease-label")
if (!isTRUE(opt$`validate-only`)) {
  required <- c(required, "out-pdf", "out-png")
}

missing <- required[vapply(required, function(x) {
  is.null(opt[[x]]) || is.na(opt[[x]]) || !nzchar(opt[[x]])
}, logical(1))]

if (length(missing) > 0) {
  stop("Missing required option(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

if (!file.exists(opt$`summary-tsv`)) {
  stop("Summary TSV does not exist: ", opt$`summary-tsv`, call. = FALSE)
}

if (!file.exists(opt$`winners-tsv`)) {
  stop("Winners TSV does not exist: ", opt$`winners-tsv`, call. = FALSE)
}


# ------------------------------------------------------------------------------
# Input loading and schema checks
# ------------------------------------------------------------------------------

summary_tbl <- read_tsv(opt$`summary-tsv`, show_col_types = FALSE, progress = FALSE)
winners_tbl <- read_tsv(opt$`winners-tsv`, show_col_types = FALSE, progress = FALSE)

summary_required <- c("direction", "cell_line", "frac_ge_thr")
summary_missing <- setdiff(summary_required, names(summary_tbl))

if (length(summary_missing) > 0) {
  stop("Summary TSV missing required column(s): ",
       paste(summary_missing, collapse = ", "), call. = FALSE)
}

winners_required <- c("cell_line", "best_frac_ge_thr", "best_dir_frac_ge_thr")
winners_missing <- setdiff(winners_required, names(winners_tbl))

if (length(winners_missing) > 0) {
  stop("Winners TSV missing required column(s): ",
       paste(winners_missing, collapse = ", "), call. = FALSE)
}


# ------------------------------------------------------------------------------
# Direction parsing
# ------------------------------------------------------------------------------

# Direction names encode both the feature-selection method and the distance
# metric using the suffix pattern <feature>_<euc|corr>.
parse_direction <- function(x) {
  x <- as.character(x)
  tibble(
    direction = x,
    feature_raw = str_replace(x, "_(euc|corr)$", ""),
    distance_raw = str_match(x, "_(euc|corr)$")[, 2]
  )
}


# ------------------------------------------------------------------------------
# Feature and distance display labels
# ------------------------------------------------------------------------------

# Feature labels preserve the terminology used in the thesis methods and figures.
# Spearman and WGCNA retain distinct connectivity labels because they represent
# different connectivity-based feature-selection procedures.
feature_labels <- c(
  Variance = "Variance",
  MAD = "Median absolute deviation",
  MeanAbsDev = "Mean absolute deviation",
  Entropy = "Entropy",
  PCA = "Principal component loadings",
  Spearman = "Spearman connectivity",
  MX = "MX",
  kTotal = "WGCNA total connectivity",
  HVG = "HVG residual variance",
  pam50 = "PAM50 gene set",
  PAM50 = "PAM50 gene set"
)

# Feature order controls legend order and the deterministic visual offsets.
# The order separates visually similar colours where possible.
feature_order <- c(
  "Variance",
  "Median absolute deviation",
  "Entropy",
  "Principal component loadings",
  "Spearman connectivity",
  "WGCNA total connectivity",
  "Mean absolute deviation",
  "MX",
  "HVG residual variance",
  "PAM50 gene set"
)

feature_palette <- c(
  "Variance" = "#0072B2",
  "Median absolute deviation" = "#E69F00",
  "Mean absolute deviation" = "#CC79A7",
  "Entropy" = "#6A3D9A",
  "Principal component loadings" = "#999999",
  "Spearman connectivity" = "#56B4E9",
  "MX" = "#000000",
  "WGCNA total connectivity" = "#D55E00",
  "HVG residual variance" = "#332288",
  "PAM50 gene set" = "#AA4499"
)

distance_labels <- c(
  euc = "Euclidean distance",
  corr = "Pearson-correlation distance"
)


# ------------------------------------------------------------------------------
# Direction metadata table
# ------------------------------------------------------------------------------

direction_info <- parse_direction(unique(summary_tbl$direction)) %>%
  mutate(
    feature_label = unname(feature_labels[feature_raw]),
    feature_label = if_else(is.na(feature_label), feature_raw, feature_label),
    distance_label = unname(distance_labels[distance_raw]),
    distance_label = if_else(is.na(distance_label), distance_raw, distance_label)
  )

unknown_dist <- direction_info %>%
  filter(is.na(distance_raw) | !(distance_raw %in% names(distance_labels)))

if (nrow(unknown_dist) > 0) {
  stop("Unsupported direction suffix in: ",
       paste(unique(unknown_dist$direction), collapse = ", "), call. = FALSE)
}

# Only feature-selection methods present in the input table are included in the
# plot and legend. This keeps PAM50 cohort-specific without special-case logic.
active_feature_order <- c(
  feature_order[feature_order %in% direction_info$feature_label],
  setdiff(sort(unique(direction_info$feature_label)), feature_order)
)

n_features <- length(active_feature_order)

if (n_features < length(feature_order) &&
    "PAM50 gene set" %in% feature_order &&
    !("PAM50 gene set" %in% active_feature_order)) {
  message("[INFO] PAM50 gene set absent from this cohort's directions.")
}


# ------------------------------------------------------------------------------
# Winner-direction expansion
# ------------------------------------------------------------------------------

# Winner tables may contain multiple tied best directions separated by semicolons.
# Each winning feature--distance representation is expanded to one row.
winner_long <- winners_tbl %>%
  mutate(
    best_dir_frac_ge_thr = as.character(best_dir_frac_ge_thr),
    best_frac_ge_thr = as.numeric(best_frac_ge_thr)
  ) %>%
  separate_rows(best_dir_frac_ge_thr, sep = ";") %>%
  mutate(best_dir_frac_ge_thr = str_trim(best_dir_frac_ge_thr)) %>%
  filter(nzchar(best_dir_frac_ge_thr)) %>%
  transmute(cell_line, direction = best_dir_frac_ge_thr, is_winner = TRUE)


# ------------------------------------------------------------------------------
# Plot table construction
# ------------------------------------------------------------------------------

plot_tbl <- summary_tbl %>%
  mutate(
    direction = as.character(direction),
    cell_line = as.character(cell_line),
    frac_ge_thr = as.numeric(frac_ge_thr)
  ) %>%
  filter(!is.na(frac_ge_thr), frac_ge_thr >= 0, frac_ge_thr <= 1) %>%
  left_join(direction_info, by = "direction") %>%
  left_join(winner_long, by = c("cell_line", "direction")) %>%
  mutate(is_winner = if_else(is.na(is_winner), FALSE, is_winner))

if (nrow(plot_tbl) == 0) {
  stop("No plottable rows after filtering frac_ge_thr to [0, 1].", call. = FALSE)
}


# ------------------------------------------------------------------------------
# Soft checks
# ------------------------------------------------------------------------------

# Missing winner rows or reduced direction counts are reported as warnings because
# the plot can still be generated. Missing required files or columns remain hard
# errors earlier in the script.
cell_lines_in_summary <- unique(plot_tbl$cell_line)
cell_lines_in_winners <- unique(winners_tbl$cell_line)

no_winner_row <- setdiff(cell_lines_in_summary, cell_lines_in_winners)

if (length(no_winner_row) > 0) {
  warning("Cell line(s) with no entry in the winners table: ",
          paste(no_winner_row, collapse = ", "),
          call. = FALSE, immediate. = TRUE)
}

dir_counts <- plot_tbl %>%
  count(cell_line, name = "n_directions")

thin_cell_lines <- dir_counts %>%
  filter(n_directions < 2 * n_features - 2)

if (nrow(thin_cell_lines) > 0) {
  warning("Cell line(s) with fewer than expected directions: ",
          paste(thin_cell_lines$cell_line, collapse = ", "),
          call. = FALSE, immediate. = TRUE)
}


# ------------------------------------------------------------------------------
# Cell-line ordering and per-cell-line range calculation
# ------------------------------------------------------------------------------

# Cell lines are ordered by their largest p-consensus fraction, then by
# the number of directions passing the threshold. The grey range segment spans
# the minimum and maximum fraction observed across all feature-selection methods
# and both distance metrics for the same cell line.
order_tbl <- plot_tbl %>%
  group_by(cell_line) %>%
  summarise(
    best_frac_ge_thr = max(frac_ge_thr, na.rm = TRUE),
    worst_frac_ge_thr = min(frac_ge_thr, na.rm = TRUE),
    n_recovered_directions = sum(frac_ge_thr >= opt$threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(best_frac_ge_thr), desc(n_recovered_directions), cell_line)

cell_levels <- rev(order_tbl$cell_line)


# ------------------------------------------------------------------------------
# Deterministic visual offsets
# ------------------------------------------------------------------------------

# Small deterministic offsets separate overlapping feature-selection points
# within each cell-line row. The x-offset is intentionally small and only improves
# visibility of tied or near-tied values; it is not used for any calculation.
offset_tbl <- tibble(
  feature_label = active_feature_order,
  feature_offset = if (n_features <= 1) {
    0
  } else {
    ((seq_along(active_feature_order) - (n_features + 1) / 2) /
       (n_features - 1)) * 0.42
  }
) %>%
  mutate(x_nudge = feature_offset * 0.012)

plot_tbl <- plot_tbl %>%
  mutate(feature_label = factor(feature_label, levels = active_feature_order)) %>%
  left_join(offset_tbl, by = "feature_label") %>%
  mutate(
    cell_line_factor = factor(cell_line, levels = cell_levels),
    y_base = as.numeric(cell_line_factor),
    y_plot = y_base + feature_offset,
    x_plot = pmin(1, pmax(0, frac_ge_thr + x_nudge)),
    distance_label = factor(distance_label, levels = distance_labels),
    point_role = factor(
      if_else(is_winner, "Winning representation", "Other representations"),
      levels = c("Other representations", "Winning representation")
    )
  )

missing_palette <- setdiff(levels(plot_tbl$feature_label), names(feature_palette))

if (length(missing_palette) > 0) {
  stop("No colour mapping for feature-selection method(s): ",
       paste(missing_palette, collapse = ", "), call. = FALSE)
}

range_tbl <- order_tbl %>%
  transmute(
    cell_line,
    cell_line_factor = factor(cell_line, levels = cell_levels),
    y_base = as.numeric(cell_line_factor),
    x_min = worst_frac_ge_thr,
    x_max = best_frac_ge_thr
  )


# ------------------------------------------------------------------------------
# Plot dimensions and labels
# ------------------------------------------------------------------------------

n_cell_lines <- length(cell_levels)
plot_height <- max(5.8, 0.34 * n_cell_lines + 3.4)
plot_width <- 12.6

threshold_text <- format(opt$threshold, trim = TRUE, scientific = FALSE)

subtitle_text <- sprintf(
  "%s; threshold: p_consensus \u2265 %s",
  opt$`disease-label`,
  threshold_text
)

x_lab <- sprintf(
  "Fraction of cell-line\u2013tumour neighbourhood pairs with p_consensus \u2265 %s",
  threshold_text
)


# ------------------------------------------------------------------------------
# Validate-only check output
# ------------------------------------------------------------------------------

if (isTRUE(opt$`validate-only`)) {
  cat("[VALIDATE-ONLY] ", opt$`disease-label`, "\n", sep = "")
  cat("  cell lines               : ", n_cell_lines, "\n", sep = "")
  cat("  feature-selection methods: ", n_features, " (",
      paste(active_feature_order, collapse = ", "), ")\n", sep = "")
  cat("  distance metrics         : ",
      paste(levels(plot_tbl$distance_label), collapse = ", "), "\n", sep = "")
  cat("  rows in summary_tbl      : ", nrow(summary_tbl), "\n", sep = "")
  cat("  winner rows matched      : ", sum(plot_tbl$is_winner),
      " of ", n_cell_lines, " cell lines\n", sep = "")
  cat("  planned output height    : ", round(plot_height, 2),
      "in; width: ", plot_width, "in\n", sep = "")
  cat("[OK] Validation complete -- no files written.\n")
  quit(status = 0)
}


# ------------------------------------------------------------------------------
# Plot construction
# ------------------------------------------------------------------------------

p <- ggplot(plot_tbl, aes(x = x_plot, y = y_plot)) +
  geom_vline(
    xintercept = opt$threshold,
    colour = "#6F6F6F",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_segment(
    data = range_tbl,
    aes(x = x_min, xend = x_max, y = y_base, yend = y_base),
    inherit.aes = FALSE,
    colour = "#BFBFBF",
    linewidth = 2.1,
    alpha = 0.55,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = feature_label, size = point_role, colour = point_role),
    shape = 21,
    stroke = 0.32
  ) +
  facet_wrap(~ distance_label, ncol = 2) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_along(cell_levels),
    labels = cell_levels,
    expand = expansion(add = c(0.55, 0.75))
  ) +
  scale_fill_manual(
    name = "Feature-selection method",
    values = feature_palette[active_feature_order],
    drop = FALSE,
    guide = guide_legend(
      ncol = 2,
      order = 1,
      override.aes = list(shape = 21, size = 3.0, colour = "#4D4D4D")
    )
  ) +
  scale_size_manual(
    name = "Point type",
    values = c(
      "Other representations" = 2.15,
      "Winning representation" = 3.25
    ),
    guide = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "grey75",
        colour = c("#4D4D4D", "black"),
        alpha = c(0.85, 1)
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "Other representations" = "#4D4D4D",
      "Winning representation" = "black"
    ),
    guide = "none"
  ) +
  scale_alpha_identity() +
  labs(
    title = "P-consensus fraction across feature\u2013distance representations",
    subtitle = subtitle_text,
    x = x_lab,
    y = "Cell line",
    caption = paste0(
      "Grey bar: range of the p-consensus fraction across all feature-selection methods ",
      "and both distance metrics for that cell line. Larger, black-outlined point: winning feature\u2013distance ",
      "representation for that cell line. Small deterministic offsets separate overlapping points."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 7.6, colour = "#595959", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text.y = element_text(size = 8.8),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.25),
    panel.spacing.x = unit(14, "pt"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.6),
    legend.key.height = unit(0.42, "cm"),
    plot.margin = margin(10, 22, 10, 16)
  )


# ------------------------------------------------------------------------------
# Output writing
# ------------------------------------------------------------------------------

dir.create(dirname(opt$`out-pdf`), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`out-png`), recursive = TRUE, showWarnings = FALSE)

ggsave(opt$`out-pdf`, p, width = plot_width, height = plot_height, device = cairo_pdf)
ggsave(opt$`out-png`, p, width = plot_width, height = plot_height, dpi = 320)

cat("[OK] Wrote ", opt$`out-pdf`, "\n", sep = "")
cat("[OK] Wrote ", opt$`out-png`, "\n", sep = "")
cat("[INFO] Cell lines: ", n_cell_lines, "\n", sep = "")
cat("[INFO] Feature-selection methods: ", n_features, "\n", sep = "")
cat("[INFO] Winner rows highlighted: ", sum(plot_tbl$is_winner), "\n", sep = "")