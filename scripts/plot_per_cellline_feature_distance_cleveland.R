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

option_list <- list(
  make_option("--summary-tsv", type = "character",
              help = "p_consensus_cellline_direction_summary.tsv"),
  make_option("--winners-tsv", type = "character",
              help = "p_consensus_winners_by_frac_ge_thr.tsv"),
  make_option("--out-pdf", type = "character", help = "Output PDF path"),
  make_option("--out-png", type = "character", help = "Output PNG path"),
  make_option("--disease-label", type = "character",
              help = "Full disease label for visible plot text"),
  make_option("--threshold", type = "double", default = 0.7,
              help = "p_consensus threshold used for recovery score")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("summary-tsv", "winners-tsv", "out-pdf", "out-png", "disease-label")
missing <- required[vapply(required, function(x) {
  is.null(opt[[x]]) || !nzchar(opt[[x]])
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

parse_direction <- function(x) {
  x <- as.character(x)
  tibble(
    direction = x,
    feature_raw = str_replace(x, "_(euc|corr)$", ""),
    distance_raw = str_match(x, "_(euc|corr)$")[, 2]
  )
}

feature_labels <- c(
  Variance = "Variance",
  MAD = "Median absolute deviation",
  MeanAbsDev = "Mean absolute deviation",
  Entropy = "Entropy",
  PCA = "Principal component loadings",
  Spearman = "Mean absolute Spearman correlation",
  MX = "MX score",
  kTotal = "WGCNA total connectivity",
  HVG = "Highly variable genes",
  pam50 = "PAM50 gene set",
  PAM50 = "PAM50 gene set"
)

feature_palette <- c(
  "Variance" = "#0072B2",
  "Median absolute deviation" = "#E69F00",
  "Mean absolute deviation" = "#CC79A7",
  "Entropy" = "#6A3D9A",
  "Principal component loadings" = "#999999",
  "Mean absolute Spearman correlation" = "#56B4E9",
  "MX score" = "#000000",
  "WGCNA total connectivity" = "#D55E00",
  "Highly variable genes" = "#332288",
  "PAM50 gene set" = "#AA4499"
)

distance_labels <- c(
  euc = "Euclidean distance",
  corr = "Pearson correlation distance"
)

distance_shapes <- c(
  "Euclidean distance" = 21,
  "Pearson correlation distance" = 24
)

direction_info <- parse_direction(unique(summary_tbl$direction)) %>%
  mutate(
    feature_label = unname(feature_labels[feature_raw]),
    feature_label = if_else(is.na(feature_label), feature_raw, feature_label),
    distance_label = unname(distance_labels[distance_raw]),
    distance_label = if_else(is.na(distance_label), distance_raw, distance_label)
  )

unknown_dist <- direction_info %>% filter(is.na(distance_raw) | !(distance_raw %in% names(distance_labels)))
if (nrow(unknown_dist) > 0) {
  stop("Unsupported direction suffix in: ",
       paste(unique(unknown_dist$direction), collapse = ", "),
       call. = FALSE)
}

feature_order <- c(
  "Variance",
  "Median absolute deviation",
  "Mean absolute deviation",
  "Entropy",
  "Principal component loadings",
  "Mean absolute Spearman correlation",
  "MX score",
  "WGCNA total connectivity",
  "Highly variable genes",
  "PAM50 gene set"
)
feature_order <- c(feature_order[feature_order %in% direction_info$feature_label],
                   setdiff(sort(unique(direction_info$feature_label)), feature_order))

direction_order <- direction_info %>%
  mutate(
    feature_label = factor(feature_label, levels = feature_order),
    distance_label = factor(distance_label, levels = distance_labels)
  ) %>%
  arrange(feature_label, distance_label, direction) %>%
  pull(direction)

winner_long <- winners_tbl %>%
  mutate(
    best_dir_frac_ge_thr = as.character(best_dir_frac_ge_thr),
    best_frac_ge_thr = as.numeric(best_frac_ge_thr)
  ) %>%
  separate_rows(best_dir_frac_ge_thr, sep = ";") %>%
  mutate(best_dir_frac_ge_thr = str_trim(best_dir_frac_ge_thr)) %>%
  filter(nzchar(best_dir_frac_ge_thr)) %>%
  transmute(
    cell_line,
    direction = best_dir_frac_ge_thr,
    is_winner = TRUE,
    winner_score = best_frac_ge_thr
  )

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

order_tbl <- plot_tbl %>%
  group_by(cell_line) %>%
  summarise(
    best_frac_ge_thr = max(frac_ge_thr, na.rm = TRUE),
    n_recovered_directions = sum(frac_ge_thr >= opt$threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(best_frac_ge_thr), desc(n_recovered_directions), cell_line)

n_directions <- length(direction_order)
offset_tbl <- tibble(
  direction = direction_order,
  direction_index = seq_along(direction_order),
  direction_offset = if (n_directions <= 1) {
    0
  } else {
    ((seq_along(direction_order) - (n_directions + 1) / 2) / (n_directions - 1)) * 0.52
  }
)

cell_levels <- rev(order_tbl$cell_line)
plot_tbl <- plot_tbl %>%
  left_join(offset_tbl, by = "direction") %>%
  mutate(
    cell_line_factor = factor(cell_line, levels = cell_levels),
    y_base = as.numeric(cell_line_factor),
    y_plot = y_base + direction_offset,
    feature_label = factor(feature_label, levels = feature_order),
    distance_label = factor(distance_label, levels = distance_labels)
  )

missing_palette <- setdiff(levels(plot_tbl$feature_label), names(feature_palette))
if (length(missing_palette) > 0) {
  stop("No colour mapping for feature representation(s): ",
       paste(missing_palette, collapse = ", "), call. = FALSE)
}

plot_height <- max(5.8, 0.34 * length(cell_levels) + 3.4)
plot_width <- 11.2

threshold_text <- format(opt$threshold, trim = TRUE, scientific = FALSE)
subtitle_text <- sprintf("%s; threshold: p_consensus \u2265 %s", opt$`disease-label`, threshold_text)

p <- ggplot(plot_tbl, aes(x = frac_ge_thr, y = y_plot)) +
  geom_vline(xintercept = opt$threshold, colour = "#6F6F6F",
             linetype = "dashed", linewidth = 0.35) +
  geom_point(
    aes(fill = feature_label, shape = distance_label),
    colour = "#4D4D4D",
    size = 2.05,
    stroke = 0.28,
    alpha = 0.80
  ) +
  geom_point(
    data = plot_tbl %>% filter(is_winner),
    aes(fill = feature_label, shape = distance_label),
    colour = "black",
    size = 3.15,
    stroke = 0.90,
    alpha = 1
  ) +
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
    name = "Feature representation",
    values = feature_palette[feature_order],
    drop = FALSE
  ) +
  scale_shape_manual(
    name = "Distance metric",
    values = distance_shapes,
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      override.aes = list(shape = 21, size = 3.2, colour = "#4D4D4D", alpha = 1)
    ),
    shape = guide_legend(
      nrow = 1,
      override.aes = list(fill = "white", size = 3.2, colour = "#4D4D4D", alpha = 1)
    )
  ) +
  labs(
    title = "Consensus recovery across feature\u2013distance representations",
    subtitle = subtitle_text,
    x = "Fraction of tumour-neighbour pairs with p_consensus \u2265 threshold",
    y = "Cell line"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 10.5),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text.y = element_text(size = 8.8),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.6),
    legend.key.height = unit(0.42, "cm"),
    legend.spacing.y = unit(0.12, "cm"),
    plot.margin = margin(10, 20, 10, 16)
  )

dir.create(dirname(opt$`out-pdf`), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`out-png`), recursive = TRUE, showWarnings = FALSE)

ggsave(opt$`out-pdf`, p, width = plot_width, height = plot_height, device = cairo_pdf)
ggsave(opt$`out-png`, p, width = plot_width, height = plot_height, dpi = 320)

cat("[OK] Wrote ", opt$`out-pdf`, "\n", sep = "")
cat("[OK] Wrote ", opt$`out-png`, "\n", sep = "")
cat("[INFO] Cell lines: ", length(cell_levels), "\n", sep = "")
cat("[INFO] Directions: ", n_directions, "\n", sep = "")
cat("[INFO] Winner rows highlighted: ", sum(plot_tbl$is_winner), "\n", sep = "")
