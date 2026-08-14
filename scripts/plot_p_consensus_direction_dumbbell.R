#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(stringr)
  library(scales)
})

option_list <- list(
  make_option("--brca", type = "character", help = "BRCA p_consensus_direction_summary.tsv"),
  make_option("--nbl", type = "character", help = "NBL p_consensus_direction_summary.tsv"),
  make_option("--rbl", type = "character", help = "RBL p_consensus_direction_summary.tsv"),
  make_option("--out_prefix", type = "character", help = "Output prefix for PDF, PNG, and SVG"),
  make_option("--caption", type = "character", help = "Caption text output path"),
  make_option("--threshold", type = "double", default = 0.7,
              help = "Strong-consensus p_consensus threshold [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
required_args <- c("brca", "nbl", "rbl", "out_prefix", "caption")
missing_args <- required_args[vapply(required_args, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required option(s): ", paste(missing_args, collapse = ", "))
}

metric_palette <- c(
  Euclidean = "#0072B2",
  Correlation = "#E69F00"
)
metric_shapes <- c(
  Euclidean = 16,
  Correlation = 17
)
segment_colour <- "#BDBDBD"
reference_colour <- "#6F6F6F"
label_close_threshold <- 5
cohort_label_map <- c(
  BRCA = "(a) Breast cancer",
  NBL = "(b) Neuroblastoma",
  RBL = "(c) Retinoblastoma"
)

cohort_inputs <- tibble::tibble(
  cohort = factor(c("BRCA", "NBL", "RBL"), levels = c("BRCA", "NBL", "RBL")),
  path = c(opt$brca, opt$nbl, opt$rbl)
)

read_summary <- function(cohort, path) {
  if (!file.exists(path)) {
    stop("Input TSV does not exist: ", path)
  }

  x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  required_cols <- c("direction", "frac_ge_thr")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("Input TSV missing required column(s) ", paste(missing_cols, collapse = ", "),
         ": ", path)
  }

  x %>%
    transmute(
      cohort = as.character(cohort),
      source_tsv = normalizePath(path),
      direction = as.character(direction),
      frac_ge_thr = as.numeric(frac_ge_thr)
    )
}

summary_tbl <- bind_rows(lapply(seq_len(nrow(cohort_inputs)), function(i) {
  read_summary(cohort_inputs$cohort[i], cohort_inputs$path[i])
}))

plot_tbl <- summary_tbl %>%
  filter(str_detect(direction, "_(euc|corr)$")) %>%
  mutate(
    feature = str_replace(direction, "_(euc|corr)$", ""),
    feature = if_else(feature == "pam50", "PAM50", feature),
    distance = if_else(str_ends(direction, "_euc"), "Euclidean", "Correlation"),
    value_pct = 100 * frac_ge_thr
  )

if (nrow(plot_tbl) == 0) {
  stop("No Euclidean or correlation p_consensus directions found in input TSVs.")
}
if (any(is.na(plot_tbl$frac_ge_thr))) {
  bad <- plot_tbl %>%
    filter(is.na(frac_ge_thr)) %>%
    distinct(cohort, direction)
  stop("Non-numeric or missing frac_ge_thr values for: ",
       paste(paste(bad$cohort, bad$direction, sep = ":"), collapse = ", "))
}

duplicates <- plot_tbl %>%
  count(cohort, feature, distance, name = "n") %>%
  filter(n > 1)
if (nrow(duplicates) > 0) {
  stop("Duplicate cohort/feature/distance rows found in summary TSVs: ",
       paste(paste(duplicates$cohort, duplicates$feature, duplicates$distance, sep = ":"),
             collapse = ", "))
}

feature_order <- plot_tbl %>%
  filter(feature != "PAM50") %>%
  group_by(feature) %>%
  summarise(grand_mean_strong_consensus_fraction = mean(frac_ge_thr), .groups = "drop") %>%
  arrange(desc(grand_mean_strong_consensus_fraction), feature) %>%
  pull(feature)
if ("PAM50" %in% unique(plot_tbl$feature)) {
  feature_order <- c(feature_order, "PAM50")
}

plot_tbl <- plot_tbl %>%
  mutate(
    cohort = factor(cohort, levels = c("BRCA", "NBL", "RBL")),
    feature = factor(feature, levels = rev(feature_order)),
    distance = factor(distance, levels = c("Euclidean", "Correlation"))
  )

segment_tbl <- plot_tbl %>%
  select(cohort, feature, distance, value_pct) %>%
  pivot_wider(names_from = distance, values_from = value_pct) %>%
  filter(!is.na(Euclidean), !is.na(Correlation))

label_tbl <- plot_tbl %>%
  select(cohort, feature, distance, value_pct) %>%
  pivot_wider(names_from = distance, values_from = value_pct) %>%
  filter(!is.na(Euclidean), !is.na(Correlation)) %>%
  mutate(
    euc_value = Euclidean,
    corr_value = Correlation,
    close_pair = abs(euc_value - corr_value) < label_close_threshold
  ) %>%
  pivot_longer(
    cols = c(Euclidean, Correlation),
    names_to = "distance",
    values_to = "value_pct"
  ) %>%
  mutate(
    distance = factor(distance, levels = c("Euclidean", "Correlation")),
    label = paste0(round(value_pct), "%"),
    raw_side = case_when(
      distance == "Euclidean" & euc_value <= corr_value ~ "left",
      distance == "Euclidean" & euc_value > corr_value ~ "right",
      distance == "Correlation" & euc_value <= corr_value ~ "right",
      TRUE ~ "left"
    ),
    x_nudge = case_when(
      raw_side == "left" & value_pct <= 3 ~ 0.6,
      raw_side == "right" & value_pct >= 97 ~ -0.6,
      raw_side == "left" ~ -1.4,
      TRUE ~ 1.4
    ),
    label_x = pmin(100, pmax(0, value_pct + x_nudge)),
    hjust = case_when(
      raw_side == "left" & value_pct <= 3 ~ 0,
      raw_side == "right" & value_pct >= 97 ~ 1,
      raw_side == "left" ~ 1,
      TRUE ~ 0
    ),
    vjust = case_when(
      close_pair & distance == "Euclidean" ~ 1.45,
      close_pair & distance == "Correlation" ~ -0.55,
      TRUE ~ 0.5
    )
  )

p <- ggplot(plot_tbl, aes(x = value_pct, y = feature)) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    linewidth = 0.45,
    colour = reference_colour
  ) +
  geom_segment(
    data = segment_tbl,
    aes(x = Euclidean, xend = Correlation, y = feature, yend = feature),
    inherit.aes = FALSE,
    linewidth = 0.65,
    colour = segment_colour,
    lineend = "round"
  ) +
  geom_point(
    aes(colour = distance, shape = distance),
    size = 2.8,
    stroke = 0.75
  ) +
  geom_text(
    data = label_tbl,
    aes(x = label_x, y = feature, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = 2.45,
    colour = "#222222"
  ) +
  facet_wrap(~ cohort, nrow = 1, labeller = labeller(cohort = cohort_label_map),
             scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_colour_manual(
    name = "Distance metric",
    values = metric_palette,
    breaks = c("Euclidean", "Correlation"),
    drop = FALSE
  ) +
  scale_shape_manual(
    name = "Distance metric",
    values = metric_shapes,
    breaks = c("Euclidean", "Correlation"),
    drop = FALSE
  ) +
  labs(
    x = "Mean p-consensus fraction across tumour-neighbourhood pairs (%)",
    y = "Feature set"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9.5),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#BDBDBD"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.3),
    axis.text.y = element_text(size = 9.5, colour = "#222222"),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.spacing.x = unit(1.1, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )

out_dir <- dirname(opt$out_prefix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$caption), recursive = TRUE, showWarnings = FALSE)

pdf_path <- paste0(opt$out_prefix, ".pdf")
png_path <- paste0(opt$out_prefix, ".png")
svg_path <- paste0(opt$out_prefix, ".svg")

ggsave(pdf_path, p, width = 10.5, height = 5.7, units = "in")
ggsave(png_path, p, width = 10.5, height = 5.7, units = "in", dpi = 300)
ggsave(svg_path, p, width = 10.5, height = 5.7, units = "in", device = grDevices::svg)

caption_text <- sprintf(
  paste(
    "Strong consensus was defined as p_consensus >= %.2f.",
    "Here p_consensus is the fraction of clustering methods that placed a",
    "cell-line--tumour pair in the same tumour-neighbourhood relation.",
    "Plotted values are 100 x frac_ge_thr from the validated",
    "p_consensus_direction_summary.tsv files; feature rows are ordered by",
    "the grand mean of frac_ge_thr across BRCA, NBL, RBL, Euclidean, and",
    "correlation entries present in the plotted data."
  ),
  opt$threshold
)

caption_lines <- c(
  caption_text,
  "",
  "Source TSVs:",
  paste0("BRCA: ", normalizePath(opt$brca)),
  paste0("NBL: ", normalizePath(opt$nbl)),
  paste0("RBL: ", normalizePath(opt$rbl)),
  "",
  "Palette:",
  paste0("Euclidean: ", metric_palette[["Euclidean"]], " circle"),
  paste0("Correlation: ", metric_palette[["Correlation"]], " triangle"),
  paste0("Connector/reference grey: ", segment_colour, " / ", reference_colour),
  "",
  "Feature order:",
  paste(feature_order, collapse = ", ")
)
writeLines(caption_lines, opt$caption)

cat("[INFO] Wrote ", pdf_path, "\n", sep = "")
cat("[INFO] Wrote ", png_path, "\n", sep = "")
cat("[INFO] Wrote ", svg_path, "\n", sep = "")
cat("[INFO] Wrote ", opt$caption, "\n", sep = "")
