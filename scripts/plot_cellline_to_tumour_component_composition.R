#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
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

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)))
  } else {
    fread(path)
  }
}

same_cancer_component_plot_config <- list(
  output_pdf = "Fig_cellline_to_tumour_same_cancer_top50_component_composition.pdf",
  output_png = "Fig_cellline_to_tumour_same_cancer_top50_component_composition.png",
  plot_width_in = 12.0,
  plot_height_in = 14.5,
  png_dpi = 300L,
  base_font_size = 15,
  axis_title_font_size = 17,
  axis_text_font_size = 13,
  legend_title_font_size = 14,
  legend_text_font_size = 13,
  facet_label_font_size = 15,
  title_font_size = 18,
  subtitle_font_size = 14,
  bar_width = 0.86,
  legend_columns = 2L,
  top_n = 50L
)

all_top50_component_plot_config <- list(
  output_pdf = "Fig_cellline_to_tumour_all_top50_component_composition.pdf",
  output_png = "Fig_cellline_to_tumour_all_top50_component_composition.png",
  plot_width_in = 12.0,
  plot_height_in = 14.5,
  png_dpi = 300L,
  base_font_size = 15,
  axis_title_font_size = 17,
  axis_text_font_size = 13,
  legend_title_font_size = 14,
  legend_text_font_size = 13,
  facet_label_font_size = 15,
  title_font_size = 18,
  subtitle_font_size = 14,
  bar_width = 0.86,
  legend_columns = 2L,
  top_n = 50L
)

lineage_levels <- c("BRCA", "NBL", "RBL")

normalise_scores <- function(scores) {
  required <- c("cell_line_group", "cell_lineage", "tumour", "score")
  if (!"cell_line_group" %in% names(scores) && "cell_line" %in% names(scores)) {
    scores[, cell_line_group := cell_line]
  }
  if (!"tum_lineage" %in% names(scores) && "tumour_lineage" %in% names(scores)) {
    setnames(scores, "tumour_lineage", "tum_lineage")
  }
  required <- c(required, "tum_lineage")
  missing_cols <- setdiff(required, names(scores))
  if (length(missing_cols)) {
    stop("Score table is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  scores[, cell_line_group := as.character(cell_line_group)]
  scores[, cell_lineage := as.character(cell_lineage)]
  scores[, tumour := as.character(tumour)]
  scores[, tum_lineage := as.character(tum_lineage)]
  scores[, score := as.numeric(score)]
  scores <- scores[is.finite(score)]
  scores
}

normalise_components <- function(components) {
  if ("sample" %in% names(components) && !"tumour" %in% names(components)) {
    setnames(components, "sample", "tumour")
  }
  if ("component" %in% names(components) && !"component_id" %in% names(components)) {
    setnames(components, "component", "component_id")
  }
  if ("comp_size" %in% names(components) && !"component_size" %in% names(components)) {
    setnames(components, "comp_size", "component_size")
  }
  missing_cols <- setdiff(c("tumour", "component_id"), names(components))
  if (length(missing_cols)) {
    stop("Component table is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!"component_size" %in% names(components)) components[, component_size := NA_integer_]
  components[, tumour := as.character(tumour)]
  components[, component_id := as.character(component_id)]
  components
}

rank_cellline_to_tumour_scores <- function(scores, components) {
  # Load cell-line-to-tumour scores and assign cell-line-to-tumour ranks
  ranked <- copy(scores)
  setorder(ranked, cell_line_group, -score, tumour)
  ranked[, rank_cellline_to_tumour := seq_len(.N), by = cell_line_group]
  ranked[, same_cancer_type := cell_lineage == tum_lineage]
  ranked <- merge(
    ranked,
    unique(components[, .(tumour, component_id, component_size)]),
    by = "tumour",
    all.x = TRUE
  )
  if (anyNA(ranked$component_id)) {
    missing_n <- uniqueN(ranked[is.na(component_id)]$tumour)
    stop("Missing component assignment for ", missing_n, " ranked tumours.")
  }
  setnames(ranked, c("cell_lineage", "tum_lineage"),
           c("cell_line_cancer_type", "tumour_cancer_type"))
  ranked
}

build_component_metadata <- function(ranked) {
  # Summarise patient consensus-network component metadata from ranked tumours
  tumour_meta <- unique(ranked[, .(
    tumour,
    tumour_cancer_type,
    component_id,
    component_size
  )])
  component_meta <- tumour_meta[, {
    lineage_counts <- table(tumour_cancer_type)
    dominant <- names(sort(lineage_counts, decreasing = TRUE))[1]
    purity <- max(as.numeric(lineage_counts)) / sum(as.numeric(lineage_counts))
    list(
      component_label = paste0("Component ", component_id[1]),
      component_dominant_lineage = dominant,
      component_size = component_size[which(!is.na(component_size))[1]],
      component_tumour_count = .N,
      component_purity = purity
    )
  }, by = component_id]
  component_meta[is.na(component_size), component_size := component_tumour_count]
  component_meta
}

build_same_cancer_component_composition <- function(ranked, component_meta, top_n) {
  # Identify same-cancer-type ranked tumours before selecting retained top-ranked tumours
  same_candidates <- ranked[same_cancer_type == TRUE]
  available <- same_candidates[, .(
    same_cancer_type_tumours_available = .N
  ), by = .(cell_line_group, cell_line_cancer_type)]

  same_top <- same_candidates[order(cell_line_group, rank_cellline_to_tumour)][
    ,
    head(.SD, min(top_n, .N)),
    by = .(cell_line_group, cell_line_cancer_type)
  ]
  retained <- same_top[, .(
    n_same_cancer_type_tumours_retained = .N
  ), by = .(cell_line_group, cell_line_cancer_type)]

  comp <- same_top[, .(
    n_retained_tumours_in_component = .N
  ), by = .(cell_line_group, cell_line_cancer_type, component_id)]
  comp <- merge(comp, available, by = c("cell_line_group", "cell_line_cancer_type"))
  comp <- merge(comp, retained, by = c("cell_line_group", "cell_line_cancer_type"))
  comp <- merge(comp, component_meta, by = "component_id", all.x = TRUE)
  comp[, top_n_cutoff_used := top_n]
  comp[, tumour_cancer_type_used_for_filtering := cell_line_cancer_type]
  comp[, fraction_retained_tumours_in_component :=
         n_retained_tumours_in_component / n_same_cancer_type_tumours_retained]
  setcolorder(comp, c(
    "cell_line_group",
    "cell_line_cancer_type",
    "tumour_cancer_type_used_for_filtering",
    "same_cancer_type_tumours_available",
    "n_same_cancer_type_tumours_retained",
    "top_n_cutoff_used",
    "component_id",
    "component_label",
    "component_dominant_lineage",
    "component_size",
    "component_tumour_count",
    "component_purity",
    "n_retained_tumours_in_component",
    "fraction_retained_tumours_in_component"
  ))
  setorder(comp, cell_line_cancer_type, cell_line_group, -fraction_retained_tumours_in_component, component_id)
  comp
}

build_all_top50_component_composition <- function(ranked, component_meta, top_n) {
  # Select the all-top-50 cell-line-to-tumour ranks without cancer-type filtering
  all_top <- ranked[order(cell_line_group, rank_cellline_to_tumour)][
    ,
    head(.SD, min(top_n, .N)),
    by = .(cell_line_group, cell_line_cancer_type)
  ]
  retained <- all_top[, .(
    n_top50_tumours_retained = .N,
    n_same_cancer_type_tumours_in_top50 = sum(same_cancer_type),
    fraction_top50_tumours_matching_cell_line_cancer_type = mean(same_cancer_type)
  ), by = .(cell_line_group, cell_line_cancer_type)]

  comp <- all_top[, .(
    n_top50_tumours_in_component = .N
  ), by = .(cell_line_group, cell_line_cancer_type, component_id)]
  comp <- merge(comp, retained, by = c("cell_line_group", "cell_line_cancer_type"))
  comp <- merge(comp, component_meta, by = "component_id", all.x = TRUE)
  comp[, top_n_cutoff_used := top_n]
  comp[, fraction_top50_tumours_in_component :=
         n_top50_tumours_in_component / n_top50_tumours_retained]
  setcolorder(comp, c(
    "cell_line_group",
    "cell_line_cancer_type",
    "top_n_cutoff_used",
    "component_id",
    "component_label",
    "component_dominant_lineage",
    "component_size",
    "component_tumour_count",
    "component_purity",
    "n_top50_tumours_retained",
    "n_top50_tumours_in_component",
    "fraction_top50_tumours_in_component",
    "n_same_cancer_type_tumours_in_top50",
    "fraction_top50_tumours_matching_cell_line_cancer_type"
  ))
  setorder(comp, cell_line_cancer_type, cell_line_group, -fraction_top50_tumours_in_component, component_id)
  comp
}

build_component_composition_summary <- function(same_comp, all_comp) {
  # Summarise component-composition denominators and largest component fractions
  same_summary <- same_comp[, .(
    n_same_cancer_type_tumours_retained =
      unique(n_same_cancer_type_tumours_retained)[1],
    n_components_same_cancer_type = .N,
    largest_same_cancer_type_component_fraction =
      max(fraction_retained_tumours_in_component)
  ), by = .(cell_line_group, cell_line_cancer_type)]

  all_summary <- all_comp[, .(
    n_components_all_top50 = .N,
    largest_all_top50_component_fraction =
      max(fraction_top50_tumours_in_component),
    fraction_all_top50_tumours_matching_cell_line_cancer_type =
      unique(fraction_top50_tumours_matching_cell_line_cancer_type)[1]
  ), by = .(cell_line_group, cell_line_cancer_type)]

  summary <- merge(same_summary, all_summary,
                   by = c("cell_line_group", "cell_line_cancer_type"),
                   all = TRUE)
  setorder(summary, cell_line_cancer_type, cell_line_group)
  summary
}

check_fraction_sums <- function(dt, fraction_col, label) {
  sums <- dt[, .(fraction_sum = sum(get(fraction_col))),
             by = .(cell_line_group, cell_line_cancer_type)]
  bad <- sums[abs(fraction_sum - 1) > 1e-8]
  if (nrow(bad)) {
    stop(label, " component fractions do not sum to 1 for: ",
         paste(head(bad$cell_line_group, 10L), collapse = ", "))
  }
  invisible(sums)
}

component_palette <- function(labels) {
  labels <- sort(unique(as.character(labels)))
  base <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Set2")
  )
  pal <- grDevices::colorRampPalette(base)(length(labels))
  setNames(pal, labels)
}

theme_component_plot <- function(config) {
  theme_bw(base_size = config$base_font_size) +
    theme(
      plot.title = element_text(size = config$title_font_size, face = "bold"),
      plot.subtitle = element_text(size = config$subtitle_font_size),
      axis.title = element_text(size = config$axis_title_font_size, face = "bold"),
      axis.text = element_text(size = config$axis_text_font_size),
      strip.text = element_text(size = config$facet_label_font_size, face = "bold"),
      legend.title = element_text(size = config$legend_title_font_size, face = "bold"),
      legend.text = element_text(size = config$legend_text_font_size),
      legend.key.size = unit(0.42, "cm"),
      plot.margin = margin(14, 18, 14, 18)
    )
}

prepare_plot_data <- function(dt, fraction_col, count_col) {
  plot_dt <- copy(dt)
  group_order <- unique(plot_dt[
    order(match(cell_line_cancer_type, lineage_levels), cell_line_group)
  ]$cell_line_group)
  plot_dt[, cell_line_group := factor(cell_line_group, levels = rev(group_order))]
  component_order <- unique(plot_dt[order(as.integer(component_id))]$component_label)
  plot_dt[, component_label := factor(component_label, levels = component_order)]
  plot_dt[, plot_fraction := get(fraction_col)]
  plot_dt[, plot_count := get(count_col)]
  plot_dt
}

plot_same_cancer_component_composition <- function(same_comp, output_pdf, output_png) {
  # Plot same-cancer-type retained tumour component fractions by cell-line group
  config <- same_cancer_component_plot_config
  plot_dt <- prepare_plot_data(
    same_comp,
    "fraction_retained_tumours_in_component",
    "n_retained_tumours_in_component"
  )
  pal <- component_palette(plot_dt$component_label)
  p <- ggplot(plot_dt, aes(x = plot_fraction, y = cell_line_group,
                           fill = component_label)) +
    geom_col(width = config$bar_width) +
    facet_grid(cell_line_cancer_type ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1.0001), expand = c(0, 0)) +
    scale_y_discrete(drop = TRUE) +
    scale_fill_manual(values = pal, name = "Patient consensus-network component") +
    guides(fill = guide_legend(ncol = config$legend_columns, byrow = TRUE)) +
    labs(
      title = "Same-cancer-type top-ranked tumours by component",
      subtitle = paste0("Top ", config$top_n,
                        " same-cancer-type tumours retained per cell-line group where available"),
      x = "Fraction of retained tumours",
      y = "Cell-line group"
    ) +
    theme_component_plot(config)
  ggsave(output_pdf, p, width = config$plot_width_in,
         height = config$plot_height_in, units = "in", device = cairo_pdf)
  ggsave(output_png, p, width = config$plot_width_in,
         height = config$plot_height_in, units = "in",
         dpi = config$png_dpi, bg = "white")
}

plot_all_top50_component_composition <- function(all_comp, output_pdf, output_png) {
  # Plot all top-50 retained tumour component fractions by cell-line group
  config <- all_top50_component_plot_config
  plot_dt <- prepare_plot_data(
    all_comp,
    "fraction_top50_tumours_in_component",
    "n_top50_tumours_in_component"
  )
  pal <- component_palette(plot_dt$component_label)
  p <- ggplot(plot_dt, aes(x = plot_fraction, y = cell_line_group,
                           fill = component_label)) +
    geom_col(width = config$bar_width) +
    facet_grid(cell_line_cancer_type ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1.0001), expand = c(0, 0)) +
    scale_y_discrete(drop = TRUE) +
    scale_fill_manual(values = pal, name = "Patient consensus-network component") +
    guides(fill = guide_legend(ncol = config$legend_columns, byrow = TRUE)) +
    labs(
      title = "All top-50 ranked tumours by component",
      subtitle = paste0("Top ", config$top_n,
                        " cell-line-to-tumour ranks retained per cell-line group"),
      x = "Fraction of top-50 tumours",
      y = "Cell-line group"
    ) +
    theme_component_plot(config)
  ggsave(output_pdf, p, width = config$plot_width_in,
         height = config$plot_height_in, units = "in", device = cairo_pdf)
  ggsave(output_png, p, width = config$plot_width_in,
         height = config$plot_height_in, units = "in",
         dpi = config$png_dpi, bg = "white")
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
mode <- required_arg(opts, "mode")
if (!mode %in% c("same", "all")) stop("--mode must be one of: same, all")

scores_path <- required_arg(opts, "scores-long")
components_path <- required_arg(opts, "components")
top_n <- if (!is.null(opts[["top-n"]])) as.integer(opts[["top-n"]]) else 50L
same_cancer_component_plot_config$top_n <- top_n
all_top50_component_plot_config$top_n <- top_n

scores <- normalise_scores(read_table_auto(scores_path))
components <- normalise_components(read_table_auto(components_path))
ranked <- rank_cellline_to_tumour_scores(scores, components)
component_meta <- build_component_metadata(ranked)

cat("Input score table:", scores_path, "\n")
cat("Input component table:", components_path, "\n")
cat("Cell-line-to-tumour ranking direction: ranks assigned within cell_line_group by descending score\n")
cat("Cell-line groups:", uniqueN(ranked$cell_line_group), "\n")
cat("Ranked tumours:", uniqueN(ranked$tumour), "\n")
cat("Tumour components represented:", uniqueN(component_meta$component_id), "\n")
cat("Top-N cutoff:", top_n, "\n")
cat("Cell-line groups per cancer type:\n")
print(ranked[, .N, by = .(cell_line_group, cell_line_cancer_type)][
  , .(n_cell_line_groups = .N), by = cell_line_cancer_type][
  order(match(cell_line_cancer_type, lineage_levels))
])

if (mode == "same") {
  same_table <- required_arg(opts, "same-table")
  same_pdf <- required_arg(opts, "same-pdf")
  same_png <- required_arg(opts, "same-png")
  dir.create(dirname(same_table), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(same_pdf), recursive = TRUE, showWarnings = FALSE)

  same_comp <- build_same_cancer_component_composition(ranked, component_meta, top_n)
  check_fraction_sums(same_comp, "fraction_retained_tumours_in_component",
                      "Same-cancer-type")
  denom_check <- unique(same_comp[, .(
    cell_line_group,
    same_cancer_type_tumours_available,
    n_same_cancer_type_tumours_retained,
    top_n_cutoff_used
  )])
  denom_check[, expected_retained :=
                pmin(top_n_cutoff_used, same_cancer_type_tumours_available)]
  if (any(denom_check$n_same_cancer_type_tumours_retained != denom_check$expected_retained)) {
    stop("Same-cancer-type retained tumour denominators are incorrect.")
  }
  fwrite(same_comp, same_table, sep = "\t")
  plot_same_cancer_component_composition(same_comp, same_pdf, same_png)

  cat("Same-cancer-type component table:", same_table, "\n")
  cat("Same-cancer-type plot PDF:", same_pdf, "\n")
  cat("Same-cancer-type plot PNG:", same_png, "\n")
  cat("Rows:", nrow(same_comp), "\n")
  cat("Components represented:", uniqueN(same_comp$component_id), "\n")
  cat("Fraction-sum check: PASS\n")
  cat("Denominator check: PASS\n")
}

if (mode == "all") {
  same_table <- required_arg(opts, "same-table")
  all_table <- required_arg(opts, "all-table")
  summary_table <- required_arg(opts, "summary-table")
  all_pdf <- required_arg(opts, "all-pdf")
  all_png <- required_arg(opts, "all-png")
  dir.create(dirname(all_table), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(all_pdf), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(same_table)) stop("Missing same-cancer-type table: ", same_table)
  same_comp <- fread(same_table)
  all_comp <- build_all_top50_component_composition(ranked, component_meta, top_n)
  check_fraction_sums(all_comp, "fraction_top50_tumours_in_component",
                      "All-top-50")
  denom_check <- unique(all_comp[, .(
    cell_line_group,
    n_top50_tumours_retained,
    top_n_cutoff_used
  )])
  denom_check[, expected_retained := pmin(top_n_cutoff_used, uniqueN(ranked$tumour))]
  if (any(denom_check$n_top50_tumours_retained != denom_check$expected_retained)) {
    stop("All-top-50 retained tumour denominators are incorrect.")
  }
  composition_summary <- build_component_composition_summary(same_comp, all_comp)

  fwrite(all_comp, all_table, sep = "\t")
  fwrite(composition_summary, summary_table, sep = "\t")
  plot_all_top50_component_composition(all_comp, all_pdf, all_png)

  cat("All-top-50 component table:", all_table, "\n")
  cat("Component-composition summary table:", summary_table, "\n")
  cat("All-top-50 plot PDF:", all_pdf, "\n")
  cat("All-top-50 plot PNG:", all_png, "\n")
  cat("Rows:", nrow(all_comp), "\n")
  cat("Components represented:", uniqueN(all_comp$component_id), "\n")
  cat("Fraction-sum check: PASS\n")
  cat("Denominator check: PASS\n")
}
