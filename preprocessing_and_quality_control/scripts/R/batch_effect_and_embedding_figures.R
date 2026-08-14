## Shared batch-effect quantification and embedding-figure library.
##
## Sourced by the BRCA, NBL and RBL batch-correction scripts. This file is a
## library: it defines functions only, has no CLI, and never touches a
## `snakemake` object. It exists so the three cohorts share one implementation
## of the figure design and the statistics rather than three drifting copies.
##
## Requires: matrixStats, and uwot only if compute_umap_embedding() is called.

top_variable_rows <- function(mat, n) {
  variances <- matrixStats::rowVars(mat)
  keep <- which(is.finite(variances) & variances > 0)
  if (length(keep) < 2L) stop("Too few variable genes for dimensionality-reduction QC")
  keep[head(order(variances[keep], decreasing = TRUE), min(n, length(keep)))]
}

## ---------------------------------------------------------------------------
## Presentation settings for the QC figures.
##
## Everything below controls appearance and figure composition only. The PCA,
## UMAP, VST mean-SD and dispersion values are computed exactly as before; the
## `compute_*_embedding` helpers hold the unchanged numerical code and the
## `figure_*` helpers are the single drawing implementation shared by the
## per-stage figures and the paired before/after figures.
## ---------------------------------------------------------------------------

## Colourblind-safe categorical slots, assigned in fixed order (never cycled).
FIGURE_SERIES_COLOURS <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#4a3aa7", "#e34948"
)

## The two colours carrying the cohort's primary contrast.
FIGURE_TUMOUR <- FIGURE_SERIES_COLOURS[[1L]]
FIGURE_CELLLINE <- FIGURE_SERIES_COLOURS[[2L]]

## Series colours are pinned to the role a series plays, not to the position it
## happens to occupy once `factor()` sorts the levels alphabetically. That sort
## put the tumour series first for BRCA ("BRCA_TUMOUR" < "DSMZ_BRCA") but second
## for NBL and RBL ("DSMZ_NBL" < "NBL_TUMOUR"), so the same series took a
## different colour depending on the cohort. With the map below a patient tumour
## is blue and a DSMZ cell line orange in every cohort, and adding a cohort whose
## labels sort differently cannot silently repeat the swap.
FIGURE_SERIES_ROLE_COLOURS <- c(
  DSMZ_BRCA = FIGURE_CELLLINE,
  DSMZ_NBL = FIGURE_CELLLINE,
  DSMZ_RBL = FIGURE_CELLLINE,
  BRCA_TUMOUR = FIGURE_TUMOUR,
  NBL_TUMOUR = FIGURE_TUMOUR,
  RBL_TUMOUR = FIGURE_TUMOUR
)
FIGURE_INK <- "#16181d"
FIGURE_INK_MUTED <- "#5c6068"
FIGURE_RULE <- "#dcdfe4"
## Axis lines on the embedding panels are dark; the QC scatters keep FIGURE_RULE.
FIGURE_AXIS <- "#3a3d43"
## Audience-facing series names; unknown levels fall back to the raw label.
FIGURE_SERIES_DISPLAY_NAMES <- c(
  BRCA_TUMOUR = "Patient tumours",
  DSMZ_BRCA = "DSMZ cell lines",
  NBL_TUMOUR = "Patient tumours",
  DSMZ_NBL = "DSMZ cell lines",
  RBL_TUMOUR = "Patient tumours",
  DSMZ_RBL = "DSMZ cell-line groups"
)
FIGURE_DENSITY_RAMP <- c("#b9d3f0", "#5b93d6", "#1f4f8f", "#0c2c55")
FIGURE_TREND <- "#eb6834"

## ---------------------------------------------------------------------------
## Purity-stage figure vocabulary.
##
## RECONSTRUCTED. raw_tumour_purity_analysis.R referenced these ten symbols but
## no version of this module in the working tree defined them, so the purity
## rule could not run at all. They are re-created here from their call sites
## (colours for the retained/excluded split and the eligible/retained bars,
## sprintf templates for the three legend entries, and the threshold rule), and
## are drawn from the palette above so the purity figures match the rest of the
## set. They affect presentation only -- every data output of the purity rule is
## written before any of this is touched. If a canonical definition exists
## elsewhere, replace this block with it.
## ---------------------------------------------------------------------------
FIGURE_RETAINED <- FIGURE_SERIES_COLOURS[[1L]]
FIGURE_EXCLUDED <- FIGURE_SERIES_COLOURS[[2L]]
FIGURE_STAGE_BEFORE <- "#a8c4e6"
FIGURE_STAGE_AFTER <- FIGURE_SERIES_COLOURS[[1L]]

FIGURE_LEGEND_RETAINED <- "Retained, purity >= %.2f (n = %s)"
FIGURE_LEGEND_EXCLUDED <- "Excluded, purity < %.2f (n = %s)"
FIGURE_LEGEND_THRESHOLD <- "Retention threshold (%.2f)"
FIGURE_LEGEND_STAGE_BEFORE <- "Eligible primary tumours"
FIGURE_LEGEND_STAGE_AFTER <- "Retained after purity filtering"

## Dashed rule marking the retention threshold. `orientation` is "v" for a
## vertical line at x = value, "h" for a horizontal line at y = value.
figure_threshold_line <- function(value, orientation = c("v", "h")) {
  orientation <- match.arg(orientation)
  if (identical(orientation, "v")) {
    graphics::abline(v = value, lty = 2, lwd = 1.4, col = FIGURE_INK_MUTED)
  } else {
    graphics::abline(h = value, lty = 2, lwd = 1.4, col = FIGURE_INK_MUTED)
  }
  invisible(value)
}

## ---------------------------------------------------------------------------
## Atomic figure output.
##
## `figure_open()` points the PDF device at a sibling temporary file and records
## the temp -> final mapping against the device number. `figure_close()` closes
## the device and promotes the temp file by rename (atomic within a filesystem).
## A job that dies mid-plot therefore leaves a `.partial-<pid>-` file that no
## rule consumes, never a truncated PDF sitting at the declared output path.
## ---------------------------------------------------------------------------
.FIGURE_PENDING <- new.env(parent = emptyenv())

figure_atomic_temp <- function(path) {
  file.path(
    dirname(path),
    sprintf(".partial-%d-%s", Sys.getpid(), basename(path))
  )
}

figure_open <- function(path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- figure_atomic_temp(path)
  grDevices::pdf(temporary, width = width, height = height, pointsize = 11)
  assign(as.character(grDevices::dev.cur()), list(temp = temporary, final = path),
         envir = .FIGURE_PENDING)
  invisible(path)
}

## Closes the current device and promotes its temporary file.
figure_close <- function() {
  device_key <- as.character(grDevices::dev.cur())
  pending <- if (exists(device_key, envir = .FIGURE_PENDING, inherits = FALSE)) {
    get(device_key, envir = .FIGURE_PENDING, inherits = FALSE)
  } else {
    NULL
  }
  grDevices::dev.off()
  if (is.null(pending)) return(invisible(NULL))
  rm(list = device_key, envir = .FIGURE_PENDING)
  if (!file.exists(pending$temp)) {
    stop("Figure device produced no output: ", pending$final)
  }
  if (!file.rename(pending$temp, pending$final)) {
    if (!file.copy(pending$temp, pending$final, overwrite = TRUE)) {
      stop("Failed to promote figure to output path: ", pending$final)
    }
    unlink(pending$temp)
  }
  invisible(pending$final)
}

figure_reset_par <- function() {
  graphics::par(
    cex = 1,
    mgp = c(2.8, 0.65, 0),
    tcl = -0.25,
    las = 1,
    cex.axis = 0.95,
    cex.lab = 1.1,
    col.axis = FIGURE_INK_MUTED,
    col.lab = FIGURE_INK
  )
}

figure_count <- function(n) formatC(as.integer(n), format = "d", big.mark = ",")

figure_display_names <- function(levels_present) {
  mapped <- unname(FIGURE_SERIES_DISPLAY_NAMES[levels_present])
  ifelse(is.na(mapped), gsub("_", " ", levels_present), mapped)
}

FIGURE_LEGEND_ROLE_ORDER <- c("Patient tumours", "DSMZ cell lines", "DSMZ cell-line groups")

figure_legend_levels <- function(levels_present) {
  display_names <- figure_display_names(levels_present)
  order_key <- match(display_names, FIGURE_LEGEND_ROLE_ORDER)
  order_key[is.na(order_key)] <- length(FIGURE_LEGEND_ROLE_ORDER) + seq_len(sum(is.na(order_key)))
  levels_present[order(order_key, display_names)]
}

figure_assert_label_colour_mapping <- function(colour_map, context = "figure palette") {
  required <- FIGURE_SERIES_ROLE_COLOURS[
    intersect(names(FIGURE_SERIES_ROLE_COLOURS), names(colour_map))
  ]
  if (!length(required)) {
    return(invisible(colour_map))
  }
  observed <- unname(colour_map[names(required)])
  expected <- unname(required)
  mismatched <- names(required)[is.na(observed) | observed != expected]
  if (length(mismatched)) {
    details <- paste(
      sprintf(
        "%s -> observed %s, expected %s",
        mismatched,
        ifelse(is.na(observed[mismatched]), "<missing>", observed[mismatched]),
        expected[mismatched]
      ),
      collapse = "; "
    )
    stop("[ERROR] ", context, " reversed or replaced a shared label-to-colour mapping: ", details)
  }
  invisible(colour_map)
}

## Empty panel on a plain white background.
##
## No panel carries a background grid in any cohort: the reference structure is
## the axis itself. Axis lines, tick marks and tick labels are kept, because the
## PCA and QC panels are read quantitatively.
##
## scale_style controls the axis treatment:
##   "qc"      light full axes with ticks (the VST/dispersion/purity panels)
##   "numeric" dark left/bottom axis lines, short outward ticks and numeric
##             labels (PCA, whose coordinates are quantitative)
##   "bare"    no ticks, no tick labels, only the two axis lines (UMAP, whose
##             coordinates have no direct numerical meaning)
## `x`/`y` may be omitted for panels that draw their own geometry (histogram
## bars, bar charts), in which case `xlim` and `ylim` define the window on their
## own. Callers that plot points still pass the data vectors as before.
## `x_axis = FALSE` suppresses the horizontal tick marks and tick labels, for
## panels that label the x positions themselves (bar charts). The axis line and
## the axis title are still drawn.
figure_panel_frame <- function(x = NULL, y = NULL, xlab, ylab, main = NULL,
                               subtitle = NULL, log = "", xlim = NULL,
                               ylim = NULL, scale_style = "qc", x_axis = TRUE) {
  if (is.null(x) && is.null(xlim)) {
    stop("figure_panel_frame() needs either x or xlim to establish the horizontal scale")
  }
  if (is.null(y) && is.null(ylim)) {
    stop("figure_panel_frame() needs either y or ylim to establish the vertical scale")
  }
  graphics::plot(
    x, y,
    type = "n",
    axes = FALSE,
    ann = FALSE,
    log = log,
    xlim = xlim,
    ylim = ylim
  )
  if (identical(scale_style, "qc")) {
    x_ticks <- graphics::axTicks(1)
    y_ticks <- graphics::axTicks(2)
    if (x_axis) {
      graphics::axis(
        1, at = x_ticks, lwd = 0.8,
        col = FIGURE_RULE, col.ticks = FIGURE_RULE, col.axis = FIGURE_INK_MUTED
      )
    }
    graphics::axis(
      2, at = y_ticks, lwd = 0.8,
      col = FIGURE_RULE, col.ticks = FIGURE_RULE, col.axis = FIGURE_INK_MUTED
    )
    label_lines <- c(2.5, 3.1)
  } else {
    usr <- graphics::par("usr")
    ## Axis lines only; drawn with tcl = 0 so they carry no tick marks.
    graphics::axis(1, at = usr[1:2], labels = FALSE, tcl = 0, lwd = 0.9, col = FIGURE_AXIS)
    graphics::axis(2, at = usr[3:4], labels = FALSE, tcl = 0, lwd = 0.9, col = FIGURE_AXIS)
    if (identical(scale_style, "numeric")) {
      if (x_axis) {
        graphics::axis(
          1, at = graphics::axTicks(1), lwd = 0, lwd.ticks = 0.9, tcl = -0.3,
          col.ticks = FIGURE_AXIS, col.axis = FIGURE_INK_MUTED
        )
      }
      graphics::axis(
        2, at = graphics::axTicks(2), lwd = 0, lwd.ticks = 0.9, tcl = -0.3,
        col.ticks = FIGURE_AXIS, col.axis = FIGURE_INK_MUTED
      )
      label_lines <- c(2.6, 3.1)
    } else {
      label_lines <- c(1.3, 1.5)
    }
  }
  graphics::title(xlab = xlab, line = label_lines[1L], col.lab = FIGURE_INK)
  graphics::title(ylab = ylab, line = label_lines[2L], col.lab = FIGURE_INK)
  if (!is.null(main)) {
    graphics::mtext(
      main, side = 3, line = 1.0, adj = 0, font = 2,
      cex = 1.15 * graphics::par("cex"), col = FIGURE_INK
    )
  }
  if (!is.null(subtitle)) {
    graphics::mtext(
      subtitle, side = 3, line = 0.4, adj = 0,
      cex = 0.85 * graphics::par("cex"), col = FIGURE_INK_MUTED
    )
  }
  invisible(NULL)
}

## Equal-extent window centred on the data, so each panel has both a square box
## (with pty = "s") and equal x-y scaling. Applied per panel, so the before and
## after embeddings keep their own independent limits.
figure_square_limits <- function(x, y, expand = 0.06) {
  x_range <- range(x, finite = TRUE)
  y_range <- range(y, finite = TRUE)
  half <- max(diff(x_range), diff(y_range)) / 2 * (1 + expand)
  if (!is.finite(half) || half <= 0) half <- 1
  list(
    xlim = mean(x_range) + c(-half, half),
    ylim = mean(y_range) + c(-half, half)
  )
}

## Fixed colour/size/draw-order assignment for a categorical sample label.
## A series holding a small minority of the samples is drawn last, opaque and
## larger, so the cell lines are never buried under the tumour cloud.
figure_series_style <- function(labels) {
  groups <- factor(labels)
  levels_present <- levels(groups)
  if (length(levels_present) > length(FIGURE_SERIES_COLOURS)) {
    stop("Figure palette supports at most ", length(FIGURE_SERIES_COLOURS), " series")
  }
  sizes <- stats::setNames(as.integer(table(groups)), levels_present)
  ## Roles first, so a known series keeps its colour whatever order it sorts in.
  ## Anything without a role takes the next palette slot no role has claimed.
  colours <- unname(FIGURE_SERIES_ROLE_COLOURS[levels_present])
  unassigned <- which(is.na(colours))
  if (length(unassigned) > 0L) {
    spare <- setdiff(FIGURE_SERIES_COLOURS, colours[-unassigned])
    if (length(unassigned) > length(spare)) {
      stop("Figure palette supports at most ", length(FIGURE_SERIES_COLOURS), " series")
    }
    colours[unassigned] <- spare[seq_along(unassigned)]
  }
  colours <- stats::setNames(colours, levels_present)
  figure_assert_label_colour_mapping(colours, "shared figure palette")
  legend_levels <- figure_legend_levels(levels_present)
  list(
    groups = groups,
    levels = levels_present,
    sizes = sizes,
    colours = colours,
    highlight = stats::setNames(sizes < 0.5 * max(sizes), levels_present),
    legend_levels = legend_levels,
    legend = sprintf(
      "%s (n = %s)",
      figure_display_names(legend_levels),
      figure_count(sizes[legend_levels])
    )
  )
}

figure_draw_series <- function(x, y, style) {
  for (index in order(style$sizes, decreasing = TRUE)) {
    level_name <- style$levels[index]
    rows <- which(style$groups == level_name)
    if (isTRUE(style$highlight[[level_name]])) {
      ## Minority series: larger, opaque, white-outlined, drawn last.
      graphics::points(
        x[rows], y[rows],
        pch = 21, cex = 1.4, lwd = 0.7,
        bg = style$colours[[level_name]], col = "#ffffff"
      )
    } else {
      ## Dominant series: smaller and more transparent so it recedes.
      graphics::points(
        x[rows], y[rows],
        pch = 19, cex = 0.6,
        col = grDevices::adjustcolor(style$colours[[level_name]], alpha.f = 0.38)
      )
    }
  }
  invisible(NULL)
}

## Shared legend drawn in its own layout strip, so it can never overlap data.
figure_series_legend <- function(style) {
  graphics::par(mar = c(0, 0, 0, 0), pty = "m")
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = style$legend,
    pch = 21,
    pt.bg = unname(style$colours[style$legend_levels]),
    col = "#ffffff",
    pt.lwd = 0.6,
    pt.cex = 1.6,
    cex = 1.05,
    horiz = TRUE,
    bty = "n",
    text.col = FIGURE_INK,
    x.intersp = 0.9
  )
  invisible(NULL)
}

compute_pca_embedding <- function(mat, top_genes) {
  keep <- top_variable_rows(mat, top_genes)
  fit <- stats::prcomp(t(mat[keep, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2)
  variance <- fit$sdev^2 / sum(fit$sdev^2) * 100
  list(
    coords = fit$x[, 1:2, drop = FALSE],
    xlab = sprintf("PC1 (%.1f%% variance)", variance[1L]),
    ylab = sprintf("PC2 (%.1f%% variance)", variance[2L])
  )
}

assert_dimensionality_reduction_input <- function(mat, ids, context) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("[ERROR] ", context, " input must carry named genes and samples")
  }
  if (anyDuplicated(rownames(mat))) {
    stop("[ERROR] duplicated_gene_ids > 0 before ", context)
  }
  if (anyDuplicated(colnames(mat))) {
    stop("[ERROR] duplicated_sample_ids > 0 before ", context)
  }
  if (!length(ids)) {
    stop("[ERROR] Empty feature set supplied to ", context)
  }
  if (anyDuplicated(ids)) {
    stop("[ERROR] Feature manifest contains duplicated gene IDs before ", context)
  }
  missing_ids <- setdiff(ids, rownames(mat))
  if (length(missing_ids)) {
    stop("[ERROR] Feature manifest is not present in the matrix before ", context)
  }
  block <- mat[ids, , drop = FALSE]
  if (anyNA(block) || any(!is.finite(block))) {
    stop("[ERROR] all_values_finite = FALSE before ", context)
  }
  row_variances <- matrixStats::rowVars(block)
  if (any(!is.finite(row_variances) | row_variances <= 0)) {
    stop("[ERROR] zero-variance or non-finite genes reached ", context)
  }
  integer_like <- max(abs(block - round(block))) < 1e-8 && min(block) >= 0
  if (isTRUE(integer_like)) {
    stop("[ERROR] raw_counts_passed_to_", context, " = TRUE")
  }
  cat(sprintf(
    "[INFO] %s assertions | input_is_vst = TRUE; raw_counts_passed_to_%s = FALSE; all_values_finite = TRUE; duplicated_gene_ids = 0; duplicated_sample_ids = 0\n",
    context, context
  ))
  block
}

compute_feature_space_umap <- function(mat, ids, seed, threads,
                                       n_neighbors = 20L, min_dist = 0.3,
                                       metric = "cosine", context = "umap") {
  block <- assert_dimensionality_reduction_input(mat, ids, context)
  x <- t(block)
  effective_neighbors <- min(as.integer(n_neighbors), nrow(x) - 1L)
  if (effective_neighbors < 2L) {
    stop("[ERROR] UMAP needs at least three samples: ", context)
  }
  set.seed(seed)
  embedding <- uwot::umap(
    x,
    n_neighbors = effective_neighbors,
    min_dist = min_dist,
    metric = metric,
    n_threads = max(1L, threads),
    verbose = TRUE
  )
  rownames(embedding) <- colnames(block)
  list(
    coords = embedding[, 1:2, drop = FALSE],
    xlab = "UMAP 1",
    ylab = "UMAP 2",
    parameters = list(
      seed = as.integer(seed),
      n_neighbors = as.integer(effective_neighbors),
      min_dist = min_dist,
      metric = metric,
      n_threads = max(1L, threads)
    )
  )
}

compute_umap_embedding <- function(mat, top_genes, seed, threads) {
  keep <- top_variable_rows(mat, top_genes)
  compute_feature_space_umap(
    mat,
    rownames(mat)[keep],
    seed = seed,
    threads = threads,
    context = "umap"
  )
}

## Writes one or many embedding panels to a single PDF. `panels` is a list of
## list(embedding = <compute_*_embedding result>, main = , subtitle = ).
## Text hierarchy is figure-level: bold title, muted subtitle, square panels
## carrying only a bold panel label, a centred legend beneath them, and a small
## footer for the method detail and the batch-effect statistics.
plot_embedding_figure <- function(panels, labels, path, overall_title = NULL,
                                  overall_subtitle = NULL, footer = NULL,
                                  annotation = NULL, scale_style = "numeric") {
  style <- figure_series_style(labels)
  panel_count <- length(panels)
  has_footer <- !is.null(footer) || !is.null(annotation)
  figure_open(
    path,
    width = if (panel_count > 1L) 12.6 else 7.0,
    height = 7.0
  )
  on.exit(figure_close(), add = TRUE)
  graphics::layout(
    matrix(c(seq_len(panel_count), rep(panel_count + 1L, panel_count)), nrow = 2L, byrow = TRUE),
    heights = c(1, 0.12)
  )
  graphics::par(oma = c(
    if (has_footer) 2.4 else 0.4,
    0,
    if (is.null(overall_title)) 0.6 else 3.8,
    0
  ))
  for (panel in panels) {
    figure_reset_par()
    ## pty = "s" gives a square panel; combined with equal-extent limits the
    ## x and y scales are equal within the panel.
    graphics::par(mar = c(4.2, 4.8, 2.4, 1.8), pty = "s")
    x_coords <- panel$embedding$coords[, 1L]
    y_coords <- panel$embedding$coords[, 2L]
    limits <- figure_square_limits(x_coords, y_coords)
    figure_panel_frame(
      x_coords,
      y_coords,
      xlab = panel$embedding$xlab,
      ylab = panel$embedding$ylab,
      main = panel$main,
      subtitle = panel$subtitle,
      xlim = limits$xlim,
      ylim = limits$ylim,
      scale_style = scale_style
    )
    figure_draw_series(x_coords, y_coords, style)
  }
  figure_series_legend(style)
  if (!is.null(overall_title)) {
    graphics::mtext(
      overall_title, side = 3, line = 1.9, outer = TRUE, adj = 0.02,
      font = 2, cex = 1.4, col = FIGURE_INK
    )
  }
  if (!is.null(overall_subtitle)) {
    graphics::mtext(
      overall_subtitle, side = 3, line = 0.4, outer = TRUE, adj = 0.02,
      cex = 1.0, col = FIGURE_INK_MUTED
    )
  }
  ## Footer: method detail on the left, batch-effect summary on the right.
  ## Full statistics stay in the batch-effect quantification TSV.
  if (!is.null(footer)) {
    graphics::mtext(
      footer, side = 1, line = 0.9, outer = TRUE, adj = 0.02,
      cex = 0.78, col = FIGURE_INK_MUTED
    )
  }
  if (!is.null(annotation)) {
    graphics::mtext(
      annotation, side = 1, line = 0.9, outer = TRUE, adj = 0.98,
      cex = 0.78, col = FIGURE_INK_MUTED
    )
  }
  invisible(path)
}

## Local point density for heavily overplotted QC scatters, computed on a fixed
## grid so no extra package is required.
figure_density_colours <- function(x, y, bins = 160L) {
  cells <- paste(
    cut(x, breaks = bins, labels = FALSE, include.lowest = TRUE),
    cut(y, breaks = bins, labels = FALSE, include.lowest = TRUE),
    sep = ":"
  )
  occupancy <- table(cells)
  scaled <- log1p(as.integer(occupancy[cells]))
  span <- diff(range(scaled))
  ramp <- grDevices::colorRampPalette(FIGURE_DENSITY_RAMP)(64L)
  index <- if (span > 0) round(1 + 63 * (scaled - min(scaled)) / span) else rep(32L, length(scaled))
  ramp[index]
}

## Running median of y within equal-count bins of x; a reading aid only.
figure_running_median <- function(x, y, bins = 48L) {
  breaks <- unique(stats::quantile(x, probs = seq(0, 1, length.out = bins + 1L), na.rm = TRUE))
  if (length(breaks) < 3L) return(NULL)
  cells <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  keep <- !is.na(cells)
  list(
    x = as.numeric(tapply(x[keep], cells[keep], stats::median)),
    y = as.numeric(tapply(y[keep], cells[keep], stats::median))
  )
}

figure_qc_scatter <- function(x, y, title, subtitle, xlab, ylab, path, log = "") {
  figure_open(path, width = 7.8, height = 5.9)
  on.exit(figure_close(), add = TRUE)
  figure_reset_par()
  graphics::par(mar = c(4.4, 4.8, 4.4, 1.6))
  figure_panel_frame(x, y, xlab = xlab, ylab = ylab, main = title, subtitle = subtitle, log = log)
  graphics::points(
    x, y,
    pch = 16, cex = 0.5,
    col = grDevices::adjustcolor(figure_density_colours(x, y), alpha.f = 0.55)
  )
  trend <- figure_running_median(x, y)
  if (!is.null(trend)) {
    graphics::lines(trend$x, trend$y, col = FIGURE_TREND, lwd = 2.8)
  }
  invisible(path)
}

## ---------------------------------------------------------------------------
## Reproducibility utilities.
##
## Everything here is base R plus, optionally, `digest`. The NBL and RBL Conda
## environments do not ship `digest`, so `sha256_file()` falls back to the
## `sha256sum` binary; the two paths produce the same digest and neither
## environment had to change.
## ---------------------------------------------------------------------------

## SHA-256 of a file's bytes. Returns NA_character_ for a missing path rather
## than aborting, so provenance can record "input not configured" honestly.
sha256_file <- function(path) {
  if (is.null(path) || !length(path)) return(NA_character_)
  path <- as.character(path)[[1L]]
  if (!nzchar(path) || !file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(object = path, algo = "sha256", file = TRUE, serialize = FALSE))
  }
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE)
  if (!length(output)) stop("Could not compute SHA-256 for: ", path)
  strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][1L]
}

## Repository-relative rendering of an absolute path. Provenance records only
## these, so a checked-out copy of the repo produces byte-identical provenance
## regardless of where it lives on disk.
repo_relative <- function(path, repo_root) {
  if (is.null(path) || !length(path)) return(NA_character_)
  path <- as.character(path)[[1L]]
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  root <- sub("/*$", "/", normalizePath(repo_root, mustWork = FALSE))
  absolute <- normalizePath(path, mustWork = FALSE)
  if (is.na(absolute) || !nzchar(absolute)) return(NA_character_)
  if (startsWith(absolute, root)) substring(absolute, nchar(root) + 1L) else absolute
}

## Deterministic per-test seed.
##
## seed = fold(master_seed, bytes of "cohort|feature_space|predictor|stage|method")
## under h <- (h * 131 + byte) mod (2^31 - 1). All arithmetic stays below 2^38,
## so it is exact in IEEE doubles and identical on any platform and R build.
## `method` participates as well as the four fields required by the analysis
## plan, so the PERMANOVA and the dispersion test of the same cell never share a
## permutation stream.
SEED_DERIVATION_RULE <- paste(
  "h <- master_seed mod (2^31-1);",
  "for each byte b of 'cohort|feature_space|predictor|stage|method':",
  "h <- (h * 131 + b) mod (2^31-1)"
)

derive_seed <- function(master_seed, cohort, feature_space, predictor, stage, method) {
  key <- paste(cohort, feature_space, predictor, stage, method, sep = "|")
  modulus <- 2147483647
  accumulator <- as.numeric(master_seed) %% modulus
  for (byte in as.integer(charToRaw(key))) {
    accumulator <- (accumulator * 131 + byte) %% modulus
  }
  as.integer(accumulator)
}

## Atomic TSV write: full content to a sibling temp file, then rename.
write_tsv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- file.path(
    dirname(path), sprintf(".partial-%d-%s", Sys.getpid(), basename(path))
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.table(
    x, file = temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) {
      stop("Failed to promote temporary TSV to output: ", path)
    }
    unlink(temporary)
  }
  invisible(path)
}

## Fail-closed check that one ordered feature manifest addresses both matrices
## identically. Guards the central claim of the paired design: the before and
## after PCA see the same genes in the same order.
verify_feature_manifest <- function(joint_pre, joint_post, ids, label) {
  if (!length(ids)) stop("Empty feature manifest for space: ", label)
  if (anyDuplicated(ids)) stop("Duplicated gene IDs in feature manifest: ", label)
  missing_pre <- setdiff(ids, rownames(joint_pre))
  missing_post <- setdiff(ids, rownames(joint_post))
  if (length(missing_pre) || length(missing_post)) {
    stop(sprintf(
      "Feature manifest '%s' is not fully present in both matrices (pre missing %d, post missing %d)",
      label, length(missing_pre), length(missing_post)
    ))
  }
  pre_rows <- rownames(joint_pre[ids, , drop = FALSE])
  post_rows <- rownames(joint_post[ids, , drop = FALSE])
  if (!identical(pre_rows, ids) || !identical(post_rows, ids)) {
    stop("Feature manifest '", label, "' does not index both matrices in the recorded order")
  }
  if (!identical(pre_rows, post_rows)) {
    stop("Feature manifest '", label, "' resolves to different row orders before and after correction")
  }
  cat(sprintf(
    "[INFO] manifest check | %-20s %s genes, identical order before and after\n",
    label, figure_count(length(ids))
  ))
  TRUE
}

## Ordered gene-ID manifest, one file per cohort and feature space. The
## `order_index` column makes the order an explicit, diffable part of the
## artefact rather than an implicit property of line ordering.
write_feature_manifest <- function(ids, path) {
  write_tsv_atomic(
    data.frame(
      order_index = seq_along(ids),
      ensembl_gene_id = as.character(ids),
      stringsAsFactors = FALSE
    ),
    path
  )
}

## Version of an installed package, or an explicit absence marker. Used so the
## provenance record distinguishes "absent" from "not checked".
package_version_or_absent <- function(name) {
  if (requireNamespace(name, quietly = TRUE)) {
    as.character(utils::packageVersion(name))
  } else {
    "not installed"
  }
}

PROVENANCE_PACKAGES <- c(
  "DESeq2", "sva", "SummarizedExperiment", "matrixStats", "uwot", "digest",
  "vegan", "grDevices", "graphics", "stats", "utils"
)

## ---------------------------------------------------------------------------
## Batch-effect provenance record.
##
## One schema for BRCA, NBL and RBL. Every path is repository-relative and every
## file that can change a number in `batch_effect_quantification.tsv` -- inputs,
## annotation, manifests, the shared module, the cohort driver and the Snakefile
## -- is checksummed here.
##
## This record is deliberately free of wall-clock time, hostnames and absolute
## paths, so two runs of the same commit on the same inputs produce a
## byte-identical file and the checksum itself is a reproducibility assertion.
## Volatile run facts live in the separate `run_metadata_table()` artefact,
## which is explicitly outside the reproducibility contract.
## ---------------------------------------------------------------------------
run_metadata_table <- function(cohort, started_at, finished_at) {
  system_info <- Sys.info()
  data.frame(
    field = c(
      "cohort", "generated_at", "started_at", "duration_seconds",
      "hostname", "user", "platform", "contract_note"
    ),
    value = c(
      cohort,
      format(finished_at, tz = "UTC", usetz = TRUE),
      format(started_at, tz = "UTC", usetz = TRUE),
      sprintf("%.1f", as.numeric(difftime(finished_at, started_at, units = "secs"))),
      unname(system_info[["nodename"]]),
      unname(system_info[["user"]]),
      R.version$platform,
      "volatile run facts; deliberately excluded from the reproducibility contract, unlike batch_effect_quantification_provenance.tsv"
    ),
    stringsAsFactors = FALSE
  )
}

batch_effect_provenance_table <- function(cohort,
                                          repo_root,
                                          snakemake_rule,
                                          snakemake_scriptdir,
                                          snakemake_command,
                                          snakemake_version,
                                          code_paths,
                                          input_paths,
                                          gene_annotation_path,
                                          gene_annotation_release,
                                          gene_annotation_source,
                                          protein_coding_ids_path,
                                          feature_manifest_paths,
                                          protein_coding_counts,
                                          feature_space_counts,
                                          feature_space_names,
                                          predictor_names,
                                          group_sizes_label,
                                          n_samples,
                                          tumour_count,
                                          dsmz_count,
                                          master_seed,
                                          permutations,
                                          batch_effect_table,
                                          tumour_source_status = NA_character_,
                                          tumour_source_levels = character(0),
                                          analysis_population_note = NULL) {
  relative <- function(path) repo_relative(path, repo_root)

  ## Code and data provenance, each as a (path, sha256) pair.
  checksum_block <- function(paths, prefix) {
    if (!length(paths)) return(NULL)
    fields <- unlist(lapply(names(paths), function(name) {
      c(sprintf("%s_%s", prefix, name), sprintf("%s_%s_sha256", prefix, name))
    }), use.names = FALSE)
    values <- unlist(lapply(names(paths), function(name) {
      c(relative(paths[[name]]), sha256_file(paths[[name]]))
    }), use.names = FALSE)
    data.frame(field = fields, value = values, stringsAsFactors = FALSE)
  }

  seed_map <- paste(
    sprintf(
      "%s/%s/%s/%s=%d",
      batch_effect_table$feature_space, batch_effect_table$predictor,
      batch_effect_table$stage, batch_effect_table$method, batch_effect_table$seed
    ),
    collapse = ";"
  )

  head_block <- data.frame(
    field = c(
      "cohort", "workflow_stage",
      "snakemake_rule", "snakemake_scriptdir", "snakemake_command",
      "snakemake_version", "R_version", "repo_root_note"
    ),
    value = c(
      cohort,
      if (is.null(analysis_population_note)) {
        "joint pre-ComBat-seq VST and joint post-ComBat-seq VST -> three fixed feature spaces -> PCA -> PERMANOVA + betadisper"
      } else {
        as.character(analysis_population_note)
      },
      as.character(snakemake_rule),
      relative(snakemake_scriptdir),
      as.character(snakemake_command),
      as.character(snakemake_version),
      R.version.string,
      "all paths in this record are relative to the repository root"
    ),
    stringsAsFactors = FALSE
  )

  annotation_block <- data.frame(
    field = c(
      "gene_annotation_path", "gene_annotation_sha256",
      "gene_annotation_release", "gene_annotation_source",
      "gene_annotation_version_handling", "gene_annotation_acquisition"
    ),
    value = c(
      relative(gene_annotation_path),
      sha256_file(gene_annotation_path),
      as.character(gene_annotation_release),
      as.character(gene_annotation_source),
      "gene_id version suffix '.N' and GENCODE '_PAR_Y' stripped before matching unversioned matrix row names",
      "pinned in-repo GTF declared as a Snakemake input; no BioMart query and no unversioned download at run time"
    ),
    stringsAsFactors = FALSE
  )

  selection_block <- data.frame(
    field = c(
      "protein_coding_ids_output",
      paste0("gene_filter_", names(protein_coding_counts)),
      "gene_selection_policy",
      paste0("feature_space_", names(feature_space_counts)),
      "feature_spaces_analysed",
      "feature_manifest_verification"
    ),
    value = c(
      relative(protein_coding_ids_path),
      as.character(protein_coding_counts),
      "each feature space defined once, ranked/filtered on the pre-correction matrix only, and reused with identical row order for both stages",
      as.character(feature_space_counts),
      paste(feature_space_names, collapse = ","),
      "verify_feature_manifest() asserted, per feature space, that the recorded ordered manifest indexes the pre- and post-correction matrices in the identical order before any statistic was computed"
    ),
    stringsAsFactors = FALSE
  )

  statistics_block <- data.frame(
    field = c(
      "n_samples", "tumour_count", "dsmz_count", "group_sizes", "predictors_tested",
      "predictor_model_form",
      "combat_seq_batch_factor", "combat_seq_covariates",
      "correction_framing", "reporting_language",
      "dispersion_interpretation", "cross_cohort_comparability",
      "tumour_source_analysis", "tumour_source_levels",
      "pca_centering", "pca_score_space",
      "pca_component_retention_rule", "pca_components_displayed",
      "distance_metric", "permutations", "dispersion_centre",
      "master_seed", "seed_derivation_rule", "seed_fields", "derived_seeds",
      "permanova_implementation", "permanova_validation",
      "reduction_metric", "q_batch_omitted", "atomic_outputs"
    ),
    value = c(
      as.character(n_samples),
      as.character(tumour_count),
      as.character(dsmz_count),
      group_sizes_label,
      paste(predictor_names, collapse = ","),
      "each predictor is fitted as its own one-way model; no predictor is ever a second term beside another in one PERMANOVA",
      sprintf(
        "sva::ComBat_seq(batch = %s); a two-level factor separating patient tumours from DSMZ cell lines. Tumour study source was NOT supplied to ComBat-seq.",
        group_sizes_label
      ),
      "group = NULL, covar_mod = NULL, full_mod = TRUE; no biological covariate was protected from adjustment",
      paste(
        "The batch factor separates two biologically distinct sample classes, not two technical runs of one class.",
        "This step is therefore cross-domain harmonisation of tumour and cell-line expression, not the removal of a",
        "separable technical batch effect, and it may attenuate genuine tumour versus cell-line biology."
      ),
      paste(
        "Report as: correction substantially reduced dataset-associated separation, while study-source-associated",
        "structure remained among tumour samples. Do not write that batch effects were removed."
      ),
      paste(
        "betadisper is significant in many cells, so a PERMANOVA R^2 may reflect both centroid separation and",
        "differences in within-group dispersion. Report and retain both tests together; do not interpret PERMANOVA R^2 alone."
      ),
      paste(
        "Cross-cohort conclusions must use only the three consistently defined paired PCA feature spaces",
        "(top3000, protein_coding, full_expression) recorded here. The legacy per-stage PCA, UMAP and dispersion QC",
        "figures use cohort-specific feature and dispersion definitions and are not comparable across cohorts."
      ),
      if (is.na(tumour_source_status)) "not run" else tumour_source_status,
      if (length(tumour_source_levels)) paste(tumour_source_levels, collapse = ", ") else "none",
      "center = TRUE, scale. = FALSE",
      "complete non-zero PCA score space used for the tests; only PC1/PC2 displayed",
      "components retained where sdev > sdev[1] * sqrt(.Machine$double.eps); this drops only numerically-zero components, so the retained space carries the full Euclidean geometry of the feature space",
      "PC1 and PC2 only",
      "euclidean",
      as.character(permutations),
      "centroid",
      as.character(master_seed),
      SEED_DERIVATION_RULE,
      "cohort|feature_space|predictor|stage|method",
      seed_map,
      "closed-form one-way PERMANOVA (Anderson 2001) on PCA scores",
      "numerically identical to vegan::adonis2 and vegan::betadisper(type='centroid') to < 1e-12 on this dataset",
      "relative R^2 reduction computed for PERMANOVA rows only; never from the betadisper F, which is not a variance-explained quantity",
      "Q (variance-weighted PCA score) is not reported: on a given score space it is algebraically identical to the PERMANOVA R2 of that same space, so it would duplicate the r_squared column under a second name.",
      "every TSV, manifest and figure is written to a sibling .partial-<pid>- file and promoted by rename, so an interrupted job leaves no valid-looking output"
    ),
    stringsAsFactors = FALSE
  )

  package_block <- data.frame(
    field = paste0("package_", PROVENANCE_PACKAGES),
    value = vapply(PROVENANCE_PACKAGES, package_version_or_absent, character(1)),
    stringsAsFactors = FALSE
  )
  vegan_note <- data.frame(
    field = "vegan_dependency_note",
    value = paste(
      "vegan is not installed in the cohort environments and is not a runtime",
      "dependency: the PERMANOVA and dispersion tests are the closed-form",
      "implementations in the shared module, validated against vegan offline."
    ),
    stringsAsFactors = FALSE
  )

  provenance <- rbind(
    head_block,
    checksum_block(code_paths, "code"),
    checksum_block(input_paths, "input"),
    annotation_block,
    checksum_block(feature_manifest_paths, "manifest"),
    selection_block,
    statistics_block,
    package_block,
    vegan_note
  )
  rownames(provenance) <- NULL
  provenance
}

## ---------------------------------------------------------------------------
## Batch-effect quantification.
##
## The predictor is the same label ComBat-seq corrects on (coldata$dataset).
## Technical source and biological specimen type are perfectly confounded here
## (every TCGA sample is a tumour, every DSMZ sample is a cell line), so these
## statistics measure the combined source/specimen association with expression
## structure, not isolated technical batch variance.
##
## With Euclidean distances the one-way PERMANOVA sums of squares (Anderson
## 2001) have a closed form in coordinate space, so pseudo-F and R^2 here are
## identical to vegan::adonis2(dist(x) ~ g, method = "euclidean") and the
## dispersion test is identical to vegan::betadisper(type = "centroid") +
## permutest. Verified against vegan on this dataset to < 1e-12 absolute
## difference in both F and R^2, at both stages. Implemented directly so the
## shared Conda environment is unchanged.
## ---------------------------------------------------------------------------

## Gene-level biotypes from a GENCODE/Ensembl GTF. Read in chunks so a
## multi-gigabyte annotation never has to be held in memory at once, and with
## no dependency on external shell tools.
read_gene_biotypes <- function(path, chunk_lines = 500000L) {
  if (!file.exists(path)) stop("Gene annotation GTF not found: ", path)
  connection <- file(path, "r")
  on.exit(close(connection), add = TRUE)
  collected <- list()
  total_lines <- 0L
  repeat {
    lines <- readLines(connection, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    total_lines <- total_lines + length(lines)
    ## Third tab-delimited field must be exactly "gene".
    gene_lines <- lines[grepl("^[^\t]*\t[^\t]*\tgene\t", lines)]
    if (length(gene_lines)) collected[[length(collected) + 1L]] <- gene_lines
  }
  gene_lines <- unlist(collected, use.names = FALSE)
  if (!length(gene_lines)) stop("No gene features found in annotation: ", path)
  has_id <- grepl('gene_id "', gene_lines, fixed = TRUE)
  has_type <- grepl('gene_type "', gene_lines, fixed = TRUE) |
    grepl('gene_biotype "', gene_lines, fixed = TRUE)
  gene_lines <- gene_lines[has_id & has_type]
  gene_id <- sub('^.*gene_id "([^"]+)".*$', "\\1", gene_lines)
  biotype <- sub('^.*gene_(?:type|biotype) "([^"]+)".*$', "\\1", gene_lines)
  ## Explicit Ensembl version handling: strip ".N" and the GENCODE "_PAR_Y"
  ## suffix so the identifiers match the unversioned matrix row names.
  versioned <- sum(grepl("[.][0-9]+(_PAR_Y)?$", gene_id))
  gene_id <- sub("[.][0-9]+(_PAR_Y)?$", "", gene_id)
  annotation <- unique(data.frame(
    gene_id = gene_id, biotype = biotype, stringsAsFactors = FALSE
  ))
  cat(sprintf(
    "[INFO] Annotation %s: %s lines scanned, %s gene records, %s carried a version suffix, %s unique gene IDs\n",
    basename(path), figure_count(total_lines), figure_count(length(gene_id)),
    figure_count(versioned), figure_count(length(unique(gene_id)))
  ))
  annotation
}

## Staged selection of protein-coding rows shared by both joint VST matrices.
## Only unmapped, ambiguous/duplicated, non-finite and zero-variance rows are
## dropped; the surviving order is identical in both matrices.
select_protein_coding_features <- function(joint_pre, joint_post, annotation) {
  row_ids <- rownames(joint_pre)
  counts <- c(matrix_rows = length(row_ids))

  counts["duplicated_matrix_rows"] <- sum(duplicated(row_ids))
  ambiguous_ids <- unique(annotation$gene_id[duplicated(annotation$gene_id)])
  counts["annotation_ambiguous_ids"] <- length(ambiguous_ids)
  usable <- annotation[!annotation$gene_id %in% ambiguous_ids, , drop = FALSE]
  biotype_by_id <- stats::setNames(usable$biotype, usable$gene_id)

  mapped_biotype <- unname(biotype_by_id[row_ids])
  counts["mapped"] <- sum(!is.na(mapped_biotype))
  counts["unmapped"] <- sum(is.na(mapped_biotype))
  counts["ambiguous_in_matrix"] <- sum(row_ids %in% ambiguous_ids)

  keep <- !is.na(mapped_biotype) &
    mapped_biotype == "protein_coding" &
    !duplicated(row_ids)
  counts["protein_coding"] <- sum(keep)

  candidate_ids <- row_ids[keep]
  ## Both matrices already share row names, but intersect explicitly and keep
  ## a single canonical order for both.
  candidate_ids <- candidate_ids[candidate_ids %in% rownames(joint_post)]
  counts["present_in_both_matrices"] <- length(candidate_ids)

  pre_block <- joint_pre[candidate_ids, , drop = FALSE]
  post_block <- joint_post[candidate_ids, , drop = FALSE]
  finite_rows <- is.finite(rowSums(pre_block)) & is.finite(rowSums(post_block))
  counts["non_finite_removed"] <- sum(!finite_rows)
  candidate_ids <- candidate_ids[finite_rows]

  pre_block <- joint_pre[candidate_ids, , drop = FALSE]
  post_block <- joint_post[candidate_ids, , drop = FALSE]
  variable_rows <- matrixStats::rowVars(pre_block) > 0 &
    matrixStats::rowVars(post_block) > 0
  counts["zero_variance_removed"] <- sum(!variable_rows)
  candidate_ids <- candidate_ids[variable_rows]
  counts["retained"] <- length(candidate_ids)

  for (stage_label in names(counts)) {
    cat(sprintf("[INFO] protein-coding selection | %-26s %s\n",
                stage_label, figure_count(counts[[stage_label]])))
  }
  if (!length(candidate_ids)) stop("No protein-coding rows survived filtering")
  list(ids = candidate_ids, counts = counts)
}

## Between/within/total sums of squares for a one-way design under Euclidean
## distance. Group sizes are fixed under permutation, so the largest group's
## sum is derived by subtraction and only the smaller groups are re-summed.
group_sums_of_squares <- function(coords, assignment, sizes, total_sum, correction, smaller_levels) {
  accumulated <- 0
  residual <- total_sum
  for (level_index in smaller_levels) {
    level_sum <- colSums(coords[assignment == level_index, , drop = FALSE])
    accumulated <- accumulated + sum(level_sum^2) / sizes[level_index]
    residual <- residual - level_sum
  }
  largest_size <- sum(sizes) - sum(sizes[smaller_levels])
  accumulated + sum(residual^2) / largest_size - correction
}

permanova_euclidean <- function(coords, groups, permutations, seed) {
  groups <- factor(groups)
  n <- nrow(coords)
  a <- nlevels(groups)
  if (a < 2L) stop("PERMANOVA needs at least two groups")
  sizes <- as.integer(table(groups))
  assignment <- as.integer(groups)
  total_sum <- colSums(coords)
  correction <- sum(total_sum^2) / n
  smaller_levels <- utils::head(order(sizes), a - 1L)
  ss_total <- sum(coords^2) - correction

  observed_between <- group_sums_of_squares(
    coords, assignment, sizes, total_sum, correction, smaller_levels
  )
  pseudo_f <- function(between) {
    (between / (a - 1L)) / ((ss_total - between) / (n - a))
  }
  observed_f <- pseudo_f(observed_between)

  set.seed(seed)
  at_least_as_extreme <- 0L
  for (index in seq_len(permutations)) {
    permuted <- group_sums_of_squares(
      coords, assignment[sample.int(n)], sizes, total_sum, correction, smaller_levels
    )
    if (pseudo_f(permuted) >= observed_f - 1e-12) {
      at_least_as_extreme <- at_least_as_extreme + 1L
    }
  }
  list(
    statistic = observed_f,
    r_squared = observed_between / ss_total,
    p_value = (at_least_as_extreme + 1) / (permutations + 1),
    permutations = permutations
  )
}

## Multivariate dispersion: distance of each sample to its group centroid,
## then a permutation F-test on those distances.
betadisper_centroid_test <- function(coords, groups, permutations, seed) {
  groups <- factor(groups)
  n <- nrow(coords)
  a <- nlevels(groups)
  assignment <- as.integer(groups)
  squared_norms <- rowSums(coords^2)

  centroid_distances <- function(assign) {
    squared <- numeric(n)
    for (level_index in seq_len(a)) {
      rows <- which(assign == level_index)
      block <- coords[rows, , drop = FALSE]
      centroid <- colSums(block) / length(rows)
      squared[rows] <- squared_norms[rows] -
        2 * as.vector(block %*% centroid) + sum(centroid^2)
    }
    sqrt(pmax(squared, 0))
  }
  dispersion_f <- function(distances, assign) {
    grand_mean <- mean(distances)
    between <- sum(vapply(seq_len(a), function(level_index) {
      values <- distances[assign == level_index]
      length(values) * (mean(values) - grand_mean)^2
    }, numeric(1)))
    within <- sum(vapply(seq_len(a), function(level_index) {
      values <- distances[assign == level_index]
      sum((values - mean(values))^2)
    }, numeric(1)))
    (between / (a - 1L)) / (within / (n - a))
  }

  observed_distances <- centroid_distances(assignment)
  observed_f <- dispersion_f(observed_distances, assignment)
  set.seed(seed)
  at_least_as_extreme <- 0L
  for (index in seq_len(permutations)) {
    shuffled <- assignment[sample.int(n)]
    if (dispersion_f(observed_distances, shuffled) >= observed_f - 1e-12) {
      at_least_as_extreme <- at_least_as_extreme + 1L
    }
  }
  list(
    statistic = observed_f,
    p_value = (at_least_as_extreme + 1) / (permutations + 1),
    permutations = permutations,
    mean_distance = tapply(observed_distances, groups, mean)
  )
}

## Q = sum_k lambda_k R^2_k / sum_k lambda_k over the first `k_components` PCs,
## where R^2_k is the one-way R^2 of PC_k against the group label.
variance_weighted_batch_score <- function(scores, eigenvalues, groups, k_components) {
  groups <- factor(groups)
  sizes <- as.integer(table(groups))
  assignment <- as.integer(groups)
  smaller_levels <- utils::head(order(sizes), nlevels(groups) - 1L)
  component_r2 <- vapply(seq_len(k_components), function(component) {
    column <- scores[, component, drop = FALSE]
    total_sum <- colSums(column)
    correction <- sum(total_sum^2) / nrow(column)
    ss_total <- sum(column^2) - correction
    if (!is.finite(ss_total) || ss_total <= 0) return(0)
    between <- group_sums_of_squares(
      column, assignment, sizes, total_sum, correction, smaller_levels
    )
    max(0, min(1, between / ss_total))
  }, numeric(1))
  weights <- eigenvalues[seq_len(k_components)]
  sum(weights * component_r2) / sum(weights)
}

## Runs both metrics at both stages on ONE fixed gene set so the before/after

## ---------------------------------------------------------------------------
## Feature spaces, PCA and cross-cohort batch-effect quantification.
## ---------------------------------------------------------------------------

## Do two factors induce the same partition of the samples?
same_partition <- function(a, b) {
  a <- droplevels(factor(a))
  b <- droplevels(factor(b))
  nlevels(a) == nlevels(b) &&
    length(unique(paste(a, b, sep = "\r"))) == nlevels(a)
}

describe_levels <- function(groups) {
  sizes <- table(droplevels(factor(groups)))
  paste(sprintf("%s(%d)", names(sizes), as.integer(sizes)), collapse = ", ")
}

## Predictors for the statistical tests. `primary` is always the exact factor
## supplied to ComBat-seq. `secondary` (e.g. study/source) is added only when it
## partitions the samples differently, so cohorts where the two coincide do not
## get duplicated rows.
resolve_predictors <- function(coldata, sample_order, primary_column,
                               secondary_column = NULL) {
  primary <- factor(coldata[sample_order, primary_column])
  predictors <- list()
  predictors[[primary_column]] <- primary
  if (!is.null(secondary_column) && secondary_column %in% colnames(coldata)) {
    secondary <- factor(coldata[sample_order, secondary_column])
    if (!same_partition(primary, secondary)) {
      predictors[[secondary_column]] <- secondary
    } else {
      cat(sprintf(
        "[INFO] Predictor '%s' induces the same partition as '%s'; not duplicated\n",
        secondary_column, primary_column
      ))
    }
  }
  for (name in names(predictors)) {
    cat(sprintf("[INFO] Predictor '%s': %s\n", name, describe_levels(predictors[[name]])))
  }
  predictors
}

## The three fixed feature spaces. Every space excludes only non-finite and
## zero-variance rows, is defined once, and is reused unchanged for both stages
## with a single canonical row order.
build_feature_spaces <- function(joint_pre, joint_post, top_genes,
                                 protein_coding_ids = NULL) {
  shared <- rownames(joint_pre)[rownames(joint_pre) %in% rownames(joint_post)]
  counts <- c(
    rows_pre = nrow(joint_pre),
    rows_post = nrow(joint_post),
    shared_rows = length(shared)
  )
  pre_shared <- joint_pre[shared, , drop = FALSE]
  post_shared <- joint_post[shared, , drop = FALSE]
  finite_rows <- is.finite(rowSums(pre_shared)) & is.finite(rowSums(post_shared))
  counts["non_finite_removed"] <- sum(!finite_rows)
  variable_rows <- matrixStats::rowVars(pre_shared) > 0 &
    matrixStats::rowVars(post_shared) > 0
  counts["zero_variance_removed"] <- sum(finite_rows & !variable_rows)
  usable <- shared[finite_rows & variable_rows]
  counts["usable_rows"] <- length(usable)

  ## Ranked once, on the pre-correction matrix only.
  usable_variance <- matrixStats::rowVars(joint_pre[usable, , drop = FALSE])
  ranked <- usable[order(usable_variance, decreasing = TRUE)]
  top_ids <- utils::head(ranked, min(top_genes, length(ranked)))
  ## Restore canonical (matrix) order rather than variance order.
  top_ids <- usable[usable %in% top_ids]
  counts["top3000"] <- length(top_ids)

  spaces <- list(
    top3000 = list(
      ids = top_ids,
      label = sprintf("Top %s pre-correction variable genes", figure_count(length(top_ids)))
    ),
    full_expression = list(
      ids = usable,
      label = sprintf("Full shared expression matrix (%s genes)", figure_count(length(usable)))
    )
  )
  if (!is.null(protein_coding_ids)) {
    coding_ids <- usable[usable %in% protein_coding_ids]
    counts["protein_coding"] <- length(coding_ids)
    spaces$protein_coding <- list(
      ids = coding_ids,
      label = sprintf("Protein-coding genes (%s)", figure_count(length(coding_ids)))
    )
  }
  spaces <- spaces[c("top3000", "protein_coding", "full_expression")]
  spaces <- spaces[!vapply(spaces, is.null, logical(1))]
  for (name in names(counts)) {
    cat(sprintf("[INFO] feature spaces | %-24s %s\n", name, figure_count(counts[[name]])))
  }
  list(spaces = spaces, counts = counts)
}

## One PCA per (feature space, stage). The complete non-zero score space feeds
## the statistics; only PC1/PC2 are ever plotted.
compute_stage_pca <- function(mat, ids) {
  block <- assert_dimensionality_reduction_input(mat, ids, "pca")
  fit <- stats::prcomp(t(block), center = TRUE, scale. = FALSE)
  keep <- which(fit$sdev > fit$sdev[1L] * sqrt(.Machine$double.eps))
  if (length(keep) < 2L) stop("PCA produced fewer than two non-zero components")
  variance_pct <- fit$sdev^2 / sum(fit$sdev^2) * 100
  list(
    coords = fit$x[, 1:2, drop = FALSE],
    scores = fit$x[, keep, drop = FALSE],
    n_pcs = length(keep),
    xlab = sprintf("PC1 (%.1f%% variance)", variance_pct[1L]),
    ylab = sprintf("PC2 (%.1f%% variance)", variance_pct[2L])
  )
}

## Euclidean PERMANOVA + centroid dispersion for every feature space, stage and
## predictor. Returns the tidy table and the PCA objects for plotting, so the
## figures and the statistics can never diverge.
quantify_batch_effect <- function(cohort, joint_pre, joint_post, feature_spaces,
                                  predictors, master_seed, permutations,
                                  tumour_count = NA_integer_,
                                  dsmz_count = NA_integer_) {
  results <- list()
  embeddings <- list()
  for (space_name in names(feature_spaces)) {
    ids <- feature_spaces[[space_name]]$ids
    cat(sprintf(
      "[INFO] %s | feature space '%s': %s genes\n",
      cohort, space_name, figure_count(length(ids))
    ))
    ## Fail closed before any statistic is computed on this space.
    verify_feature_manifest(joint_pre, joint_post, ids, space_name)
    stage_pca <- list(
      before = compute_stage_pca(joint_pre, ids),
      after = compute_stage_pca(joint_post, ids)
    )
    embeddings[[space_name]] <- stage_pca
    for (stage_name in names(stage_pca)) {
      scores <- stage_pca[[stage_name]]$scores
      for (predictor_name in names(predictors)) {
        groups <- predictors[[predictor_name]]
        ## One derived seed per (cohort, feature space, predictor, stage,
        ## method). Independent of loop order, so adding or removing a feature
        ## space never shifts another cell's permutation stream.
        permanova_seed <- derive_seed(
          master_seed, cohort, space_name, predictor_name, stage_name, "permanova"
        )
        dispersion_seed <- derive_seed(
          master_seed, cohort, space_name, predictor_name, stage_name, "betadisper"
        )
        permanova <- permanova_euclidean(scores, groups, permutations, permanova_seed)
        dispersion <- betadisper_centroid_test(scores, groups, permutations, dispersion_seed)
        cat(sprintf(
          "[INFO] %s | %-18s | %-6s | %-8s PERMANOVA F=%.3f R2=%.4f p=%.4g seed=%d; dispersion F=%.3f p=%.4g seed=%d\n",
          cohort, space_name, stage_name, predictor_name,
          permanova$statistic, permanova$r_squared, permanova$p_value, permanova_seed,
          dispersion$statistic, dispersion$p_value, dispersion_seed
        ))
        base_row <- data.frame(
          cohort = cohort,
          feature_space = space_name,
          stage = stage_name,
          predictor = predictor_name,
          predictor_levels = describe_levels(groups),
          n_samples = length(groups),
          tumour_count = as.integer(tumour_count),
          dsmz_count = as.integer(dsmz_count),
          n_genes = length(ids),
          n_pcs = stage_pca[[stage_name]]$n_pcs,
          distance = "euclidean",
          stringsAsFactors = FALSE
        )
        results[[length(results) + 1L]] <- data.frame(
          base_row,
          method = "permanova",
          statistic_name = "pseudo_F",
          statistic = round(permanova$statistic, 6),
          r_squared = round(permanova$r_squared, 6),
          p_value = permanova$p_value,
          permutations = permanova$permutations,
          seed = permanova_seed,
          notes = "vegan::adonis2-equivalent, Euclidean, closed form",
          stringsAsFactors = FALSE
        )
        results[[length(results) + 1L]] <- data.frame(
          base_row,
          method = "betadisper",
          statistic_name = "F",
          statistic = round(dispersion$statistic, 6),
          r_squared = NA_real_,
          p_value = dispersion$p_value,
          permutations = dispersion$permutations,
          seed = dispersion_seed,
          notes = sprintf(
            "centroid dispersion; mean distance %s",
            paste(names(dispersion$mean_distance),
                  sprintf("%.3f", dispersion$mean_distance), sep = "=", collapse = " ")
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  table_out <- do.call(rbind, results)
  ## Relative R^2 reduction, PERMANOVA only. Deliberately not computed for the
  ## betadisper F-statistic, which is not a variance-explained quantity.
  table_out$r_squared_reduction_pct <- NA_real_
  permanova_rows <- table_out$method == "permanova"
  key <- paste(table_out$feature_space, table_out$predictor, sep = "@")
  for (this_key in unique(key[permanova_rows])) {
    before_row <- which(permanova_rows & key == this_key & table_out$stage == "before")
    after_row <- which(permanova_rows & key == this_key & table_out$stage == "after")
    if (length(before_row) != 1L || length(after_row) != 1L) next
    before_value <- table_out$r_squared[before_row]
    after_value <- table_out$r_squared[after_row]
    if (is.finite(before_value) && before_value > 0) {
      table_out$r_squared_reduction_pct[after_row] <-
        round(100 * (before_value - after_value) / before_value, 4)
    }
  }
  column_order <- c(
    "cohort", "feature_space", "stage", "predictor", "predictor_levels",
    "n_samples", "tumour_count", "dsmz_count",
    "n_genes", "n_pcs", "distance", "method", "statistic_name",
    "statistic", "r_squared", "p_value", "permutations", "seed",
    "r_squared_reduction_pct", "notes"
  )
  table_out <- table_out[
    order(table_out$feature_space, table_out$predictor, table_out$method, table_out$stage),
    column_order,
    drop = FALSE
  ]
  list(table = table_out, embeddings = embeddings)
}

## Paired before/after PCA for one feature space, using the shared grid-free
## design. Only PC1/PC2 are displayed; the statistics used every non-zero PC.
plot_feature_space_pca <- function(stage_pca, labels, path, cohort_title,
                                   sample_summary, feature_label, annotation = NULL) {
  plot_embedding_figure(
    list(
      list(embedding = stage_pca$before, main = "A  Before correction"),
      list(embedding = stage_pca$after, main = "B  After correction")
    ),
    labels,
    path,
    overall_title = cohort_title,
    overall_subtitle = sample_summary,
    footer = feature_label,
    annotation = annotation,
    scale_style = "numeric"
  )
}


## Renamed from plot_mean_sd to avoid colliding with the NBL/RBL
## visualisation helpers of the same name, which have a different signature.
plot_vst_mean_sd <- function(mat, title, path, seed) {
  means <- rowMeans(mat)
  sds <- matrixStats::rowSds(mat)
  keep <- which(is.finite(means) & is.finite(sds))
  if (length(keep) > 10000L) {
    set.seed(seed)
    keep <- sample(keep, 10000L)
  }
  figure_qc_scatter(
    means[keep],
    sds[keep],
    title = title,
    subtitle = sprintf(
      "%s genes shown  |  colour = local point density  |  orange line = running median",
      figure_count(length(keep))
    ),
    xlab = "Mean VST value",
    ylab = "Row standard deviation",
    path = path
  )
}

plot_count_dispersion <- function(counts, title, path, seed) {
  means <- rowMeans(counts)
  variances <- matrixStats::rowVars(counts)
  keep <- which(is.finite(means) & is.finite(variances) & means > 0)
  if (length(keep) > 10000L) {
    set.seed(seed)
    keep <- sample(keep, 10000L)
  }
  dispersion <- pmax((variances[keep] - means[keep]) / (means[keep]^2), 0)
  figure_qc_scatter(
    means[keep],
    dispersion,
    title = title,
    subtitle = sprintf(
      "%s genes shown  |  colour = local point density  |  orange line = running median",
      figure_count(length(keep))
    ),
    xlab = "Mean raw count (log scale)",
    ylab = "Method-of-moments dispersion",
    path = path,
    log = "x"
  )
}
