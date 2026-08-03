#!/usr/bin/env Rscript

# plot_pan_cancer_tumour_cell_line_alignment_umap.R
# ==============================================================================
# PAN-CANCER TUMOUR/CELL-LINE ALIGNMENT UMAP
# ==============================================================================
#
# SCRIPT OVERVIEW
# ---------------
# This script generates UMAP (Uniform Manifold Approximation and Projection)
# embeddings of pan-cancer transcriptomic data using the current graph-informed
# DESeq2 marker-derived feature panel. The visualisations assess whether
# cell-line profiles occupy the broad transcriptomic regions of their annotated
# tumour lineages.
#
# BIOLOGICAL CONTEXT
# ------------------
# Pan-cancer analysis examines molecular patterns that transcend individual
# cancer types, revealing both shared oncogenic programmes and lineage-specific
# biology. This alignment visualisation plots breast cancer, neuroblastoma,
# and retinoblastoma samples; HEME profiles may be present in upstream inputs
# but are not included in this figure.
#
#   1. Do cell lines cluster with tumours of their annotated cancer type?
#   2. Are there cross-cancer molecular similarities suggesting shared biology?
#   3. Which cell lines may have diverged from their tissue of origin?
#
# THE FEATURE PANEL
# -----------------
# Rather than using all expressed genes or generic highly variable genes, this
# script uses the graph-informed DESeq2 marker-derived pan-cancer feature panel
# assembled upstream from recurrent genes, accepted singleton candidates, and
# accepted non-recurrent candidates.
#
# UMAP METHODOLOGY
# ----------------
# UMAP is a non-linear dimensionality reduction technique that preserves
# both local neighbourhood structure and global topology. Key advantages
# for cancer transcriptomics:
#
#   - Preserves cluster structure better than t-SNE for well-separated groups
#   - More computationally efficient than t-SNE for large datasets
#   - Produces more reproducible embeddings with fixed random seeds
#   - Better preserves global relationships between distant clusters
#
# The script evaluates multiple distance metrics (cosine, Euclidean) to
# assess robustness of the observed patterns to methodological choices.
#
# VISUALISATION DESIGN
# --------------------
# Hue encodes lineage only. Sample type is encoded by marker shape/fill/stroke:
# tumours are smaller semi-transparent filled circles; cell lines are larger
# hollow diamonds with lineage-coloured outlines, drawn last.
#
# OUTPUT FILES
# ------------
# For each distance metric, the script generates:
#   - Coordinate table (TSV): UMAP coordinates with sample metadata
#   - Visualisation (PDF/SVG/PNG): Pan-cancer overview plot
#   - Summary table (TSV): Analysis parameters and output file inventory
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The following R packages provide essential functionality for data processing,
# dimensionality reduction, and visualisation.

suppressPackageStartupMessages({
  
  # readr: Fast and consistent functions for reading and writing delimited files.
  # Provides better defaults than base R for character encoding and type inference.
  library(readr)
  
  # dplyr: A grammar of data manipulation with intuitive verbs for transforming
  # data frames. Essential for metadata processing and result aggregation.
  library(dplyr)
  
  # tibble: Modern reimplementation of data frames with more stringent subsetting
  # rules and improved printing. Tibbles never convert strings to factors.
  library(tibble)
  
  # ggplot2: The grammar of graphics implementation for R. Provides a layered
  # approach to building statistical visualisations.
  library(ggplot2)
  
  # stringr: Consistent string manipulation functions built on stringi.
  # Used for sample identifier parsing and cancer type normalisation.
  library(stringr)
  
  # optparse: Systematic command-line argument parsing following GNU conventions.
  # Enables flexible parameterisation without source code modification.
  library(optparse)
  
  # uwot: Fast implementation of UMAP for dimensionality reduction.
  # Supports multiple distance metrics and parallel computation.
  library(uwot)

  # ggrepel: Optional non-overlapping labels for any metadata-supported
  # future cluster annotation. The default plot leaves outlier clusters
  # unannotated unless explicit annotation data are supplied.
  library(ggrepel)
  
  library(scales)
})

# ------------------------------------------------------------------------------
# NULL-COALESCING OPERATOR
# ------------------------------------------------------------------------------
# The %||% operator returns the left operand if it is not NULL, empty, or
# entirely NA; otherwise returns the right operand. This extended version
# handles edge cases common in data processing pipelines.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# ==============================================================================
# SECTION 2: THREAD CONFIGURATION FOR HPC ENVIRONMENTS
# ==============================================================================
# High-performance computing (HPC) clusters allocate specific CPU resources
# to each job. This section configures thread limits to respect these
# allocations, preventing resource contention and ensuring efficient execution.
#
# SLURM (Simple Linux Utility for Resource Management) is a common HPC
# workload manager. The SLURM_CPUS_PER_TASK environment variable indicates
# the number of CPUs allocated to the current job.
#
# Multiple numerical libraries may spawn threads independently. Setting
# consistent limits across all threading backends prevents over-subscription
# where more threads are created than available CPUs.

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
if (is.na(n_threads) || n_threads < 1) n_threads <- 8

# Configure thread limits for various numerical backends.
Sys.setenv(
  OMP_NUM_THREADS = n_threads,         # OpenMP (general parallel computing)
  OPENBLAS_NUM_THREADS = n_threads,    # OpenBLAS (linear algebra)
  MKL_NUM_THREADS = n_threads,         # Intel MKL (linear algebra)
  VECLIB_MAXIMUM_THREADS = n_threads,  # Apple Accelerate framework
  NUMEXPR_NUM_THREADS = n_threads      # NumExpr (numerical expressions)
)

# Configure R-specific parallelisation.
options(mc.cores = n_threads)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = n_threads)

# ==============================================================================
# SECTION 3: COMMAND-LINE ARGUMENT PARSING
# ==============================================================================
# Command-line arguments enable flexible parameterisation of the analysis
# without modifying source code. This is essential for reproducible research
# and automated pipeline integration.

option_list <- list(
  
  # --pipe_root: Root directory of the analysis pipeline.
  # All relative paths are resolved relative to this directory.
  make_option("--pipe_root", type="character", default=".",
              help="Pipeline root [default: %default]"),

  # --expr_rds: Path to the expression matrix containing integrated
  # pan-cancer data. The matrix may be in either genes x samples or
  # samples x genes orientation; the script handles both.
  make_option("--expr_rds", type="character",
              default="results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds",
              help="Joint expression matrix RDS [default: %default]"),

  # --meta_tsv: Sample metadata including cancer type and sample type
  # annotations required for visualisation aesthetics.
  make_option("--meta_tsv", type="character",
              default="results/unsupervised/multicohort_cancer/inputs/joint_metadata.tsv",
              help="Joint metadata TSV [default: %default]"),

  # --source_meta_tsv: Optional metadata table used only to recover cohort /
  # source labels when --meta_tsv does not already contain a cohort column.
  make_option("--source_meta_tsv", type="character",
              default="results/unsupervised/multicohort_cancer/inputs/joint_metadata.tsv",
              help="Metadata TSV used to join cohort/source labels by sample_id when absent from --meta_tsv [default: %default]"),

  # --feature_list: Path to the marker-panel feature list containing one
  # Ensembl gene identifier per line. Required when --feature_mode is
  # feature_list; ignored when --feature_mode is all_genes.
  make_option("--feature_list", type="character",
              default="results/unsupervised/pan_cancer/feature_space/pan_cancer_features_clean.txt",
              help="Feature list file (one ENSG per line); used when --feature_mode=feature_list [default: %default]"),

  # --feature_mode: Whether to subset the expression matrix to a
  # supplied gene list ('feature_list') or use all genes present in the
  # matrix after missing-value and zero-variance filtering ('all_genes').
  make_option("--feature_mode", type="character", default="feature_list",
              help="feature_list (default; marker-panel feature space) or all_genes (use all genes after NA + zero-var filtering)"),

  # --feature_label: Label embedded in output figure / coordinate
  # filenames so a single output directory can hold runs from different
  # feature sets.
  make_option("--feature_label", type="character", default=NULL,
              help="Feature-set label inserted into output filenames (default: PAN_CANCER_MARKER_PANEL when feature_mode=feature_list, ALL_GENES when feature_mode=all_genes)"),

  # --out_stem: Filename stem (without metric / feature_label / .ext)
  # used to compose figure, coordinate, and summary output names.
  # Default preserves the existing pan-cancer alignment figure naming.
  make_option("--out_stem", type="character",
              default="pan_cancer_tumour_cell_line_alignment_umap",
              help="Output filename stem (figures/coords use {stem}_{feature_label}_{metric}; summary uses {stem} or --summary_basename) [default: %default]"),

  # --summary_basename: Explicit basename for the summary TSV. When
  # unset, defaults to summary_{out_stem}.tsv (preserves existing
  # pan-cancer summary filename).
  make_option("--summary_basename", type="character", default=NULL,
              help="Summary TSV basename (default: summary_{out_stem}.tsv)"),

  # --outdir: Output directory for coordinates, plots, and summary files.
  make_option("--outdir", type="character",
              default="results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap",
              help="Output directory [default: %default]"),

  # --dist_metrics: Distance metrics for UMAP computation.
  # Multiple metrics can be specified to assess result robustness.
  make_option("--dist_metrics", type="character",
              default="cosine,euclidean",
              help="Comma-separated uwot metrics [default: %default]"),

  # --n_neighbors: UMAP n_neighbors parameter controlling the balance
  # between local and global structure preservation.
  # Higher values emphasise global topology; lower values capture local detail.
  make_option("--n_neighbors", type="integer", default=30,
              help="UMAP n_neighbors [default: %default]"),
  
  # --min_dist: UMAP min_dist parameter controlling point packing density.
  # Lower values create tighter clusters; higher values spread points evenly.
  make_option("--min_dist", type="double", default=0.3,
              help="UMAP min_dist [default: %default]"),
  
  # --seed: Random seed for reproducibility.
  # UMAP involves stochastic optimisation; fixing the seed ensures
  # identical embeddings across runs.
  make_option("--seed", type="integer", default=42,
              help="Seed [default: %default]"),

  # --page: Output size preset for different presentation contexts.
  # Alignment presets use the established full-width figure size.
  make_option("--page", type="character", default="pan_cancer_alignment",
              help="Output size preset: pan_cancer_alignment, all_gene_alignment, slide (16:9), a4, a4_landscape [default: %default]"),

  make_option("--width", type="double", default=9.0,
              help="Width in inches for alignment page presets [default: %default]"),
  make_option("--height", type="double", default=7.0,
              help="Height in inches for alignment page presets [default: %default]"),
  make_option("--slide_width", type="double", default=13.33,
              help="Width in inches for slide page preset [default: %default]"),
  make_option("--slide_height", type="double", default=7.5,
              help="Height in inches for slide page preset [default: %default]"),
  make_option("--legend_position", type="character", default="bottom",
              help="ggplot legend position: bottom, top, left, right, none [default: %default]"),
  make_option("--source_cancer_slide_threshold", type="integer", default=6,
              help="Use slide page preset for source-cancer UMAPs above this number of legend entries [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
legend_position <- tolower(opt$legend_position)
if (!legend_position %in% c("bottom", "top", "left", "right", "none")) {
  stop("--legend_position must be one of bottom, top, left, right, none", call. = FALSE)
}
legend_direction <- if (legend_position %in% c("left", "right")) "vertical" else "horizontal"
alignment_page_presets <- c(
  "pan_cancer_alignment",
  "all_gene_alignment",
  "thesis"
)

# ------------------------------------------------------------------------------
# PATH RESOLUTION
# ------------------------------------------------------------------------------
# Resolve all paths relative to the pipeline root directory to ensure
# correct file access regardless of working directory.

pipe_root <- normalizePath(opt$pipe_root, mustWork = TRUE)

abs_from_root <- function(p) {
  # abs_from_root(): Converts relative paths to absolute paths.
  #
  # If the path is already absolute, it is returned unchanged.
  # Otherwise, it is resolved relative to the pipeline root.
  
  if (is.null(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(pipe_root, p)
}

# Resolve feature_mode and the feature-list path.
feature_mode <- tolower(trimws(opt$feature_mode %||% "feature_list"))
if (!(feature_mode %in% c("feature_list", "all_genes"))) {
  stop("Invalid --feature_mode: ", feature_mode,
       " (must be 'feature_list' or 'all_genes')")
}

feature_list_path <- opt$feature_list

# Resolve all input and output paths.
expr_rds <- abs_from_root(opt$expr_rds)
meta_tsv <- abs_from_root(opt$meta_tsv)
source_meta_tsv <- abs_from_root(opt$source_meta_tsv)
outdir   <- abs_from_root(opt$outdir)
if (feature_mode == "feature_list") {
  if (is.null(feature_list_path) || !nzchar(feature_list_path)) {
    stop("--feature_mode=feature_list requires --feature_list to be set.")
  }
  feature_list_path <- abs_from_root(feature_list_path)
} else {
  # all_genes: feature-list file is not required and is not consulted.
  feature_list_path <- NA_character_
}

# Resolve feature_label, out_stem, summary_basename with sensible defaults
feature_label <- opt$feature_label
if (is.null(feature_label) || !nzchar(feature_label)) {
  feature_label <- if (feature_mode == "all_genes") "ALL_GENES" else "PAN_CANCER_MARKER_PANEL"
}
out_stem <- opt$out_stem %||% "pan_cancer_tumour_cell_line_alignment_umap"
summary_basename <- opt$summary_basename
if (is.null(summary_basename) || !nzchar(summary_basename)) {
  summary_basename <- paste0("summary_", out_stem, ".tsv")
}

# Create output directory.
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Parse distance metrics from comma-separated string.
dist_metrics <- strsplit(opt$dist_metrics, ",")[[1]] %>% trimws()
dist_metrics <- dist_metrics[nzchar(dist_metrics)]

# Log configuration for reproducibility documentation.
cat("[INFO] expr_rds        :", expr_rds, "\n")
cat("[INFO] meta_tsv        :", meta_tsv, "\n")
cat("[INFO] source_meta_tsv :", source_meta_tsv, "\n")
cat("[INFO] feature_mode    :", feature_mode, "\n")
cat("[INFO] feature_list    :", ifelse(is.na(feature_list_path), "(not used)", feature_list_path), "\n")
cat("[INFO] feature_label   :", feature_label, "\n")
cat("[INFO] out_stem        :", out_stem, "\n")
cat("[INFO] summary_basename:", summary_basename, "\n")
cat("[INFO] outdir          :", outdir, "\n")
cat("[INFO] metrics         :", paste(dist_metrics, collapse=", "), "\n")
cat(sprintf("[INFO] threads         : %d\n", n_threads))

# ==============================================================================
# SECTION 4: UTILITY FUNCTIONS
# ==============================================================================

clean_ensg <- function(ids) {
  # clean_ensg(): Removes version suffixes from Ensembl gene identifiers.
  #
  # Ensembl identifiers may include version numbers (e.g., ENSG00000141510.16).
  # These versions can cause matching failures when different data sources
  # use different annotation versions. This function strips the version
  # suffix to enable robust identifier matching.
  #
  # Parameters:
  #   ids: A character vector of Ensembl gene identifiers
  #
  # Returns:
  #   The identifiers with version suffixes removed
  #
  # Example:
  #   "ENSG00000141510.16" -> "ENSG00000141510"
  
  sub("\\.\\d+$", "", ids)
}

read_feature_list <- function(path) {
  # read_feature_list(): Loads and validates a gene list file.
  #
  # The function reads a text file containing one gene identifier per line,
  # removes comments (lines starting with #), strips whitespace, removes
  # Ensembl version suffixes, and returns unique identifiers.
  #
  # Parameters:
  #   path: Path to the gene list file
  #
  # Returns:
  #   A character vector of cleaned, unique Ensembl gene identifiers
  
  if (!file.exists(path)) stop("feature list not found: ", path)
  
  g <- readLines(path, warn = FALSE)
  g <- g[!grepl("^\\s*#", g)]  # Remove comment lines
  g <- trimws(g)                # Strip whitespace
  g <- g[nzchar(g)]             # Remove empty lines
  g <- clean_ensg(g)            # Remove version suffixes
  g <- unique(g)                # Deduplicate
  
  if (length(g) == 0) stop("feature list empty after cleaning: ", path)
  g
}

as_samples_x_genes <- function(x) {
  # as_samples_x_genes(): Standardises matrix orientation to samples x genes.
  #
  # Expression matrices may be stored in either genes x samples (common in
  # Bioconductor) or samples x genes (required by uwot::umap) orientation.
  # This function detects the current orientation by examining whether
  # row or column names contain Ensembl gene identifiers, and transposes
  # if necessary.
  #
  # Parameters:
  #   x: An expression matrix or data frame
  #
  # Returns:
  #   A matrix with samples as rows and genes as columns
  
  X <- as.matrix(x)
  rn <- rownames(X) %||% character(0)
  cn <- colnames(X) %||% character(0)

  # Detect orientation by checking for ENSG prefix in dimension names.
  rn_has_ensg <- length(rn) > 0 && mean(grepl("^ENSG", rn)) > 0.5
  cn_has_ensg <- length(cn) > 0 && mean(grepl("^ENSG", cn)) > 0.5

  # If genes are in rows, transpose to put genes in columns.
  if (rn_has_ensg && !cn_has_ensg) X <- t(X)

  # Clean gene identifiers in column names.
  if (!is.null(colnames(X))) colnames(X) <- clean_ensg(colnames(X))
  
  X
}

# ==============================================================================
# SECTION 5: VISUALISATION THEME
# ==============================================================================
# A consistent visual theme ensures comparability across all pipeline outputs.
# This theme is optimised for UMAP plots with the following design principles:
#
#   - Clean white background for maximum point visibility
#   - No gridlines (UMAP coordinates are relative, not absolute)
#   - Bottom legends so the UMAP panel occupies most of the figure width
#   - Clear typographic hierarchy for title, subtitle, and caption

theme_umap_main <- function() {
  # theme_umap_main(): Custom ggplot2 theme for pan-cancer UMAP plots.
  #
  # This theme extends theme_minimal() with modifications suited for
  # multi-disease visualisations where many groups must be distinguished.
  #
  # Returns:
  #   A ggplot2 theme object
  
  theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(linewidth = 0.4, colour = "grey20"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text         = element_text(size = 10, colour = "grey20"),
      axis.title        = element_text(size = 11, colour = "grey10"),
      legend.position   = legend_position,
      legend.direction  = legend_direction,
      legend.box        = legend_direction,
      legend.title      = element_text(face = "bold", size = 9),
      legend.text       = element_text(size = 8),
      legend.key.size   = grid::unit(0.72, "lines"),
      legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
      plot.title        = element_text(size = 14, face = "bold", colour = "grey10"),
      plot.subtitle     = element_text(size = 10, colour = "grey40"),
      plot.caption      = element_text(size = 8, colour = "grey50", hjust = 0),
      plot.margin       = margin(14, 16, 10, 10),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

# ==============================================================================
# SECTION 6: DATA PREPARATION FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
lineage_display_labels <- c(
  BRCA = "Breast Cancer",
  NBL  = "Neuroblastoma",
  RBL  = "Retinoblastoma"
)

sample_type_display_label <- function(x) {
  x0 <- as.character(x)
  x_low <- tolower(trimws(x0))
  dplyr::case_when(
    x_low %in% c("tumour", "tumor", "tumours", "tumors") ~ "Tumour",
    x_low %in% c("cell_line", "cell line", "cellline", "cell", "cells") ~ "Cell line",
    TRUE ~ x0
  )
}

sample_type_code <- function(x) {
  x0 <- as.character(x)
  x_low <- tolower(trimws(x0))
  dplyr::case_when(
    x_low %in% c("tumour", "tumor", "tumours", "tumors") ~ "tumour",
    x_low %in% c("cell_line", "cell line", "cellline", "cell", "cells") ~ "cell_line",
    TRUE ~ x0
  )
}

prepare_plot_data <- function(df) {
  df %>%
    mutate(
      type = factor(sample_type_label, levels = c("Tumour", "Cell line")),
      disease = cancer_type_label
    ) %>%
    filter(!is.na(disease)) %>%
    mutate(
      disease = factor(disease, levels = c("Breast Cancer","Neuroblastoma","Retinoblastoma"))
    )
}

# ==============================================================================
# SECTION 7: PLOTTING FUNCTION
# ==============================================================================
# Hue = lineage; marker shape/fill/stroke = sample type.

pan_cancer_alignment_plot <- function(alignment, palette, subtitle_text,
                                      caption_text = NULL,
                                      show_lineage_legend = TRUE,
                                      point_alpha = 0.58,
                                      size_tumour = 1.05, size_cl = 3.0,
                                      stroke_cl = 1.0,
                                      label_df = NULL, annot_df = NULL) {

  df <- alignment %>% mutate(type = factor(type, levels = c("Tumour", "Cell line")))
  type_legend_df <- tibble(
    UMAP1 = rep(min(df$UMAP1, na.rm = TRUE), 2),
    UMAP2 = rep(min(df$UMAP2, na.rm = TRUE), 2),
    stype = factor(c("Tumour", "Cell line"), levels = c("Tumour", "Cell line"))
  )

  gg <- ggplot(df, aes(UMAP1, UMAP2)) +

    # Tumours: filled lineage circles, semi-transparent, no heavy outline.
    geom_point(
      data = df %>% filter(type == "Tumour"),
      aes(fill = disease),
      shape = 21, colour = "transparent", alpha = point_alpha,
      size = size_tumour, stroke = 0
    ) +

    # Cell lines: hollow lineage-coloured diamonds, drawn last.
    geom_point(
      data = df %>% filter(type == "Cell line"),
      aes(colour = disease),
      shape = 23, fill = "white", alpha = 1, size = size_cl, stroke = stroke_cl
    ) +

    scale_fill_manual(
      values = palette,
      name = "Lineage",
      guide = if (show_lineage_legend) "legend" else "none"
    ) +
    scale_colour_manual(values = palette, guide = "none") +

    # Invisible points used only to produce a sample-type legend whose symbols
    # match the plotted tumour and cell-line markers.
    geom_point(
      data = type_legend_df,
      aes(UMAP1, UMAP2, shape = stype),
      inherit.aes = FALSE, colour = "grey35", fill = "white",
      size = 2.8, stroke = 0.9, alpha = 0, show.legend = TRUE
    ) +
    scale_shape_manual(
      values = c(`Tumour` = 21, `Cell line` = 23),
      name = "Sample type",
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          colour = c("grey45", "grey25"),
          fill = c("grey45", "white"),
          size = c(2.5, 3.0),
          stroke = c(0.25, 1.0)
        )
      )
    ) +

    theme_umap_main() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    labs(subtitle = subtitle_text, caption = caption_text)

  if (!is.null(label_df) && nrow(label_df) > 0) {
    gg <- gg + geom_text(
      data = label_df,
      aes(x = x, y = y, label = label, colour = disease, hjust = hjust),
      inherit.aes = FALSE,
      fontface = "bold", size = 4.0, show.legend = FALSE
    )
  }

  if (!is.null(annot_df) && nrow(annot_df) > 0) {
    gg <- gg + ggrepel::geom_text_repel(
      data = annot_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.0, colour = "grey25", segment.colour = "grey60",
      min.segment.length = 0, box.padding = 0.6, max.overlaps = Inf,
      show.legend = FALSE
    )
  }

  if (show_lineage_legend) {
    gg <- gg + guides(
      fill = guide_legend(
        title = "Lineage",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = 21,
          colour = "grey30",
          fill = unname(palette),
          alpha = 1,
          size = 3.2,
          stroke = 0.25
        )
      )
    )
  }

  gg
}

source_diagnostic_plot <- function(df, colour_col, colour_title, subtitle_text,
                                   caption_text = NULL, legend_rows = 2) {
  ggplot(
    df,
    aes(UMAP1, UMAP2, colour = .data[[colour_col]], shape = sample_type_label)
  ) +
    geom_point(alpha = 0.72, size = 1.55, stroke = 0.65) +
    scale_shape_manual(
      values = c(`Tumour` = 16, `Cell line` = 4),
      name = "Sample type",
      guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.8, alpha = 1))
    ) +
    labs(
      colour = colour_title,
      subtitle = subtitle_text,
      caption = caption_text,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_umap_main() +
    guides(
      colour = guide_legend(
        nrow = legend_rows,
        byrow = TRUE,
        override.aes = list(size = 3.0, alpha = 1)
      )
    )
}

# ==============================================================================
# SECTION 8: OUTPUT SAVING FUNCTION
# ==============================================================================

save_plot <- function(plot, path_base, page = "pan_cancer_alignment") {
  # save_plot(): Saves plot in PDF, SVG, and PNG formats.
  #
  # The function supports multiple size presets optimised for different
  # presentation contexts:
  #
  #   pan_cancer_alignment / all_gene_alignment:
  #       full-width alignment figure (9.0" x 7.0")
  #   slide: 16:9 aspect ratio for PowerPoint/Google Slides (13.33" x 7.5")
  #   a4: A4 portrait for print documents (8.27" x 11.69")
  #   a4_landscape: A4 landscape for wide figures (11.69" x 8.27")
  #
  # PDF output uses Cairo for better font rendering and transparency support.
  # PNG output uses 400 DPI for high-resolution presentations.
  #
  # Parameters:
  #   plot:      A ggplot object
  #   path_base: Output path without extension (extensions added automatically)
  #   page:      Size preset name
  
  page <- tolower(page)

  if (page %in% alignment_page_presets) {
    width <- opt$width
    height <- opt$height
  } else if (page == "slide") {
    width <- opt$slide_width
    height <- opt$slide_height
  } else if (page == "a4") {
    width <- 8.27
    height <- 11.69
  } else if (page == "a4_landscape") {
    width <- 11.69
    height <- 8.27
  } else {
    stop(
      "Unknown --page: ", page,
      " (use pan_cancer_alignment, all_gene_alignment, slide, a4, a4_landscape)"
    )
  }

  # Save vector (PDF/SVG) and raster (PNG) versions.
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, 
         units = "in", device = grDevices::cairo_pdf)
  ggsave(paste0(path_base, ".svg"), plot, width = width, height = height,
         units = "in", device = grDevices::svg)
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, 
         units = "in", dpi = 400)
}

page_dimensions <- function(page) {
  page <- tolower(page)
  if (page %in% alignment_page_presets) {
    return(c(width = opt$width, height = opt$height))
  }
  if (page == "slide") {
    return(c(width = opt$slide_width, height = opt$slide_height))
  }
  if (page == "a4") {
    return(c(width = 8.27, height = 11.69))
  }
  if (page == "a4_landscape") {
    return(c(width = 11.69, height = 8.27))
  }
  c(width = NA_real_, height = NA_real_)
}

# ==============================================================================
# SECTION 9: DATA LOADING AND VALIDATION
# ==============================================================================
# This section loads the expression matrix and metadata, performs extensive
# validation, and aligns the two data structures.

# Validate input file existence.
if (!file.exists(expr_rds)) stop("Missing expr_rds: ", expr_rds)
if (!file.exists(meta_tsv)) stop("Missing meta_tsv: ", meta_tsv)

# Load expression matrix.
X0 <- readRDS(expr_rds)

# Handle case where RDS contains a list structure.
# Some pipeline outputs wrap the matrix in a list with an 'expr' element.
if (is.list(X0) && "expr" %in% names(X0)) {
  X0 <- X0$expr
  cat("[INFO] Extracted 'expr' element from list structure\n")
}

# Standardise to samples x genes orientation.
X  <- as_samples_x_genes(X0)
rm(X0)  # Free memory

if (is.null(rownames(X))) stop("Expression matrix has no rownames (sample IDs).")
if (nrow(X) <= 2 || ncol(X) <= 2) {
  stop("Expression matrix looks like a placeholder: ", nrow(X), " samples x ",
       ncol(X), " genes after orientation.")
}
cat("[INFO] Expression matrix after orientation:", nrow(X), "samples x", ncol(X), "genes\n")

# Load sample metadata.
meta <- read_tsv(meta_tsv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# METADATA COLUMN NORMALISATION
# ------------------------------------------------------------------------------
# Different data sources may use different column names for the same concept.
# This section normalises column names to a consistent schema.

if (!("sample_id" %in% names(meta))) stop("meta_tsv must contain sample_id")

# Accept 'lineage' as alternative to 'cancer_type'.
if (!("cancer_type" %in% names(meta)) && ("lineage" %in% names(meta))) {
  meta <- meta %>% rename(cancer_type = lineage)
}

# Accept 'type' as alternative to 'sample_type'.
if (!("sample_type" %in% names(meta)) && ("type" %in% names(meta))) {
  meta <- meta %>% rename(sample_type = type)
}

# Validate required columns.
need <- c("sample_id","cancer_type","sample_type")
miss <- setdiff(need, names(meta))
if (length(miss) > 0) stop("meta_tsv missing required columns: ", paste(miss, collapse=", "))

# ------------------------------------------------------------------------------
# SAMPLE TYPE STANDARDISATION
# ------------------------------------------------------------------------------
# Standardise sample type to code-safe internal values while keeping a
# separate display label for plot legends and coordinate tables.

meta <- meta %>%
  mutate(
    sample_type = sample_type_code(sample_type),
    sample_type_label = sample_type_display_label(sample_type)
  )

# ------------------------------------------------------------------------------
# CANCER TYPE STANDARDISATION
# ------------------------------------------------------------------------------
# Normalise cancer type labels to consistent uppercase abbreviations.
# This handles variations in naming conventions across data sources.

meta <- meta %>%
  mutate(
    cancer_type = stringr::str_trim(cancer_type),
    cancer_type = toupper(cancer_type),
    cancer_type = dplyr::recode(
      cancer_type,
      "NEUROBLASTOMA" = "NBL",
      "RETINOBLASTOMA" = "RBL",
      "BREAST CANCER" = "BRCA",
      "BREAST" = "BRCA",
      .default = cancer_type
    )
  )

# ------------------------------------------------------------------------------
# DATASET SOURCE CORRECTION
# ------------------------------------------------------------------------------
# Recover cohort/source labels when the primary metadata does not carry them.
# The corrected dataset_source keeps SRP409177/PRJNA904244 distinct from TCGA.

if (!("cohort" %in% names(meta))) {
  meta$cohort <- NA_character_
}
meta <- meta %>% mutate(cohort = as.character(cohort))

if (file.exists(source_meta_tsv)) {
  source_meta <- read_tsv(source_meta_tsv, show_col_types = FALSE)
  if (all(c("sample_id", "cohort") %in% names(source_meta))) {
    source_meta <- source_meta %>%
      transmute(sample_id = as.character(sample_id),
                cohort_source_join = as.character(cohort)) %>%
      distinct(sample_id, .keep_all = TRUE)
    meta <- meta %>%
      left_join(source_meta, by = "sample_id") %>%
      mutate(
        cohort = dplyr::coalesce(
          dplyr::na_if(cohort, ""),
          dplyr::na_if(cohort_source_join, "")
        )
      ) %>%
      select(-cohort_source_join)
  } else {
    warning("[WARN] source_meta_tsv lacks sample_id/cohort columns: ", source_meta_tsv)
  }
} else {
  warning("[WARN] source_meta_tsv not found; dataset_source will use cohort when available: ",
          source_meta_tsv)
}

meta <- meta %>%
  mutate(
    cancer_type_label = dplyr::recode(
      cancer_type,
      BRCA = lineage_display_labels[["BRCA"]],
      NBL  = lineage_display_labels[["NBL"]],
      RBL  = lineage_display_labels[["RBL"]],
      .default = cancer_type
    ),
    dataset_source = case_when(
      grepl("^SRP409177_SRR", sample_id) ~ "SRA_SRP409177",
      !is.na(cohort) & nzchar(cohort) ~ cohort,
      TRUE ~ "unknown"
    ),
    source_cancer_label = paste(dataset_source, cancer_type, sep = " | "),
    source_sample_label = paste(dataset_source, sample_type_label, sep = " | ")
  )

# Diagnostic output: verify metadata structure.
cat("\n[CHECK] Raw meta counts (before filtering to X):\n")
print(table(meta$cancer_type, meta$sample_type))
cat("[CHECK] Raw source counts (before filtering to X):\n")
print(meta %>% dplyr::count(dataset_source, cancer_type, sample_type_label))
cat("[CHECK] n meta sample_id:", nrow(meta), "\n")
cat("[CHECK] n unique meta sample_id:", dplyr::n_distinct(meta$sample_id), "\n")

# ------------------------------------------------------------------------------
# EXPRESSION-METADATA ALIGNMENT
# ------------------------------------------------------------------------------
# Ensure that metadata and expression matrix contain the same samples
# in the same order.

# Identify samples in metadata but missing from expression matrix.
missing_in_X <- setdiff(meta$sample_id, rownames(X))
if (length(missing_in_X) > 0) {
  cat("\n[MISSING] In meta but not in X:", length(missing_in_X), "\n")
  print(missing_in_X)
}

# Filter metadata to samples present in expression matrix.
meta <- meta %>% filter(sample_id %in% rownames(X))

# Subset and reorder expression matrix to match metadata.
X <- X[meta$sample_id, , drop = FALSE]

# Ensure perfect alignment.
if (!all(meta$sample_id == rownames(X))) {
  meta <- meta %>% slice(match(rownames(X), sample_id))
}
stopifnot(all(meta$sample_id == rownames(X)))

# Diagnostic output: verify alignment.
cat("\n[CHECK] After filtering meta to X rownames:\n")
print(table(meta$cancer_type, meta$sample_type))
cat("[CHECK] n samples in X:", nrow(X), "\n")

# Additional diagnostics for cell line coverage.
cat("\n[CHECK] Cell lines per cancer_type:\n")
print(meta %>% dplyr::filter(sample_type == "cell_line") %>% 
        dplyr::count(cancer_type, sort=TRUE))

cat("\n[CHECK] Cell lines with technical replicates (collapsed by _lib suffix):\n")
meta %>%
  filter(sample_type == "cell_line", cancer_type %in% c("NBL","RBL")) %>%
  mutate(cell_line_base = stringr::str_replace(sample_id, "_lib.*$", "")) %>%
  group_by(cancer_type) %>%
  summarise(
    n_samples = n(),
    n_unique_cell_lines = n_distinct(cell_line_base),
    .groups = "drop"
  ) %>%
  print()

# ==============================================================================
# SECTION 10: FEATURE-SET RESOLUTION (mode-aware)
# ==============================================================================
# Mode 1 (feature_list, default): subset the expression matrix to the supplied
# gene list (existing pan-cancer DEG-set behaviour; requires every listed gene
# to be present in the matrix).
# Mode 2 (all_genes): use every column of the matrix, after removing genes
# with any missing values and genes with zero variance across samples. No
# external gene list is consulted; no class label, lineage, or sample-type
# annotation is used to select features.

n_genes_dropped_na      <- 0L
n_genes_dropped_zerovar <- 0L

if (feature_mode == "feature_list") {
  genes_fixed   <- read_feature_list(feature_list_path)
  genes_present <- intersect(colnames(X), genes_fixed)

  cat("[INFO] Feature-list mode: requested =", length(genes_fixed), "\n")
  cat("[INFO]                    present in expr =", length(genes_present), "\n")

  missing <- setdiff(genes_fixed, colnames(X))
  if (length(missing) > 0) {
    cat("[WARN] Feature-list genes missing from expr (showing up to 20):\n")
    print(head(missing, 20))
  }

  if (length(genes_present) != length(genes_fixed)) {
    stop("Not all feature-list genes are present in expr. Present=",
         length(genes_present), " Expected=", length(genes_fixed),
         ". Fix upstream gene harmonization or feature-list IDs.")
  }

  X_sub <- X[, genes_present, drop = FALSE]

} else {
  # all_genes mode: take every column of the matrix and drop NA-bearing /
  # zero-variance genes. UMAP fitting is purely on numeric expression values.
  cat("[INFO] all_genes mode: starting from", ncol(X), "genes\n")

  na_cols <- colSums(is.na(X)) > 0
  n_genes_dropped_na <- sum(na_cols)
  if (n_genes_dropped_na > 0) {
    X <- X[, !na_cols, drop = FALSE]
    cat(sprintf("[INFO] Dropped %d genes with any missing values\n",
                n_genes_dropped_na))
  }

  # matrixStats::colVars would be faster but adds a dependency. apply() on
  # the surviving columns is acceptable here because gene counts are typically
  # in the tens of thousands and this is a one-off filter.
  col_vars <- apply(X, 2, stats::var, na.rm = FALSE)
  zv_cols  <- !is.finite(col_vars) | col_vars == 0
  n_genes_dropped_zerovar <- sum(zv_cols)
  if (n_genes_dropped_zerovar > 0) {
    X <- X[, !zv_cols, drop = FALSE]
    cat(sprintf("[INFO] Dropped %d genes with zero variance\n",
                n_genes_dropped_zerovar))
  }

  if (ncol(X) <= 2) {
    stop("After NA + zero-var filtering only ", ncol(X),
         " genes remain in all_genes mode; cannot run UMAP.")
  }

  X_sub <- X
  cat("[INFO] all_genes mode: retained", ncol(X_sub), "genes after filtering\n")
}

cat("[INFO] X_sub:", nrow(X_sub), "samples x", ncol(X_sub), "genes\n")
cat("[INFO] Sample summary (pre-lineage filter):\n")
print(table(meta$cancer_type, meta$sample_type))

# limit to BRCA/NBL/RBL for UMAP and plotting.
keep <- meta$cancer_type %in% c("BRCA", "NBL", "RBL")
meta <- meta[keep, ]
X_sub <- X_sub[meta$sample_id, , drop = FALSE]
cat("[INFO] limited to BRCA/NBL/RBL:", nrow(X_sub), "samples\n")
print(table(meta$cancer_type, meta$sample_type))

non_tcga_ids_labelled_tcga <- meta %>%
  dplyr::filter(dataset_source == "TCGA", !grepl("^TCGA-", sample_id))
if (nrow(non_tcga_ids_labelled_tcga) > 0) {
  warning("[WARN] Non-TCGA sample IDs labelled as TCGA in dataset_source: ",
          nrow(non_tcga_ids_labelled_tcga),
          " (showing up to 20 in the log)")
  print(non_tcga_ids_labelled_tcga %>%
          dplyr::select(sample_id, dataset_source, cohort, cancer_type, sample_type_label) %>%
          head(20))
} else {
  cat("[CHECK] No non-TCGA sample IDs are labelled TCGA in dataset_source.\n")
}

cat("[CHECK] Final source table:\n")
final_source_table <- meta %>%
  dplyr::count(dataset_source, cancer_type, sample_type_label, name = "n") %>%
  dplyr::arrange(dataset_source, cancer_type, sample_type_label)
print(final_source_table)

# ==============================================================================
# SECTION 11: COLOUR PALETTE DEFINITION
# ==============================================================================
# The colour palette uses distinct hues for each cancer type, following
# colour-blind-friendly principles where possible. The colours are selected
# to maximise perceptual distinguishability.

disease_palette <- c(
  "Breast Cancer"     = "#CC79A7",
  "Neuroblastoma"     = "#0072B2",
  "Retinoblastoma"    = "#E69F00"
)

# ==============================================================================
# SECTION 12: UMAP COMPUTATION AND VISUALISATION
# ==============================================================================
# For each distance metric, compute UMAP embedding and generate visualisation.
# This loop structure enables systematic comparison across metrics.

summaries <- list()

for (dm in dist_metrics) {
  metric_label <- gsub("[^A-Za-z0-9]+", "_", dm)
  tag <- paste0(feature_label, "_", metric_label)
  file_stem <- paste0(out_stem, "_", tag)
  cat("\n[INFO] Running UMAP:", tag, "\n")

  # Set seed for reproducibility.
  # UMAP involves stochastic gradient descent; fixing the seed ensures
  # identical embeddings across pipeline runs.
  set.seed(opt$seed)
  
  # Compute UMAP embedding.
  # Parameters:
  #   X: samples x genes matrix (subset to DEG features)
  #   n_neighbors: number of neighbours for local structure (30)
  #   min_dist: minimum distance between embedded points (0.3)
  #   metric: distance metric for neighbour calculations
  #   n_threads: parallel threads for computation
  #   scale: FALSE to preserve VST-normalised scale
  emb <- uwot::umap(
    X           = X_sub,
    n_neighbors = opt$n_neighbors,
    min_dist    = opt$min_dist,
    metric      = dm,
    n_threads   = n_threads,
    scale       = FALSE,
    verbose     = TRUE
  )

  # Format coordinates as data frame.
  coords <- as.data.frame(emb)
  colnames(coords) <- c("UMAP1","UMAP2")
  coords$sample_id <- rownames(X_sub)

  # Merge with metadata for plotting.
  coords <- coords %>%
    left_join(meta, by = "sample_id") %>%
    mutate(
      metric = dm,
      feature_label = feature_label,
      feature_mode = feature_mode
    ) %>%
    select(
      sample_id, UMAP1, UMAP2,
      cancer_type, cancer_type_label,
      sample_type, sample_type_label,
      any_of("cohort"),
      dataset_source, source_cancer_label, source_sample_label,
      metric, feature_label, feature_mode,
      everything()
    )
  
  # Save coordinate table.
  coords_path <- file.path(outdir, paste0("coords_", file_stem, ".tsv"))
  write_tsv(coords, coords_path)

  source_composition_path <- file.path(outdir, paste0("source_composition_", file_stem, ".tsv"))
  source_composition <- bind_rows(
    coords %>%
      dplyr::count(dataset_source, cancer_type, sample_type_label, name = "n") %>%
      mutate(
        summary_scope = "dataset_source_cancer_sample_type",
        source_cancer_label = NA_character_
      ) %>%
      select(summary_scope, dataset_source, cancer_type, source_cancer_label,
             sample_type_label, n),
    coords %>%
      dplyr::count(source_cancer_label, sample_type_label, name = "n") %>%
      mutate(
        summary_scope = "source_cancer_sample_type",
        dataset_source = NA_character_,
        cancer_type = NA_character_
      ) %>%
      select(summary_scope, dataset_source, cancer_type, source_cancer_label,
             sample_type_label, n)
  ) %>%
    mutate(metric = dm, feature_label = feature_label, feature_mode = feature_mode) %>%
    select(metric, feature_label, feature_mode, everything())
  write_tsv(source_composition, source_composition_path)

  # Generate visualisation.
  plot_df <- prepare_plot_data(coords)
  x_limits <- range(plot_df$UMAP1, na.rm = TRUE)
  x_span <- diff(x_limits)
  x_mid <- mean(x_limits)

  label_df <- plot_df %>%
    dplyr::filter(type == "Tumour") %>%
    dplyr::group_by(disease) %>%
    dplyr::summarise(
      x = stats::median(UMAP1),
      y = stats::quantile(UMAP2, 0.95) + 1.0,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = dplyr::recode(
        as.character(disease),
        "Breast Cancer" = "Breast Cancer",
        "Neuroblastoma" = "Neuroblastoma",
        "Retinoblastoma" = "Retinoblastoma"
      ),
      disease = factor(disease, levels = levels(plot_df$disease))
    ) %>%
    dplyr::mutate(
      hjust = dplyr::case_when(
        x <= x_mid - 0.25 * x_span ~ 0,
        x >= x_mid + 0.25 * x_span ~ 1,
        TRUE ~ 0.5
      ),
      x = dplyr::case_when(
        hjust == 0 ~ pmax(x, x_limits[1] + 0.4),
        hjust == 1 ~ pmin(x, x_limits[2] - 0.4),
        TRUE ~ pmin(pmax(x, x_limits[1] + 0.8), x_limits[2] - 0.8)
      )
    )

  # Do not annotate detached clusters by default. Cluster labels should only be
  # enabled after coordinate-table inspection supports a metadata-backed label.
  annot_df <- NULL

  subtitle <- sprintf("%s (n=%d genes); metric=%s; nn=%d; min_dist=%.2f",
                      feature_label, ncol(X_sub), dm,
                      opt$n_neighbors, opt$min_dist)

  p <- pan_cancer_alignment_plot(
    alignment     = plot_df,
    palette       = disease_palette,
    subtitle_text = subtitle,
    caption_text  = NULL,
    show_lineage_legend = nrow(label_df) == 0,
    label_df      = label_df,
    annot_df      = annot_df
  )
  plot_base <- file.path(outdir, paste0("Fig_", file_stem))
  save_plot(p, plot_base, page = opt$page)

  source_plot_base <- file.path(outdir, paste0("Fig_", out_stem, "_", feature_label, "_SOURCE_", metric_label))
  p_source <- source_diagnostic_plot(
    df = coords,
    colour_col = "dataset_source",
    colour_title = "Dataset source",
    subtitle_text = subtitle,
    legend_rows = 2
  )
  save_plot(p_source, source_plot_base, page = opt$page)

  source_cancer_plot_base <- file.path(outdir, paste0("Fig_", out_stem, "_", feature_label, "_SOURCE_CANCER_", metric_label))
  source_cancer_page <- if (dplyr::n_distinct(coords$source_cancer_label) > opt$source_cancer_slide_threshold) "slide" else opt$page
  p_source_cancer <- source_diagnostic_plot(
    df = coords,
    colour_col = "source_cancer_label",
    colour_title = "Dataset source | lineage",
    subtitle_text = subtitle,
    legend_rows = 3
  )
  save_plot(p_source_cancer, source_cancer_plot_base, page = source_cancer_page)

  # Record summary for this metric.
  main_dims <- page_dimensions(opt$page)
  source_cancer_dims <- page_dimensions(source_cancer_page)
  summaries[[dm]] <- tibble(
    feature_set            = feature_label,
    feature_mode           = feature_mode,
    dist_metric            = dm,
    n_genes                = ncol(X_sub),
    n_samples              = nrow(X_sub),
    n_tumours              = sum(meta$sample_type == "tumour"),
    n_cell_lines           = sum(meta$sample_type == "cell_line"),
    n_dataset_sources      = dplyr::n_distinct(meta$dataset_source),
    n_non_tcga_ids_labelled_tcga = nrow(non_tcga_ids_labelled_tcga),
    n_genes_dropped_na     = n_genes_dropped_na,
    n_genes_dropped_zerovar= n_genes_dropped_zerovar,
    coords_tsv             = basename(coords_path),
    source_composition_tsv = basename(source_composition_path),
    plot_pdf               = basename(paste0(plot_base, ".pdf")),
    plot_svg               = basename(paste0(plot_base, ".svg")),
    plot_png               = basename(paste0(plot_base, ".png")),
    source_plot_pdf        = basename(paste0(source_plot_base, ".pdf")),
    source_plot_svg        = basename(paste0(source_plot_base, ".svg")),
    source_plot_png        = basename(paste0(source_plot_base, ".png")),
    source_cancer_plot_pdf = basename(paste0(source_cancer_plot_base, ".pdf")),
    source_cancer_plot_svg = basename(paste0(source_cancer_plot_base, ".svg")),
    source_cancer_plot_png = basename(paste0(source_cancer_plot_base, ".png")),
    width_in               = unname(main_dims[["width"]]),
    height_in              = unname(main_dims[["height"]]),
    source_cancer_width_in = unname(source_cancer_dims[["width"]]),
    source_cancer_height_in= unname(source_cancer_dims[["height"]])
  )

  cat("[SAVED] ", coords_path, "\n", sep="")
  cat("[SAVED] ", source_composition_path, "\n", sep="")
  cat("[SAVED] ", paste0(plot_base, ".pdf/.svg/.png"), "\n", sep="")
  cat("[SAVED] ", paste0(source_plot_base, ".pdf/.svg/.png"), "\n", sep="")
  cat("[SAVED] ", paste0(source_cancer_plot_base, ".pdf/.svg/.png"), "\n", sep="")
}

# ==============================================================================
# SECTION 13: SUMMARY OUTPUT
# ==============================================================================
# Compile and save a summary table documenting all analyses performed.

summary_tbl <- bind_rows(summaries)
summary_path <- file.path(outdir, summary_basename)
write_tsv(summary_tbl, summary_path)

cat("\n[SAVED] ", summary_path, "\n", sep="")
cat("\n[SUCCESS] UMAP plotting complete.\n")
