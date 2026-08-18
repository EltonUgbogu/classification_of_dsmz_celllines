#!/usr/bin/env Rscript

# unsup_plot_helpers.R

# ==============================================================================
# HELPER FUNCTIONS FOR UNSUPERVISED ANALYSIS AND DIMENSIONALITY REDUCTION
# ==============================================================================
#
# LIBRARY OVERVIEW
# ----------------
# This script provides a collection of helper functions that support the
# unsupervised clustering and tumour neighbourhood analysis pipeline. The
# functions fall into three categories:
#
#   1. Method Configuration: Functions for parsing and managing analytical
#      method identifiers (feature sets, distance metrics) that define
#      the multi-directional ensemble approach.
#
#   2. Visualisation Theme: A consistent ggplot2 theme for dimensionality
#      reduction plots, ensuring visual coherence across all pipeline outputs.
#
#   3. Dimensionality Reduction Helpers: Functions for generating PCA and
#      UMAP visualisations that overlay cell lines onto tumour-derived
#      coordinate systems.
#
# DESIGN PHILOSOPHY
# -----------------
# The pipeline evaluates multiple combinations of feature selection methods
# and distance metrics to identify robust biological patterns. This library
# centralises the logic for handling these "method" or "direction" identifiers,
# ensuring consistent naming conventions and file path generation throughout
# the analysis.
#
# The visualisation functions implement a standardised approach to displaying
# high-dimensional transcriptomic data in two dimensions, with particular
# attention to the visual distinction between primary tumour samples and
# cell line overlays.
#
# COLOUR CONVENTIONS
# ------------------
# The library defines two colour palettes:
#
#   PAM50 colours: For breast cancer molecular subtypes (Basal, Her2, LumA,
#   LumB, Normal). These colours follow conventions established in the
#   breast cancer literature for immediate recognition by domain experts.
#
#   Cluster colours: For numerically-labelled clusters from unsupervised
#   analysis. These use a divergent palette to maximise visual distinction
#   between adjacent cluster numbers.
#
# USAGE CONTEXT
# -------------
# This library is sourced by multiple R scripts in the pipeline, including:
#   - Consensus clustering scripts
#   - Tumour neighbourhood computation
#   - Quality control and visualisation scripts
#   - Pan-cancer analysis modules
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: PACKAGE LOADING
# ------------------------------------------------------------------------------
# The following packages provide essential functionality for configuration
# parsing, data manipulation, and visualisation.

suppressPackageStartupMessages({
  
  # yaml: Provides functions for reading YAML configuration files.
  # YAML format is used throughout the pipeline to specify analysis parameters,
  # file paths, and method definitions in a human-readable format.
  library(yaml)
  
  # dplyr: Provides a grammar of data manipulation with intuitive verbs
  # (filter, select, mutate, summarise, arrange, join) for transforming

  # data frames. Essential for metadata handling and result aggregation.
  library(dplyr)
  
  # ggplot2: The grammar of graphics implementation for R. Provides a
  # layered approach to building statistical visualisations where each
  # component (data, aesthetics, geometries, scales, themes) is specified
  # independently.
  library(ggplot2)
  
  # uwot: A fast implementation of the Uniform Manifold Approximation and
  # Projection (UMAP) algorithm for dimensionality reduction. UMAP preserves
  # both local and global structure better than t-SNE while being computationally
  # more efficient for large datasets.
  library(uwot)
  
 # ggrepel: Extends ggplot2 with text and label geoms that automatically
  # repel overlapping labels. This is essential for creating readable
  # scatter plots where many points require annotation.
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# NULL-COALESCING OPERATOR
# ------------------------------------------------------------------------------
# The %||% operator provides a concise syntax for handling NULL values.
# It returns the left operand if it is not NULL, otherwise returns the
# right operand. This pattern is particularly useful for providing default
# values when optional parameters are not specified.
#
# Example usage:
#   dataset_type <- meta$dataset_type %||% "Tumour"
# If meta$dataset_type is NULL, "Tumour" is used as the default.

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# DATASET TYPE NORMALISATION
# ------------------------------------------------------------------------------
# Normalises dataset type labels to a standardised format ("Tumour" or "Cell line").
# This function handles variations in metadata labelling across different cohorts
# (e.g., "TARGET", "TCGA", "Tumor", "Primary" -> "Tumour"; "DSMZ", "cellline" -> "Cell line").
#
# This ensures that plotting functions work correctly regardless of how the
# metadata labels tumour samples and cell lines, making the pipeline cohort-independent.

normalise_dataset_type <- function(x) {
  # normalise_dataset_type(): Standardises dataset type labels.
  #
  # Converts various dataset type labels to a standardised format:
  #   - Tumour variants: "tumour", "tumor", "primary", "tcga", "target", etc. -> "Tumour"
  #   - Cell line variants: "cell line", "cell_line", "cellline", "dsmz", "cell" -> "Cell line"
  #
  # Parameters:
  #   x: A character vector of dataset type labels (may include NA or NULL)
  #
  # Returns:
  #   A character vector with normalised labels ("Tumour" or "Cell line")
  
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("tumour", "tumor", "primary", "tcga", "target", "tumour sample", "tumor sample")] <- "Tumour"
  x[x %in% c("cell line", "cell_line", "cellline", "dsmz", "cell")] <- "Cell line"
  x
}

# ==============================================================================
# SECTION 2: METHOD CONFIGURATION FUNCTIONS
# ==============================================================================
# The pipeline employs a multi-directional ensemble strategy, systematically
# evaluating combinations of feature selection methods and distance metrics.
# Each combination is identified by a "method_id" or "direction" string
# (e.g., "pam50_euc", "HVG_corr", "Variance_euclidean").
#
# These functions provide utilities for:
#   - Parsing method identifiers from command-line arguments
#   - Generating standardised filename tags
#   - Retrieving full method context from configuration files

get_method_id <- function(default = "pam50_euc") {
  # get_method_id(): Extracts the method identifier from command-line arguments.
  #
  # This function searches the command-line arguments for a --method flag
  # and extracts its value. If no --method argument is found, the specified

  # default value is returned.
  #
  # The method_id encodes both the feature set and distance metric in a
  # single identifier, following the convention: {feature}_{distance}
  #
  # Parameters:
 #   default: The method_id to use if none is specified (default: "pam50_euc")
  #
  # Returns:
  #   A character string containing the method identifier
  #
  # Example command-line usage:
  #   Rscript analysis.R --method=HVG_corr
  #   Rscript analysis.R --method HVG_corr
  
  args <- commandArgs(trailingOnly = TRUE)
  
  # Search for arguments matching the --method pattern.
  m_arg <- args[grepl("^--method", args)]
  
  if (length(m_arg) == 0) return(default)
  
  # Extract the value, handling both --method=value and --method value formats.
  m <- sub("^--method(=)?", "", m_arg[1])
  
  if (m == "") default else m
}

make_method_tags <- function(feature, distance) {
  # make_method_tags(): Generates standardised filename tags from method components.
  #
  # To ensure consistent file naming across the pipeline, this function
  # converts feature set and distance metric names to uppercase tag format.
  # These tags are embedded in output filenames to enable unambiguous
  # identification of analytical provenance.
  #
  # Parameters:
  #   feature:  The feature set name (e.g., "pam50", "HVG", "Variance")
  #   distance: The distance metric name (e.g., "correlation", "euclidean")
  #
  # Returns:
  #   A list with two elements:
  #     - feature_tag: Uppercase feature identifier (e.g., "PAM50", "HVG")
  #     - dist_tag: Abbreviated distance tag ("CORR" or "EUC")
  #
  # NAMING CONVENTIONS:
  #   Feature tags preserve the original capitalisation in uppercase form.
  #   Distance tags use standardised abbreviations:
  #     - "correlation" / "corr" -> "CORR"
  #     - "euclidean" / "euc" -> "EUC"
  
  feature_tag <- toupper(feature)
  
  # Map distance metric names to standardised abbreviations.
  dist_tag <- switch(
    tolower(distance),
    "correlation" = "CORR", "corr" = "CORR",
    "euclidean" = "EUC", "euc" = "EUC",
    toupper(distance)  # Fallback: use uppercase of input
  )
  
  list(feature_tag = feature_tag, dist_tag = dist_tag)
}

get_method_context <- function(
  cfg, method_id,
  tumour_prefix  = NULL,
  overlay_prefix = NULL
) {
  # get_method_context(): Retrieves complete method configuration from YAML config.
  #
  # This function looks up a method_id in the configuration file and returns
  # a comprehensive context object containing all information needed for
  # file naming and method-aware processing.
  #
  # The function validates that the specified method_id exists in the
  # configuration, providing early failure with a clear error message if
  # the method is not defined.
  #
  # Parameters:
  #   cfg:            The parsed YAML configuration object
  #   method_id:      The method identifier to look up (e.g., "pam50_euc")
  #   tumour_prefix:  Prefix for tumour-only output files (optional; falls back to config or generic)
  #   overlay_prefix: Prefix for joint tumour-cell line files (optional; falls back to config or generic)
  #
  # Returns:
  #   A list containing:
  #     - method_id: The original method identifier
  #     - feature: The feature set name (e.g., "PAM50", "HVG")
  #     - distance: The distance metric (e.g., "correlation", "euclidean")
  #     - feature_tag: Uppercase feature tag for filenames
  #     - dist_tag: Abbreviated distance tag for filenames
  #     - cluster_prefix: Full prefix for tumour clustering outputs
  #     - overlay_prefix: Full prefix for joint analysis outputs
  #
  # FILENAME CONVENTION:
  # Output files follow the pattern:
  #   {prefix}_{FEATURE}_{DIST}_*.{ext}
  # For example:
  #   TUMOUR_PAM50_CORR_clusters.rds
  #   TUMOUR-DSMZ_HVG_EUC_umap.pdf
  #
  # The prefix is cohort-independent and can be overridden via config or function arguments.
  
  # Validate that the method_id exists in configuration.
  if (is.null(cfg$methods[[method_id]])) {
    stop(sprintf("Unknown method_id '%s' – check config.yaml$methods", method_id))
  }
  
  # --------------------------------------------------------------------------
  # Cohort-independent filename prefixes
  # --------------------------------------------------------------------------
  # If prefixes are not supplied, fall back to config or generic defaults.
  # This keeps filenames informative without hard-coding BRCA/TCGA.
  tumour_prefix  <- tumour_prefix  %||% cfg$labels$tumour_prefix  %||% "TUMOUR"
  overlay_prefix <- overlay_prefix %||% cfg$labels$overlay_prefix %||% "TUMOUR-DSMZ"
  
  # Extract method configuration.
  m_cfg <- cfg$methods[[method_id]]
  feature <- m_cfg$feature
  distance <- m_cfg$distance
  
  # Generate filename tags.
  tags <- make_method_tags(feature, distance)
  
  # Construct full filename prefixes.
  cluster_prefix <- sprintf("%s_%s_%s", tumour_prefix, tags$feature_tag, tags$dist_tag)
  overlay_prefix_full <- sprintf("%s_%s_%s", overlay_prefix, tags$feature_tag, tags$dist_tag)
  
  list(
    method_id = method_id,
    feature = feature,
    distance = distance,
    feature_tag = tags$feature_tag,
    dist_tag = tags$dist_tag,
    cluster_prefix = cluster_prefix,
    overlay_prefix = overlay_prefix_full
  )
}

# ==============================================================================
# SECTION 3: SAMPLE IDENTIFIER UTILITIES
# ==============================================================================
# Cell line sample identifiers often contain technical prefixes and suffixes
# that clutter visualisations. This function extracts the core cell line
# name for display purposes while preserving the full identifier for
# data integrity.

shorten_cell_line_name <- function(x) {
  # shorten_cell_line_name(): Extracts display-friendly cell line names.
  #
  # DSMZ cell line identifiers may include sequencing batch prefixes
  # (e.g., "NG-12345_") and library suffixes (e.g., "_lib001_002") that
  # are essential for data tracking but obscure the biological identity
  # in visualisations.
  #
  # This function removes these technical components to produce clean
  # labels suitable for plot annotations.
  #
  # Parameters:
  #   x: A character vector of cell line identifiers
  #
  # Returns:
  #   A character vector with technical prefixes and suffixes removed
  #
  # Example:
  #   "NG-12345_MCF7_lib001_002" -> "MCF7"
  
  x %>%
    gsub("^NG-[0-9]+_", "", .) %>%   # Remove sequencing batch prefix
    gsub("_lib[0-9_]+$", "", .)      # Remove library suffix
}

# ==============================================================================
# SECTION 4: COLOUR PALETTES
# ==============================================================================
# Consistent colour coding is essential for interpretability across multiple
# visualisations. The pipeline defines two colour palettes:
#
#   1. PAM50 colours: For breast cancer molecular subtypes
#   2. Cluster colours: For numerically-labelled unsupervised clusters

# ------------------------------------------------------------------------------
# PAM50 MOLECULAR SUBTYPE COLOURS
# ------------------------------------------------------------------------------
# Note: PAM50 palette is only used when plotting PAM50 subtype annotations.
# It does not imply the pipeline is BRCA-only; this palette is available
# for any analysis that uses PAM50 subtype classification.
#
# The PAM50 classifier identifies five intrinsic molecular subtypes of breast
# cancer, each with distinct biological characteristics and clinical implications:
#
#   Basal-like (green): Triple-negative phenotype, high proliferation,
#   aggressive behaviour, poor prognosis but may respond to chemotherapy.
#
#   HER2-enriched (orange): HER2 amplification, targetable with anti-HER2
#   therapies (trastuzumab, pertuzumab).
#
#   Luminal A (blue): ER-positive, low proliferation, best prognosis,
#   typically treated with endocrine therapy alone.
#
#   Luminal B (red): ER-positive, higher proliferation than Luminal A,
#   may benefit from chemotherapy in addition to endocrine therapy.
#
#   Normal-like (purple): Expression profile similar to normal breast tissue,
#   relatively rare, prognosis intermediate between Luminal A and Luminal B.
#
# These colours follow conventions established in the breast cancer literature
# to facilitate immediate recognition by researchers familiar with PAM50.

pam50_colours <- c(
  "Basal"  = "#4DAF4A",
  "Her2"   = "#FF7F00",
  "LumA"   = "#377EB8",
  "LumB"   = "#E41A1C",
  "Normal" = "#984EA3"
)

# ------------------------------------------------------------------------------
# UNSUPERVISED CLUSTER COLOURS
# ------------------------------------------------------------------------------
# For clusters identified through unsupervised analysis (where biological
# identity is not predetermined), a palette of visually distinct colours
# is provided. The palette supports up to 8 clusters; additional clusters
# will be assigned colours from the default ggplot2 hue scale.
#
# The colours are selected to maximise perceptual distinguishability,
# following principles of colour theory for categorical data visualisation.

cluster_colours <- c(
  "1" = "#E41A1C",
  "2" = "#377EB8",
  "3" = "#4DAF4A",
  "4" = "#FF7F00",
  "5" = "#984EA3",
  "6" = "#FFFF33",
  "7" = "#A65628",
  "8" = "#F781BF"
)

# ==============================================================================
# SECTION 5: VISUALISATION THEME
# ==============================================================================
# A consistent visual theme across all pipeline outputs facilitates comparison
# and interpretation. This custom ggplot2 theme is optimised for dimensionality
# reduction plots (PCA, UMAP) with the following design principles:
#
#   - Clean white background to maximise data-ink ratio
#   - Minimal gridlines to avoid visual clutter
#   - Clear axis labels and titles with appropriate hierarchy
#   - Consistent legend positioning and formatting
#   - Adequate margins for figure integration

theme_dimred <- function(base_size = 11) {
  # theme_dimred(): Custom ggplot2 theme for dimensionality reduction plots.
  #
  # This theme extends theme_classic() with modifications specifically
  # suited for scatter plots of PCA and UMAP coordinates. Key features:
  #
  #   - White background: Provides maximum contrast for coloured points
  #   - No panel grid: Coordinates in reduced space are relative, so
  #     gridlines add visual noise without aiding interpretation
  #   - Bold title: Establishes clear visual hierarchy
  #   - Right-aligned legend: Preserves square aspect ratio of plot area
  #   - Subtle axis styling: Present but not distracting
  #
  # Parameters:
  #   base_size: Base font size in points (default: 11)
  #
  # Returns:
  #   A ggplot2 theme object that can be added to any ggplot
  #
  # Usage:
  #   ggplot(data, aes(x, y)) + geom_point() + theme_dimred(base_size = 12)
  
  theme_classic(base_size = base_size) +
    theme(
      # Background: Clean white canvas
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      
      # Title: Bold, dark grey, left-aligned
      plot.title = element_text(
        size = base_size + 3,
        face = "bold",
        colour = "grey10",
        hjust = 0,
        margin = margin(b = 5)
      ),
      
      # Subtitle: Lighter grey, provides methodological context
      plot.subtitle = element_text(
        size = base_size,
        colour = "grey40",
        hjust = 0,
        margin = margin(b = 10)
      ),
      
      # Axes: Subtle but present
      axis.line = element_line(linewidth = 0.4, colour = "grey30"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text = element_text(size = base_size - 1, colour = "grey20"),
      axis.title = element_text(size = base_size, colour = "grey10"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      
      # Legend: Right-aligned, boxed for clarity
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_rect(
        fill = "white",
        colour = "grey80",
        linewidth = 0.3
      ),
      legend.margin = margin(6, 8, 6, 8),
      legend.title = element_text(
        size = base_size,
        face = "bold",
        colour = "grey20"
      ),
      legend.text = element_text(
        size = base_size - 1,
        colour = "grey30"
      ),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      
      # Margins: Adequate space for figure integration
      plot.margin = margin(15, 15, 10, 10)
    )
}

# ==============================================================================
# SECTION 6: PCA VISUALISATION
# ==============================================================================
# Principal Component Analysis (PCA) is a linear dimensionality reduction
# technique that identifies orthogonal axes of maximum variance in the data.
# PCA is particularly useful for:
#
#   - Visualising global structure in high-dimensional expression data
#   - Identifying batch effects or technical artefacts
#   - Assessing whether cell lines occupy similar regions of expression
#     space as primary tumours
#
# The first two principal components (PC1, PC2) typically capture the
# dominant biological variation, though the proportion of variance explained
# should always be reported to contextualise the visualisation.

save_pca_by_cluster <- function(mat,
                                meta_df,
                                out_pdf,
                                title = "PCA: samples",
                                subtitle = NULL,
                                label_dsmz = FALSE) {
  # save_pca_by_cluster(): Generates PCA scatter plot with cluster colouring.
  #
  # This function performs PCA on an expression matrix and generates a
  # scatter plot of the first two principal components. Samples are coloured
  # by cluster assignment, with optional overlay of DSMZ cell lines.
  #
  # BIOLOGICAL CONTEXT:
  # PCA reveals the dominant axes of transcriptomic variation. When cell
  # lines are overlaid onto tumour-derived PCA coordinates, their position
  # indicates which aspects of tumour biology they capture:
  #
  #   - Cell lines near tumour clusters: Likely represent that molecular subtype
  #   - Cell lines in sparse regions: May have diverged from primary tumour biology
  #   - Cell lines forming separate clusters: Potentially culture-adapted
  #
  # Parameters:
  #   mat:        Expression matrix (genes x samples)
  #   meta_df:    Sample metadata with columns: sample, dataset_type, cluster
  #   out_pdf:    Output file path for PDF
  #   title:      Plot title (character)
  #   subtitle:   Plot subtitle; if NULL, variance explained is shown
  #   label_dsmz: If TRUE, cell line points are labelled with shortened names
  #
  # Returns:
  #   The ggplot object (invisibly); also saves PDF and PNG files
  #
  # OUTPUT FILES:
  #   - {out_pdf}: Vector PDF for scalable graphics
  #   - {out_pdf%.pdf}.png: Raster PNG for presentations
  
  # Validate that all matrix columns have corresponding metadata.
  stopifnot(all(colnames(mat) %in% meta_df$sample))
  
  # Ensure matrix column order matches metadata row order.
  mat <- mat[, meta_df$sample, drop = FALSE]
  
  # Perform PCA on samples (columns become rows for prcomp).
  # Scaling is not applied because VST-transformed data already has
  # stabilised variance across the expression range.
  pca <- prcomp(t(mat), scale. = FALSE)
  
  # Calculate proportion of variance explained by each component.
  # This contextualises the 2D visualisation within the full data structure.
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  # Extract PC scores for plotting.
  pc_df <- as.data.frame(pca$x[, 1:2])
  colnames(pc_df) <- c("PC1", "PC2")
  pc_df$sample <- rownames(pca$x)
  
  # Merge with metadata for aesthetic mappings.
  plot_df <- pc_df %>%
    left_join(meta_df, by = "sample") %>%
    mutate(
      dataset_type = normalise_dataset_type(dataset_type %||% "Tumour")
    )
  
  # Generate default subtitle showing variance explained.
  if (is.null(subtitle)) {
    subtitle <- sprintf(
      "PC1: %.1f%%, PC2: %.1f%% variance explained",
      100 * var_expl[1], 100 * var_expl[2]
    )
  }
  
  # Separate tumour and cell line samples for layered plotting.
  # Tumours are plotted first (background), cell lines overlaid (foreground).
  tumour_df <- plot_df %>% filter(dataset_type == "Tumour")
  dsmz_df   <- plot_df %>% filter(dataset_type == "Cell line")
  
  # Safety check: ensure tumour samples are present.
  if (nrow(tumour_df) == 0) {
    stop("No tumour samples detected after dataset_type normalisation. Check meta_df$dataset_type values.")
  }
  
  # Build the plot with tumour samples as the base layer.
  p <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = PC1, y = PC2, colour = cluster, shape = dataset_type),
      size = 1.8, alpha = 0.8
    )
  
  # Add cell line overlay if present.
  if (nrow(dsmz_df) > 0) {
    p <- p +
      geom_point(
        data = dsmz_df,
        aes(x = PC1, y = PC2, shape = dataset_type),
        size = 2.0, stroke = 0.7, colour = "black"
      )
    
    # Add labels for cell lines if requested.
    if (label_dsmz) {
      dsmz_df <- dsmz_df %>%
        mutate(cell_line_short = shorten_cell_line_name(sample))
      
      p <- p +
        ggrepel::geom_label_repel(
          data = dsmz_df,
          aes(x = PC1, y = PC2, label = cell_line_short),
          size = 2.4,
          fill = "white",
          label.size = 0.1,
          label.r = unit(0.1, "lines"),
          colour = "black",
          max.overlaps = Inf,
          segment.size = 0.2,
          segment.alpha = 0.6
        )
    }
  }
  
  # Determine colours for present clusters.
  # Only clusters actually in the data receive colour assignments.
  present_clusters <- as.character(sort(unique(as.integer(as.character(tumour_df$cluster)))))
  plot_colours <- cluster_colours[present_clusters]
  
  # Generate additional colours if more than 8 clusters are present.
  if (any(is.na(plot_colours))) {
    extra <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra
  }
  
  # Apply scales and theme.
  p <- p +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name   = "Sample type"
    ) +
    theme_dimred(base_size = 11) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "PC1",
      y        = "PC2"
    )
  
  # Save outputs in both vector (PDF) and raster (PNG) formats.
  ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
  ggsave(sub("\\.pdf$", ".png", out_pdf), p,
         width = 7, height = 5.5, units = "in", dpi = 300)
  
  invisible(p)
}

# ==============================================================================
# SECTION 7: UMAP VISUALISATION
# ==============================================================================
# Uniform Manifold Approximation and Projection (UMAP) is a non-linear
# dimensionality reduction technique that preserves both local and global
# structure in high-dimensional data. Compared to t-SNE, UMAP:
#
#   - Better preserves global relationships between clusters
#   - Is more computationally efficient
#   - Produces more reproducible embeddings (with fixed random seed)
#
# UMAP is particularly valuable for visualising cancer transcriptomics because
# molecular subtypes often form discrete clusters in expression space, and
# UMAP's neighbourhood-preserving properties highlight these groupings.
#
# KEY PARAMETERS:
#   n_neighbors: Controls the balance between local and global structure.
#   Higher values emphasise global topology; lower values capture finer
#   local detail. Default of 20 is suitable for datasets with clear clusters.
#
#   min_dist: Controls point packing in the embedding. Lower values create
#   tighter clusters; higher values spread points more evenly. Default of
#   0.3 provides clear cluster separation without excessive compression.
#
#   metric: Distance metric for neighbour calculations. Cosine distance
#   (default) is appropriate for gene expression data where relative
#   expression patterns are more informative than absolute levels.

save_umap_by_cluster <- function(mat,
                                 meta_df,
                                 out_pdf,
                                 title = "UMAP: samples",
                                 subtitle = "Coloured by cluster",
                                 n_neighbors = 20,
                                 min_dist = 0.3,
                                 metric = "cosine",
                                 label_dsmz = FALSE) {
  # save_umap_by_cluster(): Generates UMAP embedding with cluster colouring.
  #
  # This function computes a UMAP embedding of the expression data and
  # generates a scatter plot with samples coloured by cluster assignment.
  # Cell lines can be overlaid to visualise their position relative to
  # tumour-defined structure.
  #
  # BIOLOGICAL INTERPRETATION:
  # UMAP embeddings of cancer transcriptomes typically reveal:
  #   - Distinct clusters corresponding to molecular subtypes
  #   - Gradients reflecting continuous biological processes (e.g., EMT)
  #   - Outliers representing unusual biology or technical issues
  #
  # Cell line positions in the embedding indicate their transcriptomic
  # relationship to primary tumours. Cell lines embedded within tumour
  # clusters are likely to recapitulate that subtype's biology.
  #
  # Parameters:
  #   mat:         Expression matrix (genes x samples)
  #   meta_df:     Sample metadata with columns: sample, dataset_type, cluster
  #   out_pdf:     Output file path for PDF
  #   title:       Plot title (character)
  #   subtitle:    Plot subtitle (character)
  #   n_neighbors: UMAP n_neighbors parameter (default: 20)
  #   min_dist:    UMAP min_dist parameter (default: 0.3)
  #   metric:      Distance metric for UMAP (default: "cosine")
  #   label_dsmz:  If TRUE, cell line points are labelled
  #
  # Returns:
  #   The ggplot object (invisibly); also saves PDF and PNG files
  
  # Validate metadata coverage.
  stopifnot(all(colnames(mat) %in% meta_df$sample))
  
  # Ensure consistent sample ordering.
  mat <- mat[, meta_df$sample, drop = FALSE]
  
  # Compute UMAP embedding with fixed seed for reproducibility.
  # The seed ensures identical embeddings across pipeline runs,
  # which is essential for comparing results across methods.
  set.seed(42)
  um <- uwot::umap(
    t(mat),
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = metric
  )
  
  # Format embedding results.
  colnames(um) <- c("UMAP1", "UMAP2")
  rownames(um) <- colnames(mat)
  
  # Merge with metadata for aesthetic mappings.
  um_df <- as.data.frame(um) %>%
    mutate(sample = rownames(um)) %>%
    left_join(meta_df, by = "sample") %>%
    mutate(
      dataset_type = normalise_dataset_type(dataset_type %||% "Tumour")
    )
  
  # Separate sample types for layered plotting.
  tumour_df <- um_df %>% filter(dataset_type == "Tumour")
  dsmz_df   <- um_df %>% filter(dataset_type == "Cell line")
  
  # Safety check: ensure tumour samples are present.
  if (nrow(tumour_df) == 0) {
    stop("No tumour samples detected after dataset_type normalisation. Check meta_df$dataset_type values.")
  }
  
  # Determine colours for present clusters.
  present_clusters <- as.character(sort(unique(as.integer(as.character(tumour_df$cluster)))))
  plot_colours <- cluster_colours[present_clusters]
  if (any(is.na(plot_colours))) {
    extra <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra
  }
  
  # Build the plot with tumours as the base layer.
  p <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = UMAP1, y = UMAP2, colour = cluster, shape = dataset_type),
      size = 1.8,
      alpha = 0.75
    )
  
  # Add cell line overlay if present.
  if (nrow(dsmz_df) > 0) {
    p <- p +
      geom_point(
        data = dsmz_df,
        aes(x = UMAP1, y = UMAP2, shape = dataset_type),
        size = 2,
        stroke = 0.7,
        colour = "black"
      )
    
    # Add labels for cell lines if requested.
    if (label_dsmz) {
      dsmz_df <- dsmz_df %>%
        mutate(cell_line_short = shorten_cell_line_name(sample))
      
      p <- p +
        ggrepel::geom_label_repel(
          data = dsmz_df,
          aes(x = UMAP1, y = UMAP2, label = cell_line_short),
          size = 2.4,
          fill = "white",
          label.size = 0.1,
          label.r = unit(0.1, "lines"),
          colour = "black",
          max.overlaps = Inf,
          segment.size = 0.2,
          segment.alpha = 0.6
        )
    }
  }
  
  # Apply scales, theme, and annotations.
  p <- p +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name   = "Sample type"
    ) +
    coord_fixed(ratio = 1) +  # Preserve aspect ratio for accurate distance perception
    theme_dimred(base_size = 11) +
    guides(
      colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
      shape  = guide_legend(override.aes = list(size = 3))
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "UMAP 1",
      y        = "UMAP 2"
    )
  
  # Save outputs in both vector and raster formats.
  ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
  ggsave(sub("\\.pdf$", ".png", out_pdf), p,
         width = 7, height = 5.5, units = "in", dpi = 300)
  
  invisible(p)
}