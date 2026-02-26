#!/usr/bin/env Rscript

# plot_pan_cancer_umap_deg_set.R
# ==============================================================================
# PAN-CANCER UMAP VISUALISATION USING DIFFERENTIALLY EXPRESSED GENE SET
# ==============================================================================
#
# SCRIPT OVERVIEW
# ---------------
# This script generates UMAP (Uniform Manifold Approximation and Projection)
# embeddings of pan-cancer transcriptomic data using a curated set of
# differentially expressed genes (DEGs). The visualisations enable assessment
# of tumour-cell line alignment across multiple cancer types simultaneously.
#
# BIOLOGICAL CONTEXT
# ------------------
# Pan-cancer analysis examines molecular patterns that transcend individual
# cancer types, revealing both shared oncogenic programmes and lineage-specific
# biology. By integrating data from multiple diseases (breast cancer,
# neuroblastoma, retinoblastoma, haematological malignancies), this analysis
# addresses fundamental questions:
#
#   1. Do cell lines cluster with tumours of their annotated cancer type?
#   2. Are there cross-cancer molecular similarities suggesting shared biology?
#   3. Which cell lines may have diverged from their tissue of origin?
#
# THE DEG FEATURE SET
# -------------------
# Rather than using all expressed genes or generic highly variable genes,
# this script employs a curated set of differentially expressed genes
# identified through the pan-cancer validation pipeline. These genes were
# selected because they:
#
#   - Show significant expression differences between cancer types
#   - Discriminate between tumours and cell lines within cancer types
#   - Are robustly expressed across the integrated dataset
#
# Using this focused gene set (typically 200-300 genes) reduces noise from
# uninformative genes while preserving biologically meaningful variation.
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
# Fill = cancer type (disease), outline = sample type (tumour white, cell line black).
#
#   - Fill: Cancer type (disease origin)
#   - Outline: Tumour = white, cell line = black (pch 21)
#
# Cell lines (black outline) remain visible against tumour clusters.
#
# OUTPUT FILES
# ------------
# For each distance metric, the script generates:
#   - Coordinate table (TSV): UMAP coordinates with sample metadata
#   - Visualisation (PDF/PNG): Pan-cancer overview plot
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
  
  # tibble: Modern reimplementation of data frames with stricter subsetting
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
              default="results/unsupervised/pan_cancer/pan_cancer_expr.rds",
              help="Joint expression matrix RDS [default: %default]"),

  # --meta_tsv: Sample metadata including cancer type and sample type
  # annotations required for visualisation aesthetics.
  make_option("--meta_tsv", type="character",
              default="results/unsupervised/pan_cancer/inputs/joint_metadata.tsv",
              help="Joint metadata TSV [default: %default]"),

  # --deg_set: Path to the DEG feature set file containing one Ensembl
  # gene identifier per line. These genes define the feature space for
  # UMAP computation.
  make_option("--deg_set", type="character",
              default="results/unsupervised/pan_cancer/pan_cancer_features_clean.txt",
              help="DEG feature set (one ENSG per line) [default: %default]"),

  # --outdir: Output directory for coordinates, plots, and summary files.
  make_option("--outdir", type="character",
              default="results/unsupervised/pan_cancer/umap_deg_set",
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
  # 'slide' optimises for 16:9 PowerPoint slides; 'a4' variants for print.
  make_option("--page", type="character", default="slide",
              help="Output size preset: slide (16:9), a4, a4_landscape [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

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

# Resolve all input and output paths.
expr_rds <- abs_from_root(opt$expr_rds)
meta_tsv <- abs_from_root(opt$meta_tsv)
deg_set  <- abs_from_root(opt$deg_set)
outdir   <- abs_from_root(opt$outdir)

# Create output directory.
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Parse distance metrics from comma-separated string.
dist_metrics <- strsplit(opt$dist_metrics, ",")[[1]] %>% trimws()
dist_metrics <- dist_metrics[nzchar(dist_metrics)]

# Log configuration for reproducibility documentation.
cat("[INFO] expr_rds:", expr_rds, "\n")
cat("[INFO] meta_tsv:", meta_tsv, "\n")
cat("[INFO] deg_set :", deg_set, "\n")
cat("[INFO] outdir  :", outdir, "\n")
cat("[INFO] metrics :", paste(dist_metrics, collapse=", "), "\n")
cat(sprintf("[INFO] threads: %d\n", n_threads))

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
  
  if (!file.exists(path)) stop("deg_set not found: ", path)
  
  g <- readLines(path, warn = FALSE)
  g <- g[!grepl("^\\s*#", g)]  # Remove comment lines
  g <- trimws(g)                # Strip whitespace
  g <- g[nzchar(g)]             # Remove empty lines
  g <- clean_ensg(g)            # Remove version suffixes
  g <- unique(g)                # Deduplicate
  
  if (length(g) == 0) stop("deg_set empty after cleaning: ", path)
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
#   - Bottom legend position to preserve square plot aspect ratio
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
      axis.text         = element_text(size = 9, colour = "grey20"),
      axis.title        = element_text(size = 10, colour = "grey10"),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9),
      legend.key.size   = grid::unit(0.8, "lines"),
      legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
      plot.title        = element_text(size = 13, face = "bold", colour = "grey10"),
      plot.subtitle     = element_text(size = 10, colour = "grey40"),
      plot.caption      = element_text(size = 7, colour = "grey50", hjust = 0),
      plot.margin       = margin(12, 12, 8, 8),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

# ==============================================================================
# SECTION 6: DATA PREPARATION FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
prepare_plot_data <- function(df) {
  df %>%
    mutate(
      type = ifelse(sample_type == "Cell Line", "CL", "tumour"),
      type = factor(type, levels = c("tumour", "CL")),
      disease = case_when(
        cancer_type == "BRCA" ~ "Breast Cancer",
        cancer_type == "NBL"  ~ "Neuroblastoma",
        cancer_type == "RBL"  ~ "Retinoblastoma",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(disease)) %>%
    mutate(
      disease = factor(disease, levels = c("Breast Cancer","Neuroblastoma","Retinoblastoma"))
    )
}

# ==============================================================================
# SECTION 7: PLOTTING FUNCTION
# ==============================================================================
# Fill = disease; outline = sample type (tumour = white, cell line = black).

pan_cancer_alignment_plot <- function(alignment, palette, subtitle_text,
                                      caption_text = NULL, show_legend = TRUE,
                                      point_alpha = 0.7,
                                      stroke_tumour = 0.35, stroke_cl = 0.55,
                                      size_tumour = 0.85, size_cl = 1.15) {

  df <- alignment %>% mutate(type = factor(type, levels = c("tumour", "CL")))

  gg <- ggplot(df, aes(UMAP1, UMAP2, fill = disease)) +

    geom_point(
      data = df %>% filter(type == "tumour"),
      aes(colour = type),
      pch = 21, alpha = point_alpha, size = size_tumour, stroke = stroke_tumour
    ) +

    geom_point(
      data = df %>% filter(type == "CL"),
      aes(colour = type),
      pch = 21, alpha = point_alpha, size = size_cl, stroke = stroke_cl
    ) +

    scale_colour_manual(values = c(`CL` = "black", `tumour` = "white"), guide = "none") +
    scale_fill_manual(values = palette) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    labs(subtitle = subtitle_text, caption = caption_text)

  if (!show_legend) gg <- gg + guides(fill = "none")
  else gg <- gg + guides(fill = guide_legend(title = "Cancer type"))

  gg
}

# ==============================================================================
# SECTION 8: OUTPUT SAVING FUNCTION
# ==============================================================================

save_plot <- function(plot, path_base, page = "slide") {
  # save_plot(): Saves plot in both PDF and PNG formats.
  #
  # The function supports multiple size presets optimised for different
  # presentation contexts:
  #
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

  if (page == "slide") {
    width <- 13.33
    height <- 7.5
  } else if (page == "a4") {
    width <- 8.27
    height <- 11.69
  } else if (page == "a4_landscape") {
    width <- 11.69
    height <- 8.27
  } else {
    stop("Unknown --page: ", page, " (use slide, a4, a4_landscape)")
  }

  # Save vector (PDF) and raster (PNG) versions.
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, 
         units = "in", device = grDevices::cairo_pdf)
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, 
         units = "in", dpi = 400)
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
# Standardise sample type labels to consistent values for filtering and plotting.

meta <- meta %>%
  mutate(
    sample_type = case_when(
      sample_type %in% c("tumour","Tumour","TUMOUR") ~ "Tumour",
      sample_type %in% c("cell_line","Cell Line","CELL_LINE","cellline") ~ "Cell Line",
      TRUE ~ sample_type
    )
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

# Diagnostic output: verify metadata structure.
cat("\n[CHECK] Raw meta counts (before filtering to X):\n")
print(table(meta$cancer_type, meta$sample_type))
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
print(meta %>% dplyr::filter(sample_type == "Cell Line") %>% 
        dplyr::count(cancer_type, sort=TRUE))

cat("\n[CHECK] Cell lines with technical replicates (collapsed by _lib suffix):\n")
meta %>%
  filter(sample_type == "Cell Line", cancer_type %in% c("NBL","RBL")) %>%
  mutate(cell_line_base = stringr::str_replace(sample_id, "_lib.*$", "")) %>%
  group_by(cancer_type) %>%
  summarise(
    n_samples = n(),
    n_unique_cell_lines = n_distinct(cell_line_base),
    .groups = "drop"
  ) %>%
  print()

# ==============================================================================
# SECTION 10: DEG FEATURE SET LOADING AND VALIDATION
# ==============================================================================
# Load the differentially expressed gene set and validate that all genes
# are present in the expression matrix.

genes_fixed <- read_feature_list(deg_set)
genes_present <- intersect(colnames(X), genes_fixed)

cat("[INFO] DEG set genes:", length(genes_fixed), "\n")
cat("[INFO] Present in expr:", length(genes_present), "\n")

# Report missing genes (indicates annotation version mismatch).
missing <- setdiff(genes_fixed, colnames(X))
if (length(missing) > 0) {
  cat("[WARN] DEG genes missing from expr (showing up to 20):\n")
  print(head(missing, 20))
}

# Require complete gene coverage for reproducibility.
if (length(genes_present) != length(genes_fixed)) {
  stop("Not all DEG genes are present in expr. Present=", length(genes_present),
       " Expected=", length(genes_fixed),
       ". Fix upstream gene harmonization or DEG IDs.")
}

# Subset expression matrix to DEG genes only.
X_sub <- X[, genes_present, drop = FALSE]

cat("[INFO] X_sub:", nrow(X_sub), "samples x", ncol(X_sub), "genes\n")
cat("[INFO] Sample summary:\n")
print(table(meta$cancer_type, meta$sample_type))

# Restrict to BRCA/NBL/RBL for UMAP and plotting (match pan_cancer_umap_main.R).
keep <- meta$cancer_type %in% c("BRCA", "NBL", "RBL")
meta <- meta[keep, ]
X_sub <- X_sub[meta$sample_id, , drop = FALSE]
cat("[INFO] Restricted to BRCA/NBL/RBL:", nrow(X_sub), "samples\n")

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
  tag <- paste0("DEG_SET_", dm)
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
  coords <- coords %>% left_join(meta, by = "sample_id")
  
  # Save coordinate table.
  coords_path <- file.path(outdir, paste0("coords_", tag, ".tsv"))
  write_tsv(coords, coords_path)

  # Generate visualisation.
  plot_df <- prepare_plot_data(coords)
  subtitle <- sprintf("DEG set (n=%d genes); metric=%s; nn=%d; min_dist=%.2f",
                      ncol(X_sub), dm, opt$n_neighbors, opt$min_dist)

  p <- pan_cancer_alignment_plot(
    alignment     = plot_df,
    palette       = disease_palette,
    subtitle_text = subtitle,
    caption_text  = NULL,
    show_legend   = TRUE
  )
  plot_base <- file.path(outdir, paste0("Fig_umap_", tag))
  save_plot(p, plot_base, page = opt$page)

  # Record summary for this metric.
  summaries[[dm]] <- tibble(
    feature_set = "DEG_SET",
    dist_metric = dm,
    n_genes     = ncol(X_sub),
    n_samples   = nrow(X_sub),
    n_tumours   = sum(meta$sample_type == "Tumour"),
    n_cell_lines= sum(meta$sample_type == "Cell Line"),
    coords_tsv  = basename(coords_path),
    plot_pdf    = basename(paste0(plot_base, ".pdf")),
    plot_png    = basename(paste0(plot_base, ".png"))
  )

  cat("[SAVED] ", coords_path, "\n", sep="")
  cat("[SAVED] ", paste0(plot_base, ".pdf/.png"), "\n", sep="")
}

# ==============================================================================
# SECTION 13: SUMMARY OUTPUT
# ==============================================================================
# Compile and save a summary table documenting all analyses performed.

summary_tbl <- bind_rows(summaries)
summary_path <- file.path(outdir, "summary_umap_DEG_SET.tsv")
write_tsv(summary_tbl, summary_path)

cat("\n[SAVED] ", summary_path, "\n", sep="")
cat("\n[SUCCESS] UMAP plotting complete.\n")