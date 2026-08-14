#!/usr/bin/env Rscript
# =============================================================================
# feature_selection_unsupervised.R
# =============================================================================
#
# PURPOSE:
# This script performs unsupervised, multi-method feature selection from a
# joint VST-normalised expression matrix. The selected genes define alternative
# reduced expression spaces for downstream clustering, tumour-neighbourhood
# analysis, similarity estimation and consensus graph construction.
#
# The script does not treat any single statistic as the only definition of an
# informative gene. Instead, it generates multiple method-specific gene lists
# that capture complementary properties of the expression matrix.
#
# APPROACH:
# The workflow has two feature-selection tiers:
#
#   1. First-pass selection criteria
#      These methods each retain a broad top-N gene list, usually 3,000 genes.
#      They score genes using dispersion, mean-variance adjusted variability,
#      distributional complexity, or contribution to leading multivariate axes:
#
#        - Variance
#        - HVG / mean-variance trend-adjusted variance
#        - Median Absolute Deviation (MAD)
#        - Mean Absolute Deviation (MeanAbsDev)
#        - Shannon entropy
#        - PCA loading score
#
#   2. Connectivity-based criteria
#      These methods operate on the union of the first-pass selected genes.
#      They retain smaller feature sets, usually 500 genes, by scoring genes
#      according to co-expression centrality and, for MX, centrality combined
#      with scaled variance:
#
#        - Spearman mean absolute connectivity
#        - MX score: Spearman connectivity × scaled variance
#        - WGCNA kTotal
#
# Each selected gene list is written separately and used downstream as a
# method-specific representation. This allows later consensus steps to favour
# tumour-cell-line relationships that remain stable across several definitions
# of gene informativeness.
#
# INPUTS:
#   --config   : Path to config.yaml containing paths and parameters
#   --vst_rds  : Path to VST-normalised expression matrix (.rds file)
#
# OUTPUTS:
#   - BRCA_MX_genes_top500.txt : Final selected MX gene list
#   - BRCA_MX_joint_ranks_kTotal_vs_MX.tsv : Full ranking table
#   - UpSet plots comparing method overlap
#   - Various intermediate files for QC and method comparison
#
# USAGE:
#   Rscript feature_selection_unsupervised.R \
#     --config config/config.yaml \
#     --vst_rds results/vst_joint.rds
#
# =============================================================================

# =============================================================================
# LOAD REQUIRED PACKAGES
# =============================================================================
# Note: All packages should be pre-installed in the conda/mamba environment
# (e.g., fs_wgcna). Packages are not installed within this script to ensure
# reproducibility and avoid permission issues on HPC systems.
#
# The suppressPackageStartupMessages() function silences the verbose loading
# messages that packages such as dplyr and ggplot2 typically print, keeping
# the console output clean and focused on the script's progress messages.

suppressPackageStartupMessages({
  library(optparse)      # Provides command-line argument parsing functionality
  library(yaml)          # Enables reading of YAML configuration files
  library(matrixStats)   # Provides fast row/column statistics (rowVars, rowSds, rowRanks)
                         # These functions are optimised in C and much faster than apply()
  library(stats)         # Contains base R statistics functions (cor, prcomp, scale)
  library(dplyr)         # Facilitates data manipulation (tibble, mutate, arrange, filter)
  library(WGCNA)         # Implements Weighted Gene Co-expression Network Analysis
                         # Used for soft-thresholded network construction and hub detection
  library(readr)         # Provides fast file I/O operations (write_tsv, read_tsv)
  library(ggplot2)       # Enables publication-quality data visualisation
  library(UpSetR)        # Creates UpSet plots for visualising set intersections
                         # More effective than Venn diagrams for >3 sets
})

# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================
# The expected command-line arguments are defined using optparse. This approach
# allows the script to be called from Snakemake or the command line with
# flexible parameter specification.

option_list <- list(
  make_option(
    c("--config"),
    type = "character",
    default = NULL,
    help = "Path to config.yaml (optional if other args provided)",
    metavar = "FILE"
  ),
  make_option(
    c("--vst_rds"),
    type = "character",
    help = "Path to joint VST .rds",
    metavar = "FILE"
  ),
  make_option(
    c("--outdir"),
    type = "character",
    default = NULL,
    help = "Output directory for feature selection results",
    metavar = "DIR"
  ),
  make_option(
    c("--scores_tsv"),
    type = "character",
    default = NULL,
    help = "Output path for joint HVG scores table (TSV)"
  ),
  make_option(
    c("--mx_list"),
    type = "character",
    default = NULL,
    help = "Output path for final MX gene list (TXT)"
  ),
  make_option(
    c("--profile"),
    type = "character",
    default = NULL,
    help = "Profile name under config$profiles (e.g. brca/nbl/rbl/heme/pan_cancer)"
  ),
  make_option(
    c("--scripts_dir"),
    type = "character",
    default = ".",
    help = "Base scripts directory (so we can source helper files robustly)"
  ),
  make_option(
    c("--top_n_method"),
    type = "integer",
    default = NULL,
    help = "Number of top genes per method (default: 3000)"
  ),
  make_option(
    c("--final_top"),
    type = "integer",
    default = NULL,
    help = "Final number of genes to select per method (default: 500)"
  ),
  make_option(
    c("--method_topn"),
    type = "character",
    default = NULL,
    help = paste0(
      "Comma-separated method:N pairs specifying the top-N for each method ",
      "(e.g. 'Variance:3000,MAD:3000,MX:500,HVG:3000'). ",
      "Overrides --top_n_method for named methods. ",
      "Used by the Snakefile to pass per-method top-N values."
    )
  )
)

# Arguments are parsed from the command line.
opt <- parse_args(OptionParser(option_list = option_list))

# The script validates that required arguments have been provided.
if (is.null(opt$vst_rds)) {
  stop("--vst_rds must be provided.")
}

# Shared plotting helpers (theme_dimred) are sourced from an external file, if available.
helper_path <- file.path(opt$scripts_dir, "R", "brca_unsup_methods.R")
if (file.exists(helper_path)) {
  source(helper_path)
} else {
  message("[WARN] helper file not found: ", helper_path, " (using basic ggplot2 theme)")
  theme_dimred <- function(base_size = 12) ggplot2::theme_bw(base_size = base_size)
}

# =============================================================================
# LOAD CONFIGURATION
# =============================================================================
# The config.yaml file contains paths and parameters that may vary between
# runs or environments. If not provided, use command-line arguments directly.

cfg <- NULL
profile_cfg <- NULL

if (!is.null(opt$config)) {
  cfg <- yaml::read_yaml(opt$config)
}

if (!is.null(opt$profile) && !is.null(cfg) && !is.null(cfg$profiles[[opt$profile]])) {
  profile_cfg <- cfg$profiles[[opt$profile]]
}

# Output directory: use --outdir if provided, otherwise from profile config, or error
if (!is.null(opt$outdir)) {
  fs_subdir <- opt$outdir
} else if (!is.null(profile_cfg) && !is.null(profile_cfg$paths$unsup_root)) {
  fs_subdir <- file.path(profile_cfg$paths$unsup_root, "feature_selection_unsupervised")
} else {
  stop("Provide --outdir or (--config and --profile with profiles.<profile>.paths.unsup_root).")
}

# The output directory is created if it does not already exist.
if (!dir.exists(fs_subdir)) {
  dir.create(fs_subdir, recursive = TRUE, showWarnings = FALSE)
}

# Configuration details are logged for debugging and reproducibility purposes.
message("[FEATURE_SELECTION] outdir: ", fs_subdir)
message("[FEATURE_SELECTION] vst_rds: ", opt$vst_rds)


# =============================================================================
# MAIN FEATURE SELECTION FUNCTION
# =============================================================================

run_unsupervised_feature_selection <- function(
  vst_rds,
  outdir = fs_subdir,
  top_n_method = 3000,
  final_top = 500,
  method_topn_map = list(),
  seed = 123,
  quiet = FALSE
) {
  # ---------------------------------------------------------------------------
  # FUNCTION: run_unsupervised_feature_selection
  # ---------------------------------------------------------------------------
  # This function implements a two-tier feature-selection workflow:
  #
  #   Tier 1: Six first-pass criteria each select a broad top-N list.
  #           These criteria measure different properties of single-gene
  #           expression behaviour or multivariate contribution.
  #
  #   Tier 2: Three connectivity-based criteria are computed on the union of the
  #           first-pass gene lists. These criteria rank genes by co-expression
  #           centrality, or by centrality combined with scaled variance.
  #
  # The canonical outputs preserve each method-specific gene list separately.
  # Downstream rules treat these lists as alternative representations rather than
  # merging them into one final gene list.
  #
  # PARAMETERS:
  # -----------
  #   vst_rds : character
  #       Path to VST-normalised expression matrix stored as an RDS file.
  #       The matrix should be oriented as genes (rows) × samples (columns)
  #       with gene identifiers (e.g., Ensembl IDs) as rownames.
  #
  #   outdir : character
  #       Output directory for all results. The directory is created if it
  #       does not exist.
  #
  #   top_n_method : integer, default = 3000
#       Number of top genes to select from each individual method.
#       These genes form the candidate pool for network-based refinement.
#       A value of 3000 typically includes ~15-20% of expressed genes.
  #
  #   final_top : integer, default = 500
  #       Final number of genes to select after network-based ranking.
  #       500 genes is a common choice that balances information content
  #       with computational tractability for clustering.
  #
  #   seed : integer, default = 123
  #       Random seed for reproducibility. Although most operations are
  #       deterministic, this ensures consistent results if any stochastic
  #       operations are introduced.
  #
  #   quiet : logical, default = FALSE
  #       If TRUE, progress messages and debug information are suppressed.
  #       Useful for batch processing or when embedding in pipelines.
  #
  # RETURNS:
  # --------
  #   A list (invisible) containing:
  #     - candidates: Character vector of all candidate genes (union of methods)
  #     - mx: Character vector of final top genes by MX score
  #     - kTotal: Character vector of final top genes by WGCNA connectivity
  #     - overlap: Character vector of genes in both mx and kTotal
  #     - scores_df: Data frame with all scores and rankings for each gene
  #     - joint_tsv: Path to the saved ranking table
  #     - genes_file: Path to the saved final gene list
  # ---------------------------------------------------------------------------
  
  # The random seed is set to ensure reproducibility of any stochastic operations.
  # While most methods in this script are deterministic, setting a seed is
  # good practice and ensures consistent results if the code is modified.
  set.seed(seed)
  
  # The output directory is created if it does not exist.
  # recursive = TRUE creates parent directories if needed.
  # showWarnings = FALSE suppresses warnings if the directory already exists.
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  subdir <- outdir
  
  # ===========================================================================
  # SECTION 1: LOAD AND VALIDATE INPUT DATA
  # ===========================================================================
  
  if (!quiet) message("[INFO] Loading VST matrix from: ", vst_rds)
  counts <- readRDS(vst_rds)
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  # The input must be a numeric matrix with genes as rows and samples as columns.
  # Gene identifiers (e.g., ENSG00000141510) must be present as rownames.
  # Sample identifiers should be present as colnames (though optional).
  #
  # VST (Variance Stabilising Transformation) normalisation is assumed to have
  # been applied upstream (e.g., via DESeq2::vst()). VST transforms count data
  # to approximately homoskedastic values, making variance-based methods more
  # reliable.
  
  if (!is.matrix(counts)) stop("VST object must be a matrix")
  if (is.null(rownames(counts))) stop("Matrix must have gene rownames")

  ensembl_fraction <- function(ids) {
    if (is.null(ids) || length(ids) == 0) {
      return(0)
    }
    mean(grepl("^ENS[A-Z]*G[0-9]+", ids))
  }

  row_gene_frac <- ensembl_fraction(rownames(counts))
  col_gene_frac <- ensembl_fraction(colnames(counts))
  if (col_gene_frac > row_gene_frac && col_gene_frac >= 0.5) {
    if (!quiet) {
      message("[WARN] Input matrix appears to be samples x genes; transposing to genes x samples")
    }
    counts <- t(counts)
  }
  
  n_genes <- nrow(counts)
  n_samples <- ncol(counts)
  if (!quiet) message("[INFO] Matrix: ", n_genes, " genes × ", n_samples, " samples")
  
  # ===========================================================================
  # SECTION 2: DATA CLEANING - REMOVE UNINFORMATIVE GENES
  # ===========================================================================
  # Before feature selection, genes that cannot contribute to distinguishing
  # samples are removed. This improves computational efficiency and prevents
  # edge cases in downstream calculations.
  
  # ---------------------------------------------------------------------------
  # STEP 2.1: Remove genes with zero expression across all samples
  # ---------------------------------------------------------------------------
  # Genes with zero counts in every sample provide no information and would
  # cause division-by-zero or undefined variance calculations.
  #
  # In VST-normalised data, true zeros are rare (VST applies a pseudocount),
  # but this check catches any edge cases.
  
  zero_counts <- rowSums(counts == 0) == n_samples
  if (sum(zero_counts) > 0) {
    if (!quiet) message("[INFO] Removing ", sum(zero_counts), " genes with zero counts")
    counts <- counts[!zero_counts, , drop = FALSE]
    n_genes <- nrow(counts)
    if (!quiet) message("[INFO] Matrix after removing zero counts: ",
                        n_genes, " genes × ", n_samples, " samples")
  }
  
  # ---------------------------------------------------------------------------
  # STEP 2.2: Remove genes with near-zero variance (constant expression)
  # ---------------------------------------------------------------------------
  # Genes with constant or near-constant expression across samples cannot
  # distinguish between cell types, conditions, or clusters.
  #
  # A threshold of 1e-8 is used rather than exactly zero to handle
  # floating-point precision issues. Standard deviation is used rather than
  # variance because it is in the same units as the data.
  #
  # matrixStats::rowSds() is used instead of apply(counts, 1, sd) because
  # it is implemented in C and is approximately 10× faster for large matrices.
  
  tz_row_vars <- matrixStats::rowSds(counts, na.rm = TRUE)
  tz_row_vars[is.na(tz_row_vars)] <- 0  # Replace NA with 0 (treated as constant)
  tz_keep_genes <- tz_row_vars > 1e-8
  tz_n_removed <- sum(!tz_keep_genes)
  
  if (tz_n_removed > 0) {
    if (!quiet) message("[INFO] Removing ", tz_n_removed, " constant/zero-variance genes")
    counts <- counts[tz_keep_genes, , drop = FALSE]
    n_genes <- nrow(counts)
    if (!quiet) message("[INFO] Matrix after variance filter: ",
                        n_genes, " genes × ", n_samples, " samples")
  }
  
  # A safety check ensures that genes remain after filtering.
  # This would only fail if the input data is severely problematic.
  if (nrow(counts) == 0) stop("No variable genes left after filtering")
  
  # ===========================================================================
  # HELPER FUNCTIONS
  # ===========================================================================
  
  # ---------------------------------------------------------------------------
  # HELPER: Write gene lists to text files
  # ---------------------------------------------------------------------------
  # This helper function writes a character vector of gene names to a text file,
  # one gene per line. This format is easy to read in R, Python, or shell scripts.
  
  write_vec <- function(vec, filename) {
    writeLines(vec, file.path(subdir, filename))
  }
  
  # ---------------------------------------------------------------------------
  # HELPER: Save UpSet plot with explanatory legend
  # ---------------------------------------------------------------------------
  # UpSet plots are useful for visualising intersections
  # between more than 3 sets. This function creates an UpSet plot with an

  # additional legend panel explaining what each method measures.
  #
  # IMPORTANT: The upset() call must be wrapped in print() for proper rendering
  # on headless systems (HPC clusters, Docker containers) where there is no
  # active graphics device.
  #
  # PARAMETERS:
  #   sets : Named list of character vectors, each containing gene IDs
  #   outfile_png : Path for the output PNG file
  #   title : Optional title for the legend panel
  
  save_upset_with_legend <- function(sets, outfile_png, title = "") {
    # Ensure all sets contain unique elements (duplicates would skew counts)
    sets <- lapply(sets, unique)
    
    # Remove empty sets (UpSetR requires at least 2 non-empty sets)
    sets <- sets[vapply(sets, length, integer(1)) > 0]
    stopifnot(length(sets) >= 2)
    
    # Convert list to UpSetR input format (binary membership matrix)
    upset_input <- UpSetR::fromList(sets)
    
    # Open PNG device with high resolution for publication quality
    # Width/height in pixels, resolution in DPI
    png(outfile_png, width = 2200, height = 1400, res = 180)
    
    # Create a two-panel layout: main UpSet plot (4 parts) + legend (1 part)
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
    par(mar = c(2, 2, 3, 1))  # Set margins for the main plot
    
    # MAIN UPSET PLOT
    # CRITICAL: print() is required for UpSetR to render on headless systems.
    # Without print(), the plot may not appear in the saved file.
    print(
      UpSetR::upset(
        upset_input,
        order.by = "freq",           # Order intersections by frequency (most common first)
        mb.ratio = c(0.65, 0.35),    # Ratio of main bar plot to set size bars
        mainbar.y.label = "Intersection size",
        sets.x.label = "Set size"
      )
    )
    
    # LEGEND PANEL
    # This panel provides context for interpreting the UpSet plot by
    # explaining what each feature selection method measures.
    par(mar = c(2, 2, 3, 2))
    plot.new()
    title(main = "Legend / Notes", cex.main = 1.1)
    
    legend_text <- c(
      "Variance: highest gene expression variance across samples",
      "HVG: variance above mean-variance LOESS trend (trend-corrected)",
      "MAD: median absolute deviation across samples",
      "MeanAbsDev: mean absolute deviation across samples",
      "Entropy: Shannon entropy of gene expression distribution",
      "PCA: highest absolute loadings on PC1-PC5",
      "Spearman: highest mean absolute Spearman correlation (hub genes)",
      "MX: Spearman connectivity × scaled gene variance",
      "kTotal: WGCNA total network connectivity (hubness)"
    )
    
    # Position text in the legend panel
    text(x = 0, y = 1, adj = c(0, 1), labels = paste(legend_text, collapse = "\n"),
         cex = 0.9)
    
    # Add optional title
    if (nzchar(title)) {
      mtext(title, side = 3, line = -1, adj = 0, cex = 1.0, font = 2)
    }
    
    # Close the PNG device (required to save the file)
    dev.off()
  }
  
  # ===========================================================================
  # SECTION 3: FIRST-PASS FEATURE SELECTION CRITERIA
  # ===========================================================================
  # These six criteria each score genes before any gene-gene connectivity is
  # calculated. They are grouped here by pipeline role rather than by one shared
  # statistical property:
  #
  #   - Variance, MAD and MeanAbsDev describe expression spread using different
  #     deviation measures.
  #   - HVG scores variance relative to the mean-variance trend.
  #   - Entropy scores the shape and occupancy of each gene's expression
  #     distribution after binning.
  #   - PCA loadings score contribution to the leading multivariate axes.
  #
  # The union of these six top-N lists becomes the candidate pool for the
  # connectivity-based criteria.
  
  if (!quiet) message("[INFO] Running six first-pass feature selection criteria...")

  take_top_names <- function(scores, n) {
    genes <- names(sort(scores, decreasing = TRUE))
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0 || n <= 0) {
      return(character(0))
    }
    genes[seq_len(min(n, length(genes)))]
  }
  
  # ---------------------------------------------------------------------------
  # METHOD 1: VARIANCE
  # ---------------------------------------------------------------------------
  # Genes with the highest variance across samples are selected. The rationale
  # is that genes exhibiting greater variability are more likely to be
  # biologically informative (e.g., differentially expressed between subtypes).
  #
  # FORMULA: Var(x) = (1/n) × Σ(xᵢ - μ)²
  #
  # ADVANTAGES:
  #   - Simple and intuitive interpretation
  #   - Computationally efficient
  #
  # LIMITATIONS:
  #   - Can be dominated by outliers (single extreme values inflate variance)
  #   - Biased toward highly expressed genes (variance scales with mean)
  #
  # The rowVars() function from matrixStats is used for efficiency (~10× faster
  # than apply(counts, 1, var)).
  
  top_var <- take_top_names(rowVars(counts, na.rm = TRUE), top_n_method)

  # ---------------------------------------------------------------------------
  # METHOD 1b: HVG (MEAN-VARIANCE TREND ADJUSTED)
  # ---------------------------------------------------------------------------
  # Raw variance is biased toward highly expressed genes because variance
  # scales with mean. A true HVG method removes this bias by selecting genes
  # with higher variance than expected given their average expression level.
  #
  # FORMULA: HVG score = log10(var(x) + ε) - E[log10(var) | mean(x)]
  #
  # WHERE:
  #   - E[log10(var) | mean(x)] is the expected log-variance at that mean,
  #     estimated by a LOESS smooth of log10(variance) ~ mean expression
  #   - A positive residual means the gene is more variable than expected
  #
  # ADVANTAGES:
  #   - Mean-variance trend corrected (unlike raw Variance)
  #   - Identifies biologically variable genes regardless of expression level
  #   - Equivalent to the HVG selection used in Seurat/scran pipelines

  gene_mean <- rowMeans(counts, na.rm = TRUE)
  gene_var  <- matrixStats::rowVars(counts, na.rm = TRUE)
  names(gene_var) <- rownames(counts)

  eps <- 1e-8
  ok  <- is.finite(gene_mean) & is.finite(gene_var) & gene_var > 0

  hvg_fit <- loess(
    log10(gene_var[ok] + eps) ~ gene_mean[ok],
    span   = 0.3,
    degree = 2
  )

  expected_logvar <- rep(NA_real_, length(gene_var))
  names(expected_logvar) <- rownames(counts)
  expected_logvar[ok] <- predict(hvg_fit, newdata = gene_mean[ok])

  hvg_score <- log10(gene_var + eps) - expected_logvar
  names(hvg_score) <- rownames(counts)
  hvg_score[!is.finite(hvg_score)] <- -Inf

  top_hvg <- take_top_names(hvg_score, top_n_method)

  # ---------------------------------------------------------------------------
  # METHOD 2: MEDIAN ABSOLUTE DEVIATION (MAD)
  # ---------------------------------------------------------------------------
  # MAD provides a robust measure of variability that is less sensitive to
  # outliers than variance. This metric is particularly useful when the data
  # may contain extreme values or non-normal distributions.
  #
  # FORMULA: MAD = median(|xᵢ - median(x)|)
  #
  # INTERPRETATION:
  # MAD measures how far typical observations are from the median. Unlike
  # variance, a single extreme value has minimal impact on MAD.
  #
  # ADVANTAGES:
  #   - Robust to outliers (breakdown point of 50%)
  #   - Works well for non-normal distributions
  #
  # LIMITATIONS:
  #   - Less efficient than variance for normal distributions
  #   - May miss genes with important outlier expression patterns
  
  top_mad <- take_top_names(apply(counts, 1, mad, na.rm = TRUE), top_n_method)
  
  # ---------------------------------------------------------------------------
  # METHOD 3: MEAN ABSOLUTE DEVIATION FROM MEAN
  # ---------------------------------------------------------------------------
  # This method is similar to MAD but uses the mean instead of the median as
  # the centre point, providing a different perspective on dispersion.
  #
  # FORMULA: MeanAbsDev = (1/n) × Σ|xᵢ - mean(x)|
  #
  # COMPARISON TO MAD:
  #   - Uses mean instead of median (more sensitive to outliers)
  #   - Uses mean of deviations instead of median (less robust)
  #   - Falls between variance (sensitive) and MAD (robust) in outlier handling
  #
  # The sweep() function efficiently subtracts row means from each row,
  # which is faster than a loop-based approach.
  
  mean_vals <- rowMeans(counts, na.rm = TRUE)
  absdev <- rowMeans(abs(sweep(counts, 1, mean_vals, "-")), na.rm = TRUE)
  top_mean_absdev <- take_top_names(absdev, top_n_method)
  
  # ---------------------------------------------------------------------------
  # METHOD 4: SHANNON ENTROPY
  # ---------------------------------------------------------------------------
  # Shannon entropy measures the "complexity" or "uncertainty" of a gene's
  # expression distribution. Genes with high entropy have expression values
  # spread across many levels, while low-entropy genes have expression
  # concentrated at few values.
  #
  # FORMULA: H(x) = -Σ p(xᵢ) × log₂(p(xᵢ))
  #
  # WHERE:
  #   - p(xᵢ) is the probability of expression level xᵢ
  #   - Expression values are binned into 20 discrete bins
  #   - log₂ is used so entropy is measured in bits
  #
  # INTERPRETATION:
  #   - H = 0: All values in one bin (constant expression)
  #   - H = log₂(20) ≈ 4.32: Uniform distribution across all bins
  #
# WHY USE ENTROPY?
# High-entropy genes may reflect more nuanced biological variation.
  # For example, a gene expressed at three distinct levels (low/medium/high)
  # corresponding to three cell states would have higher entropy than a
  # gene with only two states.
  #
  # ADVANTAGES:
  #   - Captures distributional complexity, not just spread
  #   - Can identify genes with multi-modal expression patterns
  #
  # LIMITATIONS:
  #   - Sensitive to binning strategy (number of bins, bin edges)
  #   - Computationally more expensive than variance-based methods
  
  entropy_score <- apply(counts, 1, function(x) {
    x <- x[!is.na(x)]  # Remove missing values
    
    # Genes with no variation (constant expression) have zero entropy
    if (length(x) < 2 || var(x) == 0) return(0)
    
    # Bin expression values into 20 bins and compute frequencies
    # cut() assigns each value to a bin; table() counts values per bin
    p <- table(cut(x, breaks = 20)) / length(x)
    
    # Remove zero-probability bins to avoid log(0), which is undefined
    p <- p[p > 0]
    
    # Apply Shannon entropy formula: H = -Σ p × log₂(p)
    -sum(p * log2(p))
  })
  top_entropy <- take_top_names(entropy_score, top_n_method)
  
  # ---------------------------------------------------------------------------
  # METHOD 5: PCA LOADINGS
  # ---------------------------------------------------------------------------
  # Principal Component Analysis (PCA) identifies the directions of maximum
  # variance in the data. Genes with high loadings on the top principal
  # components are those that contribute most to explaining the overall
  # structure of the data.
  #
  # APPROACH:
  # 1. Centre the data (subtract mean from each gene)
  # 2. Perform PCA on samples (genes as features)
  # 3. Extract loadings (gene contributions) for top 5 PCs
  # 4. Sum absolute loadings across PCs as the gene importance score
  #
# WHY SUM ABSOLUTE LOADINGS?
# A gene might load heavily on PC2 but not PC1. By summing across multiple
# PCs, genes important for any major axis of variation are included.
  # Absolute values are used because the sign of a loading is arbitrary
  # (PC directions can be flipped).
  #
  # ADVANTAGES:
  #   - Captures multivariate structure (gene-gene relationships)
  #   - Identifies genes important for major data patterns
  #
  # LIMITATIONS:
  #   - Assumes linear relationships
  #   - May miss genes important for minor (but biologically relevant) variation
  
  # Centre the data by subtracting the mean from each gene (row)
  # scale(center=TRUE, scale=FALSE) mean-centres without dividing by SD
  X_scaled <- t(scale(t(counts), center = TRUE, scale = FALSE))
  
  # Perform PCA on the transposed matrix
  # prcomp expects samples as rows and features (genes) as columns
  # center=FALSE because the data is already centred
  pca <- prcomp(t(X_scaled), center = FALSE, scale. = FALSE)
  
  # Extract loadings (rotation matrix) for the top 5 PCs
  # pca$rotation has genes as rows and PCs as columns
  # Each entry represents how much that gene contributes to that PC
  n_pcs <- min(5, ncol(pca$rotation))
  loadings <- abs(pca$rotation[, 1:n_pcs, drop = FALSE])
  
  # Sum absolute loadings across PCs to get an overall importance score
  gene_scores <- rowSums(loadings)
  top_pca <- take_top_names(gene_scores, top_n_method)
  
  # ===========================================================================
# SECTION 4: COMBINE CANDIDATES FROM ALL METHODS
# ===========================================================================
# The union of genes selected by any method forms the candidate pool for
# network-based refinement. This ensures that genes identified by only one
# method are not missed.
  #
  # EXAMPLE:
  # If each method selects 3000 genes with 50% pairwise overlap, the union
  # might contain ~5000-8000 unique genes, depending on method agreement.
  
  raw_candidates <- unique(c(top_var, top_hvg, top_mad, top_mean_absdev, top_entropy, top_pca))
  all_candidates <- raw_candidates[!is.na(raw_candidates) & raw_candidates %in% rownames(counts)]
  dropped_candidates <- length(raw_candidates) - length(all_candidates)
  if (!quiet && dropped_candidates > 0) {
    message("[WARN] Dropped ", dropped_candidates,
            " candidate genes absent from the expression matrix or missing IDs")
  }
  if (length(all_candidates) == 0) {
    stop("No candidate genes remain after filtering to expression matrix rownames.")
  }
  # The connectivity-based criteria are computed on the union of the first-pass
  # selected genes. This is intentionally permissive: a gene selected by any one
  # first-pass criterion remains eligible for connectivity scoring.
  #
  # This is a union, not an intersection, and it is not limited to variance-selected
  # genes. The union preserves method-specific candidates before the smaller
  # connectivity-based lists are derived.
  if (!quiet) message("[INFO] Union of candidates: ", length(all_candidates), " genes")
  
  # ===========================================================================
  # SECTION 5: CONNECTIVITY-BASED FEATURE SELECTION CRITERIA
  # ===========================================================================
  # These criteria operate on the first-pass candidate union. Unlike the
  # first-pass criteria, they use gene-gene relationships within the candidate
  # pool.
  #
  # In this script, a hub gene means a co-expression-central gene: a gene with
  # high total or average absolute connectivity to other candidate genes. This is
  # a network-position measure, not evidence that the gene is a causal regulator.
  #
  # Spearman connectivity ranks genes by mean absolute rank correlation to the
  # candidate pool. MX combines this connectivity with scaled variance. WGCNA
  # kTotal ranks genes by total unsigned soft-thresholded network connectivity.
  
  # ---------------------------------------------------------------------------
  # METHOD 6: SPEARMAN CONNECTIVITY
  # ---------------------------------------------------------------------------
  # Spearman correlation measures monotonic relationships between genes.
  # For each gene, the average absolute correlation with all other genes
  # is computed. Genes with high average correlation are "hub" genes.
  #
  # WHY SPEARMAN INSTEAD OF PEARSON?
  #   - Spearman is rank-based, making it robust to outliers
  #   - Captures monotonic non-linear relationships
  #   - More appropriate for gene expression data, which often has outliers
  #
  # MATHEMATICAL NOTE:
  # Spearman correlation is equivalent to Pearson correlation computed on
  # ranked data. This script exploits this fact for computational efficiency:
  #   1. Rank all expression values once (O(n log n) per gene)
  #   2. Compute Pearson correlations on ranks (faster than Spearman directly)
  
  if (!quiet) message("[INFO] Computing Spearman correlation on candidates...")
  
  # Subset the expression matrix to candidate genes only
  # This reduces the correlation matrix from ~20000² to ~5000² entries,
  # dramatically reducing memory usage and computation time
  X_cand <- counts[all_candidates, , drop = FALSE]
  
  if (!quiet) message("[DEBUG] X_cand dimensions: ", nrow(X_cand), " × ", ncol(X_cand))
  if (!quiet) message("[DEBUG] X_cand rownames sample (first 5): ",
                      paste(head(rownames(X_cand), 5), collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # MEMORY-EFFICIENT BLOCK-WISE CORRELATION COMPUTATION
  # ---------------------------------------------------------------------------
  # Computing a full correlation matrix for 5000+ genes would require
  # storing a 5000 × 5000 matrix (~200MB) plus temporary memory during
  # computation. This can cause memory issues on systems with limited RAM.
  #
  # SOLUTION: Compute correlations in blocks and only store the running
  # sum needed for the mean. This keeps memory usage bounded.
  #
  # ALGORITHM:
  # 1. Pre-rank all genes once (Spearman = Pearson on ranks)
  # 2. For each block of 500 genes:
  #    a. Compute correlations between block and all genes
  #    b. Store only the row means (average correlation per gene)
  #    c. Free the correlation block immediately
  # 3. Result: average absolute correlation for each gene
  
  if (!quiet) message("[INFO] Computing Spearman connectivity (block-wise, memory-efficient)...")
  
  X <- X_cand
  g <- nrow(X)               # Number of genes
  block <- 500               # Block size (adjustable based on available RAM)
  avg_corr <- numeric(g)     # Will store average absolute correlation per gene
  names(avg_corr) <- rownames(X)
  
  # PRE-RANK ROWS
  # matrixStats::rowRanks() is much faster than apply(X, 1, rank)
  # ties.method = "average" assigns tied values the average of their ranks
  # preserveShape = TRUE maintains the matrix dimensions
  if (!quiet) message("[INFO] Ranking expression values for Spearman correlation...")
  Xr <- matrixStats::rowRanks(X, ties.method = "average", preserveShape = TRUE)
  
  # BLOCK-WISE CORRELATION COMPUTATION
  if (!quiet) message("[INFO] Computing block-wise correlations (block size: ", block, ")...")
  for (i in seq(1, g, by = block)) {
    i2 <- min(g, i + block - 1)  # End index for this block
    Xi <- Xr[i:i2, , drop = FALSE]  # Extract block of ranked genes
    
    # Compute Pearson correlation between block and all genes
    # cor(t(Xi), t(Xr)) gives a (block_size × g) matrix
    # Using "pearson" on ranks is equivalent to Spearman
    C <- cor(t(Xi), t(Xr), method = "pearson", use = "pairwise.complete.obs")
    
    # Store only the mean of absolute correlations for each gene in this block
    avg_corr[i:i2] <- rowMeans(abs(C), na.rm = TRUE)
    
    if (!quiet) message(sprintf("[INFO] Spearman blocks done: %d-%d / %d", i, i2, g))
    
    # Free memory immediately after each block
    # rm() removes the object; gc() triggers garbage collection
    rm(C)
    gc()
  }
  
  # Free the ranked matrix (no longer needed)
  rm(Xr)
  gc()
  
  if (!quiet) message("[DEBUG] avg_corr length: ", length(avg_corr))
  if (!quiet) message("[DEBUG] avg_corr names sample (first 5): ",
                      paste(head(names(avg_corr), 5), collapse = ", "))
  
  # Select top genes by Spearman connectivity
  final_top500_spearman <- take_top_names(avg_corr, final_top)
  final_top3000_spearman <- take_top_names(avg_corr, top_n_method)
  
  # ---------------------------------------------------------------------------
  # METHOD 7: MX SCORE (SPEARMAN × VARIANCE)
  # ---------------------------------------------------------------------------
  # The MX score combines network connectivity with variability.
  #
  # RATIONALE:
  # A gene might be highly connected but have low variance (e.g., a stably
  # expressed housekeeping gene that co-varies with many others). Such genes
  # are not useful for distinguishing cell types or conditions.
  #
  # Conversely, a gene might be highly variable but not well-connected
  # (potentially technical noise rather than biological signal).
  #
  # By multiplying connectivity by variance, genes that are BOTH variable
  # AND connected are prioritised. These genes are likely to be biologically
  # meaningful regulators that distinguish different cellular states.
  #
  # FORMULA: MX = avg_correlation × (variance / max_variance)
  #
  # The variance is scaled to [0, 1] to make it comparable to correlation
  # (which is already bounded by [-1, 1]).
  
  if (!quiet) message("[INFO] Computing Spearman × Variance (MX)...")
  
  # Compute variance for candidate genes
  gene_var <- rowVars(X_cand, na.rm = TRUE)
  
  # CRITICAL: Force names onto gene_var
  # rowVars() may return an unnamed vector, which would break downstream
  # operations that rely on gene names for matching
  names(gene_var) <- rownames(X_cand)
  
  if (!quiet) {
    message("[DEBUG] gene_var length: ", length(gene_var))
    message("[DEBUG] gene_var names sample (first 5): ",
            paste(head(names(gene_var), 5), collapse = ", "))
    message("[DEBUG] gene_var and avg_corr name match: ",
            identical(names(gene_var), names(avg_corr)))
  }
  
  # Scale variance to [0, 1] range
  max_var <- max(gene_var, na.rm = TRUE)
  if (!quiet) message("[DEBUG] max_var: ", max_var)
  gene_var_scaled <- gene_var / max_var
  
  # Compute MX Score: product of connectivity and scaled variance
  mx_score <- avg_corr * gene_var_scaled
  
  # CRITICAL: Ensure mx_score has names
  # Element-wise operations in R can sometimes drop names
  names(mx_score) <- names(avg_corr)
  
  if (!quiet) message("[DEBUG] mx_score length: ", length(mx_score))
  if (!quiet) message("[DEBUG] mx_score names sample (first 5): ",
                      paste(head(names(mx_score), 5), collapse = ", "))
  
  # Select top genes by MX score
  final_top500_mx <- take_top_names(mx_score, final_top)
  final_top3000_mx <- take_top_names(mx_score, top_n_method)
  
  if (!quiet) message("[DEBUG] final_top500_mx length: ", length(final_top500_mx))
  
  # Free memory before WGCNA step
  if (!quiet) message("[INFO] Freeing memory before WGCNA step...")
  gc()
  
  # ===========================================================================
  # SECTION 6: WGCNA kTOTAL (METHOD 8)
  # ===========================================================================
  # WGCNA (Weighted Gene Co-expression Network Analysis) provides a more
  # sophisticated network-based approach using "soft thresholding."
  #
  # KEY CONCEPT: SOFT THRESHOLDING
  # Instead of using a hard cutoff to define whether genes are connected
  # (correlation > 0.7 → connected, otherwise not), soft thresholding raises
  # correlations to a power:
  #
  #   connection_weight = |correlation|^power
  #
  # This preserves the continuous nature of correlations while emphasising
  # stronger connections. Higher powers create sparser networks with only
  # the strongest connections retained.
  #
  # kTOTAL: TOTAL CONNECTIVITY
  # For each gene, kTotal is the sum of its connection weights to all other
  # genes. High-kTotal genes are network "hubs" - they are strongly correlated
  # with many other genes and likely play central regulatory roles.
  
  if (!quiet) message("[INFO] Running WGCNA soft threshold and kTotal...")
  
  # Enable WGCNA multi-threading for faster computation
  # The number of threads is read from SLURM environment variable if available,
  # otherwise defaults to 8
  n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
  WGCNA::allowWGCNAThreads(nThreads = n_threads)
  if (!quiet) message("[INFO] WGCNA using ", n_threads, " threads")
  
  # ---------------------------------------------------------------------------
  # DATA PREPARATION FOR WGCNA
  # ---------------------------------------------------------------------------
  # CRITICAL: WGCNA expects data in the opposite orientation from typical
  # expression matrices:
  #   - Standard format: genes (rows) × samples (columns)
  #   - WGCNA format: samples (rows) × genes (columns)
  #
  # Additionally, WGCNA works best with z-scored data (mean=0, sd=1 per gene).
  
  if (!quiet) message("[DEBUG] Using candidate genes for WGCNA (", nrow(X_cand), " genes)")
  if (!quiet) message("[DEBUG] X_cand dimensions (genes×samples): ",
                      nrow(X_cand), " × ", ncol(X_cand))
  
  # Transpose: genes×samples → samples×genes
  datExpr <- t(X_cand)
  
  # Z-score per gene (per column in the transposed matrix)
  # scale() with default arguments centres and scales each column
  datExpr <- scale(datExpr, center = TRUE, scale = TRUE)
  
  # Convert to data.frame (WGCNA often prefers this format)
  datExpr <- as.data.frame(datExpr)
  
  if (!quiet) message("[DEBUG] datExpr dimensions (samples×genes): ",
                      nrow(datExpr), " × ", ncol(datExpr))
  if (!quiet) message("[DEBUG] datExpr colnames sample (first 5 genes): ",
                      paste(head(colnames(datExpr), 5), collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # SOFT THRESHOLD SELECTION
  # ---------------------------------------------------------------------------
  # The soft threshold power is chosen to achieve approximate "scale-free"
  # network topology, which is characteristic of biological networks.
  #
  # SCALE-FREE NETWORKS:
  # In a scale-free network, the degree distribution follows a power law:
  #   P(k) ∝ k^(-γ)
  # Most nodes have few connections, but a few "hub" nodes have many.
  # This is typical of biological networks like protein-protein interactions
  # and gene regulatory networks.
  #
  # SELECTION CRITERIA:
  #   1. Scale-free topology fit R² ≥ 0.8 (good fit to power law)
  #   2. Slope between -1.5 and -2.5 (characteristic of biological networks)
  #
  # If no power meets these criteria, a default of 5 is used.
  
  powers <- c(1:10, seq(12, 30, 2))  # Candidate powers to test
  sft <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = powers,
    verbose = 0,
    networkType = "unsigned"  # Consider both positive and negative correlations
  )
  
  sft_table <- sft$fitIndices
  
  # Identify powers meeting the criteria
  candidates_sft <- sft_table$SFT.R.sq >= 0.8 &
                    sft_table$slope <= -1.5 &
                    sft_table$slope >= -2.5
  
  # Select the power with highest R² among valid candidates
  softPower <- if (sum(candidates_sft) == 0) {
    5  # Default if no candidate meets criteria
  } else {
    sft_table$Power[candidates_sft][which.max(sft_table$SFT.R.sq[candidates_sft])]
  }
  if (!quiet) message("[INFO] Using soft power: ", softPower)
  
  # ---------------------------------------------------------------------------
  # COMPUTE ADJACENCY MATRIX AND kTOTAL
  # ---------------------------------------------------------------------------
  # The adjacency matrix A has entries:
  #   A_ij = |cor(gene_i, gene_j)|^power
  #
  # kTotal for gene i is the sum of its adjacencies:
  #   kTotal_i = Σⱼ A_ij
  #
  # This measures total network connectivity - how strongly a gene is
  # connected to all other genes in the network.
  
  if (!quiet) message("[DEBUG] Computing adjacency matrix...")
  adjacency <- WGCNA::adjacency(datExpr, power = softPower, type = "unsigned")
  
  if (!quiet) message("[DEBUG] Adjacency matrix computed: ",
                      nrow(adjacency), " × ", ncol(adjacency))
  if (!quiet) message("[DEBUG] Adjacency colnames sample (first 5 genes): ",
                      paste(head(colnames(adjacency), 5), collapse = ", "))
  
  # Compute kTotal per gene
  # Since genes are columns in datExpr, they are also columns in adjacency
  # colSums gives the sum of each column (total connectivity per gene)
  kTotal <- colSums(adjacency, na.rm = TRUE)
  names(kTotal) <- colnames(adjacency)
  
  # Free adjacency matrix immediately (saves ~200MB for 5000 genes)
  rm(adjacency)
  gc()
  
  if (!quiet) {
    message("[DEBUG] kTotal length: ", length(kTotal))
    message("[DEBUG] kTotal has NAs: ", sum(is.na(kTotal)))
    message("[DEBUG] kTotal has zeros: ", sum(kTotal == 0, na.rm = TRUE))
    message("[DEBUG] kTotal has names: ", !is.null(names(kTotal)))
    message("[DEBUG] kTotal names sample (first 5): ",
            paste(head(names(kTotal), 5), collapse = ", "))
  }
  
  # ===========================================================================
  # SECTION 7: JOINT RANKING - COMBINE kTOTAL AND MX SCORES
  # ===========================================================================
  # A comprehensive ranking is created that considers both WGCNA connectivity
  # and the MX score, providing multiple perspectives on gene importance.
  
  if (!quiet) {
    message("[DEBUG] Before intersection:")
    message("[DEBUG]   kTotal length: ", length(kTotal))
    message("[DEBUG]   kTotal names sample (first 5): ",
            paste(head(names(kTotal), 5), collapse = ", "))
    message("[DEBUG]   mx_score length: ", length(mx_score))
    message("[DEBUG]   mx_score names sample (first 5): ",
            paste(head(names(mx_score), 5), collapse = ", "))
  }
  
  # SANITY CHECKS
  # Both vectors must have names for intersection to work correctly
  stopifnot(!is.null(names(kTotal)),
            !is.null(names(mx_score)),
            length(names(kTotal)) > 0,
            length(names(mx_score)) > 0)
  
  # Identify genes present in both scoring systems
  genes <- intersect(names(kTotal), names(mx_score))
  
  if (!quiet) {
    message("[DEBUG] After intersection:")
    message("[DEBUG]   Intersection length: ", length(genes))
    message("[DEBUG]   Intersection sample (first 5): ",
            paste(head(genes, 5), collapse = ", "))
  }
  
  # Validate intersection size
  # Both methods were computed on the same candidate set, so the intersection
  # should be substantial (close to the full candidate set size)
  if (length(genes) < 100) {
    warning("[WARNING] Intersection is very small (", length(genes), " genes)!")
    warning("[WARNING] Expected ~", length(all_candidates), " genes (candidate set size)")
    
    # Debug: show which genes are missing from each set
    kTotal_only <- setdiff(names(kTotal), names(mx_score))
    mx_only <- setdiff(names(mx_score), names(kTotal))
    if (!quiet) {
      message("[DEBUG] Genes only in kTotal (first 10): ",
              paste(head(kTotal_only, 10), collapse = ", "))
      message("[DEBUG] Genes only in mx_score (first 10): ",
              paste(head(mx_only, 10), collapse = ", "))
      message("[DEBUG] Total kTotal_only: ", length(kTotal_only))
      message("[DEBUG] Total mx_only: ", length(mx_only))
    }
  } else {
    if (!quiet) message("[DEBUG] Intersection size looks good: ", length(genes), " genes")
  }
  
  # Final validation
  stopifnot(length(genes) > 100)
  
  # ---------------------------------------------------------------------------
  # BUILD COMPREHENSIVE RANKING TABLE
  # ---------------------------------------------------------------------------
  # This table contains all scores and multiple ranking approaches for each
  # gene, enabling detailed comparison of different ranking strategies.
  
  if (length(genes) == 0) {
    # Handle edge case of empty intersection (should not happen with proper data)
    if (!quiet) message("[WARNING] Creating empty scores_df because no genes in intersection")
    scores_df <- dplyr::tibble(
      Gene = character(0),
      WGCNA_kTotal = numeric(0),
      Spearman_meanAbs = numeric(0),
      Variance = numeric(0),
      Variance_scaled = numeric(0),
      MX_Score = numeric(0),
      Rank_kTotal = integer(0),
      Rank_MX = integer(0),
      Pctl_kTotal = numeric(0),
      Pctl_MX = numeric(0),
      Rank_Sum = integer(0),
      Rank_Prod = integer(0),
      Rank_Combo = integer(0)
    )
  } else {
    if (!quiet) message("[DEBUG] Creating scores_df with ", length(genes), " genes")
    
    scores_df <- dplyr::tibble(
      Gene = genes,
      WGCNA_kTotal = kTotal[genes],           # WGCNA total connectivity
      Spearman_meanAbs = avg_corr[genes],     # Mean absolute Spearman correlation
      Variance = gene_var[genes],              # Raw variance
      Variance_scaled = gene_var_scaled[genes], # Variance scaled to [0, 1]
      MX_Score = mx_score[genes]               # Combined connectivity × variance
    ) %>%
      dplyr::mutate(
        # RANK BY EACH METHOD
        # rank() with negative values so that highest scores get rank 1
        # ties.method = "min" assigns the same rank to tied values
        Rank_kTotal = rank(-WGCNA_kTotal, ties.method = "min"),
        Rank_MX = rank(-MX_Score, ties.method = "min"),
        
        # PERCENTILE SCORES
        # Transform ranks to percentiles (1 = top, 0 = bottom)
        # This allows comparison across different ranking methods
        Pctl_kTotal = if (n() > 1) 1 - (Rank_kTotal - 1) / (n() - 1) else numeric(0),
        Pctl_MX = if (n() > 1) 1 - (Rank_MX - 1) / (n() - 1) else numeric(0),
        
        # COMBINED RANKINGS
        # Rank_Sum: Sum of ranks (low = good in both methods)
        # Rank_Prod: Product of ranks (strongly penalises genes bad in either)
        Rank_Sum = Rank_kTotal + Rank_MX,
        Rank_Prod = Rank_kTotal * Rank_MX,
        
        # FINAL COMBINED RANK
        # Primarily by sum, with product as a tiebreaker
        # The 0.001 factor ensures product only affects ties
        Rank_Combo = rank(Rank_Sum + 0.001 * sqrt(Rank_Prod))
      ) %>%
      dplyr::arrange(Rank_Combo)  # Sort by combined rank (best first)
  }
  
  if (!quiet) message("[DEBUG] scores_df dimensions: ",
                      nrow(scores_df), " rows × ", ncol(scores_df), " cols")
  
  # ---------------------------------------------------------------------------
  # EXTRACT TOP GENES BY EACH METHOD
  # ---------------------------------------------------------------------------
  
  if (nrow(scores_df) > 0) {
    # Top genes by kTotal (WGCNA connectivity)
    top500_kTotal <- dplyr::arrange(scores_df, Rank_kTotal) %>%
      dplyr::slice_head(n = min(final_top, nrow(scores_df))) %>%
      dplyr::pull(Gene)
    
    top3000_kTotal <- dplyr::arrange(scores_df, Rank_kTotal) %>%
      dplyr::slice_head(n = min(top_n_method, nrow(scores_df))) %>%
      dplyr::pull(Gene)
    
    # Top genes by MX score
    top500_MX <- dplyr::arrange(scores_df, Rank_MX) %>%
      dplyr::slice_head(n = min(final_top, nrow(scores_df))) %>%
      dplyr::pull(Gene)
    
    # Overlap between methods (genes selected by both kTotal and MX)
    overlap_500 <- intersect(top500_kTotal, top500_MX)
    
    if (!quiet) {
      message("[DEBUG] Top genes extracted:")
      message("[DEBUG]   top500_kTotal: ", length(top500_kTotal))
      message("[DEBUG]   top500_MX: ", length(top500_MX))
      message("[DEBUG]   overlap_500: ", length(overlap_500))
    }
  } else {
    # Handle empty results
    if (!quiet) message("[WARNING] scores_df is empty, creating empty top gene lists")
    top500_kTotal <- character(0)
    top3000_kTotal <- character(0)
    top500_MX <- character(0)
    overlap_500 <- character(0)
  }
  
  # ===========================================================================
  # SECTION 8: UPSET PLOTS - VISUALISE METHOD AGREEMENT
  # ===========================================================================
  # UpSet plots effectively visualise which genes are selected by which
  # methods and the degree of overlap between methods.
  #
  # Two plots are created:
  #   A. Top 3000: Shows broader agreement at the candidate selection stage
  #   B. Top 500: Shows agreement at the final selection stage
  
  if (!quiet) message("[INFO] Creating UpSet plots comparing all 9 methods...")
  
  # Safety checks for required vectors
  stopifnot(!is.null(names(avg_corr)), length(names(avg_corr)) > 0)
  stopifnot(!is.null(names(mx_score)), length(names(mx_score)) > 0)
  stopifnot(!is.null(names(kTotal)), length(names(kTotal)) > 0)
  
  # ---------------------------------------------------------------------------
  # UPSET PLOT A: TOP 3000 (Candidate Selection Stage)
  # ---------------------------------------------------------------------------
  # This plot compares all methods at the same threshold (top 3000 genes),
  # showing which methods agree on gene selection at the broader level.
  
  sets_3000 <- list(
    Variance   = top_var[1:min(top_n_method, length(top_var))],
    HVG        = top_hvg[1:min(top_n_method, length(top_hvg))],
    MAD        = top_mad[1:min(top_n_method, length(top_mad))],
    MeanAbsDev = top_mean_absdev[1:min(top_n_method, length(top_mean_absdev))],
    Entropy    = top_entropy[1:min(top_n_method, length(top_entropy))],
    PCA        = top_pca[1:min(top_n_method, length(top_pca))],
    Spearman   = final_top3000_spearman,
    MX         = final_top3000_mx,
    kTotal     = top3000_kTotal
  )
  
  # Ensure uniqueness and remove empty sets
  sets_3000 <- lapply(sets_3000, unique)
  nonempty_3000 <- vapply(sets_3000, length, integer(1)) > 0
  sets_3000_use <- sets_3000[nonempty_3000]
  
  if (length(sets_3000_use) >= 2) {
    # Write genes unique to each method
    # These are genes selected by one method but not by any other
    for (nm in names(sets_3000_use)) {
      others_union <- unique(unlist(sets_3000_use[names(sets_3000_use) != nm], use.names = FALSE))
      only_nm <- setdiff(sets_3000_use[[nm]], others_union)
      if (length(only_nm) > 0) {
        writeLines(only_nm, file.path(subdir, paste0("unique_to_", nm, "_top3000.txt")))
      }
    }
    
    # Create UpSet plot with legend
    save_upset_with_legend(
      sets_3000_use,
      file.path(subdir, "upset_all8_top3000_with_legend.png"),
      title = "UpSet: Top 3000 genes per method (9 methods)"
    )
    
    if (!quiet) message("[INFO] UpSet plot (top3000) saved: upset_all8_top3000_with_legend.png")
  } else {
    warning("[WARNING] Insufficient method sets (top3000) for UpSet plot, skipping")
  }
  
  # ---------------------------------------------------------------------------
  # UPSET PLOT B: TOP 500 (Final Selection Stage)
  # ---------------------------------------------------------------------------
  # This plot compares final shortlists to identify genes consistently
  # selected across multiple methods.
  
  sets_500 <- list(
    Variance   = top_var[1:min(final_top, length(top_var))],
    HVG        = top_hvg[1:min(final_top, length(top_hvg))],
    MAD        = top_mad[1:min(final_top, length(top_mad))],
    MeanAbsDev = top_mean_absdev[1:min(final_top, length(top_mean_absdev))],
    Entropy    = top_entropy[1:min(final_top, length(top_entropy))],
    PCA        = top_pca[1:min(final_top, length(top_pca))],
    Spearman   = final_top500_spearman,
    MX         = final_top500_mx,
    kTotal     = top500_kTotal
  )
  
  sets_500 <- lapply(sets_500, unique)
  nonempty_500 <- vapply(sets_500, length, integer(1)) > 0
  sets_500_use <- sets_500[nonempty_500]
  
  # ---------------------------------------------------------------------------
  # EMIT ALL METHOD GENE LISTS FOR DOWNSTREAM PIPELINE FAN-OUT
  # ---------------------------------------------------------------------------
  # Two output directories are maintained:
  #
  #   feature_sets_top{final_top}/ — legacy flat directory; all methods at
  #     final_top (500 by default).  Kept for backward compatibility with
  #     non-Snakemake callers and the multicohort benchmark script.
  #
  #   feature_sets/ — canonical Snakemake output directory.  Each method is
  #     written with its own per-method top-N value derived from method_topn_map
  #     (passed via --method_topn).  File naming: genes_top{N}_{method}.txt.
  #     This is what the Snakefile rule output: declarations expect.

  # --- Legacy flat directory (all methods at final_top) ---
  feature_sets_dir <- file.path(subdir, paste0("feature_sets_top", final_top))
  dir.create(feature_sets_dir, showWarnings = FALSE, recursive = TRUE)

  if (!quiet) {
    message("[INFO] Emitting gene lists (legacy flat) for ", length(sets_500_use),
            " methods to: ", feature_sets_dir)
  }

  # Write one gene list file per method (legacy, all at final_top)
  for (nm in names(sets_500_use)) {
    if (length(sets_500_use[[nm]]) > 0) {
      fn <- file.path(feature_sets_dir, paste0("genes_top", final_top, "_", nm, ".txt"))
      writeLines(sets_500_use[[nm]], fn)
      if (!quiet) {
        message("[INFO]   [legacy] ", nm, ": ", length(sets_500_use[[nm]]),
                " genes -> ", basename(fn))
      }
    }
  }

  # --- Canonical Snakemake output directory (per-method top-N) ---
  # Directory name is 'feature_sets' (no top-N suffix) because different methods
  # use different N values.  File names encode the method-specific N.
  canonical_sets_dir <- file.path(subdir, "feature_sets")
  dir.create(canonical_sets_dir, showWarnings = FALSE, recursive = TRUE)

  # Build a combined gene-vector map: top3000 for broad methods, top500 for
  # network methods, using the per-gene vectors already computed above.
  all_method_vectors <- list(
    Variance   = top_var,
    HVG        = top_hvg,
    MAD        = top_mad,
    MeanAbsDev = top_mean_absdev,
    Entropy    = top_entropy,
    PCA        = top_pca,
    Spearman   = final_top3000_spearman,   # full ranked vector; sliced below
    MX         = c(final_top500_mx,        # MX canonical = top-500
                   names(sort(mx_score, decreasing = TRUE))),  # extend for safety
    kTotal     = c(top500_kTotal,
                   names(sort(kTotal,    decreasing = TRUE)))
  )
  # Re-rank MX and kTotal from their score vectors (already sorted) for safety.
  all_method_vectors[["MX"]]     <- names(sort(mx_score, decreasing = TRUE))
  all_method_vectors[["kTotal"]] <- names(sort(kTotal,   decreasing = TRUE))

  # Default per-method top-N (matches Snakefile defaults when --method_topn absent).
  default_topn <- list(
    Variance = 3000, HVG = 3000, MAD = 3000, MeanAbsDev = 3000,
    Entropy  = 3000, PCA = 3000, Spearman = 500, MX = 500, kTotal = 500
  )

  if (!quiet) {
    message("[INFO] Emitting canonical gene lists (per-method top-N) to: ", canonical_sets_dir)
  }

  for (nm in names(all_method_vectors)) {
    vec <- all_method_vectors[[nm]]
    if (is.null(vec) || length(vec) == 0) next
    vec <- unique(vec[!is.na(vec)])
    # Resolve per-method N: --method_topn takes precedence, then default.
    n_out_genes <- if (!is.null(method_topn_map[[nm]])) {
      as.integer(method_topn_map[[nm]])
    } else {
      as.integer(default_topn[[nm]] %||% final_top)
    }
    top_genes <- head(vec, n_out_genes)
    fn <- file.path(canonical_sets_dir, paste0("genes_top", n_out_genes, "_", nm, ".txt"))
    writeLines(top_genes, fn)
    if (!quiet) {
      message("[INFO]   [canonical] ", nm, " (top-", n_out_genes, "): ",
              length(top_genes), " genes -> ", basename(fn))
    }
  }
  
  if (length(sets_500_use) >= 2) {
    # Save method sizes
    method_sizes_500 <- dplyr::tibble(
      Method = names(sets_500_use),
      Size = vapply(sets_500_use, length, integer(1))
    )
    readr::write_tsv(method_sizes_500, file.path(subdir, "method_set_sizes_top500.tsv"))
    
    # Write unique genes per method
    for (nm in names(sets_500_use)) {
      others_union <- unique(unlist(sets_500_use[names(sets_500_use) != nm], use.names = FALSE))
      only_nm <- setdiff(sets_500_use[[nm]], others_union)
      if (length(only_nm) > 0) {
        writeLines(only_nm, file.path(subdir, paste0("unique_to_", nm, "_top500.txt")))
      }
    }
    
    save_upset_with_legend(
      sets_500_use,
      file.path(subdir, "upset_all8_top500_with_legend.png"),
      title = "UpSet: Top 500 genes per method (9 methods)"
    )
    
    if (!quiet) message("[INFO] UpSet plot (top500) saved: upset_all8_top500_with_legend.png")
  } else {
    warning("[WARNING] Insufficient method sets (top500) for UpSet plot, skipping")
  }
  
  # ---------------------------------------------------------------------------
  # MX vs kTotal DETAILED COMPARISON
  # ---------------------------------------------------------------------------
  # Generate detailed comparison between the two primary network-based methods.
  
  if (length(final_top500_mx) > 0 && length(top500_kTotal) > 0) {
    only_in_mx     <- setdiff(final_top500_mx, top500_kTotal)
    only_in_kTotal <- setdiff(top500_kTotal, final_top500_mx)
    in_both        <- intersect(final_top500_mx, top500_kTotal)
    
    writeLines(only_in_mx,     file.path(subdir, "only_in_MX_top500.txt"))
    writeLines(only_in_kTotal, file.path(subdir, "only_in_kTotal_top500.txt"))
    writeLines(in_both,        file.path(subdir, "overlap_MX_kTotal_top500.txt"))
    
    if (!quiet) {
      message("[INFO] MX vs kTotal comparison:")
      message("[INFO]   Only in MX: ", length(only_in_mx), " genes")
      message("[INFO]   Only in kTotal: ", length(only_in_kTotal), " genes")
      message("[INFO]   In both: ", length(in_both), " genes")
    }
  }
  
  # ===========================================================================
  # SECTION 9: SAVE ALL RESULTS
  # ===========================================================================
  
  if (!quiet) message("[INFO] Saving results to: ", subdir)
  
  # ---------------------------------------------------------------------------
  # INTERMEDIATE FILES: Individual method results
  # ---------------------------------------------------------------------------
  # These files enable exploration and QC of individual methods
  
  write_vec(top_var, "top3000_variance.txt")
  write_vec(top_hvg, "top3000_HVG.txt")
  write_vec(top_mad, "top3000_MAD.txt")
  write_vec(top_mean_absdev, "top3000_meanAbsDev.txt")
  write_vec(top_entropy, "top3000_entropy.txt")
  write_vec(top_pca, "top3000_PCAloadings.txt")
  
  # ---------------------------------------------------------------------------
  # PRIMARY OUTPUT FILES
  # ---------------------------------------------------------------------------
  # These files are referenced by Snakefile and downstream scripts
  
  profile_tag <- if (!is.null(opt$profile)) toupper(opt$profile) else "PROFILE"
  
  # If Snakemake (or another caller) provides explicit output paths, honour those
  # exactly; otherwise fall back to profile-tagged filenames under subdir.
  joint_tsv  <- opt$scores_tsv %||%
    file.path(subdir, sprintf("%s_MX_joint_ranks_kTotal_vs_MX.tsv", profile_tag))
  genes_file <- opt$mx_list %||%
    file.path(subdir, sprintf("%s_MX_genes_top%d.txt", profile_tag, final_top))
  
  if (!quiet) {
    message("[DEBUG] Writing output files:")
    message("[DEBUG]   joint_tsv: ", joint_tsv, " (", nrow(scores_df), " rows)")
    message("[DEBUG]   genes_file: ", genes_file, " (", length(final_top500_mx), " genes)")
  }
  
  readr::write_tsv(scores_df, joint_tsv)
  writeLines(final_top500_mx, genes_file)
  
  # Individual method top500 lists for debugging
  write_vec(top500_kTotal, "top500_WGCNA_kTotal.txt")
  write_vec(top500_MX, "top500_MX.txt")
  write_vec(overlap_500, "overlap_top500_kTotal_and_MX.txt")
  
  # ---------------------------------------------------------------------------
  # OVERLAP SUMMARY STATISTICS
  # ---------------------------------------------------------------------------
  # Jaccard Index measures similarity between two sets:
  #   J(A, B) = |A ∩ B| / |A ∪ B|
  # Values range from 0 (no overlap) to 1 (identical sets)
  
  union_500 <- length(union(top500_kTotal, top500_MX))
  jaccard <- if (union_500 > 0) length(overlap_500) / union_500 else 0
  
  overlap_summary <- dplyr::tibble(
    Metric = c("Top500_kTotal", "Top500_MX", "Overlap"),
    Count = c(length(top500_kTotal), length(top500_MX), length(overlap_500)),
    Jaccard_Index = c(NA, NA, jaccard)
  )
  readr::write_tsv(overlap_summary, file.path(subdir, "overlap_summary.tsv"))
  
  # ---------------------------------------------------------------------------
  # TOP 200 ANALYSIS (Higher Confidence Subset)
  # ---------------------------------------------------------------------------
  top200_spearman <- names(sort(avg_corr, decreasing = TRUE))[1:200]
  top200_mx <- names(sort(mx_score, decreasing = TRUE))[1:200]
  intersection_200 <- intersect(top200_spearman, top200_mx)
  
  write_vec(top200_spearman, "top200_spearman.txt")
  write_vec(top200_mx, "top200_spearman_variance.txt")
  write_vec(intersection_200, "intersection_top200_spearman_and_spearman_variance.txt")
  
  # ---------------------------------------------------------------------------
  # SCATTER PLOT: kTotal vs MX Score
  # ---------------------------------------------------------------------------
  # This visualisation helps assess agreement between the two main ranking methods.
  # A strong positive correlation indicates the methods identify similar hub genes.
  
  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = WGCNA_kTotal, y = MX_Score)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggplot2::labs(
      x = "WGCNA kTotal",
      y = "MX Score (Spearman × Variance)",
      title = "kTotal vs MX across genes"
    ) +
    theme_dimred(base_size = 12)  # Custom theme from brca_unsup_methods.R
  
  ggplot2::ggsave(
    file.path(subdir, "kTotal_vs_MX_scatter.png"),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # ---------------------------------------------------------------------------
  # CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  cat("\n=== UNSUPERVISED FEATURE SELECTION SUMMARY ===\n")
  cat(sprintf("Input matrix: %d genes × %d samples\n", n_genes, n_samples))
  cat(sprintf("Top %d per method\n", top_n_method))
  cat(sprintf("Union candidates: %d\n", length(all_candidates)))
  cat(sprintf("Final top %d (MX): %d\n", final_top, length(final_top500_mx)))
  cat(sprintf("Top %d kTotal: %d\n", final_top, length(top500_kTotal)))
  cat(sprintf("Overlap kTotal & MX: %d (Jaccard = %.3f)\n",
              length(overlap_500), overlap_summary$Jaccard_Index[3]))
  cat(sprintf("Top 200 intersection: %d\n", length(intersection_200)))
  cat(sprintf("All results saved in: %s\n", subdir))
  
  # ---------------------------------------------------------------------------
  # RETURN RESULTS
  # ---------------------------------------------------------------------------
  # A list of results is returned for programmatic access.
  # invisible() suppresses automatic printing when called from the command line.
  
  invisible(list(
    candidates = all_candidates,
    mx = final_top500_mx,
    kTotal = top500_kTotal,
    overlap = overlap_500,
    intersection_200 = intersection_200,
    scores_df = scores_df,
    joint_tsv = joint_tsv,
    genes_file = genes_file
  ))
}


# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# The null-coalescing operator is defined, returning y if x is NULL.
# This operator provides default values when configuration parameters are missing.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Extract parameters: command-line args take precedence over config
top_n_method <- opt$top_n_method %||%
                if (!is.null(cfg) && !is.null(cfg$feature_selection$top_n_method)) {
                  cfg$feature_selection$top_n_method
                } else {
                  3000
                }

final_top <- opt$final_top %||%
             if (!is.null(cfg) && !is.null(cfg$feature_selection$final_top)) {
               cfg$feature_selection$final_top
             } else {
               500
             }

# ---------------------------------------------------------------------------
# Parse --method_topn into a named integer vector.
# Format: "Variance:3000,MAD:3000,MX:500,HVG:3000,..."
# This overrides top_n_method and final_top for named methods, enabling the
# Snakefile to request per-method top-N values that match the declared outputs.
# ---------------------------------------------------------------------------
method_topn_map <- list()
if (!is.null(opt$method_topn) && nzchar(opt$method_topn)) {
  pairs <- strsplit(opt$method_topn, ",")[[1]]
  for (p in pairs) {
    kv <- strsplit(trimws(p), ":")[[1]]
    if (length(kv) == 2) {
      method_topn_map[[trimws(kv[1])]] <- as.integer(trimws(kv[2]))
    }
  }
  message("[FEATURE_SELECTION] Per-method top-N from --method_topn: ",
          paste(names(method_topn_map), unlist(method_topn_map), sep = "=", collapse = ", "))
}

# The main function is executed with parameters from command-line or config.
# Default values (3000 candidates, 500 final) are used if not specified.
run_unsupervised_feature_selection(
  vst_rds = opt$vst_rds,
  outdir  = fs_subdir,
  top_n_method = top_n_method,
  final_top    = final_top,
  method_topn_map = method_topn_map,
  seed         = 123,
  quiet        = FALSE
)
