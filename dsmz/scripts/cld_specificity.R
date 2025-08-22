#!/usr/bin/env Rscript

# Cancer/Lineage/Disease (CLD) Gene Specificity Analysis
# Identifies disease-specific, enhanced, and ubiquitous genes in DSMZ cell lines

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

options(stringsAsFactors = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(edgeR)
  library(pheatmap)
})

# Define null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(42)

# Configuration parameters
config <- list(
  # Input file paths
  dsmz_rds = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata_with_sample_name.csv",
  outdir = "/home/chu25/dsmz/results/cld_specificity",
  
  # Disease/lineage definition
  disease_col = NULL,  # Column name for disease; if NULL, derive from Origin
  
  # Expression filtering thresholds
  min_cpm_gene = 0.5,              # Minimum CPM for gene expression
  min_samples_per_disease = 2,     # Minimum samples per disease group
  
  # Gene specificity classification thresholds (applied to linear CPM)
  fc_enriched = 4,           # Fold-change for disease-enriched genes
  fc_group = 4,              # Fold-change for group-enriched genes  
  max_group_size = 5,        # Maximum size for enriched groups
  fc_enhanced = 2,           # Fold-change for disease-enhanced genes
  ubiquitous_min_cpm = 5,    # Minimum median CPM for ubiquitous genes
  ubiquitous_cv_max = 0.6    # Maximum CV for ubiquitous genes
)

# Create output directory
dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Map origin descriptions to broad organ categories
#' @param x Character vector of origin descriptions
#' @return Character vector of organ categories
derive_organ_categories <- function(x) {
  case_when(
    str_detect(x, regex("breast", TRUE)) ~ "Breast",
    str_detect(x, regex("lung", TRUE)) ~ "Lung",
    str_detect(x, regex("liver|hepato|bile", TRUE)) ~ "Liver/Biliary",
    str_detect(x, regex("adrenal", TRUE)) ~ "Adrenal",
    str_detect(x, regex("brain|gli", TRUE)) ~ "Brain/CNS",
    str_detect(x, regex(paste0("lymph|mantle|hodgkin|effusion|",
                              "b[- ]?cell|t[- ]?cell|nk|",
                              "natural\\s*killer"), TRUE)) ~ "Lymphoid",
    str_detect(x, regex("leukemia", TRUE)) ~ "Leukemia",
    str_detect(x, regex("myeloid|megakaryocytic|erythroid|monocytic", 
                       TRUE)) ~ "Myeloid",
    str_detect(x, regex("bone\\s*marrow", TRUE)) ~ "Bone Marrow",
    str_detect(x, regex("abdomen|pleural", TRUE)) ~ "Body Cavity",
    is.na(x) | x == "NONE" ~ "Unknown",
    TRUE ~ "Other"
  )
}

#' Calculate coefficient of variation
#' @param v Numeric vector
#' @return Coefficient of variation (CV)
calculate_cv <- function(v) {
  mu <- mean(v, na.rm = TRUE)
  if (!is.finite(mu) || mu == 0) return(Inf)
  sd(v, na.rm = TRUE) / mu
}

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load and validate DSMZ count and metadata files
#' @param counts_file Path to RDS file with count data
#' @param meta_file Path to CSV file with metadata
#' @return List with raw counts and metadata
load_dsmz_data <- function(counts_file, meta_file) {
  cat("[INFO] Loading DSMZ counts and metadata...\n")
  
  # Load data files
  dsmz_raw <- readRDS(counts_file)
  dsmz_meta <- read.csv(meta_file)
  
  # Validate required columns
  required_count_cols <- c("Ensembl_ID", "gene_name")
  stopifnot(all(required_count_cols %in% colnames(dsmz_raw)))
  stopifnot("sample_name" %in% names(dsmz_meta))
  
  cat(sprintf("[INFO] Raw data: %d genes x %d columns\n", 
              nrow(dsmz_raw), ncol(dsmz_raw)))
  cat(sprintf("[INFO] Metadata: %d samples\n", nrow(dsmz_meta)))
  
  list(raw = dsmz_raw, meta = dsmz_meta)
}

#' Build count matrix from raw DSMZ data
#' @param dsmz_raw Raw DSMZ data frame
#' @return Count matrix (genes x samples)
build_count_matrix <- function(dsmz_raw) {
  cat("[INFO] Building count matrix...\n")
  
  # Identify sample columns (exclude annotation columns)
  annotation_cols <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annotation_cols)
  
  # Convert sample columns to numeric
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  # Create count matrix
  counts <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(counts) <- dsmz_raw$Ensembl_ID
  
  # Remove columns that are entirely NA
  keep_col <- vapply(as.data.frame(counts), 
                    function(v) any(!is.na(v)), logical(1))
  
  if (!all(keep_col)) {
    n_dropped <- sum(!keep_col)
    message(sprintf("[WARN] Dropping %d columns that are all NA", 
                   n_dropped))
    counts <- counts[, keep_col, drop = FALSE]
  }
  
  counts
}

#' Align metadata with count matrix columns
#' @param dsmz_meta Metadata data frame
#' @param counts Count matrix
#' @return List with aligned metadata and count matrix
align_metadata_counts <- function(dsmz_meta, counts) {
  cat("[INFO] Aligning metadata with count matrix...\n")
  
  # Create sample_id column for matching
  dsmz_meta <- dsmz_meta %>% 
    mutate(sample_id = sample_name)
  
  # Find matched samples
  matched_samples <- intersect(colnames(counts), dsmz_meta$sample_id)
  
  # Filter to matched samples only
  counts <- counts[, matched_samples, drop = FALSE]
  dsmz_meta <- dsmz_meta %>% 
    filter(sample_id %in% matched_samples) %>% 
    distinct(sample_id, .keep_all = TRUE)
  
  # Ensure order matches
  stopifnot(identical(colnames(counts), dsmz_meta$sample_id))
  
  cat(sprintf("[INFO] Final matrix: %d genes x %d samples\n", 
              nrow(counts), ncol(counts)))
  
  list(meta = dsmz_meta, counts = counts)
}

# =============================================================================
# DISEASE CLASSIFICATION FUNCTIONS
# =============================================================================

#' Define disease/lineage labels for samples
#' @param dsmz_meta Metadata data frame
#' @param disease_col Column name for disease (if available)
#' @return Updated metadata with disease column
define_disease_labels <- function(dsmz_meta, disease_col = NULL) {
  cat("[INFO] Defining disease/lineage labels...\n")
  
  if (!is.null(disease_col) && disease_col %in% names(dsmz_meta)) {
    # Use existing disease column
    disease <- dsmz_meta[[disease_col]]
    cat(sprintf("[INFO] Using existing disease column: %s\n", disease_col))
  } else {
    # Derive from Origin column
    origin_col <- dsmz_meta$Origin %||% NA_character_
    disease <- derive_organ_categories(origin_col)
    cat("[INFO] Derived disease labels from Origin column\n")
  }
  
  # Clean disease names for R compatibility
  dsmz_meta$disease <- make.names(disease)
  
  dsmz_meta
}

#' Filter samples by minimum disease group size
#' @param dsmz_meta Metadata with disease labels
#' @param counts Count matrix
#' @param min_samples Minimum samples per disease
#' @return List with filtered metadata and counts
filter_by_disease_size <- function(dsmz_meta, counts, min_samples) {
  cat("[INFO] Filtering by disease group size...\n")
  
  # Count samples per disease
  disease_counts <- table(dsmz_meta$disease)
  keep_diseases <- names(disease_counts)[disease_counts >= min_samples]
  
  if (!length(keep_diseases)) {
    stop(sprintf("No diseases with >= %d samples. Lower the threshold.", 
                min_samples))
  }
  
  # Filter to diseases with sufficient samples
  keep_samples <- dsmz_meta$disease %in% keep_diseases
  counts_filtered <- counts[, keep_samples, drop = FALSE]
  meta_filtered <- dsmz_meta[keep_samples, , drop = FALSE]
  
  cat(sprintf("[INFO] Retained %d diseases with >= %d samples each\n", 
              length(keep_diseases), min_samples))
  
  # Print disease summary
  final_counts <- table(meta_filtered$disease)
  cat("[INFO] Final disease distribution:\n")
  print(final_counts)
  
  list(meta = meta_filtered, counts = counts_filtered)
}

# =============================================================================
# EXPRESSION ANALYSIS FUNCTIONS
# =============================================================================

#' Normalize counts to CPM (Counts Per Million)
#' @param counts Raw count matrix
#' @return CPM-normalized matrix
normalize_to_cpm <- function(counts) {
  cat("[INFO] Performing CPM normalization...\n")
  
  y <- DGEList(counts = counts)
  cpm_mat <- edgeR::cpm(y, log = FALSE, prior.count = 0)
  
  cat(sprintf("[INFO] CPM matrix dimensions: %d genes x %d samples\n",
              nrow(cpm_mat), ncol(cpm_mat)))
  
  cpm_mat
}

#' Calculate mean expression per disease
#' @param cpm_mat CPM-normalized expression matrix
#' @param dsmz_meta Metadata with disease labels
#' @return Matrix of disease mean expressions
calculate_disease_means <- function(cpm_mat, dsmz_meta) {
  cat("[INFO] Calculating per-disease mean expressions...\n")
  
  disease_levels <- sort(unique(dsmz_meta$disease))
  
  disease_means <- sapply(disease_levels, function(disease) {
    sample_ids <- dsmz_meta$sample_id[dsmz_meta$disease == disease]
    
    if (length(sample_ids) == 1L) {
      # Single sample - return as is
      cpm_mat[, sample_ids]
    } else {
      # Multiple samples - calculate row means
      rowMeans(cpm_mat[, sample_ids, drop = FALSE])
    }
  })
  
  colnames(disease_means) <- disease_levels
  
  cat(sprintf("[INFO] Disease means matrix: %d genes x %d diseases\n",
              nrow(disease_means), ncol(disease_means)))
  
  disease_means
}

#' Filter genes by minimum expression threshold
#' @param disease_means Matrix of disease mean expressions
#' @param min_cpm Minimum CPM threshold
#' @return Filtered disease means matrix
filter_expressed_genes <- function(disease_means, min_cpm) {
  cat("[INFO] Filtering genes by expression threshold...\n")
  
  # Keep genes with at least min_cpm in at least one disease
  expressed_genes <- rowMaxs(disease_means) >= min_cpm
  filtered_means <- disease_means[expressed_genes, , drop = FALSE]
  
  cat(sprintf("[INFO] Genes retained after CPM filter: %d (was %d)\n", 
              nrow(filtered_means), nrow(disease_means)))
  
  filtered_means
}

# =============================================================================
# GENE SPECIFICITY CLASSIFICATION
# =============================================================================

#' Classify a single gene's specificity pattern
#' @param expression_vector Numeric vector of disease mean CPM for one gene
#' @param thresholds List of classification thresholds
#' @return List with classification and associated disease label
classify_gene_specificity <- function(expression_vector, thresholds) {
  small_value <- 1e-6  # Small value to avoid division by zero
  
  # Sort diseases by expression level
  expression_order <- order(expression_vector, decreasing = TRUE)
  top_disease <- expression_order[1]
  second_disease <- if (length(expression_vector) >= 2) {
    expression_order[2]
  } else {
    expression_order[1]
  }
  
  top_expression <- expression_vector[top_disease]
  second_expression <- expression_vector[second_disease]
  remaining_expression <- expression_vector[-top_disease]
  
  # Count diseases with meaningful expression
  n_expressed <- sum(expression_vector >= thresholds$min_cpm_gene)
  median_remaining <- if (length(remaining_expression)) {
    median(remaining_expression)
  } else {
    0
  }
  
  # Classification 1: Disease-enriched
  # Only one disease expressed AND strong fold-change over second
  if (n_expressed == 1 && 
      (top_expression / (second_expression + small_value)) >= 
      thresholds$fc_enriched) {
    return(list(
      class = "enriched", 
      label = names(expression_vector)[top_disease]
    ))
  }
  
  # Classification 2: Group-enriched  
  # Small group (2-5) of diseases much higher than others
  max_group_size <- min(thresholds$max_group_size, length(expression_vector))
  
  for (group_size in 2:max_group_size) {
    group_diseases <- expression_order[1:group_size]
    group_expression <- expression_vector[group_diseases]
    non_group_expression <- expression_vector[-group_diseases]
    
    if (!length(non_group_expression)) next
    
    # All group members must be >= fc_group * median(non-group)
    threshold_expression <- thresholds$fc_group * 
                           (median(non_group_expression) + small_value)
    
    if (all(group_expression >= threshold_expression)) {
      group_label <- paste(names(expression_vector)[group_diseases], 
                          collapse = ",")
      return(list(
        class = paste0("group_enriched_", group_size), 
        label = group_label
      ))
    }
  }
  
  # Classification 3: Disease-enhanced
  # Top disease >= fc_enhanced * median of remaining
  if ((top_expression / (median_remaining + small_value)) >= 
      thresholds$fc_enhanced) {
    return(list(
      class = "enhanced", 
      label = names(expression_vector)[top_disease]
    ))
  }
  
  # Classification 4: Ubiquitous
  # Broadly expressed with low variability
  median_expression <- median(expression_vector)
  cv_expression <- calculate_cv(expression_vector)
  
  if (median_expression >= thresholds$ubiquitous_min_cpm && 
      cv_expression <= thresholds$ubiquitous_cv_max) {
    return(list(class = "ubiquitous", label = NA_character_))
  }
  
  # Classification 5: Mixed (default)
  list(class = "mixed", label = NA_character_)
}

#' Classify specificity patterns for all genes
#' @param disease_means Matrix of disease mean expressions
#' @param config Configuration list with thresholds
#' @return Data frame with gene classifications
classify_all_genes <- function(disease_means, config) {
  cat("[INFO] Classifying gene specificity patterns...\n")
  
  # Apply classification to each gene (row)
  classification_results <- apply(disease_means, 1, function(gene_expression) {
    classify_gene_specificity(gene_expression, config)
  })
  
  # Extract classification results
  specificity_classes <- vapply(classification_results, 
                               function(result) result$class, 
                               character(1))
  
  specificity_labels <- vapply(classification_results, 
                              function(result) result$label, 
                              character(1))
  
  # Create results data frame
  results_df <- tibble(
    gene = rownames(disease_means),
    class = specificity_classes,
    label = specificity_labels,
    top_disease = colnames(disease_means)[max.col(disease_means, 
                                                 ties.method = "first")],
    max_expression_cpm = rowMaxs(disease_means)
  )
  
  results_df
}

# =============================================================================
# CORRELATION ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate disease-disease correlations
#' @param disease_means Matrix of disease mean expressions
#' @return Correlation matrix
calculate_disease_correlations <- function(disease_means) {
  cat("[INFO] Calculating disease-disease correlations...\n")
  
  correlation_matrix <- cor(disease_means, method = "spearman", 
                           use = "pairwise.complete.obs")
  
  correlation_matrix
}

#' Create correlation heatmap
#' @param correlation_matrix Disease correlation matrix
#' @param output_file Output PDF file path
create_correlation_heatmap <- function(correlation_matrix, output_file) {
  cat("[INFO] Creating disease correlation heatmap...\n")
  
  # Convert correlation to distance for clustering
  # Distance = 1 - correlation, clamped to [0, 2]
  correlation_distance <- as.dist(pmax(0, 1 - correlation_matrix))
  
  # Create heatmap
  pdf(output_file, width = 8, height = 7)
  pheatmap(
    correlation_matrix,
    clustering_distance_rows = correlation_distance,
    clustering_distance_cols = correlation_distance,
    clustering_method = "complete",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    main = "Disease-Disease Spearman Correlation (Mean CPM)",
    fontsize = 10
  )
  dev.off()
  
  cat(sprintf("[INFO] Heatmap saved to: %s\n", output_file))
}

# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

#' Save analysis results to files
#' @param disease_means Matrix of disease mean expressions
#' @param specificity_results Data frame with gene classifications
#' @param correlation_matrix Disease correlation matrix
#' @param config Configuration list
save_analysis_results <- function(disease_means, specificity_results, 
                                 correlation_matrix, config) {
  cat("[INFO] Saving analysis results...\n")
  
  outdir <- config$outdir
  
  # Save disease means
  disease_means_df <- as.data.frame(disease_means) %>% 
    rownames_to_column("gene")
  write.csv(disease_means_df, 
           file.path(outdir, "disease_means_cpm.csv"), 
           row.names = FALSE)
  
  # Save gene specificity classifications
  write.csv(specificity_results, 
           file.path(outdir, "gene_specificity_calls.csv"), 
           row.names = FALSE)
  
  # Save classification summary
  classification_counts <- specificity_results %>% 
    count(class, name = "n_genes") %>%
    arrange(desc(n_genes))
  
  write.csv(classification_counts, 
           file.path(outdir, "gene_specificity_counts.csv"), 
           row.names = FALSE)
  
  # Save disease correlations
  correlation_df <- as.data.frame(correlation_matrix) %>% 
    rownames_to_column("disease")
  write.csv(correlation_df, 
           file.path(outdir, "disease_spearman_cor.csv"), 
           row.names = FALSE)
  
  # Print classification summary
  cat("[INFO] Gene specificity classification summary:\n")
  print(classification_counts)
  
  cat(sprintf("[INFO] Results saved to: %s\n", outdir))
}

#' Save diagnostic information
#' @param disease_means Matrix of disease mean expressions
#' @param config Configuration list
save_diagnostic_info <- function(disease_means, config) {
  cat("[INFO] Saving diagnostic information...\n")
  
  outdir <- config$outdir
  
  # Save sample of gene and disease names for debugging
  writeLines(utils::capture.output(head(rownames(disease_means))), 
            file.path(outdir, "head_genes.txt"))
  
  writeLines(utils::capture.output(head(colnames(disease_means))), 
            file.path(outdir, "head_diseases.txt"))
  
  # Save session information for reproducibility
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

#' Main analysis function
main_analysis <- function() {
  cat("[INFO] Starting CLD gene specificity analysis...\n")
  
  # 1. Load and prepare data
  dsmz_data <- load_dsmz_data(config$dsmz_rds, config$dsmz_meta_csv)
  counts_matrix <- build_count_matrix(dsmz_data$raw)
  
  # 2. Align metadata with count matrix
  aligned_data <- align_metadata_counts(dsmz_data$meta, counts_matrix)
  
  # 3. Define disease labels and filter by group size
  meta_with_disease <- define_disease_labels(aligned_data$meta, 
                                           config$disease_col)
  
  filtered_data <- filter_by_disease_size(meta_with_disease, 
                                         aligned_data$counts,
                                         config$min_samples_per_disease)
  
  # 4. Normalize expression data
  cpm_matrix <- normalize_to_cpm(filtered_data$counts)
  
  # 5. Calculate disease means and filter genes
  disease_means <- calculate_disease_means(cpm_matrix, filtered_data$meta)
  
  expressed_disease_means <- filter_expressed_genes(disease_means, 
                                                   config$min_cpm_gene)
  
  # 6. Classify gene specificity patterns
  specificity_results <- classify_all_genes(expressed_disease_means, config)
  
  # 7. Analyze disease-disease relationships
  correlation_matrix <- calculate_disease_correlations(expressed_disease_means)
  
  # 8. Create visualizations
  heatmap_file <- file.path(config$outdir, "disease_correlation_heatmap.pdf")
  create_correlation_heatmap(correlation_matrix, heatmap_file)
  
  # 9. Save results
  save_analysis_results(expressed_disease_means, specificity_results, 
                       correlation_matrix, config)
  
  save_diagnostic_info(expressed_disease_means, config)
  
  cat(sprintf("[INFO] Analysis completed successfully!\n"))
  cat(sprintf("[INFO] All outputs saved to: %s\n", config$outdir))
}

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Run analysis if script is executed (not sourced)
if (!interactive()) {
  main_analysis()
}