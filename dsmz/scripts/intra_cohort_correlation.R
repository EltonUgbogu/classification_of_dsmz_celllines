#!/usr/bin/env Rscript

# ===============================================================
# Intra-Cohort Correlation Analysis Pipeline (Refactored)
# 
# This pipeline performs comprehensive correlation analysis between
# DSMZ cell lines and TCGA samples, focusing on:
# 1. Within-group sample correlations (intra-cohort consistency)
# 2. Centroid-based correlations between groups
# 3. Cross-platform comparison of internal consistency
# 4. Organ-to-TCGA cohort mapping analysis
# ===============================================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(matrixStats)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
})

# Custom operators and utilities
`%||%` <- function(x, y) if (is.null(x)) y else x
set.seed(42)

# Configuration parameters
CONFIG <- list(
  # File paths
  tcga_se_rds   = Sys.getenv("TCGA_SE_RDS", "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds"),
  dsmz_rds      = Sys.getenv("DSMZ_RDS", "/home/chu25/data/dsmz/DSMZ_count_gene.rds"),
  dsmz_meta_csv = Sys.getenv("DSMZ_META_CSV", "/home/chu25/data/dsmz/DSMZ_metadata_with_sample_name.csv"),
  outdir        = Sys.getenv("OUTDIR", "/home/chu25/dsmz/results/intra_cohort_corr"),
  
  # Analysis parameters
  var_gene_set  = Sys.getenv("VAR_GENE_SET", "all"),  # "tcga_5k", "tcga_10k", or "all"
  min_shared    = as.numeric(Sys.getenv("MIN_SHARED", "1000")),
  min_group_n   = as.numeric(Sys.getenv("MIN_GROUP_N", "3")),  # min samples per group
  correlation_method = Sys.getenv("CORR_METHOD", "spearman")
)

# Logging function
logmsg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(...)))
}

# ==============================================================================
# ORGAN TO TCGA COHORT MAPPING
# ==============================================================================

#' Comprehensive mapping between organ systems and TCGA cohorts
ORGAN_TO_TCGA_MAPPING <- list(
  "Adrenal"       = c("ACC", "PCPG"),
  "Bladder"       = "BLCA",
  "Brain/CNS"     = c("GBM", "LGG"),
  "Breast"        = "BRCA",
  "Cervix"        = "CESC",
  "Liver/Biliary" = c("LIHC", "CHOL"),
  "Colon/Rectum"  = c("COAD", "READ"),
  "Esophagus"     = "ESCA",
  "Head/Neck"     = "HNSC",
  "Kidney"        = c("KICH", "KIRC", "KIRP"),
  "Lung"          = c("LUAD", "LUSC"),
  "Mesothelium"   = "MESO",
  "Ovary"         = "OV",
  "Pancreas"      = "PAAD",
  "Prostate"      = "PRAD",
  "Sarcoma"       = "SARC",
  "Skin"          = "SKCM",
  "Stomach"       = "STAD",
  "Testis"        = "TGCT",
  "Thyroid"       = "THCA",
  "Thymus"        = "THYM",
  "Uterus"        = c("UCEC", "UCS"),
  "Uveal/Eye"     = "UVM",
  "Lymphoid"      = "DLBC",
  "Leukemia"      = "LAML",
  "Myeloid"       = "LAML",
  "Bone Marrow"   = "LAML",
  "Body Cavity"   = "MESO",
  "Other"         = character(0),
  "Unknown"       = character(0)
)

# ==============================================================================
# DATA LOADING AND PREPROCESSING FUNCTIONS
# ==============================================================================

#' Load and preprocess TCGA data
#' @param filepath Path to TCGA SummarizedExperiment RDS file
#' @return List containing SummarizedExperiment object, count matrix, and sample labels
load_tcga_data <- function(filepath) {
  logmsg("Loading TCGA data from: %s", filepath)
  
  tcga_se <- readRDS(filepath)
  count_matrix <- assay(tcga_se)
  
  # Handle Ensembl ID versioning
  gene_ids_no_version <- sub("\\..*$", "", rownames(count_matrix))
  
  if (any(duplicated(gene_ids_no_version))) {
    logmsg("Collapsing %d duplicated gene IDs after version removal", 
           sum(duplicated(gene_ids_no_version)))
    count_matrix <- rowsum(count_matrix, gene_ids_no_version, reorder = TRUE)
    rownames(count_matrix) <- sort(unique(gene_ids_no_version))
  } else {
    rownames(count_matrix) <- gene_ids_no_version
  }
  
  # Extract cohort labels
  sample_labels <- extract_tcga_cohort_labels(tcga_se)
  
  # Ensure data types
  storage.mode(count_matrix) <- "double"
  
  list(
    se = tcga_se,
    counts = count_matrix,
    labels = sample_labels
  )
}

#' Extract TCGA cohort labels from SummarizedExperiment
#' @param tcga_se TCGA SummarizedExperiment object
#' @return Named character vector of cohort labels
extract_tcga_cohort_labels <- function(tcga_se) {
  # Try different possible column names for cohort information
  possible_columns <- c("project_id", "study", "disease_type")
  
  cohort_labels <- NULL
  for (col_name in possible_columns) {
    if (col_name %in% colnames(colData(tcga_se))) {
      cohort_labels <- colData(tcga_se)[[col_name]]
      logmsg("Using '%s' column for TCGA cohort labels", col_name)
      break
    }
  }
  
  if (is.null(cohort_labels)) {
    stop("No suitable cohort label column found in TCGA data")
  }
  
  # Clean and standardize labels
  cohort_labels <- as.character(cohort_labels)
  cohort_labels <- sub("^TCGA-", "", cohort_labels)
  names(cohort_labels) <- colnames(assay(tcga_se))
  
  return(cohort_labels)
}

#' Load and preprocess DSMZ data
#' @param count_filepath Path to DSMZ count RDS file
#' @param metadata_filepath Path to DSMZ metadata CSV file
#' @return List containing raw data, count matrix, and annotated metadata
load_dsmz_data <- function(count_filepath, metadata_filepath) {
  logmsg("Loading DSMZ data...")
  
  # Load count data
  raw_data <- readRDS(count_filepath)
  metadata <- read.csv(metadata_filepath)
  
  # Validate required columns
  required_count_cols <- c("Ensembl_ID", "gene_name")
  required_meta_cols <- c("sample_name")
  
  if (!all(required_count_cols %in% colnames(raw_data))) {
    stop("Missing required columns in DSMZ count data: ", 
         paste(setdiff(required_count_cols, colnames(raw_data)), collapse = ", "))
  }
  
  if (!all(required_meta_cols %in% names(metadata))) {
    stop("Missing required columns in DSMZ metadata: ", 
         paste(setdiff(required_meta_cols, names(metadata)), collapse = ", "))
  }
  
  # Process count matrix
  count_matrix <- build_dsmz_count_matrix(raw_data)
  
  # Process and align metadata
  processed_metadata <- process_dsmz_metadata(metadata, colnames(count_matrix))
  
  # Filter count matrix to match metadata
  count_matrix <- count_matrix[, processed_metadata$sample_id, drop = FALSE]
  
  # Ensure alignment
  stopifnot(identical(colnames(count_matrix), processed_metadata$sample_id))
  
  # Ensure data types
  storage.mode(count_matrix) <- "double"
  
  list(
    raw = raw_data,
    counts = count_matrix,
    metadata = processed_metadata
  )
}

#' Build DSMZ count matrix from raw data
#' @param raw_data Raw DSMZ data with gene annotations and sample columns
#' @return Processed count matrix
build_dsmz_count_matrix <- function(raw_data) {
  # Identify annotation and sample columns
  annotation_columns <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_columns <- setdiff(colnames(raw_data), annotation_columns)
  
  # Convert sample columns to numeric
  raw_data[sample_columns] <- lapply(raw_data[sample_columns], function(x) {
    as.numeric(as.character(x))
  })
  
  # Create count matrix
  count_matrix <- as.matrix(raw_data[, sample_columns, drop = FALSE])
  rownames(count_matrix) <- raw_data$Ensembl_ID
  
  # Handle duplicate gene IDs
  if (any(duplicated(rownames(count_matrix)))) {
    logmsg("Collapsing %d duplicated gene IDs in DSMZ data", 
           sum(duplicated(rownames(count_matrix))))
    count_matrix <- rowsum(count_matrix, rownames(count_matrix), reorder = TRUE)
  }
  
  # Remove columns with all NA values
  valid_columns <- vapply(as.data.frame(count_matrix), 
                         function(col) any(!is.na(col)), 
                         logical(1))
  
  if (!all(valid_columns)) {
    n_removed <- sum(!valid_columns)
    logmsg("Removing %d DSMZ samples with all NA values", n_removed)
    count_matrix <- count_matrix[, valid_columns, drop = FALSE]
  }
  
  return(count_matrix)
}

#' Process DSMZ metadata and create organ annotations
#' @param metadata Raw DSMZ metadata
#' @param available_samples Available sample IDs from count matrix
#' @return Processed metadata with organ annotations
process_dsmz_metadata <- function(metadata, available_samples) {
  processed_metadata <- metadata %>%
    mutate(sample_id = sample_name) %>%
    filter(sample_id %in% available_samples) %>%
    mutate(organ = classify_organ_from_origin(Origin))
  
  return(processed_metadata)
}

#' Classify organs from DSMZ origin descriptions
#' @param origin_descriptions Vector of origin descriptions
#' @return Vector of standardized organ classifications
classify_organ_from_origin <- function(origin_descriptions) {
  organ_classifications <- case_when(
    str_detect(origin_descriptions, regex("breast", TRUE)) ~ "Breast",
    str_detect(origin_descriptions, regex("lung", TRUE)) ~ "Lung",
    str_detect(origin_descriptions, regex("liver|hepato|bile", TRUE)) ~ "Liver/Biliary",
    str_detect(origin_descriptions, regex("adrenal", TRUE)) ~ "Adrenal",
    str_detect(origin_descriptions, regex("brain|gli", TRUE)) ~ "Brain/CNS",
    str_detect(origin_descriptions, regex("lymph|mantle|hodgkin|effusion|b[- ]?cell|t[- ]?cell|nk|natural\\s*killer", TRUE)) ~ "Lymphoid",
    str_detect(origin_descriptions, regex("leukemia", TRUE)) ~ "Leukemia",
    str_detect(origin_descriptions, regex("myeloid|megakaryocytic|erythroid|monocytic", TRUE)) ~ "Myeloid",
    str_detect(origin_descriptions, regex("bone\\s*marrow", TRUE)) ~ "Bone Marrow",
    str_detect(origin_descriptions, regex("abdomen|pleural", TRUE)) ~ "Body Cavity",
    is.na(origin_descriptions) | origin_descriptions == "NONE" ~ "Unknown",
    TRUE ~ "Other"
  )
  
  return(organ_classifications)
}

# ==============================================================================
# GENE SELECTION FUNCTIONS
# ==============================================================================

#' Select variable genes based on TCGA cohort expression patterns
#' @param tcga_counts TCGA count matrix
#' @param tcga_labels TCGA cohort labels
#' @param selection_mode Gene selection strategy ("all", "tcga_5k", "tcga_10k")
#' @return List containing selected genes and TCGA cohort means
select_variable_genes <- function(tcga_counts, tcga_labels, selection_mode = "all") {
  # Calculate cohort-wise mean expression
  cohorts <- sort(unique(tcga_labels))
  cohort_means <- sapply(cohorts, function(cohort) {
    cohort_samples <- names(tcga_labels)[tcga_labels == cohort]
    rowMeans(tcga_counts[, cohort_samples, drop = FALSE])
  })
  colnames(cohort_means) <- cohorts
  
  # Select genes based on mode
  selected_genes <- rownames(cohort_means)
  selection_mode <- tolower(selection_mode)
  
  if (selection_mode %in% c("tcga_5k", "tcga_10k")) {
    # Calculate inter-quartile range across cohorts for each gene
    gene_iqr <- matrixStats::rowIQRs(cohort_means)
    names(gene_iqr) <- rownames(cohort_means)
    
    # Select top variable genes
    n_genes <- ifelse(selection_mode == "tcga_5k", 5000, 10000)
    n_genes <- min(n_genes, length(gene_iqr))
    
    selected_genes <- names(sort(gene_iqr, decreasing = TRUE)[1:n_genes])
    
    logmsg("Selected top %d variable genes based on TCGA cohort IQR", n_genes)
  } else {
    logmsg("Using all %d available genes", length(selected_genes))
  }
  
  list(
    genes = selected_genes,
    tcga_cohort_means = cohort_means
  )
}

# ==============================================================================
# CORRELATION ANALYSIS FUNCTIONS
# ==============================================================================

#' Calculate within-group correlations for all groups
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param group_labels Named vector of group labels for samples
#' @param min_samples Minimum samples per group to include
#' @param method Correlation method ("spearman", "pearson")
#' @return Tibble with group, correlation values in long format
calculate_within_group_correlations <- function(expression_matrix, group_labels, 
                                               min_samples = 3, method = "spearman") {
  # Validate inputs
  stopifnot(identical(names(group_labels), colnames(expression_matrix)))
  
  groups <- sort(unique(group_labels))
  correlation_results <- vector("list", length(groups))
  
  logmsg("Calculating within-group correlations for %d groups", length(groups))
  
  for (i in seq_along(groups)) {
    group <- groups[i]
    group_samples <- names(group_labels)[group_labels == group]
    
    if (length(group_samples) < min_samples) {
      logmsg("  Skipping group '%s' (only %d samples, minimum %d required)", 
             group, length(group_samples), min_samples)
      next
    }
    
    # Calculate correlation matrix for group
    group_expression <- expression_matrix[, group_samples, drop = FALSE]
    correlation_matrix <- suppressWarnings(
      cor(group_expression, method = method, use = "pairwise.complete.obs")
    )
    
    # Extract upper triangle (unique pairs)
    upper_triangle_correlations <- extract_upper_triangle(correlation_matrix)
    
    correlation_results[[i]] <- tibble(
      group = group,
      rho = upper_triangle_correlations
    )
  }
  
  # Combine results and filter finite values
  combined_results <- bind_rows(correlation_results) %>%
    filter(is.finite(rho))
  
  return(combined_results)
}

#' Extract upper triangle values from a correlation matrix
#' @param correlation_matrix Square correlation matrix
#' @return Vector of upper triangle values (excluding diagonal)
extract_upper_triangle <- function(correlation_matrix) {
  if (ncol(correlation_matrix) < 2) {
    return(numeric(0))
  }
  return(correlation_matrix[upper.tri(correlation_matrix, diag = FALSE)])
}

#' Calculate group centroids (mean expression profiles)
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param group_labels Named vector of group labels for samples
#' @return Matrix of group centroids (genes x groups)
calculate_group_centroids <- function(expression_matrix, group_labels) {
  stopifnot(identical(names(group_labels), colnames(expression_matrix)))
  
  groups <- sort(unique(group_labels))
  centroids <- sapply(groups, function(group) {
    group_samples <- names(group_labels)[group_labels == group]
    rowMeans(expression_matrix[, group_samples, drop = FALSE])
  })
  
  colnames(centroids) <- groups
  return(centroids)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create violin plot of within-group correlations
#' @param correlation_data Tibble with group and rho columns
#' @param title Plot title
#' @param x_label X-axis label
#' @param output_path Output file path
#' @param order_by_median Whether to order groups by median correlation
create_correlation_violin_plot <- function(correlation_data, title, x_label, 
                                         output_path, order_by_median = TRUE) {
  if (order_by_median) {
    # Order groups by median correlation
    group_medians <- correlation_data %>%
      group_by(group) %>%
      summarise(median_rho = median(rho, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(median_rho))
    
    correlation_data$group <- factor(correlation_data$group, levels = group_medians$group)
    
    # Save median summary
    median_output_path <- sub("\\.pdf$", "_medians.csv", output_path)
    write.csv(group_medians, median_output_path, row.names = FALSE)
  }
  
  # Create plot
  p <- ggplot(correlation_data, aes(group, rho)) +
    geom_violin(fill = "grey85", alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(x = x_label, 
         y = "Within-group sample–sample Spearman ρ",
         title = title) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  # Adjust text size for many groups
  if (nlevels(correlation_data$group) > 25) {
    p <- p + theme(axis.text.x = element_text(size = 8))
  }
  
  # Save plot
  pdf(output_path, width = max(6, nlevels(correlation_data$group) * 0.3), height = 5.5)
  print(p)
  dev.off()
  
  logmsg("Saved correlation plot: %s", output_path)
}

#' Create centroid correlation heatmap
#' @param centroids Matrix of group centroids
#' @param title Plot title
#' @param output_path Output file path
#' @param correlation_method Correlation method for centroids
create_centroid_heatmap <- function(centroids, title, output_path, 
                                   correlation_method = "spearman") {
  # Calculate centroid correlations
  centroid_correlations <- suppressWarnings(
    cor(centroids, method = correlation_method, use = "pairwise.complete.obs")
  )
  
  # Create heatmap
  pdf(output_path, width = max(7, ncol(centroids) * 0.4), 
      height = max(6, nrow(centroids) * 0.4))
  
  pheatmap(
    centroid_correlations,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = title,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
    fontsize = max(8, min(12, 200 / ncol(centroids)))
  )
  
  dev.off()
  
  logmsg("Saved centroid heatmap: %s", output_path)
  return(centroid_correlations)
}

#' Create DSMZ vs TCGA comparison scatter plot
#' @param comparison_data Data frame with comparison results
#' @param output_path Output file path
create_comparison_scatter_plot <- function(comparison_data, output_path) {
  if (nrow(comparison_data) == 0) {
    logmsg("No data available for DSMZ vs TCGA comparison scatter plot")
    return()
  }
  
  p <- ggplot(comparison_data, aes(tcga_median_rho, dsmz_median_rho, label = organ)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60", size = 0.8) +
    geom_point(size = 2.5, alpha = 0.8, color = "steelblue") +
    geom_text_repel(size = 3, max.overlaps = 50, 
                    box.padding = 0.3, point.padding = 0.2) +
    theme_bw(base_size = 11) +
    labs(
      x = "TCGA median within-cohort ρ",
      y = "DSMZ median within-organ ρ",
      title = "Internal consistency: DSMZ vs TCGA\n(mapped by organ → cohort)",
      subtitle = "Diagonal line represents perfect agreement"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Add correlation coefficient if enough points
  if (nrow(comparison_data) >= 3) {
    correlation_coef <- cor(comparison_data$tcga_median_rho, 
                           comparison_data$dsmz_median_rho, 
                           use = "complete.obs")
    
    p <- p + annotate("text", 
                     x = min(comparison_data$tcga_median_rho, na.rm = TRUE) + 0.02,
                     y = max(comparison_data$dsmz_median_rho, na.rm = TRUE) - 0.02,
                     label = sprintf("r = %.3f", correlation_coef),
                     hjust = 0, vjust = 1, size = 4, color = "darkred")
  }
  
  pdf(output_path, width = 7, height = 6)
  print(p)
  dev.off()
  
  logmsg("Saved comparison scatter plot: %s", output_path)
}

# ==============================================================================
# ORGAN-TCGA MAPPING FUNCTIONS
# ==============================================================================

#' Create organ to TCGA cohort mapping table
#' @param dsmz_medians Data frame with DSMZ organ medians
#' @param tcga_medians Data frame with TCGA cohort medians
#' @return Data frame with mapped comparisons
create_organ_tcga_mapping <- function(dsmz_medians, tcga_medians) {
  # Expand organ-to-TCGA mapping into pairs
  mapping_pairs <- map_dfr(names(ORGAN_TO_TCGA_MAPPING), function(organ) {
    tcga_cohorts <- ORGAN_TO_TCGA_MAPPING[[organ]]
    if (length(tcga_cohorts) == 0) {
      return(tibble())
    }
    tibble(organ = organ, cohort = tcga_cohorts)
  })
  
  # Filter to available data
  available_pairs <- mapping_pairs %>%
    filter(
      cohort %in% tcga_medians[[1]], 
      organ %in% dsmz_medians[[1]]
    )
  
  # Join with median data
  comparison_data <- available_pairs %>%
    left_join(dsmz_medians, by = setNames(names(dsmz_medians)[1], "organ")) %>%
    rename(dsmz_median_rho = !!names(dsmz_medians)[2]) %>%
    left_join(tcga_medians, by = setNames(names(tcga_medians)[1], "cohort")) %>%
    rename(tcga_median_rho = !!names(tcga_medians)[2])
  
  return(comparison_data)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTIONS
# ==============================================================================

#' Perform DSMZ within-organ correlation analysis
#' @param dsmz_data DSMZ data list
#' @param selected_genes Vector of genes to use
#' @param config Configuration list
#' @param output_dir Output directory
#' @return List with correlation results and centroids
analyze_dsmz_correlations <- function(dsmz_data, selected_genes, config, output_dir) {
  logmsg("Analyzing DSMZ within-organ correlations...")
  
  # Create sample labels
  dsmz_labels <- setNames(dsmz_data$metadata$organ, dsmz_data$metadata$sample_id)
  dsmz_labels <- dsmz_labels[colnames(dsmz_data$counts)]
  
  # Calculate within-organ correlations
  dsmz_correlations <- calculate_within_group_correlations(
    expression_matrix = dsmz_data$counts[selected_genes, , drop = FALSE],
    group_labels = dsmz_labels,
    min_samples = config$min_group_n,
    method = config$correlation_method
  )
  
  # Create violin plot
  create_correlation_violin_plot(
    correlation_data = dsmz_correlations,
    title = "DSMZ internal consistency by organ",
    x_label = "DSMZ organ",
    output_path = file.path(output_dir, "dsmz_within_group_correlations.pdf")
  )
  
  # Calculate and visualize centroids
  dsmz_centroids <- calculate_group_centroids(
    dsmz_data$counts[selected_genes, , drop = FALSE], 
    dsmz_labels
  )
  
  create_centroid_heatmap(
    centroids = dsmz_centroids,
    title = "DSMZ organ centroids: Spearman correlation",
    output_path = file.path(output_dir, "dsmz_centroid_correlation_heatmap.pdf")
  )
  
  return(list(
    correlations = dsmz_correlations,
    centroids = dsmz_centroids
  ))
}

#' Perform TCGA within-cohort correlation analysis
#' @param tcga_data TCGA data list
#' @param selected_genes Vector of genes to use
#' @param config Configuration list
#' @param output_dir Output directory
#' @return List with correlation results and centroids
analyze_tcga_correlations <- function(tcga_data, selected_genes, config, output_dir) {
  logmsg("Analyzing TCGA within-cohort correlations...")
  
  # Calculate within-cohort correlations
  tcga_correlations <- calculate_within_group_correlations(
    expression_matrix = tcga_data$counts[selected_genes, , drop = FALSE],
    group_labels = tcga_data$labels,
    min_samples = config$min_group_n,
    method = config$correlation_method
  )
  
  # Create violin plot
  create_correlation_violin_plot(
    correlation_data = tcga_correlations,
    title = "TCGA internal consistency by cohort",
    x_label = "TCGA cohort",
    output_path = file.path(output_dir, "tcga_within_group_correlations.pdf")
  )
  
  # Use pre-calculated cohort means as centroids
  tcga_centroids <- select_variable_genes(tcga_data$counts, tcga_data$labels, "all")$tcga_cohort_means
  tcga_centroids <- tcga_centroids[selected_genes, , drop = FALSE]
  
  create_centroid_heatmap(
    centroids = tcga_centroids,
    title = "TCGA cohort centroids: Spearman correlation", 
    output_path = file.path(output_dir, "tcga_centroid_correlation_heatmap.pdf")
  )
  
  return(list(
    correlations = tcga_correlations,
    centroids = tcga_centroids
  ))
}

#' Compare DSMZ and TCGA median correlations
#' @param dsmz_results DSMZ analysis results
#' @param tcga_results TCGA analysis results
#' @param output_dir Output directory
#' @return Comparison data frame
compare_dsmz_tcga_medians <- function(dsmz_results, tcga_results, output_dir) {
  logmsg("Comparing DSMZ vs TCGA median within-group correlations...")
  
  # Calculate median correlations
  dsmz_medians <- dsmz_results$correlations %>%
    group_by(group) %>%
    summarise(median_within_rho = median(rho, na.rm = TRUE), .groups = 'drop') %>%
    rename(organ = group)
  
  tcga_medians <- tcga_results$correlations %>%
    group_by(group) %>%
    summarise(median_within_rho = median(rho, na.rm = TRUE), .groups = 'drop') %>%
    rename(cohort = group)
  
  # Save individual median files
  write.csv(dsmz_medians, file.path(output_dir, "dsmz_within_median_by_organ.csv"), row.names = FALSE)
  write.csv(tcga_medians, file.path(output_dir, "tcga_within_median_by_cohort.csv"), row.names = FALSE)
  
  # Create organ-to-TCGA mapping
  comparison_data <- create_organ_tcga_mapping(dsmz_medians, tcga_medians)
  
  # Save comparison results
  write.csv(comparison_data, 
           file.path(output_dir, "dsmz_vs_tcga_median_within_comparison.csv"), 
           row.names = FALSE)
  
  # Create scatter plot
  create_comparison_scatter_plot(
    comparison_data = comparison_data,
    output_path = file.path(output_dir, "dsmz_vs_tcga_median_within_scatter.pdf")
  )
  
  return(comparison_data)
}

# ==============================================================================
# SUMMARY AND REPORTING FUNCTIONS
# ==============================================================================

#' Generate comprehensive analysis summary
#' @param config Configuration parameters
#' @param gene_info Gene selection information
#' @param dsmz_results DSMZ analysis results
#' @param tcga_results TCGA analysis results
#' @param comparison_data Cross-platform comparison data
#' @param output_dir Output directory
generate_analysis_summary <- function(config, gene_info, dsmz_results, tcga_results, 
                                     comparison_data, output_dir) {
  # Extract group information
  dsmz_groups <- sort(unique(dsmz_results$correlations$group))
  tcga_groups <- sort(unique(tcga_results$correlations$group))
  
  # Calculate summary statistics
  dsmz_overall_median <- median(dsmz_results$correlations$rho, na.rm = TRUE)
  tcga_overall_median <- median(tcga_results$correlations$rho, na.rm = TRUE)
  
  # Create summary lines
  summary_lines <- c(
    "Intra-Cohort Correlation Analysis Summary",
    "==========================================",
    "",
    sprintf("Analysis Date: %s", Sys.time()),
    sprintf("Gene Selection: %s (%d genes)", config$var_gene_set, length(gene_info$genes)),
    sprintf("Correlation Method: %s", config$correlation_method),
    sprintf("Minimum Group Size: %d samples", config$min_group_n),
    "",
    "DSMZ Analysis:",
    sprintf("  - Organs analyzed (n≥%d): %d", config$min_group_n, length(dsmz_groups)),
    sprintf("  - Organ list: %s", paste(dsmz_groups, collapse = ", ")),
    sprintf("  - Overall median within-organ correlation: %.3f", dsmz_overall_median),
    sprintf("  - Total sample pairs analyzed: %d", nrow(dsmz_results$correlations)),
    "",
    "TCGA Analysis:",
    sprintf("  - Cohorts analyzed (n≥%d): %d", config$min_group_n, length(tcga_groups)),
    sprintf("  - Cohort list: %s", paste(tcga_groups, collapse = ", ")),
    sprintf("  - Overall median within-cohort correlation: %.3f", tcga_overall_median),
    sprintf("  - Total sample pairs analyzed: %d", nrow(tcga_results$correlations)),
    "",
    "Cross-Platform Comparison:",
    sprintf("  - Organ-cohort mappings: %d", nrow(comparison_data)),
    "",
    "Output Files Generated:",
    "  Main Analyses:",
    "    - dsmz_within_group_correlations.pdf",
    "    - tcga_within_group_correlations.pdf", 
    "    - dsmz_centroid_correlation_heatmap.pdf",
    "    - tcga_centroid_correlation_heatmap.pdf",
    "    - dsmz_vs_tcga_median_within_scatter.pdf",
    "  Data Files:",
    "    - dsmz_within_median_by_organ.csv",
    "    - tcga_within_median_by_cohort.csv",
    "    - dsmz_vs_tcga_median_within_comparison.csv",
    "    - SUMMARY.txt (this file)"
  )
  
  # Add top/bottom performing groups
  if (nrow(dsmz_results$correlations) > 0) {
    dsmz_group_medians <- dsmz_results$correlations %>%
      group_by(group) %>%
      summarise(median_rho = median(rho, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(median_rho))
    
    summary_lines <- c(summary_lines,
      "",
      "Top 5 DSMZ Organs (by median within-organ correlation):",
      sprintf("  %s", paste(head(sprintf("%-15s: %.3f", 
                                        dsmz_group_medians$group, 
                                        dsmz_group_medians$median_rho), 5), 
                            collapse = "\n  "))
    )
  }
  
  if (nrow(tcga_results$correlations) > 0) {
    tcga_group_medians <- tcga_results$correlations %>%
      group_by(group) %>%
      summarise(median_rho = median(rho, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(median_rho))
    
    summary_lines <- c(summary_lines,
      "",
      "Top 5 TCGA Cohorts (by median within-cohort correlation):",
      sprintf("  %s", paste(head(sprintf("%-15s: %.3f", 
                                        tcga_group_medians$group, 
                                        tcga_group_medians$median_rho), 5), 
                            collapse = "\n  "))
    )
  }
  
  # Write summary to console and file
  summary_text <- paste(summary_lines, collapse = "\n")
  cat(summary_text, "\n")
  writeLines(summary_lines, file.path(output_dir, "SUMMARY.txt"))
  
  logmsg("Analysis summary saved to: %s", file.path(output_dir, "SUMMARY.txt"))
}

#' Validate analysis inputs and configuration
#' @param config Configuration list
validate_analysis_inputs <- function(config) {
  # Check required files exist
  required_files <- c(config$tcga_se_rds, config$dsmz_rds, config$dsmz_meta_csv)
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
  }
  
  # Validate parameters
  if (!config$var_gene_set %in% c("all", "tcga_5k", "tcga_10k")) {
    stop("var_gene_set must be one of: 'all', 'tcga_5k', 'tcga_10k'")
  }
  
  if (!config$correlation_method %in% c("spearman", "pearson")) {
    stop("correlation_method must be 'spearman' or 'pearson'")
  }
  
  if (config$min_group_n < 2) {
    stop("min_group_n must be at least 2")
  }
  
  if (config$min_shared < 100) {
    warning("min_shared is very low (", config$min_shared, "), results may be unreliable")
  }
  
  logmsg("Input validation passed")
}

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

#' Main analysis pipeline
main_analysis <- function() {
  logmsg("Starting Intra-Cohort Correlation Analysis Pipeline")
  
  # Create output directory
  dir.create(CONFIG$outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Print configuration
  logmsg("Configuration:")
  logmsg("  Output directory: %s", CONFIG$outdir)
  logmsg("  Gene selection: %s", CONFIG$var_gene_set)
  logmsg("  Correlation method: %s", CONFIG$correlation_method)
  logmsg("  Minimum group size: %d", CONFIG$min_group_n)
  logmsg("  Minimum shared genes: %d", CONFIG$min_shared)
  
  # Validate inputs
  validate_analysis_inputs(CONFIG)
  
  # === LOAD DATA ===
  tcga_data <- load_tcga_data(CONFIG$tcga_se_rds)
  dsmz_data <- load_dsmz_data(CONFIG$dsmz_rds, CONFIG$dsmz_meta_csv)
  
  # === FIND COMMON GENES ===
  common_genes <- intersect(rownames(tcga_data$counts), rownames(dsmz_data$counts))
  logmsg("Shared genes between datasets: %d", length(common_genes))
  
  if (length(common_genes) < CONFIG$min_shared) {
    warning(sprintf("Only %d shared genes found (minimum %d). Check Ensembl versions & annotations.",
                   length(common_genes), CONFIG$min_shared))
  }
  
  # Filter both datasets to common genes
  tcga_data$counts <- tcga_data$counts[common_genes, , drop = FALSE]
  dsmz_data$counts <- dsmz_data$counts[common_genes, , drop = FALSE]
  
  # === SELECT VARIABLE GENES ===
  gene_selection <- select_variable_genes(
    tcga_data$counts, 
    tcga_data$labels, 
    CONFIG$var_gene_set
  )
  
  selected_genes <- intersect(gene_selection$genes, common_genes)
  logmsg("Using %d genes for correlation analysis", length(selected_genes))
  
  # === PERFORM ANALYSES ===
  
  # Analyze DSMZ within-organ correlations
  dsmz_results <- analyze_dsmz_correlations(
    dsmz_data, selected_genes, CONFIG, CONFIG$outdir
  )
  
  # Analyze TCGA within-cohort correlations  
  tcga_results <- analyze_tcga_correlations(
    tcga_data, selected_genes, CONFIG, CONFIG$outdir
  )
  
  # Compare DSMZ vs TCGA medians
  comparison_data <- compare_dsmz_tcga_medians(
    dsmz_results, tcga_results, CONFIG$outdir
  )
  
  # === GENERATE SUMMARY ===
  generate_analysis_summary(
    CONFIG, gene_selection, dsmz_results, tcga_results, 
    comparison_data, CONFIG$outdir
  )
  
  logmsg("Analysis completed successfully!")
  logmsg("Results saved to: %s", CONFIG$outdir)
  
  # Return results for potential further analysis
  invisible(list(
    tcga = tcga_results,
    dsmz = dsmz_results,
    comparison = comparison_data,
    genes_used = selected_genes
  ))
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Check if required packages are installed
check_required_packages <- function() {
  required_packages <- c(
    "tidyverse", "SummarizedExperiment", "matrixStats", 
    "pheatmap", "viridis", "ggrepel"
  )
  
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  
  if (length(missing_packages) > 0) {
    stop(sprintf("Missing required packages: %s\nInstall with: install.packages(c(%s))",
                paste(missing_packages, collapse = ", "),
                paste(sprintf("'%s'", missing_packages), collapse = ", ")))
  }
  
  logmsg("All required packages are available")
}

# ==============================================================================
# SCRIPT EXECUTION
# ==============================================================================

# Execute main analysis if script is run directly
if (!interactive()) {
  tryCatch({
    # Check dependencies
    check_required_packages()
    
    # Run main analysis
    main_analysis()
    
  }, error = function(e) {
    logmsg("FATAL ERROR: %s", e$message)
    traceback()
    quit(status = 1)
  })
}