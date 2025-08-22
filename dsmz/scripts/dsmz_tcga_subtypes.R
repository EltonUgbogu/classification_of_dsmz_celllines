#!/usr/bin/env Rscript

# ===============================================
# DSMZ–TCGA Subtype Pipeline (Refactored)
# 
# This pipeline:
# 1. Builds subtype templates from TCGA data
# 2. Evaluates templates using 80/20 train/test split
# 3. Filters template genes for cell-line relevance
# 4. Predicts DSMZ cell line subtypes (if accuracy >= 80%)
# 5. Optionally analyzes AJCC stage signatures
# ===============================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(TCGAbiolinks)
})

# Custom operators and utilities
`%||%` <- function(x, y) if (is.null(x)) y else x
set.seed(42)

# Configuration parameters
CONFIG <- list(
  # File paths
  tcga_se_rds   = Sys.getenv("TCGA_SE_RDS",  "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds"),
  dsmz_rds      = Sys.getenv("DSMZ_COUNT_RDS", "/home/chu25/data/dsmz/DSMZ_count_gene.rds"),
  dsmz_meta_csv = Sys.getenv("DSMZ_META_CSV", "/home/chu25/data/dsmz/DSMZ_metadata_with_sample_name.csv"),
  outdir        = Sys.getenv("OUTDIR", "/home/chu25/dsmz/results/subtypes"),
  
  # Analysis thresholds
  lfc_min                 = 1.0,    # Minimum log fold change for DE genes
  fdr_max                 = 0.01,   # Maximum FDR for DE genes
  abs_lfc_tumor_cl_max    = 2.0,    # Max |logFC| between tumor and cell lines
  dsmz_top50_min_samples  = 2,      # Min DSMZ samples for gene inclusion
  min_per_subtype         = 15,     # Min samples per subtype to consider
  accuracy_pass           = 0.80,   # Min accuracy to proceed with DSMZ prediction
  stage_lfc_min           = 2.0,    # Min logFC for stage signature genes
  
  # Flags
  run_stage = TRUE
)

# Logging function
logmsg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(...)))
}

# ==============================================================================
# DATA LOADING AND PREPROCESSING FUNCTIONS
# ==============================================================================

#' Load and preprocess TCGA data
load_tcga_data <- function(tcga_se_path) {
  logmsg("Loading TCGA SE: %s", tcga_se_path)
  tcga_se <- readRDS(tcga_se_path)
  tcga_counts <- assay(tcga_se)
  
  # Strip version numbers from Ensembl IDs and handle duplicates
  tcga_genes_novers <- sub("\\..*$", "", rownames(tcga_counts))
  if (any(duplicated(tcga_genes_novers))) {
    logmsg("Collapsing duplicated TCGA Ensembl IDs after version strip.")
    tcga_counts <- rowsum(tcga_counts, tcga_genes_novers, reorder = TRUE)
    rownames(tcga_counts) <- sort(unique(tcga_genes_novers))
  } else {
    rownames(tcga_counts) <- tcga_genes_novers
  }
  
  list(se = tcga_se, counts = tcga_counts)
}

#' Load and preprocess DSMZ data
load_dsmz_data <- function(dsmz_rds_path, dsmz_meta_path) {
  logmsg("Loading DSMZ counts & metadata.")
  
  # Load count data
  dsmz_raw <- readRDS(dsmz_rds_path)
  stopifnot(all(c("Ensembl_ID", "gene_name") %in% colnames(dsmz_raw)))
  
  # Load metadata
  dsmz_meta <- read.csv(dsmz_meta_path)
  stopifnot("sample_name" %in% names(dsmz_meta))
  
  # Build count matrix
  annot_cols <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot_cols)
  
  # Convert to numeric and create matrix
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  dsmz_counts <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(dsmz_counts) <- dsmz_raw$Ensembl_ID
  
  # Handle duplicated gene IDs
  if (any(duplicated(rownames(dsmz_counts)))) {
    dsmz_counts <- rowsum(dsmz_counts, rownames(dsmz_counts), reorder = TRUE)
  }
  
  # Remove columns with all NA values
  keep_col <- vapply(as.data.frame(dsmz_counts), function(v) any(!is.na(v)), logical(1))
  if (!all(keep_col)) {
    logmsg("Dropping %d DSMZ columns with all NA.", sum(!keep_col))
    dsmz_counts <- dsmz_counts[, keep_col, drop = FALSE]
  }
  
  # Align metadata to counts
  dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)
  matched_ids <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))
  dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched_ids)
  dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop = FALSE]
  
  stopifnot(identical(colnames(dsmz_counts), dsmz_meta$sample_id))
  
  list(counts = dsmz_counts, metadata = dsmz_meta)
}

#' Extract TCGA cohort labels
extract_tcga_labels <- function(tcga_se) {
  labs <- NULL
  for (nm in c("project_id", "study", "disease_type")) {
    if (nm %in% colnames(colData(tcga_se))) {
      labs <- colData(tcga_se)[[nm]]
      break
    }
  }
  stopifnot(!is.null(labs))
  labs <- sub("^TCGA-", "", as.character(labs))
  names(labs) <- colnames(assay(tcga_se))
  return(labs)
}

# ==============================================================================
# NORMALIZATION AND TRANSFORMATION FUNCTIONS
# ==============================================================================

#' Apply voom normalization with quantile normalization
voom_quantile_norm <- function(count_mat) {
  y <- DGEList(counts = count_mat)
  y <- calcNormFactors(y)  # TMM normalization
  v <- voom(y, plot = FALSE)  # Convert to logCPM with precision weights
  E <- v$E
  # Apply quantile normalization
  normalizeBetweenArrays(E, method = "quantile")
}

#' Z-score transformation by genes (rows)
z_score_by_gene <- function(mat) {
  Z <- t(scale(t(mat)))
  Z[is.na(Z)] <- 0
  Z
}

# ==============================================================================
# SUBTYPE ANALYSIS FUNCTIONS
# ==============================================================================

#' Get subtype information for a cohort
get_subtype_info <- function(cohort, subtab, tcga_counts, min_per_subtype) {
  if (is.null(subtab)) {
    return(list(success = FALSE, message = "Subtype table not available"))
  }
  
  # Filter subtype table for current cohort
  st_cohort <- subtab %>%
    dplyr::filter(cohort == paste0("TCGA-", cohort)) %>%
    dplyr::select(sample, Subtype_Selected)
  
  # Find common samples between TCGA counts and subtype data
  common_samples <- intersect(colnames(tcga_counts), st_cohort$sample)
  
  if (length(common_samples) < min_per_subtype * 2) {
    return(list(success = FALSE, message = "Too few samples with subtype information"))
  }
  
  st_cohort <- st_cohort %>% filter(sample %in% common_samples)
  subtype_labels <- setNames(st_cohort$Subtype_Selected, st_cohort$sample)
  
  # Keep only subtypes with sufficient samples
  subtype_counts <- sort(table(subtype_labels), decreasing = TRUE)
  keep_subtypes <- names(subtype_counts[subtype_counts >= min_per_subtype])
  
  if (length(keep_subtypes) < 2) {
    return(list(success = FALSE, message = "Insufficient samples per subtype"))
  }
  
  # Filter to keep only samples with retained subtypes
  subtype_labels <- subtype_labels[subtype_labels %in% keep_subtypes]
  filtered_counts <- tcga_counts[, names(subtype_labels), drop = FALSE]
  
  list(
    success = TRUE,
    counts = filtered_counts,
    labels = subtype_labels,
    subtypes = keep_subtypes
  )
}

#' Build subtype templates using differential expression analysis
build_subtype_templates <- function(train_counts, train_labels, lfc_min, fdr_max) {
  # Create design matrix
  design <- model.matrix(~ 0 + factor(train_labels))
  colnames(design) <- levels(factor(train_labels))
  
  # Fit linear model
  fit <- lmFit(train_counts, design)
  
  # Create contrasts (each subtype vs. average of others)
  subtypes <- levels(factor(train_labels))
  contrast_list <- lapply(subtypes, function(s) {
    others <- setdiff(subtypes, s)
    paste0(s, " - (", paste(others, collapse = " + "), ")/", length(others))
  })
  names(contrast_list) <- subtypes
  
  contrasts <- makeContrasts(contrasts = unlist(contrast_list), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrasts))
  
  # Extract up-regulated genes for each subtype
  templates <- lapply(subtypes, function(s) {
    tt <- topTable(fit2, coef = s, number = Inf, sort.by = "none")
    rownames(tt)[which(tt$logFC >= lfc_min & tt$adj.P.Val <= fdr_max)]
  })
  names(templates) <- subtypes
  
  return(templates)
}

#' Filter template genes for cell-line relevance
filter_templates_for_cell_lines <- function(templates, tcga_expr, dsmz_expr, 
                                           abs_lfc_max, min_samples) {
  # Filter 1: Remove genes with large expression difference between TCGA and DSMZ
  tcga_mean <- rowMeans(tcga_expr)
  dsmz_mean <- rowMeans(dsmz_expr, na.rm = TRUE)
  lfc_tumor_dsmz <- tcga_mean - dsmz_mean
  
  genes_to_drop <- names(which(abs(lfc_tumor_dsmz) > abs_lfc_max))
  
  # Filter 2: Keep genes expressed above median in sufficient DSMZ samples
  dsmz_medians <- apply(dsmz_expr, 1, median, na.rm = TRUE)
  above_median <- sweep(dsmz_expr, 1, dsmz_medians, FUN = ">")
  genes_to_keep <- names(which(rowSums(above_median, na.rm = TRUE) >= min_samples))
  
  # Apply filters to templates
  filtered_templates <- lapply(templates, function(genes) {
    setdiff(intersect(genes, genes_to_keep), genes_to_drop)
  })
  
  return(filtered_templates)
}

#' NTP-style prediction using template matching
ntp_predict <- function(expr_z, template_list) {
  subtypes <- names(template_list)
  score_matrix <- matrix(NA_real_, nrow = length(subtypes), ncol = ncol(expr_z),
                        dimnames = list(subtypes, colnames(expr_z)))
  
  for (subtype in subtypes) {
    template_genes <- intersect(template_list[[subtype]], rownames(expr_z))
    if (length(template_genes) == 0) next
    
    # Calculate mean z-score across template genes
    score_matrix[subtype, ] <- colMeans(expr_z[template_genes, , drop = FALSE])
  }
  
  predictions <- apply(score_matrix, 2, function(v) names(which.max(v)))
  scores <- apply(score_matrix, 2, max, na.rm = TRUE)
  
  list(pred = predictions, score = scores, score_mat = score_matrix)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create confusion matrix plot
plot_confusion_matrix <- function(true_labels, predicted_labels, cohort, output_path) {
  cm <- table(True = true_labels, Pred = predicted_labels)
  accuracy <- mean(predicted_labels == true_labels)
  
  pdf(output_path, width = 6, height = 5)
  pheatmap(cm, cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, 
           main = sprintf("%s NTP hold-out (acc=%.1f%%)", cohort, 100 * accuracy))
  dev.off()
  
  write.csv(as.data.frame(cm), 
           sub("\\.pdf$", ".csv", output_path), 
           row.names = FALSE)
  
  return(accuracy)
}

#' Plot template sizes
plot_template_sizes <- function(templates, cohort, output_path) {
  template_sizes <- data.frame(
    subtype = names(templates),
    n_genes = vapply(templates, length, integer(1))
  ) %>%
    arrange(desc(n_genes))
  
  write.csv(template_sizes, sub("\\.pdf$", ".csv", output_path), row.names = FALSE)
  
  p <- ggplot(template_sizes, aes(reorder(subtype, -n_genes), n_genes)) +
    geom_col(fill = "grey60") +
    geom_text(aes(label = n_genes), vjust = -0.2, size = 3) +
    theme_bw(base_size = 11) +
    labs(x = "Subtype", y = "# template genes",
         title = sprintf("%s: Subtype template sizes (post cell-line filters)", cohort)) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  ggsave(output_path, p, width = 7, height = 4)
  
  return(template_sizes)
}

#' Plot DSMZ prediction scores
plot_dsmz_scores <- function(predictions, cohort, output_path) {
  call_df <- tibble(
    sample_id = names(predictions$pred),
    subtype = unname(predictions$pred),
    score = unname(predictions$score)
  ) %>% 
    arrange(desc(score))
  
  write.csv(call_df, sub("\\.pdf$", ".csv", output_path), row.names = FALSE)
  
  p <- ggplot(call_df, aes(subtype, score)) +
    geom_violin(fill = "grey85") +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    theme_bw(base_size = 11) +
    labs(x = "Predicted subtype", y = "Template score (mean Z)",
         title = sprintf("%s: DSMZ NTP scores by predicted subtype", cohort)) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  ggsave(output_path, p, width = 7, height = 4)
  
  return(call_df)
}

# ==============================================================================
# STAGE ANALYSIS FUNCTIONS
# ==============================================================================

#' Analyze AJCC stage signatures
analyze_stage_signatures <- function(tcga_se, tcga_counts, cohort_indices, 
                                    dsmz_expr, stage_lfc_min, output_dir) {
  # Extract stage information
  stage_cols <- intersect(
    c("ajcc_pathologic_stage", "ajcc_pathologic_stage_system_edition"),
    colnames(colData(tcga_se))
  )
  
  if (length(stage_cols) == 0) {
    logmsg("  No AJCC stage column found; skipping stage analysis for cohort.")
    return(invisible())
  }
  
  stage_data <- as.character(colData(tcga_se)[cohort_indices, stage_cols[1]])
  names(stage_data) <- colnames(tcga_counts)[cohort_indices]
  
  # Simplify stage labels
  stage_simple <- toupper(gsub("[^IVX]", "", stage_data))
  stage_simple[stage_simple == ""] <- NA
  
  valid_stages <- !is.na(stage_simple)
  if (sum(valid_stages) < 40) {
    logmsg("  Too few staged samples; skipping stage analysis.")
    return(invisible())
  }
  
  # Filter to samples with valid stages
  staged_counts <- tcga_counts[, valid_stages, drop = FALSE]
  stage_labels <- stage_simple[valid_stages]
  
  # Normalize expression data
  staged_expr <- voom_quantile_norm(staged_counts)
  
  # Calculate stage mean profiles
  stages <- names(sort(table(stage_labels), decreasing = TRUE))
  stage_means <- sapply(stages, function(s) {
    rowMeans(staged_expr[, stage_labels == s, drop = FALSE])
  })
  colnames(stage_means) <- stages
  
  # Generate stage signatures
  stage_signatures <- generate_stage_signatures(stage_means, stage_lfc_min)
  
  # Analyze DSMZ vs stages
  analyze_dsmz_stage_correlation(dsmz_expr, stage_means, stage_signatures, 
                                stages, output_dir)
}

#' Generate stage-specific gene signatures
generate_stage_signatures <- function(stage_means, lfc_min) {
  stages <- colnames(stage_means)
  
  # For each stage, find genes with high expression vs other stages
  best_other_lfc <- apply(stage_means, 1, function(gene_expr) {
    sapply(seq_along(gene_expr), function(i) {
      gene_expr[i] - max(gene_expr[-i])
    })
  })
  
  # Extract signature genes for each stage
  signatures <- lapply(stages, function(s) {
    stage_idx <- which(stages == s)
    signature_genes <- rownames(stage_means)[which(best_other_lfc[stage_idx, ] >= lfc_min)]
    signature_genes
  })
  names(signatures) <- stages
  
  return(signatures)
}

#' Analyze correlation between DSMZ samples and stage profiles
analyze_dsmz_stage_correlation <- function(dsmz_expr, stage_means, stage_signatures, 
                                         stages, output_dir) {
  # Find common genes
  common_genes <- intersect(rownames(stage_means), rownames(dsmz_expr))
  
  # Calculate Spearman correlations
  correlation_matrix <- outer(colnames(dsmz_expr), colnames(stage_means), 
                             Vectorize(function(dsmz_sample, stage) {
    suppressWarnings(
      cor(dsmz_expr[common_genes, dsmz_sample], 
          stage_means[common_genes, stage], 
          method = "spearman")
    )
  }))
  
  rownames(correlation_matrix) <- colnames(dsmz_expr)
  colnames(correlation_matrix) <- colnames(stage_means)
  
  # Save correlation results
  write.csv(
    as.data.frame(correlation_matrix) %>% rownames_to_column("DSMZ_sample"),
    file.path(output_dir, "stage_spearman_dsmz_vs_stage.csv"),
    row.names = FALSE
  )
  
  # Create correlation heatmap
  pdf(file.path(output_dir, "stage_correlation_heatmap.pdf"), width = 7.5, height = 5.5)
  pheatmap(correlation_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
           main = "DSMZ vs AJCC stage (Spearman ρ)")
  dev.off()
  
  # Calculate enrichment scores and create prioritization
  calculate_stage_prioritization(dsmz_expr, stage_signatures, correlation_matrix, 
                               stages, output_dir)
}

#' Calculate stage prioritization for DSMZ samples
calculate_stage_prioritization <- function(dsmz_expr, stage_signatures, 
                                         correlation_matrix, stages, output_dir) {
  # Z-score DSMZ expression
  dsmz_z <- z_score_by_gene(dsmz_expr)
  
  # Calculate enrichment scores
  enrichment_scores <- sapply(stages, function(stage) {
    sig_genes <- intersect(stage_signatures[[stage]], rownames(dsmz_z))
    if (length(sig_genes) > 0) {
      colMeans(dsmz_z[sig_genes, , drop = FALSE])
    } else {
      rep(NA_real_, ncol(dsmz_z))
    }
  })
  
  enrichment_scores <- as.data.frame(enrichment_scores)
  rownames(enrichment_scores) <- colnames(dsmz_expr)
  
  # Save enrichment scores
  write.csv(
    enrichment_scores %>% rownames_to_column("DSMZ_sample"),
    file.path(output_dir, "stage_signature_scores_dsmz.csv"),
    row.names = FALSE
  )
  
  # Generate prioritization for each stage
  for (stage in stages) {
    create_stage_prioritization(correlation_matrix, enrichment_scores, stage, output_dir)
  }
}

#' Create prioritization ranking for a specific stage
create_stage_prioritization <- function(correlation_matrix, enrichment_scores, 
                                       stage, output_dir) {
  prioritization_df <- tibble(
    DSMZ_sample = rownames(correlation_matrix),
    rho = correlation_matrix[, stage],
    score = enrichment_scores[rownames(correlation_matrix), stage]
  ) %>%
    mutate(
      rank_rho = rank(-rho, ties.method = "average"),
      rank_score = rank(-score, ties.method = "average", na.last = "keep"),
      rank_mean = (rank_rho + rank_score) / 2
    ) %>%
    arrange(rank_mean)
  
  # Save prioritization results
  write.csv(
    prioritization_df,
    file.path(output_dir, sprintf("stage_%s_prioritisation.csv", stage)),
    row.names = FALSE
  )
  
  # Plot top 20 prioritized samples
  top_samples <- head(prioritization_df, 20)
  
  p <- ggplot(top_samples, aes(x = reorder(DSMZ_sample, rank_mean), y = rank_mean)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(x = "DSMZ sample", y = "Mean rank (lower = better)",
         title = sprintf("DSMZ prioritisation for stage %s", stage))
  
  ggsave(
    file.path(output_dir, sprintf("stage_%s_prioritisation_top20.pdf", stage)),
    p, width = 7, height = 5
  )
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Process a single cohort
process_cohort <- function(cohort, tcga_data, dsmz_data, subtab, config) {
  logmsg("=== Processing Cohort: %s ===", cohort)
  
  # Create output directory
  output_dir <- file.path(config$outdir, cohort)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Extract cohort-specific TCGA data
  cohort_labels <- extract_tcga_labels(tcga_data$se)
  cohort_indices <- which(cohort_labels == cohort)
  
  if (length(cohort_indices) < 40) {
    logmsg("  Skipping %s (too few samples: %d).", cohort, length(cohort_indices))
    return(invisible())
  }
  
  tcga_cohort_counts <- tcga_data$counts[, cohort_indices, drop = FALSE]
  
  # Find common genes between TCGA and DSMZ
  common_genes <- intersect(rownames(tcga_cohort_counts), rownames(dsmz_data$counts))
  logmsg("  Shared genes: %d", length(common_genes))
  
  tcga_cohort_counts <- tcga_cohort_counts[common_genes, , drop = FALSE]
  dsmz_cohort_counts <- dsmz_data$counts[common_genes, , drop = FALSE]
  
  # Normalize expression data
  tcga_lcpm <- voom_quantile_norm(tcga_cohort_counts)
  dsmz_lcpm <- voom_quantile_norm(dsmz_cohort_counts)
  
  # === SUBTYPE ANALYSIS ===
  subtype_info <- get_subtype_info(cohort, subtab, tcga_cohort_counts, config$min_per_subtype)
  
  if (subtype_info$success) {
    process_subtype_analysis(subtype_info, tcga_lcpm, dsmz_lcpm, config, output_dir, cohort)
  } else {
    logmsg("  %s", subtype_info$message)
  }
  
  # === STAGE ANALYSIS ===
  if (config$run_stage) {
    analyze_stage_signatures(tcga_data$se, tcga_data$counts, cohort_indices, 
                            dsmz_lcpm, config$stage_lfc_min, output_dir)
  }
  
  logmsg("=== Completed: %s ===", cohort)
}

#' Process subtype analysis for a cohort
process_subtype_analysis <- function(subtype_info, tcga_lcpm, dsmz_lcpm, config, 
                                   output_dir, cohort) {
  logmsg("  Processing subtypes: %s", paste(subtype_info$subtypes, collapse = ", "))
  
  # Perform train/test split (80/20)
  set.seed(42)
  train_indices <- unlist(tapply(names(subtype_info$labels), subtype_info$labels, 
                                function(samples) {
    sample(samples, floor(0.8 * length(samples)))
  }))
  
  test_indices <- setdiff(names(subtype_info$labels), train_indices)
  
  # Split data
  train_labels <- subtype_info$labels[train_indices]
  test_labels <- subtype_info$labels[test_indices]
  
  train_counts <- subtype_info$counts[, train_indices, drop = FALSE]
  test_counts <- subtype_info$counts[, test_indices, drop = FALSE]
  
  # Normalize training and test data
  train_lcpm <- voom_quantile_norm(train_counts)
  test_lcpm <- voom_quantile_norm(test_counts)
  
  # Build subtype templates
  templates <- build_subtype_templates(train_lcpm, train_labels, 
                                     config$lfc_min, config$fdr_max)
  
  # Filter templates for cell-line relevance
  filtered_templates <- filter_templates_for_cell_lines(
    templates, train_lcpm, dsmz_lcpm, 
    config$abs_lfc_tumor_cl_max, config$dsmz_top50_min_samples
  )
  
  # Save and plot template sizes
  template_sizes <- plot_template_sizes(filtered_templates, cohort, 
                                       file.path(output_dir, "template_sizes_barplot.pdf"))
  
  # Evaluate on test set
  combined_expr <- cbind(train_lcpm, test_lcpm)
  combined_z <- z_score_by_gene(combined_expr)
  test_z <- combined_z[, colnames(test_lcpm), drop = FALSE]
  
  test_predictions <- ntp_predict(test_z, filtered_templates)
  accuracy <- plot_confusion_matrix(test_labels, test_predictions$pred, cohort,
                                  file.path(output_dir, "confusion_matrix_test.pdf"))
  
  logmsg("  %s hold-out accuracy: %.1f%%", cohort, 100 * accuracy)
  
  # If accuracy is sufficient, predict DSMZ subtypes
  if (accuracy >= config$accuracy_pass) {
    logmsg("  Accuracy >= %.0f%%; classifying DSMZ.", 100 * config$accuracy_pass)
    
    # Z-score combined TCGA + DSMZ data
    all_expr_z <- z_score_by_gene(cbind(tcga_lcpm, dsmz_lcpm))
    dsmz_z <- all_expr_z[, colnames(dsmz_lcpm), drop = FALSE]
    
    # Predict DSMZ subtypes
    dsmz_predictions <- ntp_predict(dsmz_z, filtered_templates)
    
    # Save and plot results
    plot_dsmz_scores(dsmz_predictions, cohort, 
                    file.path(output_dir, "dsmz_ntp_scores_by_subtype.pdf"))
  } else {
    logmsg("  Accuracy below threshold; skipping DSMZ classification.")
  }
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

#' Main pipeline execution
main <- function() {
  # Create output directory
  dir.create(CONFIG$outdir, showWarnings = FALSE, recursive = TRUE)
  
  logmsg("Starting DSMZ-TCGA Subtype Pipeline")
  logmsg("Configuration:")
  logmsg("  Output directory: %s", CONFIG$outdir)
  logmsg("  LFC threshold: %.2f", CONFIG$lfc_min)
  logmsg("  FDR threshold: %.3f", CONFIG$fdr_max)
  logmsg("  Accuracy threshold: %.0f%%", 100 * CONFIG$accuracy_pass)
  
  # === LOAD DATA ===
  tcga_data <- load_tcga_data(CONFIG$tcga_se_rds)
  dsmz_data <- load_dsmz_data(CONFIG$dsmz_rds, CONFIG$dsmz_meta_csv)
  
  # Find common genes between datasets
  common_genes <- intersect(rownames(tcga_data$counts), rownames(dsmz_data$counts))
  logmsg("Total shared genes between TCGA and DSMZ: %d", length(common_genes))
  
  # Filter both datasets to common genes
  tcga_data$counts <- tcga_data$counts[common_genes, , drop = FALSE]
  dsmz_data$counts <- dsmz_data$counts[common_genes, , drop = FALSE]
  
  # === LOAD SUBTYPE INFORMATION ===
  subtype_table <- tryCatch({
    TCGAbiolinks::PanCancerAtlas_subtypes()
  }, error = function(e) {
    logmsg("Warning: PanCancerAtlas_subtypes() not available. Subtype analysis will be limited.")
    NULL
  })
  
  # === DETERMINE COHORTS TO PROCESS ===
  tcga_labels <- extract_tcga_labels(tcga_data$se)
  available_cohorts <- sort(unique(tcga_labels))
  
  # Check for user-specified target cohorts
  target_cohorts_env <- Sys.getenv("TARGET_COHORTS", "")
  if (target_cohorts_env != "") {
    target_cohorts <- strsplit(target_cohorts_env, ",")[[1]]
    target_cohorts <- trimws(target_cohorts)
    cohorts_to_process <- intersect(target_cohorts, available_cohorts)
    
    if (length(cohorts_to_process) == 0) {
      stop("No valid cohorts found in TARGET_COHORTS environment variable")
    }
    
    logmsg("Processing user-specified cohorts: %s", paste(cohorts_to_process, collapse = ", "))
    
    missing_cohorts <- setdiff(target_cohorts, available_cohorts)
    if (length(missing_cohorts) > 0) {
      logmsg("Warning: Requested cohorts not found: %s", paste(missing_cohorts, collapse = ", "))
    }
  } else {
    cohorts_to_process <- available_cohorts
    logmsg("Processing all available cohorts: %s", paste(cohorts_to_process, collapse = ", "))
  }
  
  # === PROCESS EACH COHORT ===
  results_summary <- list()
  
  for (cohort in cohorts_to_process) {
    tryCatch({
      process_cohort(cohort, tcga_data, dsmz_data, subtype_table, CONFIG)
      results_summary[[cohort]] <- "SUCCESS"
    }, error = function(e) {
      logmsg("ERROR processing cohort %s: %s", cohort, e$message)
      results_summary[[cohort]] <- paste("ERROR:", e$message)
    })
  }
  
  # === GENERATE SUMMARY REPORT ===
  generate_summary_report(results_summary, CONFIG$outdir)
  
  logmsg("Pipeline completed successfully!")
  logmsg("Results saved to: %s", CONFIG$outdir)
}

#' Generate a summary report of all processed cohorts
generate_summary_report <- function(results_summary, output_dir) {
  summary_file <- file.path(output_dir, "pipeline_summary.txt")
  
  cat("DSMZ-TCGA Subtype Pipeline Summary\n", file = summary_file)
  cat("==================================\n\n", file = summary_file, append = TRUE)
  cat(sprintf("Run date: %s\n", Sys.time()), file = summary_file, append = TRUE)
  cat(sprintf("Total cohorts processed: %d\n\n", length(results_summary)), 
      file = summary_file, append = TRUE)
  
  # Count successes and failures
  successes <- sum(results_summary == "SUCCESS")
  failures <- length(results_summary) - successes
  
  cat(sprintf("Successful: %d\n", successes), file = summary_file, append = TRUE)
  cat(sprintf("Failed: %d\n\n", failures), file = summary_file, append = TRUE)
  
  # List results for each cohort
  cat("Cohort Results:\n", file = summary_file, append = TRUE)
  cat("---------------\n", file = summary_file, append = TRUE)
  
  for (cohort in names(results_summary)) {
    cat(sprintf("%-15s: %s\n", cohort, results_summary[[cohort]]), 
        file = summary_file, append = TRUE)
  }
  
  # Add configuration information
  cat("\nConfiguration Used:\n", file = summary_file, append = TRUE)
  cat("-------------------\n", file = summary_file, append = TRUE)
  cat(sprintf("LFC minimum: %.2f\n", CONFIG$lfc_min), file = summary_file, append = TRUE)
  cat(sprintf("FDR maximum: %.3f\n", CONFIG$fdr_max), file = summary_file, append = TRUE)
  cat(sprintf("Accuracy threshold: %.0f%%\n", 100 * CONFIG$accuracy_pass), 
      file = summary_file, append = TRUE)
  cat(sprintf("Min samples per subtype: %d\n", CONFIG$min_per_subtype), 
      file = summary_file, append = TRUE)
  cat(sprintf("Stage analysis: %s\n", ifelse(CONFIG$run_stage, "Enabled", "Disabled")), 
      file = summary_file, append = TRUE)
  
  logmsg("Summary report saved to: %s", summary_file)
}

#' Utility function to check if all required packages are installed
check_dependencies <- function() {
  required_packages <- c(
    "tidyverse", "SummarizedExperiment", "edgeR", "limma", 
    "pheatmap", "ggrepel", "TCGAbiolinks"
  )
  
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  
  if (length(missing_packages) > 0) {
    stop(sprintf("Missing required packages: %s\nPlease install with: install.packages(c(%s))", 
                paste(missing_packages, collapse = ", "),
                paste(sprintf("'%s'", missing_packages), collapse = ", ")))
  }
  
  logmsg("All required packages are available")
}

#' Utility function to validate input files exist
validate_input_files <- function(config) {
  required_files <- c(
    config$tcga_se_rds,
    config$dsmz_rds,
    config$dsmz_meta_csv
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    stop(sprintf("Missing required input files:\n%s", paste(missing_files, collapse = "\n")))
  }
  
  logmsg("All required input files found")
}

# ==============================================================================
# SCRIPT EXECUTION
# ==============================================================================

# Run the pipeline if this script is executed directly
if (!interactive()) {
  tryCatch({
    # Check dependencies and input files
    check_dependencies()
    validate_input_files(CONFIG)
    
    # Run main pipeline
    main()
    
  }, error = function(e) {
    logmsg("FATAL ERROR: %s", e$message)
    quit(status = 1)
  })
}