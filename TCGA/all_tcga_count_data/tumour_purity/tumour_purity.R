## Load required libraries ---------------------------------------------------
#
suppressPackageStartupMessages({
  library(tidyestimate)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tibble)
})

map_ensembl_to_hgnc <- function(expr_matrix, gtf_path) {
  # Skip if rownames are already gene symbols (not Ensembl IDs)
  if (!any(grepl("^ENSG\\d", rownames(expr_matrix), perl = TRUE))) {
    message("[INFO] Row names appear to be gene symbols. Skipping ID mapping.")
    return(expr_matrix)
  }
  
  message("[INFO] Mapping Ensembl IDs to HGNC symbols using GTF: ", gtf_path)
  
  # Remove version suffix from Ensembl IDs (e.g., ENSG000001.1 → ENSG000001)
  ensembl_ids <- sub("\\.\\d+$", "", rownames(expr_matrix))
  
  # Read GTF and extract gene_id ↔ gene_name mapping
  gtf <- data.table::fread(
    file = gtf_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    data.table = TRUE,
    col.names = c("seqname", "source", "feature", "start", "end",
                  "score", "strand", "frame", "attributes"),
    showProgress = FALSE
  )
  
  # Extract gene_id and gene_name from attributes
  gene_map <- gtf[feature == "gene", .(
    gene_id = sub('.*gene_id "([^"]+)".*', '\\1', attributes),
    gene_name = sub('.*gene_name "([^"]+)".*', '\\1', attributes)
  )]
  
  # Remove duplicates and empty names
  gene_map <- gene_map[!duplicated(gene_id) & nzchar(gene_name)]
  
  # Create named vector: Ensembl ID → HGNC symbol
  id_to_symbol <- setNames(gene_map$gene_name, gene_map$gene_id)
  
  # Map and filter
  mapped_symbols <- id_to_symbol[ensembl_ids]
  valid <- !is.na(mapped_symbols) & nzchar(mapped_symbols)
  
  message(sprintf(
    "[INFO] Successfully mapped %d / %d Ensembl IDs (%.1f%%).",
    sum(valid), length(ensembl_ids), 100 * mean(valid)
  ))
  
  # Subset matrix and assign uppercase HGNC symbols
  expr_mapped <- expr_matrix[valid, , drop = FALSE]
  rownames(expr_mapped) <- toupper(mapped_symbols[valid])
  
  # Collapse duplicate gene symbols by summing
  before <- nrow(expr_mapped)
  expr_mapped <- as.data.frame(expr_mapped) %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(dplyr::across(everything(), sum), .groups = "drop") %>%
    tibble::column_to_rownames("symbol")
  expr_mapped <- as.matrix(expr_mapped)
  
  if (nrow(expr_mapped) < before) {
    message(sprintf("[INFO] Collapsed duplicate HGNC symbols: %d → %d genes.", before, nrow(expr_mapped)))
  }
  
  return(expr_mapped)
}


## Helper: Select appropriate tidyestimate scoring function ------------------
get_estimate_scoring_function <- function() {
  if ("estimate_score" %in% getNamespaceExports("tidyestimate")) {
    return(tidyestimate::estimate_score)
  }
  if ("estimate_scores" %in% getNamespaceExports("tidyestimate")) {
    return(tidyestimate::estimate_scores)
  }
  stop("[ERROR] Neither 'estimate_score' nor 'estimate_scores' found in tidyestimate package.")
}


## Helper: Standardize and enrich ESTIMATE output ----------------------------
standardize_estimate_scores <- function(scores_df) {
  # Convert to data.frame with sample rownames
  if ("sample" %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df$sample <- NULL
  }
  if (is.matrix(scores_df)) {
    scores_df <- as.data.frame(scores_df)
  }
  
  # Assign default column names if missing
  if (ncol(scores_df) == 3 && (is.null(colnames(scores_df)) || all(colnames(scores_df) == ""))) {
    colnames(scores_df) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
    message("[INFO] Assigned default ESTIMATE column names.")
  }
  
  # Normalize column names using pattern matching
  col_lower <- tolower(colnames(scores_df))
  new_names <- colnames(scores_df)
  new_names[grepl("strom", col_lower)] <- "StromalScore"
  new_names[grepl("imm",   col_lower)] <- "ImmuneScore"
  new_names[grepl("est",   col_lower)] <- "ESTIMATEScore"
  colnames(scores_df) <- new_names
  
  # Enforce standard names for 3-column output
  standard_cols <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
  if (ncol(scores_df) == 3 && length(intersect(colnames(scores_df), standard_cols)) < 3) {
    colnames(scores_df) <- standard_cols
    message("[INFO] Enforced standard column names by position.")
  }
  
  # Derive tumourPurity using Yoshihara et al. formula if missing
  if (!"tumourPurity" %in% colnames(scores_df) && "ESTIMATEScore" %in% colnames(scores_df)) {
    scores_df$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)
    scores_df$tumourPurity <- pmin(1, pmax(0, scores_df$tumourPurity))  # Clamp to [0,1]
    message("[INFO] Derived tumourPurity from ESTIMATEScore.")
  }
  
  # Final validation
  required_cols <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "tumourPurity")
  missing_cols <- setdiff(required_cols, colnames(scores_df))
  if (length(missing_cols) > 0) {
    warning("[WARN] Missing expected columns: ", paste(missing_cols, collapse = ", "))
  } else {
    message("[INFO] All required ESTIMATE columns are present.")
  }
  
  return(scores_df)
}


## Helper: Filter samples by tumour purity threshold -------------------------
filter_samples_by_purity <- function(se_object, purity_threshold = 0.7) {
  message("[INFO] Applying tumour purity filtering (threshold >= ", purity_threshold, ")...")
  
  if (!"purity" %in% colnames(colData(se_object))) {
    warning("[WARN] 'purity' column not found in colData. Skipping filtering.")
    return(se_object)
  }
  
  purity_vals <- colData(se_object)$purity
  message("[INFO] tumour purity summary:")
  print(summary(purity_vals))
  
  # Define keep mask
  keep <- !is.na(purity_vals) & purity_vals >= purity_threshold
  
  # Apply filter
  se_filtered <- se_object[, keep]
  
  # Report results
  n_total <- ncol(se_object)
  n_na <- sum(is.na(purity_vals))
  n_low <- sum(!is.na(purity_vals) & purity_vals < purity_threshold)
  n_kept <- sum(keep)
  
  message(sprintf(
    "[INFO] Filtering summary:\n  Total samples: %d\n  NA purity: %d\n  Low purity (< %.2f): %d\n  Kept (>= %.2f): %d (%.1f%%)",
    n_total, n_na, purity_threshold, n_low, purity_threshold, n_kept, 100 * n_kept / n_total
  ))
  
  return(se_filtered)
}


## Main Function: Run full tidyestimate + purity analysis --------------------
run_tidyestimate_purity_analysis <- function(
  tcga_se_path = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds",
  tcga_se_updated_path = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_update.rds",
  tcga_se_filtered_path = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds",
  output_dir = "/home/chu25/dsmz/results/tumor_purity_results",
  gtf_path = "/home/chu25/data/GTF/Homo_sapiens.GRCh38.107.gtf",
  assay_name = NULL,
  save_updated = TRUE,
  purity_threshold = 0.7,
  apply_purity_filter = TRUE
) {
  message("[INFO] Starting tidyestimate purity analysis pipeline...")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load SummarizedExperiment
  if (!file.exists(tcga_se_path)) {
    stop(sprintf("[ERROR] Input file not found: %s", tcga_se_path))
  }
  tcga_se <- readRDS(tcga_se_path)
  
  # Extract expression matrix
  expr_matrix <- if (is.null(assay_name)) SummarizedExperiment::assay(tcga_se) else SummarizedExperiment::assay(tcga_se, assay_name)
  message(sprintf("[INFO] Input expression matrix: %d genes × %d samples", nrow(expr_matrix), ncol(expr_matrix)))
  
  # Map to HGNC symbols
  expr_hgnc <- map_ensembl_to_hgnc(expr_matrix, gtf_path)
  
  # Filter to ESTIMATE signature genes
  message("[INFO] Filtering to ESTIMATE signature genes...")
  expr_sig <- suppressMessages(
    tidyestimate::filter_common_genes(
      as.data.frame(expr_hgnc),
      id = "hgnc_symbol",
      tidy = FALSE,
      tell_missing = TRUE,
      find_alias = TRUE
    )
  )
  
  # Run ESTIMATE scoring
  score_fn <- get_estimate_scoring_function()
  message("[INFO] Computing ESTIMATE scores...")
  raw_scores <- score_fn(expr_sig, is_affymetrix = FALSE)
  
  # Standardize and enrich scores
  estimate_scores <- standardize_estimate_scores(raw_scores)
  message(sprintf("[INFO] Final scores: %d samples × %d metrics", nrow(estimate_scores), ncol(estimate_scores)))
  
  # Save scores to CSV
  scores_csv <- file.path(output_dir, "tcga_tidyestimate_scores.csv")
  write.csv(estimate_scores, scores_csv, row.names = TRUE)
  message("[INFO] ESTIMATE scores saved to: ", scores_csv)
  
  # Update SummarizedExperiment with tumourPurity
  result <- list(
    scores = estimate_scores,
    output_csv = scores_csv,
    tcga_se_updated = NULL,
    tcga_se_filtered = NULL,
    filtered_counts = NULL,
    filtered_output_path = NULL
  )
  
  if (save_updated && "tumourPurity" %in% colnames(estimate_scores)) {
    common_samples <- intersect(rownames(estimate_scores), colnames(tcga_se))
    colData(tcga_se)$purity <- NA_real_
    colData(tcga_se)[common_samples, "purity"] <- estimate_scores[common_samples, "tumourPurity"]
    
    saveRDS(tcga_se, tcga_se_updated_path)
    message(sprintf("[INFO] Updated colData with tumourPurity (%d/%d samples). Saved to: %s",
                    length(common_samples), ncol(tcga_se), tcga_se_updated_path))
    
    result$tcga_se_updated <- tcga_se
  }
  
  # Apply purity filtering
  if (apply_purity_filter && "tumourPurity" %in% colnames(estimate_scores)) {
    tcga_se_filtered <- filter_samples_by_purity(tcga_se, purity_threshold)
    saveRDS(tcga_se_filtered, tcga_se_filtered_path)
    message("[INFO] Filtered SummarizedExperiment saved to: ", tcga_se_filtered_path)
    
    filtered_counts <- SummarizedExperiment::assay(tcga_se_filtered)
    message(sprintf("[INFO] Filtered expression matrix: %d genes × %d samples",
                    nrow(filtered_counts), ncol(filtered_counts)))
    
    result$tcga_se_filtered <- tcga_se_filtered
    result$filtered_counts <- filtered_counts
    result$filtered_output_path <- tcga_se_filtered_path
  }
  
  # Generate diagnostic plots
  if ("tumourPurity" %in% colnames(estimate_scores)) {
    message("[INFO] Generating diagnostic correlation plots...")
    
    plot_data <- estimate_scores %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample")
    
    p1 <- ggplot(plot_data, aes(x = StromalScore, y = tumourPurity)) +
      geom_point(alpha = 0.6, size = 1.2) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggpubr::stat_cor(method = "spearman", label.y = 0.1) +
      ggpubr::stat_cor(method = "pearson",  label.y = 0.05) +
      theme_bw(base_size = 12) +
      labs(title = "tumour Purity vs Stromal Score", x = "Stromal Score", y = "tumour Purity")
    
    p2 <- ggplot(plot_data, aes(x = ImmuneScore, y = tumourPurity)) +
      geom_point(alpha = 0.6, size = 1.2) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggpubr::stat_cor(method = "spearman", label.y = 0.1) +
      ggpubr::stat_cor(method = "pearson",  label.y = 0.05) +
      theme_bw(base_size = 12) +
      labs(title = "tumour Purity vs Immune Score", x = "Immune Score", y = "tumour Purity")
    
    p3 <- ggplot(plot_data, aes(x = tumourPurity)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8, color = "black") +
      geom_vline(xintercept = purity_threshold, color = "red", linetype = "dashed", size = 1) +
      annotate("text", x = purity_threshold + 0.05, y = Inf, vjust = 1.5, hjust = 0,
               label = sprintf("Threshold = %.2f", purity_threshold), color = "red") +
      theme_bw(base_size = 12) +
      labs(title = "Distribution of tumour Purity", x = "tumour Purity", y = "Sample Count")
    
    # Save combined plot
    plot_file <- file.path(output_dir, "purity_correlation_diagnostics.pdf")
    combined_plot <- ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = "AUTO")
    ggsave(plot_file, combined_plot, width = 14, height = 10, dpi = 300)
    message("[INFO] Diagnostic plots saved to: ", plot_file)
  }
  
  invisible(result)
}

## Execute Pipeline ----------------------------------------------------------
tryCatch({
  analysis_result <- run_tidyestimate_purity_analysis(
    apply_purity_filter = TRUE,
    purity_threshold = 0.7
  )
  
  message("[SUCCESS] tidyestimate purity analysis completed successfully.")
  
  if (!is.null(analysis_result$tcga_se_filtered)) {
    message(sprintf(
      "[INFO] Filtered dataset ready: %d genes × %d samples",
      nrow(analysis_result$filtered_counts), ncol(analysis_result$filtered_counts)
    ))
  }
}, error = function(e) {
  message("[ERROR] Pipeline failed: ", e$message)
})
