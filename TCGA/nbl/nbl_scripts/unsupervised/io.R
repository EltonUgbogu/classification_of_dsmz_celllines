---
title: "NBL Tumour Purity Estimation with tidyestimate"
author: "chu25"
date: "`r format(Sys.time(), '%Y-%m-%d')`"
output:
  pdf_document:
    latex_engine: xelatex
    toc: true
    toc_depth: 2
    fig_width: 7
    fig_height: 5
    fig_caption: true
---

```{r setup, include=FALSE}
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  error = FALSE
)

# Required packages
required_packages <- c(
  "tidyestimate", "SummarizedExperiment", "dplyr", "ggplot2", "ggpubr",
  "data.table", "tibble", "arrow", "rtracklayer"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[INFO] Installing package: %s", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Load libraries
suppressPackageStartupMessages({
  library(tidyestimate)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tibble)
  library(arrow)
  library(rtracklayer)
})
```
```{r}
map_ensembl_to_hgnc <- function(expr_matrix, gtf_path) {
  # Skip if already HGNC symbols
  if (!any(grepl("^ENSG\\d", rownames(expr_matrix), perl = TRUE))) {
    message("[INFO] Row names are already gene symbols. Skipping mapping.")
    return(expr_matrix)
  }

  if (!file.exists(gtf_path)) {
    stop(sprintf("[ERROR] GTF file not found: %s", gtf_path))
  }

  message("[INFO] Loading GTF for Ensembl → HGNC mapping: ", gtf_path)
  gtf <- rtracklayer::import(gtf_path)
  gtf_genes <- gtf[gtf$type == "gene"]

  gene_map <- data.frame(
    gene_id = gtf_genes$gene_id,
    gene_name = gtf_genes$gene_name,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(gene_name), gene_name != "", !duplicated(gene_id))

  # Remove version from Ensembl IDs
  ensembl_ids <- sub("\\.\\d+$", "", rownames(expr_matrix))
  id_to_symbol <- setNames(gene_map$gene_name, gene_map$gene_id)

  mapped_symbols <- id_to_symbol[ensembl_ids]
  valid <- !is.na(mapped_symbols) & nzchar(mapped_symbols)

  message(sprintf(
    "[INFO] Mapped %d / %d Ensembl IDs (%.1f%%)",
    sum(valid), length(ensembl_ids), 100 * mean(valid)
  ))

  expr_mapped <- expr_matrix[valid, , drop = FALSE]
  rownames(expr_mapped) <- toupper(mapped_symbols[valid])

  # Collapse duplicates by summing
  before <- nrow(expr_mapped)
  expr_mapped <- expr_mapped %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(across(everything(), sum), .groups = "drop") %>%
    tibble::column_to_rownames("symbol") %>%
    as.matrix()

  if (nrow(expr_mapped) < before) {
    message(sprintf("[INFO] Collapsed duplicates: %d → %d genes", before, nrow(expr_mapped)))
  }

  return(expr_mapped)
}
```
```{r}
get_estimate_scoring_function <- function() {
  if (exists("estimate_score", where = asNamespace("tidyestimate"))) {
    return(tidyestimate::estimate_score)
  }
  if (exists("estimate_scores", where = asNamespace("tidyestimate"))) {
    return(tidyestimate::estimate_scores)
  }
  stop("[ERROR] No valid scoring function found in tidyestimate.")
}
```
```{r}
standardize_estimate_scores <- function(scores_df) {
  if (inherits(scores_df, "data.frame") && "sample" %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df$sample <- NULL
  }
  scores_df <- as.data.frame(scores_df)

  # Normalize column names
  col_lower <- tolower(colnames(scores_df))
  new_names <- colnames(scores_df)
  new_names[grepl("stromal", col_lower)] <- "StromalScore"
  new_names[grepl("immune", col_lower)] <- "ImmuneScore"
  new_names[grepl("estimate", col_lower)] <- "ESTIMATEScore"
  colnames(scores_df) <- new_names

  # Derive tumourPurity if missing
  if (!"tumourPurity" %in% colnames(scores_df) && "ESTIMATEScore" %in% colnames(scores_df)) {
    scores_df$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)
    scores_df$tumourPurity <- pmin(1, pmax(0, scores_df$tumourPurity))
    message("[INFO] Derived tumourPurity from ESTIMATEScore.")
  }

  required <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "tumourPurity")
  missing <- setdiff(required, colnames(scores_df))
  if (length(missing) > 0) {
    warning("[WARN] Missing columns: ", paste(missing, collapse = ", "))
  } else {
    message("[INFO] All required ESTIMATE columns present.")
  }

  return(scores_df)
}
```
```{r}
filter_samples_by_purity <- function(se_object, purity_threshold = 0.7) {
  if (!"purity" %in% colnames(colData(se_object))) {
    warning("[WARN] 'purity' column missing in colData. Skipping filter.")
    return(se_object)
  }

  purity_vals <- colData(se_object)$purity
  keep <- !is.na(purity_vals) & purity_vals >= purity_threshold
  se_filtered <- se_object[, keep]

  message(sprintf(
    "[INFO] Purity filter (≥ %.2f): %d → %d samples (%.1f%% kept)",
    purity_threshold, ncol(se_object), ncol(se_filtered), 100 * mean(keep)
  ))
  return(se_filtered)
}
```
```{r}
run_tidyestimate_purity_analysis <- function(
  se_input_path,
  gtf_path,
  se_updated_path = NULL,
  se_filtered_path = NULL,
  output_dir = getwd(),
  assay_name = NULL,
  save_updated = TRUE,
  purity_threshold = 0.7,
  apply_purity_filter = TRUE
) {
  message("[INFO] Starting tidyestimate tumour purity analysis...")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  tumour_name <- "nbl"

  # --- Load input: Parquet or RDS ---
  if (grepl("\\.parquet$", se_input_path, ignore.case = TRUE)) {
    message("[INFO] Reading Parquet file: ", se_input_path)
    count_matrix <- arrow::read_parquet(se_input_path)
    if (!is.matrix(count_matrix)) {
      stop("Parquet must contain a matrix with samples as columns, genes as rows.")
    }
    # Create minimal SummarizedExperiment
    se_object <- SummarizedExperiment(assays = list(counts = count_matrix))
  } else if (grepl("\\.rds$", se_input_path, ignore.case = TRUE)) {
    message("[INFO] Reading RDS SummarizedExperiment: ", se_input_path)
    se_object <- readRDS(se_input_path)
  } else {
    stop("Input must be .parquet or .rds")
  }

  # Extract assay
  expr_matrix <- if (is.null(assay_name)) assay(se_object) else assay(se_object, assay_name)
  message(sprintf("[INFO] Input matrix: %d genes × %d samples", nrow(expr_matrix), ncol(expr_matrix)))

  # Map to HGNC
  expr_hgnc <- map_ensembl_to_hgnc(expr_matrix, gtf_path)

  # Filter to ESTIMATE genes
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

  # Run scoring
  score_fn <- get_estimate_scoring_function()
  message("[INFO] Running ESTIMATE scoring...")
  raw_scores <- score_fn(expr_sig, is_affymetrix = FALSE)
  estimate_scores <- standardize_estimate_scores(raw_scores)

  # Save scores
  scores_csv <- file.path(output_dir, paste0(tumour_name, "_tidyestimate_scores.csv"))
  write.csv(estimate_scores, scores_csv, row.names = TRUE)
  message("[INFO] Scores saved: ", scores_csv)

  result <- list(scores = estimate_scores, output_csv = scores_csv)

  # Update SE with purity
  if (save_updated && !is.null(se_updated_path) && "tumourPurity" %in% colnames(estimate_scores)) {
    common_samples <- intersect(rownames(estimate_scores), colnames(se_object))
    colData(se_object)$purity <- NA_real_
    colData(se_object)[common_samples, "purity"] <- estimate_scores[common_samples, "tumourPurity"]
    saveRDS(se_object, se_updated_path)
    message(sprintf("[INFO] Updated SE saved: %s (%d samples)", se_updated_path, length(common_samples)))
    result$se_updated <- se_object
  }

  # Apply filtering
  if (apply_purity_filter && "purity" %in% colnames(colData(se_object)) && !is.null(se_filtered_path)) {
    se_filtered <- filter_samples_by_purity(se_object, purity_threshold)
    saveRDS(se_filtered, se_filtered_path)
    message("[INFO] Filtered SE saved: ", se_filtered_path)
    result$se_filtered <- se_filtered
  }

  # Diagnostic plots
  if ("tumourPurity" %in% colnames(estimate_scores)) {
    plot_data <- estimate_scores %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample")

    p1 <- ggplot(plot_data, aes(x = StromalScore, y = tumourPurity)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
      stat_cor(method = "spearman", label.y = 0.15) +
      theme_bw() +
      labs(title = "Purity vs Stromal Score")

    p2 <- ggplot(plot_data, aes(x = ImmuneScore, y = tumourPurity)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
      stat_cor(method = "spearman", label.y = 0.15) +
      theme_bw() +
      labs(title = "Purity vs Immune Score")

    p3 <- ggplot(plot_data, aes(x = tumourPurity)) +
      geom_histogram(bins = 40, fill = "skyblue", color = "black", alpha = 0.8) +
      geom_vline(xintercept = purity_threshold, color = "red", linetype = "dashed") +
      annotate("text", x = purity_threshold + 0.05, y = Inf, vjust = 1.5,
               label = sprintf("Threshold = %.2f", purity_threshold), color = "red") +
      theme_bw() +
      labs(title = "Tumour Purity Distribution")

    plot_file <- file.path(output_dir, paste0(tumour_name, "_purity_diagnostics.pdf"))
    ggsave(plot_file, ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = "AUTO"),
           width = 14, height = 10, dpi = 300)
    message("[INFO] Plots saved: ", plot_file)
  }

  invisible(result)
}
```

```{r}
# === CONFIGURATION ===
NBL_INPUT_PARQUET <- "/home/chu25/data/nbl/final_nbl_count_matrix_hgnc.parquet"
GTF_FILE_PATH     <- "/home/chu25/data/GTF/Homo_sapiens.GRCh38.107.gtf"
OUTPUT_DIR        <- "/home/chu25/nbl/results/tumour_purity_results"
UPDATED_SE_PATH   <- "/home/chu25/data/nbl/ALL_nbl_STAR_Counts_SE_updated.rds"
FILTERED_SE_PATH  <- "/home/chu25/data/nbl/ALL_nbl_STAR_Counts_SE_filtered.rds"

# === RUN ANALYSIS ===
tryCatch({
  result <- run_tidyestimate_purity_analysis(
    se_input_path = NBL_INPUT_PARQUET,
    gtf_path = GTF_FILE_PATH,
    se_updated_path = UPDATED_SE_PATH,
    se_filtered_path = FILTERED_SE_PATH,
    output_dir = OUTPUT_DIR,
    purity_threshold = 0.7,
    apply_purity_filter = TRUE,
    save_updated = TRUE
  )
  message("[SUCCESS] Pipeline completed. See: ", OUTPUT_DIR)
}, error = function(e) {
  message("[ERROR] Pipeline failed: ", e$message)
})
```

































