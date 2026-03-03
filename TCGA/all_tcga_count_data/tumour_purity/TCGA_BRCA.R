#!/usr/bin/env Rscript
#TCGA_BRCA.R
## Load required libraries ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyestimate)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tibble)
  library(RColorBrewer)
  library(tidyr)
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


## Helper: Subset SummarizedExperiment to one TCGA project ------------------
subset_se_to_tcga_project <- function(se, project = "TCGA-BRCA") {
  cd <- as.data.frame(SummarizedExperiment::colData(se))

  # Common column names seen in TCGA-derived objects
  candidate_cols <- c(
    "project_id", "project", "tcga_project", "cohort", "disease", "study",
    "cancer_type", "primary_site", "tumor_type", "tumour_type"
  )

  hit_cols <- intersect(candidate_cols, colnames(cd))

  # If none of the common columns exist, try any column containing "project"
  if (length(hit_cols) == 0) {
    hit_cols <- grep("project", colnames(cd), value = TRUE, ignore.case = TRUE)
  }

  if (length(hit_cols) == 0) {
    stop("[ERROR] Could not find a project/cohort column in colData() to filter TCGA-BRCA.")
  }

  # Use the first matching column that actually contains the requested project string
  chosen <- NA_character_
  for (cc in hit_cols) {
    vals <- as.character(cd[[cc]])
    if (any(grepl(project, vals, fixed = TRUE))) {
      chosen <- cc
      break
    }
  }
  if (is.na(chosen)) {
    # fallback: try case-insensitive match
    for (cc in hit_cols) {
      vals <- as.character(cd[[cc]])
      if (any(grepl(project, vals, ignore.case = TRUE))) {
        chosen <- cc
        break
      }
    }
  }

  if (is.na(chosen)) {
    # Helpful debug output
    msg <- paste0(
      "[ERROR] None of the candidate columns contained '", project, "'.\n",
      "Checked columns: ", paste(hit_cols, collapse = ", "), "\n",
      "Example values (first 5) from first hit column '", hit_cols[1], "': ",
      paste(head(as.character(cd[[hit_cols[1]]]), 5), collapse = " | ")
    )
    stop(msg)
  }

  keep <- grepl(project, as.character(cd[[chosen]]), ignore.case = TRUE)
  message(sprintf("[INFO] Subsetting to %s using colData$%s: kept %d / %d samples.",
                  project, chosen, sum(keep), ncol(se)))

  se[, keep]
}


## Helper: Collapse duplicate sample IDs by summing counts --------------------
collapse_duplicate_samples <- function(mat) {
  stopifnot(!is.null(colnames(mat)))
  if (anyDuplicated(colnames(mat)) == 0) return(mat)

  message("[INFO] Collapsing duplicated sample IDs by summing columns...")
  df <- as.data.frame(mat)
  df$gene_id <- rownames(mat)

  df_long <- tidyr::pivot_longer(df, -gene_id, names_to = "sample", values_to = "count")
  df_sum  <- df_long %>%
    dplyr::group_by(gene_id, sample) %>%
    dplyr::summarise(count = sum(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample, values_from = count)

  out <- as.matrix(df_sum[, -1])
  rownames(out) <- df_sum$gene_id
  storage.mode(out) <- "numeric"
  out
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
  output_dir = "/home/chu25/TCGA/all_tcga_count_data/tumour_purity/results/TCGA_BRCA",
  gtf_path = "/home/chu25/data/GTF/Homo_sapiens.GRCh38.107.gtf",
  assay_name = NULL,
  save_updated = TRUE,
  purity_threshold = 0.7,
  apply_purity_filter = TRUE,
  focus_project = "TCGA-BRCA"
) {
  message("[INFO] Starting tidyestimate purity analysis pipeline...")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set up log file
  log_file <- file.path(output_dir, paste0("tcga_purity_pipeline_", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                            ".log"))
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "output", split = TRUE)
  sink(log_conn, type = "message")
  
  message(sprintf("[INFO] Log file created: %s", log_file))
  message(sprintf("[INFO] Pipeline start time: %s", Sys.time()))
  
  # Wrap pipeline in tryCatch to ensure log cleanup
  tryCatch({
  
    # Load SummarizedExperiment
    if (!file.exists(tcga_se_path)) {
      stop(sprintf("[ERROR] Input file not found: %s", tcga_se_path))
    }
    tcga_se <- readRDS(tcga_se_path)
    
    # Subset to one TCGA project (BRCA) if requested
    if (!is.null(focus_project) && nzchar(focus_project)) {
      tcga_se <- subset_se_to_tcga_project(tcga_se, project = focus_project)
      
      # Make outputs project-specific so you don't overwrite pan-TCGA results
      output_dir <- file.path(output_dir, focus_project)
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      tcga_se_updated_path  <- sub("\\.rds$", paste0("_", focus_project, ".rds"), tcga_se_updated_path)
      tcga_se_filtered_path <- sub("\\.rds$", paste0("_", focus_project, ".rds"), tcga_se_filtered_path)
    }
    
    # Define tag early for use throughout
    tag <- if (!is.null(focus_project) && nzchar(focus_project)) focus_project else "TCGA_ALL"
    
    # Extract expression matrix
    expr_matrix <- if (is.null(assay_name)) SummarizedExperiment::assay(tcga_se) else SummarizedExperiment::assay(tcga_se, assay_name)
    message(sprintf("[INFO] Input expression matrix: %d genes × %d samples", nrow(expr_matrix), ncol(expr_matrix)))
    
    # Ensure unique sample IDs (critical for filter_common_genes)
    n_dup <- anyDuplicated(colnames(expr_matrix))
    if (n_dup > 0) {
      dup <- colnames(expr_matrix)[duplicated(colnames(expr_matrix))]
      message(sprintf("[WARN] Found %d duplicated sample IDs. Collapsing by summing counts.", length(unique(dup))))
      expr_matrix <- collapse_duplicate_samples(expr_matrix)
      message(sprintf("[INFO] After collapsing duplicates: %d genes × %d samples", nrow(expr_matrix), ncol(expr_matrix)))
    } else {
      message("[INFO] All sample IDs are unique.")
    }
    
    # Map to HGNC symbols
    expr_hgnc <- map_ensembl_to_hgnc(expr_matrix, gtf_path)
    
    # Filter to ESTIMATE signature genes
    message("[INFO] Filtering to ESTIMATE signature genes...")
    
    common_symbols <- unique(tidyestimate::common_genes$hgnc_symbol)
    keep_genes <- intersect(rownames(expr_hgnc), common_symbols)
    expr_sig_mat <- expr_hgnc[keep_genes, , drop = FALSE]
    message(sprintf("[INFO] Kept %d / %d ESTIMATE signature genes", length(keep_genes), length(common_symbols)))
    
    expr_sig <- as.data.frame(expr_sig_mat, check.names = FALSE) %>%
      tibble::rownames_to_column("hgnc_symbol")
    
    # Run ESTIMATE scoring
    score_fn <- get_estimate_scoring_function()
    message("[INFO] Computing ESTIMATE scores...")
    raw_scores <- score_fn(expr_sig, is_affymetrix = FALSE)
    
    # Standardize and enrich scores
    estimate_scores <- standardize_estimate_scores(raw_scores)
    message(sprintf("[INFO] Final scores: %d samples × %d metrics", nrow(estimate_scores), ncol(estimate_scores)))
    
    # Ensure score rownames match sample IDs (important for common_samples)
    if (nrow(estimate_scores) == ncol(tcga_se) &&
        !setequal(rownames(estimate_scores), colnames(tcga_se))) {
      message("[INFO] Forcing ESTIMATE score rownames to tcga_se sample IDs.")
      rownames(estimate_scores) <- colnames(tcga_se)
    }
    
    # Save scores to CSV
    scores_csv <- file.path(output_dir, sprintf("%s_tidyestimate_scores.csv", tag))
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
      # Record retention counts BEFORE filtering
      n_before <- ncol(tcga_se)
      
      tcga_se_filtered <- filter_samples_by_purity(tcga_se, purity_threshold)
      saveRDS(tcga_se_filtered, tcga_se_filtered_path)
      message("[INFO] Filtered SummarizedExperiment saved to: ", tcga_se_filtered_path)
      
      filtered_counts <- SummarizedExperiment::assay(tcga_se_filtered)
      message(sprintf("[INFO] Filtered expression matrix: %d genes × %d samples",
                      nrow(filtered_counts), ncol(filtered_counts)))
      
      # Record retention counts AFTER filtering
      n_after <- ncol(tcga_se_filtered)
      
      message(sprintf(
        "[RETENTION] %s: %d -> %d (%.1f%% kept)",
        tag, n_before, n_after, 100 * n_after / n_before
      ))
      
      # Simple retention plot (2 bars: Before vs After)
      df_ret <- data.frame(
        Stage = factor(c("Before", "After"), levels = c("Before", "After")),
        Count = c(n_before, n_after)
      )
      
      
      before_grey <- "#D3D3D3"
      after_col   <- RColorBrewer::brewer.pal(3, "Set2")[2]  # stable "Set2" colour
      
      p_ret <- ggplot(df_ret, aes(x = Stage, y = Count, fill = Stage)) +
        geom_col(width = 0.7, color = NA) +
        geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 5) +
        scale_fill_manual(values = c("Before" = before_grey, "After" = after_col), guide = "none") +
        labs(
          title = paste0(tag, ": Sample retention after purity filtering"),
          subtitle = sprintf("Purity ≥ %.1f | %d → %d (%.1f%% kept)",
                             purity_threshold, n_before, n_after, 100 * n_after / n_before),
          x = NULL,
          y = "Number of samples"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position  = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
        )
      
      ret_pdf <- file.path(output_dir, sprintf("%s_retention_before_after.pdf", tag))
      ggsave(ret_pdf, p_ret, width = 6.5, height = 5, dpi = 300)
      message("[SUCCESS] Retention plot saved: ", ret_pdf)
      
      result$tcga_se_filtered <- tcga_se_filtered
      result$filtered_counts <- filtered_counts
      result$filtered_output_path <- tcga_se_filtered_path
    }
    
    # Generate diagnostic plots
    if ("tumourPurity" %in% colnames(estimate_scores)) {
      message("[INFO] Generating diagnostic correlation plots...")
      
      
      pal <- RColorBrewer::brewer.pal(3, "Set2")[c(2, 1)]  # two colours
      
      plot_data <- estimate_scores %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::mutate(
          PurityGroup = cut(
            tumourPurity,
            breaks = c(-Inf, purity_threshold, Inf),
            labels = c(paste0("<", purity_threshold), paste0(">=", purity_threshold)),
            include.lowest = TRUE
          )
        )
      
      # A: Purity vs StromalScore 
      p1 <- ggplot(plot_data, aes(x = StromalScore, y = tumourPurity, colour = PurityGroup)) +
        geom_point(size = 2.2, alpha = 0.85, shape = 16) +
        geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
        scale_colour_manual(values = pal, name = "Purity") +
        ggpubr::stat_cor(method = "spearman",
                         label.y = 0.95, size = 4,
                         label.x.npc = "left", colour = "black") +
        ggpubr::stat_cor(method = "pearson",
                         label.y = 0.88, size = 4,
                         label.x.npc = "left", colour = "black") +
        labs(
          title = "TCGA-BRCA: Tumour Purity vs Stromal Score",
          x     = "Stromal Score",
          y     = "Tumour Purity"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5),
          legend.position  = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # B: Purity vs ImmuneScore 
      p2 <- ggplot(plot_data, aes(x = ImmuneScore, y = tumourPurity, colour = PurityGroup)) +
        geom_point(size = 2.2, alpha = 0.85, shape = 16) +
        geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
        scale_colour_manual(values = pal, name = "Purity") +
        ggpubr::stat_cor(method = "spearman",
                         label.y = 0.95, size = 4,
                         label.x.npc = "left", colour = "black") +
        ggpubr::stat_cor(method = "pearson",
                         label.y = 0.88, size = 4,
                         label.x.npc = "left", colour = "black") +
        labs(
          title = "TCGA-BRCA: Tumour Purity vs Immune Score",
          x     = "Immune Score",
          y     = "Tumour Purity"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5),
          legend.position  = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # C: Purity distribution 
      # Compute max count for stable annotation positioning
      y_top <- max(ggplot_build(
        ggplot(plot_data, aes(x = tumourPurity)) + geom_histogram(bins = 40)
      )$data[[1]]$count, na.rm = TRUE)
      
      p3 <- ggplot(plot_data, aes(x = tumourPurity)) +
        geom_histogram(bins = 40, alpha = 0.85, fill = pal[2], colour = NA) +
        geom_vline(xintercept = purity_threshold, linetype = "dashed",
                   colour = "firebrick", linewidth = 1) +
        annotate(
          "text",
          x = purity_threshold + 0.02,
          y = y_top * 0.95,
          vjust = 1,
          hjust = 0,
          label = paste("Threshold =", purity_threshold),
          colour = "firebrick",
          size = 4
        ) +
        labs(
          title = "TCGA-BRCA: Distribution of Tumour Purity",
          x     = "Tumour Purity",
          y     = "Number of Samples"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5),
          legend.position  = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
        )
      
      tag <- if (!is.null(focus_project) && nzchar(focus_project)) focus_project else "TCGA_ALL"
      combined_plot <- ggpubr::ggarrange(
        p1, p2, p3,
        ncol = 2, nrow = 2,
        labels = c("A", "B", "C"),
        font.label = list(size = 14, face = "bold")
      )
      
      plot_file <- file.path(output_dir, sprintf("%s_purity_diagnostics.pdf", tag))
      ggsave(plot_file, combined_plot, width = 14, height = 10, dpi = 300)
      message("[INFO] Diagnostic plots saved to: ", plot_file)
    }
    
    message(sprintf("[INFO] Pipeline end time: %s", Sys.time()))
    
    result
    
  }, error = function(e) {
    message("[ERROR] Pipeline failed with error: ", e$message)
    message("[ERROR] Traceback:")
    print(sys.calls())
    stop(e)
  }, finally = {
    # Close log file
    sink(type = "message")
    sink(type = "output")
    close(log_conn)
    message(sprintf("[INFO] Log file closed: %s", log_file))
  })
}

## Execute Pipeline ----------------------------------------------------------
tryCatch({
  analysis_result <- run_tidyestimate_purity_analysis(
    focus_project = "TCGA-BRCA",
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