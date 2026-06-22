#!/usr/bin/env Rscript

# rbl_tumour_purity_analysis.R
# Tumour purity pipeline (tidyestimate) for RBL cohort
#
# INPUT DEFAULTS:
#   data/rbl/count_data/rbl_tumour_count.rds
#   data/rbl/count_data/rbl_tumour_sample_metadata.csv
#   data/rbl/count_data/rbl_ensembl_to_hgnc.tsv
#
# OUTPUT DIR:
#   results/tumour_purity_analysis/rbl

suppressPackageStartupMessages({
  library(tidyestimate)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tibble)
  library(RColorBrewer)
})

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE))
  }
  normalizePath(
    "preprocessing_and_quality_control/rbl/tumour_purity_analysis/rbl_tumour_purity_analysis.R",
    mustWork = FALSE
  )
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIR <- dirname(SCRIPT_PATH)
WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
PIPELINE_DIR <- normalizePath(file.path(WORKFLOW_DIR, "..", ".."), mustWork = FALSE)
RBL_DATA_ROOT <- file.path(PIPELINE_DIR, "data", "rbl")
RBL_COUNT_DIR <- file.path(RBL_DATA_ROOT, "count_data")

default_count_matrix <- file.path(RBL_COUNT_DIR, "rbl_tumour_count.rds")
default_metadata <- file.path(RBL_COUNT_DIR, "rbl_tumour_sample_metadata.csv")
default_map <- file.path(RBL_COUNT_DIR, "rbl_ensembl_to_hgnc.tsv")
default_output_dir <- file.path(PIPELINE_DIR, "results", "tumour_purity_analysis", "rbl")

# --------------------------------------------------------------
# 1. Map Ensembl → HGNC (robust, deduplicates symbols)
# --------------------------------------------------------------
map_ensembl_to_hgnc <- function(expr_matrix, map_tsv_path) {
  # expr_matrix: genes x samples, rownames = Ensembl IDs or HGNC

  # If it already looks like HGNC, just upper-case and return
  if (!any(grepl("^ENSG\\d", rownames(expr_matrix), perl = TRUE))) {
    message("[INFO] Row names appear to be HGNC symbols. Skipping ID mapping.")
    rownames(expr_matrix) <- toupper(rownames(expr_matrix))
    return(expr_matrix)
  }

  message("[INFO] Loading Ensembl → HGNC map from: ", map_tsv_path)
  if (!file.exists(map_tsv_path)) {
    stop("[ERROR] Mapping file not found: ", map_tsv_path)
  }

  gene_map <- data.table::fread(
    map_tsv_path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE
  )

  if (!all(c("Ensembl_ID", "HGNC_Symbol") %in% colnames(gene_map))) {
    stop("[ERROR] Mapping file must have columns: Ensembl_ID, HGNC_Symbol")
  }

  # Strip version number
  ensembl_ids <- sub("\\.\\d+$", "", rownames(expr_matrix))
  id_to_symbol <- setNames(gene_map$HGNC_Symbol, gene_map$Ensembl_ID)

  mapped_symbols <- id_to_symbol[ensembl_ids]
  valid <- !is.na(mapped_symbols) & nzchar(mapped_symbols)

  message(sprintf(
    "[INFO] Successfully mapped %d / %d Ensembl IDs (%.1f%%).",
    sum(valid), length(ensembl_ids), 100 * mean(valid)
  ))

  if (sum(valid) == 0) {
    stop("[ERROR] No Ensembl IDs could be mapped to HGNC symbols.")
  }

  # Build explicit data.frame with symbol column
  expr_df <- as.data.frame(expr_matrix[valid, , drop = FALSE])
  expr_df$symbol <- toupper(mapped_symbols[valid])

  # Move symbol to first column
  expr_df <- expr_df %>%
    dplyr::relocate(symbol)

  before <- nrow(expr_df)

  # Collapse duplicate HGNC symbols by summing numeric columns
  expr_mapped <- expr_df %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = where(is.numeric),
        .fns  = ~ sum(.x),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    tibble::column_to_rownames("symbol") %>%
    as.matrix()

  if (nrow(expr_mapped) < before) {
    message(sprintf(
      "[INFO] Collapsed %d → %d genes (duplicate HGNC symbols)",
      before, nrow(expr_mapped)
    ))
  } else {
    message("[INFO] No duplicate HGNC symbols to collapse.")
  }

  expr_mapped
}

# --------------------------------------------------------------
# 2. ESTIMATE scoring helper (estimate_score vs estimate_scores)
# --------------------------------------------------------------
get_estimate_scoring_function <- function() {
  exports <- try(getNamespaceExports("tidyestimate"), silent = TRUE)
  if (inherits(exports, "try-error")) {
    stop("[ERROR] tidyestimate not available.")
  }
  if ("estimate_score" %in% exports)  return(tidyestimate::estimate_score)
  if ("estimate_scores" %in% exports) return(tidyestimate::estimate_scores)
  stop("[ERROR] Neither 'estimate_score' nor 'estimate_scores' found in tidyestimate.")
}

# --------------------------------------------------------------
# 3. Standardize ESTIMATE score columns (+ purity)
# --------------------------------------------------------------
standardize_estimate_scores <- function(scores_df) {
  scores_df <- as.data.frame(scores_df)

  if ("sample" %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df$sample <- NULL
  }
  if (is.matrix(scores_df)) scores_df <- as.data.frame(scores_df)

  col_lower <- tolower(colnames(scores_df))
  new_names <- colnames(scores_df)
  new_names[grepl("strom", col_lower)] <- "StromalScore"
  new_names[grepl("\\bimm", col_lower)] <- "ImmuneScore"
  new_names[grepl("estimate", col_lower) | grepl("\\best\\b", col_lower)] <- "ESTIMATEScore"
  colnames(scores_df) <- new_names

  # Handle 3-col outputs with no names
  if (ncol(scores_df) == 3 && (is.null(colnames(scores_df)) || all(colnames(scores_df) == ""))) {
    colnames(scores_df) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
    message("[INFO] Assigned default ESTIMATE column names.")
  }

  # Derive purity if missing
  if (!"tumourPurity" %in% colnames(scores_df) && "ESTIMATEScore" %in% colnames(scores_df)) {
    scores_df$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)
    scores_df$tumourPurity <- pmin(1, pmax(0, scores_df$tumourPurity))
    message("[INFO] Derived tumourPurity from ESTIMATEScore.")
  }

  req <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "tumourPurity")
  miss <- setdiff(req, colnames(scores_df))
  if (length(miss) > 0) {
    warning("[WARN] Missing expected columns: ", paste(miss, collapse = ", "))
  }

  scores_df
}

# --------------------------------------------------------------
# 4. Sample metadata and cohort helpers
# --------------------------------------------------------------
load_sample_metadata <- function(meta_csv_path) {
  if (!file.exists(meta_csv_path)) {
    stop("[ERROR] Sample metadata file not found: ", meta_csv_path,
         "\nRun the RBL count merge first.")
  }

  meta <- data.table::fread(meta_csv_path, sep = ",", header = TRUE, data.table = FALSE)

  if (!"cohort" %in% colnames(meta)) {
    stop("[ERROR] Metadata must contain a cohort column.")
  }

  if (!"sample" %in% colnames(meta)) {
    alternatives <- intersect(c("sample_id", "run", "Run", "srr", "SRR"), colnames(meta))
    if (length(alternatives) == 0) {
      stop("[ERROR] Metadata must contain a sample column or one of: sample_id, run, Run, srr, SRR.")
    }
    meta$sample <- meta[[alternatives[1]]]
    message("[INFO] Using metadata column as sample ID: ", alternatives[1])
  }

  meta <- meta[!is.na(meta$sample) & meta$sample != "" &
               !is.na(meta$cohort) & meta$cohort != "", , drop = FALSE]

  meta$sample <- as.character(meta$sample)
  meta$cohort <- sub("_primary_tumours$", "", as.character(meta$cohort))

  meta
}

cohort_for_samples <- function(samples, metadata) {
  sample_to_cohort <- setNames(metadata$cohort, metadata$sample)
  cohort <- sample_to_cohort[samples]

  missing <- is.na(cohort)
  if (any(missing)) {
    stop(
      "[ERROR] Some expression samples have no cohort mapping in metadata:\n  ",
      paste(unique(samples[missing]), collapse = ", "),
      "\nCheck that rbl_tumour_sample_metadata.csv matches the count matrix columns."
    )
  }

  unname(cohort)
}

make_retention_summary <- function(all_samples, kept_samples, metadata) {
  before <- data.frame(
    Sample = all_samples,
    Cohort = cohort_for_samples(all_samples, metadata),
    Stage = rep("Before Filtering", length(all_samples)),
    stringsAsFactors = FALSE
  )

  after <- data.frame(
    Sample = kept_samples,
    Cohort = cohort_for_samples(kept_samples, metadata),
    Stage = rep("After Filtering", length(kept_samples)),
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(before, after) %>%
    dplyr::count(Stage, Cohort, name = "Count")
}

# --------------------------------------------------------------
# 5. Filter samples by purity  (ROBUST VERSION)
# --------------------------------------------------------------
filter_samples_by_purity <- function(expr_matrix, purity_scores, threshold = 0.7) {
  purity_scores <- as.data.frame(purity_scores)

  if (!"tumourPurity" %in% colnames(purity_scores)) {
    message("[WARN] tumourPurity not found in scores, skipping filtering.")
    return(expr_matrix)
  }

  common <- intersect(colnames(expr_matrix), rownames(purity_scores))
  message(sprintf("[INFO] %d overlapping samples between expression and purity scores.", length(common)))

  if (length(common) == 0) {
    message("[WARN] No overlap between expression samples and purity scores; skipping filtering.")
    return(expr_matrix)
  }

  # Pull tumourPurity for overlapping samples and keep names = sample IDs
  purity_vec <- purity_scores[common, "tumourPurity", drop = TRUE]
  names(purity_vec) <- common

  message("[INFO] Tumour purity summary (overlapping samples):")
  print(summary(purity_vec))

  keep <- !is.na(purity_vec) & purity_vec >= threshold
  kept <- common[keep]

  n_total <- length(common)
  n_na    <- sum(is.na(purity_vec))
  n_low   <- sum(!is.na(purity_vec) & purity_vec < threshold)
  n_kept  <- length(kept)
  pct_kept <- if (n_total > 0) 100 * n_kept / n_total else NA_real_

  msg <- paste0(
    "[INFO] Filtering summary:\n",
    "  Overlapping samples: ", n_total, "\n",
    "  NA purity: ", n_na, "\n",
    "  Low purity (< ", sprintf("%.2f", threshold), "): ", n_low, "\n",
    "  Kept (>= ", sprintf("%.2f", threshold), "): ", n_kept,
    " (", ifelse(is.na(pct_kept), "NA", sprintf("%.1f", pct_kept)), "%)"
  )
  message(msg)

  expr_matrix[, kept, drop = FALSE]
}

# --------------------------------------------------------------
# 6. Read RBL RDS and ensure genes x samples orientation
# --------------------------------------------------------------
read_rbl_expr_matrix <- function(path, metadata = NULL) {
  if (!file.exists(path)) {
    stop("[ERROR] RBL RDS not found: ", path)
  }
  message("[INFO] Loading RBL RDS: ", path)
  obj <- readRDS(path)

  if (inherits(obj, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("[ERROR] RBL RDS is a SummarizedExperiment, but the optional SummarizedExperiment package is not installed.")
    }
    mat <- SummarizedExperiment::assay(obj)
  } else if (is.matrix(obj)) {
    mat <- obj
  } else if (is.data.frame(obj)) {
    mat <- as.matrix(obj)
  } else {
    stop("[ERROR] RBL RDS must be a matrix, data.frame, or optional SummarizedExperiment object.")
  }

  rn <- rownames(mat)
  cn <- colnames(mat)

  message(sprintf("[INFO] Raw matrix dims: %d x %d", nrow(mat), ncol(mat)))

  metadata_samples <- if (!is.null(metadata) && "sample" %in% colnames(metadata)) {
    as.character(metadata$sample)
  } else {
    character(0)
  }
  row_matches_metadata <- !is.null(rn) && any(rn %in% metadata_samples)
  col_matches_metadata <- !is.null(cn) && any(cn %in% metadata_samples)
  row_is_ensg <- !is.null(rn) && any(grepl("^ENSG", rn))
  col_is_ensg <- !is.null(cn) && any(grepl("^ENSG", cn))

  if (row_matches_metadata && !col_matches_metadata) {
    message("[INFO] Detected metadata samples on rows; transposing to genes x samples.")
    mat <- t(mat)
  } else if (col_matches_metadata) {
    message("[INFO] Detected metadata samples on columns; using matrix as-is.")
  } else if (!row_is_ensg && col_is_ensg) {
    message("[WARN] ENSG IDs appear on columns; transposing as fallback.")
    mat <- t(mat)
  } else {
    message("[INFO] Using matrix as-is (assuming genes on rows).")
  }

  message(sprintf("[INFO] Final orientation: %d genes x %d samples",
                  nrow(mat), ncol(mat)))
  mat
}

# --------------------------------------------------------------
# 7. CORE PIPELINE: ESTIMATE + purity filter + plots + outputs
# --------------------------------------------------------------
run_rbl_tidyestimate_purity <- function(
  expr_matrix_path = Sys.getenv("SE_PATH", default_count_matrix),
  meta_csv_path    = Sys.getenv("META_CSV", default_metadata),
  map_tsv_path     = Sys.getenv("MAP_TSV", default_map),
  output_dir       = Sys.getenv("OUTPUT_DIR", default_output_dir),
  purity_threshold = as.numeric(Sys.getenv("PURITY_THRESHOLD", "0.7"))
) {
  message("[INFO] Starting RBL tidyestimate purity pipeline...")
  if (is.na(purity_threshold)) {
    stop("[ERROR] PURITY_THRESHOLD must be numeric.")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  
  # Set up log file
  log_file <- file.path(output_dir, paste0("rbl_purity_pipeline_", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                            ".log"))
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "output", split = TRUE)
  sink(log_conn, type = "message")
  
  message(sprintf("[INFO] Log file created: %s", log_file))
  message(sprintf("[INFO] Pipeline start time: %s", Sys.time()))
  
  # Wrap pipeline in tryCatch to ensure log cleanup
  tryCatch({

    # 1) Load metadata, then load and orient expression matrix
    metadata <- load_sample_metadata(meta_csv_path)
    message(sprintf("[INFO] Loaded sample metadata: %d rows", nrow(metadata)))

    counts <- read_rbl_expr_matrix(expr_matrix_path, metadata)
    stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)))
    message("[INFO] Input matrix: ", nrow(counts), " genes × ", ncol(counts), " samples")

    cohort_for_samples(colnames(counts), metadata)
    message("[INFO] Metadata cohort mapping found for all expression samples.")

    # 3) Ensembl → HGNC
    counts_hgnc <- map_ensembl_to_hgnc(counts, map_tsv_path)

    # 4) Common-gene filter using the exported tidyestimate API.
    message("[INFO] Selecting common ESTIMATE genes with tidyestimate::filter_common_genes...")
    counts_for_filter <- as.data.frame(counts_hgnc) %>%
      tibble::rownames_to_column("hgnc_symbol")

    counts_sig_df <- tidyestimate::filter_common_genes(
      counts_for_filter,
      id = "hgnc_symbol",
      tidy = TRUE,
      tell_missing = FALSE,
      find_alias = FALSE
    )

    gene_col <- intersect(c("hgnc_symbol", "entrezgene_id", "symbol"), colnames(counts_sig_df))[1]
    if (is.na(gene_col)) {
      gene_col <- colnames(counts_sig_df)[1]
    }

    counts_sig <- counts_sig_df %>%
      as.data.frame() %>%
      tibble::column_to_rownames(gene_col) %>%
      as.matrix()

    if (nrow(counts_sig) < 50) {
      warning(sprintf(
        "[WARN] Only %d common ESTIMATE genes found; proceeding with ALL HGNC genes.",
        nrow(counts_sig)
      ))
      counts_sig <- counts_hgnc
    } else {
      message(sprintf("[INFO] Using %d common genes for ESTIMATE scoring.", nrow(counts_sig)))
    }

    # 5) Prepare for tidyestimate (HGNC symbols in first column)
    df_for_estimate <- as.data.frame(counts_sig) %>%
      tibble::rownames_to_column("hgnc_symbol")

    # 6) Run ESTIMATE
    score_fn <- get_estimate_scoring_function()
    message("[INFO] Running ESTIMATE (tidyestimate)...")
    scores_raw <- score_fn(df_for_estimate, is_affymetrix = FALSE)
    scores     <- standardize_estimate_scores(scores_raw)

    # --- FORCE SAMPLE NAME ALIGNMENT (critical for filtering) ---
    if (nrow(scores) == ncol(counts_sig)) {
      message("[INFO] Forcing ESTIMATE score rownames to expression sample IDs.")
      rownames(scores) <- colnames(counts_sig)
    } else {
      warning(sprintf(
        "[WARN] nrow(scores) = %d but ncol(counts_sig) = %d; skipping forced alignment.",
        nrow(scores), ncol(counts_sig)
      ))
    }

    overlap <- intersect(colnames(counts_hgnc), rownames(scores))
    message(sprintf("[INFO] Overlap between counts and ESTIMATE scores: %d samples", length(overlap)))

    message(sprintf("[INFO] Scores: %d samples × %d metrics",
                    nrow(scores), ncol(scores)))

    # 7) Save scores
    scores_csv <- file.path(output_dir, "rbl_estimate_scores.csv")
    write.csv(scores, scores_csv, row.names = TRUE)
    message("[INFO] Scores saved to: ", scores_csv)

    # 8) Purity-based filtering (HGNC + Ensembl)
    counts_filtered_hgnc <- filter_samples_by_purity(counts_hgnc, scores, purity_threshold)
    kept_samples <- colnames(counts_filtered_hgnc)

    filtered_hgnc_rds <- file.path(
      output_dir,
      sprintf("rbl_filtered_hgnc_counts_purity%.1f.rds", purity_threshold)
    )
    saveRDS(counts_filtered_hgnc, filtered_hgnc_rds)
    message("[INFO] Filtered HGNC counts saved to: ", filtered_hgnc_rds)

    # Filter original Ensembl matrix with the same sample list
    ensembl_filtered <- counts[, kept_samples, drop = FALSE]
    ensembl_out <- file.path(
      output_dir,
      sprintf("rbl_filtered_count_matrix_ensembl_purity%.1f.rds", purity_threshold)
    )
    saveRDS(ensembl_filtered, ensembl_out)
    message("[INFO] Filtered Ensembl matrix saved to: ", ensembl_out)

    # 9) Cohort retention summary
    retention_summary <- make_retention_summary(colnames(counts), kept_samples, metadata)
    message("[INFO] Cohort retention counts:")
    print(retention_summary)

    # ----------------------------------------------------------
    # 10a. Bar plot BEFORE vs AFTER (same axis)
    # ----------------------------------------------------------
    message("[INFO] Generating bar plot (sample retention)...")
    total_before <- sum(retention_summary$Count[retention_summary$Stage == "Before Filtering"])
    total_after <- sum(retention_summary$Count[retention_summary$Stage == "After Filtering"])

    cohort_levels <- sort(unique(retention_summary$Cohort))
    stage_levels <- c("Before Filtering", "After Filtering")
    retention_plot <- merge(
      expand.grid(Stage = stage_levels, Cohort = cohort_levels, stringsAsFactors = FALSE),
      retention_summary,
      by = c("Stage", "Cohort"),
      all.x = TRUE,
      sort = FALSE
    )
    retention_plot$Count[is.na(retention_plot$Count)] <- 0
    retention_plot <- retention_plot %>%
      mutate(
        Cohort = factor(Cohort, levels = cohort_levels),
        Stage = factor(Stage, levels = stage_levels)
      )

    n_cohorts <- length(cohort_levels)
    cohort_palette <- RColorBrewer::brewer.pal(max(3, min(8, n_cohorts)), "Set2")
    if (n_cohorts > length(cohort_palette)) {
      cohort_palette <- grDevices::colorRampPalette(cohort_palette)(n_cohorts)
    }
    cohort_colors <- setNames(cohort_palette[seq_len(n_cohorts)], cohort_levels)

    p_bar <- ggplot(retention_plot, aes(x = Cohort, y = Count, fill = Stage)) +
      geom_col(position = position_dodge(width = 0.8),
               width = 0.7, color = "white", linewidth = 0.6) +
      geom_text(
        aes(label = Count, group = Stage),
        position = position_dodge(width = 0.8),
        vjust    = -0.6,
        size     = 4.2,
        fontface = "bold",
        color    = "black"
      ) +
      scale_fill_manual(
        values = c("Before Filtering" = "#D3D3D3", "After Filtering" = "#66C2A5"),
        name = NULL
      ) +
      labs(
        title    = "RBL Sample Retention by Cohort and Purity Filtering",
        subtitle = sprintf("Purity >= %.1f | n = %d -> %d",
                           purity_threshold, total_before, total_after),
        x = "Cohort",
        y = "Number of Samples"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title        = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle     = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.title.x      = element_text(face = "bold"),
        axis.title.y      = element_text(face = "bold"),
        legend.position   = "top",
        legend.text       = element_text(size = 11),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.border      = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )

    bar_pdf <- file.path(output_dir, "rbl_sample_retention_by_cohort.pdf")
    ggsave(bar_pdf, p_bar, width = 10, height = 7, dpi = 300, device = cairo_pdf)
    message("[SUCCESS] Bar plot saved: ", bar_pdf)

    # ----------------------------------------------------------
    # 10b. Purity diagnostics (Stromal/Immune + histogram)
    # ----------------------------------------------------------
    message("[INFO] Generating purity diagnostics plots...")

    plot_data <- scores %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      mutate(
        Cohort      = cohort_for_samples(Sample, metadata),
        PurityGroup = cut(
          tumourPurity,
          breaks = c(-Inf, purity_threshold, Inf),
          labels = c(paste("<", purity_threshold), paste(">=", purity_threshold)),
          include.lowest = TRUE
        )
      )

    pal <- RColorBrewer::brewer.pal(3, "Set2")[c(2, 1)]  # two colours

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
        title = "Tumour Purity vs Stromal Score",
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
        title = "Tumour Purity vs Immune Score",
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

    # C: Purity distribution by cohort
    p3 <- ggplot(plot_data, aes(x = tumourPurity, fill = Cohort)) +
      geom_histogram(bins = 40, alpha = 0.8, position = "stack", colour = NA) +
      geom_vline(xintercept = purity_threshold, linetype = "dashed",
                 colour = "firebrick", linewidth = 1) +
      annotate(
        "text",
        x = purity_threshold + 0.02,
        y = Inf,
        vjust = 1.5,
        hjust = 0,
        label = paste("Threshold =", purity_threshold),
        colour = "firebrick",
        size = 4
      ) +
      scale_fill_manual(values = cohort_colors, name = "Cohort") +
      labs(
        title = "Distribution of Tumour Purity by Cohort",
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

    combined <- ggpubr::ggarrange(
      p1, p2, p3,
      ncol = 2, nrow = 2,
      labels = c("A", "B", "C"),
      font.label = list(size = 14, face = "bold")
    )

    diag_pdf <- file.path(output_dir, "rbl_purity_diagnostics.pdf")
    ggsave(diag_pdf, combined, width = 14, height = 10, dpi = 300, device = cairo_pdf)
    message("[SUCCESS] Diagnostics saved: ", diag_pdf)

    message(sprintf("[INFO] Pipeline end time: %s", Sys.time()))
    message("[SUCCESS] RBL tidyestimate purity pipeline finished.")

    result <- list(
      scores                  = scores,
      counts_filtered_hgnc    = counts_filtered_hgnc,
      counts_filtered_ensembl = ensembl_filtered
    )
    
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

# --------------------------------------------------------------
# 8. MAIN (for direct Rscript use)
# --------------------------------------------------------------
if (sys.nframe() == 0) {
  run_rbl_tidyestimate_purity()
}
