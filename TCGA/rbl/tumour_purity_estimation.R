#!/usr/bin/env Rscript

# tumour_purity_estimation.R
# Tumour purity pipeline (tidyestimate) for RBL cohort
#
# INPUT:
#   /home/chu25/data/rbl/rbl_raw_count_matrix.rds   (66 samples x 60660 genes)
#   /home/chu25/data/rbl/final_rbl_ensembl_to_hgnc.tsv
#
# OUTPUT DIR:
#   /home/chu25/data/rbl/results/tumour_purity_results

suppressPackageStartupMessages({
  library(tidyestimate)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tibble)
  library(RColorBrewer)
  library(SummarizedExperiment)
})

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
# 4. Extract project from RBL sample names (GSE/TARGET-style)
# --------------------------------------------------------------
extract_project <- function(samples) {
  # Take everything before first underscore
  prefix <- sub("^([^_]+).*", "\\1", samples)

  project <- dplyr::case_when(
    grepl("^GSE[0-9]+$", prefix) ~ prefix,        # GSE111168, GSE196420, GSE268136, ...
    grepl("^SRP[0-9]+$", prefix) ~ prefix,
    grepl("^TARGET", prefix)      ~ "TARGET",
    TRUE                          ~ "Unknown"
  )

  if (any(project == "Unknown")) {
    bad <- unique(samples[project == "Unknown"])
    stop(
      "[ERROR] Some samples could not be mapped to a project:\n  ",
      paste(bad, collapse = ", "),
      "\nPlease adjust extract_project() to cover these patterns."
    )
  }

  project
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
read_rbl_expr_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("[ERROR] RBL RDS not found: ", path)
  }
  message("[INFO] Loading RBL RDS: ", path)
  obj <- readRDS(path)

  if (inherits(obj, "SummarizedExperiment")) {
    mat <- SummarizedExperiment::assay(obj)
  } else if (is.matrix(obj)) {
    mat <- obj
  } else if (is.data.frame(obj)) {
    mat <- as.matrix(obj)
  } else {
    stop("[ERROR] RBL RDS must be matrix/data.frame/SummarizedExperiment.")
  }

  rn <- rownames(mat)
  cn <- colnames(mat)

  message(sprintf("[INFO] Raw matrix dims: %d x %d", nrow(mat), ncol(mat)))

  # From conversion log: 66 samples x 60660 genes
  # If rows look like GSE* and columns like ENSG*, transpose
  row_is_sample <- !is.null(rn) && any(grepl("^GSE[0-9]+_", rn))
  col_is_ensg   <- !is.null(cn) && any(grepl("^ENSG", cn))

  if (row_is_sample && col_is_ensg) {
    message("[INFO] Detected samples on rows and Ensembl genes on columns → transposing.")
    mat <- t(mat)
  } else if (!any(grepl("^ENSG", rownames(mat))) &&
             any(grepl("^ENSG", colnames(mat)))) {
    message("[WARN] ENSG IDs appear on columns → transposing as fallback.")
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
  expr_matrix_path = "/home/chu25/data/rbl/rbl_raw_count_matrix.rds",
  output_dir       = "/home/chu25/data/rbl/results/tumour_purity_results",
  map_tsv_path     = "/home/chu25/data/rbl/final_rbl_ensembl_to_hgnc.tsv",
  purity_threshold = 0.7
) {
  message("[INFO] Starting RBL tidyestimate purity pipeline...")
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

    # 1) Load and orient expression matrix
    counts <- read_rbl_expr_matrix(expr_matrix_path)
    stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)))
    message("[INFO] Input matrix: ", nrow(counts), " genes × ", ncol(counts), " samples")

    # 2) Project assignment BEFORE filtering
    projects_before <- extract_project(colnames(counts))
    project_before_tab <- table(projects_before)
    message("[INFO] Project counts BEFORE filtering:")
    print(project_before_tab)

    # 3) Ensembl → HGNC
    counts_hgnc <- map_ensembl_to_hgnc(counts, map_tsv_path)

    # 4) Manual common-gene filter (avoid filter_common_genes bug)
    message("[INFO] Selecting common ESTIMATE genes (manual filter)...")
    common_df      <- tidyestimate::common_genes
    common_symbols <- unique(common_df$hgnc_symbol)
    keep_genes     <- intersect(rownames(counts_hgnc), common_symbols)

    if (length(keep_genes) < 50) {
      warning(sprintf(
        "[WARN] Only %d common ESTIMATE genes found; proceeding with ALL HGNC genes.",
        length(keep_genes)
      ))
      counts_sig <- counts_hgnc
    } else {
      message(sprintf("[INFO] Using %d common genes for ESTIMATE scoring.", length(keep_genes)))
      counts_sig <- counts_hgnc[keep_genes, , drop = FALSE]
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

    # 9) Project after filtering
    project_after_tab <- table(extract_project(kept_samples))
    message("[INFO] Project counts AFTER filtering:")
    print(project_after_tab)

    # ----------------------------------------------------------
    # 10a. Bar plot BEFORE vs AFTER (same axis)
    # ----------------------------------------------------------
    message("[INFO] Generating bar plot (sample retention)...")
    total_before <- sum(project_before_tab)
    total_after  <- sum(project_after_tab)

    proj_levels <- sort(unique(c(names(project_before_tab), names(project_after_tab))))
    project_before_tab <- project_before_tab[proj_levels]
    project_after_tab  <- project_after_tab[proj_levels]
    project_before_tab[is.na(project_before_tab)] <- 0
    project_after_tab[is.na(project_after_tab)]   <- 0

    df_plot <- data.frame(
      Project = rep(proj_levels, 2),
      Count   = c(as.numeric(project_before_tab), as.numeric(project_after_tab)),
      Stage   = rep(c("Before Filtering", "After Filtering"), each = length(proj_levels))
    ) %>%
      mutate(
        Project = factor(Project, levels = proj_levels),
        Stage   = factor(Stage, levels = c("Before Filtering", "After Filtering"))
      )

    n_proj <- length(proj_levels)
    palette <- RColorBrewer::brewer.pal(max(3, min(8, n_proj)), "Set2")[seq_len(n_proj)]
    project_colors <- setNames(palette, proj_levels)

    fill_vals <- c(
      setNames(rep("#D3D3D3", n_proj),
               paste("Before Filtering", proj_levels, sep = ".")),
      setNames(project_colors,
               paste("After Filtering", proj_levels, sep = "."))
    )

    p_bar <- ggplot(df_plot, aes(x = Project, y = Count,
                                 fill = interaction(Stage, Project))) +
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
      scale_fill_manual(values = fill_vals, name = NULL) +
      labs(
        title    = "RBL: Sample Retention by Project and Purity Filtering",
        subtitle = sprintf("Purity >= %.1f | n = %d -> %d",
                           purity_threshold, total_before, total_after),
        x = "Project",
        y = "Number of Samples"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title        = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle     = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.title.x      = element_text(face = "bold"),
        axis.title.y      = element_text(face = "bold"),
        legend.position   = "top",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 11),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.border      = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )

    bar_pdf <- file.path(output_dir, "rbl_sample_retention_by_project_same_axis.pdf")
    ggsave(bar_pdf, p_bar, width = 10, height = 7, dpi = 300)
    message("[SUCCESS] Bar plot saved: ", bar_pdf)

    # ----------------------------------------------------------
    # 10b. Purity diagnostics (Stromal/Immune + histogram)
    # ----------------------------------------------------------
    message("[INFO] Generating purity diagnostics plots...")

    plot_data <- scores %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      mutate(
        Project     = extract_project(Sample),
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
        title = "RBL: Tumour Purity vs Stromal Score",
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
        title = "RBL: Tumour Purity vs Immune Score",
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

    # C: Purity distribution by project
    p3 <- ggplot(plot_data, aes(x = tumourPurity, fill = Project)) +
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
      scale_fill_manual(values = project_colors, name = "Project") +
      labs(
        title = "RBL: Distribution of Tumour Purity by Project",
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
    ggsave(diag_pdf, combined, width = 14, height = 10, dpi = 300)
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
  run_rbl_tidyestimate_purity(
    expr_matrix_path = "/home/chu25/data/rbl/rbl_raw_count_matrix.rds",
    output_dir       = "/home/chu25/data/rbl/results/tumour_purity_results",
    map_tsv_path     = "/home/chu25/data/rbl/final_rbl_ensembl_to_hgnc.tsv",
    purity_threshold = 0.7
  )
}