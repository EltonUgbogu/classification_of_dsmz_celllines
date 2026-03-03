#!/usr/bin/env Rscript
# hema_tumour_purity_estimation.R
# Robust tidyestimate-based tumour purity pipeline for the HEMA cohort.
# 
# Key fixes (v2):
# - Correctly strips Ensembl version suffixes from BOTH expression matrix AND mapping file
# - Adds sanity checks for ESTIMATE signature gene overlap
# - Fails early with clear error if mapping is insufficient

suppressPackageStartupMessages({
  library(tidyestimate)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  library(data.table)
  library(tibble)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

# Ensure matrix is genes × samples (required for map_ensembl_to_hgnc and ESTIMATE).
# Detects orientation by ENSG pattern on row/col names; transposes if needed.
ensure_genes_by_samples <- function(mat) {
  stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)))

  row_is_ensg <- any(grepl("^ENSG", rownames(mat)))
  col_is_ensg <- any(grepl("^ENSG", colnames(mat)))

  if (!row_is_ensg && col_is_ensg) {
    message("[INFO] Detected samples×genes; transposing to genes×samples.")
    mat <- t(mat)
  }
  mat
}

# Helper to strip Ensembl version suffix: ENSG00000141510.13 → ENSG00000141510
strip_ensembl_version <- function(ids) {

  sub("\\.\\d+$", "", ids)
}

# Map Ensembl → HGNC and collapse duplicates
map_ensembl_to_hgnc <- function(expr_matrix, map_tsv_path) {
  gene_ids <- rownames(expr_matrix)
  
  # Check if already HGNC symbols (no ENSG pattern)
  if (!any(grepl("^ENSG\\d", gene_ids, perl = TRUE))) {
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

  # Strip version from BOTH the expression matrix AND the mapping file
  ensembl_ids_expr <- strip_ensembl_version(gene_ids)
  ensembl_ids_map  <- strip_ensembl_version(gene_map$Ensembl_ID)
  
  message("[INFO] Stripped Ensembl version suffixes from expression matrix and mapping file.")
  message("[INFO] Expression matrix gene IDs (first 5): ", 
          paste(head(gene_ids, 5), collapse = ", "))
  message("[INFO] After stripping (first 5): ", 
          paste(head(ensembl_ids_expr, 5), collapse = ", "))
  
  # Build lookup: bare Ensembl ID → HGNC symbol
  # Handle duplicates by keeping first occurrence
  gene_map$ensembl_bare <- ensembl_ids_map
  gene_map_unique <- gene_map[!duplicated(gene_map$ensembl_bare), ]
  
  id_to_symbol <- setNames(gene_map_unique$HGNC_Symbol, gene_map_unique$ensembl_bare)
  
  # Map expression matrix genes
  mapped_symbols <- id_to_symbol[ensembl_ids_expr]
  valid <- !is.na(mapped_symbols) & nzchar(mapped_symbols)

  n_mapped <- sum(valid)
  n_total <- length(ensembl_ids_expr)
  pct_mapped <- round(100 * n_mapped / n_total, 1)
  
  message(sprintf(
    "[INFO] Successfully mapped %d / %d Ensembl IDs (%.1f%%).",
    n_mapped, n_total, pct_mapped
  ))

  if (n_mapped == 0) {
    # Debug: show some IDs that failed to map
    sample_expr <- head(ensembl_ids_expr, 10)
    sample_map <- head(names(id_to_symbol), 10)
    message("[DEBUG] Sample expression IDs: ", paste(sample_expr, collapse = ", "))
    message("[DEBUG] Sample mapping IDs: ", paste(sample_map, collapse = ", "))
    stop("[ERROR] No Ensembl IDs could be mapped to HGNC symbols. Check ID formats.")
  }
  
  if (pct_mapped < 10) {
    warning(sprintf("[WARN] Only %.1f%% of genes mapped. Results may be unreliable.", pct_mapped))
  }

  expr_df <- as.data.frame(expr_matrix[valid, , drop = FALSE])
  expr_df$symbol <- toupper(mapped_symbols[valid])
  expr_df <- expr_df %>% dplyr::relocate(symbol)

  before <- nrow(expr_df)

  # Collapse duplicate HGNC symbols by summing counts
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
  }
  
  message(sprintf("[INFO] Final HGNC matrix: %d genes × %d samples", 
                  nrow(expr_mapped), ncol(expr_mapped)))

  expr_mapped
}

# Check overlap with ESTIMATE common genes
check_estimate_gene_overlap <- function(expr_matrix, min_genes = 1000) {
  hgnc_genes <- rownames(expr_matrix)
  
  # Get ESTIMATE common genes (the gene universe used by tidyestimate)
  common_genes_df <- tidyestimate::common_genes
  estimate_genes <- unique(common_genes_df$hgnc_symbol)
  
  n_overlap <- sum(hgnc_genes %in% estimate_genes)
  pct_overlap <- round(100 * n_overlap / length(estimate_genes), 1)
  
  message(sprintf(
    "[INFO] ESTIMATE gene overlap: %d / %d common genes found (%.1f%%)",
    n_overlap, length(estimate_genes), pct_overlap
  ))
  
  if (n_overlap < min_genes) {
    # Show some genes that are present vs expected
    present <- intersect(hgnc_genes, estimate_genes)
    message("[DEBUG] Sample present genes (first 20): ", 
            paste(head(present, 20), collapse = ", "))
    missing <- setdiff(estimate_genes, hgnc_genes)
    message("[DEBUG] Sample missing genes (first 20): ", 
            paste(head(missing, 20), collapse = ", "))
    
    stop(sprintf(
      "[ERROR] Insufficient ESTIMATE gene overlap (need >= %d, found %d).\n" %+%
      "  Check your Ensembl → HGNC mapping or use an HGNC-based expression matrix.",
      min_genes, n_overlap
    ))
  }
  
  invisible(list(overlap = n_overlap, pct = pct_overlap))
}

# String concatenation helper
`%+%` <- function(a, b) paste0(a, b)

get_estimate_scoring_function <- function() {
  exports <- try(getNamespaceExports("tidyestimate"), silent = TRUE)
  if (inherits(exports, "try-error")) {
    stop("[ERROR] tidyestimate not available.")
  }
  if ("estimate_score" %in% exports)  return(tidyestimate::estimate_score)
  if ("estimate_scores" %in% exports) return(tidyestimate::estimate_scores)
  stop("[ERROR] Neither 'estimate_score' nor 'estimate_scores' found in tidyestimate.")
}

standardize_estimate_scores <- function(scores_df) {
  scores_df <- as.data.frame(scores_df)
  if ("sample" %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df$sample <- NULL
  }
  col_lower <- tolower(colnames(scores_df))
  new_names <- colnames(scores_df)
  new_names[grepl("strom", col_lower)] <- "StromalScore"
  new_names[grepl("\\bimm", col_lower)] <- "ImmuneScore"
  new_names[grepl("estimate", col_lower) | grepl("\\best\\b", col_lower)] <- "ESTIMATEScore"
  colnames(scores_df) <- new_names

  if (ncol(scores_df) == 3 && (is.null(colnames(scores_df)) || all(colnames(scores_df) == ""))) {
    colnames(scores_df) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
    message("[INFO] Assigned default ESTIMATE column names.")
  }

  if (!"tumourPurity" %in% colnames(scores_df) && "ESTIMATEScore" %in% colnames(scores_df)) {
    scores_df$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)
    scores_df$tumourPurity <- pmin(1, pmax(0, scores_df$tumourPurity))
    message("[INFO] Derived tumourPurity from ESTIMATEScore.")
  }
  scores_df
}

filter_samples_by_purity <- function(expr_matrix, scores_df, threshold = 0.7) {
  scores_df <- as.data.frame(scores_df)
  if (!"tumourPurity" %in% colnames(scores_df)) {
    message("[WARN] tumourPurity not found; skipping filtering.")
    return(expr_matrix)
  }

  common <- intersect(colnames(expr_matrix), rownames(scores_df))
  message(sprintf("[INFO] %d overlapping samples between expression and purity scores.", length(common)))
  if (length(common) == 0) {
    message("[WARN] No overlap; skipping filtering.")
    return(expr_matrix)
  }

  purity_vec <- scores_df[common, "tumourPurity", drop = TRUE]
  names(purity_vec) <- common

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
    "  Low purity (< ", sprintf('%.2f', threshold), "): ", n_low, "\n",
    "  Kept (>= ", sprintf('%.2f', threshold), "): ", n_kept,
    " (", ifelse(is.na(pct_kept), "NA", sprintf('%.1f', pct_kept)), "%)"
  )
  message(msg)

  expr_matrix[, kept, drop = FALSE]
}

read_hema_expr <- function(path) {
  if (!file.exists(path)) stop("[ERROR] HEMA RDS not found: ", path)
  message("[INFO] Loading HEMA RDS: ", path)
  mat <- readRDS(path)
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      mat <- as.matrix(mat)
    } else {
      stop("[ERROR] RDS must be a matrix or data.frame.")
    }
  }
  message(sprintf("[INFO] Loaded matrix: %d × %d", nrow(mat), ncol(mat)))

  mat <- ensure_genes_by_samples(mat)
  message(sprintf("[INFO] Matrix dims now: %d genes × %d samples", nrow(mat), ncol(mat)))

  mat
}

run_hema_purity <- function(expr_matrix_path, map_tsv_path, output_dir, purity_threshold = 0.7) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  
  # Set up log file
  log_file <- file.path(output_dir, paste0("hema_purity_pipeline_", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                            ".log"))
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "output", split = TRUE)
  sink(log_conn, type = "message")
  
  message(sprintf("[INFO] Log file created: %s", log_file))
  message(sprintf("[INFO] Pipeline start time: %s", Sys.time()))
  
  # Wrap pipeline in tryCatch to ensure log cleanup
  tryCatch({

    counts <- read_hema_expr(expr_matrix_path)
    stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)))

    message("[INFO] Input matrix: ", nrow(counts), " genes × ", ncol(counts), " samples")

    # Map Ensembl → HGNC (with version stripping fix)
    counts_hgnc <- map_ensembl_to_hgnc(counts, map_tsv_path)
    
    # Sanity check: ensure sufficient ESTIMATE common gene overlap
    message("[INFO] Checking ESTIMATE gene overlap...")
    check_estimate_gene_overlap(counts_hgnc, min_genes = 1000)

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

    df_for_estimate <- as.data.frame(counts_sig) %>%
      tibble::rownames_to_column("hgnc_symbol")

    score_fn <- get_estimate_scoring_function()
    message("[INFO] Running ESTIMATE (tidyestimate)...")
    scores_raw <- score_fn(df_for_estimate, is_affymetrix = FALSE)
    scores     <- standardize_estimate_scores(scores_raw)

    if (nrow(scores) == ncol(counts_sig)) {
      rownames(scores) <- colnames(counts_sig)
    }
    
    # Check for all-NA purity (should not happen now with sanity check)
    if (all(is.na(scores$tumourPurity))) {
      stop("[ERROR] All tumourPurity values are NA. This should not happen after the sanity check.")
    }

    scores_csv <- file.path(output_dir, "estimate_scores.csv")
    write.csv(scores, scores_csv, row.names = TRUE)
    message("[INFO] Scores saved to: ", scores_csv)

    counts_filtered_hgnc <- filter_samples_by_purity(counts_hgnc, scores, purity_threshold)
    kept_samples <- colnames(counts_filtered_hgnc)

    filtered_hgnc_rds <- file.path(
      output_dir,
      sprintf("filtered_hgnc_purity%.1f.rds", purity_threshold)
    )
    saveRDS(counts_filtered_hgnc, filtered_hgnc_rds)
    message("[INFO] Filtered HGNC counts saved to: ", filtered_hgnc_rds)

    ensembl_filtered <- counts[, kept_samples, drop = FALSE]
    ensembl_out <- file.path(
      output_dir,
      sprintf("filtered_ensembl_purity%.1f.rds", purity_threshold)
    )
    saveRDS(ensembl_filtered, ensembl_out)
    message("[INFO] Filtered Ensembl matrix saved to: ", ensembl_out)

    # ------------------------------------------------------------------
    # Purity histogram
    # ------------------------------------------------------------------
    if ("tumourPurity" %in% colnames(scores)) {
      hist_path <- file.path(output_dir, "purity_histogram.pdf")
      valid_purity <- scores$tumourPurity
      valid_purity <- valid_purity[!is.na(valid_purity)]
      if (length(valid_purity) == 0L) {
        message("[WARN] No non-NA tumourPurity values; writing placeholder histogram.")
        pdf(hist_path, width = 7, height = 5)
        plot.new()
        text(
          0.5, 0.5,
          labels = sprintf(
            "No valid tumour purity estimates\n(threshold = %.2f)", purity_threshold
          ),
          cex = 1.1
        )
        dev.off()
      } else {
        pdf(hist_path, width = 7, height = 5)
        hist(
          valid_purity,
          breaks = 40,
          main = sprintf("HEMA: Tumour purity (threshold %.2f)", purity_threshold),
          xlab  = "Tumour purity",
          col = "steelblue",
          border = "white"
        )
        abline(v = purity_threshold, col = "red", lwd = 2, lty = 2)
        abline(v = median(valid_purity), col = "darkgreen", lwd = 2, lty = 3)
        legend("topright", 
               legend = c(sprintf("Threshold (%.2f)", purity_threshold),
                         sprintf("Median (%.2f)", median(valid_purity))),
               col = c("red", "darkgreen"),
               lty = c(2, 3),
               lwd = 2)
        dev.off()
        message("[INFO] Purity histogram saved to: ", hist_path)
      }
    }
    
    # N_before → N_after (raw matrix is genes × samples, so samples = ncol)
    N_before <- ncol(counts)
    N_after  <- length(kept_samples)
    message("\n[INFO] === PURITY SUMMARY ===")
    message(sprintf("  [RETENTION] N_before → N_after: %d → %d (purity >= %.2f)", N_before, N_after, purity_threshold))
    message(sprintf("  Total samples (N_before): %d", N_before))
    message(sprintf("  Valid purity estimates: %d", sum(!is.na(scores$tumourPurity))))
    message(sprintf("  Samples >= %.2f purity (N_after): %d", purity_threshold, N_after))
    message(sprintf("  Median purity: %.3f", median(scores$tumourPurity, na.rm = TRUE)))
    message(sprintf("  Mean purity: %.3f", mean(scores$tumourPurity, na.rm = TRUE)))

    message(sprintf("[INFO] Pipeline end time: %s", Sys.time()))

    result <- list(
      scores = scores,
      counts_filtered_hgnc = counts_filtered_hgnc,
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

# ------------------------------------------------------------------
# HEMA-specific cohort-level plots (CLL / BL / MPN)
# ------------------------------------------------------------------
make_hema_cohort_plots <- function(meta_csv, scores_df, purity_threshold, outdir) {
  if (!file.exists(meta_csv)) {
    message("[WARN] Metadata file not found, skipping HEMA cohort plots: ", meta_csv)
    return(invisible(NULL))
  }

  message("[INFO] Generating cohort-level plots from: ", meta_csv)

  meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
  if (!("sample" %in% colnames(meta))) stop("[ERROR] Metadata needs a 'sample' column.")

  # Derive cohort (CLL / BL / MPN)
  if ("cohort" %in% colnames(meta)) {
    meta$Cohort <- meta$cohort
  } else if ("Cohort" %in% colnames(meta)) {
    meta$Cohort <- meta$Cohort
  } else if ("project" %in% colnames(meta)) {
    meta$Cohort <- dplyr::case_when(
      meta$project == "GSE237907" ~ "CLL",
      meta$project == "GSE199413" ~ "BL",
      meta$project == "GSE228758" ~ "MPN",
      TRUE                        ~ NA_character_
    )
  } else {
    stop("[ERROR] Could not infer cohort; add 'cohort' or 'project' to metadata.")
  }

  meta$Cohort <- factor(meta$Cohort, levels = c("CLL", "BL", "MPN"))
  meta <- meta[!is.na(meta$Cohort), , drop = FALSE]

  # Join scores + meta
  scores_df <- as.data.frame(scores_df)
  scores_df <- tibble::rownames_to_column(scores_df, "sample")

  merged <- dplyr::inner_join(scores_df, meta, by = "sample")

  # Basic sanity
  needed <- c("tumourPurity", "StromalScore", "ImmuneScore")
  missing <- setdiff(needed, colnames(merged))
  if (length(missing) > 0) {
    warning("[WARN] Missing columns in merged data: ", paste(missing, collapse = ", "),
            "\nPurity diagnostics will be partial.")
  }

  # Use project as x-axis grouping if available; otherwise fall back to Cohort
  has_project <- "project" %in% colnames(merged)
  group_var <- if (has_project) "project" else "Cohort"

  # Cohort loop
  cohort_levels <- levels(meta$Cohort)

  # Style helpers
  before_grey <- "#D3D3D3"
  after_blue  <- "#1E88E5"
  pal_groups  <- RColorBrewer::brewer.pal(3, "Set2")[c(2, 1)]

  for (coh in cohort_levels) {
    message("[INFO] --- Cohort: ", coh, " ---")

    meta_c <- meta[meta$Cohort == coh, , drop = FALSE]
    merged_c <- merged[merged$Cohort == coh, , drop = FALSE]

    if (nrow(meta_c) == 0) {
      message("[WARN] No metadata samples for cohort: ", coh, " (skipping).")
      next
    }
    if (!"tumourPurity" %in% colnames(merged_c)) {
      message("[WARN] tumourPurity missing for cohort: ", coh, " (skipping purity diagnostics).")
    }

    # ------------------------------------------------------------
    # A) Retention by project (same-axis) WITHIN THIS COHORT
    # ------------------------------------------------------------
    # BEFORE counts are based on meta (all samples), AFTER counts are purity-kept
    before_tab <- if (has_project && "project" %in% colnames(meta_c)) table(meta_c$project) else table(meta_c$Cohort)
    grp_levels <- names(before_tab)

    # AFTER (only samples in merged_c with valid purity and passing threshold)
    after_tab <- integer(0)
    if ("tumourPurity" %in% colnames(merged_c)) {
      kept <- !is.na(merged_c$tumourPurity) & merged_c$tumourPurity >= purity_threshold
      if (has_project && "project" %in% colnames(merged_c)) {
        after_tab <- table(merged_c$project[kept])
      } else {
        after_tab <- table(merged_c$Cohort[kept])
      }
    }

    # unify levels
    grp_levels <- sort(unique(c(names(before_tab), names(after_tab))))
    before_counts <- as.numeric(before_tab[grp_levels]); before_counts[is.na(before_counts)] <- 0
    after_counts  <- as.numeric(after_tab[grp_levels]);  after_counts[is.na(after_counts)]  <- 0

    df_ret <- data.frame(
      Group = rep(grp_levels, 2),
      Count = c(before_counts, after_counts),
      Stage = rep(c("Before", "After"), each = length(grp_levels))
    )
    df_ret$Group <- factor(df_ret$Group, levels = grp_levels)
    df_ret$Stage <- factor(df_ret$Stage, levels = c("Before", "After"))

    total_before <- sum(before_counts)
    total_after  <- sum(after_counts)

    p_ret <- ggplot(df_ret, aes(x = Group, y = Count, fill = Stage)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = NA) +
      geom_text(
        aes(label = Count, group = Stage),
        position = position_dodge(width = 0.8),
        vjust = -0.4, size = 3.8
      ) +
      scale_fill_manual(values = c("Before" = before_grey, "After" = after_blue), name = NULL) +
      labs(
        title    = sprintf("%s: Sample retention by %s", coh, if (has_project) "project" else "cohort"),
        subtitle = sprintf("Purity ≥ %.1f | n = %d → %d", purity_threshold, total_before, total_after),
        x        = if (has_project) "Project" else "Cohort",
        y        = "Number of samples"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, colour = "grey40"),
        axis.title.x     = element_text(face = "bold"),
        axis.title.y     = element_text(face = "bold"),
        legend.position  = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
      )

    ret_pdf <- file.path(outdir, sprintf("%s_sample_retention_by_project_same_axis.pdf", coh))
    ggsave(ret_pdf, p_ret, width = 7.5, height = 5.5, dpi = 300)
    message("[INFO] Saved: ", ret_pdf)

    # ------------------------------------------------------------
    # B) Purity diagnostics (per cohort): Stromal/Immune correlations + histogram
    # ------------------------------------------------------------
    if (!all(c("tumourPurity", "StromalScore", "ImmuneScore") %in% colnames(merged_c))) {
      message("[WARN] Missing tumourPurity/StromalScore/ImmuneScore for cohort ", coh, "; skipping purity diagnostics.")
      next
    }

    plot_data <- merged_c %>%
      dplyr::filter(!is.na(tumourPurity)) %>%
      dplyr::mutate(
        PurityGroup = cut(
          tumourPurity,
          breaks = c(-Inf, purity_threshold, Inf),
          labels = c(paste0("<", purity_threshold), paste0(">=", purity_threshold)),
          include.lowest = TRUE
        )
      )

    if (nrow(plot_data) == 0) {
      message("[WARN] No non-NA tumourPurity values for cohort ", coh, "; skipping purity diagnostics.")
      next
    }

    p1 <- ggplot(plot_data, aes(x = StromalScore, y = tumourPurity, colour = PurityGroup)) +
      geom_point(size = 2.2, alpha = 0.85, shape = 16) +
      geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
      scale_colour_manual(values = pal_groups, name = "Purity") +
      ggpubr::stat_cor(method = "spearman", label.y.npc = "top", label.x.npc = "left", size = 4, colour = "black") +
      ggpubr::stat_cor(method = "pearson",  label.y.npc = "top", label.x.npc = "left", size = 4, colour = "black", vjust = 2) +
      labs(title = sprintf("%s: Purity vs Stromal score", coh),
           x = "Stromal Score", y = "Tumour Purity") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
      )

    p2 <- ggplot(plot_data, aes(x = ImmuneScore, y = tumourPurity, colour = PurityGroup)) +
      geom_point(size = 2.2, alpha = 0.85, shape = 16) +
      geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
      scale_colour_manual(values = pal_groups, name = "Purity") +
      ggpubr::stat_cor(method = "spearman", label.y.npc = "top", label.x.npc = "left", size = 4, colour = "black") +
      ggpubr::stat_cor(method = "pearson",  label.y.npc = "top", label.x.npc = "left", size = 4, colour = "black", vjust = 2) +
      labs(title = sprintf("%s: Purity vs Immune score", coh),
           x = "Immune Score", y = "Tumour Purity") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
      )

    p3 <- ggplot(plot_data, aes(x = tumourPurity)) +
      geom_histogram(bins = 30, alpha = 0.85, fill = pal_groups[2], colour = NA) +
      geom_vline(xintercept = purity_threshold, linetype = "dashed", colour = "firebrick", linewidth = 1) +
      labs(title = sprintf("%s: Tumour purity distribution", coh),
           x = "Tumour Purity", y = "Number of samples") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
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

    diag_pdf <- file.path(outdir, sprintf("%s_purity_diagnostics.pdf", coh))
    ggsave(diag_pdf, combined, width = 14, height = 10, dpi = 300)
    message("[INFO] Saved: ", diag_pdf)
  }

  invisible(NULL)
}

# MAIN
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_cli <- function(args) {
    res <- list()
    i <- 1
    while (i <= length(args)) {
      key <- args[i]
      if (!startsWith(key, "--")) stop("Expected named arg starting with '--', got: ", key)
      if (i == length(args)) stop("No value provided for ", key)
      val <- args[i + 1]
      res[[sub("^--", "", key)]] <- val
      i <- i + 2
    }
    res
  }
  opts <- parse_cli(args)

  expr_path <- opts[["expr"]]     %||% stop("Missing --expr")
  map_path  <- opts[["map"]]      %||% stop("Missing --map")
  out_dir   <- opts[["outdir"]]   %||% stop("Missing --outdir")
  thr       <- if (!is.null(opts[["threshold"]])) as.numeric(opts[["threshold"]]) else 0.7

  # NEW: optional metadata path for cohort plots
  meta_path <- opts[["meta"]]

  message("[MAIN] Running HEMA tumour purity pipeline with:")
  message("   expr_matrix_path = ", expr_path)
  message("   map_tsv_path     = ", map_path)
  message("   output_dir       = ", out_dir)
  message("   purity_threshold = ", thr)
  if (!is.null(meta_path)) {
    message("   meta_csv         = ", meta_path)
  }

  res <- run_hema_purity(
    expr_matrix_path = expr_path,
    map_tsv_path     = map_path,
    output_dir       = out_dir,
    purity_threshold = thr
  )

  # HEMA-specific cohort plots, if metadata is available
  if (!is.null(meta_path)) {
    make_hema_cohort_plots(
      meta_csv        = meta_path,
      scores_df       = res$scores,
      purity_threshold = thr,
      outdir          = out_dir
    )
  } else {
    message("[INFO] No --meta provided; skipping HEMA cohort plots.")
  }
}