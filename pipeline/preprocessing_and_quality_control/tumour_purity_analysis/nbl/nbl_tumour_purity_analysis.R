#!/usr/bin/env Rscript
# nbl_tumour_purity_analysis.R
# =============================================================================
# Purpose (bioinformatics context):
#   This script estimates tumour microenvironment signals and tumour purity from
#   bulk RNA-seq expression data, then filters samples based on a purity threshold.
#   It is designed for cohort QC in cancer transcriptomics, where low-purity
#   samples can distort downstream analyses (clustering, correlation, DE, etc.).
#
# High-level workflow:
#   1) Load a SummarizedExperiment containing a gene-by-sample count matrix.
#   2) Ensure gene identifiers are HGNC symbols (map from Ensembl if needed).
#   3) Run ESTIMATE scoring (stromal, immune, overall ESTIMATE score).
#   4) Derive tumour purity from the ESTIMATE score when not provided.
#   5) Filter samples by a chosen tumour purity threshold.
#   6) Save filtered matrices and generate cohort-level diagnostic plots.
#
# Exam emphasis:
#   The script enforces identifier consistency, makes sample filtering traceable,
#   and provides quantitative diagnostics showing how purity relates to stromal
#   and immune signatures and how filtering changes cohort composition.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyestimate)           # ESTIMATE-based scoring using HGNC gene symbols
  library(dplyr)                  # data manipulation
  library(ggplot2)                # plotting
  library(ggpubr)                 # correlation annotations and plot arrangement
  library(data.table)             # fast TSV reading for gene mapping
  library(tibble)                 # rownames/column helpers
  library(RColorBrewer)           # categorical palettes (projects, groups)
  library(patchwork)              # plot composition (loaded; not used in final assembly here)
  library(SummarizedExperiment)   # standard container for expression matrices + metadata
})

# --------------------------------------------------------------
# 1. Map Ensembl → HGNC
# --------------------------------------------------------------
# Bioinformatics rationale:
#   ESTIMATE scoring is defined over gene symbols, not Ensembl identifiers.
#   This function standardises gene identifiers so scoring is consistent across
#   cohorts and aligns with reference gene sets used by ESTIMATE.
map_ensembl_to_hgnc <- function(expr_matrix, map_tsv_path) {
  # expr_matrix:
  #   A genes-by-samples matrix. Rownames are expected to be gene identifiers.
  #   If the identifiers are already HGNC-like, the mapping is skipped.
  if (!any(grepl("^ENSG\\d", rownames(expr_matrix)))) {
    message("[INFO] Row names are already HGNC. Skipping mapping.")
    # Uppercasing reduces symbol mismatches due to inconsistent case.
    rownames(expr_matrix) <- toupper(rownames(expr_matrix))
    return(expr_matrix)
  }

  # The mapping table provides a stable crosswalk between Ensembl gene IDs
  # (without version suffix) and HGNC gene symbols.
  message("[INFO] Loading Ensembl → HGNC map...")
  gene_map <- fread(map_tsv_path, sep = "\t", header = TRUE, data.table = FALSE)
  stopifnot(all(c("Ensembl_ID", "HGNC_Symbol") %in% colnames(gene_map)))

  # id_to_symbol supports fast lookup from Ensembl gene ID to HGNC symbol.
  id_to_symbol <- setNames(gene_map$HGNC_Symbol, gene_map$Ensembl_ID)

  # Ensembl IDs may include version suffixes (e.g. ENSG... .12); these are removed
  # to match typical mapping tables and GTF conventions.
  ensembl_ids    <- sub("\\.\\d+$", "", rownames(expr_matrix))
  mapped_symbols <- id_to_symbol[ensembl_ids]

  # Only genes with a valid HGNC symbol are retained for symbol-based workflows.
  valid <- !is.na(mapped_symbols) & nzchar(mapped_symbols)

  message(sprintf(
    "[INFO] Mapped %d / %d genes (%.1f%%)",
    sum(valid), length(ensembl_ids), 100 * mean(valid)
  ))

  expr_mapped <- expr_matrix[valid, , drop = FALSE]
  rownames(expr_mapped) <- toupper(mapped_symbols[valid])

  # Some Ensembl IDs map to the same HGNC symbol (many-to-one mapping).
  # Collapsing duplicates prevents overcounting a single gene symbol.
  df <- as.data.frame(expr_mapped) %>% rownames_to_column("symbol")
  if (any(duplicated(df$symbol))) {
    message("[INFO] Collapsing duplicate HGNC symbols by sum...")
    expr_mapped <- df %>%
      group_by(symbol) %>%
      summarise(across(where(is.numeric), sum), .groups = "drop") %>%
      column_to_rownames("symbol") %>%
      as.matrix()
  }

  expr_mapped
}

# --------------------------------------------------------------
# 2. ESTIMATE scoring + purity
# --------------------------------------------------------------
# Bioinformatics rationale:
#   Different scoring functions and packages may return slightly different column
#   names. This helper standardises score column names and derives tumour purity
#   from ESTIMATEScore when purity is not explicitly provided.
standardize_estimate_scores <- function(scores_df) {
  scores_df <- as.data.frame(scores_df)

  # Some functions may return an explicit "sample" column instead of rownames.
  if ("sample" %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df$sample <- NULL
  }

  # Standardise score column names using robust pattern matching.
  col_lower <- tolower(colnames(scores_df))
  new_names <- colnames(scores_df)
  new_names[grepl("strom", col_lower)]              <- "StromalScore"
  new_names[grepl("\\bimm", col_lower)]             <- "ImmuneScore"
  new_names[grepl("estimate|\\best\\b", col_lower)] <- "ESTIMATEScore"
  colnames(scores_df) <- new_names

  # If tumourPurity is missing, derive it using the commonly used transformation
  # from ESTIMATEScore. Values are clamped to [0, 1] for interpretability.
  if (!"tumourPurity" %in% colnames(scores_df) && "ESTIMATEScore" %in% colnames(scores_df)) {
    scores_df$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)
    scores_df$tumourPurity <- pmin(1, pmax(0, scores_df$tumourPurity))
    message("[INFO] Derived tumourPurity from ESTIMATEScore.")
  }

  scores_df
}

# --------------------------------------------------------------
# 3. Filter by purity
# --------------------------------------------------------------
# Bioinformatics rationale:
#   Low-purity samples contain substantial non-tumour signal (stromal/immune),
#   which can dominate expression variation and bias tumour-centric analyses.
#   Filtering retains samples likely to reflect tumour biology more directly.
filter_samples_by_purity <- function(expr_matrix, purity_scores, threshold = 0.7) {
  if (!"tumourPurity" %in% colnames(purity_scores)) {
    message("[WARN] tumourPurity not found in scores, skipping filtering.")
    return(expr_matrix)
  }

  # Filtering requires consistent sample identifiers between expression columns
  # and purity score rownames.
  common <- intersect(colnames(expr_matrix), rownames(purity_scores))
  if (length(common) == 0) {
    message("[WARN] No overlap between expression samples and purity scores; skipping filtering.")
    return(expr_matrix)
  }

  purity_vec <- purity_scores[common, "tumourPurity", drop = TRUE]
  keep <- !is.na(purity_vec) & purity_vec >= threshold
  kept <- common[keep]

  message(sprintf(
    "[INFO] Kept %d / %d samples (purity ≥ %.2f)",
    length(kept), length(common), threshold
  ))

  expr_matrix[, kept, drop = FALSE]
}

# --------------------------------------------------------------
# 4. Extract project from sample names
# --------------------------------------------------------------
# Bioinformatics rationale:
#   Cancer cohorts often combine samples from multiple sources (GEO series,
#   SRA projects, TARGET, etc.). This function infers cohort provenance from
#   sample naming conventions to enable stratified QC and retention reporting.
extract_project <- function(samples) {
  # The prefix heuristic assumes sample IDs begin with a project token
  # followed by an underscore and additional identifiers.
  prefix <- sub("^([^_]+).*", "\\1", samples)

  project <- dplyr::case_when(
    grepl("^GSE[0-9]+$", prefix) ~ prefix,
    grepl("^SRP[0-9]+$", prefix) ~ prefix,
    grepl("^TARGET", prefix)     ~ "TARGET",
    TRUE                         ~ "Unknown"
  )

  # Failing early here prevents silent mislabelling of samples, which would
  # weaken any project-level QC interpretation.
  if (any(project == "Unknown")) {
    bad <- unique(samples[project == "Unknown"])
    stop(
      "[ERROR] Some samples could not be mapped to a project: \n  ",
      paste(bad, collapse = ", "),
      "\nPlease check their naming / patterns in extract_project()."
    )
  }

  project
}

# --------------------------------------------------------------
# 5. MAIN PIPELINE
# --------------------------------------------------------------
# This function orchestrates the full tumour purity workflow:
#   - loads data from disk
#   - assigns project labels
#   - maps gene IDs to HGNC
#   - scores tumour microenvironment signatures
#   - filters by purity and saves outputs
#   - generates diagnostics to justify the filtering choice
run_tidyestimate_with_barplot <- function(
  se_path          = "/home/chu25/data/nbl/nbl_combined_count_matrix_ensembl.rds",
  output_dir       = "/work/ugbogu/pipeline/results/tumour_purity_analysis/nbl",
  map_tsv          = "/home/chu25/data/nbl/final_nbl_ensembl_to_hgnc.tsv",
  purity_threshold = 0.7
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set up log file
  log_file <- file.path(output_dir, paste0("nbl_purity_pipeline_", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                            ".log"))
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "output", split = TRUE)
  sink(log_conn, type = "message")
  
  message(sprintf("[INFO] Log file created: %s", log_file))
  message(sprintf("[INFO] Pipeline start time: %s", Sys.time()))
  
  # Wrap pipeline in tryCatch to ensure log cleanup
  tryCatch({

    # Load SummarizedExperiment and extract the expression assay.
    # The assay is treated as a genes-by-samples matrix for scoring and filtering.
    message("[INFO] Loading SummarizedExperiment...")
    se <- readRDS(se_path)
    counts <- SummarizedExperiment::assay(se)
    stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)))
    message("[INFO] Loaded: ", nrow(se), " genes × ", ncol(se), " samples")

    # Project assignment before filtering supports an explicit retention summary.
    colData(se)$project <- extract_project(colnames(se))
    project_before <- table(colData(se)$project)
    message("[INFO] Project counts BEFORE filtering:")
    print(project_before)

    # Map Ensembl identifiers to HGNC symbols for ESTIMATE compatibility.
    counts_hgnc <- map_ensembl_to_hgnc(counts, map_tsv)

    # ----------------------------------------------------------
    # Manual common-gene filter
    # ----------------------------------------------------------
    # ESTIMATE scoring relies on a defined set of genes. The script selects the
    # intersection of the cohort's HGNC symbols with the package's common gene set.
    # A safeguard allows scoring to proceed even if the overlap is unexpectedly low.
    message("[INFO] Selecting common ESTIMATE genes (manual filter)...")
    common_df      <- tidyestimate::common_genes
    common_symbols <- unique(common_df$hgnc_symbol)
    keep_genes     <- intersect(rownames(counts_hgnc), common_symbols)

    if (length(keep_genes) < 50) {
      warning(sprintf(
        "[WARN] Only %d common genes found; proceeding with ALL HGNC genes.",
        length(keep_genes)
      ))
      counts_sig <- counts_hgnc
    } else {
      message(sprintf("[INFO] Using %d common genes for ESTIMATE scoring.", length(keep_genes)))
      counts_sig <- counts_hgnc[keep_genes, , drop = FALSE]
    }

    # Convert to the input format expected by tidyestimate: a data frame with an
    # explicit gene symbol column and sample columns containing expression values.
    df_for_estimate <- as.data.frame(counts_sig) %>%
      tibble::rownames_to_column("hgnc_symbol")

    # Run ESTIMATE scoring and standardise outputs for consistent downstream use.
    message("[INFO] Running ESTIMATE...")
    scores_raw <- tidyestimate::estimate_score(df_for_estimate, is_affymetrix = FALSE)
    scores     <- standardize_estimate_scores(scores_raw)
    message(sprintf("[INFO] Scores: %d samples × %d metrics", nrow(scores), ncol(scores)))

    # Save scores for reproducibility and later reporting.
    scores_csv <- file.path(output_dir, "estimate_scores.csv")
    write.csv(scores, scores_csv, row.names = TRUE)
    message("[INFO] Scores saved to: ", scores_csv)

    # Filter by tumour purity and save the filtered expression matrix in HGNC space.
    counts_filtered <- filter_samples_by_purity(counts_hgnc, scores, purity_threshold)
    filtered_rds <- file.path(output_dir, sprintf("filtered_hgnc_counts_purity%.1f.rds", purity_threshold))
    saveRDS(counts_filtered, filtered_rds)
    message("[INFO] Filtered counts (HGNC) saved to: ", filtered_rds)

    # Also save the same filtered sample subset in the original Ensembl space.
    # This is practical when downstream workflows expect Ensembl identifiers.
    kept_samples <- colnames(counts_filtered)
    ensembl_filtered <- counts[, kept_samples, drop = FALSE]
    ensembl_out <- file.path(
      "/home/chu25/data/nbl",
      sprintf("nbl_filtered_count_matrix_ensembl_purity%.1f.rds", purity_threshold)
    )
    saveRDS(ensembl_filtered, ensembl_out)
    message("[INFO] Filtered Ensembl matrix saved to: ", ensembl_out)

    # Project counts after filtering quantify sample retention by study source.
    project_after <- table(extract_project(kept_samples))
    message("[INFO] Project counts AFTER filtering:")
    print(project_after)

    # --------------------------------------------------------------
    # 6. Bar plot: project retention before vs after filtering
    # --------------------------------------------------------------
    # This plot provides an audit trail showing whether purity filtering
    # disproportionately removes samples from specific data sources.
    message("[INFO] Generating bar plot...")
    total_before <- sum(project_before)
    total_after  <- sum(project_after)

    proj_levels <- sort(unique(c(names(project_before), names(project_after))))
    project_before <- project_before[proj_levels]
    project_after  <- project_after[proj_levels]
    project_before[is.na(project_before)] <- 0
    project_after[is.na(project_after)]   <- 0

    df_plot <- data.frame(
      Project = rep(proj_levels, 2),
      Count   = c(as.numeric(project_before), as.numeric(project_after)),
      Stage   = rep(c("Before Filtering", "After Filtering"), each = length(proj_levels))
    ) %>%
      mutate(
        Project = factor(Project, levels = proj_levels),
        Stage   = factor(Stage, levels = c("Before Filtering", "After Filtering"))
      )

    # A consistent project colour palette improves interpretability across plots.
    n_proj <- length(proj_levels)
    palette <- RColorBrewer::brewer.pal(max(3, min(8, n_proj)), "Set2")[seq_len(n_proj)]
    project_colors <- setNames(palette, proj_levels)

    # The "Before" bars are grey to visually distinguish baseline cohort composition.
    fill_vals <- c(
      setNames(rep("#D3D3D3", n_proj), paste("Before Filtering", proj_levels, sep = ".")),
      setNames(project_colors, paste("After Filtering", proj_levels, sep = "."))
    )

    p_bar <- ggplot(df_plot, aes(x = Project, y = Count, fill = interaction(Stage, Project))) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white", linewidth = 0.6) +
      geom_text(
        aes(label = Count, group = Stage),
        position = position_dodge(width = 0.8),
        vjust = -0.6, size = 4.2, fontface = "bold", color = "black"
      ) +
      scale_fill_manual(values = fill_vals, name = NULL) +
      labs(
        title    = "Sample Retention by Project and Purity Filtering",
        subtitle = sprintf("Purity >= %.1f | n = %d -> %d", purity_threshold, total_before, total_after),
        x = "Project",
        y = "Number of Samples"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle      = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.title.x       = element_text(face = "bold"),
        axis.title.y       = element_text(face = "bold"),
        legend.position    = "top",
        legend.text        = element_text(size = 11),
        panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_rect(color = "gray80", fill = NA, size = 0.5)
      )

    bar_pdf <- file.path(output_dir, "sample_retention_by_project_same_axis.pdf")
    ggsave(bar_pdf, p_bar, width = 10, height = 7, dpi = 300, device = cairo_pdf)
    message("[SUCCESS] Bar plot saved: ", bar_pdf)

    # --------------------------------------------------------------
    # 7. Purity diagnostics
    # --------------------------------------------------------------
    # These plots justify the purity threshold by showing relationships between
    # estimated purity and tumour microenvironment scores, and by visualising the
    # purity distribution across projects.
    plot_data <- scores %>%
      as.data.frame() %>%
      rownames_to_column("Sample") %>%
      mutate(
        Project     = extract_project(Sample),
        PurityGroup = cut(
          tumourPurity,
          breaks = c(-Inf, purity_threshold, Inf),
          labels = c(paste("<", purity_threshold), paste(">=", purity_threshold)),
          include.lowest = TRUE
        )
      )

    pal <- RColorBrewer::brewer.pal(3, "Set2")[c(2, 1)]

    # A: Purity vs stromal score (expected negative association in many cohorts)
    p1 <- ggplot(plot_data, aes(x = StromalScore, y = tumourPurity, colour = PurityGroup)) +
      geom_point(size = 2.2, alpha = 0.85, shape = 16) +
      geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
      scale_colour_manual(values = pal, name = "Purity") +
      ggpubr::stat_cor(method = "spearman", label.y = 0.95, size = 4,
                       label.x.npc = "left", colour = "black") +
      ggpubr::stat_cor(method = "pearson", label.y = 0.88, size = 4,
                       label.x.npc = "left", colour = "black") +
      labs(title = "Tumour Purity vs Stromal Score", x = "Stromal Score", y = "Tumour Purity") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
      )

    # B: Purity vs immune score (often negative association, cohort-dependent)
    p2 <- ggplot(plot_data, aes(x = ImmuneScore, y = tumourPurity, colour = PurityGroup)) +
      geom_point(size = 2.2, alpha = 0.85, shape = 16) +
      geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
      scale_colour_manual(values = pal, name = "Purity") +
      ggpubr::stat_cor(method = "spearman", label.y = 0.95, size = 4,
                       label.x.npc = "left", colour = "black") +
      ggpubr::stat_cor(method = "pearson", label.y = 0.88, size = 4,
                       label.x.npc = "left", colour = "black") +
      labs(title = "Tumour Purity vs Immune Score", x = "Immune Score", y = "Tumour Purity") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
      )

    # C: Purity distribution by project (checks for project-specific biases)
    p3 <- ggplot(plot_data, aes(x = tumourPurity, fill = Project)) +
      geom_histogram(bins = 40, alpha = 0.8, position = "stack", colour = NA) +
      geom_vline(xintercept = purity_threshold, linetype = "dashed", colour = "firebrick", linewidth = 1) +
      annotate(
        "text", x = purity_threshold + 0.02, y = Inf, vjust = 1.5, hjust = 0,
        label = paste("Threshold =", purity_threshold),
        colour = "firebrick", size = 4
      ) +
      scale_fill_manual(values = project_colors, name = NULL) +
      labs(title = "Distribution of Tumour Purity by Project",
           x = "Tumour Purity", y = "Number of Samples") +
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

    diag_pdf <- file.path(output_dir, "purity_diagnostics.pdf")
    ggsave(diag_pdf, combined, width = 14, height = 10, dpi = 300, device = cairo_pdf)
    message("[SUCCESS] Diagnostics saved: ", diag_pdf)

    message(sprintf("[INFO] Pipeline end time: %s", Sys.time()))
    message("[SUCCESS] All done! Output in: ", output_dir)
    
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
# 8. Run with cohort-specific paths
# --------------------------------------------------------------
# The run call defines:
#   - the input SummarizedExperiment file (counts in Ensembl space),
#   - the mapping table to HGNC symbols,
#   - the output directory for scores, filtered matrices, and plots,
#   - the purity threshold used for sample retention.
run_tidyestimate_with_barplot(
  se_path    = "/home/chu25/data/nbl/nbl_combined_count_matrix_ensembl.rds",
  output_dir = "/work/ugbogu/pipeline/results/tumour_purity_analysis/nbl",
  map_tsv    = "/home/chu25/data/nbl/final_nbl_ensembl_to_hgnc.tsv",
  purity_threshold = 0.7
)