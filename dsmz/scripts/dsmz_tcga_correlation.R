#!/usr/bin/env Rscript

# TCGA-DSMZ Correlation Analysis Pipeline
# Refactored for improved readability and maintainability

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

options(stringsAsFactors = FALSE)

# Load required packages with suppressed startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(matrixStats)
  library(limma)
  library(pheatmap)
  library(uwot)
  library(ggrepel)
})

# Define null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(42)

# Configuration parameters
config <- list(
  # Input file paths
  tcga_se_rds = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds",
  dsmz_rds = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",
  purity_tsv = NULL,  # optional external purity file
  
  # Output directory
  outdir = "/home/chu25/dsmz/results/correlation",
  
  # Analysis parameters
  var_gene_set = "tcga_5k",  # "tcga_5k", "tcga_10k", or "all"
  min_shared = 1000          # minimum shared genes threshold
)

# Create output directory
dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load and validate TCGA data
load_tcga_data <- function(file_path) {
  cat("[INFO] Loading TCGA data...\n")
  tcga_se <- readRDS(file_path)
  tcga_counts <- assay(tcga_se)
  
  cat(sprintf("[INFO] TCGA dims: %d genes x %d samples\n", 
              nrow(tcga_counts), ncol(tcga_counts)))
  
  list(se = tcga_se, counts = tcga_counts)
}

#' Load and validate DSMZ data
load_dsmz_data <- function(counts_path, meta_path) {
  cat("[INFO] Loading DSMZ data...\n")
  
  # Load counts and metadata
  dsmz_raw <- readRDS(counts_path)
  dsmz_meta <- read.csv(meta_path)
  
  # Validate required columns
  required_cols <- c("Ensembl_ID", "gene_name")
  stopifnot(all(required_cols %in% colnames(dsmz_raw)))
  stopifnot("sample_name" %in% names(dsmz_meta))
  
  cat(sprintf("[INFO] DSMZ table: %d rows x %d cols\n", 
              nrow(dsmz_raw), ncol(dsmz_raw)))
  
  list(raw = dsmz_raw, meta = dsmz_meta)
}

# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

#' Harmonize gene IDs between datasets
harmonize_gene_ids <- function(tcga_counts) {
  cat("[INFO] Harmonizing gene IDs...\n")
  
  # Strip version from TCGA: ENSG...#.## -> ENSG...#
  tcga_genes_novers <- sub("\\..*$", "", rownames(tcga_counts))
  
  if (any(duplicated(tcga_genes_novers))) {
    cat("[INFO] Collapsing TCGA duplicate genes after removing version...\n")
    tcga_counts <- rowsum(tcga_counts, tcga_genes_novers, reorder = TRUE)
    rownames(tcga_counts) <- sort(unique(tcga_genes_novers))
  } else {
    rownames(tcga_counts) <- tcga_genes_novers
  }
  
  tcga_counts
}

#' Build DSMZ count matrix
build_dsmz_matrix <- function(dsmz_raw) {
  cat("[INFO] Building DSMZ count matrix...\n")
  
  # Identify sample columns
  annotation_cols <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annotation_cols)
  
  # Convert to numeric
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  # Build matrix
  dsmz_counts <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(dsmz_counts) <- dsmz_raw$Ensembl_ID
  
  # Collapse duplicate Ensembl_IDs
  if (any(duplicated(rownames(dsmz_counts)))) {
    dsmz_counts <- rowsum(dsmz_counts, rownames(dsmz_counts), 
                         reorder = TRUE)
  }
  
  # Remove all-NA columns
  keep_col <- vapply(as.data.frame(dsmz_counts), 
                    function(v) any(!is.na(v)), logical(1))
  
  if (!all(keep_col)) {
    n_dropped <- sum(!keep_col)
    message(sprintf("[WARN] Dropping %d DSMZ columns that are all NA", 
                   n_dropped))
    dsmz_counts <- dsmz_counts[, keep_col, drop = FALSE]
  }
  
  dsmz_counts
}

#' Align metadata with count matrix
align_metadata <- function(dsmz_meta, dsmz_counts, outdir) {
  cat("[INFO] Aligning metadata with count matrix...\n")
  
  dsmz_meta <- dsmz_meta %>% 
    mutate(sample_id = sample_name)
  
  matched_ids <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))
  unmatched_md <- setdiff(dsmz_meta$sample_id, colnames(dsmz_counts))
  unmatched_cnt <- setdiff(colnames(dsmz_counts), dsmz_meta$sample_id)
  
  # Save unmatched samples for debugging
  if (length(unmatched_md)) {
    write.csv(tibble(sample_name = unmatched_md),
             file.path(outdir, "unmatched_metadata_samples.csv"), 
             row.names = FALSE)
  }
  
  if (length(unmatched_cnt)) {
    write.csv(tibble(sample_name = unmatched_cnt),
             file.path(outdir, "unmatched_count_columns.csv"), 
             row.names = FALSE)
  }
  
  # Filter to matched samples
  dsmz_meta <- dsmz_meta %>% 
    filter(sample_id %in% matched_ids)
  dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop = FALSE]
  
  stopifnot(identical(colnames(dsmz_counts), dsmz_meta$sample_id))
  
  cat(sprintf("[INFO] Matched: %d | Unmatched(meta): %d | Unmatched(counts): %d\n",
              length(matched_ids), length(unmatched_md), length(unmatched_cnt)))
  
  list(meta = dsmz_meta, counts = dsmz_counts)
}

#' Map DSMZ samples to broad organ categories
map_to_organs <- function(dsmz_meta) {
  cat("[INFO] Mapping samples to organ categories...\n")
  
  dsmz_meta %>%
    mutate(
      organ = case_when(
        str_detect(Origin, regex("breast", TRUE)) ~ "Breast",
        str_detect(Origin, regex("lung", TRUE)) ~ "Lung",
        str_detect(Origin, regex("liver|hepato|bile", TRUE)) ~ 
          "Liver/Biliary",
        str_detect(Origin, regex("adrenal", TRUE)) ~ "Adrenal",
        str_detect(Origin, regex("brain|gli", TRUE)) ~ "Brain/CNS",
        str_detect(Origin, regex(paste0("lymph|mantle|hodgkin|effusion|",
                                       "b[- ]?cell|t[- ]?cell|nk|",
                                       "natural\\s*killer"), TRUE)) ~ 
          "Lymphoid",
        str_detect(Origin, regex("leukemia", TRUE)) ~ "Leukemia",
        str_detect(Origin, regex("myeloid|megakaryocytic|erythroid|monocytic", 
                                TRUE)) ~ "Myeloid",
        str_detect(Origin, regex("bone\\s*marrow", TRUE)) ~ "Bone Marrow",
        str_detect(Origin, regex("abdomen|pleural", TRUE)) ~ "Body Cavity",
        is.na(Origin) | Origin == "NONE" ~ "Unknown",
        TRUE ~ "Other"
      )
    )
}

# =============================================================================
# PURITY ADJUSTMENT FUNCTIONS
# =============================================================================

#' Load tumor purity data
load_purity_data <- function(tcga_se, purity_file = NULL) {
  purity <- NULL
  
  # Try to get purity from SummarizedExperiment
  if ("purity" %in% colnames(colData(tcga_se))) {
    purity <- as.numeric(colData(tcga_se)$purity)
    names(purity) <- colnames(assay(tcga_se))
    cat("[INFO] Using purity from SummarizedExperiment\n")
  } else if (!is.null(purity_file) && file.exists(purity_file)) {
    # Load from external file
    ptab <- read.delim(purity_file, stringsAsFactors = FALSE)
    stopifnot(all(c("sample", "purity") %in% colnames(ptab)))
    purity <- setNames(ptab$purity, ptab$sample)
    cat("[INFO] Using external purity file\n")
  } else {
    cat("[INFO] No purity data available\n")
  }
  
  purity
}

#' Apply tumor purity adjustment to TCGA data
apply_purity_adjustment <- function(tcga_counts, purity) {
  if (is.null(purity)) {
    return(tcga_counts)
  }
  
  cat("[INFO] Applying tumor purity adjustment...\n")
  
  # Filter to samples with purity data
  keep_samp <- colnames(tcga_counts)[!is.na(purity)]
  purity <- purity[keep_samp]
  tcga_tmp <- tcga_counts[, keep_samp, drop = FALSE]
  
  # Identify purity-associated genes
  ct <- apply(tcga_tmp, 1, function(v) {
    suppressWarnings(cor.test(as.numeric(v), purity, method = "spearman"))
  })
  
  pv <- vapply(ct, function(x) x$p.value, numeric(1))
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))
  padj <- p.adjust(pv, method = "BH")
  
  # Remove strongly purity-associated genes
  drop_genes <- names(which(padj < 0.01 & rh < -0.4))
  
  if (length(drop_genes)) {
    cat(sprintf("[INFO] Removing %d purity-associated genes\n", 
                length(drop_genes)))
    tcga_tmp <- tcga_tmp[setdiff(rownames(tcga_tmp), drop_genes), , 
                        drop = FALSE]
  }
  
  # Apply linear adjustment
  infilt <- 1 - purity
  design <- model.matrix(~ infilt)
  fit <- lmFit(tcga_tmp, design)
  beta <- fit$coefficients[, 2, drop = FALSE]
  tcga_adj_tmp <- as.matrix(tcga_tmp) - beta %*% t(infilt)
  
  # Update original matrix
  tcga_adj <- tcga_counts
  common_adj <- intersect(rownames(tcga_adj_tmp), rownames(tcga_adj))
  tcga_adj[common_adj, colnames(tcga_adj_tmp)] <- 
    tcga_adj_tmp[common_adj, ]
  
  tcga_adj
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Extract TCGA labels and compute cohort means
compute_tcga_means <- function(tcga_se, tcga_adj) {
  cat("[INFO] Computing TCGA cohort means...\n")
  
  # Extract labels
  labs <- NULL
  label_cols <- c("project_id", "study", "disease_type")
  
  for (nm in label_cols) {
    if (nm %in% colnames(colData(tcga_se))) {
      labs <- colData(tcga_se)[[nm]]
      break
    }
  }
  
  stopifnot(!is.null(labs))
  labs <- as.character(labs)
  labs <- sub("^TCGA-", "", labs)
  names(labs) <- colnames(tcga_adj)
  
  # Compute means by cohort
  tcga_levels <- sort(unique(labs))
  tcga_means <- sapply(tcga_levels, function(ct) {
    rowMeans(tcga_adj[, labs == ct, drop = FALSE])
  })
  colnames(tcga_means) <- tcga_levels
  
  cat(sprintf("[INFO] TCGA cohorts available: %d\n", length(tcga_levels)))
  
  list(means = tcga_means, labels = labs)
}

#' Select variable genes based on configuration
select_variable_genes <- function(tcga_means, var_gene_set) {
  sel_genes <- rownames(tcga_means)
  
  if (tolower(var_gene_set) == "tcga_5k") {
    iq <- matrixStats::rowIQRs(tcga_means)
    n_genes <- min(5000, length(iq))
    sel_genes <- names(sort(setNames(iq, rownames(tcga_means)), 
                           decreasing = TRUE))[seq_len(n_genes)]
  } else if (tolower(var_gene_set) == "tcga_10k") {
    iq <- matrixStats::rowIQRs(tcga_means)
    n_genes <- min(10000, length(iq))
    sel_genes <- names(sort(setNames(iq, rownames(tcga_means)), 
                           decreasing = TRUE))[seq_len(n_genes)]
  }
  
  cat(sprintf("[INFO] Variable-gene setting: %s (%d genes used)\n", 
              var_gene_set, length(sel_genes)))
  
  sel_genes
}

#' Compute Spearman correlations between DSMZ and TCGA
compute_correlations <- function(dsmz_counts, tcga_means, dsmz_meta, 
                                sel_genes) {
  cat("[INFO] Computing Spearman correlations...\n")
  cat(sprintf("[INFO] Using %d genes across %d TCGA cohorts and %d DSMZ samples\n",
              length(sel_genes), ncol(tcga_means), ncol(dsmz_counts)))
  
  spearman_scores <- outer(
    dsmz_meta$sample_id, colnames(tcga_means),
    Vectorize(function(cell, cohort) {
      suppressWarnings(cor(dsmz_counts[sel_genes, cell], 
                          tcga_means[sel_genes, cohort],
                          method = "spearman", 
                          use = "pairwise.complete.obs"))
    })
  )
  
  dimnames(spearman_scores) <- list(dsmz_meta$sample_id, 
                                   colnames(tcga_means))
  
  spearman_scores
}

# =============================================================================
# ORGAN MAPPING FUNCTIONS
# =============================================================================

#' Define organ to TCGA mapping
get_organ_tcga_mapping <- function() {
  list(
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
}

#' Perform organ-constrained TCGA matching
perform_organ_mapping <- function(dsmz_meta, spearman_scores, tcga_means) {
  cat("[INFO] Performing organ-constrained mapping...\n")
  
  organ2tcga <- get_organ_tcga_mapping()
  
  dsmz_map <- dsmz_meta %>%
    rowwise() %>%
    mutate(
      allowed = list(intersect(organ2tcga[[organ]] %||% character(0), 
                              colnames(tcga_means))),
      tcga_code = {
        cand <- allowed[[1]]
        if (length(cand) == 0 || 
            !(sample_id %in% rownames(spearman_scores))) {
          NA_character_
        } else {
          cand[which.max(spearman_scores[sample_id, cand])]
        }
      }
    ) %>% 
    ungroup()
  
  dsmz_map
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Create diagnostic plots
create_diagnostic_plots <- function(summary_tbl, spearman_scores, dsmz_map, 
                                   outdir) {
  cat("[INFO] Creating diagnostic plots...\n")
  
  # Violin plot of correlations by organ
  p1 <- ggplot(summary_tbl, aes(organ, best_rho_overall)) +
    geom_violin(fill = "grey85") +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    theme_bw(base_size = 11) + 
    labs(x = NULL, y = "Best overall Spearman ρ") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  ggsave(file.path(outdir, "violin_rho_by_organ.pdf"), 
         p1, width = 7, height = 4)
  
  # Heatmap of top DSMZ samples
  topN <- head(arrange(summary_tbl, desc(best_rho_overall))$sample_id, 25)
  hm <- spearman_scores[topN, , drop = FALSE]
  
  pdf(file.path(outdir, "heatmap_top25_dsmz_vs_tcga_cohorts.pdf"), 
      width = 8, height = 6)
  pheatmap(hm, cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  dev.off()
  
  # Organ-TCGA assignment counts
  counts_long <- dsmz_map %>%
    mutate(TCGA = ifelse(is.na(tcga_code), "Unassigned", tcga_code)) %>%
    count(Organ = organ, TCGA, name = "Freq")
  
  heat_mat <- counts_long %>%
    tidyr::pivot_wider(names_from = TCGA, values_from = Freq, 
                      values_fill = 0) %>%
    column_to_rownames("Organ") %>% 
    as.matrix()
  
  pdf(file.path(outdir, "organ_tcga_assignment_counts.pdf"), 
      width = 8, height = 6)
  pheatmap(heat_mat, cluster_rows = FALSE, cluster_cols = FALSE, 
           display_numbers = TRUE)
  dev.off()
}

#' Create UMAP visualization
create_umap_plots <- function(tcga_means, dsmz_counts, dsmz_meta, 
                             sel_genes, outdir) {
  cat("[INFO] Creating UMAP visualization...\n")
  
  # Combine data
  X <- cbind(tcga_means[sel_genes, , drop = FALSE], 
             dsmz_counts[sel_genes, , drop = FALSE])
  Xz <- t(scale(t(X)))
  Xz[is.na(Xz)] <- 0
  
  # Create labels
  lab <- c(paste0("TCGA_", colnames(tcga_means)), 
           paste0("DSMZ_", colnames(dsmz_counts)))
  src <- c(rep("TCGA", ncol(tcga_means)), 
           rep("DSMZ", ncol(dsmz_counts)))
  group <- c(colnames(tcga_means), dsmz_meta$organ)
  
  # Compute UMAP
  um <- umap(Xz, n_neighbors = 15, min_dist = 0.2, metric = "cosine", 
             verbose = FALSE)
  
  # Create data frame
  um_df <- as.data.frame(um)
  colnames(um_df) <- c("UMAP1", "UMAP2")
  um_df$label <- lab
  um_df$source <- src
  um_df$group <- group
  
  # Save UMAP coordinates
  write.csv(um_df, file.path(outdir, "umap_tcga_means_plus_dsmz.csv"), 
           row.names = FALSE)
  
  # Plot 1: Source comparison
  p1 <- ggplot(um_df, aes(UMAP1, UMAP2, color = source, shape = source)) +
    geom_point(size = 1.8, alpha = 0.8) +
    theme_bw(base_size = 11) + 
    labs(color = NULL, shape = NULL) +
    ggtitle(sprintf("UMAP (genes=%d): DSMZ samples vs. TCGA cohort means", 
                   length(sel_genes)))
  
  ggsave(file.path(outdir, "umap_tcga_means_plus_dsmz.pdf"), 
         p1, width = 7, height = 5)
  
  # Plot 2: DSMZ organs
  p2 <- ggplot(um_df, aes(UMAP1, UMAP2)) +
    geom_point(data = subset(um_df, source == "TCGA"),
               color = "grey80", size = 1.2, alpha = 0.8) +
    geom_point(data = subset(um_df, source == "DSMZ"),
               aes(color = group), size = 1.8, alpha = 0.9) +
    theme_bw(base_size = 11) + 
    labs(color = "DSMZ organ") +
    ggtitle("UMAP: DSMZ colored by organ; TCGA cohort means in gray")
  
  ggsave(file.path(outdir, "umap_dsmz_by_organ.pdf"), 
         p2, width = 7, height = 5)
}

# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

#' Save correlation results
save_results <- function(spearman_scores, summary_tbl, dsmz_map, outdir) {
  cat("[INFO] Saving results...\n")
  
  # Wide and long format correlations
  all_wide <- as.data.frame(spearman_scores) %>% 
    rownames_to_column("sample_id")
  all_long <- all_wide %>% 
    pivot_longer(-sample_id, names_to = "tcga_cohort", values_to = "rho")
  
  write.csv(all_wide, file.path(outdir, "spearman_all_wide.csv"), 
           row.names = FALSE)
  write.csv(all_long, file.path(outdir, "spearman_all_long.csv"), 
           row.names = FALSE)
  saveRDS(spearman_scores, file.path(outdir, "spearman_all.rds"))
  
  # Summary tables
  write.csv(summary_tbl, file.path(outdir, "spearman_best_by_sample.csv"), 
           row.names = FALSE)
  write.csv(dsmz_map %>% select(sample_id, Origin, organ, tcga_code),
           file.path(outdir, "dsmz_tcga_mapping.csv"), 
           row.names = FALSE)
}

#' Save debug information
save_debug_info <- function(tcga_counts, dsmz_raw, common_genes, outdir) {
  cat("[INFO] Saving debug information...\n")
  
  # Gene lists for debugging
  writeLines(utils::capture.output(head(rownames(tcga_counts))),
            file.path(outdir, "head_tcga_genes.txt"))
  writeLines(utils::capture.output(head(dsmz_raw$Ensembl_ID)),
            file.path(outdir, "head_dsmz_genes.txt"))
  writeLines(common_genes[1:min(100, length(common_genes))],
            file.path(outdir, "common_genes_head.txt"))
  
  # Session info
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

main <- function() {
  cat("[INFO] Starting TCGA-DSMZ correlation analysis pipeline...\n")
  
  # 1. Load data
  tcga_data <- load_tcga_data(config$tcga_se_rds)
  dsmz_data <- load_dsmz_data(config$dsmz_rds, config$dsmz_meta_csv)
  
  # 2. Process TCGA data
  tcga_counts <- harmonize_gene_ids(tcga_data$counts)
  storage.mode(tcga_counts) <- "double"
  
  # 3. Process DSMZ data
  dsmz_counts <- build_dsmz_matrix(dsmz_data$raw)
  storage.mode(dsmz_counts) <- "double"
  
  aligned_data <- align_metadata(dsmz_data$meta, dsmz_counts, config$outdir)
  dsmz_meta <- map_to_organs(aligned_data$meta)
  dsmz_counts <- aligned_data$counts
  
  # 4. Find common genes
  common_genes <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
  cat(sprintf("[INFO] Shared genes: %d\n", length(common_genes)))
  
  if (length(common_genes) < config$min_shared) {
    cat("[WARN] Few shared genes; check Ensembl ID consistency\n")
  }
  
  tcga_counts <- tcga_counts[common_genes, , drop = FALSE]
  dsmz_counts <- dsmz_counts[common_genes, , drop = FALSE]
  
  # 5. Apply purity adjustment
  purity <- load_purity_data(tcga_data$se, config$purity_tsv)
  tcga_adj <- apply_purity_adjustment(tcga_counts, purity)
  
  # 6. Compute TCGA cohort means
  tcga_results <- compute_tcga_means(tcga_data$se, tcga_adj)
  tcga_means <- tcga_results$means
  
  # Ensure matrices are aligned
  stopifnot(identical(rownames(tcga_means), rownames(dsmz_counts)))
  
  # 7. Select variable genes
  sel_genes <- select_variable_genes(tcga_means, config$var_gene_set)
  
  # 8. Compute correlations
  spearman_scores <- compute_correlations(dsmz_counts, tcga_means, 
                                        dsmz_meta, sel_genes)
  
  # 9. Perform organ mapping
  dsmz_map <- perform_organ_mapping(dsmz_meta, spearman_scores, tcga_means)
  
  # 10. Create summary table
  best_idx <- max.col(spearman_scores, ties.method = "first")
  best_overall <- tibble(
    sample_id = rownames(spearman_scores),
    best_tcga_overall = colnames(spearman_scores)[best_idx],
    best_rho_overall = spearman_scores[cbind(seq_len(nrow(spearman_scores)), 
                                            best_idx)]
  )
  
  rho_allowed <- rep(NA_real_, nrow(dsmz_map))
  ok <- !is.na(dsmz_map$tcga_code)
  if (any(ok)) {
    rho_allowed[ok] <- spearman_scores[cbind(dsmz_map$sample_id[ok], 
                                            dsmz_map$tcga_code[ok])]
  }
  
  summary_tbl <- dsmz_map %>%
    select(sample_id, Origin, organ, tcga_code) %>%
    left_join(best_overall, by = "sample_id") %>%
    mutate(rho_allowed = rho_allowed)
  
  # 11. Save results
  save_results(spearman_scores, summary_tbl, dsmz_map, config$outdir)
  save_debug_info(tcga_counts, dsmz_data$raw, common_genes, config$outdir)
  
  # 12. Create visualizations
  create_diagnostic_plots(summary_tbl, spearman_scores, dsmz_map, 
                         config$outdir)
  create_umap_plots(tcga_means, dsmz_counts, dsmz_meta, sel_genes, 
                   config$outdir)
  
  cat(sprintf("[INFO] Pipeline completed. Results saved to: %s\n", 
              config$outdir))
}

# Run the pipeline
if (!interactive()) {
  main()
}