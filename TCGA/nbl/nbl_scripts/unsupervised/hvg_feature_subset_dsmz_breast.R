suppressPackageStartupMessages({
library(matrixStats)
library(stats)
library(dplyr)
library(WGCNA)
library(matrixStats)
library(dplyr)
library(readr)
library(ggplot2)
})

#' @title Unsupervised Feature Selection (5 Methods + WGCNA)
#'
#' @description
#' Performs robust, multi-method feature selection on a **VST-transformed count matrix**:
#' - Variance, MAD, MeanAbsDev, Entropy, PCA loadings → top 3000 each
#' - Union → candidates
#' - Rank by Spearman connectivity
#' - Rank by Spearman × Variance (MX)
#' - Compare to WGCNA kTotal
#' - Save all results + scatter plot
#'
#' @param vst_rds      Path to `.rds` file containing **gene × sample** VST matrix.
#' @param outdir       Output directory (created if missing).
#' @param top_n_method Number ofsbatch genes per method (default 3000).
#' @param final_top    Final number of selected genes (default 500).
#' @param seed         Random seed for reproducibility.
#' @param quiet        Suppress messages? Default `FALSE`.
#'
#' @return Invisibly a list with all results.
#' @export
run_unsupervised_feature_selection <- function(
  vst_rds,
  outdir = "dsmz_breast",
  top_n_method = 3000,
  final_top = 500,
  seed = 123,
  quiet = FALSE
) {
  ## ------------------------------------------------------------------ ##
  ## 0. SETUP
  ## ------------------------------------------------------------------ ##
  set.seed(seed)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  subdir <- file.path(outdir, "feature_selection_unsupervised")
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)

  if (!quiet) message("[INFO] Loading VST matrix from: ", vst_rds)
  counts <- readRDS(vst_rds)
  if (!is.matrix(counts)) stop("VST object must be a matrix")
  if (is.null(rownames(counts))) stop("Matrix must have gene rownames")

  n_genes <- nrow(counts); n_samples <- ncol(counts)
  if (!quiet) message("[INFO] Matrix: ", n_genes, " genes × ", n_samples, " samples")

  ## ------------------------------------------------------------------ ##
  ## 1. HELPER: write vector to file
  ## ------------------------------------------------------------------ ##
  write_vec <- function(vec, filename) {
    writeLines(vec, file.path(subdir, filename))
  }

  ## ------------------------------------------------------------------ ##
  ## 2. FEATURE SELECTION METHODS
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Running 5 feature selection methods...")

  # Variance
  top_var <- names(sort(rowVars(counts, na.rm = TRUE), decreasing = TRUE))[1:top_n_method]

  # MAD
  top_mad <- names(sort(apply(counts, 1, mad, na.rm = TRUE), decreasing = TRUE))[1:top_n_method]

  # Mean Absolute Deviation from Mean
  mean_vals <- rowMeans(counts, na.rm = TRUE)
  absdev <- rowMeans(abs(sweep(counts, 1, mean_vals, "-")), na.rm = TRUE)
  top_mean_absdev <- names(sort(absdev, decreasing = TRUE))[1:top_n_method]

  # Entropy
  top_entropy <- names(sort(apply(counts, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2 || var(x) == 0) return(0)
    p <- table(cut(x, breaks = 20)) / length(x)
    -sum(p * log2(p + 1e-8))
  }), decreasing = TRUE))[1:top_n_method]

  # PCA Loadings (top 5 PCs)
 
  tz_row_vars <- matrixStats::rowSds(counts, na.rm = TRUE)
  tz_row_vars[is.na(tz_row_vars)] <- 0
  tz_keep_genes <- tz_row_vars > 1e-8
  tz_n_removed <- sum(!tz_keep_genes)

  if (tz_n_removed > 0) {
    if (!quiet) message("[INFO] Removing ", tz_n_removed, " constant/zero-variance genes")
    counts <- counts[tz_keep_genes, , drop = FALSE]
  }

  if (nrow(counts) == 0) stop("No variable genes left after filtering")
   X_scaled <- t(scale(t(counts), center = TRUE, scale = FALSE))

  pca <- prcomp(t(X_scaled), center = FALSE, scale. = FALSE)
  loadings <- abs(pca$rotation[, 1:min(5, ncol(pca$rotation)), drop = FALSE])
  gene_scores <- rowSums(loadings)
  top_pca <- names(sort(gene_scores, decreasing = TRUE))[1:top_n_method]

  ## ------------------------------------------------------------------ ##
  ## 3. UNION OF CANDIDATES
  ## ------------------------------------------------------------------ ##
  all_candidates <- unique(c(top_var, top_mad, top_mean_absdev, top_entropy, top_pca))
  if (!quiet) message("[INFO] Union of candidates: ", length(all_candidates), " genes")

  ## ------------------------------------------------------------------ ##
  ## 4. SPEARMAN CONNECTIVITY RANKING
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Computing Spearman correlation on candidates...")
  X_cand <- counts[all_candidates, , drop = FALSE]
  cor_mat <- cor(t(X_cand), method = "spearman", use = "pairwise.complete.obs")
  avg_corr <- rowMeans(abs(cor_mat), na.rm = TRUE)
  final_top500_spearman <- names(sort(avg_corr, decreasing = TRUE))[1:final_top]

  ## ------------------------------------------------------------------ ##
  ## 5. SPEARMAN × VARIANCE (MX)
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Computing Spearman × Variance (MX)...")
  gene_var <- rowVars(X_cand, na.rm = TRUE)
  gene_var_scaled <- gene_var / max(gene_var, na.rm = TRUE)
  mx_score <- avg_corr * gene_var_scaled
  final_top500_mx <- names(sort(mx_score, decreasing = TRUE))[1:final_top]

  ## ------------------------------------------------------------------ ##
  ## 6. WGCNA kTotal
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Running WGCNA soft threshold and kTotal...")
  if (!requireNamespace("WGCNA", quietly = TRUE)) stop("Install WGCNA")
  WGCNA::allowWGCNAThreads()

  datExpr <- t(scale(t(counts), center = TRUE, scale = TRUE))
  powers <- c(1:10, seq(12, 30, 2))
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "unsigned")
  sft_table <- sft$fitIndices
  # Pick power with highest SFT.R.sq and slope in [-2.5, -1.5]
  candidates <- sft_table$SFT.R.sq >= 0.8 & sft_table$slope <= -1.5 & sft_table$slope >= -2.5
  if (sum(candidates) == 0) {
    softPower <- 5  # safe default
  } else {
    softPower <- sft_table$Power[candidates][which.max(sft_table$SFT.R.sq[candidates])]
  }
  message("[INFO] Using soft power: ", softPower)

  adjacency <- WGCNA::adjacency(datExpr, power = softPower, type = "unsigned")
  kTotal <- rowSums(adjacency, na.rm = TRUE)

  ## ------------------------------------------------------------------ ##
  ## 7. JOINT RANKING: kTotal vs MX
  ## ------------------------------------------------------------------ ##
  genes <- intersect(names(kTotal), names(mx_score))
  scores_df <- dplyr::tibble(
    Gene = genes,
    WGCNA_kTotal = kTotal[genes],
    Spearman_meanAbs = avg_corr[genes],
    Variance = gene_var[genes],
    Variance_scaled = gene_var_scaled[genes],
    MX_Score = mx_score[genes]
  ) %>%
    dplyr::mutate(
      Rank_kTotal = rank(-WGCNA_kTotal, ties.method = "min"),
      Rank_MX = rank(-MX_Score, ties.method = "min"),
      Pctl_kTotal = 1 - (Rank_kTotal - 1) / (n() - 1),
      Pctl_MX = 1 - (Rank_MX - 1) / (n() - 1),
      Rank_Sum = Rank_kTotal + Rank_MX,
      Rank_Prod = Rank_kTotal * Rank_MX,
      Rank_Combo = rank(Rank_Sum + 0.001 * sqrt(Rank_Prod))
    ) %>%
    dplyr::arrange(Rank_Combo)

  top500_kTotal <- dplyr::arrange(scores_df, Rank_kTotal) %>% dplyr::slice_head(n = final_top) %>% dplyr::pull(Gene)
  top500_MX <- dplyr::arrange(scores_df, Rank_MX) %>% dplyr::slice_head(n = final_top) %>% dplyr::pull(Gene)
  overlap_500 <- intersect(top500_kTotal, top500_MX)

  ## ------------------------------------------------------------------ ##
  ## 8. SAVE ALL RESULTS
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Saving results to: ", subdir)

  write_vec(top_var, "top3000_variance.txt")
  write_vec(top_mad, "top3000_MAD.txt")
  write_vec(top_mean_absdev, "top3000_meanAbsDev.txt")
  write_vec(top_entropy, "top3000_entropy.txt")
  write_vec(top_pca, "top3000_PCAloadings.txt")
  write_vec(final_top500_spearman, "FINAL_top500_spearman.txt")
  write_vec(final_top500_mx, "FINAL_top500_spearman_variance.txt")
  write_vec(top500_kTotal, "top500_WGCNA_kTotal.txt")
  write_vec(top500_MX, "top500_MX.txt")
  write_vec(overlap_500, "overlap_top500_kTotal_and_MX.txt")

  readr::write_tsv(scores_df, file.path(subdir, "joint_ranks_kTotal_vs_MX.tsv"))

  # Overlap summary
  overlap_summary <- dplyr::tibble(
    Metric = c("Top500_kTotal", "Top500_MX", "Overlap"),
    Count = c(length(top500_kTotal), length(top500_MX), length(overlap_500))
  )
  overlap_summary$Jaccard_Index <- length(overlap_500) / length(union(top500_kTotal, top500_MX))
  readr::write_tsv(overlap_summary, file.path(subdir, "overlap_summary.tsv"))

  # Intersection top 200
  top200_spearman <- names(sort(avg_corr, decreasing = TRUE))[1:200]
  top200_mx <- names(sort(mx_score, decreasing = TRUE))[1:200]
  intersection_200 <- intersect(top200_spearman, top200_mx)
  write_vec(top200_spearman, "top200_spearman.txt")
  write_vec(top200_mx, "top200_spearman_variance.txt")
  write_vec(intersection_200, "intersection_top200_spearman_and_spearman_variance.txt")

  ## ------------------------------------------------------------------ ##
  ## 9. SCATTER PLOT: kTotal vs MX
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Saving scatter plot...")
  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = WGCNA_kTotal, y = MX_Score)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggplot2::labs(
      x = "WGCNA kTotal",
      y = "MX Score (Spearman × Variance)",
      title = "kTotal vs MX across genes"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave(file.path(subdir, "kTotal_vs_MX_scatter.png"), p, width = 7, height = 5, dpi = 300)

  ## ------------------------------------------------------------------ ##
  ## 10. FINAL SUMMARY
  ## ------------------------------------------------------------------ ##
  cat("\n=== UNSUPERVISED FEATURE SELECTION SUMMARY ===\n")
  cat(sprintf("Input matrix: %d genes × %d samples\n", n_genes, n_samples))
  cat(sprintf("Top %d per method: %d each\n", top_n_method, top_n_method))
  cat(sprintf("Union candidates: %d\n", length(all_candidates)))
  cat(sprintf("Final top %d (Spearman): %d\n", final_top, length(final_top500_spearman)))
  cat(sprintf("Final top %d (MX): %d\n", final_top, length(final_top500_mx)))
  cat(sprintf("Top %d kTotal: %d\n", final_top, length(top500_kTotal)))
  cat(sprintf("Overlap kTotal & MX: %d (Jaccard = %.3f)\n",
              length(overlap_500), overlap_summary$Jaccard_Index[3]))
  cat(sprintf("Top 200 intersection: %d\n", length(intersection_200)))
  cat(sprintf("All results saved in: %s\n", subdir))

  invisible(list(
    candidates = all_candidates,
    spearman = final_top500_spearman,
    mx = final_top500_mx,
    kTotal = top500_kTotal,
    overlap = overlap_500,
    intersection_200 = intersection_200,
    scores_df = scores_df
  ))
}



run_unsupervised_feature_selection(
  vst_rds   = "/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/results/dsmz_vst_nbl_batch_corrected.rds",
  outdir    = "/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/results/feature_selection_unsupervised",
  top_n_method = 3000,
  final_top = 500
)