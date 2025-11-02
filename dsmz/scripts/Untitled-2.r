# =============================================================================
# PAM50-Based Subtype Reference Building & Label Transfer Pipeline
# Author: [Your Name]
# Date: 2025-10-31
# =============================================================================

library(dplyr)
library(ggplot2)
library(uwot)
library(tidyr)

# --------------------------------------------------------------------------- #
# (A) Build Balanced TCGA Reference (Top 20 per PAM50 subtype)
# --------------------------------------------------------------------------- #
build_tcga_reference <- function(
  M_pam50,              # Full PAM50 matrix (genes x all samples)
  V_tcga,               # Vector of TCGA sample names
  colData_tcga,         # colData from SummarizedExperiment (or data.frame)
  subtype_labels = c("HER2-high", "HER2-low", "LumA", "LumB", "Basal"),
  k_per_subtype = 20,
  seed = 42,
  outdir = "output"
) {
  set.seed(seed)
  ref_dir <- file.path(outdir, "tcga_reference_selection")
  dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------- 1. Subset & scale TCGA PAM50 -------------------
  cat(sprintf("[A.1] Subsetting PAM50 for %d TCGA samples\n", length(V_tcga)))
  M_pam50_tcga <- M_pam50[, colnames(V_tcga), drop = FALSE]
  X_tcga <- t(scale(M_pam50_tcga))  # samples x genes
  cat(sprintf("[A.1] Scaled: %d samples x %d genes\n", nrow(X_tcga), ncol(X_tcga)))
  
  if (nrow(X_tcga) < 5) stop("[FATAL] Too few TCGA samples.")
  
  # ------------------- 2. Hierarchical clustering (k=5) -------------------
  cat("[A.2] Hierarchical clustering (k=5, euclidean, average)...\n")
  d_tcga <- dist(X_tcga, method = "euclidean")
  hc_tcga <- hclust(d_tcga, method = "average")
  cl_tcga <- cutree(hc_tcga, k = 5)
  
  # ------------------- 3. Map clusters to PAM50 subtypes -------------------
  cat("[A.3] Mapping clusters to PAM50 subtypes...\n")
  tcga_sub_col <- intersect(c("PAM50", "Subtype", "BRCA_Subtype", "molecular_subtype"), colnames(colData_tcga))[1]
  
  if (!is.na(tcga_sub_col)) {
    lab_tcga_full <- as.character(colData_tcga[colnames(V_tcga), tcga_sub_col])
    sub_lvls <- sort(unique(na.omit(lab_tcga_full)))
    cat(sprintf("[A.3] Found annotation '%s' with %d subtypes\n", tcga_sub_col, length(sub_lvls)))
    
    sub_means <- sapply(sub_lvls, function(s) rowMeans(M_pam50_tcga[, lab_tcga_full == s, drop = FALSE], na.rm = TRUE))
    
    cl_means <- sapply(1:5, function(k) rowMeans(M_pam50_tcga[, cl_tcga == k, drop = FALSE], na.rm = TRUE))
    
    g_common <- intersect(rownames(sub_means), rownames(cl_means))
    corr_mat <- cor(cl_means[g_common, ], sub_means[g_common, ], method = "spearman")
    cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])
    cl_to_sub_named <- setNames(subtype_labels[match(cl_to_sub, sub_lvls)], 1:5)
  } else {
    cat("[A.3] [WARN] No annotation → using k-means fallback\n")
    km <- kmeans(X_tcga, centers = 5, nstart = 50)
    km_cent <- t(sapply(1:5, function(k) colMeans(X_tcga[km$cluster == k, , drop = FALSE])))
    ord <- order(prcomp(km_cent)$x[,1])
    cl_to_sub_named <- setNames(subtype_labels[ord], 1:5)
  }
  
  map_df <- data.frame(Cluster = names(cl_to_sub_named), Mapped_Subtype = cl_to_sub_named)
  cat("[A.3] Cluster → Subtype Mapping:\n"); print(map_df, row.names = FALSE)
  
  # ------------------- 4. Final subtype assignment -------------------
  subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)], levels = subtype_labels)
  cat("[A.4] Subtype counts:\n"); print(table(subtype_tcga))
  
  # ------------------- 5. Select top-k nearest to centroid -------------------
  get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
  get_k_nearest <- function(X, cen, k) {
    d <- sqrt(rowSums((X - matrix(cen, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
    ord <- order(d)[seq_len(min(k, length(d)))]
    data.frame(sample = rownames(X)[ord], dist = d[ord], stringsAsFactors = FALSE)
  }
  
  cat(sprintf("[A.5] Selecting up to %d samples per subtype...\n", k_per_subtype))
  sel_list <- lapply(subtype_labels, function(sub) {
    idx <- which(subtype_tcga == sub)
    if (length(idx) == 0) return(NULL)
    cen <- get_centroid(X_tcga, idx)
    nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = k_per_subtype)
    nearest$subtype <- sub
    nearest
  })
  sel_df <- bind_rows(sel_list) %>% relocate(subtype, .before = sample)
  sel_df <- sel_df %>% group_by(subtype) %>% slice_head(n = k_per_subtype) %>% ungroup()
  
  cat("[A.5] Final reference counts:\n"); print(table(sel_df$subtype))
  stopifnot(all(table(sel_df$subtype) <= k_per_subtype))
  stopifnot(length(unique(sel_df$sample)) == nrow(sel_df))
  
  # ------------------- 6. Save reference -------------------
  ref_samples <- sel_df$sample
  M_pam50_tcga_ref <- M_pam50_tcga[, ref_samples, drop = FALSE]
  
  out_csv <- file.path(ref_dir, "tcga_top20_per_subtype.csv")
  write.csv(sel_df, out_csv, row.names = FALSE)
  cat(sprintf("[A.6] Saved selection → %s\n", out_csv))
  
  out_rds <- file.path(ref_dir, "PAM50_tcga_top20_per_subtype.rds")
  saveRDS(M_pam50_tcga_ref, out_rds)
  cat(sprintf("[A.6] Saved matrix → %s\n", out_rds))
  
  # ------------------- 7. UMAP visualization -------------------
  cat("[A.7] Generating UMAP...\n")
  set.seed(seed)
  emb <- <- uwot::umap(X_tcga, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample <- rownames(X_tcga)
  emb_df$subtype <- subtype_tcga
  emb_df$chosen <- ifelse(emb_df$sample %in% ref_samples, "Top20", "Other")
  
  p_ref <- ggplot(emb_df[emb_df$chosen == "Top20", ], aes(UMAP1, UMAP2, color = subtype)) +
    geom_point(alpha = 0.95, size = 2) + theme_bw() +
    labs(title = "TCGA-BRCA: Top20 per PAM50 Subtype (PAM50 Space)")
  
  out_pdf <- file.path(ref_dir, "umap_tcga_reference_top20.pdf")
  ggsave(out_pdf, p_ref, width = 7, height = 6)
  cat(sprintf("[A.7] Saved UMAP → %s\n", out_pdf))
  
  # ------------------- 8. Export plain list -------------------
  dumpfile <- file.path(ref_dir, "tcga_top20_samples_by_subtype.txt")
  con <- file(dumpfile, "wt"); on.exit(close(con), add = TRUE)
  for (nm in names(split(sel_df$sample, sel_df$subtype))) {
    writeLines(sprintf("[%s]", nm), con)
    writeLines(split(sel_df$sample, sel_df$subtype)[[nm]], con)
    writeLines("", con)
  }
  cat(sprintf("[A.7] Saved list → %s\n", dumpfile))
  
  invisible(list(
    selection = sel_df,
    matrix = M_pam50_tcga_ref,
    subtype_map = subtype_tcga,
    umap = emb_df,
    plot = p_ref
  ))
}


# --------------------------------------------------------------------------- #
# (B) DSMZ Subtype Assignment via TCGA Centroids
# --------------------------------------------------------------------------- #
assign_dsmz_subtypes <- function(
  M_pam50,               # Full PAM50 matrix
  V_dsmz,                # DSMZ sample names
  tcga_ref_matrix,       # From build_tcga_reference()
  tcga_ref_selection,    # sel_df from build_tcga_reference()
  subtype_labels = c("HER2-high", "HER2-low", "LumA", "LumB", "Basal"),
  seed = 42,
  outdir = "output"
) {
  dsmz_dir <- file.path(outdir, "dsmz_from_tcga_reference")
  dir.create(dsmz_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------- 1. Validate reference -------------------
  cat(sprintf("[B.1] Loading TCGA reference: %d genes x %d samples\n",
              nrow(tcga_ref_matrix), ncol(tcga_ref_matrix)))
  M_ref <- tcga_ref_matrix
  sel_info <- tcga_ref_selection[tcga_ref_selection$sample %in% colnames(M_ref), ]
  stopifnot(nrow(sel_info) > 0)
  
  # ------------------- 2. Build centroids -------------------
  centroids <- sapply(subtype_labels, function(s) {
    ss <- sel_info$sample[sel_info$subtype == s]
    if (length(ss) == 0) return(rep(NA_real_, nrow(M_ref)))
    rowMeans(M_ref[, ss, drop = FALSE], na.rm = TRUE)
  })
  rownames(centroids) <- rownames(M_ref)
  cat(sprintf("[B.2] Built %d centroids\n", ncol(centroids)))
  
  # ------------------- 3. Prepare DSMZ -------------------
  M_pam50_dsmz <- M_pam50[, colnames(V_dsmz), drop = FALSE]
  cat(sprintf("[B.3] DSMZ matrix: %d genes x %d samples\n", nrow(M_pam50_dsmz), ncol(M_pam50_dsmz)))
  
  # ------------------- 4. Common genes -------------------
  g_common <- intersect(rownames(centroids), rownames(M_pam50_dsmz))
  cat(sprintf("[B.4] %d common genes\n", length(g_common)))
  if (length(g_common) < 10) stop("[FATAL] Too few overlapping genes.")
  
  # ------------------- 5. Assign via Spearman -------------------
  assign_sample <- function(x, C) {
    cor_vec <- apply(C, 2, function(cen) {
      suppressWarnings(cor(x, cen, method = "spearman", use = "pairwise.complete.obs"))
    })
    best <- which.max(cor_vec)
    c(label = names(best), score = cor_vec[best])
  }
  
  cat(sprintf("[B.5] Assigning %d DSMZ samples...\n", ncol(M_pam50_dsmz)))
  assign_mat <- t(apply(M_pam50_dsmz[g_common, , drop = FALSE], 2, assign_sample, C = centroids[g_common, , drop = FALSE]))
  dsmz_assign <- as.data.frame(assign_mat, stringsAsFactors = FALSE)
  dsmz_assign$sample <- rownames(assign_mat)
  dsmz_assign$label <- factor(dsmz_assign$label, levels = subtype_labels)
  dsmz_assign$score <- as.numeric(dsmz_assign$score)
  
  cat("[B.5] Predicted counts:\n"); print(table(dsmz_assign$label))
  
  # ------------------- 6. Save assignment -------------------
  out_csv <- file.path(dsmz_dir, "dsmz_pred_subtypes_from_tcga_top20.csv")
  write.csv(dsmz_assign, out_csv, row.names = FALSE)
  cat(sprintf("[B.6] Saved → %s\n", out_csv))
  
  # ------------------- 7. UMAP visualization -------------------
  cat("[B.7] Generating UMAP...\n")
  set.seed(seed)
  X_dsmz <- t(scale(M_pam50_dsmz[g_common, , drop = FALSE]))
  emb <- uwot::umap(X_dsmz, n_neighbors = 10, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample <- rownames(X_dsmz)
  emb_df$pred <- dsmz_assign$label[match(emb_df$sample, dsmz_assign$sample)]
  
  p_pred <- ggplot(emb_df, aes(UMAP1, UMAP2, color = pred)) +
    geom_point(alpha = 0.95, size = 2.1) + theme_bw() +
    labs(title = "DSMZ: Predicted PAM50 Subtype (TCGA Centroid)")
  
  out_pdf <- file.path(dsmz_dir, "umap_dsmz_predicted_labels.pdf")
  ggsave(out_pdf, p_pred, width = 6.5, height = 5.5)
  cat(sprintf("[B.7] Saved UMAP → %s\n", out_pdf))
  
  out_tsv <- file.path(dsmz_dir, "dsmz_predicted_label_map.tsv")
  write.table(dsmz_assign[, c("sample", "label")], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[B.7] Saved map → %s\n", out_tsv))
  
  invisible(list(
    assignment = dsmz_assign,
    centroids = centroids[g_common, ],
    umap = emb_df,
    plot = p_pred
  ))
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# result_A <- build_tcga_reference(
#   M_pam50 = M_pam50,
#   V_tcga = V_tcga,
#   colData_tcga = colData(tcga_se),
#   outdir = "output"
# )
#
# result_B <- assign_dsmz_subtypes(
#   M_pam50 = M_pam50,
#   V_dsmz = V_dsmz,
#   tcga_ref_matrix = result_A$matrix,
#   tcga_ref_selection = result_A$selection,
#   outdir = "output"
# )
#
# head(result_B$assignment)
# =============================================================================