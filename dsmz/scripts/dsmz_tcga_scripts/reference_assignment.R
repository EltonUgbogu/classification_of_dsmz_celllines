# PAM50-Based Subtype Reference Building & Label Transfer Pipeline

# --------------------------------------------------------------------------- #
# (A) Build Balanced TCGA Reference (Top 20 per PAM50 subtype)
# --------------------------------------------------------------------------- #
build_tcga_reference <- function(
  M_pam50,              # PAM50 matrix (genes x all samples)
  V_tcga,               # TCGA matrix/data.frame (genes x TCGA samples) OR a matrix-like with these columns
  colData_tcga,         # colData (data.frame) with project_id and PAM50 subtype column
  k_per_subtype = 20,
  seed = 42,
  outdir = "output"
) {
  set.seed(seed)
  ref_dir <- file.path(outdir, "tcga_reference_selection")
  dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------- 1. Harmonize samples & subset to TCGA-BRCA -------------------
  stopifnot("project_id" %in% colnames(colData_tcga))
  stopifnot(!is.null(colnames(V_tcga)), !is.null(rownames(colData_tcga)))
  # Align colData and V_tcga by shared sample names
  sample_names <- intersect(colnames(V_tcga), rownames(colData_tcga))
  if (length(sample_names) == 0) stop("No overlap between V_tcga columns and colData_tcga rownames")
  V_tcga <- V_tcga[, sample_names, drop = FALSE]
  colData_tcga <- colData_tcga[sample_names, , drop = FALSE]
  # BRCA-only filter
  tcga_is_brca <- as.character(colData_tcga[, "project_id"]) == "TCGA-BRCA"
  if (!all(tcga_is_brca, na.rm = TRUE)) {
    keep_names <- rownames(colData_tcga)[tcga_is_brca]
    cat(sprintf("[A.1] Filtering to BRCA-only: %d → %d samples\n", ncol(V_tcga), length(keep_names)))
    V_tcga <- V_tcga[, keep_names, drop = FALSE]
    colData_tcga <- colData_tcga[keep_names, , drop = FALSE]
  }
  # Subset to PAM50 rows already implied by M_pam50
  cat(sprintf("[A.1] Using %d PAM50 genes x %d TCGA-BRCA samples\n", nrow(M_pam50), ncol(V_tcga)))
  M_pam50_tcga <- M_pam50[, colnames(V_tcga), drop = FALSE]
  
  PC_tcga <- make_pcs_matrix(M_pam50_tcga, n_hvg = nrow(M_pam50_tcga), max_pc = 30)
  X_tcga <- PC_tcga
  cat(sprintf("[A.1] PC space: %d samples x %d PCs\n", nrow(X_tcga), ncol(X_tcga)))
  
  if (nrow(X_tcga) < 5) stop("[FATAL] Too few TCGA samples.")
  


  # ------------------- 2. Hierarchical clustering (k=4) -------------------
  cat("[A.2] Hierarchical clustering (k=4, euclidean, average)...\n")
  hc_res <- run_fixed_hc(X_tcga, k = 4, dist_method = "euclidean", linkage = "average")
  hc_tcga <- hc_res$hc
  cl_tcga <- hc_res$clusters
  
  # ------------------- 3. Map clusters to PAM50 subtypes -------------------
  cat("[A.3] Mapping clusters to PAM50 subtypes...\n")
  tcga_sub_col <- intersect(c("paper_BRCA_Subtype_PAM50"), colnames(colData_tcga))[1]
  if (is.na(tcga_sub_col)) stop("'paper_BRCA_Subtype_PAM50' not found in colData_tcga")
  # Build subtype and cluster centroids (genes x categories) using calc_centroids for consistency
  lab_tcga_full <- as.character(colData_tcga[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga_full)))
  hvgs_use <- rownames(M_pam50_tcga)  # all PAM50 genes present
  sub_means <- calc_centroids(M_pam50_tcga, hvgs_use, lab_tcga_full)
  cl_means  <- calc_centroids(M_pam50_tcga, hvgs_use, cl_tcga)
  g_common <- intersect(rownames(sub_means), rownames(cl_means))
  corr_mat <- cor(cl_means[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])
  # Preserve cluster keys from cl_means matrix
  cl_to_sub_named <- setNames(unname(cl_to_sub), colnames(cl_means))
  map_df <- data.frame(Cluster = names(cl_to_sub_named), Mapped_Subtype = cl_to_sub_named)
  cat("[A.3] Cluster → Subtype Mapping:\n"); print(map_df, row.names = FALSE)
  
  # ------------------- 4. Final subtype assignment -------------------
  # Map each TCGA sample's numeric cluster to its mapped subtype label
  cl_keys <- names(cl_to_sub_named)
  # Ensure we can index by cluster labels; coerce clusters to character to match keys
  subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)], levels = sort(unique(unname(cl_to_sub_named))))
  cat("[A.4] Subtype counts:\n"); print(table(subtype_tcga))
  
  # ------------------- 5. Select top-k nearest to centroid -------------------
  
  cat(sprintf("[A.5] Selecting up to %d samples per subtype...\n", k_per_subtype))
  sel_list <- lapply(levels(subtype_tcga), function(sub) {
    idx <- which(subtype_tcga == sub)
    if (length(idx) == 0) return(NULL)
    cen <- get_centroid(X_tcga, idx)
    nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = k_per_subtype)
    nearest$subtype <- sub
    nearest
  })
  sel_df <- dplyr::bind_rows(sel_list)
  sel_df <- dplyr::relocate(sel_df, subtype, .before = sample)
  sel_df <- dplyr::group_by(sel_df, subtype)
  sel_df <- dplyr::slice_head(sel_df, n = k_per_subtype)
  sel_df <- dplyr::ungroup(sel_df)
  
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
  emb <- make_umap(X_tcga, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample <- rownames(X_tcga)
  emb_df$subtype <- subtype_tcga
  emb_df$chosen <- ifelse(emb_df$sample %in% ref_samples, "Top20", "Other")
  
  p_ref <- ggplot2::ggplot(emb_df[emb_df$chosen == "Top20", ], ggplot2::aes(UMAP1, UMAP2, color = subtype)) +
    ggplot2::geom_point(alpha = 0.95, size = 2) + ggplot2::theme_bw() +
    ggplot2::labs(title = "TCGA-BRCA: Top20 per PAM50 Subtype (PAM50 Space)")
  
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
  tcga_ref_matrix,       # From build_tcga_reference() (genes x selected samples)
  tcga_ref_selection,    # sel_df with sample and subtype
  seed = 42,
  outdir = "output",
  knn_k = 5              # K for KNN assignment
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
  subtype_levels <- sort(unique(as.character(tcga_ref_selection$subtype)))
  centroids <- sapply(subtype_levels, function(s) {
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
  
  # ------------------- 5. Assign via KNN on reference samples -------------------
  # Build reference matrix and labels
  X_ref <- t(scale(tcga_ref_matrix[g_common, tcga_ref_selection$sample, drop = FALSE]))
  y_ref <- factor(tcga_ref_selection$subtype, levels = subtype_levels)
  # DSMZ matrix
  X_dsmz <- t(scale(M_pam50_dsmz[g_common, , drop = FALSE]))
  # Cosine distance KNN: convert to cosine similarity => distance = 1 - cos
  unit_norm <- function(A) A / sqrt(rowSums(A^2) + 1e-8)
  Xr <- unit_norm(X_ref); Xq <- unit_norm(X_dsmz)
  sim <- Xq %*% t(Xr)                      # samples_dsmz x samples_ref
  assign_knn <- function(sim_row) {
    ord <- order(sim_row, decreasing = TRUE)[seq_len(min(knn_k, length(sim_row)))]
    labs <- y_ref[ord]
    tab <- sort(table(labs), decreasing = TRUE)
    as.character(names(tab)[1])
  }
  preds <- apply(sim, 1, assign_knn)
  scores <- apply(sim, 1, function(v) max(v, na.rm = TRUE))
  dsmz_assign <- data.frame(sample = rownames(X_dsmz), label = factor(preds, levels = subtype_levels), score = scores, stringsAsFactors = FALSE)
  
  cat("[B.5] Predicted counts:\n"); print(table(dsmz_assign$label))
  
  # ------------------- 6. Save assignment -------------------
  out_csv <- file.path(dsmz_dir, "dsmz_pred_subtypes_from_tcga_top20.csv")
  write.csv(dsmz_assign, out_csv, row.names = FALSE)
  cat(sprintf("[B.6] Saved → %s\n", out_csv))
  
  # ------------------- 7. UMAP visualization -------------------
  cat("[B.7] Generating UMAP...\n")
  set.seed(seed)
  # Combine TCGA reference (Top20 per subtype) and DSMZ for a joint UMAP
  X_ref_umap <- X_ref
  X_joint <- rbind(X_ref_umap, X_dsmz)
  emb <- make_umap(X_joint, n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample <- rownames(X_joint)
  emb_df$set <- ifelse(emb_df$sample %in% rownames(X_ref_umap), "TCGA_BRCA_Ref", "DSMZ")
  # annotate labels
  lab_map <- setNames(as.character(tcga_ref_selection$subtype), tcga_ref_selection$sample)
  emb_df$label <- NA_character_
  emb_df$label[emb_df$set == "TCGA_BRCA_Ref"] <- lab_map[emb_df$sample[emb_df$set == "TCGA_BRCA_Ref"]]
  emb_df$label[emb_df$set == "DSMZ"] <- as.character(dsmz_assign$label[match(emb_df$sample[emb_df$set == "DSMZ"], dsmz_assign$sample)])
  
  p_pred <- ggplot2::ggplot(emb_df, ggplot2::aes(UMAP1, UMAP2, color = label, shape = set)) +
    ggplot2::geom_point(alpha = 0.9, size = 1.9) + ggplot2::theme_bw() +
    ggplot2::labs(title = "TCGA-BRCA (Top20/subtype) + DSMZ (PAM50)", subtitle = sprintf("KNN (k=%d) assignment for DSMZ in PAM50 space", knn_k))
  
  out_pdf <- file.path(dsmz_dir, "umap_tcga_brca_ref_plus_dsmz.pdf")
  ggplot2::ggsave(out_pdf, p_pred, width = 7.5, height = 6)
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