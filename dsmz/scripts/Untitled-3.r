# =============================================================================
# PAM50 Reference Builder & DSMZ Label Transfer (Refactor aligned to script2)
# Author: You
# Date: 2025-11-01
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(uwot)
})

# --- Helpers assumed from your script2 environment ---------------------------
# - log_info(), log_warn(), log_debug()
# - safe_pdf(file, expr)
# If they don't exist in your session, define minimal fallbacks:
if (!exists("log_info")) log_info <- function(fmt, ...) cat(sprintf(paste0("[INFO] ", fmt, "\n"), ...))
if (!exists("log_warn")) log_warn <- function(fmt, ...) cat(sprintf(paste0("[WARN] ", fmt, "\n"), ...))
if (!exists("log_debug")) log_debug <- function(fmt, ...) cat(sprintf(paste0("[DEBUG] ", fmt, "\n"), ...))
if (!exists("safe_pdf")) {
  safe_pdf <- function(file, expr) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file); on.exit(grDevices::dev.off(), add = TRUE); force(expr)
  }
}

# =============================================================================
# (A) Build Balanced TCGA Reference (Top k per PAM50 subtype)
# =============================================================================
build_tcga_reference <- function(
  M_pam50,                      # PAM50 matrix (genes x samples)
  V_tcga,                       # character vector of TCGA sample IDs (colnames of M_pam50)
  colData_tcga,                 # data.frame or DataFrame with sample annotations (rownames/row ids = sample IDs)
  subtype_labels = c("HER2-high", "HER2-low", "LumA", "LumB", "Basal"),
  k_per_subtype = 20,
  tcga_sub_col_opts = c("paper_BRCA_Subtype_PAM50", "PAM50", "Subtype", "BRCA_Subtype", "molecular_subtype"),
  scale_genes = TRUE,           # z-score genes across TCGA subset before distance/UMAP
  seed = 42,
  outdir = "output"
) {
  set.seed(seed)
  ref_dir <- file.path(outdir, "tcga_reference_selection")
  dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------- 1. Validate inputs -------------------
  stopifnot(is.matrix(M_pam50), is.character(V_tcga))
  if (is.null(colnames(M_pam50))) stop("[FATAL] M_pam50 must have colnames = sample IDs.")
  if (is.null(rownames(M_pam50))) stop("[FATAL] M_pam50 must have rownames = gene IDs.")
  V_tcga <- intersect(V_tcga, colnames(M_pam50))
  if (length(V_tcga) < 5) stop("[FATAL] Too few TCGA samples after intersection with M_pam50.")

  # Coerce colData to data.frame with rownames = sample IDs
  if (is.null(colData_tcga)) stop("[FATAL] colData_tcga cannot be NULL.")
  if (!is.data.frame(colData_tcga)) colData_tcga <- as.data.frame(colData_tcga)
  rn_ok <- rownames(colData_tcga)
  if (is.null(rn_ok)) {
    stop("[FATAL] colData_tcga must have rownames equal to TCGA sample IDs.")
  }

  # ------------------- 2. Subset & (optionally) scale -------------------
  log_info("[A.1] Subsetting PAM50 for %d TCGA samples", length(V_tcga))
  M_pam50_tcga <- M_pam50[, V_tcga, drop = FALSE]
  X_tcga <- t(M_pam50_tcga)  # samples x genes
  if (scale_genes) {
    log_info("[A.1] Z-scoring genes across TCGA subset")
    X_tcga <- scale(X_tcga)                   # center/scale per gene (column)
  }
  X_tcga <- as.matrix(X_tcga)
  log_info("[A.1] Matrix ready: %d samples x %d genes", nrow(X_tcga), ncol(X_tcga))

  # ------------------- 3. Cluster TCGA (k=5) -------------------
  log_info("[A.2] Hierarchical clustering (k=5, euclidean/average)")
  d_tcga <- dist(X_tcga, method = "euclidean")
  hc_tcga <- hclust(d_tcga, method = "average")
  cl_tcga <- cutree(hc_tcga, k = 5)

  # ------------------- 4. Map clusters to PAM50 subtypes -------------------
  log_info("[A.3] Mapping clusters to PAM50 subtypes")
  tcga_sub_col <- intersect(tcga_sub_col_opts, colnames(colData_tcga))[1]
  if (!is.na(tcga_sub_col)) {
    lab_tcga_full <- as.character(colData_tcga[rownames(X_tcga), tcga_sub_col])
    sub_lvls <- sort(unique(na.omit(lab_tcga_full)))
    log_info("[A.3] Using annotation '%s' with %d subtype levels", tcga_sub_col, length(sub_lvls))

    # Mean profiles by annotated subtype and by cluster (use *unscaled* gene space for biological centroids)
    sub_means <- sapply(sub_lvls, function(s) rowMeans(M_pam50_tcga[, lab_tcga_full == s, drop = FALSE], na.rm = TRUE))
    cl_means  <- sapply(1:5,     function(k) rowMeans(M_pam50_tcga[, cl_tcga == k,          drop = FALSE], na.rm = TRUE))

    g_common <- intersect(rownames(sub_means), rownames(cl_means))
    if (length(g_common) < 10) stop("[FATAL] Too few overlapping genes when mapping clusters to subtypes.")
    corr_mat <- cor(cl_means[g_common, ], sub_means[g_common, ], method = "spearman", use = "pairwise.complete.obs")
    cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])

    # Map to desired output label order if possible; otherwise keep discovered names
    cl_to_sub_named <- setNames(
      ifelse(cl_to_sub %in% subtype_labels, cl_to_sub, cl_to_sub),
      nm = as.character(1:5)
    )
  } else {
    log_warn("[A.3] No subtype annotation found → fallback to k-means template ordering")
    km <- kmeans(X_tcga, centers = 5, nstart = 50)
    km_cent <- t(sapply(1:5, function(k) colMeans(X_tcga[km$cluster == k, , drop = FALSE])))
    ord <- order(prcomp(km_cent)$x[, 1])
    cl_to_sub_named <- setNames(subtype_labels[ord], nm = as.character(1:5))
  }

  map_df <- data.frame(Cluster = names(cl_to_sub_named), Mapped_Subtype = cl_to_sub_named, row.names = NULL)
  log_info("[A.3] Cluster → Subtype mapping:\n%s", paste(utils::capture.output(print(map_df, row.names = FALSE)), collapse = "\n"))

  # ------------------- 5. Final subtype per TCGA sample -------------------
  subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)],
                         levels = union(subtype_labels, unique(unname(cl_to_sub_named))))
  log_info("[A.4] Subtype counts:\n%s", paste(utils::capture.output(print(table(subtype_tcga))), collapse = "\n"))

  # ------------------- 6. Select top-k nearest to subtype centroid -------------------
  get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
  get_k_nearest <- function(X, cen, k) {
    d <- sqrt(rowSums((X - matrix(cen, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
    ord <- order(d)[seq_len(min(k, length(d)))]
    data.frame(sample = rownames(X)[ord], dist = d[ord], stringsAsFactors = FALSE)
  }

  log_info("[A.5] Selecting up to %d samples per subtype", k_per_subtype)
  sel_list <- lapply(levels(subtype_tcga), function(sub) {
    idx <- which(subtype_tcga == sub)
    if (length(idx) == 0) return(NULL)
    cen <- get_centroid(X_tcga, idx)
    nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = k_per_subtype)
    nearest$subtype <- sub
    nearest
  })
  sel_df <- bind_rows(sel_list) %>% relocate(subtype, .before = sample) %>%
    group_by(subtype) %>% slice_head(n = k_per_subtype) %>% ungroup()

  log_info("[A.5] Final reference counts:\n%s", paste(utils::capture.output(print(table(sel_df$subtype))), collapse = "\n"))
  stopifnot(length(unique(sel_df$sample)) == nrow(sel_df))

  # ------------------- 7. Save reference (CSV + RDS) -------------------
  ref_samples <- sel_df$sample
  M_pam50_tcga_ref <- M_pam50_tcga[, ref_samples, drop = FALSE]

  out_csv <- file.path(ref_dir, sprintf("tcga_top%d_per_subtype.csv", k_per_subtype))
  write.csv(sel_df, out_csv, row.names = FALSE)
  log_info("[A.6] Saved selection → %s", out_csv)

  out_rds <- file.path(ref_dir, sprintf("PAM50_tcga_top%d_per_subtype.rds", k_per_subtype))
  saveRDS(M_pam50_tcga_ref, out_rds)
  log_info("[A.6] Saved matrix → %s", out_rds)

  # ------------------- 8. UMAP visualization (highlight selected) -------------------
  log_info("[A.7] Generating UMAP for TCGA subset")
  set.seed(seed)
  emb <- uwot::umap(X_tcga, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample  <- rownames(X_tcga)
  emb_df$subtype <- subtype_tcga
  emb_df$chosen  <- ifelse(emb_df$sample %in% ref_samples, "TopK", "Other")

  p_ref <- ggplot(emb_df[emb_df$chosen == "TopK", ], aes(UMAP1, UMAP2, color = subtype)) +
    geom_point(alpha = 0.95, size = 2) + theme_bw() +
    labs(title = sprintf("TCGA-BRCA: Top%d per PAM50 Subtype (PAM50 space)", k_per_subtype))

  out_pdf <- file.path(ref_dir, sprintf("umap_tcga_reference_top%d.pdf", k_per_subtype))
  ggsave(out_pdf, p_ref, width = 7, height = 6)
  log_info("[A.7] Saved UMAP → %s", out_pdf)

  # ------------------- 9. Export plain list (per subtype) -------------------
  dumpfile <- file.path(ref_dir, sprintf("tcga_top%d_samples_by_subtype.txt", k_per_subtype))
  con <- file(dumpfile, "wt"); on.exit(try(close(con), silent = TRUE), add = TRUE)
  by_sub <- split(sel_df$sample, sel_df$subtype)
  for (nm in names(by_sub)) {
    writeLines(sprintf("[%s]", nm), con)
    writeLines(by_sub[[nm]], con)
    writeLines("", con)
  }
  log_info("[A.7] Saved list → %s", dumpfile)

  invisible(list(
    selection   = sel_df,
    matrix      = M_pam50_tcga_ref,   # genes x selected samples
    subtype_map = subtype_tcga,       # factor for all TCGA subset samples
    umap        = emb_df,
    plot        = p_ref,
    hc          = hc_tcga,
    clusters    = cl_tcga
  ))
}

# =============================================================================
# (B) DSMZ Subtype Assignment via TCGA Centroids
# =============================================================================
assign_dsmz_subtypes <- function(
  M_pam50,                 # PAM50 matrix (genes x samples)
  V_dsmz,                  # character vector of DSMZ sample IDs (colnames of M_pam50)
  tcga_ref_matrix,         # genes x selected TCGA reference samples (from build_tcga_reference()$matrix)
  tcga_ref_selection,      # data.frame with columns: subtype, sample (from build_tcga_reference()$selection)
  subtype_labels = c("HER2-high", "HER2-low", "LumA", "LumB", "Basal"),
  scale_genes = TRUE,      # z-score within DSMZ subset before UMAP (Spearman is rank-based; scaling is benign)
  seed = 42,
  outdir = "output"
) {
  set.seed(seed)
  dsmz_dir <- file.path(outdir, "dsmz_from_tcga_reference")
  dir.create(dsmz_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------- 1. Validate reference & DSMZ set -------------------
  stopifnot(is.matrix(M_pam50), is.matrix(tcga_ref_matrix))
  if (is.null(colnames(M_pam50)) || is.null(rownames(M_pam50)))
    stop("[FATAL] M_pam50 must have row/colnames.")
  if (is.null(rownames(tcga_ref_matrix))) stop("[FATAL] tcga_ref_matrix must have rownames (genes).")

  V_dsmz <- intersect(V_dsmz, colnames(M_pam50))
  if (length(V_dsmz) == 0) stop("[FATAL] No DSMZ samples in M_pam50 colnames.")
  log_info("[B.1] Loading TCGA reference: %d genes x %d samples", nrow(tcga_ref_matrix), ncol(tcga_ref_matrix))

  # ------------------- 2. Build centroids (in raw PAM50 space) -------------------
  req_cols <- c("subtype", "sample")
  if (!all(req_cols %in% colnames(tcga_ref_selection))) {
    stop(sprintf("[FATAL] tcga_ref_selection must contain columns: %s", paste(req_cols, collapse = ", ")))
  }
  sel_info <- tcga_ref_selection[tcga_ref_selection$sample %in% colnames(tcga_ref_matrix), ]
  if (!nrow(sel_info)) stop("[FATAL] Reference selection is empty after intersecting with tcga_ref_matrix colnames.")

  centroids <- sapply(subtype_labels, function(s) {
    ss <- sel_info$sample[sel_info$subtype == s]
    if (!length(ss)) return(rep(NA_real_, nrow(tcga_ref_matrix)))
    rowMeans(tcga_ref_matrix[, ss, drop = FALSE], na.rm = TRUE)
  })
  rownames(centroids) <- rownames(tcga_ref_matrix)
  log_info("[B.2] Built %d centroids", ncol(centroids))

  # ------------------- 3. Prepare DSMZ matrix & genes -------------------
  M_pam50_dsmz <- M_pam50[, V_dsmz, drop = FALSE]
  log_info("[B.3] DSMZ matrix: %d genes x %d samples", nrow(M_pam50_dsmz), ncol(M_pam50_dsmz))

  g_common <- intersect(rownames(centroids), rownames(M_pam50_dsmz))
  if (length(g_common) < 10) stop("[FATAL] Too few overlapping genes between centroids and DSMZ.")
  log_info("[B.4] %d common genes", length(g_common))

  # ------------------- 4. Assign via Spearman correlation -------------------
  assign_sample <- function(x, C) {
    cor_vec <- apply(C, 2, function(cen) suppressWarnings(
      cor(x, cen, method = "spearman", use = "pairwise.complete.obs")
    ))
    best <- which.max(cor_vec)
    c(label = colnames(C)[best], score = cor_vec[best])
  }

  log_info("[B.5] Assigning %d DSMZ samples", ncol(M_pam50_dsmz))
  assign_mat <- t(apply(M_pam50_dsmz[g_common, , drop = FALSE], 2, assign_sample,
                        C = centroids[g_common, , drop = FALSE]))
  dsmz_assign <- as.data.frame(assign_mat, stringsAsFactors = FALSE)
  dsmz_assign$sample <- rownames(assign_mat)
  dsmz_assign$label  <- factor(dsmz_assign$label, levels = subtype_labels)
  dsmz_assign$score  <- as.numeric(dsmz_assign$score)

  log_info("[B.5] Predicted counts:\n%s", paste(utils::capture.output(print(table(dsmz_assign$label))), collapse = "\n"))

  # ------------------- 5. Save assignments -------------------
  out_csv <- file.path(dsmz_dir, "dsmz_pred_subtypes_from_tcga_topK.csv")
  write.csv(dsmz_assign, out_csv, row.names = FALSE)
  log_info("[B.6] Saved → %s", out_csv)

  # ------------------- 6. UMAP visualization (optional scaling for display) -------------------
  log_info("[B.7] Generating UMAP for DSMZ subset")
  set.seed(seed)
  X_dsmz <- t(M_pam50_dsmz[g_common, , drop = FALSE])   # samples x genes
  if (scale_genes) {
    X_dsmz <- scale(X_dsmz)
  }
  X_dsmz <- as.matrix(X_dsmz)
  emb <- uwot::umap(X_dsmz, n_neighbors = 10, min_dist = 0.3, metric = "cosine")
  emb_df <- as.data.frame(emb); colnames(emb_df) <- c("UMAP1", "UMAP2")
  emb_df$sample <- rownames(X_dsmz)
  emb_df$pred   <- dsmz_assign$label[match(emb_df$sample, dsmz_assign$sample)]

  p_pred <- ggplot(emb_df, aes(UMAP1, UMAP2, color = pred)) +
    geom_point(alpha = 0.95, size = 2.1) + theme_bw() +
    labs(title = "DSMZ: Predicted PAM50 Subtype (TCGA centroids)")

  out_pdf <- file.path(dsmz_dir, "umap_dsmz_predicted_labels.pdf")
  ggsave(out_pdf, p_pred, width = 6.5, height = 5.5)
  log_info("[B.7] Saved UMAP → %s", out_pdf)

  # Quick 2-column map for easy joins
  out_tsv <- file.path(dsmz_dir, "dsmz_predicted_label_map.tsv")
  write.table(dsmz_assign[, c("sample", "label")], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  log_info("[B.7] Saved map → %s", out_tsv)

  invisible(list(
    assignment = dsmz_assign,
    centroids  = centroids[g_common, , drop = FALSE],
    umap       = emb_df,
    plot       = p_pred
  ))
}

# =============================================================================
# USAGE EXAMPLE (consistent with script2 style)
# =============================================================================
# ref_tcga <- build_tcga_reference(
#   M_pam50         = M_pam50,
#   V_tcga          = V_tcga,                   # character vector of TCGA sample IDs
#   colData_tcga    = colData(tcga_se),         # or a data.frame with rownames = TCGA IDs
#   k_per_subtype   = 20,
#   outdir          = "output"
# )
#
# dsmz_pred <- assign_dsmz_subtypes(
#   M_pam50            = M_pam50,
#   V_dsmz             = V_dsmz,                # character vector of DSMZ sample IDs
#   tcga_ref_matrix    = ref_tcga$matrix,
#   tcga_ref_selection = ref_tcga$selection,
#   outdir             = "output"
# )
#
# head(dsmz_pred$assignment)
# =============================================================================
