# TCGA reference building and DSMZ assignment
# Place in R/reference_assignment.R

build_tcga_reference <- function(M_pam50, V_tcga, tcga_se, tcga_sub_col_opts = c("PAM50","Subtype","BRCA_Subtype","molecular_subtype"),
                                 subtype_labels = c("HER2-high","HER2-low","LumA","LumB","Basal"), k_per_subtype = 20, outdir = NULL) {
  tcga_sub_col <- tcga_sub_col_opts[tcga_sub_col_opts %in% colnames(colData(tcga_se))][1]
  M_pam50_tcga <- M_pam50[, colnames(V_tcga), drop = FALSE]
  X_tcga <- t(scale(M_pam50_tcga))
  d_tcga  <- dist(X_tcga, method = "euclidean")
  hc_tcga <- hclust(d_tcga, method = "average")
  cl_tcga <- cutree(hc_tcga, k = 5)
  if (!is.na(tcga_sub_col)) {
    lab_tcga_full <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
    sub_lvls <- sort(unique(na.omit(lab_tcga_full)))
    sub_means_tcga <- sapply(sub_lvls, function(s) rowMeans(M_pam50_tcga[, lab_tcga_full == s, drop = FALSE], na.rm = TRUE))
    cl_means_tcga <- sapply(1:5, function(k) rowMeans(M_pam50_tcga[, cl_tcga == k, drop = FALSE], na.rm = TRUE))
    g_common <- intersect(rownames(sub_means_tcga), rownames(cl_means_tcga))
    corr_mat <- cor(cl_means_tcga[g_common, , drop = FALSE],
                    sub_means_tcga[g_common, , drop = FALSE], method = "spearman")
    cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])
    cl_to_sub_named <- setNames(subtype_labels[match(cl_to_sub, sub_lvls)], 1:5)
  } else {
    km5 <- kmeans(X_tcga, centers = 5, nstart = 50)
    km_cent <- t(sapply(1:5, function(k) colMeans(X_tcga[km5$cluster == k, , drop = FALSE], na.rm = TRUE)))
    ord <- order(prcomp(km_cent)$x[,1])
    cl_to_sub_named <- setNames(subtype_labels[ord], 1:5)
  }
  subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)], levels = subtype_labels)
  sel_list <- list(); summary_counts <- c()
  for (sub in subtype_labels) {
    idx <- which(subtype_tcga == sub)
    if (length(idx) == 0) next
    cen <- get_centroid(X_tcga, idx)
    nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = min(k_per_subtype, length(idx)))
    nearest$subtype <- sub
    sel_list[[sub]] <- nearest
    summary_counts[sub] <- nrow(nearest)
  }
  sel_df <- dplyr::bind_rows(sel_list) %>% dplyr::relocate(subtype, .before = sample)
  sel_df <- sel_df %>% group_by(subtype) %>% slice_head(n = k_per_subtype) %>% ungroup()
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    write.csv(sel_df, file.path(outdir, "tcga_top20_per_subtype.csv"), row.names = FALSE)
    saveRDS(M_pam50_tcga[, sel_df$sample, drop = FALSE], file.path(outdir, "PAM50_tcga_top20_per_subtype.rds"))
  }
  list(sel_df = sel_df, subtype_tcga = subtype_tcga, cl_to_sub_named = cl_to_sub_named)
}

assign_dsmz_from_tcga <- function(M_pam50, sel_df, subtype_labels, outdir = NULL) {
  ref_samples <- sel_df$sample
  M_ref <- M_pam50[, ref_samples, drop = FALSE]
  centroids <- sapply(subtype_labels, function(s) {
    ss <- sel_df$sample[sel_df$subtype == s]
    if (!length(ss)) return(rep(NA_real_, nrow(M_ref)))
    rowMeans(M_ref[, ss, drop = FALSE], na.rm = TRUE)
  })
  rownames(centroids) <- rownames(M_ref)
  M_pam50_dsmz <- M_pam50[, !colnames(M_pam50) %in% ref_samples, drop = FALSE]
  g <- intersect(rownames(centroids), rownames(M_pam50_dsmz))
  if (length(g) < 10) stop("[FATAL] Too few overlapping PAM50 genes between reference and DSMZ.")
  assign_one <- function(x, C) {
    cs <- apply(C, 2, function(cen) suppressWarnings(cor(x, cen, method = "spearman", use = "pairwise.complete.obs")))
    lab <- names(which.max(cs))
    score <- max(cs, na.rm = TRUE)
    c(label = lab, score = score)
  }
  dsmz_assign <- t(apply(M_pam50_dsmz[g, , drop = FALSE], 2, assign_one, C = centroids[g, , drop = FALSE]))
  dsmz_assign <- as.data.frame(dsmz_assign, stringsAsFactors = FALSE)
  dsmz_assign$sample <- rownames(dsmz_assign)
  dsmz_assign$label  <- factor(dsmz_assign$label, levels = subtype_labels)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    write.csv(dsmz_assign, file.path(outdir, "dsmz_pred_subtypes_from_tcga_top20.csv"), row.names = FALSE)
  }
  dsmz_assign
}