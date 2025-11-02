Awesome—here’s a clean, self-contained add-on that does exactly what you asked:

1. **Build a TCGA reference** in PAM50 space by:

   * clustering TCGA (PAM50, k=5 → mapped to PAM50 labels),
   * computing a centroid per subtype,
   * selecting the **20 nearest TCGA samples** to each subtype centroid,
   * saving both the list of samples and the PAM50 matrix of these references.

2. **Assign DSMZ** subtypes **independently** using **only** those TCGA references:

   * compute one **reference centroid per subtype** from the top-20 samples,
   * **label each DSMZ sample** by **maximum Spearman correlation** to those five centroids,
   * save a table and a quick UMAP plot.

> Drop this block **after** your script has already created `M_pam50`, `V_tcga`, `V_dsmz`, and `tcga_se` (you already have those earlier in your pipeline).

---

```r
## =====================================================================
## (A) TCGA-only balanced reference: 20 per PAM50 subtype (PAM50 space)
## =====================================================================
set.seed(42)

# PAM50 matrix for TCGA samples (genes x samples)
M_pam50_tcga <- M_pam50[, colnames(V_tcga), drop = FALSE]

# Scale samples in PAM50 space for clustering/nearest calculations
X_tcga <- t(scale(M_pam50_tcga))   # samples x genes (scaled)
if (nrow(X_tcga) < 5) {
  stop(sprintf("[FATAL] Too few TCGA samples in PAM50 space (%d) for clustering.", nrow(X_tcga)))
}

# Hierarchical clustering (k=5)
d_tcga  <- dist(X_tcga, method = "euclidean")
hc_tcga <- hclust(d_tcga, method = "average")
cl_tcga <- cutree(hc_tcga, k = 5)

# Find a TCGA subtype column if present
tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]
subtype_labels <- c("HER2-high","HER2-low","LumA","LumB","Basal")

# Map HC clusters to PAM50 names using correlation to TCGA subtype centroids (if annotation exists)
if (!is.na(tcga_sub_col)) {
  lab_tcga_full <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga_full)))

  # centroids in PAM50 space (genes as rows)
  sub_means_tcga <- sapply(sub_lvls, function(s) {
    rowMeans(M_pam50_tcga[, lab_tcga_full == s, drop = FALSE], na.rm = TRUE)
  })

  cl_means_tcga <- sapply(1:5, function(k) {
    rowMeans(M_pam50_tcga[, cl_tcga == k, drop = FALSE], na.rm = TRUE)
  })

  g_common <- intersect(rownames(sub_means_tcga), rownames(cl_means_tcga))
  corr_mat <- cor(cl_means_tcga[g_common, , drop = FALSE],
                  sub_means_tcga[g_common, , drop = FALSE],
                  method = "spearman")

  cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])
  cl_to_sub_named <- setNames(subtype_labels[match(cl_to_sub, sub_lvls)], 1:5)
} else {
  # If no annotation exists, give stable names via simple ordering heuristic
  km5 <- kmeans(X_tcga, centers = 5, nstart = 50)
  km_cent <- t(sapply(1:5, function(k) colMeans(X_tcga[km5$cluster == k, , drop = FALSE], na.rm = TRUE)))
  ord <- order(prcomp(km_cent)$x[,1])
  cl_to_sub_named <- setNames(subtype_labels[ord], 1:5)
}

# Final subtype per TCGA sample (factor with fixed order)
subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)], levels = subtype_labels)
cat("[TCGA ref] Cluster→subtype size:\n"); print(table(subtype_tcga))

# ---- Helpers for centroid + nearest selection in scaled space ----
get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
get_k_nearest <- function(X, centroid, k = 20) {
  d <- sqrt(rowSums((X - matrix(centroid, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
  ord <- order(d, decreasing = FALSE)
  data.frame(sample = rownames(X)[ord], dist = d[ord], stringsAsFactors = FALSE)
}

# ---- Select up to top-20 per subtype (closest to the subtype centroid) ----
sel_list <- list(); summary_counts <- c()
for (sub in subtype_labels) {
  idx <- which(subtype_tcga == sub)
  if (length(idx) == 0) {
    cat(sprintf("[WARN] No TCGA samples for subtype '%s' — skipping.\n", sub))
    next
  }
  cen <- get_centroid(X_tcga, idx)                               # centroid in *scaled* space
  nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = min(20, length(idx)))
  nearest$subtype <- sub
  sel_list[[sub]] <- nearest
  summary_counts[sub] <- nrow(nearest)
}
sel_df <- dplyr::bind_rows(sel_list) %>% dplyr::relocate(subtype, .before = sample)

# Save selection
dir.create(file.path(outdir, "tcga_reference_selection"), showWarnings = FALSE, recursive = TRUE)
write.csv(sel_df, file.path(outdir, "tcga_reference_selection", "tcga_top20_per_subtype.csv"), row.names = FALSE)

# Save PAM50-only expression matrix for these references (genes x selected samples)
ref_samples <- sel_df$sample
M_pam50_tcga_ref <- M_pam50_tcga[, ref_samples, drop = FALSE]
saveRDS(M_pam50_tcga_ref, file.path(outdir, "tcga_reference_selection", "PAM50_tcga_top20_per_subtype.rds"))

# Quick UMAP to visualize chosen references
set.seed(42)
emb_ref <- uwot::umap(X_tcga, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb_ref <- as.data.frame(emb_ref); colnames(emb_ref) <- c("UMAP1","UMAP2")
emb_ref$sample  <- rownames(X_tcga)
emb_ref$subtype <- subtype_tcga
emb_ref$chosen  <- ifelse(emb_ref$sample %in% ref_samples, "Top20", "Other")
p_ref <- ggplot(emb_ref, aes(UMAP1, UMAP2, color = subtype, shape = chosen)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("TCGA BRCA (PAM50) — clusters & selected Top20 per subtype")
ggsave(file.path(outdir, "tcga_reference_selection", "umap_tcga_reference_top20.pdf"), p_ref, width = 7, height = 6)

# Export plain list by subtype
by_sub <- split(sel_df$sample, sel_df$subtype)
dumpfile <- file.path(outdir, "tcga_reference_selection", "tcga_top20_samples_by_subtype.txt")
con_ref <- file(dumpfile, "wt"); on.exit(close(con_ref), add = TRUE)
for (nm in names(by_sub)) {
  writeLines(sprintf("[%s]", nm), con_ref)
  writeLines(by_sub[[nm]], con_ref)
  writeLines("", con_ref)
}
cat(sprintf("[TCGA ref] Wrote selection to: %s\n", file.path(outdir, "tcga_reference_selection")))

## =====================================================================
## (B) DSMZ assignment from TCGA references (independent label transfer)
## =====================================================================

# Load references (genes x selected TCGA samples) and their grouping
M_ref <- readRDS(file.path(outdir, "tcga_reference_selection", "PAM50_tcga_top20_per_subtype.rds"))
sel_info <- read.csv(file.path(outdir, "tcga_reference_selection", "tcga_top20_per_subtype.csv"), stringsAsFactors = FALSE)

# Keep only samples that really exist in the matrix
sel_info <- sel_info[sel_info$sample %in% colnames(M_ref), ]
stopifnot(nrow(sel_info) > 0)

# Build one centroid per subtype in PAM50 space (genes as rows)
centroids <- sapply(subtype_labels, function(s) {
  ss <- sel_info$sample[sel_info$subtype == s]
  if (!length(ss)) return(rep(NA_real_, nrow(M_ref)))
  rowMeans(M_ref[, ss, drop = FALSE], na.rm = TRUE)
})
rownames(centroids) <- rownames(M_ref)  # genes
colnames(centroids) <- subtype_labels

# Prepare DSMZ PAM50 (genes x DSMZ samples)
M_pam50_dsmz <- M_pam50[, colnames(V_dsmz), drop = FALSE]

# Use only common PAM50 genes
g <- intersect(rownames(centroids), rownames(M_pam50_dsmz))
if (length(g) < 10) stop("[FATAL] Too few overlapping PAM50 genes between reference and DSMZ.")

# Spearman correlation for robustness (rank-based, scale-free)
assign_one <- function(x, C) {
  cs <- apply(C, 2, function(cen) suppressWarnings(cor(x, cen, method = "spearman", use = "pairwise.complete.obs")))
  lab <- names(which.max(cs))
  score <- max(cs, na.rm = TRUE)
  c(label = lab, score = score)
}

# Compute assignment per DSMZ sample
dsmz_assign <- t(apply(M_pam50_dsmz[g, , drop = FALSE], 2, assign_one, C = centroids[g, , drop = FALSE]))
dsmz_assign <- as.data.frame(dsmz_assign, stringsAsFactors = FALSE)
dsmz_assign$sample <- rownames(dsmz_assign)
dsmz_assign$label  <- factor(dsmz_assign$label, levels = subtype_labels)

# Save table
dir.create(file.path(outdir, "dsmz_from_tcga_reference"), showWarnings = FALSE, recursive = TRUE)
write.csv(dsmz_assign, file.path(outdir, "dsmz_from_tcga_reference", "dsmz_pred_subtypes_from_tcga_top20.csv"), row.names = FALSE)

# Quick DSMZ-only UMAP (PAM50 space), color by predicted subtype
set.seed(42)
X_dsmz <- t(scale(M_pam50_dsmz[g, , drop = FALSE]))  # samples x genes (scaled)
emb_d <- uwot::umap(X_dsmz, n_neighbors = 10, min_dist = 0.3, metric = "cosine")
emb_d <- as.data.frame(emb_d); colnames(emb_d) <- c("UMAP1","UMAP2")
emb_d$sample <- rownames(X_dsmz)
emb_d$pred   <- factor(dsmz_assign$label[match(emb_d$sample, dsmz_assign$sample)], levels = subtype_labels)

p_pred <- ggplot(emb_d, aes(UMAP1, UMAP2, color = pred)) +
  geom_point(alpha = 0.95, size = 2.1) + theme_bw() +
  ggtitle("DSMZ (PAM50) — predicted PAM50 subtype (nearest TCGA centroid)")
ggsave(file.path(outdir, "dsmz_from_tcga_reference", "umap_dsmz_predicted_labels.pdf"), p_pred, width = 6.5, height = 5.5)

# Also export a simple 2-column mapping for downstream use
write.table(
  dsmz_assign[, c("sample","label")],
  file = file.path(outdir, "dsmz_from_tcga_reference", "dsmz_predicted_label_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("[DSMZ assignment] Done. Results in: ", file.path(outdir, "dsmz_from_tcga_reference"), "\n")
```

### Notes / why this is robust

* Everything is computed in **PAM50 space**, which is the cleanest axis for breast subtypes.
* DSMZ labels are assigned **only** from **TCGA top-20 centroids**—no leakage from the rest of TCGA.
* **Spearman correlation** makes the transfer stable to scale and mild batch effects.
* If a subtype has <20 TCGA samples available, the code just takes as many as exist and moves on (with a warning).

If you’d like, we can wrap parts (A) and (B) into two functions (`build_tcga_reference()` and `assign_dsmz_from_reference()`) and wire them behind a flag so your Snakemake rule can flip between “build reference” and “label DSMZ” in one step.
