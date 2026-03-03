ens2sym_union <- map_ensembl_to_symbol(expr_mat = DSMZ_counts)
pam50_DSMZ <- subset_to_pam50(DSMZ_counts, ens2sym_union, pam50_symbols, "mean")
saveRDS(pam50_DSMZ, file.path(config$outdir, "DSMZ-BRCA_PAM50_matrix.rds"))

# -------------------------------------------------------------------------
# 1. GENE CLUSTERING (genes x samples)
# -------------------------------------------------------------------------
cat("\n[1/2] Gene-level Dynamic Clustering clustering...\n")
hc_genes <- run_hierarchical_clustering(
  Exp_mat = pam50_DSMZ,
  title   = "DSMZ-BRCA_PAM50_genes",
  outdir  = config$outdir
)

# Heatmap: Genes (rows) – Dynamic clusters
save_heatmap_pub(
  mat = pam50_DSMZ,
  annotation_row = data.frame(Cluster = hc_genes$dynamic$clusters_labeled,
                              row.names = rownames(pam50_DSMZ)),
  out_pdf = file.path(config$outdir, "DSMZ-BRCA_PAM50_heatmap_genes_dynamic.pdf"),
  main = "PAM50 Genes – Unsupervised Dynamic Clustering",
  show_rownames = TRUE,
  cluster_rows = hc_genes$dynamic$hc,
  cluster_cols = TRUE,
  width = 10, height = 8
)

cat("\n[1/2] Gene-level optimal clustering...\n")
hc_genes <- run_hierarchical_clustering(
  Exp_mat = pam50_DSMZ,
  title = "DSMZ-BRCA_PAM50_genes",
  max_k = 5,
  use_optimal_k = TRUE,
  outdir = config$outdir
)

# Heatmap: Optimal k
save_heatmap_pub(
  mat = pam50_DSMZ,
  annotation_row = data.frame(Cluster = hc_genes$optimal$clusters_labeled, row.names = rownames(pam50_DSMZ)),
  out_pdf = file.path(config$outdir, "DSMZ-BRCA_PAM50_heatmap_genes_optimal.pdf"),
  main = sprintf("PAM50 Genes – Optimal k=%d", hc_genes$optimal$k),
  show_rownames = TRUE,
  cluster_rows = hc_genes$optimal$hc
)

# === GENE CLUSTERING (force k=4 for biology) ===
cat("\n[1/2] Gene-level clustering (forced k=4)...\n")
hc_genes_fixed <- run_fixed_hc(pam50_DSMZ, k = 4, seed = 42)
clusters_gene_k4 <- hc_genes_fixed$clusters

save_heatmap_pub(
  mat = pam50_DSMZ,
  annotation_row = data.frame(Cluster = clusters_gene_k4, row.names = rownames(pam50_DSMZ)),
  out_pdf = file.path(config$outdir, "DSMZ-BRCA_PAM50_heatmap_genes_fixed_k4.pdf"),
  main = "PAM50 Genes – Biologically Informed k=4",
  show_rownames = TRUE,
  cluster_rows = hc_genes_fixed$hc
)



# -------------------------------------------------------------------------
# 2. PATIENT CLUSTERING (samples x genes)
# -------------------------------------------------------------------------
cat("\n[2/2] Patient-level clustering...\n")
pam50_DSMZ_t <- t(pam50_DSMZ)

hc_patients_dy <- run_hierarchical_clustering(
  Exp_mat = pam50_DSMZ_t,
  title   = "DSMZ-BRCA_PAM50_patients",
  outdir  = config$outdir
)

# Heatmap: Patients (rows) – Dynamic clusters
save_heatmap_pub(
  mat = pam50_DSMZ_t,
  annotation_row = data.frame(Cluster = hc_patients_dy$dynamic$clusters_labeled,
                              row.names = rownames(pam50_DSMZ_t)),
  out_pdf = file.path(config$outdir, "DSMZ-BRCA_PAM50_heatmap_patients_dynamic.pdf"),
  main = "DSMZ-BRCA Patients – Unsupervised Dynamic Clustering",
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = hc_patients_dy$dynamic$hc,
  cluster_cols = TRUE,
  width = 12, height = 10
)

# optimal k
hc_patients_fixed <- run_hierarchical_clustering(
  Exp_mat = pam50_DSMZ_t,
  title = "DSMZ-BRCA_PAM50_patients",
  max_k = 6,
  use_optimal_k = TRUE,
  outdir = config$outdir
)

save_heatmap_pub(
  mat = pam50_DSMZ_t,
  annotation_row = data.frame(Cluster = hc_patients_fixed$optimal$clusters_labeled, row.names = rownames(pam50_DSMZ_t)),
  out_pdf = file.path(config$outdir, "DSMZ-BRCA_PAM50_heatmap_patients_optimal.pdf"),
  main = sprintf("DSMZ-BRCA Patients – Optimal k=%d", hc_patients_fixed$optimal$k),
  show_rownames = FALSE,
  cluster_rows = hc_patients_fixed$optimal$hc
)

# Plot silhouette curve
# SAFE CALL
if (!is.null(hc_patients_fixed$optimal) && 
    !is.null(hc_patients_fixed$optimal$all_ks) &&
    length(hc_patients$optimal$all_ks) > 0) {
  plot_silhouette_curve(
    hc_patients_fixed$optimal,
    file.path(config$outdir, "DSMZ-BRCA_PAM50_silhouette_vs_k_patients.pdf")
  )
}
# -------------------------------------------------------------------------
# Done
# -------------------------------------------------------------------------
message("\n[SUCCESS]  heatmaps saved in:\n    ", config$outdir)





