# nolint start
# Load required libraries with suppressed startup messages for cleaner output
suppressPackageStartupMessages({
  library(tidyverse)          # Data manipulation and visualization (dplyr, ggplot2, etc.)
  library(SummarizedExperiment) # Container for genomic data (TCGA data format)
  library(DESeq2)             # Differential expression analysis and VST normalization
  library(matrixStats)        # Fast matrix operations (rowMeans, rowSds, rowIQRs)
  library(pheatmap)           # Pretty heatmaps for visualization
  library(uwot)               # UMAP dimensionality reduction
  library(cluster)            # Clustering algorithms (silhouette analysis)
  library(fpc)                # Calinski-Harabasz index
  library(clusterCrit)        # Davies-Bouldin index
  library(patchwork)
  library(genefilter)
})

# ---------- CONFIG ----------
# Define input file paths
tcga_se_rds   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"  # TCGA BRCA RNA-seq data
dsmz_rds      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"                                      # DSMZ cell line RNA-seq data
dsmz_meta_csv <- "/home/chu25/data/dsmz/DSMZ_metadata.csv"                                       # DSMZ sample metadata

# Define output directory and create it if it doesn't exist
outdir <- "/home/chu25/dsmz/results/brca_vst_dual"  # Where to save all results
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)  # Create directory structure

# Subsetting filters (set to NULL to skip filtering)
tcga_project_filter <- c("TCGA-BRCA")   # TCGA projects to keep; set NULL to skip filtering
dsmz_organ_filter   <- c("Breast")      # DSMZ organ(s) to keep; set NULL to skip filtering

# PAM50 marker genes for breast cancer subtyping
# These are key biomarkers: ESR1=ER, PGR=PR, ERBB2=HER2, MKI67=proliferation, KRT*=keratins
marker_symbols <- c("ESR1","PGR","ERBB2","MKI67","KRT5","KRT14","KRT17","FOXC1","GRB7",  # Luminal/basal markers
                    "BCL2","EGFR","CCNB1","BIRC5","MYBL2","KIF2C","KRT8","KRT18","GATA3","XBP1")  # Additional subtype markers

# ---------- HELPERS ----------
# Strip Ensembl version suffixes (e.g., "ENSG000001...12" → "ENSG000001...")
# This ensures compatibility between TCGA and DSMZ gene IDs
strip_ensver <- function(x) sub("\\..*$","", x)  # Remove everything after first dot

# Helper to safely convert to data.frame
as_df <- function(x) {
  if (inherits(x, "data.frame")) return(as.data.frame(x, stringsAsFactors = FALSE))
  as.data.frame(x, stringsAsFactors = FALSE)
}

# --- helper: unique Ensembl -> one symbol row ---
make_ens2sym <- function(dsmz_raw) {
  df <- as_df(dsmz_raw)
  ens2sym <- df %>%
    transmute(
      ensembl = strip_ensver(Ensembl_ID),
      symbol  = as.character(gene_name)
    ) %>%
    filter(!is.na(ensembl) & nzchar(ensembl)) %>%
    arrange(ensembl, dplyr::desc(nchar(symbol)), symbol) %>%
    distinct(ensembl, .keep_all = TRUE) %>%
    mutate(symbol_unique = make.unique(ifelse(is.na(symbol) | symbol == "", ensembl, symbol)))
  ens2sym <- as.data.frame(ens2sym, stringsAsFactors = FALSE)  # ensure not tibble
  rownames(ens2sym) <- ens2sym$ensembl
  ens2sym
}

# robust HVGs (IQR) and safe z-scoring
top_var_genes <- function(mat, n = 3000) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}
zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}
safe_pdf <- function(path, expr) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path); on.exit(dev.off(), add = TRUE); force(expr)
}

# Convert DSMZ raw data to numeric count matrix, handle duplicates
# Takes data.frame with annotation columns and converts to matrix format
build_dsmz_matrix <- function(dsmz_raw) {
  stopifnot(all(c("Ensembl_ID","gene_name") %in% names(dsmz_raw)))  # Check required columns exist
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version")  # Define annotation column names
  sample_cols <- setdiff(colnames(dsmz_raw), annot)  # Identify sample (count) columns
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(v) as.numeric(as.character(v)))  # Convert to numeric
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])  # Convert to matrix
  rownames(M) <- strip_ensver(dsmz_raw$Ensembl_ID)  # Set gene IDs as row names
  if (any(duplicated(rownames(M)))) M <- rowsum(M, rownames(M), reorder = TRUE)  # Sum duplicate genes
  M  # Return processed matrix
}

# ---- LOGGING HELPERS ----
log_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark=","), format(ns, big.mark=",")))
}

write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con)
}

# ---------- LOAD + FILTER ----------
cat("[INFO] Starting pipeline...\n")

# ---- load
tcga_se <- readRDS(tcga_se_rds)  # Load the full TCGA dataset
stopifnot("project_id" %in% colnames(colData(tcga_se)))  # Ensure project_id column exists
dsmz_raw  <- readRDS(dsmz_rds)  # Load DSMZ count data
dsmz_meta <- read.csv(dsmz_meta_csv, check.names = FALSE)  # Load DSMZ metadata

# ---- harmonize/build
tcga_counts <- assay(tcga_se)  # Extract count matrix from SummarizedExperiment
rownames(tcga_counts) <- strip_ensver(rownames(tcga_counts))  # Clean gene IDs (remove versions)
if (any(duplicated(rownames(tcga_counts)))) {  # Check for duplicate gene IDs
  tcga_counts <- rowsum(tcga_counts, rownames(tcga_counts), reorder = TRUE)  # Sum duplicates
}
storage.mode(tcga_counts) <- "double"

dsmz_counts <- build_dsmz_matrix(dsmz_raw)  # Convert to matrix format
storage.mode(dsmz_counts) <- "double"

# Capture pre-subset snapshots for accurate dimension reporting
tcga_counts_raw    <- tcga_counts
dsmz_counts_raw    <- dsmz_counts
dsmz_meta_raw      <- dsmz_meta

# ---- initial dims
log_dims("TCGA (raw)", tcga_counts_raw)
log_dims("DSMZ (raw)", dsmz_counts_raw)

# ---- align DSMZ metadata and SUBSET to organ filter (e.g., Breast)
# Create sample_id column and align metadata with counts
dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)
matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))
dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched)
dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop = FALSE]
stopifnot("organ" %in% colnames(dsmz_meta))

if (!is.null(dsmz_organ_filter)) {
  keep_ids <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
  dsmz_meta   <- dsmz_meta[dsmz_meta$sample_id %in% keep_ids, , drop = FALSE]
  dsmz_counts <- dsmz_counts[, keep_ids, drop = FALSE]
}

log_dims("DSMZ (after organ subset)", dsmz_counts)

# Guardrails for subsetting
if (!is.null(dsmz_organ_filter) && ncol(dsmz_counts) == 0L) {
  stop(sprintf("[FATAL] DSMZ organ filter %s left 0 samples.", paste(dsmz_organ_filter, collapse=", ")))
}

# ---- OPTIONAL: subset TCGA to BRCA only (recommended for BRCA focus)
if (!is.null(tcga_project_filter)) {
  # tcga_se has colData with project_id like "TCGA-BRCA"
  keep_tcga <- colnames(tcga_counts) %in% colnames(assay(tcga_se)) &
               as.character(colData(tcga_se)[colnames(tcga_counts), "project_id"]) %in% tcga_project_filter
  tcga_counts <- tcga_counts[, keep_tcga, drop = FALSE]
}

log_dims("TCGA (after project subset)", tcga_counts)

# Guardrails for subsetting
if (!is.null(tcga_project_filter) && ncol(tcga_counts) == 0L) {
  stop(sprintf("[FATAL] TCGA project filter %s left 0 samples.", paste(tcga_project_filter, collapse=", ")))
}

# ---- find common genes
common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
cat(sprintf("[INFO] Shared genes: %d\n", length(common)))
if (length(common) < 1000) cat("[WARN] Few shared genes; check Ensembl IDs\n")

tcga_counts <- tcga_counts[common, , drop = FALSE]
dsmz_counts <- dsmz_counts[common, , drop = FALSE]

log_dims("TCGA (post-intersect)", tcga_counts)
log_dims("DSMZ (post-intersect)", dsmz_counts)

# ---- merged raw matrix
Xc_raw <- cbind(tcga_counts, dsmz_counts)

# safety: fix any duplicate gene names (rare, but safe)
if (any(duplicated(rownames(Xc_raw)))) {
  Xc_raw <- rowsum(Xc_raw, rownames(Xc_raw), reorder = TRUE)
}

log_dims("Merged (Xc_raw)", Xc_raw)

# ---- HARD STOPS if anything is empty
if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA subset has 0 samples.")
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ subset has 0 samples.")
if (nrow(Xc_raw)        == 0L) stop("[FATAL] No common genes after intersection.")
if (ncol(Xc_raw)        == 0L) stop("[FATAL] Merged matrix has 0 samples.")

# ---- write a compact DIM report to disk
dim_report <- file.path(outdir, "dimension_report.txt")
con <- file(dim_report, open = "wt")
on.exit(close(con), add = TRUE)

write_dims_line(con, "TCGA (raw)", tcga_counts_raw)
write_dims_line(con, "DSMZ (raw)", dsmz_counts_raw)

write_dims_line(con, "DSMZ (after organ subset)", dsmz_counts)
write_dims_line(con, "TCGA (after project subset)", tcga_counts)

writeLines(sprintf("Shared genes\t%d", length(common)), con)

write_dims_line(con, "TCGA (post-intersect)", tcga_counts)
write_dims_line(con, "DSMZ (post-intersect)", dsmz_counts)

write_dims_line(con, "Merged (Xc_raw)", Xc_raw)

# ---- proceed as before (pipeline continues here)
batch <- factor(c(rep("TCGA", ncol(tcga_counts)), rep("DSMZ", ncol(dsmz_counts))))
all_counts <- Xc_raw  # Use the merged matrix we already created

# ---------- JOINT VST NORMALIZATION ----------
# Combine TCGA and DSMZ counts for joint normalization
coldata <- tibble(  # Create metadata for DESeq2
  sample = colnames(all_counts),  # Sample IDs
  dataset = ifelse(colnames(all_counts) %in% colnames(tcga_counts), "TCGA", "DSMZ")  # Dataset labels
) %>% column_to_rownames("sample")  # Convert to row names format

# Apply DESeq2 VST normalization (blind=TRUE for unbiased joint normalization)
# Robustness check for duplicate gene names
if (any(duplicated(rownames(all_counts)))) {
  all_counts <- rowsum(all_counts, rownames(all_counts), reorder = TRUE)
}
dds <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ 1)  # Create DESeq2 object
vsd <- vst(dds, blind = TRUE)  # Apply variance-stabilizing transformation
V <- assay(vsd)  # Extract VST-normalized matrix (genes x samples)
V_tcga <- V[, colnames(tcga_counts), drop = FALSE]  # Split back to TCGA samples
V_dsmz <- V[, colnames(dsmz_counts), drop = FALSE]  # Split back to DSMZ samples

# Sanity check: ensure column names exist
stopifnot(!is.null(colnames(V_tcga)), !is.null(colnames(V_dsmz)))

# Save VST matrices for future use
saveRDS(V_tcga, file.path(outdir, "VST_tcga_brca.rds"))  # Save TCGA VST matrix
saveRDS(V_dsmz, file.path(outdir, "VST_dsmz_brca.rds"))  # Save DSMZ VST matrix

# ---------- (1) UNSUPERVISED DISCOVERY CLUSTERING (HVG-BASED) ----------
# Discovery first, labels later: cluster on HVGs from joint VST matrix

# 1) feature space for discovery (HVGs from joint VST)
set.seed(42)
hvgs <- top_var_genes(V, n = 3000)                            # discovery features (no PAM50)
X_discovery <- t(V[hvgs, , drop = FALSE])                     # samples x genes
X_scaled    <- scale(X_discovery)                             # standardize per gene

# Save the actual HVG feature set for reproducibility
writeLines(hvgs, file.path(outdir, "HVGs_used_for_discovery_3k.txt"))

# 2) scan k with multiple criteria + plots (2…8)
k_grid <- 2:8
res_list <- list()
metrics <- lapply(k_grid, function(k) {
  km <- kmeans(X_scaled, centers = k, nstart = 50, iter.max = 100)
  # WCSS (elbow)
  wcss <- sum(km$withinss)
  # Silhouette
  sil  <- silhouette(km$cluster, dist(X_scaled))
  sil_mean <- mean(sil[, "sil_width"])
  # CH index
  ch <- tryCatch(calinhara(X_scaled, km$cluster), error = function(e) NA_real_)
  # Davies–Bouldin (lower is better)
  db <- tryCatch({
    m <- intCriteria(as.matrix(X_scaled), as.integer(km$cluster), "davies_bouldin")
    as.numeric(m$davies_bouldin)
  }, error = function(e) NA_real_)
  res_list[[as.character(k)]] <<- list(km = km, sil = sil)
  data.frame(k = k, WCSS = wcss, Silhouette = sil_mean, CH = ch, DB = db)
})
metrics <- do.call(rbind, metrics)

# save metrics
write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), row.names = FALSE)

# elbow + metrics panel
safe_pdf(file.path(outdir, "kmeans_elbow_HVG.pdf"), {
  plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", ylab = "Within-Cluster SS", main = "Elbow (HVG space)")
})
safe_pdf(file.path(outdir, "kmeans_metrics_panel_HVG.pdf"), {
  par(mfrow = c(2,2))
  plot(metrics$k, metrics$Silhouette, type = "b", xlab = "k", ylab = "Avg silhouette", main = "Silhouette (higher is better)")
  plot(metrics$k, metrics$CH,         type = "b", xlab = "k", ylab = "Calinski-Harabasz", main = "Calinski-Harabasz (higher)")
  plot(metrics$k, metrics$DB,         type = "b", xlab = "k", ylab = "Davies-Bouldin",    main = "Davies-Bouldin (lower)")
  plot(metrics$k, metrics$WCSS,       type = "b", xlab = "k", ylab = "WCSS",              main = "Elbow (lower)")
  par(mfrow = c(1,1))
})

# choose best k: max silhouette, tie-break by CH
best_k <- with(metrics, {
  cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)]
  if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand
})
message(sprintf("[kmeans] best k by silhouette→CH = %s", best_k))

# 3) final fit, save labels, silhouette plot
km_final   <- res_list[[as.character(best_k)]]$km
sil_final  <- res_list[[as.character(best_k)]]$sil
clusters   <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k)))

# annotate which samples are DSMZ vs TCGA (consistent with your earlier code)
sample_is_tcga <- colnames(V) %in% colnames(V_tcga)
dataset_lab <- ifelse(sample_is_tcga, "TCGA", "DSMZ")

cluster_df <- tibble(sample = rownames(X_scaled),
                     dataset = dataset_lab,
                     cluster = as.character(clusters))
write.csv(cluster_df, file.path(outdir, "clusters_HVG_kmeans.csv"), row.names = FALSE)

safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_HVG.pdf", best_k)), {
  plot(sil_final, main = sprintf("Silhouette (k = %d, HVG space)", best_k), cex.names = 0.8)
})

# Save per-k silhouette plots for discovery rigor
sil_dir <- file.path(outdir, "silhouettes_per_k")
dir.create(sil_dir, showWarnings = FALSE, recursive = TRUE)
for (k in k_grid) {
  silk <- res_list[[as.character(k)]]$sil
  if (!is.null(silk)) {
    pdf(file.path(sil_dir, sprintf("silhouette_k=%d.pdf", k)), width = 8, height = 6)
    plot(silk, main = sprintf("Silhouette (k = %d, HVG space)", k), cex.names = 0.7)
    dev.off()
  }
}

# Print clustering results
cat("Optimal k:", best_k, "with mean silhouette:", round(mean(sil_final[, "sil_width"]), 3), "\n")
print(table(clusters))  # Print cluster sizes

# 4) UMAP colored by discovered clusters (discovery features)
set.seed(42)
emb <- uwot::umap(X_discovery, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2")
emb$sample  <- rownames(X_discovery)
emb$dataset <- dataset_lab
emb$cluster <- clusters

p_umap <- ggplot(emb, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (HVG discovery space)")

ggsave(file.path(outdir, "umap_discovery_HVG.pdf"), p_umap, width = 7, height = 6)
write.csv(emb, file.path(outdir, "umap_discovery_HVG_coords.csv"), row.names = FALSE)

# 5) post hoc biological interpretation (does NOT influence clustering)

# (a) marker heatmap by discovered clusters
# symbol ↔ ensembl map from DSMZ table (you already built ens2sym)
marker_symbols <- c("ESR1","PGR","ERBB2","MKI67","KRT5","KRT14","KRT17","FOXC1","GRB7",
                    "BCL2","EGFR","CCNB1","BIRC5","MYBL2","KIF2C","KRT8","KRT18","GATA3","XBP1")
ens2sym <- make_ens2sym(dsmz_raw)

# Quick sanity log (optional)
ens_vec <- strip_ensver(as_df(dsmz_raw)$Ensembl_ID)
dup_ids <- sum(duplicated(ens_vec))                          # duplicate *rows*
dup_unique <- sum(table(ens_vec) > 1, na.rm = TRUE)          # duplicate *IDs*
cat(sprintf("[INFO] DSMZ Ensembl duplicate rows: %d | duplicate IDs: %d (collapsed)\n",
            dup_ids, dup_unique))

marker_ens <- rownames(ens2sym)[tolower(ens2sym$symbol) %in% tolower(marker_symbols)]
marker_ens <- intersect(unique(marker_ens), rownames(V))

if (length(marker_ens) > 0) {
  V_mark <- V[marker_ens, , drop = FALSE]
  
  # guard 1: if any rownames still collide, make them unique to satisfy pheatmap
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  
  # guard 2: use pretty labels (unique symbols) on the heatmap if you want
  lab <- ens2sym[rownames(V_mark), "symbol_unique", drop = TRUE]
  lab[is.na(lab) | lab == ""] <- rownames(V_mark)
  rownames(V_mark) <- lab
  
  # order samples by (dataset, cluster)
  ord <- order(dataset_lab, clusters)
  ann <- data.frame(dataset = dataset_lab, cluster = clusters,
                    row.names = colnames(V_mark))
  safe_pdf(file.path(outdir, "heatmap_marker_discovery_clusters.pdf"), {
    pheatmap(V_mark[, ord, drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_col = ann[ord, , drop = FALSE],
             main = "Marker expression (post hoc, discovery clusters)")
  })
} else {
  message("[marker heatmap] No marker genes found in V — skipping.")
}

# (b) optional PAM50/TCGA subtype correlation (for LABELING clusters only)
tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]

if (!is.na(tcga_sub_col)) {
  # TCGA subtype centroids in VST space, restricted to HVGs for fairness
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga)))
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[, lab_tcga == s, drop = FALSE], na.rm = TRUE))
  # cluster centroids (in gene space)
  cl_lvls <- levels(clusters)
  cl_means <- sapply(cl_lvls, function(cn) rowMeans(V[, clusters == cn, drop = FALSE], na.rm = TRUE))
  # restrict to shared genes (use all genes in V to be robust)
  g_common <- intersect(rownames(sub_means), rownames(cl_means))
  corr <- cor(cl_means[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  write.csv(corr, file.path(outdir, "cor_cluster_centroid_vs_TCGA_subtype.csv"))
  safe_pdf(file.path(outdir, "heatmap_cluster_vs_TCGA_subtype_corr.pdf"), {
    pheatmap::pheatmap(corr, main = "Cluster centroids vs TCGA subtype centroids (ρ)")
  })
}

# ---------- (2) MAP DSMZ → TCGA ----------
# Map DSMZ samples to TCGA using Spearman correlation
# Overall BRCA centroid (mean expression across all TCGA samples)
tcga_mean <- rowMeans(V_tcga, na.rm = TRUE)  # Calculate mean expression per gene across all TCGA samples
rho_overall <- sapply(colnames(V_dsmz), function(s) {  # For each DSMZ sample
  suppressWarnings(cor(V_dsmz[, s], tcga_mean, method = "spearman", use = "pairwise.complete.obs"))  # Correlate with TCGA centroid
})

# TCGA subtype centroids (if PAM50/subtype labels available)
tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")  # Possible subtype column names
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]  # Find first available column
rho_by_sub <- NULL  # Initialize subtype correlation matrix
best_tcga_subtype <- NULL  # Initialize best subtype assignments
best_tcga_rho <- NULL  # Initialize best correlation scores

if (!is.na(tcga_sub_col)) {  # If subtype column exists
  lab <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])  # Get subtype labels
  sub_levels <- sort(unique(na.omit(lab)))  # Get unique subtypes
  sub_means <- sapply(sub_levels, function(s) {  # For each subtype
    rowMeans(V_tcga[, lab == s, drop = FALSE], na.rm = TRUE)  # Calculate mean expression
  })
  
  # Compute correlations to each TCGA subtype
  rho_by_sub <- sapply(colnames(V_dsmz), function(s) {  # For each DSMZ sample
    apply(sub_means, 2, function(ref) {  # For each TCGA subtype
      suppressWarnings(cor(V_dsmz[, s], ref, method = "spearman", use = "pairwise.complete.obs"))  # Correlate
    })
  })
  colnames(rho_by_sub) <- colnames(V_dsmz)  # Set column names
  best_tcga_subtype <- apply(rho_by_sub, 2, function(v) names(which.max(v)))  # Find best subtype for each sample
  best_tcga_rho <- apply(rho_by_sub, 2, max)  # Get best correlation score
  write.csv(rho_by_sub, file.path(outdir, "rho_dsmz_to_tcga_subtypes.csv"))  # Save correlation matrix
}

# Create comprehensive mapping table
map_tbl <- tibble(  # Create results table
  sample = colnames(V_dsmz),  # DSMZ sample names
  rho_overall = unname(rho_overall),  # Overall correlation scores
  discovery_cluster = clusters[colnames(V_dsmz)]  # Discovery cluster assignments
)

if (!is.null(best_tcga_subtype)) {  # If subtype mapping was performed
  map_tbl$best_tcga_subtype <- best_tcga_subtype[map_tbl$sample]  # Add best TCGA subtype
  map_tbl$best_tcga_subtype_rho <- best_tcga_rho[map_tbl$sample]  # Add correlation score
}

write.csv(map_tbl, file.path(outdir, "mapping_dsmz_to_tcga.csv"), row.names = FALSE)  # Save mapping table
print(head(map_tbl))  # Print first few rows

# ---------- (3) HEATMAP: marker expression ----------
# Create heatmap of PAM50 marker genes across TCGA and DSMZ samples

# Symbol ↔ Ensembl map from DSMZ (drop version, de-duplicate)
ens2sym <- make_ens2sym(dsmz_raw)

# Ensembl IDs for your marker symbols (case-insensitive), present in V
marker_ens <- rownames(ens2sym)[tolower(ens2sym$symbol) %in% tolower(marker_symbols)]
marker_ens <- intersect(unique(marker_ens), rownames(V))

if (length(marker_ens) == 0L) {
  warning("[marker heatmap] No marker genes found in the joint VST matrix; skipping heatmap.")
} else {
  V_mark <- V[marker_ens, , drop = FALSE]

  # guard 1: if any rownames still collide, make them unique to satisfy pheatmap
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  
  # guard 2: use pretty labels (unique symbols) on the heatmap if you want
  lab <- ens2sym[rownames(V_mark), "symbol_unique", drop = TRUE]
  lab[is.na(lab) | lab == ""] <- rownames(V_mark)
  rownames(V_mark) <- lab

  # Unified annotations
  tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
  tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]
  tcga_sub_vec <- rep(NA_character_, ncol(V_tcga))
  if (!is.na(tcga_sub_col)) {
    tcga_sub_vec <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  }

  ann_tcga <- data.frame(
    dataset = "TCGA",
    discovery_cluster = as.character(clusters[colnames(V_tcga)]),
    TCGA_subtype = tcga_sub_vec,
    row.names = colnames(V_tcga),
    check.names = FALSE
  )
  ann_dsmz <- data.frame(
    dataset = "DSMZ",
    discovery_cluster = as.character(clusters[colnames(V_dsmz)]),
    TCGA_subtype = NA_character_,
    row.names = colnames(V_dsmz),
    check.names = FALSE
  )

  cols_order <- c(colnames(V_tcga), colnames(V_dsmz))
  V_mark <- V_mark[, cols_order, drop = FALSE]
  ann_col <- rbind(ann_tcga, ann_dsmz)[cols_order, , drop = FALSE]
  non_empty_annot <- colSums(!is.na(ann_col)) > 0
  ann_col <- if (any(non_empty_annot)) ann_col[, non_empty_annot, drop = FALSE] else NULL

  if (nrow(V_mark) >= 2 && ncol(V_mark) >= 2) {
    pdf(file.path(outdir, "heatmap_marker_vst.pdf"), width = 12, height = 8)
    pheatmap(
      V_mark, show_rownames = TRUE, show_colnames = FALSE,
      cluster_rows = TRUE, cluster_cols = TRUE,
      annotation_col = ann_col,
      main = "Marker expression (VST) across TCGA-BRCA & DSMZ-BRCA"
    )
    dev.off()
  } else {
    message(sprintf("[marker heatmap] Skipping: matrix dim = %d x %d (need ≥ 2 x 2).",
                    nrow(V_mark), ncol(V_mark)))
  }
}

# ---------- (4) HEATMAP: median z-scores by subtype/cluster ----------
# Create heatmap showing median z-scores for each subtype/cluster
# Z-score normalize each dataset separately for fair comparison

# --- Marker diagnostics & robust Ensembl mapping ---
cat("[MARKER] Requested symbols:", paste(marker_symbols, collapse=", "), "\n")

# 1) build symbol->ensembl from BOTH sources and take the union
map_dsmz <- as_df(dsmz_raw) %>%
  transmute(ensembl = strip_ensver(Ensembl_ID), symbol = as.character(gene_name)) %>%
  distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)

# try to get TCGA gene symbols if available in rowData
map_tcga <- tryCatch({
  rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se))
  symcol <- intersect(tolower(colnames(rd)), c("gene_name","symbol","hgnc_symbol"))
  if (length(symcol)) {
    data.frame(
      ensembl = strip_ensver(rownames(rd)),
      symbol  = as.character(rd[[ colnames(rd)[match(symcol[1], tolower(colnames(rd)))] ]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(symbol) & symbol != "") %>%
      distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)
  } else NULL
}, error = function(e) NULL)

ens2sym_union <- bind_rows(filter(map_dsmz, !is.na(symbol)),
                           filter(map_tcga, !is.null(map_tcga))) %>%
  distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)

# 2) Ensembl markers that exist in V
marker_ens <- ens2sym_union$ensembl[
  tolower(ens2sym_union$symbol) %in% tolower(marker_symbols)
]
marker_ens <- intersect(unique(marker_ens), rownames(V))

cat(sprintf("[MARKER] Found %d of %d PAM50/PAM-like markers in V\n",
            length(marker_ens), length(marker_symbols)))
if (length(marker_ens)) {
  cat("[MARKER] Example matches:", paste(head(marker_ens, 5), collapse=", "), "...\n")
} else {
  stop("[FATAL] No marker genes from the requested list are present in V. Check symbol mapping.")
}

V_tcga_z <- zscore_by_gene(V_tcga)  # Z-score TCGA data
V_dsmz_z <- zscore_by_gene(V_dsmz)  # Z-score DSMZ data

# Helper function to compute median expression per group
aggregate_median <- function(M, groups) {  # M = matrix, groups = group labels
  lev <- sort(unique(groups))  # Get unique group levels
  sapply(lev, function(g) {  # For each group
    rowMedians(M[, groups == g, drop = FALSE], na.rm = TRUE)  # Calculate median per gene
  })
}

# Aggregate TCGA by subtype (if available) or overall
if (!is.na(tcga_sub_col)) {  # If subtype column exists
  tcga_groups <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])  # Get subtype labels
  TC_med <- aggregate_median(V_tcga_z, tcga_groups)  # Calculate median per subtype
} else {  # If no subtype column
  TC_med <- rowMedians(V_tcga_z, na.rm = TRUE)  # Calculate overall median
  TC_med <- matrix(TC_med, ncol = 1,  # Convert to matrix
                   dimnames = list(rownames(V_tcga_z), "TCGA_BRCA"))  # Set dimension names
}

# Aggregate DSMZ by discovery clusters
DSM_med <- aggregate_median(V_dsmz_z, clusters[colnames(V_dsmz)])  # Calculate median per cluster

## ---- safe median-z heatmap ----
common_mark <- intersect(marker_ens, intersect(rownames(TC_med), rownames(DSM_med)))
if (length(common_mark) == 0L) {
  cat("[DEBUG] marker_ens length =", length(marker_ens), "\n")
  cat("[DEBUG] nrow(TC_med) =", nrow(TC_med), " nrow(DSM_med) =", nrow(DSM_med), "\n")
  cat("[DEBUG] Example TC_med genes:", paste(head(rownames(TC_med), 5), collapse=", "), "\n")
  cat("[DEBUG] Example DSM_med genes:", paste(head(rownames(DSM_med), 5), collapse=", "), "\n")
  stop("[FATAL] No overlap between marker_ens and (TC_med ∩ DSM_med). Check ID harmonization (Ensembl, version stripping) and markers.")
} else {
  H <- cbind(TC_med[common_mark, , drop = FALSE], DSM_med[common_mark, , drop = FALSE])
  pdf(file.path(outdir, "heatmap_marker_median_z.pdf"), width = 10, height = 7)
  pheatmap(H, show_rownames = TRUE,
           main = "Median z-scores of marker genes by subtype/cluster")
  dev.off()
}

# ---------- (5) UMAP EMBEDDING ----------
# Create UMAP visualization of joint TCGA+DSMZ data
set.seed(42)  # Set random seed for reproducibility
genes_umap <- top_var_genes(V, n = 3000)  # Select most variable genes for UMAP
Xu <- t(V[genes_umap, , drop = FALSE])    # Transpose to samples x genes (required for UMAP)
emb <- umap(Xu, n_neighbors = 20, min_dist = 0.3, metric = "cosine")  # Run UMAP
emb <- as.data.frame(emb)  # Convert to data frame
colnames(emb) <- c("UMAP1","UMAP2")  # Set column names
emb$sample <- rownames(Xu)  # Add sample names
emb$dataset <- ifelse(emb$sample %in% colnames(V_tcga), "TCGA", "DSMZ")  # Add dataset labels
emb$discovery_cluster <- clusters[emb$sample]  # Add discovery clusters

# Add TCGA subtypes if available
emb$TCGA_subtype <- NA  # Initialize with NA
if (!is.na(tcga_sub_col)) {  # If subtype column exists
  emb$TCGA_subtype[match(colnames(V_tcga), emb$sample)] <-  # For TCGA samples
    as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])  # Add subtype labels
}

# Create UMAP plots
p1 <- ggplot(emb, aes(UMAP1, UMAP2, color = dataset)) +  # Plot by dataset
  geom_point(alpha = 0.8, size = 1.6) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (dataset)")  # Title

p2 <- ggplot(subset(emb, dataset=="TCGA"), aes(UMAP1, UMAP2, color = TCGA_subtype)) +  # TCGA only
  geom_point(alpha = 0.8, size = 1.6) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (TCGA BRCA subtype)")  # Title

p3 <- ggplot(subset(emb, dataset=="DSMZ"), aes(UMAP1, UMAP2, color = discovery_cluster)) +  # DSMZ only
  geom_point(alpha = 0.9, size = 1.9) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (discovery clusters)")  # Title

# Save individual plots
ggsave(file.path(outdir, "umap_dataset.pdf"), p1, width = 7, height = 6)  # Save dataset plot
ggsave(file.path(outdir, "umap_tcga_subtype.pdf"), p2, width = 7, height = 6)  # Save TCGA plot
ggsave(file.path(outdir, "umap_dsmz_cluster.pdf"), p3, width = 7, height = 6)  # Save DSMZ plot

# Save combined plot (try patchwork, fallback to separate PDF)
if (requireNamespace("patchwork", quietly = TRUE)) {  # If patchwork is available
  library(patchwork)  # Load patchwork
  combined <- (p1 | p2 | p3) + plot_layout(guides = "collect")  # Combine plots
  ggsave(file.path(outdir, "umap_combined.pdf"), combined, width = 15, height = 5)  # Save combined
} else {  # Fallback method
  pdf(file.path(outdir, "umap_combined.pdf"), width = 15, height = 5)  # Open PDF
  print(p1)  # Print plot 1
  print(p2)  # Print plot 2
  print(p3)  # Print plot 3
  dev.off()  # Close PDF
}

# Save embeddings for further analysis
write.csv(emb, file.path(outdir, "umap_embeddings.csv"), row.names = FALSE)  # Save coordinates

cat("\n[OK] Done. Results in:", outdir, "\n")  # Print completion message
# nolint end
