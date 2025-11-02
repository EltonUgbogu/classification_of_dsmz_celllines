# Helpers: low-level utility functions
# Place in R/helpers.R

# Library imports for helper code (not necessary to call library() here; main script should load them)
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

# Simple logging helpers
log_info  <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
log_warn  <- function(...) cat("[WARN] ", sprintf(...), "\n", sep = "")
log_error <- function(...) cat("[ERROR]", sprintf(...), "\n", sep = "")

#' @brief Calculates minPts for clustering based on dataset size.
#' @param n Integer, total samples.
#' @param frac Numeric, fraction of samples for minPts (default: 0.02).
#' @param min_floor Integer, minimum minPts (default: 5).
#' @param max_cap Integer, maximum minPts (default: 50).
#' @return Integer, calculated minPts.
# Determines the minPts parameter for HDBSCAN based on sample size
choose_minPts <- function(n_obs) {
  # Apply heuristic: minPts is 5 times the log10 of sample size (min 10)
  val <- as.integer(round(log10(max(10, n_obs)) * 5))
  # Clamp the value to ensure minPts is between 5 and 50 (inclusive)
  min(50L, max(5L, val))
}


# Determines minCluster size for dynamic tree cut with flexible parameters
choose_minCluster <- function(n, frac = 0.03, min_floor = 8, max_cap = 150) {
  # Calculate minCluster as fraction of sample size, respecting min_floor
  p <- max(min_floor, round(frac * n))
  # Return minCluster, capped at max_cap
  as.integer(min(p, max_cap))
}


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
  ens2sym <- as.data.frame(ens2sym, stringsAsFactors = FALSE)
  rownames(ens2sym) <- ens2sym$ensembl
  ens2sym
}

# mat - matrix of data
# n - number of top variable genes to return
top_var_genes <- function(mat, n) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}

make_pcs <- function(M, n_hvg = 3000, max_pc = 30) {
  v <- matrixStats::rowVars(M); keep <- which(is.finite(v) & v > 0)
  if (length(keep) < 2) stop("[FATAL] Too few variable genes for PCA.")
  M <- M[keep, , drop = FALSE]
  hvgs_all <- names(sort(matrixStats::rowIQRs(M), decreasing = TRUE))
  hvgs <- hvgs_all[seq_len(min(n_hvg, length(hvgs_all)))]
  X <- t(scale(M[hvgs, , drop = FALSE]))
  X <- X[, colSums(!is.na(t(X))) > 0, drop = FALSE]
  if (nrow(X) < 2 || ncol(X) < 1) stop("[FATAL] PCA input has insufficient dimension.")
  pr <- prcomp(X, center = FALSE, scale. = FALSE)
  if (is.null(pr$x) || ncol(pr$x) < 1) stop("[FATAL] PCA returned 0 components.")
  kPC <- as.integer(min(max_pc, ncol(pr$x)))
  list(PC = pr$x[, 1:kPC, drop = FALSE],
       kPC = kPC,
       pr = pr,
       hvgs = hvgs)
}

# alias to allow wrappers in other files without naming collisions
helpers_make_pcs <- function(M, n_hvg = 3000, max_pc = 30) make_pcs(M, n_hvg = n_hvg, max_pc = max_pc)

# convenience: return only the PC matrix (used by I/O and clustering code)
make_pcs_matrix <- function(V_adj, n_hvg = 3000, max_pc = 30, center = TRUE, scale. = TRUE) {
  if (!is.matrix(V_adj) && !is.data.frame(V_adj)) {
    stop("V_adj must be a matrix or data.frame")
  }
  mat <- as.matrix(V_adj)
  if (ncol(mat) < 2 || nrow(mat) < 2) stop("V_adj must have at least 2 rows and 2 columns")
  pcs_obj <- helpers_make_pcs(mat, n_hvg = n_hvg, max_pc = max_pc)
  pcs_obj$PC
}

zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

safe_pdf <- function(path, expr) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path)
  on.exit({
    if (dev.cur() != 1) dev.off()
  }, add = TRUE)
  force(expr)
}

build_dsmz_matrix <- function(dsmz_raw) {
  stopifnot(all(c("Ensembl_ID","gene_name") %in% names(dsmz_raw)))
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot)
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(v) as.numeric(as.character(v)))
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(M) <- strip_ensver(dsmz_raw$Ensembl_ID)
  if (any(duplicated(rownames(M)))) {
    M <- rowsum(M, rownames(M), reorder = TRUE)
  }
  M
}

logCPM <- function(x, prior.count = 1, lib.size = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x)
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size)
}

# Map symbols -> ensembl lookup (expects ens2sym_union data.frame with ensembl,symbol)
make_sym2ens <- function(ens2sym_union) {
  setNames(ens2sym_union$ensembl, toupper(ens2sym_union$symbol))
}

# relabel helper for HDBSCAN outputs
relabel_hdb_numeric <- function(x, keep_unassigned = TRUE) {
  x <- as.character(x)
  unass <- x == "Unassigned"
  labs  <- sort(unique(x[!unass]))
  map   <- setNames(as.character(seq_along(labs)), labs)
  y     <- ifelse(unass & keep_unassigned, NA, map[x])
  factor(y, levels = as.character(seq_len(min(5, length(labs)))))
}

aggregate_median <- function(M, groups) {
  grp <- droplevels(factor(groups))
  lev <- levels(grp)
  out <- sapply(lev, function(g) {
    rowMedians(M[, grp == g, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  rownames(out) <- rownames(M)
  colnames(out) <- lev
  out
}

# centroid helpers
get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
get_k_nearest <- function(X, centroid, k = 20) {
  d <- sqrt(rowSums((X - matrix(centroid, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
  ord <- order(d, decreasing = FALSE)
  keep <- ord[seq_len(min(k, length(ord)))]
  data.frame(sample = rownames(X)[keep], dist = d[keep], stringsAsFactors = FALSE)
}

make_umap <- function(PC, seed = 42, n_neighbors = 20, min_dist = 0.3,
                     metric = "cosine") {
  set.seed(seed)
  emb <- uwot::umap(PC, n_neighbors = n_neighbors, min_dist = min_dist,
                    metric = metric, verbose = FALSE)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample   <- rownames(PC)
  emb
}


# Fixed-k hierarchical clustering
run_fixed_hc <- function(PC, k = 4, dist_method = "euclidean",
                         linkage = "average", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(nrow(PC) >= k, "Too few samples for the requested k.")

  log_info("Calculating distance matrix (%s) on %d PCs for %d samples",
           dist_method, ncol(PC), nrow(PC))
  d <- stats::dist(PC, method = dist_method)

  log_info("Building dendrogram with '%s' linkage", linkage)
  hc <- stats::hclust(d, method = linkage)

  log_info("Cutting dendrogram into k = %d clusters", k)
  clusters <- stats::cutree(hc, k = k)

  list(hc = hc, dist = d, clusters = clusters)
}

# Dynamic Tree Cut wrapper
run_dynamic_cut <- function(PC, min_cluster_size = NULL,
                            deep_split = 2, pam_stage = TRUE,
                            dist_method = "euclidean",
                            linkage = "average") {
  d  <- stats::dist(PC, method = dist_method)
  hc <- stats::hclust(d, method = linkage)
  distM <- as.matrix(d)

  if (is.null(min_cluster_size)) {
    min_cluster_size <- choose_minCluster(nrow(PC))
  }
  log_info("Running cutreeHybrid (deepSplit=%d, pamStage=%s, minClusterSize=%d)",
           deep_split, pam_stage, min_cluster_size)

  dyn <- dynamicTreeCut::cutreeHybrid(dendro = hc, distM = distM,
                                      deepSplit = deep_split,
                                      pamStage = pam_stage,
                                      minClusterSize = min_cluster_size)

  clusters <- factor(dyn$labels)
  list(hc = hc, dist = d, clusters = clusters, raw = dyn)
}

# Legacy aliases
run_dynamic_tree_cut <- function(PC, deepSplit = 2, pamStage = TRUE) {
  run_dynamic_cut(PC, deep_split = deepSplit, pam_stage = pamStage)
}

# Centroid calculation over HVGs
calc_centroids <- function(expr_mat, hvgs, group_vec) {
  stopifnot(length(group_vec) == ncol(expr_mat))
  groups <- levels(factor(group_vec))
  sapply(groups, function(g) {
    rowMeans(expr_mat[hvgs, group_vec == g, drop = FALSE], na.rm = TRUE)
  })
}

# ----------------------------------------------------------------------
# Supervised labeling of clusters using TCGA-BRCA PAM50 subtype centroids
# Assumes BRCA-only context; V_adj may contain TCGA + DSMZ samples
# metadata_ann: CSV file path or data.frame with sample_barcode & PAM50_Subtype columns
# ----------------------------------------------------------------------
label_by_tcga <- function(clusters,          # cluster assignments from clustering algorithm (for all samples)
                          V_adj,             # genes x samples of dataset clustered
                          features,          # feature genes (HVGs or PAM50)  
                          metadata_ann,      # metadata annotation CSV 
                          V_tcga,            # genes x TCGA-BRCA samples
                          subtype = "PAM50_Subtype") {  # column name in CSV metadata
  log_info("=== Starting supervised label_by_tcga ===")
  stopifnot(length(clusters) == ncol(V_adj), !is.null(colnames(V_tcga)))
  
  # Load subtype annotations from metadata annotation CSV
  if (is.character(metadata_ann)) {
    # CSV path
    tcga_df <- read.csv(metadata_ann, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(metadata_ann)) {
    # Already a data.frame
    tcga_df <- metadata_ann
  } else {
    stop("metadata_ann must be a CSV file path (character) or data.frame")
  }
  # check if the metadata_ann has the required columns
  stopifnot("sample_barcode" %in% colnames(tcga_df), subtype %in% colnames(tcga_df))
  # match the sample_barcode in the metadata_ann to the column names of V_tcga (TCGA-BRCA expression matrix)
  matched <- tcga_df[tcga_df$sample_barcode %in% colnames(V_tcga), ]
  if (nrow(matched) == 0) {
    stop("No matches between TCGA-BRCA samples and metadata_ann sample_barcode column")
  }
  # set the names of the labs to the PAM50 Subtype and the V_tcga (TCGA-BRCA expression matrix) to the matched TCGA-BRCA sample_barcode
  labs <- setNames(as.character(matched[[subtype]]), matched$sample_barcode)
  V_tcga <- V_tcga[, matched$sample_barcode, drop = FALSE] # subset the V_tcga (TCGA-BRCA expression matrix) to the matched sample_barcode
  
  # Drop samples with missing PAM50 Subtype labels in the metadata_ann
  keep <- !is.na(labs)
  if (sum(!keep) > 0) {
    log_warn("Dropping %d TCGA-BRCA samples with missing PAM50 Subtype labels in the metadata_ann.", sum(!keep))
    V_tcga <- V_tcga[, keep, drop = FALSE]
    labs <- labs[keep]
  }
  tcga_lvls <- sort(unique(labs)) # get the unique PAM50 Subtypes in the metadata_ann
  log_info("PAM50 Subtypes found in TCGA-BRCA metadata_ann: %s", paste(tcga_lvls, collapse=", "))

  # Compute centroids on shared features
  features_use <- intersect(features, intersect(rownames(V_tcga), rownames(V_adj))) # get the shared features between the V_tcga (TCGA-BRCA expression matrix) and the V_adj (dataset clustered expression matrix)
  if (length(features_use) < 40L) {
    stop(sprintf("Check correctness of selected TCGA-BRCA features: %d ", length(features_use)))
  }
  log_info("Computing centroids on %d features…", length(features_use))
  
  # Pass labels to calc_centroids (labs is named vector from metadata_ann)
  tcga_means <- calc_centroids(V_tcga, features_use, unname(labs)) # calculate the centroids of the TCGA-BRCA expression matrix on the shared features
  disc_means <- calc_centroids(V_adj, features_use, clusters) # calculate the centroids of the dataset clustered expression matrix on the shared features

  # Correlate and map clusters to PAM50 subtypes
  common_g <- intersect(rownames(tcga_means), rownames(disc_means)) # get the shared features between the tcga_means and the disc_means
  tcga_means <- tcga_means[common_g, , drop = FALSE] # subset the tcga_means to the shared features
  disc_means <- disc_means[common_g, , drop = FALSE] # subset the disc_means to the shared features
  # correlate the disc_means and tcga_means
  corr_mat <- stats::cor(disc_means, tcga_means, method = "spearman", use = "pairwise.complete.obs")
  log_info("Correlation matrix: %d clusters x %d PAM50 Subtypes in TCGA-BRCA metadata_ann", nrow(corr_mat), ncol(corr_mat))
  print(round(corr_mat, 3))

  # Map each cluster to best-matching PAM50 subtype
  best_map <- setNames(tcga_lvls[apply(corr_mat, 1, which.max)], rownames(corr_mat)) # map the clusters to the best-matching PAM50 subtype
  map_df <- data.frame(Cluster = names(best_map), PAM50_Subtype = unname(best_map), stringsAsFactors = FALSE) # create a data frame with the clusters and the best-matching PAM50 subtype
  log_info("=== Cluster → PAM50 Mapping ===\n%s", paste(capture.output(print(map_df, row.names = FALSE)), collapse = "\n"))

  # Apply mapping to the clusters
  cl_chr <- as.character(clusters) # convert the clusters to a character vector
  if (!all(unique(cl_chr) %in% names(best_map))) {
    log_warn("Missing cluster mapping in the TCGA-BRCA metadata_ann")
  }
  labeled <- as.factor(factor(unname(best_map[cl_chr]), levels = sort(unique(unname(best_map)))) # label the clusters with the best-matching PAM50 subtype
  log_info("Labeling complete. Levels: %s", paste(levels(labeled), collapse=", ")) # log the levels of the labeled clusters
  return(labeled) # return the labeled clusters
}
plot_dendrogram_dynamic <- function(hc, clusters, out_pdf) {
  cols <- WGCNA::labels2colors(as.numeric(factor(clusters)))
  safe_pdf(out_pdf, {
    WGCNA::plotDendroAndColors(hc, cols,
                               groupLabels = "Dynamic cut",
                               main = "HC + Dynamic Tree Cut (HVG PCs)",
                               dendroLabels = FALSE, cex.dendroLabels = 0.3)
  })
}

plot_silhouette <- function(sil_obj, out_pdf) {
  safe_pdf(out_pdf, {
    plot(sil_obj, main = sprintf("Silhouette – mean = %.3f", mean(sil_obj[, "sil_width"])))
  })
}






write_dimension_report <- function(outdir, dims_list) {
  dim_report <- file.path(outdir, "dimension_report.txt")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  con <- file(dim_report, open = "wt")
  on.exit(close(con), add = TRUE)
  lapply(names(dims_list), function(nm) write_dims_line(con, nm, dims_list[[nm]]))
}

safe_rds_save <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}

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

## form clustering.R
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

safe_rds_save <- function(object, file) {
  ensure_dir(dirname(file))
  saveRDS(object, file = file)
}

run_pca_and_save <- function(V_adj, outdir, n_hvg = 3000, max_pc = 30) {
  ensure_dir(outdir)
  pcs <- make_pcs_matrix(V_adj, n_hvg = n_hvg, max_pc = max_pc)
  safe_rds_save(pcs, file.path(outdir, "pca_objects.rds"))
  pcs
}