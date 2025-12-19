#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ConsensusClusterPlus)
  library(matrixStats)
})

info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
stop_with <- function(...) stop(paste0(...), call. = FALSE)

option_list <- list(
  make_option("--config",      type = "character", default = "config/config.yaml"),
  make_option("--direction",   type = "character", default = NULL),
  make_option("--kind",        type = "character", default = NULL),
  make_option("--mode",        type = "character", default = NULL),  # expr | pca
  make_option("--alg",         type = "character", default = NULL),  # km | hc
  make_option("--outdir",      type = "character", default = NULL),
  make_option("--cluster_rds", type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$direction))   stop_with("Missing --direction")
if (is.null(opt$kind))        stop_with("Missing --kind")
if (is.null(opt$mode))        stop_with("Missing --mode (expr|pca)")
if (is.null(opt$alg))         stop_with("Missing --alg (km|hc)")
if (is.null(opt$outdir))      stop_with("Missing --outdir")
if (is.null(opt$cluster_rds)) stop_with("Missing --cluster_rds")

if (!opt$mode %in% c("expr","pca")) stop_with("--mode must be expr or pca")
if (!opt$alg  %in% c("km","hc"))    stop_with("--alg must be km or hc")

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- yaml::read_yaml(opt$config)

if (is.null(cfg$methods[[opt$direction]])) {
  stop_with("Unknown direction in config methods: ", opt$direction)
}

feature  <- cfg$methods[[opt$direction]]$feature   # PAM50 | HVG
distance <- cfg$methods[[opt$direction]]$distance  # euclidean | correlation

k_grid <- cfg$clustering$k_grid
if (is.null(k_grid) || length(k_grid) < 2) stop_with("config: clustering.k_grid must have >=2 values")
k_grid <- sort(as.integer(unlist(k_grid)))
maxK   <- max(k_grid)

seed <- 42
if (!is.null(cfg$agnostic_clustering$seed)) seed <- as.integer(cfg$agnostic_clustering$seed)

n_pcs <- 20
if (!is.null(cfg$agnostic_clustering$n_pcs)) n_pcs <- as.integer(cfg$agnostic_clustering$n_pcs)

# ─────────────────────────────────────────────────────────────
# IMPORTANT:
# We ALWAYS run CCP with euclidean distance so Ward.D2 is valid.
# For "correlation" directions, we transform columns (samples)
# so euclidean becomes correlation/cosine-like geometry.
# ─────────────────────────────────────────────────────────────
ccp_distance <- "euclidean"  # ALWAYS
innerLinkage <- "ward.D2"    # ALWAYS (prevents CCP hclust(method=NULL) even for km)
finalLinkage <- "ward.D2"    # ALWAYS

info("direction=%s | feature=%s | distance=%s | mode=%s | alg=%s | kind=%s",
     opt$direction, feature, distance, opt$mode, opt$alg, opt$kind)

build_expr_mat <- function(feature, cfg, kind) {
  # Determine scope from kind
  is_cell_tumour <- grepl("_cell_tumour", kind)
  is_cell_only   <- grepl("_cell$", kind) && !is_cell_tumour
  is_tumour_only <- grepl("_tumour$", kind) && !is_cell_tumour

  if (feature == "PAM50") {
    tcga_path <- cfg$paths$tcga_brca_pam50_expr
    dsmz_path <- cfg$paths$dsmz_bcc_pam50_expr
    if (is.null(tcga_path) || is.null(dsmz_path)) {
      stop_with("Missing pam50 paths: paths.tcga_brca_pam50_expr or paths.dsmz_bcc_pam50_expr")
    }
    tcga_mat <- readRDS(tcga_path)  # genes × samples
    dsmz_mat <- readRDS(dsmz_path)

    common_genes <- intersect(rownames(tcga_mat), rownames(dsmz_mat))
    if (length(common_genes) < 20) stop_with("Too few common genes between TCGA and DSMZ PAM50 matrices.")

    if (is_cell_tumour) {
      expr_mat <- cbind(tcga_mat[common_genes, , drop = FALSE],
                        dsmz_mat[common_genes, , drop = FALSE])
      dataset <- ifelse(colnames(expr_mat) %in% colnames(dsmz_mat), "DSMZ", "TCGA")
      names(dataset) <- colnames(expr_mat)
    } else if (is_cell_only) {
      expr_mat <- dsmz_mat[common_genes, , drop = FALSE]
      dataset <- setNames(rep("DSMZ", ncol(expr_mat)), colnames(expr_mat))
    } else if (is_tumour_only) {
      expr_mat <- tcga_mat[common_genes, , drop = FALSE]
      dataset <- setNames(rep("TCGA", ncol(expr_mat)), colnames(expr_mat))
    } else {
      stop_with("Cannot determine scope from kind: ", kind)
    }

    return(list(expr_mat = expr_mat, dataset = dataset))
  }

  if (feature == "HVG") {
    cell_path   <- cfg$paths$cell_vst_rds
    tumour_path <- cfg$paths$tumour_vst_rds
    hvg_list    <- cfg$features$hvg_final_gene_list

    if (is.null(cell_path) || is.null(tumour_path)) stop_with("Missing: paths.cell_vst_rds and/or paths.tumour_vst_rds")
    if (is.null(hvg_list) || !file.exists(hvg_list)) stop_with("Missing HVG list: features.hvg_final_gene_list")

    cell_mat   <- readRDS(cell_path)    # genes × samples (DSMZ)
    tumour_mat <- readRDS(tumour_path)  # genes × samples (TCGA)

    genes <- readLines(hvg_list); genes <- genes[nzchar(genes)]

    if (is_cell_tumour) {
      keep <- intersect(genes, intersect(rownames(cell_mat), rownames(tumour_mat)))
      if (length(keep) < 50) stop_with("Too few HVG genes found in both cell and tumour matrices.")
      expr_mat <- cbind(tumour_mat[keep, , drop = FALSE],
                        cell_mat[keep, , drop = FALSE])
      dataset <- c(setNames(rep("TCGA", ncol(tumour_mat[keep, , drop = FALSE])), colnames(tumour_mat)),
                   setNames(rep("DSMZ", ncol(cell_mat[keep, , drop = FALSE])), colnames(cell_mat)))
      dataset <- dataset[colnames(expr_mat)]
      names(dataset) <- colnames(expr_mat)
    } else if (is_cell_only) {
      keep <- intersect(genes, rownames(cell_mat))
      if (length(keep) < 50) stop_with("Too few HVG genes found in cell matrix.")
      expr_mat <- cell_mat[keep, , drop = FALSE]
      dataset <- setNames(rep("DSMZ", ncol(expr_mat)), colnames(expr_mat))
    } else if (is_tumour_only) {
      keep <- intersect(genes, rownames(tumour_mat))
      if (length(keep) < 50) stop_with("Too few HVG genes found in tumour matrix.")
      expr_mat <- tumour_mat[keep, , drop = FALSE]
      dataset <- setNames(rep("TCGA", ncol(expr_mat)), colnames(expr_mat))
    } else {
      stop_with("Cannot determine scope from kind: ", kind)
    }

    return(list(expr_mat = expr_mat, dataset = dataset))
  }

  stop_with("Unsupported feature: ", feature)
}

# Transform columns so Euclidean distance is correlation-like
# (center per sample, then L2-normalise per sample)
corr_geometry_transform <- function(X) {
  # X is features × samples
  X <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  norms <- sqrt(colSums(X^2, na.rm = TRUE))
  norms[!is.finite(norms) | norms == 0] <- 1
  X <- sweep(X, 2, norms, "/")
  X
}

built    <- build_expr_mat(feature, cfg, opt$kind)
expr_mat <- built$expr_mat
dataset  <- built$dataset

# Drop constant genes
v <- rowVars(expr_mat, na.rm = TRUE)
expr_mat <- expr_mat[v > 0 & is.finite(v), , drop = FALSE]
if (anyDuplicated(colnames(expr_mat))) stop_with("Duplicate sample IDs detected in integrated matrix.")

# CCP expects: rows = features, cols = items (samples)
if (opt$mode == "expr") {
  data_for_ccp <- expr_mat  # genes × samples
  info("CCP input = expr (genes × samples): %d × %d", nrow(data_for_ccp), ncol(data_for_ccp))
} else {
  samples_mat <- t(expr_mat)  # samples × genes
  pca <- prcomp(samples_mat, center = TRUE, scale. = FALSE)
  pcs_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]  # samples × PCs
  data_for_ccp <- t(pcs_mat)  # PCs × samples
  info("CCP input = pca (PCs × samples): %d × %d", nrow(data_for_ccp), ncol(data_for_ccp))
}

# If direction is "correlation", transform so Euclidean == correlation-like geometry
if (distance == "correlation") {
  info("Applying correlation-geometry transform (center + unit norm) so Ward.D2 uses euclidean safely.")
  data_for_ccp <- corr_geometry_transform(data_for_ccp)
}

# Isolate CCP outputs
ccp_wd <- file.path(opt$outdir, "ccp")
dir.create(ccp_wd, showWarnings = FALSE, recursive = TRUE)

old_wd <- getwd()
setwd(ccp_wd)
on.exit(setwd(old_wd), add = TRUE)

# Build args (NOTE: always pass linkage to avoid CCP hclust(method=NULL) path)
ccp_args <- list(
  d            = data_for_ccp,
  maxK         = maxK,
  reps         = 1000,
  pItem        = 0.8,
  pFeature     = 1.0,
  clusterAlg   = opt$alg,
  distance     = ccp_distance,   # ALWAYS euclidean
  innerLinkage = innerLinkage,   # ALWAYS ward.D2
  finalLinkage = finalLinkage,   # ALWAYS ward.D2
  seed         = seed,
  title        = paste0("CCP_", opt$direction, "_", opt$kind),
  plot         = "png",
  writeTable   = TRUE,
  verbose      = FALSE
)

ccp_res <- do.call(ConsensusClusterPlus, ccp_args)

auc_at_k <- function(k) {
  cm <- ccp_res[[k]]$consensusMatrix
  vals <- cm[upper.tri(cm)]
  ec <- ecdf(vals)
  x <- seq(0, 1, length.out = 500)
  sum(diff(x) * (ec(head(x, -1)) + ec(tail(x, -1))) / 2)
}

auc_values <- vapply(k_grid, auc_at_k, numeric(1))
delta_auc  <- diff(auc_values)
best_k     <- k_grid[which.max(delta_auc) + 1]

info("k_grid: %s", paste(k_grid, collapse = ","))
info("AUC:    %s", paste(round(auc_values, 4), collapse = " "))
info("ΔAUC:   %s", paste(round(delta_auc, 4), collapse = " "))
info("best_k: %d", best_k)

final_clusters <- ccp_res[[best_k]]$consensusClass
final_clusters <- setNames(paste0("C", final_clusters), names(final_clusters))

clusters_csv <- file.path(opt$outdir, sprintf("tcga_dsmz_ccp_clusters_K%d.csv", best_k))
write_csv(
  tibble(sample_id = names(final_clusters),
         cluster   = unname(final_clusters),
         dataset   = dataset[names(final_clusters)]),
  clusters_csv
)

saveRDS(
  list(
    direction      = opt$direction,
    kind           = opt$kind,
    feature        = feature,
    distance       = distance,
    mode           = opt$mode,
    alg            = opt$alg,
    ccp_distance   = ccp_distance,
    innerLinkage   = innerLinkage,
    finalLinkage   = finalLinkage,
    k_grid         = k_grid,
    auc_values     = setNames(auc_values, k_grid),
    delta_auc      = delta_auc,
    best_k         = best_k,
    clusters       = final_clusters,
    dataset        = dataset,
    ccp_results    = ccp_res,
    timestamp      = Sys.time()
  ),
  opt$cluster_rds
)

info("Wrote: %s", clusters_csv)
info("Wrote: %s", opt$cluster_rds)
info("Done.")
