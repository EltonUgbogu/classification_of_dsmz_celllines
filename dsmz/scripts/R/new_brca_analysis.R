.libPaths("/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")
Sys.setenv(R_LIBS_USER="/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")

suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
  library(uwot)
  library(cluster)
  library(fpc)
  library(clusterCrit)
  library(limma)
  library(sva)
  library(dbscan)
  library(PMCMRplus)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dynamicTreeCut)  # for cutreeDynamic / cutreeHybrid
  library(WGCNA)           # for plotDendroAndColors
  library(paran)           # for Parallel Analysis for PCs
})

message(sprintf("[INFO] Using %d threads (SLURM_CPUS_PER_TASK)", getOption("mc.cores", 1L)))
## ---------- End Prelude ----------


# ---------- CONFIG ----------
tcga_se_rds   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"
dsmz_rds      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
dsmz_meta_csv <- "/home/chu25/data/dsmz/DSMZ_metadata.csv"
outdir <- "/home/chu25/dsmz/results/brca_vst_dual3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

tcga_project_filter <- c("TCGA-BRCA")
dsmz_organ_filter   <- c("Breast")
marker_symbols <- c("ESR1","PGR","ERBB2","MKI67","KRT5","KRT14","KRT17","FOXC1","GRB7",
                    "BCL2","EGFR","CCNB1","BIRC5","MYBL2","KIF2C","KRT8","KRT18","GATA3","XBP1")
pam50_symbols <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC6","CDC20",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT5","KRT14","KRT17","KRT23","MDM2","MELK",
  "MIA","MKI67","MLPH","MMP11","MYBL2","NAT1","ORC6L","PGR","PHGDH","PTTG1",
  "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T","WHSC1L1","KRT7","KRT15"
)
config <- list(batch_adjust_method = "combat_seq")

# Purity adjustment control flag
drop_purity_genes_permanently <- TRUE  # set TRUE to actually remove rows downstream

# ---------- HELPERS ----------
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

# ------- Adaptive minPts for HDBSCAN -------
# Choose minPts proportional to subset size with sensible caps
choose_minPts <- function(n, frac = 0.02, min_floor = 5, max_cap = 50) {
  p <- max(min_floor, round(frac * n))
  p <- min(p, max_cap, n - 1L)
  return(as.integer(p))
}

# Choose min cluster size as a fraction of n, with sensible caps
choose_minCluster <- function(n, frac = 0.03, min_floor = 8, max_cap = 150) {
  p <- max(min_floor, round(frac * n))
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
top_var_genes <- function(mat, n = 3000) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}
make_pcs <- function(M, n_hvg = 3000, max_pc = 30) {
  # keep genes with non-zero variance
  v <- matrixStats::rowVars(M); keep <- which(is.finite(v) & v > 0)
  if (length(keep) < 2) stop("[FATAL] Too few variable genes for PCA.")

  M <- M[keep, , drop = FALSE]

  # HVGs by IQR
  hvgs_all <- names(sort(matrixStats::rowIQRs(M), decreasing = TRUE))
  hvgs <- hvgs_all[seq_len(min(n_hvg, length(hvgs_all)))]

  X <- t(scale(M[hvgs, , drop = FALSE]))
  # drop columns that become all-NA after scaling (rare but possible)
  X <- X[, colSums(!is.na(t(X))) > 0, drop = FALSE]

  # guard: need ≥2 samples and ≥1 variable column
  if (nrow(X) < 2 || ncol(X) < 1) stop("[FATAL] PCA input has insufficient dimension.")

  pr <- prcomp(X, center = FALSE, scale. = FALSE)
  if (is.null(pr$x) || ncol(pr$x) < 1) stop("[FATAL] PCA returned 0 components. Check input variance/NA.")

  kPC <- as.integer(min(max_pc, ncol(pr$x)))
  list(PC = pr$x[, 1:kPC, drop = FALSE],
       kPC = kPC,
       pr = pr,
       hvgs = hvgs)
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
    dup_genes <- rownames(M)[duplicated(rownames(M))]
    cat(sprintf("[INFO] Collapsing %d duplicate DSMZ genes: %s\n", 
                length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
    M <- rowsum(M, rownames(M), reorder = TRUE)
  }
  M
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
logCPM <- function(x, prior.count = 1, lib.size = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x)
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size)
}
load_purity_data <- function(tcga_se, purity_file = NULL) {
  p <- NULL
  if ("purity" %in% colnames(colData(tcga_se))) {
    p <- as.numeric(colData(tcga_se)$purity)
    names(p) <- colnames(assay(tcga_se))
    cat("[INFO] Using purity from SummarizedExperiment\n")
  } else if (!is.null(purity_file) && file.exists(purity_file)) {
    tab <- read.delim(purity_file, stringsAsFactors = FALSE)
    stopifnot(all(c("sample","purity") %in% colnames(tab)))
    p <- setNames(tab$purity, tab$sample)
    cat("[INFO] Using external purity file\n")
  } else {
    cat("[INFO] No purity data available\n")
  }
  p
}
purity_adjust_log <- function(tcga_log, purity) {
  if (is.null(purity)) return(list(V = tcga_log, dropped = character(0)))

  keep_samp <- intersect(colnames(tcga_log), names(purity)[!is.na(purity)])
  if (length(keep_samp) < 3) return(list(V = tcga_log, dropped = character(0)))

  M <- tcga_log[, keep_samp, drop = FALSE]
  pu <- purity[keep_samp]

  ct <- apply(M, 1, function(v) {
    if (var(v) == 0 || var(pu) == 0) return(list(p.value = 1, estimate = 0))
    suppressWarnings(cor.test(as.numeric(v), pu, method = "spearman"))
  })
  pv <- vapply(ct, `[[`, numeric(1), "p.value"); pv[is.na(pv)] <- 1
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))
  padj <- p.adjust(pv, "BH")

  drop <- names(which(padj < 0.01 & rh < -0.4))
  message(sprintf("[INFO] Purity-associated genes for model fitting: %d", length(drop)))

  if (length(drop)) M <- M[setdiff(rownames(M), drop), , drop = FALSE]

  # fit + adjust only on remaining genes
  infilt <- 1 - pu
  design <- model.matrix(~ infilt)
  fit <- limma::lmFit(M, design)
  beta <- fit$coefficients[, "infilt", drop = FALSE]
  Madj <- as.matrix(M) - beta %*% t(infilt)

  out <- tcga_log
  common_g <- intersect(rownames(Madj), rownames(out))
  out[common_g, colnames(Madj)] <- Madj[common_g, ]

  list(V = out, dropped = drop)
}
subset_to_pam50 <- function(mat, ens2sym, pam_symbols, collapse_fun = c("sum","mean","max")) {
  cf <- match.arg(collapse_fun)
  ensembl_ids <- strip_ensver(rownames(mat))

  # which rows are mapped?
  idx_map <- which(ensembl_ids %in% rownames(ens2sym))
  if (!length(idx_map)) stop("No Ensembl IDs mapped to symbols. Check inputs.")

  # symbols corresponding to mapped rows
  sym <- ens2sym[ensembl_ids[idx_map], "symbol"]
  hit <- tolower(sym) %in% tolower(pam_symbols)

  keep_idx <- idx_map[hit]
  if (!length(keep_idx)) stop("No PAM50 symbols matched after mapping.")

  sub <- mat[keep_idx, , drop = FALSE]
  sym_keep <- ens2sym[strip_ensver(rownames(sub)), "symbol"]

  # collapse multiple Ensembl rows per symbol
  split_idx <- split(seq_len(nrow(sub)), tolower(sym_keep))
  collapsed <- lapply(split_idx, function(ix) {
    if (cf == "sum")      colSums(sub[ix, , drop = FALSE], na.rm = TRUE)
    else if (cf == "mean")  colMeans(sub[ix, , drop = FALSE], na.rm = TRUE)
    else                    apply(sub[ix, , drop = FALSE], 2, max, na.rm = TRUE)
  })
  M <- do.call(rbind, collapsed)
  rownames(M) <- toupper(names(split_idx))
  M
}

# ---------- LOAD + FILTER ----------
cat("[INFO] Starting pipeline...\n")

tcga_se <- readRDS(tcga_se_rds)
stopifnot("project_id" %in% colnames(colData(tcga_se)))
dsmz_raw  <- readRDS(dsmz_rds)
dsmz_meta <- read.csv(dsmz_meta_csv, check.names = FALSE)

tcga_counts <- assay(tcga_se)
rownames(tcga_counts) <- strip_ensver(rownames(tcga_counts))
if (any(duplicated(rownames(tcga_counts)))) {
  dup_genes <- rownames(tcga_counts)[duplicated(rownames(tcga_counts))]
  cat(sprintf("[INFO] Collapsing %d duplicate TCGA genes: %s\n", 
              length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
  tcga_counts <- rowsum(tcga_counts, rownames(tcga_counts), reorder = TRUE)
}
storage.mode(tcga_counts) <- "double"

dsmz_counts <- build_dsmz_matrix(dsmz_raw)
storage.mode(dsmz_counts) <- "double"

tcga_counts_raw <- tcga_counts
dsmz_counts_raw <- dsmz_counts
dsmz_meta_raw <- dsmz_meta

log_dims("TCGA (raw)", tcga_counts_raw)
log_dims("DSMZ (raw)", dsmz_counts_raw)

dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)
matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))
dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched)
dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop = FALSE]
stopifnot("organ" %in% colnames(dsmz_meta))

if (!is.null(dsmz_organ_filter)) {
  keep_ids <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
  dsmz_meta <- dsmz_meta[dsmz_meta$sample_id %in% keep_ids, , drop = FALSE]
  dsmz_counts <- dsmz_counts[, keep_ids, drop = FALSE]
}

log_dims("DSMZ (after organ subset)", dsmz_counts)

if (!is.null(dsmz_organ_filter) && ncol(dsmz_counts) == 0L) {
  stop(sprintf("[FATAL] DSMZ organ filter %s left 0 samples.", paste(dsmz_organ_filter, collapse=", ")))
}

if (!is.null(tcga_project_filter)) {
  keep_tcga <- colnames(tcga_counts) %in% colnames(assay(tcga_se)) &
               as.character(colData(tcga_se)[colnames(tcga_counts), "project_id"]) %in% tcga_project_filter
  tcga_counts <- tcga_counts[, keep_tcga, drop = FALSE]
}

log_dims("TCGA (after project subset)", tcga_counts)

if (!is.null(tcga_project_filter) && ncol(tcga_counts) == 0L) {
  stop(sprintf("[FATAL] TCGA project filter %s left 0 samples.", paste(tcga_project_filter, collapse=", ")))
}

if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA-BRCA filter left 0 samples.")
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ Breast filter left 0 samples.")

common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
cat(sprintf("[INFO] Shared genes: %d\n", length(common)))
if (length(common) < 1000) cat("[WARN] Few shared genes; check Ensembl IDs\n")

tcga_counts <- tcga_counts[common, , drop = FALSE]
dsmz_counts <- dsmz_counts[common, , drop = FALSE]

log_dims("TCGA (post-intersect)", tcga_counts)
log_dims("DSMZ (post-intersect)", dsmz_counts)

Xc_raw <- cbind(tcga_counts, dsmz_counts)
if (any(duplicated(rownames(Xc_raw)))) {
  dup_genes <- rownames(Xc_raw)[duplicated(rownames(Xc_raw))]
  cat(sprintf("[INFO] Collapsing %d duplicate genes in merged matrix: %s\n", 
              length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
  Xc_raw <- rowsum(Xc_raw, rownames(Xc_raw), reorder = TRUE)
}

log_dims("Merged (Xc_raw)", Xc_raw)

if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA subset has 0 samples.")
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ subset has 0 samples.")
if (nrow(Xc_raw) == 0L) stop("[FATAL] No common genes after intersection.")
if (ncol(Xc_raw) == 0L) stop("[FATAL] Merged matrix has 0 samples.")

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

dir.create(file.path(outdir, "shared_genes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "unshared_genes"), showWarnings = FALSE, recursive = TRUE)
writeLines(common, file.path(outdir, "shared_genes", "shared_genes.txt"))
writeLines(setdiff(rownames(tcga_counts_raw), common), file.path(outdir, "unshared_genes", "tcga_unshared_genes.txt"))
writeLines(setdiff(rownames(dsmz_counts_raw), common), file.path(outdir, "unshared_genes", "dsmz_unshared_genes.txt"))


# ---- Independent VST normalization per dataset first ----
cat("[INFO] VST per-dataset (TCGA only)...\n")
dds_tcga_ind <- DESeqDataSetFromMatrix(round(tcga_counts), data.frame(row.names = colnames(tcga_counts)), design = ~ 1)
vsd_tcga_ind <- vst(dds_tcga_ind, blind = TRUE)
V_tcga_ind   <- assay(vsd_tcga_ind)
saveRDS(V_tcga_ind, file.path(outdir, "VST_tcga_only.rds"))

cat("[INFO] VST per-dataset (DSMZ only)...\n")
dds_dsmz_ind <- DESeqDataSetFromMatrix(round(dsmz_counts), data.frame(row.names = colnames(dsmz_counts)), design = ~ 1)
vsd_dsmz_ind <- vst(dds_dsmz_ind, blind = TRUE)
V_dsmz_ind   <- assay(vsd_dsmz_ind)
saveRDS(V_dsmz_ind, file.path(outdir, "VST_dsmz_only.rds"))

# Apply purity adjustment to TCGA VST data if available
purity <- load_purity_data(tcga_se, NULL)
if (!is.null(purity)) {
  cat("[INFO] Applying purity adjustment to TCGA VST data...\n")
  adj <- purity_adjust_log(V_tcga_ind, purity)
  V_tcga_ind <- adj$V
  
  # always record what was flagged
  writeLines(adj$dropped, file.path(outdir, "purity_associated_genes.txt"))
  
  if (drop_purity_genes_permanently && length(adj$dropped)) {
    V_tcga_ind <- V_tcga_ind[setdiff(rownames(V_tcga_ind), adj$dropped), , drop = FALSE]
    message(sprintf("[INFO] Permanently removed %d purity-associated genes from TCGA VST.", length(adj$dropped)))
  }
}

# Optional consistency patch: Make DSMZ data consistent if purity genes were permanently dropped
if (drop_purity_genes_permanently && exists("adj")) {
  dsmz_counts <- dsmz_counts[setdiff(rownames(dsmz_counts), adj$dropped), , drop = FALSE]
  message(sprintf("[INFO] Made DSMZ data consistent by removing %d purity-associated genes.", length(adj$dropped)))
}

# ---- Merge independently normalized datasets ----
cat("[INFO] Merging independently normalized datasets...\n")
common_genes <- intersect(rownames(V_tcga_ind), rownames(V_dsmz_ind))
V_tcga_norm <- V_tcga_ind[common_genes, , drop = FALSE]
V_dsmz_norm <- V_dsmz_ind[common_genes, , drop = FALSE]

# Merge normalized datasets
V_merged <- cbind(V_tcga_norm, V_dsmz_norm)
log_dims("Merged normalized", V_merged)

# ---- Re-normalize merged data for ComBat_seq ----
cat("[INFO] Re-normalizing merged data for ComBat_seq...\n")
batch <- factor(c(rep("TCGA", ncol(V_tcga_norm)), rep("DSMZ", ncol(V_dsmz_norm))))
coldata <- tibble(
  sample = colnames(V_merged),
  dataset = ifelse(colnames(V_merged) %in% colnames(V_tcga_norm), "TCGA", "DSMZ")
) %>% column_to_rownames("sample")

# Convert back to counts for ComBat_seq (approximate inverse VST)
# Note: This is an approximation since VST is not easily invertible
# For ComBat_seq, we need raw counts, so we'll use the original counts
all_counts <- Xc_raw[common_genes, , drop = FALSE]

if (any(duplicated(rownames(all_counts)))) {
  all_counts <- rowsum(all_counts, rownames(all_counts), reorder = TRUE)
}

# Apply ComBat_seq to raw counts
cat("[INFO] Applying ComBat_seq to raw counts...\n")
set.seed(42)
all_counts_combat <- sva::ComBat_seq(counts = as.matrix(all_counts), 
                                     batch = batch, 
                                     group = NULL,
                                     covar_mod = NULL,
                                     full_mod = TRUE,
                                     shrink = FALSE,
                                     shrink.disp = FALSE,
                                     gene.subset.n = NULL
                                    )

# Convert ComBat_seq corrected counts back to VST
cat("[INFO] Converting ComBat_seq corrected counts to VST...\n")
dds_combat <- DESeqDataSetFromMatrix(round(all_counts_combat), coldata, design = ~ 1)
vsd_combat <- vst(dds_combat, blind = TRUE)
V_adj <- assay(vsd_combat)

# Keep old variable name for later code that still expects `V`
V <- V_adj

# Split back into datasets
V_tcga <- V_adj[, colnames(V_tcga_norm), drop = FALSE]
V_dsmz <- V_adj[, colnames(V_dsmz_norm), drop = FALSE]

log_dims("TCGA (ComBat_seq corrected)", V_tcga)
log_dims("DSMZ (ComBat_seq corrected)", V_dsmz)

stopifnot(!is.null(colnames(V_tcga)), !is.null(colnames(V_dsmz)))

# PCA before vs after batch correction
cat("[INFO] PCA before and after batch correction...\n")
make_pca_plot <- function(M, lab_batch, title, path_pdf, n_top = 5000) {
  # use top variable genes to stabilize PCA
  sel <- top_var_genes(M, n = min(n_top, nrow(M)))
  X   <- t(scale(M[sel, , drop = FALSE]))
  pr  <- prcomp(X, center = FALSE, scale. = FALSE)

  # percent variance
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  df <- data.frame(PC1 = pr$x[,1], PC2 = pr$x[,2],
                   batch = lab_batch,
                   sample = rownames(pr$x),
                   stringsAsFactors = FALSE)
  safe_pdf(path_pdf, {
    layout(matrix(c(1,2), nrow = 1))
    plot(df$PC1, df$PC2, pch = 19, cex = 0.7, col = as.factor(df$batch),
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title)
    legend("topright", legend = levels(as.factor(df$batch)),
           col = seq_along(levels(as.factor(df$batch))), pch = 19, cex = 0.8)
    barplot(var_expl[1:20], las = 2, main = "Explained variance (top 20 PCs)",
            ylab = "% variance", xlab = "PC")
    layout(1)
  })
  invisible(list(scores = df, var_expl = var_expl))
}

# Create batch labels for PCA
batch_labels <- factor(ifelse(colnames(V_adj) %in% colnames(V_tcga_norm), "TCGA", "DSMZ"))

pca_pre  <- make_pca_plot(V_merged, batch_labels, "PCA (pre-ComBat_seq)",  file.path(outdir, "pca_pre_batch.pdf"))
pca_post <- make_pca_plot(V_adj,    batch_labels, "PCA (post-ComBat_seq)", file.path(outdir, "pca_post_batch.pdf"))

saveRDS(V_tcga, file.path(outdir, "VST_tcga_brca.rds"))
saveRDS(V_dsmz, file.path(outdir, "VST_dsmz_brca.rds"))

# ---------- PAM50 SUBSETTING ----------
cat("[INFO] Subsetting to PAM50 genes...\n")
ens2sym_union <- bind_rows(
  tryCatch({
    cat("[INFO] Attempting gene symbol mapping via org.Hs.eg.db...\n")
    result <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(strip_ensver(rownames(V_adj))),
                         keytype = "ENSEMBL", columns = c("SYMBOL")) %>%
      filter(!is.na(SYMBOL) & nzchar(SYMBOL)) %>%
      distinct(ENSEMBL, .keep_all = TRUE) %>%
      transmute(ensembl = ENSEMBL, symbol = SYMBOL)
    cat(sprintf("[INFO] org.Hs.eg.db mapping successful: %d mappings\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("[WARN] org.Hs.eg.db mapping failed: %s\n", e$message))
    NULL
  }),
  as_df(dsmz_raw) %>%
    transmute(ensembl = strip_ensver(Ensembl_ID), symbol = as.character(gene_name)) %>%
    filter(!is.na(symbol) & nzchar(symbol)) %>%
    distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE),
  tryCatch({
    cat("[INFO] Attempting gene symbol mapping via TCGA rowData...\n")
    rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se))
    symcol <- intersect(tolower(colnames(rd)), c("gene_name","symbol","hgnc_symbol"))
    if (length(symcol)) {
      result <- data.frame(
        ensembl = strip_ensver(rownames(rd)),
        symbol = as.character(rd[[colnames(rd)[match(symcol[1], tolower(colnames(rd)))] ]]),
        stringsAsFactors = FALSE
      ) %>% filter(!is.na(symbol) & symbol != "") %>%
        distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)
      cat(sprintf("[INFO] TCGA rowData mapping successful: %d mappings using column '%s'\n", 
                  nrow(result), symcol[1]))
      result
    } else {
      cat("[WARN] TCGA rowData mapping failed: no suitable symbol column found\n")
      NULL
    }
  }, error = function(e) {
    cat(sprintf("[WARN] TCGA rowData mapping failed: %s\n", e$message))
    NULL
  })
) %>% distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)
rownames(ens2sym_union) <- ens2sym_union$ensembl

cat(sprintf("[INFO] Final gene symbol mapping summary: %d total mappings\n", nrow(ens2sym_union)))

M_pam50 <- subset_to_pam50(V_adj, ens2sym_union, pam50_symbols, collapse_fun = "mean")
log_dims("PAM50 Matrix", M_pam50)
saveRDS(M_pam50, file.path(outdir, "PAM50_matrix_VST_adj.rds"))
missing_syms <- setdiff(toupper(pam50_symbols), rownames(M_pam50))
cat(sprintf("[PAM50] Matched %d/%d genes\n", nrow(M_pam50), length(pam50_symbols)))
if (length(missing_syms)) {
  cat("[PAM50] Missing symbols:", paste(missing_syms, collapse=", "), "\n")
  writeLines(missing_syms, file.path(outdir, "PAM50_missing_symbols.txt"))
}

marker_ens <- ens2sym_union$ensembl[tolower(ens2sym_union$symbol) %in% tolower(marker_symbols)]
marker_ens <- intersect(unique(marker_ens), rownames(V_adj))
if (length(marker_ens) == 0L) {
  stop("[FATAL] No marker genes found in VST matrix. Check gene ID mapping.")
} else {
  cat(sprintf("[INFO] Found %d/%d marker genes in VST matrix.\n", 
              length(marker_ens), length(marker_symbols)))
  missing_markers <- setdiff(tolower(marker_symbols), tolower(ens2sym_union$symbol[ens2sym_union$ensembl %in% marker_ens]))
  if (length(missing_markers) > 0) {
    cat(sprintf("[WARN] Missing %d markers: %s\n", 
                length(missing_markers), paste(missing_markers, collapse=", ")))
  }
}

# ---------- (1) UNSUPERVISED DISCOVERY CLUSTERING (K-MEANS ON PCs) ----------
set.seed(42)

# Quality control and logging for adjusted matrix
cat("V_adj dims:", nrow(V_adj), "genes x", ncol(V_adj), "samples\n")
stopifnot(all(is.finite(V_adj)))
# variance check
nzv <- sum(matrixStats::rowVars(V_adj) > 0)
cat("Genes with non-zero variance:", nzv, "\n")

# Use robust PCA function instead of manual PCA
cat("[INFO] Computing robust PCA on adjusted matrix...\n")
pcs <- make_pcs(V_adj, n_hvg = 3000, max_pc = 30)
PC <- pcs$PC
kPC <- pcs$kPC
hvgs <- pcs$hvgs
var_expl <- (pcs$pr$sdev^2) / sum(pcs$pr$sdev^2) * 100
message(sprintf("[INFO] Using %d PCs explaining %.1f%% variance", kPC, sum(var_expl[1:kPC])))

# Save HVGs used
writeLines(hvgs, file.path(outdir, "HVGs_used_for_discovery_3k.txt"))

if (nrow(PC) > 5000) {
  cat("[WARN] Large sample set (%d) for k-means; consider subsampling.\n", nrow(PC))
}

# [NEW] Add scree plots for visual inspection
safe_pdf(file.path(outdir, "PCA_scree_plot.pdf"), {
  par(mfrow = c(1, 2))
  
  # Plot 1: Standard Scree Plot (Elbow)
  plot(var_expl[1:50], type = "b", 
       xlab = "Principal Component", 
       ylab = "Percent Variance Explained",
       main = "Scree Plot (Elbow)")
  abline(v = kPC, col = "red", lty = 2) # Add a line for your chosen kPC
  
  # Plot 2: Cumulative Variance
  plot(cumsum(var_expl)[1:50], type = "b", 
       xlab = "Number of Components", 
       ylab = "Cumulative Percent Variance Explained",
       main = "Cumulative Variance")
  abline(v = kPC, col = "red", lty = 2)
  abline(h = 90, col = "blue", lty = 2) # Example: 90% threshold
  
  par(mfrow = c(1, 1))
})

# Save PCA summary for reference
pca_summary <- data.frame(
  method = "Robust PCA",
  n_hvg = 3000,
  max_pc = 30,
  significant_PCs = kPC,
  total_variance_explained = sum(var_expl[1:kPC]),
  cumulative_variance_at_kPC = cumsum(var_expl)[kPC],
  stringsAsFactors = FALSE
)
write.csv(pca_summary, file.path(outdir, "PCA_robust_summary.csv"), row.names = FALSE)

k_grid <- 2:8
res_list <- list()
metrics <- lapply(k_grid, function(k) {
  km <- kmeans(PC, centers = k, nstart = 50, iter.max = 100)
  wcss <- sum(km$withinss)
  sil <- silhouette(km$cluster, dist(PC))
  sil_mean <- mean(sil[, "sil_width"])
  ch <- tryCatch(calinhara(PC, km$cluster), error = function(e) NA_real_)
  db <- tryCatch({
    m <- intCriteria(as.matrix(PC), as.integer(km$cluster), "davies_bouldin")
    as.numeric(m$davies_bouldin)
  }, error = function(e) NA_real_)
  res_list[[as.character(k)]] <<- list(km = km, sil = sil)
  data.frame(k = k, WCSS = wcss, Silhouette = sil_mean, CH = ch, DB = db)
})
metrics <- do.call(rbind, metrics)

write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), row.names = FALSE)

safe_pdf(file.path(outdir, "kmeans_elbow_HVG.pdf"), {
  plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
       ylab = "Within-Cluster SS", main = "Elbow (HVG space)")
})
safe_pdf(file.path(outdir, "kmeans_metrics_panel_HVG.pdf"), {
  par(mfrow = c(2,2))
  plot(metrics$k, metrics$Silhouette, type = "b", xlab = "k", 
       ylab = "Avg silhouette", main = "Silhouette (higher is better)")
  plot(metrics$k, metrics$CH, type = "b", xlab = "k", 
       ylab = "Calinski-Harabasz", main = "Calinski-Harabasz (higher)")
  plot(metrics$k, metrics$DB, type = "b", xlab = "k", 
       ylab = "Davies-Bouldin", main = "Davies-Bouldin (lower)")
  plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
       ylab = "WCSS", main = "Elbow (lower)")
  par(mfrow = c(1,1))
})

best_k <- with(metrics, {
  cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)]
  if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand
})
message(sprintf("[kmeans] best k by silhouette→CH = %s", best_k))

km_final <- res_list[[as.character(best_k)]]$km
sil_final <- res_list[[as.character(best_k)]]$sil
clusters_kmeans <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k)))

sample_is_tcga <- colnames(V_adj) %in% colnames(V_tcga)
dataset_lab <- ifelse(sample_is_tcga, "TCGA", "DSMZ")

cluster_df <- tibble(sample = rownames(PC),
                     dataset = dataset_lab,
                     cluster = as.character(clusters_kmeans))
write.csv(cluster_df, file.path(outdir, "clusters_kmeans.csv"), row.names = FALSE)

safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_HVG.pdf", best_k)), {
  plot(sil_final, main = sprintf("Silhouette (k = %d, HVG space)", best_k), cex.names = 0.8)
})

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

cat("Optimal k (k-means):", best_k, "with mean silhouette:", round(mean(sil_final[, "sil_width"]), 3), "\n")
print(table(clusters_kmeans, dataset_lab))

# UMAP for k-means (using same PC space as clustering)
set.seed(42)
emb <- uwot::umap(PC, n_neighbors = 20, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2")
emb$sample <- rownames(PC)
emb$dataset <- dataset_lab
emb$cluster <- clusters_kmeans

p_umap <- ggplot(emb, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle(sprintf("UMAP on %d PCs (k-means k=%d)", kPC, best_k))

ggsave(file.path(outdir, "umap_kmeans_PCs.pdf"), p_umap, width = 7, height = 6)
write.csv(emb, file.path(outdir, "umap_kmeans_PCs_coords.csv"), row.names = FALSE)

# ---------- (2) HIERARCHICAL CLUSTERING (PCs, k=5, EUCLIDEAN, AVERAGE LINKAGE) ----------
cat("[INFO] Performing hierarchical clustering on PCs...\n")
if (nrow(PC) < 5) {
  stop(sprintf("[FATAL] Too few samples (%d) for hierarchical clustering (need at least 5).", nrow(PC)))
}

# Distance metric choice justification:
# - Euclidean distance: Most common for gene expression clustering, forms compact spherical clusters
# - Manhattan distance: More sensitive to outliers, forms hyper-rectangular clusters
# - Correlation distance: Captures expression pattern similarity regardless of magnitude
# For breast cancer subtype discovery, Euclidean distance is preferred as it groups samples
# with similar overall expression profiles, which aligns with molecular subtype definitions.

# More specific logging about the inputs and methods
cat(sprintf("[INFO] Calculating distance matrix on %d PCs for %d samples using 'euclidean' method.\n", ncol(PC), nrow(PC)))

# Alternative distance metrics (uncomment to use):
# d <- dist(PC, method = "manhattan")  # Manhattan (L1 norm) - more robust to outliers
# d <- as.dist(1 - cor(t(PC)))        # Correlation distance - pattern similarity
# d <- dist(PC, method = "maximum")   # Chebyshev distance - maximum coordinate difference

d <- dist(PC, method = "euclidean")

cat("[INFO] Building dendrogram using 'average' linkage method.\n")
hc <- hclust(d, method = "average")

cat("[INFO] Cutting dendrogram into k=5 clusters.\n")
clusters_hc <- cutree(hc, k = 5)
subtype_labels <- c("HER2-high", "HER2-low", "LumA", "LumB", "Basal")

# Label clusters by correlating with TCGA PAM50 centroids
cat("[INFO] Attempting to label k=5 clusters by correlating with TCGA PAM50 centroids...\n")
tcga_sub_col_opts <- c("PAM50", "Subtype", "BRCA_Subtype", "molecular_subtype")
tcga_sub_col <- tcga_sub_col_opts[tcga_sub_col_opts %in% colnames(colData(tcga_se))][1]

if (!is.na(tcga_sub_col)) {
  cat(sprintf("[INFO] Found TCGA subtype annotation column: '%s'.\n", tcga_sub_col))
  
  # 1. Get TCGA labels and calculate TCGA centroids
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga)))
  cat(sprintf("[INFO] Found %d TCGA subtypes: %s\n", length(sub_lvls), paste(sub_lvls, collapse = ", ")))
  
  # [NEW] Log TCGA sample counts
  sub_counts <- table(lab_tcga)
  cat(sprintf("[INFO] TCGA subtype sample counts: %s\n", paste(names(sub_counts), sub_counts, sep = "=", collapse = ", ")))
  
  cat("[INFO] Calculating mean expression centroids for TCGA subtypes...\n")
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE))
  cat(sprintf("[INFO] TCGA centroid matrix dimensions: %d genes x %d subtypes.\n", nrow(sub_means), ncol(sub_means)))
  
  # 2. Calculate centroids for our new clusters
  cat("[INFO] Calculating mean expression centroids for discovered k=5 clusters...\n")
  
  # [NEW] Log discovered cluster sample counts
  cl_counts <- table(clusters_hc)
  cat(sprintf("[INFO] Discovered cluster sample counts (k=1-5): %s\n", paste(names(cl_counts), cl_counts, sep = "=", collapse = ", ")))
  
  cl_means <- sapply(1:5, function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE))
  cat(sprintf("[INFO] Discovered cluster centroid matrix dimensions: %d genes x %d clusters.\n", nrow(cl_means), ncol(cl_means)))
  
  # 3. Correlate the centroids
  g_common <- intersect(rownames(sub_means), hvgs)
  cat(sprintf("[INFO] Found %d common HVGs between datasets for correlation.\n", length(g_common)))
  
  # [NEW] Add safety check for low gene overlap
  if (length(g_common) < 500) {
    cat(sprintf("[WARN] Only %d common genes found for correlation. Results may be unstable.\n", length(g_common)))
  }
  
  cat("[INFO] Correlating cluster centroids with TCGA centroids (Spearman's rho).\n")
  corr <- cor(cl_means[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  
  # [NEW] Print the full correlation matrix to the log
  cat("[INFO] Correlation Matrix (Rows=Clusters, Cols=TCGA Subtypes):\n")
  corr_out <- capture.output(print(corr))
  cat(paste("    ", corr_out, collapse = "\n"), "\n")
  
  # 4. Assign best-matching label
  cluster_to_subtype <- apply(corr, 1, function(v) sub_lvls[which.max(v)])
  
  # [NEW] Create and print a clear mapping data frame
  final_label_names <- subtype_labels[match(cluster_to_subtype, sub_lvls)]
  label_map_df <- data.frame(Cluster_ID = 1:5,
                           Best_TCGA_Match = cluster_to_subtype,
                           Applied_Label = final_label_names)
  cat("[INFO] Final Cluster-to-Label Mapping:\n")
  map_out <- capture.output(print(label_map_df, row.names = FALSE))
  cat(paste("    ", map_out, collapse = "\n"), "\n")
  
  # 5. Apply new labels
  cat("[INFO] Relabeling cluster factors based on TCGA correlation.\n")
  clusters_hc <- factor(clusters_hc, labels = final_label_names)
  
} else {
  # This log is critical for debugging why labels might look generic
  cat(sprintf("[WARN] No TCGA subtype column found (looked for: %s).\n", paste(tcga_sub_col_opts, collapse = ", ")))
  cat(sprintf("[WARN] Using default fallback labels: %s\n", paste(subtype_labels, collapse = ", ")))
  
  # [NEW] Clearly state the fallback mapping
  label_map_str <- paste(1:5, subtype_labels, sep = "=", collapse = ", ")
  cat(sprintf("[WARN] Applying fallback mapping (Cluster=Label): %s\n", label_map_str))
  
  clusters_hc <- factor(clusters_hc, labels = subtype_labels)
}

# Save cluster assignments
cluster_hc_df <- tibble(sample = rownames(PC),
                        dataset = dataset_lab,
                        cluster = as.character(clusters_hc))
out_path_csv <- file.path(outdir, "clusters_hierarchical_HVG.csv")
cat(sprintf("[INFO] Saving hierarchical cluster assignments to: %s\n", out_path_csv))
write.csv(cluster_hc_df, out_path_csv, row.names = FALSE)

# Save dendrogram plot
out_path_pdf <- file.path(outdir, "dendrogram_hierarchical_HVG.pdf")
cat(sprintf("[INFO] Saving dendrogram plot to: %s\n", out_path_pdf))
safe_pdf(out_path_pdf, {
  plot(hc, main = "Hierarchical Clustering (HVG, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5)
  # Add line showing where the k=5 cut was made
  abline(h = hc$height[length(hc$height) - (5 - 1)], col = "red", lty = 2) 
  legend("topright", "k=5 Cut", col = "red", lty = 2, bty = "n")
})

# ---------- DYNAMIC TREE CUT ON HVG PC SPACE ----------
cat("[INFO] Performing Dynamic Tree Cut on PCs (HVG)...\n")
cat(sprintf("[INFO] Using %d PCs for %d samples.\n", ncol(PC), nrow(PC)))

cat("[INFO] Calculating Euclidean distance matrix...\n")
d_pc      <- dist(PC, method = "euclidean")

cat("[INFO] Building dendrogram with 'average' linkage (matching static HC)...\n")
hc_pc     <- hclust(d_pc, method = "average") 
distM_pc  <- as.matrix(d_pc)

# Calculate and log the adaptive min cluster size
minCl_pc  <- choose_minCluster(nrow(PC)) # adaptive min cluster size
cat(sprintf("[INFO] Adaptive 'minClusterSize' set to: %d (based on %d samples).\n", minCl_pc, nrow(PC)))

# Log key parameters for cutreeHybrid
cat(sprintf("[INFO] Running cutreeHybrid with deepSplit=%d, pamStage=%s, minClusterSize=%d...\n", 2, TRUE, minCl_pc))
dyn_pc <- cutreeHybrid(
  dendro           = hc_pc,
  distM            = distM_pc,
  deepSplit        = 2,       # 0..4: larger -> finer splits; 2 is a solid default
  pamStage         = TRUE,    # PAM refinement stabilizes tiny modules
  minClusterSize   = minCl_pc
)

clusters_hc_dyn <- factor(dyn_pc$labels) # 0 means "unassigned" (rare)

# Validate dynamic clustering results
if (is.null(dyn_pc$labels) || length(dyn_pc$labels) == 0) {
  cat("[ERROR] Dynamic tree cut failed - no labels returned. Skipping dynamic clustering.\n")
  clusters_hc_dyn <- factor(rep("Unassigned", nrow(PC)), levels = c("Unassigned", "1", "2", "3", "4", "5"))
  clusters_hc_dyn_named <- clusters_hc_dyn
} else {

# Log the results of the cut
cl_counts_raw <- table(clusters_hc_dyn)
k_dyn <- length(cl_counts_raw)
cat(sprintf("[INFO] Dynamic cut found %d clusters (including '0' for unassigned).\n", k_dyn))
cat(sprintf("[INFO] Raw dynamic cluster counts: %s\n", paste(names(cl_counts_raw), cl_counts_raw, sep = "=", collapse = ", ")))

if ("0" %in% names(cl_counts_raw)) {
  cat(sprintf("[WARN] %d samples were unassigned (cluster '0') by the dynamic cut.\n", cl_counts_raw["0"]))
}

# --- Map dynamic clusters -> PAM50 names ---
cat(sprintf("[INFO] Mapping %d dynamic clusters to TCGA PAM50 subtypes...\n", k_dyn))
clusters_hc_dyn_named <- clusters_hc_dyn # Initialize with raw cluster numbers

if (!is.na(tcga_sub_col)) {
  # Get numeric cluster levels (excluding '0' if it exists)
  cl_levels <- levels(clusters_hc_dyn)[levels(clusters_hc_dyn) != "0"]
  
  cat(sprintf("[INFO] Calculating centroids for %d dynamic clusters (excluding '0')...\n", length(cl_levels)))
  cl_means_dyn <- sapply(cl_levels, function(cn)
    rowMeans(V_adj[hvgs, clusters_hc_dyn == cn, drop = FALSE], na.rm = TRUE))
  cat(sprintf("[INFO] Dynamic cluster centroid matrix: %d genes x %d clusters.\n", nrow(cl_means_dyn), ncol(cl_means_dyn)))
  
  g_common <- intersect(rownames(sub_means), hvgs)
  cat(sprintf("[INFO] Using %d common HVGs for correlation.\n", length(g_common)))
  
  corr_dyn <- cor(cl_means_dyn[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  
  cat("[INFO] Correlation Matrix (Rows=Dynamic Clusters, Cols=TCGA Subtypes):\n")
  corr_out_dyn <- capture.output(print(corr_dyn))
  cat(paste("    ", corr_out_dyn, collapse = "\n"), "\n")
  
  # Map names
  map_dyn  <- apply(corr_dyn, 1, function(v) sub_lvls[which.max(v)])
  remap <- setNames(subtype_labels[match(map_dyn, sub_lvls)], colnames(corr_dyn))
  
  # Add "Unassigned" for cluster '0' if it exists
  if ("0" %in% levels(clusters_hc_dyn)) {
    remap <- c(remap, "0" = "Unassigned")
  }

  label_map_dyn_df <- data.frame(Dynamic_Cluster_ID = names(remap),
                                 Applied_Label = remap)
  cat("[INFO] Final Dynamic Cluster-to-Label Mapping:\n")
  map_out_dyn <- capture.output(print(label_map_dyn_df, row.names = FALSE))
  cat(paste("    ", map_out_dyn, collapse = "\n"), "\n")

  # Apply the new factor labels
  clusters_hc_dyn_named <- factor(remap[as.character(clusters_hc_dyn)], levels = unique(c("Unassigned", subtype_labels)))
  
} else {
  cat(sprintf("[WARN] No TCGA subtype column found. Dynamic clusters will be named 1, 2, 3... (or 0 if unassigned).\n"))
  # Keep clusters_hc_dyn_named as the raw numeric factors
}
} # Close the validation else block

# --- Save results ---
cluster_hc_dyn_df <- tibble(
  sample  = rownames(PC),
  dataset = dataset_lab,
  cluster = as.character(clusters_hc_dyn_named)
)

out_csv_dyn <- file.path(outdir, "clusters_hierarchical_dynamic_HVG.csv")
cat(sprintf("[INFO] Saving dynamic cluster assignments to: %s\n", out_csv_dyn))
write.csv(cluster_hc_dyn_df, out_csv_dyn, row.names = FALSE)

# --- Save plot ---
out_pdf_dyn <- file.path(outdir, "dendrogram_hierarchical_dynamic_HVG.pdf")
cat(sprintf("[INFO] Saving dynamic dendrogram with cluster bars to: %s\n", out_pdf_dyn))

safe_pdf(out_pdf_dyn, {
  cols <- WGCNA::labels2colors(as.numeric(factor(clusters_hc_dyn_named)))
  plotDendroAndColors(
    hc_pc,
    cols,
    groupLabels = "Dynamic cut (HVG PCs)",
    main = "Hierarchical clustering + Dynamic Tree Cut (HVG PCs)",
    dendroLabels = FALSE, # Set to TRUE if you want sample labels
    cex.dendroLabels = 0.3
  )
})

# QC: silhouette for dynamic labels (internal validity)
sil_dyn <- cluster::silhouette(as.integer(factor(clusters_hc_dyn)), d_pc)
safe_pdf(file.path(outdir, "silhouette_dynamic_HVG.pdf"), {
  plot(sil_dyn, main = sprintf("Silhouette - Dynamic cut (HVG PCs); mean=%.3f", mean(sil_dyn[, "sil_width"])))
})

# UMAP for hierarchical clustering (using same PC space as clustering)
set.seed(42)
emb_hc <- uwot::umap(PC, n_neighbors = 20, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
emb_hc <- as.data.frame(emb_hc); colnames(emb_hc) <- c("UMAP1","UMAP2")
emb_hc$sample <- rownames(PC)
emb_hc$dataset <- dataset_lab
emb_hc$cluster <- clusters_hc

p_hc_umap <- ggplot(emb_hc, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle(sprintf("UMAP on %d PCs (Hierarchical Clustering)", kPC))

ggsave(file.path(outdir, "umap_hierarchical_PCs.pdf"), p_hc_umap, width = 7, height = 6)
write.csv(emb_hc, file.path(outdir, "umap_hierarchical_PCs_coords.csv"), row.names = FALSE)

# ---------- (3) HIERARCHICAL CLUSTERING (PAM50, k=5, EUCLIDEAN, AVERAGE LINKAGE) ----------
cat("[INFO] Performing hierarchical clustering (PAM50)...\n")
X_pam50 <- t(scale(M_pam50))
if (nrow(X_pam50) < 5) {
  stop(sprintf("[FATAL] Too few samples (%d) for PAM50 hierarchical clustering (need at least 5).", nrow(X_pam50)))
}

# Distance metric choice for PAM50 genes:
# - Euclidean distance: Consistent with HVG clustering above, appropriate for scaled PAM50 expression
# - PAM50 genes are already curated for breast cancer subtypes, so Euclidean distance
#   effectively groups samples with similar PAM50 expression profiles

cat(sprintf("[INFO] Calculating distance matrix on %d PAM50 genes for %d samples using 'euclidean' method.\n", ncol(X_pam50), nrow(X_pam50)))
cat("[INFO] Using Euclidean distance for PAM50 hierarchical clustering (consistent with HVG analysis)\n")
d_pam50 <- dist(X_pam50, method = "euclidean")

cat("[INFO] Building PAM50 dendrogram using 'average' linkage method.\n")
hc_pam50 <- hclust(d_pam50, method = "average")

cat("[INFO] Cutting PAM50 dendrogram into k=5 clusters.\n")
clusters_hc_pam50 <- cutree(hc_pam50, k = 5)

# Label PAM50 clusters by correlating with TCGA PAM50 centroids
cat("[INFO] Attempting to label PAM50 k=5 clusters by correlating with TCGA PAM50 centroids...\n")
if (!is.na(tcga_sub_col)) {
  cat(sprintf("[INFO] Using existing TCGA subtype annotation column: '%s'.\n", tcga_sub_col))
  
  # Calculate centroids for PAM50 clusters
  cat("[INFO] Calculating mean PAM50 expression centroids for discovered k=5 clusters...\n")
  
  # [NEW] Log PAM50 cluster sample counts
  cl_counts_pam50 <- table(clusters_hc_pam50)
  cat(sprintf("[INFO] PAM50 cluster sample counts (k=1-5): %s\n", paste(names(cl_counts_pam50), cl_counts_pam50, sep = "=", collapse = ", ")))
  
  cl_means_pam50 <- sapply(1:5, function(cn) rowMeans(M_pam50[, clusters_hc_pam50 == cn, drop = FALSE], na.rm = TRUE))
  cat(sprintf("[INFO] PAM50 cluster centroid matrix dimensions: %d genes x %d clusters.\n", nrow(cl_means_pam50), ncol(cl_means_pam50)))
  
  # Find common PAM50 genes
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50))
  cat(sprintf("[INFO] Found %d common PAM50 genes between datasets for correlation.\n", length(g_common_pam50)))
  
  # [NEW] Add safety check for low gene overlap
  if (length(g_common_pam50) < 20) {
    cat(sprintf("[WARN] Only %d common PAM50 genes found for correlation. Results may be unstable.\n", length(g_common_pam50)))
  }
  
  cat("[INFO] Correlating PAM50 cluster centroids with TCGA centroids (Spearman's rho).\n")
  corr_pam50 <- cor(cl_means_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman")
  
  # [NEW] Print the full PAM50 correlation matrix to the log
  cat("[INFO] PAM50 Correlation Matrix (Rows=Clusters, Cols=TCGA Subtypes):\n")
  corr_out_pam50 <- capture.output(print(corr_pam50))
  cat(paste("    ", corr_out_pam50, collapse = "\n"), "\n")
  
  # Assign best-matching label
  cluster_to_subtype_pam50 <- apply(corr_pam50, 1, function(v) sub_lvls[which.max(v)])
  
  # [NEW] Create and print a clear PAM50 mapping data frame
  final_label_names_pam50 <- subtype_labels[match(cluster_to_subtype_pam50, sub_lvls)]
  label_map_df_pam50 <- data.frame(Cluster_ID = 1:5,
                                 Best_TCGA_Match = cluster_to_subtype_pam50,
                                 Applied_Label = final_label_names_pam50)
  cat("[INFO] Final PAM50 Cluster-to-Label Mapping:\n")
  map_out_pam50 <- capture.output(print(label_map_df_pam50, row.names = FALSE))
  cat(paste("    ", map_out_pam50, collapse = "\n"), "\n")
  
  # Apply new labels
  cat("[INFO] Relabeling PAM50 cluster factors based on TCGA correlation.\n")
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = final_label_names_pam50)
} else {
  cat(sprintf("[WARN] No TCGA subtype column found for PAM50 labeling (looked for: %s).\n", paste(tcga_sub_col_opts, collapse = ", ")))
  cat(sprintf("[WARN] Using default fallback labels for PAM50 clusters: %s\n", paste(subtype_labels, collapse = ", ")))
  
  # [NEW] Clearly state the PAM50 fallback mapping
  label_map_str_pam50 <- paste(1:5, subtype_labels, sep = "=", collapse = ", ")
  cat(sprintf("[WARN] Applying PAM50 fallback mapping (Cluster=Label): %s\n", label_map_str_pam50))
  
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = subtype_labels)
}

# Save PAM50 cluster assignments
cluster_hc_pam50_df <- tibble(sample = rownames(X_pam50),
                              dataset = dataset_lab,
                              cluster = as.character(clusters_hc_pam50))
out_path_csv_pam50 <- file.path(outdir, "clusters_hierarchical_PAM50.csv")
cat(sprintf("[INFO] Saving PAM50 hierarchical cluster assignments to: %s\n", out_path_csv_pam50))
write.csv(cluster_hc_pam50_df, out_path_csv_pam50, row.names = FALSE)

# Save PAM50 dendrogram plot
out_path_pdf_pam50 <- file.path(outdir, "dendrogram_hierarchical_PAM50.pdf")
cat(sprintf("[INFO] Saving PAM50 dendrogram plot to: %s\n", out_path_pdf_pam50))
safe_pdf(out_path_pdf_pam50, {
  plot(hc_pam50, main = "Hierarchical Clustering (PAM50, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5)
  # Add line showing where the k=5 cut was made
  abline(h = hc_pam50$height[length(hc_pam50$height) - (5 - 1)], col = "red", lty = 2) 
  legend("topright", "k=5 Cut", col = "red", lty = 2, bty = "n")
})

# ---------- DYNAMIC TREE CUT ON PAM50 SPACE ----------
cat("[INFO] Performing Dynamic Tree Cut on PAM50 space...\n")
cat(sprintf("[INFO] Using %d PAM50 genes for %d samples.\n", ncol(X_pam50), nrow(X_pam50)))

cat("[INFO] Calculating Euclidean distance matrix on PAM50 space...\n")
X_pam50   <- t(scale(M_pam50))
d_pam50   <- dist(X_pam50, method = "euclidean")

cat("[INFO] Building PAM50 dendrogram with 'average' linkage...\n")
hc_pam50d <- hclust(d_pam50, method = "average")

# Calculate and log the adaptive min cluster size for PAM50
minCl_p50 <- choose_minCluster(nrow(X_pam50))
cat(sprintf("[INFO] Adaptive 'minClusterSize' for PAM50 set to: %d (based on %d samples).\n", minCl_p50, nrow(X_pam50)))

# Log key parameters for cutreeHybrid on PAM50
cat(sprintf("[INFO] Running cutreeHybrid on PAM50 with deepSplit=%d, pamStage=%s, minClusterSize=%d...\n", 2, TRUE, minCl_p50))
dyn_p50 <- cutreeHybrid(
  dendro         = hc_pam50d,
  distM          = as.matrix(d_pam50),
  deepSplit      = 2,
  pamStage       = TRUE,
  minClusterSize = minCl_p50
)

clusters_hc_pam50_dyn <- factor(dyn_p50$labels)

# Log the results of the PAM50 cut
cl_counts_raw_p50 <- table(clusters_hc_pam50_dyn)
k_dyn_p50 <- length(cl_counts_raw_p50)
cat(sprintf("[INFO] PAM50 dynamic cut found %d clusters (including '0' for unassigned).\n", k_dyn_p50))
cat(sprintf("[INFO] Raw PAM50 dynamic cluster counts: %s\n", paste(names(cl_counts_raw_p50), cl_counts_raw_p50, sep = "=", collapse = ", ")))

if ("0" %in% names(cl_counts_raw_p50)) {
  cat(sprintf("[WARN] %d samples were unassigned (cluster '0') by the PAM50 dynamic cut.\n", cl_counts_raw_p50["0"]))
}

# Map dynamic clusters → PAM50 names using TCGA subtypes
cat(sprintf("[INFO] Mapping %d PAM50 dynamic clusters to TCGA PAM50 subtypes...\n", k_dyn_p50))
clusters_hc_pam50_dyn_named <- clusters_hc_pam50_dyn

if (!is.na(tcga_sub_col)) {
  # Get numeric cluster levels (excluding '0' if it exists)
  cl_levels_p50 <- levels(clusters_hc_pam50_dyn)[levels(clusters_hc_pam50_dyn) != "0"]
  
  cat(sprintf("[INFO] Calculating PAM50 centroids for %d dynamic clusters (excluding '0')...\n", length(cl_levels_p50)))
  cl_means_p50_dyn <- sapply(cl_levels_p50, function(cn)
    rowMeans(M_pam50[, clusters_hc_pam50_dyn == cn, drop = FALSE], na.rm = TRUE))
  cat(sprintf("[INFO] PAM50 dynamic cluster centroid matrix: %d genes x %d clusters.\n", nrow(cl_means_p50_dyn), ncol(cl_means_p50_dyn)))
  
  g_common_p50 <- intersect(rownames(sub_means), rownames(M_pam50))
  cat(sprintf("[INFO] Using %d common PAM50 genes for correlation.\n", length(g_common_p50)))
  
  corr_p50_dyn <- cor(cl_means_p50_dyn[g_common_p50, , drop = FALSE], sub_means[g_common_p50, , drop = FALSE], method = "spearman")
  
  cat("[INFO] PAM50 Dynamic Correlation Matrix (Rows=Dynamic Clusters, Cols=TCGA Subtypes):\n")
  corr_out_p50_dyn <- capture.output(print(corr_p50_dyn))
  cat(paste("    ", corr_out_p50_dyn, collapse = "\n"), "\n")
  
  map_p50_dyn  <- apply(corr_p50_dyn, 1, function(v) sub_lvls[which.max(v)])
  remap_p50    <- setNames(subtype_labels[match(map_p50_dyn, sub_lvls)], colnames(corr_p50_dyn))
  
  # Add "Unassigned" for cluster '0' if it exists
  if ("0" %in% levels(clusters_hc_pam50_dyn)) {
    remap_p50 <- c(remap_p50, "0" = "Unassigned")
  }

  label_map_p50_dyn_df <- data.frame(Dynamic_Cluster_ID = names(remap_p50),
                                     Applied_Label = remap_p50)
  cat("[INFO] Final PAM50 Dynamic Cluster-to-Label Mapping:\n")
  map_out_p50_dyn <- capture.output(print(label_map_p50_dyn_df, row.names = FALSE))
  cat(paste("    ", map_out_p50_dyn, collapse = "\n"), "\n")
  
  clusters_hc_pam50_dyn_named <- factor(remap_p50[as.character(clusters_hc_pam50_dyn)], levels = unique(c("Unassigned", subtype_labels)))
} else {
  cat(sprintf("[WARN] No TCGA subtype column found. PAM50 dynamic clusters will be named 1, 2, 3... (or 0 if unassigned).\n"))
  # Keep clusters_hc_pam50_dyn_named as the raw numeric factors
}

# --- Save PAM50 dynamic results ---
cluster_hc_p50_dyn_df <- tibble(
  sample  = rownames(X_pam50),
  dataset = dataset_lab,
  cluster = as.character(clusters_hc_pam50_dyn_named)
)

out_csv_p50_dyn <- file.path(outdir, "clusters_hierarchical_dynamic_PAM50.csv")
cat(sprintf("[INFO] Saving PAM50 dynamic cluster assignments to: %s\n", out_csv_p50_dyn))
write.csv(cluster_hc_p50_dyn_df, out_csv_p50_dyn, row.names = FALSE)

# --- Save PAM50 dynamic plot ---
out_pdf_p50_dyn <- file.path(outdir, "dendrogram_hierarchical_dynamic_PAM50.pdf")
cat(sprintf("[INFO] Saving PAM50 dynamic dendrogram with cluster bars to: %s\n", out_pdf_p50_dyn))

safe_pdf(out_pdf_p50_dyn, {
  cols <- WGCNA::labels2colors(as.numeric(factor(clusters_hc_pam50_dyn_named)))
  plotDendroAndColors(
    hc_pam50d,
    cols,
    groupLabels = "Dynamic cut (PAM50)",
    main = "Hierarchical clustering + Dynamic Tree Cut (PAM50)",
    dendroLabels = FALSE, # Set to TRUE if you want sample labels
    cex.dendroLabels = 0.3
  )
})
sil_p50_dyn <- cluster::silhouette(as.integer(factor(clusters_hc_pam50_dyn)), d_pam50)
safe_pdf(file.path(outdir, "silhouette_dynamic_PAM50.pdf"), {
  plot(sil_p50_dyn, main = sprintf("Silhouette - Dynamic cut (PAM50); mean=%.3f", mean(sil_p50_dyn[, "sil_width"])))
})

# UMAP for hierarchical clustering (PAM50)
set.seed(42)
emb_hc_pam50 <- uwot::umap(X_pam50, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb_hc_pam50 <- as.data.frame(emb_hc_pam50); colnames(emb_hc_pam50) <- c("UMAP1","UMAP2")
emb_hc_pam50$sample <- rownames(X_pam50)
emb_hc_pam50$dataset <- dataset_lab
emb_hc_pam50$cluster <- clusters_hc_pam50

p_hc_pam50_umap <- ggplot(emb_hc_pam50, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (Hierarchical Clustering, PAM50 space)")

ggsave(file.path(outdir, "umap_hierarchical_PAM50.pdf"), p_hc_pam50_umap, width = 7, height = 6)
write.csv(emb_hc_pam50, file.path(outdir, "umap_hierarchical_PAM50_coords.csv"), row.names = FALSE)

# ---------- BIOLOGICAL SIGNATURE-BASED LABELING ----------
cat("[INFO] Computing biological signature-based labels...\n")

# --- helpers: map symbols, score signatures on V_adj (genes x samples) ---
sym2ens <- setNames(ens2sym_union$ensembl, toupper(ens2sym_union$symbol))
to_ens  <- function(x) na.omit(sym2ens[toupper(x)])
score_sig <- function(M, symbols) {
  g <- intersect(to_ens(symbols), rownames(M))
  if (!length(g)) return(rep(NA_real_, ncol(M)))
  Z <- scale(t(M[g, , drop = FALSE]))  # samples x genes z
  rowMeans(Z, na.rm = TRUE)            # sample-wise score
}

# canonical marker sets
sig_ER    <- c("ESR1","PGR","BCL2","GATA3","XBP1","SLC39A6")
sig_HER2  <- c("ERBB2","GRB7","FGFR4")
sig_BASAL <- c("KRT5","KRT14","KRT17","EGFR")
sig_PROL  <- c("MKI67","CCNB1","BIRC5","UBE2C","PTTG1","CDC20")

# per-sample biological scores
S <- data.frame(
  sample = colnames(V_adj),
  ER     = score_sig(V_adj, sig_ER),
  HER2   = score_sig(V_adj, sig_HER2),
  BASAL  = score_sig(V_adj, sig_BASAL),
  PROL   = score_sig(V_adj, sig_PROL),
  row.names = colnames(V_adj),
  check.names = FALSE
)

# optional: add PAM50 centroid correlation (strong prior, still biological)
# build centroids from TCGA within your data if annotation exists
pam_subs <- c("Basal","HER2-enriched","Luminal A","Luminal B","Normal-like")
pam_cent <- NULL
if (exists("M_pam50") && ncol(M_pam50) > 0) {
  # if TCGA subtype column is present, derive centroids from it; else use literature order as a fallback
  tcga_col_try <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
  tcga_col <- tcga_col_try[tcga_col_try %in% colnames(colData(tcga_se))][1]
  if (!is.na(tcga_col)) {
    lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_col])
    subs <- intersect(pam_subs, unique(na.omit(lab_tcga)))
    pam_cent <- sapply(subs, function(s) rowMeans(M_pam50[, colnames(V_tcga)[lab_tcga==s], drop=FALSE], na.rm=TRUE))
  }
}
pam_corr <- matrix(NA_real_, nrow = ncol(V_adj), ncol = 0,
                   dimnames = list(colnames(V_adj), NULL))
if (!is.null(pam_cent)) {
  common <- intersect(rownames(M_pam50), rownames(pam_cent))
  Xs <- scale(t(M_pam50[common, , drop = FALSE])) # samples x genes (z)
  Cs <- scale(t(pam_cent[common, , drop = FALSE]))# centroids x genes (z)
  R  <- cor(Xs, t(Cs), method = "spearman", use = "pairwise.complete.obs")
  pam_corr <- R[rownames(R) %in% colnames(V_adj), , drop = FALSE]
}

# unified per-sample biological label (with confidence & Unknown)
biologic_label <- function(i, margin = 0.4, prol_hi = 0) {
  er   <- S[i, "ER"]; her2 <- S[i, "HER2"]; bas <- S[i, "BASAL"]; prol <- S[i, "PROL"]
  best <- c(ER=er, HER2=her2, BASAL=bas)
  if (all(!is.finite(best))) return(c(label="Unknown", conf=NA))

  # primary axis: which biology dominates?
  ord <- sort(best, decreasing = TRUE)
  if (length(ord) < 2 || (ord[1] - ord[2] < margin)) {
    primary <- names(ord)[1]
    # not separable enough → Unknown unless a strong PAM50 correlation exists
    if (ncol(pam_corr) > 0) {
      pc <- pam_corr[i, , drop=TRUE]
      if (is.finite(max(pc, na.rm=TRUE)) && max(pc, na.rm=TRUE) >= 0.3) {
        return(c(label = names(which.max(pc)), conf = max(pc, na.rm=TRUE)))
      }
    }
    return(c(label="Unknown", conf=NA))
  }
  primary <- names(ord)[1]

  # refine within ER-dominant: LumA vs LumB by proliferation
  if (primary == "ER") {
    lab <- if (is.finite(prol) && prol > prol_hi) "Luminal B" else "Luminal A"
    # boost confidence with PAM50 agreement when available
    conf <- if (ncol(pam_corr)>0 && all(c("Luminal A","Luminal B") %in% colnames(pam_corr))) {
      max(pam_corr[i, c("Luminal A","Luminal B")], na.rm=TRUE)
    } else ord[1] - ord[2]
    return(c(label=lab, conf=conf))
  }

  if (primary == "HER2") {
    # split HER2-enriched vs HER2-low by ER support
    lab <- if (is.finite(er) && er > 0) "HER2-low" else "HER2-enriched"
    conf <- if (ncol(pam_corr)>0 && "HER2-enriched" %in% colnames(pam_corr)) pam_corr[i, "HER2-enriched"] else ord[1] - ord[2]
    return(c(label=lab, conf=conf))
  }

  if (primary == "BASAL") {
    lab <- "Basal"
    conf <- if (ncol(pam_corr)>0 && "Basal" %in% colnames(pam_corr)) pam_corr[i, "Basal"] else ord[1] - ord[2]
    return(c(label=lab, conf=conf))
  }

  c(label="Unknown", conf=NA)
}

# per-sample calls
bio_calls <- t(sapply(seq_len(nrow(S)), function(ii) biologic_label(rownames(S)[ii])))
bio_calls <- as.data.frame(bio_calls, stringsAsFactors = FALSE)
bio_calls$sample <- rownames(S)

# roll up to cluster-level (majority vote + purity + min confidence); keep Unknown if weak
infer_cluster_label <- function(labels, conf, min_purity = 0.6, min_conf = 0.25) {
  tab <- sort(table(labels), decreasing = TRUE)
  if (!length(tab) || names(tab)[1] == "Unknown") return(c(label="Unknown", purity=NA, min_conf=NA))
  top <- names(tab)[1]; pur <- unname(tab[1] / sum(tab))
  mc <- min(conf[labels == top], na.rm = TRUE)
  if (!is.finite(mc)) mc <- NA_real_
  if (pur >= min_purity && (is.na(mc) || mc >= min_conf)) c(label=top, purity=pur, min_conf=mc) else c(label="Unknown", purity=pur, min_conf=mc)
}

  label_cluster_set <- function(cluster_factor, name){
  df <- merge(bio_calls, data.frame(sample = names(cluster_factor),
                                    cluster = as.character(cluster_factor),
                                    stringsAsFactors=FALSE),
              by="sample", all.x=TRUE)
  byc <- split(df, df$cluster)
  cl_info <- lapply(byc, function(dd) infer_cluster_label(dd$label, as.numeric(dd$conf)))
  cl_info <- do.call(rbind, cl_info)
  
  # Handle empty cluster information
  if (is.null(cl_info) || nrow(cl_info) == 0) {
    cat(sprintf("[WARN] No valid cluster information for '%s' - skipping biological labeling.\n", name))
    cl_info <- data.frame(cluster = character(0), label = character(0), 
                         purity = numeric(0), min_conf = numeric(0), 
                         set = character(0), stringsAsFactors = FALSE)
  } else {
    # Convert to data frame and ensure consistent column types
    cl_info <- as.data.frame(cl_info, stringsAsFactors = FALSE)
    cl_info$cluster <- rownames(cl_info)
    
    # Ensure purity and min_conf are numeric
    cl_info$purity <- as.numeric(as.character(cl_info$purity))
    cl_info$min_conf <- as.numeric(as.character(cl_info$min_conf))
    
    # Add set column
    cl_info$set <- name
    
    # Reorder columns
    cl_info <- cl_info[, c("cluster", "label", "purity", "min_conf", "set"), drop = FALSE]
  }
  
  list(per_sample = df, per_cluster = cl_info)
}

# apply to all available clusterings
label_sets <- list(
  kmeans        = clusters_kmeans,
  hc_static     = clusters_hc,
  hc_dynamic    = clusters_hc_dyn_named,
  hc_pam50_dyn  = clusters_hc_pam50_dyn_named
)

# Add HDBSCAN results if they exist
if (exists("clusters_hdb_num")) {
  label_sets$hdbscan_pcs <- clusters_hdb_num
}
if (exists("clusters_hdb_pam50_num")) {
  label_sets$hdbscan_pam50 <- clusters_hdb_pam50_num
}

dir.create(file.path(outdir, "bio_labels"), FALSE, TRUE)
pcs <- list(); cls <- list()
for (nm in names(label_sets)) {
  lab <- label_sets[[nm]]
  if (is.null(lab) || length(lab) == 0) {
    cat(sprintf("[WARN] Skipping biological labeling for '%s': no valid clusters.\n", nm))
    next
  }
  res <- label_cluster_set(lab, nm)
  pcs[[nm]] <- res$per_sample
  cls[[nm]] <- res$per_cluster
  write.csv(res$per_sample,  file.path(outdir, "bio_labels", sprintf("per_sample_%s.csv", nm)),  row.names=FALSE)
  write.csv(res$per_cluster, file.path(outdir, "bio_labels", sprintf("per_cluster_%s.csv", nm)), row.names=FALSE)
}

# master tables
if (length(pcs)) write.csv(dplyr::bind_rows(pcs, .id="set"), file.path(outdir, "bio_labels", "per_sample_all.csv"), row.names=FALSE)
if (length(cls)) write.csv(dplyr::bind_rows(cls, .id="set"), file.path(outdir, "bio_labels", "per_cluster_all.csv"), row.names=FALSE)

# ---------- GENE PAIR TREND ANALYSIS ----------
cat("[INFO] Performing gene pair trend analysis...\n")

trend_pairs <- list(
  c("ESR1","PGR"),
  c("ERBB2","GRB7"),
  c("KRT5","KRT14"),
  c("MKI67","CCNB1")
)

sym2ens <- setNames(ens2sym_union$ensembl, toupper(ens2sym_union$symbol))
to_ens  <- function(sym) sym2ens[toupper(sym)]

mk_trend_df <- function(M, x_sym, y_sym, labels, label_name){
  x_ens <- to_ens(x_sym); y_ens <- to_ens(y_sym)
  ok <- !is.na(x_ens) && !is.na(y_ens) && x_ens %in% rownames(M) && y_ens %in% rownames(M)
  if(!ok) return(NULL)
  tibble(
    sample = colnames(M),
    x = as.numeric(M[x_ens, ]),
    y = as.numeric(M[y_ens, ]),
    cluster = factor(as.character(labels[ colnames(M) ])),
    label_name = label_name,
    x_sym = x_sym, y_sym = y_sym
  ) %>% filter(is.finite(x) & is.finite(y))
}

fit_trends <- function(df){
  out <- list()
  # global
  m <- lm(y ~ x, data=df)
  out[[1]] <- data.frame(scope="global", slope=coef(m)[["x"]],
                         r2=summary(m)$r.squared,
                         p=coef(summary(m))[["x","Pr(>|t|)"]],
                         n=nrow(df))
  # per-cluster
  pc <- lapply(split(df, df$cluster, drop=TRUE), function(dd){
    if (nrow(dd) < 3) return(NULL)
    m <- lm(y ~ x, data=dd)
    data.frame(scope=paste0("cluster:", as.character(dd$cluster[1])),
               slope=coef(m)[["x"]],
               r2=summary(m)$r.squared,
               p=coef(summary(m))[["x","Pr(>|t|)"]],
               n=nrow(dd))
  })
  do.call(rbind, c(list(out[[1]]), pc))
}

# Run for each label set
dir.create(file.path(outdir, "trend_checks"), showWarnings=FALSE, recursive=TRUE)

label_sets <- list(
  kmeans       = clusters_kmeans,
  hc_static    = clusters_hc,
  hc_dynamic   = clusters_hc_dyn_named
)

# Add HDBSCAN results if they exist
if (exists("clusters_hdb_num")) {
  label_sets$hdbscan_pcs <- clusters_hdb_num
}
if (exists("clusters_hdb_pam50_num")) {
  label_sets$hdbscan_p50 <- clusters_hdb_pam50_num
}

all_stats <- list()

for (labname in names(label_sets)) {
  lab <- label_sets[[labname]]
  if (is.null(lab)) next
  for (pr in trend_pairs) {
    df <- mk_trend_df(V_adj, pr[1], pr[2], lab, labname)
    if (is.null(df) || nrow(df) < 10) next

    # stats
    st <- fit_trends(df)
    st$pair <- paste(pr, collapse="~")
    st$labels <- labname
    all_stats[[paste(labname, pr[1], pr[2], sep="_")]] <- st

    # plot
    p <- ggplot(df, aes(x, y, color = cluster)) +
      geom_point(alpha=.7, size=1.6) +
      geom_smooth(method="lm", se=FALSE, linewidth=0.9) +                  # per-cluster lines
      geom_smooth(aes(group=1), method="lm", se=FALSE, linewidth=1.2,      # global line
                  color="black", linetype="dashed") +
      theme_bw() +
      labs(title = sprintf("%s: %s vs %s", labname, pr[2], pr[1]),
           x = pr[1], y = pr[2],
           subtitle = "Dashed black = global OLS; colored = per-cluster OLS")

    ggsave(file.path(outdir, "trend_checks",
                     sprintf("trend_%s_%s_vs_%s.pdf", labname, pr[2], pr[1])),
           p, width=6.8, height=5.2)
  }
}

if (length(all_stats)){
  stats_df <- dplyr::bind_rows(all_stats, .id="id")
  write.csv(stats_df, file.path(outdir, "trend_checks", "trend_summary.csv"), row.names=FALSE)
}

# ---------- (4) DIFFERENTIAL EXPRESSION ANALYSIS ----------
cat("[INFO] Performing differential expression analysis...\n")

# Make sure the factor has 'Basal' as the reference level and includes ONLY present levels
subtype_levels <- c("Basal","HER2-high","HER2-low","LumA","LumB")
cl_here <- factor(clusters_hc, levels = subtype_levels)
cl_here <- droplevels(cl_here)   # drop absent levels to avoid empty comparisons
coldata$subtype <- cl_here

# Sanitize factor levels for DESeq2 (optional hygiene)
levels(coldata$subtype) <- gsub("[^A-Za-z0-9_.]", "_", levels(coldata$subtype))

dds_de <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ subtype)
dds_de <- DESeq(dds_de)

present_lvls <- levels(coldata$subtype)
targets <- setdiff(present_lvls, "Basal")
if (!length(targets)) {
  stop("[DE] No non-Basal subtypes present to contrast against Basal.")
}

de_results <- list()
for (sub in targets) {
  res <- results(dds_de, contrast = c("subtype", sub, "Basal"), alpha = 0.05)
  res <- as.data.frame(res) %>%
    mutate(ensembl = rownames(res),
           symbol  = ens2sym_union$symbol[match(ensembl, ens2sym_union$ensembl)]) %>%
    filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange), abs(log2FoldChange) > 1) %>%
    arrange(padj)
  de_results[[sub]] <- res
}
dir.create(file.path(outdir, "de_results"), showWarnings = FALSE, recursive = TRUE)
for (sub in names(de_results)) {
  write.csv(de_results[[sub]], file.path(outdir, "de_results", sprintf("de_%s_vs_Basal.csv", sub)), row.names = FALSE)
}

# Heatmap of top DE genes
top_de_genes <- unique(unlist(lapply(de_results, function(res) head(res$ensembl, 10))))
if (length(top_de_genes) > 0) {
  V_de <- V_adj[top_de_genes, , drop = FALSE]
  if (any(duplicated(rownames(V_de)))) {
    rownames(V_de) <- make.unique(rownames(V_de))
  }
  lab <- ens2sym_union$symbol[match(rownames(V_de), ens2sym_union$ensembl)]
  lab[is.na(lab) | lab == ""] <- rownames(V_de)
  rownames(V_de) <- lab
  ann <- data.frame(dataset = dataset_lab, subtype = clusters_hc, row.names = colnames(V_de))
  safe_pdf(file.path(outdir, "heatmap_de_genes.pdf"), {
    pheatmap(V_de[, order(clusters_hc), drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_col = ann,
             main = "Top DE Genes by Subtype")
  })
}

# ---------- (5) UMAP + HDBSCAN ITERATIVE SUBTYPE ASSIGNMENT (PCs) ----------
cat("[INFO] Performing iterative UMAP + HDBSCAN on PCs...\n")
set.seed(42)
if (nrow(PC) < 20) {
  stop(sprintf("[FATAL] Too few samples (%d) for UMAP+HDBSCAN (need at least 20).", nrow(PC)))
}
if (nrow(PC) > 5000) {
  cat("[WARN] Large sample set (%d) for UMAP+HDBSCAN; consider subsampling.\n", nrow(PC))
}

emb1 <- uwot::umap(PC, n_neighbors = 20, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
minPts1 <- choose_minPts(nrow(emb1))
cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d samples\n", minPts1, nrow(emb1)))
hdb1 <- hdbscan(emb1, minPts = minPts1)
clusters_hdb <- rep("Unassigned", nrow(PC))
basal_idx <- which(hdb1$cluster == which.max(table(hdb1$cluster[hdb1$cluster != 0])))
clusters_hdb[basal_idx] <- "Basal"

non_basal <- which(clusters_hdb == "Unassigned")
if (length(non_basal) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification.\n", length(non_basal))
} else {
  emb2 <- uwot::umap(PC[non_basal, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
  minPts2 <- choose_minPts(nrow(emb2))
  cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-Basal samples\n", minPts2, nrow(emb2)))
  hdb2 <- hdbscan(emb2, minPts = minPts2)
  her2_high_idx <- non_basal[which(hdb2$cluster == which.max(table(hdb2$cluster[hdb2$cluster != 0])))]
  clusters_hdb[her2_high_idx] <- "HER2-high"
}

non_her2_high <- which(clusters_hdb == "Unassigned")
if (length(non_her2_high) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification.\n", length(non_her2_high))
} else {
  emb3 <- uwot::umap(PC[non_her2_high, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
  minPts3 <- choose_minPts(nrow(emb3))
  cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-HER2-high samples\n", minPts3, nrow(emb3)))
  hdb3 <- hdbscan(emb3, minPts = minPts3)
  her2_low_idx <- non_her2_high[which(hdb3$cluster == which.max(table(hdb3$cluster[hdb3$cluster != 0])))]
  clusters_hdb[her2_low_idx] <- "HER2-low"
}

luminal_idx <- which(clusters_hdb == "Unassigned")
if (length(luminal_idx) < 2) {
  cat("[WARN] Too few Luminal samples (%d); assigning all as LumA.\n", length(luminal_idx))
  clusters_hdb[luminal_idx] <- "LumA"
} else {
  emb_lum <- uwot::umap(PC[luminal_idx, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
  km_lum <- kmeans(emb_lum, centers = 2, nstart = 50)
  lum_clusters <- ifelse(km_lum$cluster == 1, "LumA", "LumB")
  clusters_hdb[luminal_idx] <- lum_clusters
}

clusters_hdb <- factor(clusters_hdb, levels = subtype_labels)

# ------- Generic numeric relabeling for HDBSCAN results -------
# Turn stringy subtype labels ("Basal","HER2-high",...) into 1..5,
# leaving Unassigned as NA (so it doesn't imply biology).
relabel_hdb_numeric <- function(x, keep_unassigned = TRUE) {
  x <- as.character(x)
  unass <- x == "Unassigned"
  labs  <- sort(unique(x[!unass]))              # stable, deterministic order
  map   <- setNames(as.character(seq_along(labs)), labs)
  y     <- ifelse(unass & keep_unassigned, NA, map[x])
  factor(y, levels = as.character(seq_len(min(5, length(labs)))))
}

# --- Apply to HDBSCAN on PCs ---
clusters_hdb_num <- relabel_hdb_numeric(clusters_hdb, keep_unassigned = TRUE)

hdb_df <- tibble(sample = rownames(PC),
                 dataset = dataset_lab,
                 subtype = clusters_hdb_num)      # now 1..5 (NA for Unassigned)
write.csv(hdb_df, file.path(outdir, "clusters_hdbscan_PCs.csv"), row.names = FALSE)

emb_hdb <- as.data.frame(emb1); colnames(emb_hdb) <- c("UMAP1","UMAP2")
emb_hdb$sample <- rownames(PC)
emb_hdb$dataset <- dataset_lab
emb_hdb$subtype <- clusters_hdb_num

p_hdb_umap <- ggplot(emb_hdb, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle(sprintf("UMAP on %d PCs (HDBSCAN clusters: 1–5, NA=unassigned)", kPC))

ggsave(file.path(outdir, "umap_hdbscan_PCs.pdf"), p_hdb_umap, width = 7, height = 6)
write.csv(emb_hdb, file.path(outdir, "umap_hdbscan_PCs_coords.csv"), row.names = FALSE)

# ---------- (6) UMAP + HDBSCAN ITERATIVE SUBTYPE ASSIGNMENT (PAM50) ----------
cat("[INFO] Performing iterative UMAP + HDBSCAN (PAM50)...\n")
set.seed(42)
if (nrow(X_pam50) < 20) {
  stop(sprintf("[FATAL] Too few samples (%d) for PAM50 UMAP+HDBSCAN (need at least 20).", nrow(X_pam50)))
}
if (nrow(X_pam50) > 5000) {
  cat("[WARN] Large sample set (%d) for PAM50 UMAP+HDBSCAN; consider subsampling.\n", nrow(X_pam50))
}

emb1_pam50 <- uwot::umap(X_pam50, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
minPts1_pam50 <- choose_minPts(nrow(emb1_pam50))
cat(sprintf("[HDBSCAN PAM50] Using minPts = %d for n = %d samples\n", minPts1_pam50, nrow(emb1_pam50)))
hdb1_pam50 <- hdbscan(emb1_pam50, minPts = minPts1_pam50)
clusters_hdb_pam50 <- rep("Unassigned", nrow(X_pam50))
basal_idx_pam50 <- which(hdb1_pam50$cluster == which.max(table(hdb1_pam50$cluster[hdb1_pam50$cluster != 0])))
clusters_hdb_pam50[basal_idx_pam50] <- "Basal"

non_basal_pam50 <- which(clusters_hdb_pam50 == "Unassigned")
if (length(non_basal_pam50) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification (PAM50).\n", length(non_basal_pam50))
} else {
  emb2_pam50 <- uwot::umap(X_pam50[non_basal_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  minPts2_pam50 <- choose_minPts(nrow(emb2_pam50))
  cat(sprintf("[HDBSCAN PAM50] Using minPts = %d for n = %d non-Basal samples\n", minPts2_pam50, nrow(emb2_pam50)))
  hdb2_pam50 <- hdbscan(emb2_pam50, minPts = minPts2_pam50)
  her2_high_idx_pam50 <- non_basal_pam50[which(hdb2_pam50$cluster == which.max(table(hdb2_pam50$cluster[hdb2_pam50$cluster != 0])))]
  clusters_hdb_pam50[her2_high_idx_pam50] <- "HER2-high"
}

non_her2_high_pam50 <- which(clusters_hdb_pam50 == "Unassigned")
if (length(non_her2_high_pam50) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification (PAM50).\n", length(non_her2_high_pam50))
} else {
  emb3_pam50 <- uwot::umap(X_pam50[non_her2_high_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  minPts3_pam50 <- choose_minPts(nrow(emb3_pam50))
  cat(sprintf("[HDBSCAN PAM50] Using minPts = %d for n = %d non-HER2-high samples\n", minPts3_pam50, nrow(emb3_pam50)))
  hdb3_pam50 <- hdbscan(emb3_pam50, minPts = minPts3_pam50)
  her2_low_idx_pam50 <- non_her2_high_pam50[which(hdb3_pam50$cluster == which.max(table(hdb3_pam50$cluster[hdb3_pam50$cluster != 0])))]
  clusters_hdb_pam50[her2_low_idx_pam50] <- "HER2-low"
}

luminal_idx_pam50 <- which(clusters_hdb_pam50 == "Unassigned")
if (length(luminal_idx_pam50) < 2) {
  cat("[WARN] Too few Luminal samples (%d); assigning all as LumA (PAM50).\n", length(luminal_idx_pam50))
  clusters_hdb_pam50[luminal_idx_pam50] <- "LumA"
} else {
  emb_lum_pam50 <- uwot::umap(X_pam50[luminal_idx_pam50, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "cosine")
  km_lum_pam50 <- kmeans(emb_lum_pam50, centers = 2, nstart = 50)
  lum_clusters_pam50 <- ifelse(km_lum_pam50$cluster == 1, "LumA", "LumB")
  clusters_hdb_pam50[luminal_idx_pam50] <- lum_clusters_pam50
}

clusters_hdb_pam50 <- factor(clusters_hdb_pam50, levels = subtype_labels)

# --- Apply to HDBSCAN on PAM50 ---
clusters_hdb_pam50_num <- relabel_hdb_numeric(clusters_hdb_pam50, keep_unassigned = TRUE)

hdb_pam50_df <- tibble(sample = rownames(X_pam50),
                       dataset = dataset_lab,
                       subtype = clusters_hdb_pam50_num)
write.csv(hdb_pam50_df, file.path(outdir, "clusters_hdbscan_PAM50.csv"), row.names = FALSE)

emb_hdb_pam50 <- as.data.frame(emb1_pam50); colnames(emb_hdb_pam50) <- c("UMAP1","UMAP2")
emb_hdb_pam50$sample <- rownames(X_pam50)
emb_hdb_pam50$dataset <- dataset_lab
emb_hdb_pam50$subtype <- clusters_hdb_pam50_num

p_hdb_pam50_umap <- ggplot(emb_hdb_pam50, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (HDBSCAN PAM50 clusters: 1–5, NA=unassigned)")

ggsave(file.path(outdir, "umap_hdbscan_PAM50.pdf"), p_hdb_pam50_umap, width = 7, height = 6)
write.csv(emb_hdb_pam50, file.path(outdir, "umap_hdbscan_PAM50_coords.csv"), row.names = FALSE)

# Load libraries needed for plotting (if not already loaded)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("[WARN] ggplot2 package not found. Installing...\n")
  install.packages("ggplot2")
}
library(ggplot2)

# ---------- (7) STATISTICAL ANALYSIS ----------
cat("\n---------- (7) STATISTICAL ANALYSIS ----------\n")
stats_dir <- file.path(outdir, "stats_results")
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Statistical results will be saved to: %s\n", stats_dir))

fisher_results <- list()

# --- Fisher's Exact Test: Cluster vs. Dataset ---
cat("\n[INFO] --- Fisher's Exact Tests: Cluster vs. Dataset --- \n")

# K-means vs. Dataset
cont_table_kmeans <- table(clusters_kmeans, dataset_lab)
cat(sprintf("[INFO] Testing k-means vs. dataset. Contingency table (N=%d):\n", sum(cont_table_kmeans)))
cat(paste("    ", capture.output(print(cont_table_kmeans)), collapse = "\n"), "\n")
if (all(dim(cont_table_kmeans) >= 2) && sum(cont_table_kmeans) >= 5) {
  fisher_kmeans <- fisher.test(cont_table_kmeans, simulate.p.value = nrow(cont_table_kmeans) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_kmeans$p.value))
  fisher_results[["kmeans_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_kmeans_vs_dataset",
    p_value = fisher_kmeans$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# Hierarchical (static) vs. Dataset
cont_table_hc <- table(clusters_hc, dataset_lab)
cat(sprintf("[INFO] Testing Hierarchical (static, HVG) vs. dataset. Contingency table (N=%d):\n", sum(cont_table_hc)))
cat(paste("    ", capture.output(print(cont_table_hc)), collapse = "\n"), "\n")
if (all(dim(cont_table_hc) >= 2) && sum(cont_table_hc) >= 5) {
  fisher_hc <- fisher.test(cont_table_hc, simulate.p.value = nrow(cont_table_hc) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hc$p.value))
  fisher_results[["hc_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hc_vs_dataset",
    p_value = fisher_hc$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# Hierarchical (dynamic, HVG) vs. Dataset
cont_table_hc_dyn <- table(clusters_hc_dyn_named, dataset_lab)
cat(sprintf("[INFO] Testing Hierarchical (dynamic, HVG) vs. dataset. Contingency table (N=%d):\n", sum(cont_table_hc_dyn)))
cat(paste("    ", capture.output(print(cont_table_hc_dyn)), collapse = "\n"), "\n")
if (all(dim(cont_table_hc_dyn) >= 2) && sum(cont_table_hc_dyn) >= 5) {
  fisher_hc_dyn <- fisher.test(cont_table_hc_dyn, simulate.p.value = nrow(cont_table_hc_dyn) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hc_dyn$p.value))
  fisher_results[["hc_dynamic_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hc_dynamic_vs_dataset",
    p_value = fisher_hc_dyn$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# Hierarchical (dynamic, PAM50) vs. Dataset
cont_table_hc_p50_dyn <- table(clusters_hc_pam50_dyn_named, dataset_lab)
cat(sprintf("[INFO] Testing Hierarchical (dynamic, PAM50) vs. dataset. Contingency table (N=%d):\n", sum(cont_table_hc_p50_dyn)))
cat(paste("    ", capture.output(print(cont_table_hc_p50_dyn)), collapse = "\n"), "\n")
if (all(dim(cont_table_hc_p50_dyn) >= 2) && sum(cont_table_hc_p50_dyn) >= 5) {
  fisher_hc_p50_dyn <- fisher.test(cont_table_hc_p50_dyn, simulate.p.value = nrow(cont_table_hc_p50_dyn) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hc_p50_dyn$p.value))
  fisher_results[["hc_pam50_dynamic_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hc_pam50_dynamic_vs_dataset",
    p_value = fisher_hc_p50_dyn$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# HDBSCAN (HVG) vs. Dataset
cont_table_hdb <- table(clusters_hdb_num, dataset_lab)
cat(sprintf("[INFO] Testing HDBSCAN (HVG) vs. dataset. Contingency table (N=%d):\n", sum(cont_table_hdb)))
cat(paste("    ", capture.output(print(cont_table_hdb)), collapse = "\n"), "\n")
if (all(dim(cont_table_hdb) >= 2) && sum(cont_table_hdb) >= 5) {
  fisher_hdb <- fisher.test(cont_table_hdb, simulate.p.value = nrow(cont_table_hdb) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hdb$p.value))
  fisher_results[["hdb_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_vs_dataset",
    p_value = fisher_hdb$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# HDBSCAN (PAM50) vs. Dataset
cont_table_hdb_pam50 <- table(clusters_hdb_pam50_num, dataset_lab)
cat(sprintf("[INFO] Testing HDBSCAN (PAM50) vs. dataset. Contingency table (N=%d):\n", sum(cont_table_hdb_pam50)))
cat(paste("    ", capture.output(print(cont_table_hdb_pam50)), collapse = "\n"), "\n")
if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
  fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2)
  cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hdb_pam50$p.value))
  fisher_results[["hdb_pam50_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_pam50_vs_dataset",
    p_value = fisher_hdb_pam50$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
}

# --- Fisher's Exact Test: Cluster vs. TCGA PAM50 Subtype ---
cat("\n[INFO] --- Fisher's Exact Tests: Cluster vs. TCGA PAM50 Subtype --- \n")
if (!is.na(tcga_sub_col)) {
  tcga_samples <- colnames(V_tcga)
  lab_tcga <- as.character(colData(tcga_se)[tcga_samples, tcga_sub_col])
  
  # K-means vs. PAM50
  cont_table_kmeans_pam50 <- table(clusters_kmeans[tcga_samples], lab_tcga)
  cat(sprintf("[INFO] Testing k-means vs. PAM50 (on TCGA samples). Contingency table (N=%d):\n", sum(cont_table_kmeans_pam50)))
  cat(paste("    ", capture.output(print(cont_table_kmeans_pam50)), collapse = "\n"), "\n")
  if (all(dim(cont_table_kmeans_pam50) >= 2) && sum(cont_table_kmeans_pam50) >= 5) {
    fisher_kmeans_pam50 <- fisher.test(cont_table_kmeans_pam50, simulate.p.value = nrow(cont_table_kmeans_pam50) > 2)
    cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_kmeans_pam50$p.value))
    fisher_results[["kmeans_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_kmeans_vs_PAM50",
      p_value = fisher_kmeans_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
  }
  
  # Hierarchical (static) vs. PAM50
  cont_table_hc_pam50 <- table(clusters_hc[tcga_samples], lab_tcga)
  cat(sprintf("[INFO] Testing Hierarchical (static) vs. PAM50 (on TCGA samples). Contingency table (N=%d):\n", sum(cont_table_hc_pam50)))
  cat(paste("    ", capture.output(print(cont_table_hc_pam50)), collapse = "\n"), "\n")
  if (all(dim(cont_table_hc_pam50) >= 2) && sum(cont_table_hc_pam50) >= 5) {
    fisher_hc_pam50 <- fisher.test(cont_table_hc_pam50, simulate.p.value = nrow(cont_table_hc_pam50) > 2)
    cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hc_pam50$p.value))
    fisher_results[["hc_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hc_vs_PAM50",
      p_value = fisher_hc_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
  }
  
  # HDBSCAN (HVG) vs. PAM50
  cont_table_hdb_pam50 <- table(clusters_hdb_num[tcga_samples], lab_tcga)
  cat(sprintf("[INFO] Testing HDBSCAN (HVG) vs. PAM50 (on TCGA samples). Contingency table (N=%d):\n", sum(cont_table_hdb_pam50)))
  cat(paste("    ", capture.output(print(cont_table_hdb_pam50)), collapse = "\n"), "\n")
  if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
    fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2)
    cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hdb_pam50$p.value))
    fisher_results[["hdb_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_vs_PAM50",
      p_value = fisher_hdb_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
  }
  
  # HDBSCAN (PAM50) vs. PAM50
  cont_table_hdb_pam50_sub <- table(clusters_hdb_pam50_num[tcga_samples], lab_tcga)
  cat(sprintf("[INFO] Testing HDBSCAN (PAM50) vs. PAM50 (on TCGA samples). Contingency table (N=%d):\n", sum(cont_table_hdb_pam50_sub)))
  cat(paste("    ", capture.output(print(cont_table_hdb_pam50_sub)), collapse = "\n"), "\n")
  if (all(dim(cont_table_hdb_pam50_sub) >= 2) && sum(cont_table_hdb_pam50_sub) >= 5) {
    fisher_hdb_pam50_sub <- fisher.test(cont_table_hdb_pam50_sub, simulate.p.value = nrow(cont_table_hdb_pam50_sub) > 2)
    cat(sprintf("[INFO] ... p-value = %.3e\n", fisher_hdb_pam50_sub$p.value))
    fisher_results[["hdb_pam50_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_pam50_vs_PAM50",
      p_value = fisher_hdb_pam50_sub$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] ... Skipping test: table dims < 2 or total N < 5.\n")
  }
} else {
  cat("[WARN] No TCGA subtype column found. Skipping all Cluster vs. PAM50 tests.\n")
}

fisher_df <- do.call(rbind, fisher_results)
out_file_fisher <- file.path(stats_dir, "fisher_tests.csv")
write.csv(fisher_df, out_file_fisher, row.names = FALSE)
cat(sprintf("\n[INFO] Saved %d Fisher's Exact Test results to: %s\n", nrow(fisher_df), out_file_fisher))


# --- Kruskal-Wallis and DSCF Tests: Marker gene expression across subtypes ---
cat("\n[INFO] --- Kruskal-Wallis & DSCF Tests: Marker Expression vs. Clusters --- \n")
kw_results <- list()
dscf_results <- list()
methods_to_test <- c("kmeans", "hc", "hdb", "hdb_pam50")
cat(sprintf("[INFO] Testing methods: %s\n", paste(methods_to_test, collapse = ", ")))

for (method in methods_to_test) {
  cat(sprintf("\n[INFO] == Starting Method: %s ==\n", method))
  clusters <- switch(method,
                     kmeans = clusters_kmeans,
                     hc = clusters_hc,
                     hdb = clusters_hdb_num,
                     hdb_pam50 = clusters_hdb_pam50_num)
  mat <- if (method == "hdb_pam50") M_pam50 else V_adj
  genes <- if (method == "hdb_pam50") rownames(M_pam50) else marker_ens
  
  if (length(genes) == 0) {
    cat(sprintf("[WARN] No marker genes found for method '%s'. Skipping all KW tests.\n", method))
    next
  }

  for (gene in genes) {
    gene_sym <- if (method == "hdb_pam50") gene else ens2sym_union$symbol[match(gene, ens2sym_union$ensembl)]
    gene_sym <- ifelse(is.na(gene_sym) || gene_sym == "", gene, gene_sym)
    
    cat(sprintf("[INFO] Testing Gene: %s (Method: %s)\n", gene_sym, method))
    
    if (!gene %in% rownames(mat)) {
      cat(sprintf("[WARN] ... Gene %s (%s) not found in expression matrix. Skipping.\n", gene_sym, gene))
      next
    }
    
    expr <- mat[gene, ]
    if (sum(!is.na(expr)) < length(expr) * 0.5) {
      cat(sprintf("[WARN] ... Skipping %s: >50%% NA values.\n", gene_sym))
      next
    }
    
    # Ensure there are at least 2 groups with data to compare
    valid_data <- data.frame(expr = expr, clusters = clusters)[!is.na(expr), ]
    valid_groups <- length(unique(valid_data$clusters))
    
    if (valid_groups < 2) {
      cat(sprintf("[WARN] ... Skipping %s: Only %d group(s) with valid data.\n", gene_sym, valid_groups))
      next
    }
    
    kw_test <- kruskal.test(expr ~ clusters)
    if (is.na(kw_test$p.value)) {
      cat(sprintf("[WARN] ... KW test failed for %s (p=NA), likely insufficient variation.\n", gene_sym))
      next
    }
    
    cat(sprintf("[INFO] ... KW Test: statistic = %.3f, p-value = %.3e\n", kw_test$statistic, kw_test$p.value))
    kw_results[[paste(method, gene_sym, sep = "_")]] <- data.frame(
      test = sprintf("Kruskal_Wallis_%s_%s", method, gene_sym),
      gene = gene,
      symbol = gene_sym,
      p_value = kw_test$p.value,
      statistic = kw_test$statistic
    )
    
    # --- Post-hoc DSCF and Plotting ---
    if (kw_test$p.value < 0.05) {
      cat("[INFO] ... p < 0.05, running post-hoc DSCF test...\n")
      
      # work on samples with non-NA expression only
      keep <- !is.na(expr)
      expr_ok <- expr[keep]
      grp_ok  <- droplevels(factor(clusters[keep])) # Use factor() to ensure droplevels works

      # need at least two non-empty groups
      if (nlevels(grp_ok) >= 2) {
        cat(sprintf("[INFO] ... DSCF: testing on %d samples across %d groups.\n", length(expr_ok), nlevels(grp_ok)))
        dscf_df <- NULL
        # run DSCF on present groups only
        dscf_test <- try(PMCMRplus::dscfAllPairsTest(expr_ok, grp_ok, p.adjust.method = "BH"), silent = TRUE)

        if (!inherits(dscf_test, "try-error") && !is.null(dscf_test$p.value)) {
          cat("[INFO] ... DSCF test successful. Extracting pairs...\n")
          p_mat <- try(as.matrix(dscf_test$p.value), silent = TRUE)
          if (inherits(p_mat, "try-error") || length(dim(p_mat)) != 2) next
          if (nrow(p_mat) < 1 || ncol(p_mat) < 1) next # Handle 1x1 or empty matrix

          # Ensure usable dimnames. Use the actually-present group levels in grp_ok.
          lev_present <- levels(grp_ok)
          if(is.null(rownames(p_mat)) || all(rownames(p_mat) == "")) rownames(p_mat) <- lev_present
          if(is.null(colnames(p_mat)) || all(colnames(p_mat) == "")) colnames(p_mat) <- lev_present

          # Prefer upper triangle; if empty, fall back to lower triangle.
          pairs_idx <- which(upper.tri(p_mat, diag = FALSE), arr.ind = TRUE)
          swap <- FALSE
          if (nrow(pairs_idx) == 0) {
            pairs_idx <- which(lower.tri(p_mat, diag = FALSE), arr.ind = TRUE)
            swap <- TRUE
          }
          if (nrow(pairs_idx) == 0) next  # nothing to extract

          # Index directly by position
          subtype1 <- rownames(p_mat)[pairs_idx[, 1]]
          subtype2 <- colnames(p_mat)[pairs_idx[, 2]]
          pvals    <- p_mat[cbind(pairs_idx[, 1], pairs_idx[, 2])]

          if (swap) {
            tmp <- subtype1; subtype1 <- subtype2; subtype2 <- tmp
          }

          keep <- is.finite(pvals) & !is.na(subtype1) & !is.na(subtype2)
          if (!any(keep)) next

          dscf_df <- data.frame(
            test     = sprintf("DSCF_%s_%s_%s_vs_%s", method, gene_sym, subtype1[keep], subtype2[keep]),
            gene     = gene,
            symbol   = gene_sym,
            subtype1 = subtype1[keep],
            subtype2 = subtype2[keep],
            p_value  = pvals[keep],
            stringsAsFactors = FALSE
          )

          if (is.data.frame(dscf_df) && nrow(dscf_df) > 0) {
            cat(sprintf("[INFO] ... Found %d significant DSCF pairwise comparisons.\n", nrow(dscf_df)))
            dscf_results[[paste(method, gene_sym, sep = "_")]] <- dscf_df
          }
        } else {
           cat(sprintf("[WARN] ... DSCF test failed for %s.\n", gene_sym))
        }
      } else {
        cat(sprintf("[WARN] ... Skipping DSCF for %s: < 2 groups with valid data.\n", gene_sym))
      }
      
      # --- Plotting (Improved) ---
      cat(sprintf("[INFO] ... Saving violin plot for significant gene: %s\n", gene_sym))
      
      plot_data <- data.frame(expr = mat[gene, ], subtype = clusters)
      plot_data <- plot_data[!is.na(plot_data$expr), ] # Remove NAs for plotting
      plot_data$subtype <- factor(plot_data$subtype) # Ensure it's a factor

      # Calculate annotation position
      y_pos <- max(plot_data$expr, na.rm = TRUE) + (0.05 * diff(range(plot_data$expr, na.rm = TRUE)))
      x_pos <- median(1:nlevels(plot_data$subtype)) # Center the text
      
      p_box <- ggplot(plot_data, aes(x = subtype, y = expr, fill = subtype)) +
        geom_violin(alpha = 0.6, trim = FALSE, draw_quantiles = c(0.5)) + # Violin with median
        geom_jitter(width = 0.15, alpha = 0.3, size = 0.7, height = 0) +  # Add raw data points
        theme_bw() +
        labs(title = sprintf("Expression of %s across %s clusters", gene_sym, method),
             y = "VST Expression", x = "Cluster") +
        annotate(geom = "text",
                 x = x_pos,
                 y = y_pos,
                 label = sprintf("Kruskal-Wallis p = %.3e", kw_test$p.value),
                 size = 4, fontface = "italic") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") # Remove redundant legend
              
      out_plot_file <- file.path(stats_dir, sprintf("violin_plot_%s_%s.pdf", method, gene_sym))
      ggsave(out_plot_file, p_box, width = 7, height = 5)
    }
  }
}

# --- Save Final Statistical Results ---
kw_df <- tryCatch(do.call(rbind, kw_results), error = function(e) NULL)
out_file_kw <- file.path(stats_dir, "kruskal_wallis_tests.csv")
if (!is.null(kw_df) && nrow(kw_df) > 0) {
  write.csv(kw_df, out_file_kw, row.names = FALSE)
  cat(sprintf("\n[INFO] Saved %d Kruskal-Wallis results to: %s\n", nrow(kw_df), out_file_kw))
} else {
  cat("\n[INFO] No Kruskal-Wallis results to write.\n")
}

dscf_df <- tryCatch(do.call(rbind, dscf_results), error = function(e) NULL)
out_file_dscf <- file.path(stats_dir, "dscf_tests.csv")
if (!is.null(dscf_df) && nrow(dscf_df) > 0) {
  write.csv(dscf_df, out_file_dscf, row.names = FALSE)
  cat(sprintf("[INFO] Saved %d DSCF pairwise results to: %s\n", nrow(dscf_df), out_file_dscf))
} else {
  cat("[INFO] No pairwise DSCF results to write.\n")
}

cat("\n[INFO] Statistical analysis complete.\n")

# ---------- (8) MARKER HEATMAPS & CORRELATIONS ----------
cat("\n---------- (8) MARKER HEATMAPS & CORRELATIONS ----------\n")
cat("[INFO] Mapping gene symbols for marker lists...\n")

ens2sym <- make_ens2sym(dsmz_raw)
ens_vec <- strip_ensver(as_df(dsmz_raw)$Ensembl_ID)
dup_ids <- sum(duplicated(ens_vec))
dup_unique <- sum(table(ens_vec) > 1, na.rm = TRUE)
cat(sprintf("[INFO] DSMZ Ensembl duplicate rows: %d | duplicate IDs: %d (collapsed)\n",
            dup_ids, dup_unique))

marker_ens <- rownames(ens2sym)[tolower(ens2sym$symbol) %in% tolower(marker_symbols)]
marker_ens <- intersect(unique(marker_ens), rownames(V_adj))
cat(sprintf("[INFO] Found %d marker genes present in the V_adj matrix.\n", length(marker_ens)))


# --- Generating PAM50 Gene Heatmap (All Samples, ordered) ---
cat("\n[INFO] --- Generating PAM50 Gene Heatmap (All Samples, ordered) --- \n")
if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50
  cat(sprintf("[INFO] Using PAM50 matrix with %d genes and %d samples.\n", nrow(V_mark), ncol(V_mark)))
  if (any(duplicated(rownames(V_mark)))) {
    cat("[INFO] ... Making duplicated gene row names unique.\n")
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  
  cat("[INFO] ... Building annotation table with all cluster assignments.\n")
  ann <- data.frame(dataset = dataset_lab, 
                    kmeans_cluster = clusters_kmeans,
                    hc_subtype = clusters_hc,
                    hc_dyn_subtype = clusters_hc_dyn_named,
                    hc_pam50_dyn_subtype = clusters_hc_pam50_dyn_named,
                    hdb_subtype = clusters_hdb_num,
                    hdb_pam50_subtype = clusters_hdb_pam50_num,
                    row.names = colnames(V_mark))
                    
  out_pdf_pam50_all <- file.path(outdir, "heatmap_pam50_all_clusters.pdf")
  cat(sprintf("[INFO] Saving PAM50 heatmap (ordered by HDBSCAN-PAM50 clusters) to: %s\n", out_pdf_pam50_all))
  
  safe_pdf(out_pdf_pam50_all, {
    pheatmap(V_mark[, order(clusters_hdb_pam50_num), drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE, # Note: Columns are pre-ordered
             annotation_col = ann,
             main = "PAM50 Gene Expression (k-means, HC, HDBSCAN, HDBSCAN PAM50)")
  })
} else {
  cat("[WARN] No PAM50 genes found. Skipping heatmap_pam50_all_clusters.pdf.\n")
}


# --- Correlating Cluster Centroids vs. TCGA Subtype Centroids ---
cat("\n[INFO] --- Correlating Cluster Centroids vs. TCGA Subtype Centroids --- \n")
if (!is.na(tcga_sub_col)) {
  cat(sprintf("[INFO] Using TCGA annotation column: '%s'\n", tcga_sub_col))
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga)))
  cat(sprintf("[INFO] Calculating TCGA centroids for %d subtypes: %s\n", length(sub_lvls), paste(sub_lvls, collapse = ", ")))
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE))
  
  cat("[INFO] Calculating centroids for k-means, HC, and HDBSCAN clusters...\n")
  cl_means_kmeans <- sapply(levels(clusters_kmeans), function(cn) rowMeans(V_adj[hvgs, clusters_kmeans == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hc <- sapply(levels(clusters_hc), function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hdb <- sapply(levels(clusters_hdb_num), function(cn) rowMeans(V_adj[hvgs, clusters_hdb_num == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hdb_pam50 <- sapply(levels(clusters_hdb_pam50_num), function(cn) rowMeans(M_pam50[, clusters_hdb_pam50_num == cn, drop = FALSE], na.rm = TRUE))
  
  g_common <- intersect(rownames(sub_means), hvgs)
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50))
  cat(sprintf("[INFO] ... Using %d common HVGs for HVG-based cluster correlations.\n", length(g_common)))
  cat(sprintf("[INFO] ... Using %d common PAM50 genes for PAM50-based cluster correlations.\n", length(g_common_pam50)))

  corr_kmeans <- cor(cl_means_kmeans[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hc <- cor(cl_means_hc[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hdb <- cor(cl_means_hdb[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hdb_pam50 <- cor(cl_means_hdb_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman")
  
  cat("[INFO] Saving correlation matrices (cluster vs. TCGA) to CSV files...\n")
  write.csv(corr_kmeans, file.path(outdir, "cor_kmeans_vs_TCGA_subtype.csv"))
  write.csv(corr_hc, file.path(outdir, "cor_hierarchical_vs_TCGA_subtype.csv"))
  write.csv(corr_hdb, file.path(outdir, "cor_hdbscan_vs_TCGA_subtype.csv"))
  write.csv(corr_hdb_pam50, file.path(outdir, "cor_hdbscan_pam50_vs_TCGA_subtype.csv"))
  
  out_pdf_corr_hm <- file.path(outdir, "heatmap_cluster_vs_TCGA_subtype_corr.pdf")
  cat(sprintf("[INFO] Saving correlation heatmaps to: %s\n", out_pdf_corr_hm))
  safe_pdf(out_pdf_corr_hm, {
    par(mfrow = c(2, 2))
    pheatmap::pheatmap(corr_kmeans, main = "k-means vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hc, main = "Hierarchical vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hdb, main = "HDBSCAN vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hdb_pam50, main = "HDBSCAN PAM50 vs TCGA Subtypes (ρ)")
    par(mfrow = c(1, 1))
  })
} else {
  cat("[WARN] No TCGA subtype column found. Skipping all cluster vs. TCGA correlations.\n")
}


# --- Mapping DSMZ Samples to TCGA (Spearman's Rho) ---
cat("\n[INFO] --- Mapping DSMZ Samples to TCGA (Spearman's Rho) --- \n")
# Overall correlation
cat(sprintf("[INFO] Correlating %d DSMZ samples against the overall TCGA mean expression (all %d TCGA samples).\n", ncol(V_dsmz), ncol(V_tcga)))
tcga_mean <- rowMeans(V_tcga, na.rm = TRUE)
rho_overall <- sapply(colnames(V_dsmz), function(s) {
  suppressWarnings(cor(V_dsmz[, s], tcga_mean, method = "spearman", use = "pairwise.complete.obs"))
})

rho_by_sub <- NULL
best_tcga_subtype <- NULL
best_tcga_rho <- NULL

# By-subtype correlation
if (!is.na(tcga_sub_col)) {
  lab <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_levels <- sort(unique(na.omit(lab)))
  cat(sprintf("[INFO] Correlating %d DSMZ samples against %d TCGA subtype centroids...\n", ncol(V_dsmz), length(sub_levels)))
  sub_means <- sapply(sub_levels, function(s) {
    rowMeans(V_tcga[, lab == s, drop = FALSE], na.rm = TRUE)
  })
  rho_by_sub <- sapply(colnames(V_dsmz), function(s) {
    apply(sub_means, 2, function(ref) {
      suppressWarnings(cor(V_dsmz[, s], ref, method = "spearman", use = "pairwise.complete.obs"))
    })
  })
  colnames(rho_by_sub) <- colnames(V_dsmz)
  best_tcga_subtype <- apply(rho_by_sub, 2, function(v) names(which.max(v)))
  best_tcga_rho <- apply(rho_by_sub, 2, max)
  
  out_csv_dsmz_rho <- file.path(outdir, "rho_dsmz_to_tcga_subtypes.csv")
  cat(sprintf("[INFO] ... Saving DSMZ-to-TCGA subtype correlation matrix to: %s\n", out_csv_dsmz_rho))
  write.csv(rho_by_sub, out_csv_dsmz_rho)
} else {
  cat("[WARN] No TCGA subtype column. Skipping DSMZ-to-TCGA *subtype* correlation.\n")
}

# Create and save summary table
cat("[INFO] Compiling final DSMZ mapping table (clusters + TCGA correlation).\n")
map_tbl <- tibble(
  sample = colnames(V_dsmz),
  rho_overall = unname(rho_overall),
  kmeans_cluster = clusters_kmeans[colnames(V_dsmz)],
  hc_subtype = clusters_hc[colnames(V_dsmz)],
  hdb_subtype = clusters_hdb_num[colnames(V_dsmz)],
  hdb_pam50_subtype = clusters_hdb_pam50_num[colnames(V_dsmz)]
)
if (!is.null(best_tcga_subtype)) {
  map_tbl$best_tcga_subtype <- best_tcga_subtype[map_tbl$sample]
  map_tbl$best_tcga_subtype_rho <- best_tcga_rho[map_tbl$sample]
}
out_csv_dsmz_map <- file.path(outdir, "mapping_dsmz_to_tcga.csv")
write.csv(map_tbl, out_csv_dsmz_map, row.names = FALSE)
cat(sprintf("[INFO] Saved DSMZ mapping table to: %s\n", out_csv_dsmz_map))
cat("[INFO] ... Preview of mapping table:\n")
cat(paste("    ", capture.output(print(head(map_tbl))), collapse = "\n"), "\n")


# --- Generating Joint PAM50 Heatmap (TCGA + DSMZ, clustered) ---
cat("\n[INFO] --- Generating Joint PAM50 Heatmap (TCGA + DSMZ, clustered) --- \n")
if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50
  cat(sprintf("[INFO] Using PAM50 matrix with %d genes and %d total samples (TCGA + DSMZ).\n", nrow(V_mark), ncol(V_mark)))
  if (any(duplicated(rownames(V_mark)))) {
    cat("[INFO] ... Making duplicated gene row names unique.\n")
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  
  # Build TCGA annotations
  tcga_sub_vec <- rep(NA_character_, ncol(V_tcga))
  if (!is.na(tcga_sub_col)) {
    tcga_sub_vec <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  }
  ann_tcga <- data.frame(
    dataset = "TCGA",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_tcga)]),
    hc_subtype = as.character(clusters_hc[colnames(V_tcga)]),
    hc_dyn_subtype = as.character(clusters_hc_dyn_named[colnames(V_tcga)]),
    hc_pam50_dyn_subtype = as.character(clusters_hc_pam50_dyn_named[colnames(V_tcga)]),
    hdb_subtype = as.character(clusters_hdb_num[colnames(V_tcga)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50_num[colnames(V_tcga)]),
    TCGA_subtype = tcga_sub_vec,
    row.names = colnames(V_tcga),
    check.names = FALSE
  )
  cat(sprintf("[INFO] ... Built TCGA annotation table for %d samples.\n", nrow(ann_tcga)))

  # Build DSMZ annotations
  ann_dsmz <- data.frame(
    dataset = "DSMZ",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_dsmz)]),
    hc_subtype = as.character(clusters_hc[colnames(V_dsmz)]),
    hc_dyn_subtype = as.character(clusters_hc_dyn_named[colnames(V_dsmz)]),
    hc_pam50_dyn_subtype = as.character(clusters_hc_pam50_dyn_named[colnames(V_dsmz)]),
    hdb_subtype = as.character(clusters_hdb_num[colnames(V_dsmz)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50_num[colnames(V_dsmz)]),
    TCGA_subtype = NA_character_, # DSMZ samples don't have a ground-truth TCGA subtype
    row.names = colnames(V_dsmz),
    check.names = FALSE
  )
  cat(sprintf("[INFO] ... Built DSMZ annotation table for %d samples.\n", nrow(ann_dsmz)))

  # Combine annotations and data in consistent order
  cols_order <- c(colnames(V_tcga), colnames(V_dsmz))
  V_mark <- V_mark[, cols_order, drop = FALSE]
  ann_col <- rbind(ann_tcga, ann_dsmz)[cols_order, , drop = FALSE]
  
  # Clean up annotations: remove any columns that are all NA
  non_empty_annot <- colSums(!is.na(ann_col)) > 0
  ann_col <- if (any(non_empty_annot)) ann_col[, non_empty_annot, drop = FALSE] else NULL
  cat(sprintf("[INFO] ... Final annotation table has %d columns with data.\n", ifelse(is.null(ann_col), 0, ncol(ann_col))))

  # Plot heatmap
  if (nrow(V_mark) >= 2 && ncol(V_mark) >= 2) {
    out_pdf_joint_hm <- file.path(outdir, "heatmap_pam50_vst.pdf")
    cat(sprintf("[INFO] Saving joint heatmap (all samples, clustered columns) to: %s\n", out_pdf_joint_hm))
    
    pdf(out_pdf_joint_hm, width = 12, height = 8)
    pheatmap( V_mark, show_rownames = TRUE, show_colnames = FALSE,
              cluster_rows = TRUE, cluster_cols = TRUE, # Note: Columns are clustered by pheatmap
              annotation_col = ann_col,
              main = "Marker expression (VST) across TCGA-BRCA & DSMZ-BRCA"
    )
    dev.off()
  } else {
    cat(sprintf("[WARN] Skipping joint heatmap: matrix dim = %d x %d (need >= 2 x 2).\n",
                nrow(V_mark), ncol(V_mark)))
  }
} else {
   cat("[WARN] No PAM50 genes found. Skipping joint heatmap (heatmap_pam50_vst.pdf).\n")
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
  cat("[INFO] Attempting TCGA gene symbol mapping via rowData...\n")
  rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se))
  symcol <- intersect(tolower(colnames(rd)), c("gene_name","symbol","hgnc_symbol"))
  if (length(symcol)) {
    result <- data.frame(
      ensembl = strip_ensver(rownames(rd)),
      symbol  = as.character(rd[[ colnames(rd)[match(symcol[1], tolower(colnames(rd)))] ]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(symbol) & symbol != "") %>%
      distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)
    cat(sprintf("[INFO] TCGA rowData mapping successful: %d mappings using column '%s'\n", 
                nrow(result), symcol[1]))
    result
  } else {
    cat("[WARN] TCGA rowData mapping failed: no suitable symbol column found\n")
    NULL
  }
}, error = function(e) {
  cat(sprintf("[WARN] TCGA rowData mapping failed: %s\n", e$message))
  NULL
})

ens2sym_union <- bind_rows(filter(map_dsmz, !is.na(symbol)),
                           filter(map_tcga, !is.null(map_tcga))) %>%
  distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)

cat(sprintf("[INFO] Marker gene symbol mapping summary: %d total mappings\n", nrow(ens2sym_union)))

# 2) Ensembl markers that exist in V_adj
marker_ens <- ens2sym_union$ensembl[
  tolower(ens2sym_union$symbol) %in% tolower(marker_symbols)
]
marker_ens <- intersect(unique(marker_ens), rownames(V_adj))

cat(sprintf("[MARKER] Found %d of %d PAM50/PAM-like markers in V_adj\n",
            length(marker_ens), length(marker_symbols)))
if (length(marker_ens)) {
  cat("[MARKER] Example matches:", paste(head(marker_ens, 5), collapse=", "), "...\n")
} else {
  stop("[FATAL] No marker genes from the requested list are present in V_adj. Check symbol mapping.")
}

V_tcga_z <- zscore_by_gene(V_tcga)  # Z-score TCGA data
V_dsmz_z <- zscore_by_gene(V_dsmz)  # Z-score DSMZ data

# Helper function to compute median expression per group
aggregate_median <- function(M, groups) {
  grp <- droplevels(factor(groups))
  lev <- levels(grp)

  out <- sapply(lev, function(g) {
    rowMedians(M[, grp == g, drop = FALSE], na.rm = TRUE)
  })

  # ensure 2D and set names
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  rownames(out) <- rownames(M)
  colnames(out) <- lev
  out
}

# Aggregate TCGA by subtype (if available) or overall
if (!is.na(tcga_sub_col)) {  # If subtype column exists
  tcga_groups <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])  # Get subtype labels
  TCGA_med <- aggregate_median(V_tcga_z, tcga_groups)  # Calculate median per subtype
} else {  # If no subtype column
  TCGA_med <- rowMedians(V_tcga_z, na.rm = TRUE)  # Calculate overall median
  TCGA_med <- matrix(TCGA_med, ncol = 1,  # Convert to matrix
                   dimnames = list(rownames(V_tcga_z), "TCGA_BRCA"))  # Set dimension names
}

# Aggregate DSMZ by discovery clusters
DSMZ_med <- aggregate_median(V_dsmz_z, clusters_kmeans[colnames(V_dsmz)])  # Calculate median per cluster

# Extra safety: ensure rownames are preserved
if (is.null(rownames(TCGA_med)))  rownames(TCGA_med)  <- rownames(V_tcga_z)
if (is.null(rownames(DSMZ_med))) rownames(DSMZ_med) <- rownames(V_dsmz_z)

## ---- safe median-z heatmap ----
common_mark <- intersect(marker_ens, intersect(rownames(TCGA_med), rownames(DSMZ_med)))
if (length(common_mark) == 0L) {
  cat("[DEBUG] marker_ens length =", length(marker_ens), "\n")
  cat("[DEBUG] nrow(TCGA_med) =", nrow(TCGA_med), " nrow(DSMZ_med) =", nrow(DSMZ_med), "\n")
  cat("[DEBUG] Example TCGA_med genes:", paste(head(rownames(TCGA_med), 5), collapse=", "), "\n")
  cat("[DEBUG] Example DSMZ_med genes:", paste(head(rownames(DSMZ_med), 5), collapse=", "), "\n")
  stop("[FATAL] No overlap between marker and (TCGA_med ∩ DSMZ_med). Check ID harmonization (Ensembl, version stripping) and markers.")
} else {
  H <- cbind(TCGA_med[common_mark, , drop = FALSE], DSMZ_med[common_mark, , drop = FALSE])
  pdf(file.path(outdir, "heatmap_marker_median_z.pdf"), width = 10, height = 7)
  pheatmap(H, show_rownames = TRUE,
           main = "Median z-scores of marker genes by subtype/cluster")
  dev.off()
}

# ---------- (5) TCGA REFERENCE BUILDING AND DSMZ ASSIGNMENT ----------
cat("\n---------- (5) TCGA REFERENCE BUILDING AND DSMZ ASSIGNMENT ----------\n")

## =====================================================================
## (A) Build balanced TCGA reference (PAM50 gene space, k=20/subtype)
## =====================================================================
cat("\n## =====================================================================\n")
cat("## (A) Build balanced TCGA reference (PAM50 gene space, k=20/subtype)\n")
cat("## =====================================================================\n")
set.seed(42)

# PAM50 matrix for TCGA samples (genes x samples)
cat(sprintf("[A.1] Subsetting PAM50 matrix for %d TCGA samples.\n", ncol(V_tcga)))
M_pam50_tcga <- M_pam50[, colnames(V_tcga), drop = FALSE]

# Scale samples in PAM50 space for clustering/nearest calculations
X_tcga <- t(scale(M_pam50_tcga))  # samples x genes (scaled)
cat(sprintf("[A.1] Scaling TCGA data. Input: %d genes x %d samples. Output: %d samples x %d genes.\n",
            nrow(M_pam50_tcga), ncol(M_pam50_tcga), nrow(X_tcga), ncol(X_tcga)))

if (nrow(X_tcga) < 5) {
  stop(sprintf("[FATAL] Too few TCGA samples in PAM50 space (%d) for clustering.", nrow(X_tcga)))
}

# Hierarchical clustering (k=5)
cat("[A.2] Performing hierarchical clustering (k=5, euclidean, average) on scaled TCGA samples...\n")
d_tcga  <- dist(X_tcga, method = "euclidean")
hc_tcga <- hclust(d_tcga, method = "average")
cl_tcga <- cutree(hc_tcga, k = 5)

# Find a TCGA subtype column if present
tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]
subtype_labels <- c("HER2-high","HER2-low","LumA","LumB","Basal")

# Map HC clusters to PAM50 names using correlation to TCGA subtype centroids (if annotation exists)
cat("[A.3] Mapping k=5 HC clusters to PAM50 subtype names.\n")
if (!is.na(tcga_sub_col)) {
  cat(sprintf("[A.3] Found TCGA annotation column: '%s'. Starting centroid correlation...\n", tcga_sub_col))
  lab_tcga_full <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga_full)))
  cat(sprintf("[A.3] Found %d TCGA subtypes: %s\n", length(sub_lvls), paste(sub_lvls, collapse = ", ")))
  
  # centroids in PAM50 space (genes as rows)
  sub_means_tcga <- sapply(sub_lvls, function(s) {
    rowMeans(M_pam50_tcga[, lab_tcga_full == s, drop = FALSE], na.rm = TRUE)
  })
  cat(sprintf("[A.3] Calculated TCGA centroids: %d genes x %d subtypes.\n", nrow(sub_means_tcga), ncol(sub_means_tcga)))

  cl_means_tcga <- sapply(1:5, function(k) {
    rowMeans(M_pam50_tcga[, cl_tcga == k, drop = FALSE], na.rm = TRUE)
  })
  cat(sprintf("[A.3] Calculated k=5 HC cluster centroids: %d genes x %d clusters.\n", nrow(cl_means_tcga), ncol(cl_means_tcga)))

  g_common <- intersect(rownames(sub_means_tcga), rownames(cl_means_tcga))
  cat(sprintf("[A.3] Correlating centroids using %d common PAM50 genes.\n", length(g_common)))
  corr_mat <- cor(cl_means_tcga[g_common, , drop = FALSE],
                  sub_means_tcga[g_common, , drop = FALSE],
                  method = "spearman")

  cl_to_sub <- apply(corr_mat, 1, function(v) sub_lvls[which.max(v)])
  cl_to_sub_named <- setNames(subtype_labels[match(cl_to_sub, sub_lvls)], 1:5)
  
  # [NEW] Log the mapping
  map_df <- data.frame(Cluster = names(cl_to_sub_named), Mapped_Subtype = cl_to_sub_named)
  cat("[A.3] HC Cluster to TCGA Subtype Mapping (via correlation):\n")
  out_map <- capture.output(print(map_df, row.names=FALSE))
  cat(paste("    ", out_map, collapse = "\n"), "\n")
  
} else {
  # If no annotation exists, give stable names via simple ordering heuristic
  cat(sprintf("[A.3] [WARN] No TCGA annotation column found. Using k-means heuristic to label clusters.\n"))
  km5 <- kmeans(X_tcga, centers = 5, nstart = 50)
  cat("[A.3] Fallback: Ran k-means (k=5) to get stable cluster order.\n")
  km_cent <- t(sapply(1:5, function(k) colMeans(X_tcga[km5$cluster == k, , drop = FALSE], na.rm = TRUE)))
  ord <- order(prcomp(km_cent)$x[,1])
  cl_to_sub_named <- setNames(subtype_labels[ord], 1:5)

  # [NEW] Log the fallback mapping
  map_df <- data.frame(Cluster = names(cl_to_sub_named), Mapped_Subtype = cl_to_sub_named)
  cat("[A.3] Fallback HC Cluster to Subtype Mapping (via k-means PCA order):\n")
  out_map <- capture.output(print(map_df, row.names=FALSE))
  cat(paste("    ", out_map, collapse = "\n"), "\n")
}

# Final subtype per TCGA sample (factor with fixed order)
subtype_tcga <- factor(cl_to_sub_named[as.character(cl_tcga)], levels = subtype_labels)
cat("[A.4] Final TCGA sample subtype assignments (total counts):\n")
out_table <- capture.output(print(table(subtype_tcga)))
cat(paste("    ", out_table, collapse = "\n"), "\n")

# ---- Helpers for centroid + nearest selection in scaled space ----
get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
get_k_nearest <- function(X, centroid, k = 20) {
  d <- sqrt(rowSums((X - matrix(centroid, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
  ord <- order(d, decreasing = FALSE)
  keep <- ord[seq_len(min(k, length(ord)))]
  data.frame(sample = rownames(X)[keep], dist = d[keep], stringsAsFactors = FALSE)
}

# ---- Select up to top-20 per subtype (closest to the subtype centroid) ----
cat("[A.5] Selecting up to 20 samples per subtype (nearest to centroid in scaled PAM50 space)...\n")
sel_list <- list(); summary_counts <- c()
for (sub in subtype_labels) {
  idx <- which(subtype_tcga == sub)
  if (length(idx) == 0) {
    cat(sprintf("[WARN] No TCGA samples for subtype '%s' - skipping.\n", sub))
    next
  }
  cen <- get_centroid(X_tcga, idx)                           # centroid in *scaled* space
  nearest <- get_k_nearest(X_tcga[idx, , drop = FALSE], cen, k = min(20, length(idx)))
  nearest$subtype <- sub
  sel_list[[sub]] <- nearest
  summary_counts[sub] <- nrow(nearest)
}
sel_df <- dplyr::bind_rows(sel_list) %>% dplyr::relocate(subtype, .before = sample)

# Enforce ≤20 per subtype (guard-rail in case anything changes upstream)
cat("[A.5] (Applying guard-rail: ensuring n <= 20 per subtype).\n")
sel_df <- sel_df %>%
  group_by(subtype) %>% 
  slice_head(n = 20) %>% 
  ungroup()

# [NEW] Log the final selection counts
cat("[A.5] Reference selection complete. Final counts per subtype:\n")
out_table_sel <- capture.output(print(table(sel_df$subtype)))
cat(paste("    ", out_table_sel, collapse = "\n"), "\n")

# Assert the counts and uniqueness before saving
stopifnot(all(table(sel_df$subtype) <= 20))
stopifnot(length(unique(sel_df$sample)) == nrow(sel_df))

# Save selection
ref_dir <- file.path(outdir, "tcga_reference_selection")
dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)

out_csv_ref <- file.path(ref_dir, "tcga_top20_per_subtype.csv")
cat(sprintf("[A.6] Saving reference sample list (%d total samples) to: %s\n", nrow(sel_df), out_csv_ref))
write.csv(sel_df, out_csv_ref, row.names = FALSE)

# Save PAM50-only expression matrix for these references (genes x selected samples)
ref_samples <- sel_df$sample
M_pam50_tcga_ref <- M_pam50_tcga[, ref_samples, drop = FALSE]
out_rds_ref <- file.path(ref_dir, "PAM50_tcga_top20_per_subtype.rds")
cat(sprintf("[A.6] Saving PAM50 expression matrix for reference samples (%d genes x %d samples) to: %s\n",
            nrow(M_pam50_tcga_ref), ncol(M_pam50_tcga_ref), out_rds_ref))
saveRDS(M_pam50_tcga_ref, out_rds_ref)

# Quick UMAP to visualize chosen references
cat("[A.7] Generating UMAP visualization of TCGA reference selection (using 'cosine' metric on PAM50 space).\n")
set.seed(42)
emb_ref <- uwot::umap(X_tcga, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb_ref <- as.data.frame(emb_ref); colnames(emb_ref) <- c("UMAP1","UMAP2")
emb_ref$sample  <- rownames(X_tcga)
emb_ref$subtype <- subtype_tcga
emb_ref$chosen  <- ifelse(emb_ref$sample %in% ref_samples, "Top20", "Other")

# Keep umap on all samples (for better geometry), but plot only the refs:
emb_ref_plot <- emb_ref[emb_ref$sample %in% ref_samples, ]

p_ref <- ggplot(emb_ref_plot, aes(UMAP1, UMAP2, color = subtype)) +
  geom_point(alpha = 0.95, size = 2) + theme_bw() +
  ggtitle("TCGA BRCA (PAM50) - selected Top20 per subtype only")
out_pdf_umap_ref <- file.path(ref_dir, "umap_tcga_reference_top20.pdf")
ggsave(out_pdf_umap_ref, p_ref, width = 7, height = 6)
cat(sprintf("[A.7] Saved UMAP plot to: %s\n", out_pdf_umap_ref))

# Export plain list by subtype
dumpfile <- file.path(ref_dir, "tcga_top20_samples_by_subtype.txt")
cat(sprintf("[A.7] Writing plain-text sample list to: %s\n", dumpfile))
con_ref <- file(dumpfile, "wt"); on.exit(close(con_ref), add = TRUE)
by_sub <- split(sel_df$sample, sel_df$subtype)
for (nm in names(by_sub)) {
  writeLines(sprintf("[%s]", nm), con_ref)
  writeLines(by_sub[[nm]], con_ref)
  writeLines("", con_ref)
}

## =====================================================================
## (B) DSMZ assignment from TCGA references (independent label transfer)
## =====================================================================
cat("\n## =====================================================================\n")
cat("## (B) Assign DSMZ subtypes using TCGA reference centroids\n")
cat("## =====================================================================\n")
dsmz_dir <- file.path(outdir, "dsmz_from_tcga_reference")
dir.create(dsmz_dir, showWarnings = FALSE, recursive = TRUE)

# Load references (genes x selected TCGA samples) and their grouping
cat(sprintf("[B.1] Loading TCGA reference matrix (%d genes x %d samples).\n", nrow(M_pam50_tcga_ref), ncol(M_pam50_tcga_ref)))
M_ref <- M_pam50_tcga_ref # Use object from (A)
sel_info <- sel_df         # Use object from (A)

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
cat(sprintf("[B.2] Building %d TCGA reference centroids from %d PAM50 genes.\n", ncol(centroids), nrow(centroids)))

# Prepare DSMZ PAM50 (genes x DSMZ samples)
M_pam50_dsmz <- M_pam50[, colnames(V_dsmz), drop = FALSE]
cat(sprintf("[B.3] Subsetting PAM50 matrix for %d DSMZ samples (%d genes x %d samples).\n",
            ncol(M_pam50_dsmz), nrow(M_pam50_dsmz), ncol(M_pam50_dsmz)))

# Use only common PAM50 genes
g <- intersect(rownames(centroids), rownames(M_pam50_dsmz))
cat(sprintf("[B.4] Found %d overlapping PAM50 genes for correlation.\n", length(g)))
if (length(g) < 10) stop("[FATAL] Too few overlapping PAM50 genes between reference and DSMZ.")

# Spearman correlation for robustness (rank-based, scale-free)
assign_one <- function(x, C) {
  cs <- apply(C, 2, function(cen) suppressWarnings(cor(x, cen, method = "spearman", use = "pairwise.complete.obs")))
  lab <- names(which.max(cs))
  score <- max(cs, na.rm = TRUE)
  c(label = lab, score = score)
}

# Compute assignment per DSMZ sample
cat(sprintf("[B.5] Assigning subtypes to %d DSMZ samples via Spearman correlation to centroids...\n", ncol(M_pam50_dsmz)))
dsmz_assign <- t(apply(M_pam50_dsmz[g, , drop = FALSE], 2, assign_one, C = centroids[g, , drop = FALSE]))
dsmz_assign <- as.data.frame(dsmz_assign, stringsAsFactors = FALSE)
dsmz_assign$sample <- rownames(dsmz_assign)
dsmz_assign$label  <- factor(dsmz_assign$label, levels = subtype_labels)

cat("[B.5] DSMZ predicted subtype counts:\n")
out_table_dsmz <- capture.output(print(table(dsmz_assign$label)))
cat(paste("    ", out_table_dsmz, collapse = "\n"), "\n")

# Save table
out_csv_dsmz <- file.path(dsmz_dir, "dsmz_pred_subtypes_from_tcga_top20.csv")
cat(sprintf("[B.6] Saving DSMZ assignment table to: %s\n", out_csv_dsmz))
write.csv(dsmz_assign, out_csv_dsmz, row.names = FALSE)

# Quick DSMZ-only UMAP (PAM50 space), color by predicted subtype
cat("[B.7] Generating UMAP of DSMZ samples (PAM50 space, 'cosine' metric) colored by prediction.\n")
set.seed(42)
X_dsmz <- t(scale(M_pam50_dsmz[g, , drop = FALSE]))  # samples x genes (scaled)
emb_d <- uwot::umap(X_dsmz, n_neighbors = 10, min_dist = 0.3, metric = "cosine")
emb_d <- as.data.frame(emb_d); colnames(emb_d) <- c("UMAP1","UMAP2")
emb_d$sample <- rownames(X_dsmz)
emb_d$pred   <- factor(dsmz_assign$label[match(emb_d$sample, dsmz_assign$sample)], levels = subtype_labels)

p_pred <- ggplot(emb_d, aes(UMAP1, UMAP2, color = pred)) +
  geom_point(alpha = 0.95, size = 2.1) + theme_bw() +
  ggtitle("DSMZ (PAM50) - predicted PAM50 subtype (nearest TCGA centroid)")
out_pdf_umap_dsmz <- file.path(dsmz_dir, "umap_dsmz_predicted_labels.pdf")
ggsave(out_pdf_umap_dsmz, p_pred, width = 6.5, height = 5.5)
cat(sprintf("[B.7] Saved DSMZ UMAP plot to: %s\n", out_pdf_umap_dsmz))

# Also export a simple 2-column mapping for downstream use
out_tsv_dsmz <- file.path(dsmz_dir, "dsmz_predicted_label_map.tsv")
cat(sprintf("[B.7] Saving 2-column DSMZ mapping to: %s\n", out_tsv_dsmz))
write.table(
  dsmz_assign[, c("sample","label")],
  file = out_tsv_dsmz,
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("[INFO] DSMZ assignment complete. Results in: ", dsmz_dir, "\n")

# ---------- (6) UMAP EMBEDDING ----------
cat("\n---------- (6) JOINT UMAP EMBEDDING (TCGA + DSMZ) ----------\n")
# Create UMAP visualization of joint TCGA+DSMZ data
set.seed(42)  # Set random seed for reproducibility
genes_umap <- top_var_genes(V_adj, n = 3000)  # Select most variable genes for UMAP
cat(sprintf("[INFO] Selecting top %d HVGs from joint data for UMAP.\n", length(genes_umap)))
Xu <- t(V_adj[genes_umap, , drop = FALSE])    # Transpose to samples x genes (required for UMAP)
cat(sprintf("[INFO] Preparing joint matrix: %d samples x %d HVGs.\n", nrow(Xu), ncol(Xu)))
emb <- umap(Xu, n_neighbors = 20, min_dist = 0.3, metric = "cosine")  # Run UMAP
cat("[INFO] UMAP complete. Adding annotations...\n")
emb <- as.data.frame(emb)  # Convert to data frame
colnames(emb) <- c("UMAP1","UMAP2")  # Set column names
emb$sample <- rownames(Xu)  # Add sample names
emb$dataset <- ifelse(emb$sample %in% colnames(V_tcga), "TCGA", "DSMZ")  # Add dataset labels
emb$discovery_cluster <- clusters_kmeans[emb$sample]  # Add discovery clusters

# Add TCGA subtypes if available
emb$TCGA_subtype <- NA  # Initialize with NA
if (!is.na(tcga_sub_col)) {  # If subtype column exists
  cat("[INFO] Adding TCGA subtype annotations to UMAP data.\n")
  emb$TCGA_subtype[match(colnames(V_tcga), emb$sample)] <-  # For TCGA samples
    as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])  # Add subtype labels
}

# Add DSMZ predicted subtypes
cat("[INFO] Adding DSMZ predicted subtype annotations to UMAP data.\n")
emb$DSMZ_predicted_subtype <- NA
dsmz_pred_map <- setNames(dsmz_assign$label, dsmz_assign$sample)
emb$DSMZ_predicted_subtype[match(colnames(V_dsmz), emb$sample)] <- 
  as.character(dsmz_pred_map[colnames(V_dsmz)])

# Create UMAP plots
cat("[INFO] Generating UMAP plots (by dataset, TCGA subtype, DSMZ prediction).\n")
p1 <- ggplot(emb, aes(UMAP1, UMAP2, color = dataset)) +  # Plot by dataset
  geom_point(alpha = 0.8, size = 1.6) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (dataset)")  # Title

p2 <- ggplot(subset(emb, dataset=="TCGA"), aes(UMAP1, UMAP2, color = TCGA_subtype)) +  # TCGA only
  geom_point(alpha = 0.8, size = 1.6) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (TCGA BRCA subtype)")  # Title

p3 <- ggplot(subset(emb, dataset=="DSMZ"), aes(UMAP1, UMAP2, color = DSMZ_predicted_subtype)) +  # DSMZ only
  geom_point(alpha = 0.9, size = 1.9) +  # Add points
  theme_bw() +  # Clean theme
  ggtitle("UMAP (DSMZ predicted subtype)")  # Title

# Save individual plots
ggsave(file.path(outdir, "umap_dataset.pdf"), p1, width = 7, height = 6)
ggsave(file.path(outdir, "umap_tcga_subtype.pdf"), p2, width = 7, height = 6)
ggsave(file.path(outdir, "umap_dsmz_predicted_subtype.pdf"), p3, width = 7, height = 6)
cat(sprintf("[INFO] Saved individual UMAP plots to: %s\n", outdir))

# Save combined plot (try patchwork, fallback to separate PDF)
if (requireNamespace("patchwork", quietly = TRUE)) {  # If patchwork is available
  cat("[INFO] 'patchwork' package found. Creating combined UMAP PDF.\n")
  library(patchwork)  # Load patchwork
  combined <- (p1 | p2 | p3) + plot_layout(guides = "collect")  # Combine plots
  ggsave(file.path(outdir, "umap_combined.pdf"), combined, width = 15, height = 5)  # Save combined
} else {  # Fallback method
  cat("[INFO] 'patchwork' not found. Saving UMAPs to separate pages in combined PDF.\n")
  pdf(file.path(outdir, "umap_combined.pdf"), width = 15, height = 5)  # Open PDF
  print(p1)  # Print plot 1
  print(p2)  # Print plot 2
  print(p3)  # Print plot 3
  dev.off()  # Close PDF
}

# Save embeddings for further analysis
out_csv_emb <- file.path(outdir, "umap_embeddings.csv")
cat(sprintf("[INFO] Saving full UMAP embedding coordinates to: %s\n", out_csv_emb))
write.csv(emb, out_csv_emb, row.names = FALSE)  # Save coordinates

## ---- DSMZ marker heatmap (ESR1/PGR/ERBB2/KRT5/KRT14/MKI67) ----
cat("\n---------- (7) DSMZ MARKER HEATMAP ----------\n")

# markers of interest
markers <- c("ESR1","PGR","ERBB2","KRT5","KRT14","MKI67")
cat(sprintf("[INFO] Generating heatmap for %d markers: %s\n", length(markers), paste(markers, collapse=", ")))

# 1) Map symbols -> Ensembl using your union map (ens2sym_union) and subset DSMZ matrix (V_dsmz)
stopifnot(exists("ens2sym_union"), exists("V_dsmz"))
cat("[INFO] Mapping marker symbols to Ensembl IDs...\n")
sym2ens <- setNames(ens2sym_union$ensembl, toupper(ens2sym_union$symbol))
sel_ens <- sym2ens[toupper(markers)]
sel_ens <- sel_ens[!is.na(sel_ens)]

if (length(sel_ens) == 0L) {
  cat("[WARN] None of the requested markers found in mapping.\n")
} else {
  M <- V_dsmz[intersect(sel_ens, rownames(V_dsmz)), , drop = FALSE]
  cat(sprintf("[INFO] Found %d markers present in DSMZ expression data.\n", nrow(M)))
  
  if (nrow(M) == 0L) {
    cat("[WARN] No marker rows present in V_dsmz after mapping.\n")
  } else {
    # Set display rownames back to symbols (keep order as in `markers`)
    ens2sym_lookup <- setNames(names(sel_ens), sel_ens)       # ENSG -> SYMBOL
    rownames(M) <- ens2sym_lookup[rownames(M)]
    M <- M[intersect(markers, rownames(M)), , drop = FALSE]    # order rows like markers
    
    # 2) Row z-score for readability (gene-wise)
    cat("[INFO] Applying gene-wise z-score to marker matrix.\n")
    Mz <- zscore_by_gene(M)
    
    # 3) Plot with well-defined legends
    out_pdf <- file.path(outdir, "heatmap_dsmz_markers.pdf")
    dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
    cat(sprintf("[INFO] Saving DSMZ marker heatmap to: %s\n", out_pdf))
    pdf(out_pdf, width = 9, height = 5.5)
    pheatmap(Mz,
             scale = "none",
             clustering_method = "average",
             show_colnames = FALSE,
             show_rownames = TRUE,
             main = "DSMZ marker panel (z-scored)",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-3, 3, length.out = 101),
             legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
             legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
             cluster_rows = FALSE,
             cluster_cols = TRUE)
    dev.off()
  }
}

cat("\n[OK] Pipeline complete. Results in:", outdir, "\n")  # Print completion message
