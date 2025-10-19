# Replace the path with the one you copied in the previous step
.libPaths(c("/home/chu25/miniconda3/envs/r_sva_fix/lib/R/library", .libPaths()))

# nolint start
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
})

# ---------- CONFIG ----------
tcga_se_rds   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"
dsmz_rds      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
dsmz_meta_csv <- "/home/chu25/data/dsmz/DSMZ_metadata.csv"
outdir <- "/home/chu25/dsmz/results/brca_vst_dual"
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
config <- list(batch_adjust_method = "combat")

# ---------- HELPERS ----------
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)
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
  if (is.null(purity)) return(tcga_log)
  cat("[INFO] Applying purity adjustment on VST data...\n")
  keep_samp <- intersect(colnames(tcga_log), names(purity)[!is.na(purity)])
  if (length(keep_samp) < 3) {
    cat("[WARN] Too few samples with valid purity to adjust; returning input.\n")
    return(tcga_log)
  }
  M <- tcga_log[, keep_samp, drop = FALSE]
  pu <- purity[keep_samp]
  ct <- apply(M, 1, function(v) {
    if (var(v) == 0 || var(pu) == 0) return(list(p.value = 1, estimate = 0))
    suppressWarnings(cor.test(as.numeric(v), pu, method = "spearman"))
  })
  pv <- vapply(ct, `[[`, numeric(1), "p.value")
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))
  pv[is.na(pv)] <- 1
  padj <- p.adjust(pv, "BH")
  drop <- names(which(padj < 0.01 & rh < -0.4))
  if (length(drop)) {
    cat(sprintf("[INFO] Removing %d purity-associated genes\n", length(drop)))
    M <- M[setdiff(rownames(M), drop), , drop = FALSE]
  }
  infilt <- 1 - pu
  design <- model.matrix(~ infilt)
  fit <- lmFit(M, design)
  beta <- fit$coefficients[, "infilt", drop = FALSE]
  Madj <- as.matrix(M) - beta %*% t(infilt)
  out <- tcga_log
  common_g <- intersect(rownames(Madj), rownames(out))
  out[common_g, colnames(Madj)] <- Madj[common_g, ]
  out
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

batch <- factor(c(rep("TCGA", ncol(tcga_counts)), rep("DSMZ", ncol(dsmz_counts))))
all_counts <- Xc_raw

if (nrow(Xc_raw) > 50000 || ncol(Xc_raw) > 5000) {
  cat("[WARN] Large merged matrix (%d genes x %d samples); consider subsampling.\n", 
      nrow(Xc_raw), ncol(Xc_raw))
}

coldata <- tibble(
  sample = colnames(all_counts),
  dataset = ifelse(colnames(all_counts) %in% colnames(tcga_counts), "TCGA", "DSMZ")
) %>% column_to_rownames("sample")

if (any(duplicated(rownames(all_counts)))) {
  all_counts <- rowsum(all_counts, rownames(all_counts), reorder = TRUE)
}
dds <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
V <- assay(vsd)
V_tcga <- V[, colnames(tcga_counts), drop = FALSE]
V_dsmz <- V[, colnames(dsmz_counts), drop = FALSE]

# ---- VST individually per dataset (in addition to joint VST) ----
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

purity <- load_purity_data(tcga_se, NULL)
if (!is.null(purity)) {
  cat("[INFO] Applying purity adjustment to TCGA VST data...\n")
  V_tcga <- purity_adjust_log(V_tcga, purity)
}

cat("[INFO] Batch adjustment (global)...\n")
method <- tolower(config$batch_adjust_method)
V_adj <- V
if (method == "none") {
  cat("[INFO] Skipping batch correction.\n")
} else if (method == "combat_seq") {
  cat("[WARN] ComBat_seq requires raw counts; skipping batch correction for VST data.\n")
} else if (method == "combat") {
  V_adj[is.na(V_adj)] <- 0
  set.seed(42)
  V_adj <- sva::ComBat(dat = as.matrix(V_adj), batch = batch, mod = NULL, 
                       par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE)
} else {
  stop("[ERROR] Unknown batch_adjust_method. Use 'none', 'combat_seq', or 'combat'.")
}
V_tcga <- V_adj[, colnames(tcga_counts), drop = FALSE]
V_dsmz <- V_adj[, colnames(dsmz_counts), drop = FALSE]

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

pca_pre  <- make_pca_plot(V,      batch, "PCA (pre-batch correction)",  file.path(outdir, "pca_pre_batch.pdf"))
pca_post <- make_pca_plot(V_adj,  batch, "PCA (post-batch correction)", file.path(outdir, "pca_post_batch.pdf"))

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

# ---------- (1) UNSUPERVISED DISCOVERY CLUSTERING (K-MEANS) ----------
set.seed(42)
hvgs <- top_var_genes(V_adj, n = 3000)
X_discovery <- t(V_adj[hvgs, , drop = FALSE])
if (nrow(X_discovery) > 5000) {
  cat("[WARN] Large sample set (%d) for k-means; consider subsampling.\n", nrow(X_discovery))
}
X_scaled <- scale(X_discovery)
writeLines(hvgs, file.path(outdir, "HVGs_used_for_discovery_3k.txt"))

k_grid <- 2:8
res_list <- list()
metrics <- lapply(k_grid, function(k) {
  km <- kmeans(X_scaled, centers = k, nstart = 50, iter.max = 100)
  wcss <- sum(km$withinss)
  sil <- silhouette(km$cluster, dist(X_scaled))
  sil_mean <- mean(sil[, "sil_width"])
  ch <- tryCatch(calinhara(X_scaled, km$cluster), error = function(e) NA_real_)
  db <- tryCatch({
    m <- intCriteria(as.matrix(X_scaled), as.integer(km$cluster), "davies_bouldin")
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

cluster_df <- tibble(sample = rownames(X_scaled),
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

# UMAP for k-means
set.seed(42)
emb <- uwot::umap(X_discovery, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2")
emb$sample <- rownames(X_discovery)
emb$dataset <- dataset_lab
emb$cluster <- clusters_kmeans

p_umap <- ggplot(emb, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (k-means, HVG space)")

ggsave(file.path(outdir, "umap_kmeans_HVG.pdf"), p_umap, width = 7, height = 6)
write.csv(emb, file.path(outdir, "umap_kmeans_HVG_coords.csv"), row.names = FALSE)

# ---------- (2) HIERARCHICAL CLUSTERING (HVG, k=5, EUCLIDEAN, AVERAGE LINKAGE) ----------
cat("[INFO] Performing hierarchical clustering (HVG)...\n")
if (nrow(X_scaled) < 5) {
  stop(sprintf("[FATAL] Too few samples (%d) for hierarchical clustering (need at least 5).", nrow(X_scaled)))
}

# Distance metric choice justification:
# - Euclidean distance: Most common for gene expression clustering, forms compact spherical clusters
# - Manhattan distance: More sensitive to outliers, forms hyper-rectangular clusters
# - Correlation distance: Captures expression pattern similarity regardless of magnitude
# For breast cancer subtype discovery, Euclidean distance is preferred as it groups samples
# with similar overall expression profiles, which aligns with molecular subtype definitions.
cat("[INFO] Using Euclidean distance for hierarchical clustering (standard for gene expression)\n")

# Alternative distance metrics (uncomment to use):
# d <- dist(X_scaled, method = "manhattan")  # Manhattan (L1 norm) - more robust to outliers
# d <- as.dist(1 - cor(t(X_scaled)))         # Correlation distance - pattern similarity
# d <- dist(X_scaled, method = "maximum")    # Chebyshev distance - maximum coordinate difference

d <- dist(X_scaled, method = "euclidean")
hc <- hclust(d, method = "average")
clusters_hc <- cutree(hc, k = 5)
subtype_labels <- c("HER2-high", "HER2-low", "LumA", "LumB", "Basal")

# Label clusters by correlating with TCGA PAM50 centroids
tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype")
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1]
if (!is.na(tcga_sub_col)) {
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga)))
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE))
  cl_means <- sapply(1:5, function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE))
  g_common <- intersect(rownames(sub_means), hvgs)
  corr <- cor(cl_means[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  cluster_to_subtype <- apply(corr, 1, function(v) sub_lvls[which.max(v)])
  clusters_hc <- factor(clusters_hc, labels = subtype_labels[match(cluster_to_subtype, sub_lvls)])
} else {
  clusters_hc <- factor(clusters_hc, labels = subtype_labels)
}

cluster_hc_df <- tibble(sample = rownames(X_scaled),
                        dataset = dataset_lab,
                        cluster = as.character(clusters_hc))
write.csv(cluster_hc_df, file.path(outdir, "clusters_hierarchical_HVG.csv"), row.names = FALSE)

safe_pdf(file.path(outdir, "dendrogram_hierarchical_HVG.pdf"), {
  plot(hc, main = "Hierarchical Clustering (HVG, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5)
  abline(h = hc$height[length(hc$height) - 4], col = "red", lty = 2)
})

# UMAP for hierarchical clustering (HVG)
set.seed(42)
emb_hc <- uwot::umap(X_scaled, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
emb_hc <- as.data.frame(emb_hc); colnames(emb_hc) <- c("UMAP1","UMAP2")
emb_hc$sample <- rownames(X_scaled)
emb_hc$dataset <- dataset_lab
emb_hc$cluster <- clusters_hc

p_hc_umap <- ggplot(emb_hc, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (Hierarchical Clustering, HVG space)")

ggsave(file.path(outdir, "umap_hierarchical_HVG.pdf"), p_hc_umap, width = 7, height = 6)
write.csv(emb_hc, file.path(outdir, "umap_hierarchical_HVG_coords.csv"), row.names = FALSE)

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
cat("[INFO] Using Euclidean distance for PAM50 hierarchical clustering (consistent with HVG analysis)\n")
d_pam50 <- dist(X_pam50, method = "euclidean")
hc_pam50 <- hclust(d_pam50, method = "average")
clusters_hc_pam50 <- cutree(hc_pam50, k = 5)

if (!is.na(tcga_sub_col)) {
  cl_means_pam50 <- sapply(1:5, function(cn) rowMeans(M_pam50[, clusters_hc_pam50 == cn, drop = FALSE], na.rm = TRUE))
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50))
  corr_pam50 <- cor(cl_means_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman")
  cluster_to_subtype_pam50 <- apply(corr_pam50, 1, function(v) sub_lvls[which.max(v)])
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = subtype_labels[match(cluster_to_subtype_pam50, sub_lvls)])
} else {
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = subtype_labels)
}

cluster_hc_pam50_df <- tibble(sample = rownames(X_pam50),
                              dataset = dataset_lab,
                              cluster = as.character(clusters_hc_pam50))
write.csv(cluster_hc_pam50_df, file.path(outdir, "clusters_hierarchical_PAM50.csv"), row.names = FALSE)

safe_pdf(file.path(outdir, "dendrogram_hierarchical_PAM50.pdf"), {
  plot(hc_pam50, main = "Hierarchical Clustering (PAM50, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5)
  abline(h = hc_pam50$height[length(hc_pam50$height) - 4], col = "red", lty = 2)
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

# ---------- (5) UMAP + HDBSCAN ITERATIVE SUBTYPE ASSIGNMENT (HVG) ----------
cat("[INFO] Performing iterative UMAP + HDBSCAN (HVG)...\n")
set.seed(42)
if (nrow(X_scaled) < 20) {
  stop(sprintf("[FATAL] Too few samples (%d) for UMAP+HDBSCAN (need at least 20).", nrow(X_scaled)))
}
if (nrow(X_scaled) > 5000) {
  cat("[WARN] Large sample set (%d) for UMAP+HDBSCAN; consider subsampling.\n", nrow(X_scaled))
}

emb1 <- uwot::umap(X_scaled, n_neighbors = 20, min_dist = 0.3, metric = "cosine")
hdb1 <- hdbscan(emb1, minPts = 10)
clusters_hdb <- rep("Unassigned", nrow(X_scaled))
basal_idx <- which(hdb1$cluster == which.max(table(hdb1$cluster[hdb1$cluster != 0])))
clusters_hdb[basal_idx] <- "Basal"

non_basal <- which(clusters_hdb == "Unassigned")
if (length(non_basal) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification.\n", length(non_basal))
} else {
  emb2 <- uwot::umap(X_scaled[non_basal, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  hdb2 <- hdbscan(emb2, minPts = 5)
  her2_high_idx <- non_basal[which(hdb2$cluster == which.max(table(hdb2$cluster[hdb2$cluster != 0])))]
  clusters_hdb[her2_high_idx] <- "HER2-high"
}

non_her2_high <- which(clusters_hdb == "Unassigned")
if (length(non_her2_high) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification.\n", length(non_her2_high))
} else {
  emb3 <- uwot::umap(X_scaled[non_her2_high, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  hdb3 <- hdbscan(emb3, minPts = 5)
  her2_low_idx <- non_her2_high[which(hdb3$cluster == which.max(table(hdb3$cluster[hdb3$cluster != 0])))]
  clusters_hdb[her2_low_idx] <- "HER2-low"
}

luminal_idx <- which(clusters_hdb == "Unassigned")
if (length(luminal_idx) < 2) {
  cat("[WARN] Too few Luminal samples (%d); assigning all as LumA.\n", length(luminal_idx))
  clusters_hdb[luminal_idx] <- "LumA"
} else {
  emb_lum <- uwot::umap(X_scaled[luminal_idx, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "cosine")
  km_lum <- kmeans(emb_lum, centers = 2, nstart = 50)
  lum_clusters <- ifelse(km_lum$cluster == 1, "LumA", "LumB")
  clusters_hdb[luminal_idx] <- lum_clusters
}

clusters_hdb <- factor(clusters_hdb, levels = subtype_labels)
hdb_df <- tibble(sample = rownames(X_scaled),
                 dataset = dataset_lab,
                 subtype = clusters_hdb)
write.csv(hdb_df, file.path(outdir, "clusters_hdbscan_HVG.csv"), row.names = FALSE)

emb_hdb <- as.data.frame(emb1); colnames(emb_hdb) <- c("UMAP1","UMAP2")
emb_hdb$sample <- rownames(X_scaled)
emb_hdb$dataset <- dataset_lab
emb_hdb$subtype <- clusters_hdb

p_hdb_umap <- ggplot(emb_hdb, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (HDBSCAN Subtypes, HVG space)")

ggsave(file.path(outdir, "umap_hdbscan_HVG.pdf"), p_hdb_umap, width = 7, height = 6)
write.csv(emb_hdb, file.path(outdir, "umap_hdbscan_HVG_coords.csv"), row.names = FALSE)

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
hdb1_pam50 <- hdbscan(emb1_pam50, minPts = 10)
clusters_hdb_pam50 <- rep("Unassigned", nrow(X_pam50))
basal_idx_pam50 <- which(hdb1_pam50$cluster == which.max(table(hdb1_pam50$cluster[hdb1_pam50$cluster != 0])))
clusters_hdb_pam50[basal_idx_pam50] <- "Basal"

non_basal_pam50 <- which(clusters_hdb_pam50 == "Unassigned")
if (length(non_basal_pam50) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification (PAM50).\n", length(non_basal_pam50))
} else {
  emb2_pam50 <- uwot::umap(X_pam50[non_basal_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  hdb2_pam50 <- hdbscan(emb2_pam50, minPts = 5)
  her2_high_idx_pam50 <- non_basal_pam50[which(hdb2_pam50$cluster == which.max(table(hdb2_pam50$cluster[hdb2_pam50$cluster != 0])))]
  clusters_hdb_pam50[her2_high_idx_pam50] <- "HER2-high"
}

non_her2_high_pam50 <- which(clusters_hdb_pam50 == "Unassigned")
if (length(non_her2_high_pam50) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification (PAM50).\n", length(non_her2_high_pam50))
} else {
  emb3_pam50 <- uwot::umap(X_pam50[non_her2_high_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  hdb3_pam50 <- hdbscan(emb3_pam50, minPts = 5)
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
hdb_pam50_df <- tibble(sample = rownames(X_pam50),
                       dataset = dataset_lab,
                       subtype = clusters_hdb_pam50)
write.csv(hdb_pam50_df, file.path(outdir, "clusters_hdbscan_PAM50.csv"), row.names = FALSE)

emb_hdb_pam50 <- as.data.frame(emb1_pam50); colnames(emb_hdb_pam50) <- c("UMAP1","UMAP2")
emb_hdb_pam50$sample <- rownames(X_pam50)
emb_hdb_pam50$dataset <- dataset_lab
emb_hdb_pam50$subtype <- clusters_hdb_pam50

p_hdb_pam50_umap <- ggplot(emb_hdb_pam50, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (HDBSCAN Subtypes, PAM50 space)")

ggsave(file.path(outdir, "umap_hdbscan_PAM50.pdf"), p_hdb_pam50_umap, width = 7, height = 6)
write.csv(emb_hdb_pam50, file.path(outdir, "umap_hdbscan_PAM50_coords.csv"), row.names = FALSE)

# ---------- (7) STATISTICAL ANALYSIS ----------
cat("[INFO] Performing statistical analysis...\n")
dir.create(file.path(outdir, "stats_results"), showWarnings = FALSE)

# Fisher’s Exact Test: Cluster vs. Dataset
fisher_results <- list()
cont_table_kmeans <- table(clusters_kmeans, dataset_lab)
cont_table_hc <- table(clusters_hc, dataset_lab)
cont_table_hdb <- table(clusters_hdb, dataset_lab)
cont_table_hdb_pam50 <- table(clusters_hdb_pam50, dataset_lab)

if (all(dim(cont_table_kmeans) >= 2) && sum(cont_table_kmeans) >= 5) {
  fisher_kmeans <- fisher.test(cont_table_kmeans, simulate.p.value = nrow(cont_table_kmeans) > 2)
  fisher_results[["kmeans_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_kmeans_vs_dataset",
    p_value = fisher_kmeans$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] Insufficient data for Fisher’s test (k-means vs. dataset).\n")
}

if (all(dim(cont_table_hc) >= 2) && sum(cont_table_hc) >= 5) {
  fisher_hc <- fisher.test(cont_table_hc, simulate.p.value = nrow(cont_table_hc) > 2)
  fisher_results[["hc_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hc_vs_dataset",
    p_value = fisher_hc$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] Insufficient data for Fisher’s test (hierarchical vs. dataset).\n")
}

if (all(dim(cont_table_hdb) >= 2) && sum(cont_table_hdb) >= 5) {
  fisher_hdb <- fisher.test(cont_table_hdb, simulate.p.value = nrow(cont_table_hdb) > 2)
  fisher_results[["hdb_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_vs_dataset",
    p_value = fisher_hdb$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN vs. dataset).\n")
}

if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
  fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2)
  fisher_results[["hdb_pam50_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_pam50_vs_dataset",
    p_value = fisher_hdb_pam50$p.value,
    statistic = NA_real_
  )
} else {
  cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN PAM50 vs. dataset).\n")
}

# Fisher’s Exact Test: Cluster vs. TCGA PAM50 Subtype
if (!is.na(tcga_sub_col)) {
  tcga_samples <- colnames(V_tcga)
  lab_tcga <- as.character(colData(tcga_se)[tcga_samples, tcga_sub_col])
  cont_table_kmeans_pam50 <- table(clusters_kmeans[tcga_samples], lab_tcga)
  cont_table_hc_pam50 <- table(clusters_hc[tcga_samples], lab_tcga)
  cont_table_hdb_pam50 <- table(clusters_hdb[tcga_samples], lab_tcga)
  cont_table_hdb_pam50_sub <- table(clusters_hdb_pam50[tcga_samples], lab_tcga)
  
  if (all(dim(cont_table_kmeans_pam50) >= 2) && sum(cont_table_kmeans_pam50) >= 5) {
    fisher_kmeans_pam50 <- fisher.test(cont_table_kmeans_pam50, simulate.p.value = nrow(cont_table_kmeans_pam50) > 2)
    fisher_results[["kmeans_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_kmeans_vs_PAM50",
      p_value = fisher_kmeans_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (k-means vs. PAM50).\n")
  }
  
  if (all(dim(cont_table_hc_pam50) >= 2) && sum(cont_table_hc_pam50) >= 5) {
    fisher_hc_pam50 <- fisher.test(cont_table_hc_pam50, simulate.p.value = nrow(cont_table_hc_pam50) > 2)
    fisher_results[["hc_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hc_vs_PAM50",
      p_value = fisher_hc_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (hierarchical vs. PAM50).\n")
  }
  
  if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
    fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2)
    fisher_results[["hdb_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_vs_PAM50",
      p_value = fisher_hdb_pam50$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN vs. PAM50).\n")
  }
  
  if (all(dim(cont_table_hdb_pam50_sub) >= 2) && sum(cont_table_hdb_pam50_sub) >= 5) {
    fisher_hdb_pam50_sub <- fisher.test(cont_table_hdb_pam50_sub, simulate.p.value = nrow(cont_table_hdb_pam50_sub) > 2)
    fisher_results[["hdb_pam50_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_pam50_vs_PAM50",
      p_value = fisher_hdb_pam50_sub$p.value,
      statistic = NA_real_
    )
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN PAM50 vs. PAM50).\n")
  }
}

fisher_df <- do.call(rbind, fisher_results)
write.csv(fisher_df, file.path(outdir, "stats_results", "fisher_tests.csv"), row.names = FALSE)

# Kruskal-Wallis and DSCF Tests: Marker gene expression across subtypes
kw_results <- list()
dscf_results <- list()
for (method in c("kmeans", "hc", "hdb", "hdb_pam50")) {
  clusters <- switch(method,
                     kmeans = clusters_kmeans,
                     hc = clusters_hc,
                     hdb = clusters_hdb,
                     hdb_pam50 = clusters_hdb_pam50)
  mat <- if (method == "hdb_pam50") M_pam50 else V_adj
  genes <- if (method == "hdb_pam50") rownames(M_pam50) else marker_ens
  for (gene in genes) {
    gene_sym <- if (method == "hdb_pam50") gene else ens2sym_union$symbol[match(gene, ens2sym_union$ensembl)]
    gene_sym <- ifelse(is.na(gene_sym) || gene_sym == "", gene, gene_sym)
    expr <- mat[gene, ]
    if (sum(!is.na(expr)) < length(expr) * 0.5) {
      cat(sprintf("[WARN] Skipping Kruskal-Wallis for %s (%s, %s): too many NA values.\n", gene_sym, gene, method))
      next
    }
    kw_test <- kruskal.test(expr ~ clusters)
    if (is.na(kw_test$p.value)) {
      cat(sprintf("[WARN] Kruskal-Wallis failed for %s (%s, %s): insufficient variation.\n", gene_sym, gene, method))
      next
    }
    kw_results[[paste(method, gene_sym, sep = "_")]] <- data.frame(
      test = sprintf("Kruskal_Wallis_%s_%s", method, gene_sym),
      gene = gene,
      symbol = gene_sym,
      p_value = kw_test$p.value,
      statistic = kw_test$statistic
    )
    if (kw_test$p.value < 0.05) {
      # work on samples with non-NA expression only
      keep <- !is.na(expr)
      expr_ok <- expr[keep]
      grp_ok  <- droplevels(clusters[keep])

      # need at least two non-empty groups
      if (nlevels(grp_ok) >= 2 && length(unique(grp_ok)) >= 2) {
        # run DSCF on present groups only
        dscf_test <- try(PMCMRplus::dscfAllPairsTest(expr_ok, grp_ok, p.adjust.method = "BH"), silent = TRUE)

        if (!inherits(dscf_test, "try-error") && !is.null(dscf_test$p.value)) {
          p_mat <- as.matrix(dscf_test$p.value)

          # ---- robust guards for degenerate outputs ----
          if (is.null(dim(p_mat))) next  # not a matrix
          if (min(dim(p_mat)) < 2) next  # need at least 2 groups

          gr_levels <- colnames(p_mat)
          if (is.null(gr_levels) || length(gr_levels) < 2) {
            gr_levels <- rownames(p_mat)
          }
          if (is.null(gr_levels) || length(gr_levels) < 2) {
            # fallback to actually present groups in the data
            gr_levels <- levels(grp_ok)
            if (length(gr_levels) < 2) next
            if (is.null(rownames(p_mat))) rownames(p_mat) <- gr_levels
            if (is.null(colnames(p_mat))) colnames(p_mat) <- gr_levels
          }

          # make sure names are unique, drop any NA labels
          gr_levels <- unique(gr_levels[!is.na(gr_levels)])
          if (length(gr_levels) < 2) next

          comb <- utils::combn(gr_levels, 2)

          get_p <- function(a, b) {
            if (!(a %in% rownames(p_mat) && b %in% colnames(p_mat)) &&
                !(b %in% rownames(p_mat) && a %in% colnames(p_mat))) return(NA_real_)
            if (!is.na(p_mat[a, b])) return(p_mat[a, b])
            if (!is.na(p_mat[b, a])) return(p_mat[b, a])
            NA_real_
          }
          pvals <- vapply(seq_len(ncol(comb)), function(i) get_p(comb[1, i], comb[2, i]), numeric(1))

          dscf_df <- data.frame(
            test     = sprintf("DSCF_%s_%s_%s_vs_%s", method, gene_sym, comb[1, ], comb[2, ]),
            gene     = gene,
            symbol   = gene_sym,
            subtype1 = comb[1, ],
            subtype2 = comb[2, ],
            p_value  = pvals,
            stringsAsFactors = FALSE
          )
          dscf_df <- dscf_df[is.finite(dscf_df$p_value), , drop = FALSE]
          if (nrow(dscf_df) > 0)
            dscf_results[[paste(method, gene_sym, sep = "_")]] <- dscf_df
        }
      }
    }
    if (kw_test$p.value < 0.05) {
      p_box <- ggplot(data.frame(expr = mat[gene, ], subtype = clusters), 
                      aes(x = subtype, y = expr, fill = subtype)) +
        geom_boxplot() + theme_bw() +
        labs(title = sprintf("Expression of %s (%s)", gene_sym, method), y = "VST Expression") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(outdir, "stats_results", sprintf("boxplot_%s_%s.pdf", method, gene_sym)), p_box, width = 7, height = 5)
    }
  }
}

kw_df <- do.call(rbind, kw_results)
write.csv(kw_df, file.path(outdir, "stats_results", "kruskal_wallis_tests.csv"), row.names = FALSE)
dscf_df <- do.call(rbind, dscf_results)
if (nrow(dscf_df) > 0) {
  write.csv(dscf_df, file.path(outdir, "stats_results", "dscf_tests.csv"), row.names = FALSE)
}

# ---------- (8) MARKER HEATMAPS ----------
ens2sym <- make_ens2sym(dsmz_raw)
ens_vec <- strip_ensver(as_df(dsmz_raw)$Ensembl_ID)
dup_ids <- sum(duplicated(ens_vec))
dup_unique <- sum(table(ens_vec) > 1, na.rm = TRUE)
cat(sprintf("[INFO] DSMZ Ensembl duplicate rows: %d | duplicate IDs: %d (collapsed)\n",
            dup_ids, dup_unique))

marker_ens <- rownames(ens2sym)[tolower(ens2sym$symbol) %in% tolower(marker_symbols)]
marker_ens <- intersect(unique(marker_ens), rownames(V_adj))

if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  ann <- data.frame(dataset = dataset_lab, 
                    kmeans_cluster = clusters_kmeans,
                    hc_subtype = clusters_hc,
                    hdb_subtype = clusters_hdb,
                    hdb_pam50_subtype = clusters_hdb_pam50,
                    row.names = colnames(V_mark))
  safe_pdf(file.path(outdir, "heatmap_pam50_all_clusters.pdf"), {
    pheatmap(V_mark[, order(clusters_hdb_pam50), drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_col = ann,
             main = "PAM50 Gene Expression (k-means, HC, HDBSCAN, HDBSCAN PAM50)")
  })
} else {
  message("[PAM50 heatmap] No PAM50 genes found in V_adj — skipping.")
}

# PAM50/TCGA subtype correlation
if (!is.na(tcga_sub_col)) {
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_lvls <- sort(unique(na.omit(lab_tcga)))
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE))
  cl_means_kmeans <- sapply(levels(clusters_kmeans), function(cn) rowMeans(V_adj[hvgs, clusters_kmeans == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hc <- sapply(levels(clusters_hc), function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hdb <- sapply(levels(clusters_hdb), function(cn) rowMeans(V_adj[hvgs, clusters_hdb == cn, drop = FALSE], na.rm = TRUE))
  cl_means_hdb_pam50 <- sapply(levels(clusters_hdb_pam50), function(cn) rowMeans(M_pam50[, clusters_hdb_pam50 == cn, drop = FALSE], na.rm = TRUE))
  g_common <- intersect(rownames(sub_means), hvgs)
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50))
  corr_kmeans <- cor(cl_means_kmeans[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hc <- cor(cl_means_hc[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hdb <- cor(cl_means_hdb[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman")
  corr_hdb_pam50 <- cor(cl_means_hdb_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman")
  write.csv(corr_kmeans, file.path(outdir, "cor_kmeans_vs_TCGA_subtype.csv"))
  write.csv(corr_hc, file.path(outdir, "cor_hierarchical_vs_TCGA_subtype.csv"))
  write.csv(corr_hdb, file.path(outdir, "cor_hdbscan_vs_TCGA_subtype.csv"))
  write.csv(corr_hdb_pam50, file.path(outdir, "cor_hdbscan_pam50_vs_TCGA_subtype.csv"))
  safe_pdf(file.path(outdir, "heatmap_cluster_vs_TCGA_subtype_corr.pdf"), {
    par(mfrow = c(2, 2))
    pheatmap::pheatmap(corr_kmeans, main = "k-means vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hc, main = "Hierarchical vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hdb, main = "HDBSCAN vs TCGA Subtypes (ρ)")
    pheatmap::pheatmap(corr_hdb_pam50, main = "HDBSCAN PAM50 vs TCGA Subtypes (ρ)")
    par(mfrow = c(1, 1))
  })
}

# Map DSMZ to TCGA
tcga_mean <- rowMeans(V_tcga, na.rm = TRUE)
rho_overall <- sapply(colnames(V_dsmz), function(s) {
  suppressWarnings(cor(V_dsmz[, s], tcga_mean, method = "spearman", use = "pairwise.complete.obs"))
})

rho_by_sub <- NULL
best_tcga_subtype <- NULL
best_tcga_rho <- NULL

if (!is.na(tcga_sub_col)) {
  lab <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  sub_levels <- sort(unique(na.omit(lab)))
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
  write.csv(rho_by_sub, file.path(outdir, "rho_dsmz_to_tcga_subtypes.csv"))
}

map_tbl <- tibble(
  sample = colnames(V_dsmz),
  rho_overall = unname(rho_overall),
  kmeans_cluster = clusters_kmeans[colnames(V_dsmz)],
  hc_subtype = clusters_hc[colnames(V_dsmz)],
  hdb_subtype = clusters_hdb[colnames(V_dsmz)],
  hdb_pam50_subtype = clusters_hdb_pam50[colnames(V_dsmz)]
)
if (!is.null(best_tcga_subtype)) {
  map_tbl$best_tcga_subtype <- best_tcga_subtype[map_tbl$sample]
  map_tbl$best_tcga_subtype_rho <- best_tcga_rho[map_tbl$sample]
}
write.csv(map_tbl, file.path(outdir, "mapping_dsmz_to_tcga.csv"), row.names = FALSE)
print(head(map_tbl))

# Marker expression heatmap
if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark))
  }
  tcga_sub_vec <- rep(NA_character_, ncol(V_tcga))
  if (!is.na(tcga_sub_col)) {
    tcga_sub_vec <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col])
  }
  ann_tcga <- data.frame(
    dataset = "TCGA",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_tcga)]),
    hc_subtype = as.character(clusters_hc[colnames(V_tcga)]),
    hdb_subtype = as.character(clusters_hdb[colnames(V_tcga)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50[colnames(V_tcga)]),
    TCGA_subtype = tcga_sub_vec,
    row.names = colnames(V_tcga),
    check.names = FALSE
  )
  ann_dsmz <- data.frame(
    dataset = "DSMZ",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_dsmz)]),
    hc_subtype = as.character(clusters_hc[colnames(V_dsmz)]),
    hdb_subtype = as.character(clusters_hdb[colnames(V_dsmz)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50[colnames(V_dsmz)]),
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
    pdf(file.path(outdir, "heatmap_pam50_vst.pdf"), width = 12, height = 8)
    pheatmap( V_mark, show_rownames = TRUE, show_colnames = FALSE,
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
DSM_med <- aggregate_median(V_dsmz_z, clusters_kmeans[colnames(V_dsmz)])  # Calculate median per cluster

## ---- safe median-z heatmap ----
common_mark <- intersect(marker_ens, intersect(rownames(TC_med), rownames(DSM_med)))
if (length(common_mark) == 0L) {
  cat("[DEBUG] marker_ens length =", length(marker_ens), "\n")
  cat("[DEBUG] nrow(TC_med) =", nrow(TC_med), " nrow(DSM_med) =", nrow(DSM_med), "\n")
  cat("[DEBUG] Example TC_med genes:", paste(head(rownames(TC_med), 5), collapse=", "), "\n")
  cat("[DEBUG] Example DSM_med genes:", paste(head(rownames(DSM_med), 5), collapse=", "), "\n")
  stop("[FATAL] No overlap between marker and (TC_med ∩ DSM_med). Check ID harmonization (Ensembl, version stripping) and markers.")
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
emb$discovery_cluster <- clusters_kmeans[emb$sample]  # Add discovery clusters

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

