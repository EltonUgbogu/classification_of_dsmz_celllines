## ---------- Library paths & threading ----------
.libPaths(c("/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library", .libPaths()))
Sys.setenv(R_LIBS_USER = "/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")

# Set mc.cores from SLURM (fallback to 1)
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
message(sprintf("[INFO] Using %d threads (SLURM_CPUS_PER_TASK)", getOption("mc.cores", 1L)))

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
  # edgeR not attached; we call edgeR::cpm explicitly later
})

## ---------- CONFIG ----------
tcga_se_rds   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"
dsmz_rds      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
dsmz_meta_csv <- "/home/chu25/data/dsmz/DSMZ_metadata.csv"
outdir <- "/home/chu25/dsmz/results/brca_vst_dual4"
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
drop_purity_genes_permanently <- TRUE

## ---------- HELPERS ----------
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

choose_minPts <- function(n, frac = 0.02, min_floor = 5, max_cap = 50) {
  p <- max(min_floor, round(frac * n))
  p <- min(p, max_cap, n - 1L)
  as.integer(p)
}
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
  v <- matrixStats::rowVars(M); keep <- which(is.finite(v) & v > 0)
  if (length(keep) < 2) stop("[FATAL] Too few variable genes for PCA.")
  M <- M[keep, , drop = FALSE]

  hvgs_all <- names(sort(matrixStats::rowIQRs(M), decreasing = TRUE))
  hvgs <- hvgs_all[seq_len(min(n_hvg, length(hvgs_all)))]

  X <- t(scale(M[hvgs, , drop = FALSE]))
  X <- X[, colSums(!is.na(t(X))) > 0, drop = FALSE]
  if (nrow(X) < 2 || ncol(X) < 1) stop("[FATAL] PCA input has insufficient dimension.")

  pr <- prcomp(X, center = FALSE, scale. = FALSE)
  if (is.null(pr$x) || ncol(pr$x) < 1) stop("[FATAL] PCA returned 0 components. Check input variance/NA.")
  kPC <- as.integer(min(max_pc, ncol(pr$x)))
  list(PC = pr$x[, 1:kPC, drop = FALSE], kPC = kPC, pr = pr, hvgs = hvgs)
}

zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

safe_pdf <- function(path, expr) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path)
  on.exit({ if (dev.cur() != 1) dev.off() }, add = TRUE)
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

  idx_map <- which(ensembl_ids %in% rownames(ens2sym))
  if (!length(idx_map)) stop("No Ensembl IDs mapped to symbols. Check inputs.")

  sym <- ens2sym[ensembl_ids[idx_map], "symbol"]
  hit <- tolower(sym) %in% tolower(pam_symbols)
  keep_idx <- idx_map[hit]
  if (!length(keep_idx)) stop("No PAM50 symbols matched after mapping.")

  sub <- mat[keep_idx, , drop = FALSE]
  sym_keep <- ens2sym[strip_ensver(rownames(sub)), "symbol"]
  split_idx <- split(seq_len(nrow(sub)), tolower(sym_keep))
  collapsed <- lapply(split_idx, function(ix) {
    if (cf == "sum")       colSums(sub[ix, , drop = FALSE], na.rm = TRUE)
    else if (cf == "mean") colMeans(sub[ix, , drop = FALSE], na.rm = TRUE)
    else                   apply(sub[ix, , drop = FALSE], 2, max, na.rm = TRUE)
  })
  M <- do.call(rbind, collapsed)
  rownames(M) <- toupper(names(split_idx))
  M
}

## ---------- LOAD + FILTER ----------
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

# Save snapshots AFTER project/organ filters, BEFORE prefilter
tcga_after_project <- tcga_counts
dsmz_after_organ <- dsmz_counts

common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
cat(sprintf("[INFO] Shared genes: %d\n", length(common)))
if (length(common) < 1000) cat("[WARN] Few shared genes; check Ensembl IDs\n")

tcga_counts <- tcga_counts[common, , drop = FALSE]
dsmz_counts <- dsmz_counts[common, , drop = FALSE]

log_dims("TCGA (post-intersect)", tcga_counts)
log_dims("DSMZ (post-intersect)", dsmz_counts)

# Note: Xc_raw will be created after prefilter below

# ---------- MEMORY-SAFE PREFILTERING (Reduce from ~60k to ~20k genes) ----------
cat("\n[INFO] Prefiltering low-count genes to prevent OOM...\n")
# Prefilter: keep genes with enough counts in enough samples
# TCGA: at least 10 counts in at least 10 samples
# DSMZ: at least 10 counts in at least 5 samples (smaller cohort)
keep_tcga <- rowSums(tcga_counts >= 10) >= 10
keep_dsmz <- rowSums(dsmz_counts >= 10) >= 5

cat(sprintf("[INFO] Prefilter TCGA genes: %d -> %d\n", nrow(tcga_counts_raw), sum(keep_tcga)))
cat(sprintf("[INFO] Prefilter DSMZ genes: %d -> %d\n", nrow(dsmz_counts_raw), sum(keep_dsmz)))

tcga_counts <- tcga_counts[keep_tcga, , drop = FALSE]
dsmz_counts <- dsmz_counts[keep_dsmz, , drop = FALSE]

log_dims("TCGA (after prefilter)", tcga_counts)
log_dims("DSMZ (after prefilter)", dsmz_counts)

# Re-intersect after prefilter
common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
tcga_counts <- tcga_counts[common, , drop = FALSE]
dsmz_counts <- dsmz_counts[common, , drop = FALSE]
log_dims("TCGA (prefilter+intersect)", tcga_counts)
log_dims("DSMZ (prefilter+intersect)", dsmz_counts)

if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA subset has 0 samples.")
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ subset has 0 samples.")
if (nrow(tcga_counts) == 0L) stop("[FATAL] No genes after prefilter intersection.")
if (nrow(dsmz_counts) == 0L) stop("[FATAL] No genes after prefilter intersection.")

# Recreate Xc_raw after prefilter for accurate reporting
Xc_raw_final <- cbind(tcga_counts, dsmz_counts)

dim_report <- file.path(outdir, "dimension_report.txt")
con <- file(dim_report, open = "wt")
on.exit(close(con), add = TRUE)

write_dims_line(con, "TCGA (raw)", tcga_counts_raw)
write_dims_line(con, "DSMZ (raw)", dsmz_counts_raw)
write_dims_line(con, "TCGA (after project)", tcga_after_project)
write_dims_line(con, "DSMZ (after organ)", dsmz_after_organ)
writeLines(sprintf("Shared genes (pre-PF)\t%d", length(intersect(rownames(tcga_after_project), 
                                                                 rownames(dsmz_after_organ)))), con)
write_dims_line(con, "TCGA (after prefilter+intersect)", tcga_counts)
write_dims_line(con, "DSMZ (after prefilter+intersect)", dsmz_counts)
write_dims_line(con, "Merged (final)", Xc_raw_final)

dir.create(file.path(outdir, "shared_genes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "unshared_genes"), showWarnings = FALSE, recursive = TRUE)
writeLines(common, file.path(outdir, "shared_genes", "shared_genes.txt"))
writeLines(setdiff(rownames(tcga_counts_raw), common), file.path(outdir, "unshared_genes", "tcga_unshared_genes.txt"))
writeLines(setdiff(rownames(dsmz_counts_raw), common), file.path(outdir, "unshared_genes", "dsmz_unshared_genes.txt"))

## ---- Independent VST per dataset (Memory-Safe with Parallel) ----
# Helper to check if function argument exists (version compatibility)
has_arg <- function(fun, arg) arg %in% names(formals(fun))

# Setup parallel processing
suppressPackageStartupMessages({ library(BiocParallel) })
nc <- getOption("mc.cores", 1L)
if (nc > 1L) {
  register(MulticoreParam(workers = nc), default = TRUE)
  cat(sprintf("[INFO] Using %d cores for parallel processing\n", nc))
} else {
  cat("[INFO] Using sequential processing\n")
}

cat("[INFO] VST per-dataset (TCGA only)...\n")
dds_tcga_ind <- DESeqDataSetFromMatrix(
  round(tcga_counts),
  data.frame(row.names = colnames(tcga_counts)),
  design = ~ 1
)
# Only estimate size factors; vst() will handle dispersions efficiently
dds_tcga_ind <- estimateSizeFactors(dds_tcga_ind)

# VST with parametric fit and optional nsub for memory reduction
if (has_arg(DESeq2::vst, "nsub")) {
  vsd_tcga_ind <- vst(dds_tcga_ind, blind = TRUE, fitType = "parametric", nsub = 10000L)
} else {
  vsd_tcga_ind <- vst(dds_tcga_ind, blind = TRUE, fitType = "parametric")
}
V_tcga_ind   <- assay(vsd_tcga_ind)
saveRDS(V_tcga_ind, file.path(outdir, "VST_tcga_only.rds"))
# Note: Keep dds_tcga_ind for dispersion plots below

cat("[INFO] VST per-dataset (DSMZ only)...\n")
dds_dsmz_ind <- DESeqDataSetFromMatrix(
  round(dsmz_counts),
  data.frame(row.names = colnames(dsmz_counts)),
  design = ~ 1
)
dds_dsmz_ind <- estimateSizeFactors(dds_dsmz_ind)

# VST with parametric fit and optional nsub for memory reduction
if (has_arg(DESeq2::vst, "nsub")) {
  vsd_dsmz_ind <- vst(dds_dsmz_ind, blind = TRUE, fitType = "parametric", nsub = 10000L)
} else {
  vsd_dsmz_ind <- vst(dds_dsmz_ind, blind = TRUE, fitType = "parametric")
}
V_dsmz_ind   <- assay(vsd_dsmz_ind)
saveRDS(V_dsmz_ind, file.path(outdir, "VST_dsmz_only.rds"))
# Note: Keep dds_dsmz_ind for dispersion plots below

## ---------- LIGHTWEIGHT DISPERSION PLOTS (Memory-Safe) ----------
# Function to make lightweight dispersion plot using subsampled genes
make_dispersion_plot <- function(dds, title, out_pdf, n_genes = 10000L) {
  # Subsample genes to avoid OOM during dispersion estimation/plotting
  set.seed(42)
  g <- nrow(dds)
  pick <- if (g > n_genes) sample.int(g, n_genes) else seq_len(g)
  dds_sub <- dds[pick, ]
  
  cat(sprintf("  Estimating dispersions on %d genes (subsampled from %d)...\n", 
              length(pick), g))
  
  # Pass BPPARAM only if supported (newer DESeq2 versions)
  if (has_arg(DESeq2::estimateDispersions, "BPPARAM")) {
    dds_sub <- DESeq2::estimateDispersions(dds_sub, fitType = "parametric", BPPARAM = bpparam())
  } else {
    dds_sub <- DESeq2::estimateDispersions(dds_sub, fitType = "parametric")
  }
  
  safe_pdf(out_pdf, {
    DESeq2::plotDispEsts(dds_sub,
                         xlab = "Mean of normalized counts (log scale)",
                         ylab = "Dispersion",
                         main = sprintf("%s: Dispersion (%d genes)", title, length(pick))
    )
  })
  cat(sprintf("[OUTPUT] %s\n", out_pdf))
}

cat("\n[INFO] Generating dispersion plots (subsampled to prevent OOM)...\n")
make_dispersion_plot(dds_tcga_ind, "TCGA-BRCA",
                     file.path(outdir, "TCGA_Dispersion_Plot.pdf"), n_genes = 10000L)
make_dispersion_plot(dds_dsmz_ind, "DSMZ (Breast)",
                     file.path(outdir, "DSMZ_Dispersion_Plot.pdf"), n_genes = 10000L)

# Cleanup DDS objects after dispersion plots
rm(dds_tcga_ind, dds_dsmz_ind, vsd_tcga_ind, vsd_dsmz_ind); gc()

## ---------- VST EFFECT CONFIRMATION ----------
cat("\n[INFO] Generating VST effect confirmation plots...\n")

plot_vst_effect <- function(V_mat, title) {
  mean_vst <- rowMeans(V_mat)
  sd_vst   <- matrixStats::rowSds(V_mat)
  plot(mean_vst, sd_vst,
       xlab = "Mean of VST-transformed counts",
       ylab = "Row SD (VST)",
       main = paste(title, "VST: Mean vs SD"),
       pch = 19, cex = 0.5)
  lines(lowess(mean_vst, sd_vst), lwd = 2)
}

tcga_confirm_path <- file.path(outdir, "TCGA_VST_Effect_Confirmation.pdf")
safe_pdf(tcga_confirm_path, {
  plot_vst_effect(V_tcga_ind, "TCGA-BRCA")
})
cat(sprintf("[OUTPUT] TCGA VST Effect Confirmation: %s\n", tcga_confirm_path))

dsmz_confirm_path <- file.path(outdir, "DSMZ_VST_Effect_Confirmation.pdf")
safe_pdf(dsmz_confirm_path, {
  plot_vst_effect(V_dsmz_ind, "DSMZ (Breast)")
})
cat(sprintf("[OUTPUT] DSMZ VST Effect Confirmation: %s\n", dsmz_confirm_path))

cat("\n[INFO] Mean-Variance/Dispersion plots complete!\n")

# Save session information for reproducibility
cat("\n[INFO] Saving session information...\n")
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n[INFO] Pipeline complete! Output directory: ", outdir, "\n")
