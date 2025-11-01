## ---------- Library paths & threading ----------
.libPaths(c("/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library", .libPaths()))
Sys.setenv(R_LIBS_USER = "/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")
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
  library(dynamicTreeCut)  # cutreeDynamic / cutreeHybrid
  library(WGCNA)           # plotDendroAndColors
  library(paran)           # Parallel Analysis
  # Diagnostics:
  library(vegan)           # adonis2
  library(randomForest)    # randomForest
})

## ---------- CONFIG ----------
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
tcga_sub_col <- "PAM50.Subtype"  # if present in colData(tcga_se)

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
  if (is.null(pr$x) || ncol(pr$x) < 1) stop("[FATAL] PCA returned 0 components.")
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

# ---- subtype sanity info ----
if (!is.null(tcga_sub_col) && tcga_sub_col %in% colnames(colData(tcga_se))) {
  tcga_samples <- colnames(tcga_counts)
  tcga_subtype_df <- data.frame(
    Sample = tcga_samples,
    Subtype = as.character(colData(tcga_se)[tcga_samples, tcga_sub_col])
  ) %>% filter(!is.na(Subtype))
  common_subtypes <- unique(tcga_subtype_df$Subtype)
  cat(sprintf("[INFO] TCGA samples have %d subtypes: %s\n", length(common_subtypes), paste(common_subtypes, collapse=", ")))
}

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

stopifnot(ncol(tcga_counts) > 0L, ncol(dsmz_counts) > 0L, nrow(Xc_raw) > 0L, ncol(Xc_raw) > 0L)

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

## ---- Independent VST per dataset ----
cat("[INFO] VST per-dataset (TCGA only)...\n")
dds_tcga_ind <- DESeqDataSetFromMatrix(round(tcga_counts), DataFrame(row.names = colnames(tcga_counts)), design = ~ 1)
# (Optionally estimate dispersions here if you plan to plot them)
vsd_tcga_ind <- vst(dds_tcga_ind, blind = TRUE)
V_tcga_ind   <- assay(vsd_tcga_ind)
saveRDS(V_tcga_ind, file.path(outdir, "VST_tcga_only.rds"))

cat("[INFO] VST per-dataset (DSMZ only)...\n")
dds_dsmz_ind <- DESeqDataSetFromMatrix(round(dsmz_counts), DataFrame(row.names = colnames(dsmz_counts)), design = ~ 1)
vsd_dsmz_ind <- vst(dds_dsmz_ind, blind = TRUE)
V_dsmz_ind   <- assay(vsd_dsmz_ind)
saveRDS(V_dsmz_ind, file.path(outdir, "VST_dsmz_only.rds"))

## ---- Purity adjust TCGA (on VST) ----
purity <- load_purity_data(tcga_se, NULL)
if (!is.null(purity)) {
  cat("[INFO] Applying purity adjustment to TCGA VST data...\n")
  adj <- purity_adjust_log(V_tcga_ind, purity)
  V_tcga_ind <- adj$V
  writeLines(adj$dropped, file.path(outdir, "purity_associated_genes.txt"))
  if (drop_purity_genes_permanently && length(adj$dropped)) {
    V_tcga_ind <- V_tcga_ind[setdiff(rownames(V_tcga_ind), adj$dropped), , drop = FALSE]
    dsmz_counts <- dsmz_counts[setdiff(rownames(dsmz_counts), adj$dropped), , drop = FALSE]
    message(sprintf("[INFO] Removed %d purity-associated genes from both TCGA VST and DSMZ counts.", length(adj$dropped)))
  }
}

## ---- Merge normalized (VST) prior to ComBat_seq (for diagnostics) ----
common_genes <- intersect(rownames(V_tcga_ind), rownames(V_dsmz_ind))
V_tcga_norm <- V_tcga_ind[common_genes, , drop = FALSE]
V_dsmz_norm <- V_dsmz_ind[common_genes, , drop = FALSE]
V_merged    <- cbind(V_tcga_norm, V_dsmz_norm)
log_dims("Merged normalized", V_merged)

## ---- ComBat_seq on raw counts (aligned order!) ----
cat("[INFO] Applying ComBat_seq to raw counts...\n")
all_counts <- Xc_raw[common_genes, , drop = FALSE]

# Build batch vector in the SAME column order as all_counts
batch <- factor(ifelse(colnames(all_counts) %in% colnames(V_tcga_norm), "TCGA", "DSMZ"),
                levels = c("TCGA","DSMZ"))

set.seed(42)
all_counts_combat <- sva::ComBat_seq(
  counts = as.matrix(all_counts),
  batch  = batch,
  group  = NULL,
  covar_mod = NULL,
  full_mod = TRUE,
  shrink = FALSE,
  shrink.disp = FALSE
)

# VST after ComBat_seq
cat("[INFO] Converting ComBat_seq corrected counts to VST...\n")
coldata2 <- tibble(sample = colnames(all_counts_combat),
                   dataset = as.character(batch)) %>%
            column_to_rownames("sample")
dds_combat <- DESeqDataSetFromMatrix(round(all_counts_combat), DataFrame(coldata2), design = ~ 1)
vsd_combat <- vst(dds_combat, blind = TRUE)
V_adj <- assay(vsd_combat)

# Split
V_tcga <- V_adj[, colnames(V_tcga_norm), drop = FALSE]
V_dsmz <- V_adj[, colnames(V_dsmz_norm), drop = FALSE]

log_dims("TCGA (ComBat_seq corrected)", V_tcga)
log_dims("DSMZ (ComBat_seq corrected)", V_dsmz)

## ---- Diagnostics ----
run_all_diagnostics <- function(M, V_tcga_norm, name_tag, outdir, tcga_se, tcga_sub_col) {
  message(sprintf("[DIAG] Running diagnostics for %s data...", name_tag))
  # dataset labels relative to TCGA sample list
  batch_labels <- factor(ifelse(colnames(M) %in% colnames(V_tcga_norm), "TCGA", "DSMZ"),
                         levels = c("TCGA","DSMZ"))
  meta_df <- data.frame(sample = colnames(M), dataset = batch_labels, stringsAsFactors = FALSE)

  # Add subtype for TCGA columns only
  meta_df$subtype <- NA_character_
  if (!is.null(tcga_sub_col) && tcga_sub_col %in% colnames(colData(tcga_se))) {
    tcga_samples <- meta_df$sample[meta_df$dataset == "TCGA"]
    meta_df$subtype[meta_df$dataset == "TCGA"] <- as.character(colData(tcga_se)[tcga_samples, tcga_sub_col])
  }
  meta_df$subtype <- factor(meta_df$subtype)

  # PCA on HVGs
  pcs <- make_pcs(M, n_hvg = 3000, max_pc = 30)
  kpc <- min(6, ncol(pcs$PC))

  # PC regression: fraction of variance explained by dataset / subtype
  pc_reg_res <- lapply(seq_len(kpc), function(i) {
    df <- data.frame(PC = pcs$PC[, i], dataset = meta_df$dataset, subtype = meta_df$subtype)
    fit <- lm(PC ~ dataset + subtype, data = df)
    an <- anova(fit)
    total_ss <- sum(an[, "Sum Sq"], na.rm = TRUE)
    ss_dataset <- if ("dataset" %in% rownames(an)) an["dataset", "Sum Sq"] else 0
    ss_subtype <- if ("subtype" %in% rownames(an)) an["subtype", "Sum Sq"] else 0
    data.frame(PC = paste0("PC", i),
               frac_dataset = as.numeric(ss_dataset) / total_ss,
               frac_subtype = as.numeric(ss_subtype) / total_ss)
  })
  pc_reg_df <- do.call(rbind, pc_reg_res)
  write.csv(pc_reg_df, file.path(outdir, sprintf("pc_regression_%s_batch.csv", name_tag)), row.names = FALSE)

  # PERMANOVA (Euclidean on HVGs)
  hvgs <- rownames(pcs$pr$rotation)
  mat_for_dist <- t(M[hvgs, , drop = FALSE])
  dmat <- dist(mat_for_dist, method = "euclidean")
  adonis_res <- vegan::adonis2(dmat ~ dataset + subtype, data = meta_df, permutations = 999)
  capture.output(adonis_res, file = file.path(outdir, sprintf("adonis2_%s_dataset_vs_subtype.txt", name_tag)))
  adonis_r2 <- as.numeric(adonis_res["dataset", "R2"])

  # RF classify dataset from top PCs (overall accuracy)
  pc_df_rf <- as.data.frame(pcs$PC[, 1:min(10, ncol(pcs$PC)), drop = FALSE])
  pc_df_rf$dataset <- meta_df$dataset
  set.seed(42)
  rf <- randomForest(dataset ~ ., data = pc_df_rf, ntree = 500)
  rf_acc <- sum(diag(rf$confusion)) / sum(rf$confusion)

  # Per-gene ANOVA for dataset effect
  genes_to_test <- rownames(M)
  pvals <- vapply(genes_to_test, function(g) {
    y <- as.numeric(M[g, ])
    df <- data.frame(y = y, dataset = meta_df$dataset, subtype = meta_df$subtype)
    fit <- try(lm(y ~ dataset + subtype, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    an <- try(anova(fit), silent = TRUE)
    if (inherits(an, "try-error") || !("dataset" %in% rownames(an))) return(NA_real_)
    an["dataset", "Pr(>F)"]
  }, numeric(1))
  padj <- p.adjust(pvals, method = "BH")
  n_batch_genes <- sum(padj < 0.05, na.rm = TRUE)

  list(
    pc1_frac_dataset = pc_reg_df$frac_dataset[1],
    pc1_frac_subtype = pc_reg_df$frac_subtype[1],
    pc_reg_df = pc_reg_df,
    adonis_R2 = adonis_r2,
    rf_acc = rf_acc,
    n_batch_genes = n_batch_genes
  )
}

# Pre vs Post diagnostics
tcga_sub_col_val <- if (exists("tcga_sub_col")) tcga_sub_col else NULL
diag_pre  <- run_all_diagnostics(V_merged, V_tcga_norm, "pre",  outdir, tcga_se, tcga_sub_col_val)
diag_post <- run_all_diagnostics(V_adj,    V_tcga_norm, "post", outdir, tcga_se, tcga_sub_col_val)

summary_df <- data.frame(
  metric = c("PC1_Frac_Dataset","PC1_Frac_Subtype","Adonis_R2_Dataset","RF_Accuracy_Dataset","N_Batch_Genes_FDR05"),
  pre    = c(diag_pre$pc1_frac_dataset,  diag_pre$pc1_frac_subtype,  diag_pre$adonis_R2, diag_pre$rf_acc, diag_pre$n_batch_genes),
  post   = c(diag_post$pc1_frac_dataset, diag_post$pc1_frac_subtype, diag_post$adonis_R2, diag_post$rf_acc, diag_post$n_batch_genes)
)
write.csv(summary_df, file.path(outdir, "pre_post_batch_correction_summary.csv"), row.names = FALSE)
cat(sprintf("[OUTPUT] Batch correction summary saved to: %s\n",
            file.path(outdir, "pre_post_batch_correction_summary.csv")))

## ---- PCA plots (pre vs post) ----
cat("[INFO] PCA before and after batch correction...\n")
top_var_genes <- function(mat, n = 5000) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}
make_pca_plot <- function(M, lab_batch, title, path_pdf, n_top = 5000) {
  sel <- top_var_genes(M, n = min(n_top, nrow(M)))
  X   <- t(scale(M[sel, , drop = FALSE]))
  pr  <- prcomp(X, center = FALSE, scale. = FALSE)
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  df <- data.frame(PC1 = pr$x[,1], PC2 = pr$x[,2],
                   batch = lab_batch, sample = rownames(pr$x), stringsAsFactors = FALSE)
  safe_pdf(path_pdf, {
    layout(matrix(c(1,2), nrow = 1))
    cols <- as.integer(df$batch)
    plot(df$PC1, df$PC2, pch = 19, cex = 0.7, col = cols,
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title)
    legend("topright", legend = levels(df$batch), col = seq_along(levels(df$batch)), pch = 19, cex = 0.8)
    barplot(var_expl[1:20], las = 2, main = "Explained variance (top 20 PCs)",
            ylab = "% variance", xlab = "PC")
    layout(1)
  })
  invisible(list(scores = df, var_expl = var_expl))
}

batch_labels <- factor(ifelse(colnames(V_adj) %in% colnames(V_tcga_norm), "TCGA", "DSMZ"),
                       levels = c("TCGA","DSMZ"))
pca_pre  <- make_pca_plot(V_merged, batch_labels, "PCA (pre-ComBat_seq)",  file.path(outdir, "pca_pre_batch.pdf"))
pca_post <- make_pca_plot(V_adj,    batch_labels, "PCA (post-ComBat_seq)", file.path(outdir, "pca_post_batch.pdf"))
