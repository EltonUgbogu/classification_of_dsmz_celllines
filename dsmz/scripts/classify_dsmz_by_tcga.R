#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyverse); library(SummarizedExperiment); library(matrixStats)
  library(limma); library(pheatmap); library(ggrepel); library(uwot)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
set.seed(42)

# ---------------- config (edit) ----------------
TCGA_SE_RDS   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds"
DSMZ_RDS      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
DSMZ_META_CSV <- "/home/chu25/data/dsmz/DSMZ_metadata_with_sample_name.csv"
OUTDIR        <- "/home/chu25/dsmz/results/classify"
VAR_GENE_SET  <- "tcga_5k"     # "tcga_5k", "tcga_10k", "all"
PURITY_TSV    <- NULL          # optional
dir.create(OUTDIR, FALSE, TRUE)

# --------------- load ----------------
tcga_se <- readRDS(TCGA_SE_RDS)
Xtcga   <- assay(tcga_se); rownames(Xtcga) <- sub("\\..*$","",rownames(Xtcga))
if (any(duplicated(rownames(Xtcga)))) Xtcga <- rowsum(Xtcga, rownames(Xtcga), reorder=TRUE)

dsmz_raw  <- readRDS(DSMZ_RDS)
stopifnot(all(c("Ensembl_ID","gene_name") %in% colnames(dsmz_raw)))
scols <- setdiff(colnames(dsmz_raw), c("Ensembl_ID","gene_name","Ensembl_ID_with_version"))
dsmz_raw[scols] <- lapply(dsmz_raw[scols], \(x) as.numeric(as.character(x)))
Xdsmz <- as.matrix(dsmz_raw[, scols, drop=FALSE]); rownames(Xdsmz) <- dsmz_raw$Ensembl_ID
if (any(duplicated(rownames(Xdsmz)))) Xdsmz <- rowsum(Xdsmz, rownames(Xdsmz), reorder=TRUE)

meta <- read.csv(DSMZ_META_CSV) %>% mutate(sample_id = sample_name)
Xdsmz <- Xdsmz[, intersect(colnames(Xdsmz), meta$sample_id), drop=FALSE]
meta  <- meta  %>% filter(sample_id %in% colnames(Xdsmz))
stopifnot(identical(colnames(Xdsmz), meta$sample_id))

# --------------- shared genes ----------------
storage.mode(Xtcga) <- "double"; storage.mode(Xdsmz) <- "double"
g <- intersect(rownames(Xtcga), rownames(Xdsmz))
Xtcga <- Xtcga[g, , drop=FALSE]; Xdsmz <- Xdsmz[g, , drop=FALSE]

# --------------- purity (optional) -----------
purity <- NULL
if ("purity" %in% colnames(colData(tcga_se))) {
  v <- as.numeric(colData(tcga_se)$purity); names(v) <- colnames(Xtcga); purity <- v
} else if (!is.null(PURITY_TSV) && file.exists(PURITY_TSV)) {
  p <- read.delim(PURITY_TSV); purity <- setNames(p$purity, p$sample)
  purity <- purity[colnames(Xtcga)]
}
Xtcga_adj <- Xtcga
if (!is.null(purity)) {
  keep <- names(purity)[!is.na(purity)]
  purity_k <- purity[keep]; Xt <- Xtcga[, keep, drop=FALSE]
  ct <- apply(Xt, 1, \(v) suppressWarnings(cor.test(as.numeric(v), purity_k, method="spearman")))
  padj <- p.adjust(vapply(ct, `[[`, numeric(1), "p.value"), "BH")
  rho  <- vapply(ct, function(z) unname(z$estimate), numeric(1))
  dropg <- names(which(padj < 0.01 & rho < -0.4))
  if (length(dropg)) Xt <- Xt[setdiff(rownames(Xt), dropg), , drop=FALSE]
  inf <- 1 - purity_k; fit <- lmFit(Xt, model.matrix(~inf)); beta <- fit$coefficients[,2,drop=FALSE]
  Xt_adj <- as.matrix(Xt) - beta %*% t(inf)
  Xtcga_adj <- Xtcga; cg <- intersect(rownames(Xt_adj), rownames(Xtcga_adj))
  Xtcga_adj[cg, colnames(Xt_adj)] <- Xt_adj[cg, ]
}

# --------------- labels & cohort means -------
labs <- NULL
for (nm in c("project_id","study","disease_type")) if (nm %in% colnames(colData(tcga_se))) { labs <- colData(tcga_se)[[nm]]; break }
labs <- sub("^TCGA-","", as.character(labs)); names(labs) <- colnames(Xtcga_adj)
cohorts <- sort(unique(labs))
M <- sapply(cohorts, \(ct) rowMeans(Xtcga_adj[, labs==ct, drop=FALSE])); colnames(M) <- cohorts

# --------------- variable genes --------------
sel_genes <- rownames(M)
if (tolower(VAR_GENE_SET) %in% c("tcga_5k","tcga_10k")) {
  iq <- matrixStats::rowIQRs(M); n <- if (tolower(VAR_GENE_SET)=="tcga_5k") 5000 else 10000
  sel_genes <- names(sort(setNames(iq, rownames(M)), decreasing=TRUE))[seq_len(min(n,length(iq)))]
}
guse <- sel_genes

# --------------- Spearman DSMZ x cohort ------
S <- outer(
  meta$sample_id, colnames(M),
  Vectorize(function(smp, ct) suppressWarnings(cor(Xdsmz[guse, smp], M[guse, ct], method="spearman", use="pairwise.complete.obs")))
)
dimnames(S) <- list(meta$sample_id, colnames(M))

# --------------- calls & confidence ----------
top_idx <- max.col(S, ties.method="first")
second <- apply(S, 1, function(v){ v[order(v, decreasing=TRUE)][2] %||% NA_real_ })
calls <- tibble(
  sample_id = rownames(S),
  call      = colnames(S)[top_idx],
  rho_best  = S[cbind(seq_len(nrow(S)), top_idx)],
  rho_2nd   = second,
  margin    = rho_best - rho_2nd
)
write.csv(calls, file.path(OUTDIR,"dsmz_tcga_cohort_calls.csv"), row.names=FALSE)
write.csv(as.data.frame(S) %>% rownames_to_column("sample_id"),
          file.path(OUTDIR,"dsmz_tcga_spearman_wide.csv"), row.names=FALSE)

# --------------- quick plots -----------------
# Heatmap of top 30 highest-confidence DSMZ
top30 <- calls %>% arrange(desc(margin)) %>% slice_head(n=30) %>% pull(sample_id)
pdf(file.path(OUTDIR,"heatmap_top30_dsmz_vs_tcga_cohorts.pdf"), 8,6)
pheatmap(S[top30, , drop=FALSE], cluster_rows=TRUE, cluster_cols=TRUE,
         color=colorRampPalette(c("navy","white","firebrick3"))(50),
         main="DSMZ vs TCGA cohorts (Spearman)")
dev.off()

# UMAP (optional sanity)
X <- cbind(M[guse, , drop=FALSE], Xdsmz[guse, , drop=FALSE]); Xz <- t(scale(t(X))); Xz[is.na(Xz)] <- 0
src  <- c(rep("TCGA", ncol(M)), rep("DSMZ", ncol(Xdsmz)))
lab  <- c(colnames(M), colnames(Xdsmz))
um   <- uwot::umap(Xz, n_neighbors=15, min_dist=0.2, metric="cosine")
umdf <- as.data.frame(um); colnames(umdf) <- c("UMAP1","UMAP2"); umdf$source <- src; umdf$label <- lab
pdf(file.path(OUTDIR,"umap_tcga_means_plus_dsmz.pdf"), 7,5)
ggplot(umdf, aes(UMAP1,UMAP2,color=source,shape=source)) + geom_point(size=1.8, alpha=.85) + theme_bw(11)
dev.off()

cat("[OK] Disease classification complete → ", OUTDIR, "\n")
