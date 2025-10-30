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
