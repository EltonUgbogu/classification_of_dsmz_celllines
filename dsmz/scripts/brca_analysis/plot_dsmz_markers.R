#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# ---------- CONFIG (edit paths if needed) ----------
outdir <- "/home/chu25/dsmz/results/brca_vst_dual2"
rds_dsmz_vst <- file.path(outdir, "VST_dsmz_brca.rds")
ens2sym_csv  <- NULL  # set to a CSV with columns: ensembl,symbol if you saved one
mapping_csv  <- file.path(outdir, "dsmz_from_tcga_reference/dsmz_pred_subtypes_from_tcga_top20.csv")
clusters_csv <- file.path(outdir, "clusters_kmeans.csv")  # optional (for k-means labels)

markers <- c("ESR1","PGR","ERBB2","KRT5","KRT14","MKI67")
pdf_out <- file.path(outdir, "heatmap_dsmz_markers_standalone.pdf")

# ---------- HELPERS ----------
zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

# ---------- LOAD ----------
stopifnot(file.exists(rds_dsmz_vst))
V_dsmz <- readRDS(rds_dsmz_vst)            # genes (Ensembl) x DSMZ samples

# Build symbol map
ens2sym <- NULL
if (!is.null(ens2sym_csv) && file.exists(ens2sym_csv)) {
  ens2sym <- read.csv(ens2sym_csv, stringsAsFactors = FALSE) %>%
    dplyr::rename(ensembl = 1, symbol = 2) %>%
    dplyr::filter(!is.na(ensembl), nzchar(ensembl), !is.na(symbol), nzchar(symbol)) %>%
    distinct(symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE)
} else {
  # Fallback: derive from rownames if they already look like symbols (only if you saved that way)
  ens2sym <- tibble(ensembl = rownames(V_dsmz), symbol = rownames(V_dsmz))
}
stopifnot(all(c("ensembl","symbol") %in% names(ens2sym)))

# Map markers -> Ensembl present in V_dsmz
sym2ens <- setNames(ens2sym$ensembl, toupper(ens2sym$symbol))
sel_ens <- sym2ens[toupper(markers)]
sel_ens <- sel_ens[!is.na(sel_ens)]
sel_ens <- intersect(sel_ens, rownames(V_dsmz))
if (!length(sel_ens)) stop("None of the requested markers found in the DSMZ matrix.")

# Subset and label rows by symbol (keep original marker order)
M <- V_dsmz[sel_ens, , drop = FALSE]
ens2sym_lookup <- setNames(names(sym2ens), sym2ens)  # ENSG -> SYMBOL
rownames(M) <- ens2sym_lookup[rownames(M)]
M <- M[intersect(markers, rownames(M)), , drop = FALSE]

# Z-score by gene for readability
suppressPackageStartupMessages(library(matrixStats))
Mz <- zscore_by_gene(M)

# Column annotations (optional)
ann <- data.frame(row.names = colnames(Mz))
if (file.exists(mapping_csv)) {
  dsmz_assign <- read.csv(mapping_csv, stringsAsFactors = FALSE)
  pred_map <- setNames(as.character(dsmz_assign$label), dsmz_assign$sample)
  ann$Predicted <- pred_map[colnames(Mz)]
}
if (file.exists(clusters_csv)) {
  km <- read.csv(clusters_csv, stringsAsFactors = FALSE)
  km_map <- setNames(as.character(km$cluster), km$sample)
  ann$Kmeans <- km_map[colnames(Mz)]
}
if (ncol(ann) == 0) ann <- NULL

# (nice-to-have) order samples by predicted subtype, then k-means
ord <- seq_len(ncol(Mz))
if (!is.null(ann)) {
  ord <- order(ann$Predicted, ann$Kmeans, na.last = TRUE)
}
Mz <- Mz[, ord, drop = FALSE]
if (!is.null(ann)) ann <- ann[ord, , drop = FALSE]

# ---------- PLOT ----------
dir.create(dirname(pdf_out), showWarnings = FALSE, recursive = TRUE)
pdf(pdf_out, width = 9, height = 5.5)
pheatmap(Mz,
         scale = "none",
         clustering_method = "average",
         show_colnames = FALSE,
         main = "DSMZ marker panel (z-scored)",
         annotation_col = ann)
dev.off()
cat("[plot] wrote:", pdf_out, "\n")

# Optional: also write PNG
png(sub("\\.pdf$",".png", pdf_out), width = 1100, height = 700, res = 150)
pheatmap(Mz, scale="none", clustering_method="average", show_colnames=FALSE,
         main="DSMZ marker panel (z-scored)", annotation_col=ann)
dev.off()
cat("[plot] wrote:", sub("\\.pdf$",".png", pdf_out), "\n")
