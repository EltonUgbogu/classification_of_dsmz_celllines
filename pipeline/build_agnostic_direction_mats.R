#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(readr)
})

opt_list <- list(
  make_option("--cell_rds",   type="character", help="genes x samples"),
  make_option("--tumour_rds", type="character", help="genes x samples"),
  make_option("--genes",      type="character", help="txt file: Ensembl IDs"),
  make_option("--out_cell",   type="character"),
  make_option("--out_tumour", type="character")
)

opt <- parse_args(OptionParser(option_list = opt_list))

info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep="")

read_mat <- function(path) {
  obj <- readRDS(path)
  if (is.matrix(obj) || inherits(obj, "Matrix")) return(as.matrix(obj))
  if (inherits(obj, "SummarizedExperiment") || inherits(obj, "DESeqTransform")) {
    if (!requireNamespace("SummarizedExperiment", quietly=TRUE))
      stop("SummarizedExperiment not installed.")
    return(SummarizedExperiment::assay(obj))
  }
  stop("Unsupported object type: ", paste(class(obj), collapse=", "))
}

strip_ens_version <- function(x) sub("\\.\\d+$", "", x)

cell <- read_mat(opt$cell_rds)
tum  <- read_mat(opt$tumour_rds)

rownames(cell) <- strip_ens_version(rownames(cell))
rownames(tum)  <- strip_ens_version(rownames(tum))

genes <- readr::read_lines(opt$genes)
genes <- strip_ens_version(genes)
genes <- unique(genes)

keep_cell <- intersect(genes, rownames(cell))
keep_tum  <- intersect(genes, rownames(tum))
keep <- intersect(keep_cell, keep_tum)

if (length(keep) < 5) {
  stop("Too few genes after intersection. keep=", length(keep))
}

info("Gene list size: %d", length(genes))
info("CELL keep: %d | TUMOUR keep: %d | common: %d", length(keep_cell), length(keep_tum), length(keep))

cell_f <- cell[keep, , drop=FALSE]
tum_f  <- tum[keep, , drop=FALSE]

dir.create(dirname(opt$out_cell), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(opt$out_tumour), showWarnings=FALSE, recursive=TRUE)

saveRDS(cell_f, opt$out_cell)
saveRDS(tum_f,  opt$out_tumour)

info("Saved: %s", opt$out_cell)
info("Saved: %s", opt$out_tumour)
