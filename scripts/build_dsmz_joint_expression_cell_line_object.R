#!/usr/bin/env Rscript
# Build a full-expression DSMZ cell-line object using existing VST inputs.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--dsmz-vst", type="character", default=NULL,
              help="[REQUIRED] BRCA/NBL/RBL DSMZ cell-line VST matrix RDS"),
  make_option("--dsmz-meta", type="character", default=NULL,
              help="[REQUIRED] Metadata TSV for BRCA/NBL/RBL DSMZ cell lines"),
  make_option("--heme-vst", type="character", default=NULL,
              help="[REQUIRED] HEME DSMZ cell-line VST matrix RDS"),
  make_option("--output", type="character", default=NULL,
              help="[REQUIRED] Output RDS with expr, meta, genes"),
  make_option("--output-metadata", type="character", default=NULL,
              help="Optional metadata TSV export")
)
opt <- parse_args(OptionParser(option_list=option_list))

required <- c("dsmz-vst", "dsmz-meta", "heme-vst", "output")
missing <- required[vapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required argument(s): ", paste(missing, collapse=", "))
}
for (path in c(opt[["dsmz-vst"]], opt[["dsmz-meta"]], opt[["heme-vst"]])) {
  if (!file.exists(path)) stop("Input file not found: ", path)
}
dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
if (!is.null(opt[["output-metadata"]]) && nzchar(opt[["output-metadata"]])) {
  dir.create(dirname(opt[["output-metadata"]]), recursive=TRUE, showWarnings=FALSE)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

load_matrix <- function(path, label) {
  obj <- readRDS(path)
  if (is.list(obj) && "expr" %in% names(obj)) {
    obj <- obj$expr
  }
  if (!(is.matrix(obj) || is.data.frame(obj))) {
    stop(label, " input is not a matrix/data.frame or list containing $expr: ", path)
  }
  mat <- as.matrix(obj)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- strip_version(rownames(mat))
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    stop(label, " matrix has missing sample identifiers")
  }
  mat
}

normalise_type <- function(x) {
  y <- tolower(gsub("[ -]+", "_", as.character(x)))
  y[y %in% c("cell", "cells", "cellline", "cell_line")] <- "cell_line"
  y[y %in% c("tumor", "tumour", "tumors", "tumours")] <- "tumour"
  y
}

dsmz_expr <- load_matrix(opt[["dsmz-vst"]], "DSMZ")
heme_expr <- load_matrix(opt[["heme-vst"]], "HEME")

dsmz_meta <- fread(opt[["dsmz-meta"]])
if (!all(c("sample_id", "cancer_type", "sample_type") %in% names(dsmz_meta))) {
  stop("DSMZ metadata must contain sample_id, cancer_type, and sample_type columns")
}
dsmz_meta <- unique(dsmz_meta, by="sample_id")
dsmz_meta <- dsmz_meta[sample_id %in% colnames(dsmz_expr)]
dsmz_meta[, type := normalise_type(sample_type)]
dsmz_meta <- dsmz_meta[type == "cell_line"]
if (nrow(dsmz_meta) == 0L) {
  stop("No BRCA/NBL/RBL DSMZ cell-line samples found after metadata filtering")
}
dsmz_ids <- dsmz_meta$sample_id
dsmz_expr <- dsmz_expr[, dsmz_ids, drop=FALSE]
dsmz_meta[, `:=`(
  lineage = toupper(cancer_type),
  cancer_type = toupper(cancer_type),
  sample_type = "cell_line",
  source_profile = toupper(cancer_type)
)]
dsmz_meta <- dsmz_meta[, .(
  sample_id, lineage, cancer_type, type, sample_type, source_profile
)]

heme_ids <- colnames(heme_expr)
heme_meta <- data.table(
  sample_id = heme_ids,
  lineage = "HEME",
  cancer_type = "HEME",
  type = "cell_line",
  sample_type = "cell_line",
  source_profile = "HEME"
)

common_genes <- intersect(rownames(dsmz_expr), rownames(heme_expr))
common_genes <- common_genes[nzchar(common_genes)]
if (length(common_genes) == 0L) {
  stop("No common genes found between DSMZ and HEME VST matrices")
}
common_genes <- rownames(dsmz_expr)[rownames(dsmz_expr) %in% common_genes]

expr <- cbind(
  dsmz_expr[common_genes, , drop=FALSE],
  heme_expr[common_genes, , drop=FALSE]
)
meta <- rbindlist(list(dsmz_meta, heme_meta), use.names=TRUE)
meta <- meta[match(colnames(expr), sample_id)]
if (any(is.na(meta$sample_id))) {
  stop("Metadata is missing for one or more expression columns")
}
if (anyDuplicated(meta$sample_id)) {
  stop("Duplicate sample identifiers in combined metadata")
}

obj <- list(expr=expr, meta=meta, genes=common_genes)
saveRDS(obj, opt$output)
if (!is.null(opt[["output-metadata"]]) && nzchar(opt[["output-metadata"]])) {
  fwrite(meta, opt[["output-metadata"]], sep="\t")
}

cat("[OK] Full-expression DSMZ cell-line object written to:", opt$output, "\n")
cat("  Genes:", nrow(expr), "\n")
cat("  Samples:", ncol(expr), "\n")
cat("  Cancer-type distribution:\n")
print(table(meta$cancer_type))
