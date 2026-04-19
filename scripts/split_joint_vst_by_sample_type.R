#!/usr/bin/env Rscript

parse_args_simple <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opt <- parse_args_simple()

if (is.null(opt$joint_rds) || is.null(opt$out_cell) || is.null(opt$out_tumour)) {
  stop("Required arguments: --joint_rds <file> --out_cell <file> --out_tumour <file>")
}

read_mat <- function(path) {
  obj <- readRDS(path)
  if (is.matrix(obj)) return(obj)
  if (inherits(obj, "Matrix")) return(as.matrix(obj))
  if (inherits(obj, "SummarizedExperiment") || inherits(obj, "DESeqTransform")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment not installed.")
    }
    return(SummarizedExperiment::assay(obj))
  }
  stop("Unsupported object type: ", paste(class(obj), collapse = ", "))
}

message("[INFO] Reading joint VST: ", opt$joint_rds)
mat <- read_mat(opt$joint_rds)

if (is.null(colnames(mat))) {
  stop("Joint VST matrix must have sample column names.")
}

sample_ids <- colnames(mat)
cell_idx <- grepl("^NG-", sample_ids)
tumour_idx <- grepl("^(GSE|SRP|TARGET-)", sample_ids)

if (sum(cell_idx) == 0) {
  stop("No cell-line columns matched '^NG-'.")
}
if (sum(tumour_idx) == 0) {
  stop("No tumour columns matched '^(GSE|SRP|TARGET-)'.")
}

overlap_idx <- cell_idx & tumour_idx
if (any(overlap_idx)) {
  stop("Some columns matched both cell and tumour patterns: ",
       paste(sample_ids[overlap_idx], collapse = ", "))
}

unmatched_idx <- !(cell_idx | tumour_idx)
if (any(unmatched_idx)) {
  warning(
    "Ignoring unmatched columns: ",
    paste(sample_ids[unmatched_idx], collapse = ", ")
  )
}

cell_mat <- mat[, cell_idx, drop = FALSE]
tumour_mat <- mat[, tumour_idx, drop = FALSE]

dir.create(dirname(opt$out_cell), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_tumour), recursive = TRUE, showWarnings = FALSE)

saveRDS(cell_mat, opt$out_cell)
saveRDS(tumour_mat, opt$out_tumour)

message("[INFO] Cell matrix dims: ", nrow(cell_mat), " genes x ", ncol(cell_mat), " samples")
message("[INFO] Tumour matrix dims: ", nrow(tumour_mat), " genes x ", ncol(tumour_mat), " samples")
message("[INFO] Saved: ", opt$out_cell)
message("[INFO] Saved: ", opt$out_tumour)
