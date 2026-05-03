#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args_simple <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    joint_rds = "data/brca/brca_vst_joint.rds",
    genes = "resources/pam50_ensembl_ids.txt",
    metadata = "data/dsmz/DSMZ_metadata.csv",
    out_expr = "results/unsupervised/brca/tumour_neighbourhoods_input/expr_pam50.rds",
    out_map = "results/unsupervised/brca/tumour_neighbourhoods_input/cell_line_to_original_sample_id_pam50.rds"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    out[[gsub("-", "_", substring(key, 3), fixed = TRUE)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
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

opt <- parse_args_simple()

stopifnot(file.exists(opt$joint_rds))
stopifnot(file.exists(opt$genes))
stopifnot(file.exists(opt$metadata))

message("[INFO] Reading joint VST: ", opt$joint_rds)
mat <- read_mat(opt$joint_rds)

if (is.null(rownames(mat))) stop("Joint VST must have gene rownames.")
if (is.null(colnames(mat))) stop("Joint VST must have sample colnames.")

genes <- unique(scan(opt$genes, what = character(), quiet = TRUE))
genes <- intersect(genes, rownames(mat))
message("[INFO] PAM50 genes present in VST: ", length(genes))
if (length(genes) < 40L) {
  stop("Too few PAM50 genes found in joint VST: ", length(genes))
}

expr <- t(mat[genes, , drop = FALSE])
dsmz_sample_ids <- rownames(expr)[grepl("^NG-", rownames(expr))]
tcga_sample_ids <- rownames(expr)[grepl("^TCGA-", rownames(expr))]

if (length(dsmz_sample_ids) == 0L) stop("No DSMZ NG-* samples in joint VST.")
if (length(tcga_sample_ids) == 0L) stop("No TCGA-* samples in joint VST.")

meta <- read.csv(opt$metadata, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_name", "DSMZ_Cell_line_norm")
missing <- setdiff(required, colnames(meta))
if (length(missing)) {
  stop("DSMZ metadata missing required columns: ", paste(missing, collapse = ", "))
}

cell_line_of <- setNames(as.character(meta$DSMZ_Cell_line_norm), as.character(meta$sample_name))
mapping <- cell_line_of[dsmz_sample_ids]
missing_mapping <- is.na(mapping) | !nzchar(trimws(mapping))
if (any(missing_mapping)) {
  warning("Using sample IDs for ", sum(missing_mapping), " DSMZ samples missing metadata mapping.")
  mapping[missing_mapping] <- dsmz_sample_ids[missing_mapping]
}

dir.create(dirname(opt$out_expr), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out_map), recursive = TRUE, showWarnings = FALSE)

saveRDS(expr, opt$out_expr)
saveRDS(mapping, opt$out_map)

message("[INFO] Saved expression: ", opt$out_expr)
message("[INFO] Saved mapping: ", opt$out_map)
message("[INFO] Expression dims: ", nrow(expr), " samples x ", ncol(expr), " genes")
message("[INFO] DSMZ samples: ", length(dsmz_sample_ids), " | TCGA samples: ", length(tcga_sample_ids))
message("[INFO] Mapping entries: ", length(mapping))
