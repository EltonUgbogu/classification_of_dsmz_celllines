#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# -----------------------------
# CLI
# -----------------------------
option_list <- list(
  make_option("--expr_rds", type="character", help="Joint expression RDS (genes x samples OR samples x genes)"),
  make_option("--meta_tsv", type="character", help="Joint metadata TSV (must allow cell vs tumour split)"),
  make_option("--direction", type="character", help="e.g. MAD_cosine / Variance_euclidean"),
  make_option("--kind", type="character", help="hc or kmeans (passed by Snakemake)"),
  make_option("--view", type="character", help="cell | tumour | cell_tumour"),
  make_option("--out_rds", type="character", help="Output RDS (clusters object)")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(file.exists(opt$expr_rds))
stopifnot(file.exists(opt$meta_tsv))
stopifnot(nzchar(opt$direction))
stopifnot(nzchar(opt$kind))
stopifnot(nzchar(opt$view))
stopifnot(nzchar(opt$out_rds))

opt$kind <- tolower(opt$kind)
opt$view <- tolower(opt$view)
if (opt$view == "tumor") opt$view <- "tumour"
if (!(opt$kind %in% c("hc","kmeans"))) stop("[FATAL] --kind must be hc or kmeans")
if (!(opt$view %in% c("cell","tumour","cell_tumour"))) stop("[FATAL] --view must be cell|tumour|cell_tumour")

cat("[INFO] pan_agnostic_clustering\n")
cat("[INFO] expr_rds:", opt$expr_rds, "\n")
cat("[INFO] meta_tsv:", opt$meta_tsv, "\n")
cat("[INFO] direction:", opt$direction, " kind:", opt$kind, " view:", opt$view, "\n")
cat("[INFO] out_rds:", opt$out_rds, "\n")

# -----------------------------
# Load utils (same engine as per-disease)
# -----------------------------
# Locate pipeline base dir robustly (works from run/brca, etc.)
# -----------------------------
get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) return(NA_character_)
  sub("^--file=", "", file_arg[1])
}

script_path <- get_script_path()
if (is.na(script_path)) {
  stop("[FATAL] Cannot determine script path via --file=. Are you running with Rscript?")
}

# script is scripts/pan_agnostic_clustering.R relative to the repository root
base_dir <- normalizePath(file.path(dirname(script_path), ".."))

UTILS <- file.path(base_dir, "R", "agnostic_clustering_utils.R")
if (!file.exists(UTILS)) stop("[FATAL] Missing utils file: ", UTILS)
source(UTILS)
cat("[INFO] base_dir:", base_dir, "\n")
cat("[INFO] utils:", UTILS, "\n")

# -----------------------------
# Helpers
# -----------------------------
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) return(NULL)
  hit
}

infer_is_tumour <- function(meta) {
  # preferred boolean/0-1 flag
  flag_col <- pick_col(meta, c("is_tumour","is_tumor","tumour","tumor"))
  if (!is.null(flag_col)) {
    x <- meta[[flag_col]]
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    x2 <- tolower(trimws(as.character(x)))
    return(x2 %in% c("1","true","t","yes","y","tumour","tumor","tcga","primary"))
  }

  # string-based fallbacks
  type_col <- pick_col(meta, c("sample_type","type","dataset_type","source","cohort","dataset"))
  if (!is.null(type_col)) {
    x <- tolower(trimws(as.character(meta[[type_col]])))
    # tumour-like labels
    return(x %in% c("tumour","tumor","tcga","primary","primary_tumour","primary_tumor","tumour_sample","tumor_sample"))
  }

  NULL
}

# -----------------------------
# Load meta + expr
# -----------------------------
meta <- data.table::fread(opt$meta_tsv, data.table = FALSE)
expr <- readRDS(opt$expr_rds)
if (!is.matrix(expr)) expr <- as.matrix(expr)

# Identify sample id column
id_col <- pick_col(meta, c("sample_name","sample_id","sample","id","barcode","cell_line","cellline","Sample","CellLine"))
if (is.null(id_col)) {
  stop("[FATAL] Could not find sample ID column in metadata. Tried: sample_name, sample_id, sample, id, barcode, cell_line...")
}
meta_ids <- as.character(meta[[id_col]])

rn <- rownames(expr); cn <- colnames(expr)

samples_in_cols <- !is.null(cn) && any(cn %in% meta_ids)
samples_in_rows <- !is.null(rn) && any(rn %in% meta_ids)

if (!samples_in_cols && !samples_in_rows) {
  stop("[FATAL] Cannot match metadata IDs to expr rownames/colnames.\n",
       "meta head: ", paste(head(meta_ids, 3), collapse=", "), "\n",
       "expr rownames head: ", paste(head(rn, 3), collapse=", "), "\n",
       "expr colnames head: ", paste(head(cn, 3), collapse=", "))
}

# Standardize to genes x samples (what utils expects to load)
# If samples are rows, transpose.
if (samples_in_rows && !samples_in_cols) {
  expr <- t(expr)
  cn <- colnames(expr)
}

# ---- Ensembl version hygiene (match subset_expr.R behaviour) ----
strip_ens_version <- function(x) sub("\\.\\d+$", "", x)

# expr is genes x samples here
rownames(expr) <- strip_ens_version(rownames(expr))

# After stripping, duplicates can appear (ENSG...1.1 and ENSG...1.2 -> same)
# Keep the first occurrence deterministically
dup <- duplicated(rownames(expr))
if (any(dup)) {
  cat("[INFO] Dropping duplicated genes after version stripping:", sum(dup), "\n")
  expr <- expr[!dup, , drop = FALSE]
}

# Keep matched samples only
keep_ids <- intersect(colnames(expr), meta_ids)
if (length(keep_ids) < 3) stop("[FATAL] Too few matched samples after intersect(meta, expr): ", length(keep_ids))

expr <- expr[, keep_ids, drop=FALSE]
meta <- meta[match(keep_ids, meta[[id_col]]), , drop=FALSE]

# Split cell/tumour
is_tumour <- infer_is_tumour(meta)
if (is.null(is_tumour)) {
  stop("[FATAL] Cannot infer tumour vs cell from metadata. Add one of:\n",
       " - is_tumour (0/1 or TRUE/FALSE)\n",
       " - sample_type/type/dataset_type/dataset columns with values like 'CELL'/'TUMOUR' or 'TCGA'.")
}

cell_ids   <- keep_ids[!is_tumour]
tumour_ids <- keep_ids[ is_tumour]

if (length(cell_ids) < 3)   stop("[FATAL] Too few cell samples after split: ", length(cell_ids))
if (length(tumour_ids) < 3) stop("[FATAL] Too few tumour samples after split: ", length(tumour_ids))

cat("[INFO] matched samples:", length(keep_ids), "\n")
cat("[INFO] cell:", length(cell_ids), " tumour:", length(tumour_ids), "\n")

cell_gx   <- expr[, cell_ids, drop=FALSE]    # genes x samples
tumour_gx <- expr[, tumour_ids, drop=FALSE]  # genes x samples

# -----------------------------
# Map direction -> dist_method used in utils
# utils supports: euclidean, correlation (1 - Pearson)
# 
# NOTE: Pan-cancer directions with "*_cosine" are mapped to "correlation" 
# (1 - Pearson correlation), NOT true cosine distance. This is intentional:
# - The utils only support euclidean and correlation
# - Correlation is a reasonable proxy for cosine on normalized data
# - This matches the per-disease "*_corr" directions which also use correlation
# If true cosine distance is needed, extend utils to support it explicitly.
# -----------------------------
dir_lower <- tolower(opt$direction)
dist_method <- if (grepl("_cosine$", dir_lower) || grepl("cosine$", dir_lower)) {
  "correlation"
} else {
  "euclidean"
}
cat("[INFO] dist_method:", dist_method, " (note: *_cosine -> correlation, not true cosine)\n")

# -----------------------------
# Build utils "kind" string expected by run_agnostic_clustering
# -----------------------------
# raw methods: hc_cell, kmeans_tumour, kmeans_cell_tumour, etc.
kind_full <- paste0(opt$kind, "_", opt$view)   # e.g. kmeans_cell_tumour
cat("[INFO] kind_full:", kind_full, "\n")

# -----------------------------
# Write split matrices to temp RDS (because utils loads from files)
# -----------------------------
tmpdir <- tempdir()
cell_rds   <- file.path(tmpdir, sprintf("pan_cell_%s.rds", Sys.getpid()))
tumour_rds <- file.path(tmpdir, sprintf("pan_tumour_%s.rds", Sys.getpid()))
saveRDS(cell_gx, cell_rds)
saveRDS(tumour_gx, tumour_rds)

# -----------------------------
# Run clustering using shared engine
# -----------------------------
outdir <- dirname(opt$out_rds)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

run_agnostic_clustering(
  kind            = kind_full,
  cell_rds        = cell_rds,
  tumour_rds      = tumour_rds,
  outdir          = outdir,
  dist_method     = dist_method,
  cluster_rds_path = opt$out_rds
)

cat("[SUCCESS] wrote:", opt$out_rds, "\n")
