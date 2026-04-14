#!/usr/bin/env Rscript

# build_tumour_neighbourhood_input.R

# Build joint DSMZ + tumour expression matrix + mapping for tumour neighbourhoods.
# Parameterized by gene list so it can support HVG and other feature-set directions.
# Tumours are defined as any sample that does not carry the ^NG- (DSMZ cell line) prefix,
# generalising across cohorts (TCGA, TARGET, ICGC, etc.).

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
})

option_list <- list(
  make_option("--config",    dest = "config",    type = "character",
              default = "config/config.yaml", help = "Path to config YAML"),
  make_option("--vst",       dest = "vst",       type = "character",
              help = "Path to joint VST matrix (genes x samples)"),
  make_option("--gene-list", dest = "gene.list", type = "character",
              help = "Path to gene list text file"),
  make_option("--dsmz-meta", dest = "dsmz.meta", type = "character",
              help = "Path to DSMZ metadata CSV/TSV"),
  make_option("--expr-out",  dest = "expr.out",  type = "character",
              help = "Output RDS for samples x genes expression"),
  make_option("--map-out",   dest = "map.out",   type = "character",
              help = "Output RDS mapping DSMZ samples to cell lines"),
  make_option("--direction", dest = "direction", type = "character",
              default = NA, help = "Direction label (for logging only)")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("vst", "gene.list", "dsmz.meta", "expr.out", "map.out")
missing_args <- required[sapply(required, function(k) {
  val <- opt[[k]]
  is.null(val) || !nzchar(val)
})]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

cfg <- yaml::read_yaml(opt$config)

# Resolve column names from config, with safe defaults.
# cfg$deseq2$sample_id_col  — technical sample ID column (matches VST matrix column names)
# cfg$deseq2$cell_line_col  — canonical cell line name column
sample_id_col <- cfg$deseq2$sample_id_col
if (is.null(sample_id_col) || !nzchar(sample_id_col)) sample_id_col <- "sample_id"
cell_col <- cfg$deseq2$cell_line_col
if (is.null(cell_col) || !nzchar(cell_col)) cell_col <- "cell_line"

cat("=== Building tumour neighbourhood inputs ===\n")
cat("[INFO] Direction      : ", ifelse(is.na(opt$direction), "(not set)", opt$direction), "\n", sep = "")
cat("[INFO] VST matrix     : ", opt$vst, "\n", sep = "")
cat("[INFO] Gene list      : ", opt$gene.list, "\n", sep = "")
cat("[INFO] DSMZ metadata  : ", opt$dsmz.meta, "\n", sep = "")
cat("[INFO] Expr output    : ", opt$expr.out, "\n", sep = "")
cat("[INFO] Mapping output : ", opt$map.out, "\n", sep = "")
cat("[INFO] sample_id_col  : ", sample_id_col, "\n", sep = "")
cat("[INFO] cell_line_col  : ", cell_col, "\n", sep = "")

stopifnot(file.exists(opt$vst))
stopifnot(file.exists(opt$gene.list))
stopifnot(file.exists(opt$dsmz.meta))

vst_joint <- readRDS(opt$vst)
cat("[INFO] VST dims:", nrow(vst_joint), "genes x", ncol(vst_joint), "samples\n")

# DSMZ cell line samples carry the ^NG- prefix.
# All other samples are tumours regardless of cohort (TCGA, TARGET, ICGC, etc.).
n_dsmz   <- sum(grepl("^NG-", colnames(vst_joint)))
n_tumour <- ncol(vst_joint) - n_dsmz
cat("[INFO] Detected DSMZ cell line samples (NG- prefix):", n_dsmz, "\n")
cat("[INFO] Detected tumour samples (non-NG-)            :", n_tumour, "\n")

if (n_dsmz == 0) {
  stop("No DSMZ samples found in VST matrix. Please build a joint matrix first (e.g., scripts/rebuild_joint_vst.R).")
}

# Read gene list
genes <- scan(opt$gene.list, what = character(), quiet = TRUE)
cat("[INFO] Gene list length:", length(genes), "\n")

genes <- intersect(genes, rownames(vst_joint))
cat("[INFO] Genes present in VST:", length(genes), "\n")
if (length(genes) < 50) {
  stop("Too few genes found in VST matrix: ", length(genes))
}

# Build samples x genes matrix
expr_mat <- t(vst_joint[genes, , drop = FALSE])

# Extract DSMZ sample IDs (columns in original VST, now rownames in expr_mat)
dsmz_sample_ids <- rownames(expr_mat)[grepl("^NG-", rownames(expr_mat))]

# Read DSMZ metadata.
# Files are TSV (even if named *.csv), so use read.delim.
meta <- read.delim(opt$dsmz.meta, stringsAsFactors = FALSE, check.names = FALSE)
if (!sample_id_col %in% colnames(meta)) {
  stop("sample_id_col '", sample_id_col, "' not found in metadata columns: ",
       paste(colnames(meta), collapse = ", "))
}

# Build mapping: technical sample ID -> canonical cell line name.
# Keys must be the technical IDs that appear as column names in the VST matrix
# (i.e., meta[[sample_id_col]]), so that lookups via dsmz_sample_ids resolve correctly.
if (cell_col %in% colnames(meta)) {
  cell_line_of <- setNames(meta[[cell_col]], meta[[sample_id_col]])
} else {
  # No cell_line column: derive the label from the technical sample ID.
  # e.g. NG-13640_BC3_lib... -> BC3
  # Do NOT use cohort/lineage (e.g. "HEME") — that would map every sample to the same label.
  meta$cell_line_derived <- meta[[sample_id_col]]
  dsmz_rows <- grepl("^NG-", meta[[sample_id_col]])
  if (any(dsmz_rows)) {
    meta$cell_line_derived[dsmz_rows] <- sub(
      "^NG-[0-9]+_(.+)_lib.*", "\\1", meta[[sample_id_col]][dsmz_rows]
    )
  }
  cell_line_of <- setNames(meta$cell_line_derived, meta[[sample_id_col]])
}

# Map DSMZ sample IDs to cell lines
cell_line_to_sample_id <- cell_line_of[dsmz_sample_ids]
# For samples not found in metadata, fall back to the sample ID itself
cell_line_to_sample_id[is.na(cell_line_to_sample_id)] <- dsmz_sample_ids[is.na(cell_line_to_sample_id)]

# Create output directories
dir.create(dirname(opt$expr.out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$map.out),  recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(expr_mat,              opt$expr.out)
saveRDS(cell_line_to_sample_id, opt$map.out)

cat("\n[OK] Saved expression matrix:\n  ", opt$expr.out, "\n", sep = "")
cat("[OK] Saved mapping:\n  ", opt$map.out, "\n", sep = "")
cat("[INFO] Final dims: ", nrow(expr_mat), " samples x ", ncol(expr_mat), " genes\n", sep = "")
cat("[INFO] DSMZ: ",    sum(grepl("^NG-", rownames(expr_mat))),
    " | tumour: ", sum(!grepl("^NG-", rownames(expr_mat))), "\n", sep = "")
cat("[INFO] Mapping entries:", length(cell_line_to_sample_id), " DSMZ samples\n")
