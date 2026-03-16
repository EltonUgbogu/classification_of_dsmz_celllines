#!/usr/bin/env Rscript

#build_tumour_neighbourhood_input.R

# Build joint DSMZ+TCGA expression matrix + mapping for tumour neighbourhoods
# Parameterized by gene list so it can support HVG and other feature-set directions.

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
              help = "Path to DSMZ metadata CSV"),
  make_option("--expr-out",  dest = "expr.out",  type = "character",
              help = "Output RDS for samples x genes expression"),
  make_option("--map-out",   dest = "map.out",   type = "character",
              help = "Output RDS mapping DSMZ samples to cell lines"),
  make_option("--direction", dest = "direction", type = "character",
              default = NA, help = "Direction label (for logging only)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Required arguments are keyed by their dest names defined above.
required <- c("vst", "gene.list", "dsmz.meta", "expr.out", "map.out")
missing_args <- required[sapply(required, function(k) {
  val <- opt[[k]]
  is.null(val) || !nzchar(val)
})]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

cfg <- yaml::read_yaml(opt$config)

cat("=== Building tumour neighbourhood inputs ===\n")
cat("[INFO] Direction      : ", ifelse(is.na(opt$direction), "(not set)", opt$direction), "\n", sep = "")
cat("[INFO] VST matrix     : ", opt$vst, "\n", sep = "")
cat("[INFO] Gene list      : ", opt$gene.list, "\n", sep = "")
cat("[INFO] DSMZ metadata  : ", opt$dsmz.meta, "\n", sep = "")
cat("[INFO] Expr output    : ", opt$expr.out, "\n", sep = "")
cat("[INFO] Mapping output : ", opt$map.out, "\n", sep = "")

stopifnot(file.exists(opt$vst))
stopifnot(file.exists(opt$gene.list))
stopifnot(file.exists(opt$dsmz.meta))

vst_joint <- readRDS(opt$vst)
cat("[INFO] VST dims:", nrow(vst_joint), "genes x", ncol(vst_joint), "samples\n")

# DSMZ/TCGA sanity check (informative only)
n_dsmz <- sum(grepl("^NG-", colnames(vst_joint)))
n_tcga <- sum(grepl("^TCGA-", colnames(vst_joint)))
cat("[INFO] Detected DSMZ samples (NG- prefix):", n_dsmz, "\n")
cat("[INFO] Detected TCGA samples (TCGA- prefix):", n_tcga, "\n")

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

# Build samples × genes matrix
expr_mat <- t(vst_joint[genes, , drop = FALSE])

# Extract DSMZ sample IDs (columns in original VST, now rownames in expr_mat)
dsmz_sample_ids <- rownames(expr_mat)[grepl("^NG-", rownames(expr_mat))]

# Read DSMZ metadata and build mapping (same logic as PAM50 script)
# Metadata files are TSV (even if named *.csv), so use read.delim.
meta <- read.delim(opt$dsmz.meta, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("sample_name" %in% colnames(meta))
stopifnot("sample_id" %in% colnames(meta))

# Value for mapping = canonical cell line name (key = tech sample ID).
# Prefer explicit cell line column if present.
cell_col <- intersect(
  c("Cell_line", "Cell_Line", "CellLine", "cell_line", "cellLine"),
  colnames(meta)
)[1]

if (!is.na(cell_col)) {
  cell_line_of <- setNames(meta[[cell_col]], meta$sample_name)
} else {
  # No cell_line column: derive from DSMZ tech ID (e.g. NG-13640_BC3_lib... -> BC3).
  # Do NOT use cohort/lineage (e.g. "HEME") — that would map every sample to the same label.
  meta$cell_line_derived <- meta$sample_id
  dsmz_rows <- grepl("^NG-", meta$sample_id)
  if (any(dsmz_rows)) {
    meta$cell_line_derived[dsmz_rows] <- sub("^NG-[0-9]+_(.+)_lib.*", "\\1", meta$sample_id[dsmz_rows])
  }
  cell_line_of <- setNames(meta$cell_line_derived, meta$sample_name)
}

# Map DSMZ sample IDs to cell lines
cell_line_to_sample_id <- cell_line_of[dsmz_sample_ids]
# For samples not found in metadata, use the sample ID itself
cell_line_to_sample_id[is.na(cell_line_to_sample_id)] <- dsmz_sample_ids[is.na(cell_line_to_sample_id)]

# Create output directories
dir.create(dirname(opt$expr.out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$map.out), recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(expr_mat, opt$expr.out)
saveRDS(cell_line_to_sample_id, opt$map.out)

cat("\n[OK] Saved expression matrix:\n  ", opt$expr.out, "\n", sep = "")
cat("[OK] Saved mapping:\n  ", opt$map.out, "\n", sep = "")
cat("[INFO] Final dims: ", nrow(expr_mat), " samples x ", ncol(expr_mat), " genes\n", sep = "")
cat("[INFO] DSMZ: ", sum(grepl("^NG-", rownames(expr_mat))),
    " | TCGA: ", sum(grepl("^TCGA-", rownames(expr_mat))), "\n", sep = "")
cat("[INFO] Mapping entries:", length(cell_line_to_sample_id), " DSMZ samples\n")
