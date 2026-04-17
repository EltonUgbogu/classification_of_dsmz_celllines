#!/usr/bin/env Rscript
# =============================================================================
# preflight_checks.R
# =============================================================================
#
# Pre-flight checks for pan-cancer validation pipeline
# Verifies gene overlap and metadata structure
#
# USAGE:
# Rscript scripts/preflight_checks.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("========================================\n")
cat("PAN-CANCER VALIDATION PRE-FLIGHT CHECKS\n")
cat("========================================\n\n")

# Check 1: Gene list
cat("[1] Checking pan-cancer feature gene list...\n")
gene_file <- "results/unsupervised/pan_cancer/pan_cancer_features_clean.txt"
if (!file.exists(gene_file)) {
  stop("Gene file not found: ", gene_file)
}
genes <- readLines(gene_file)
genes <- genes[genes != ""]  # Remove empty lines
cat("  ✓ Loaded", length(genes), "genes\n")
if (length(genes) != 257) {
  cat("  ⚠ WARNING: Expected 257 genes, found", length(genes), "\n")
}

# Check 2: VST file overlap
cat("\n[2] Checking gene overlap with VST matrices...\n")
vst_files <- list(
  BRCA = "results/brca/inputs/vst_joint.rds",
  NBL = "results/nbl/inputs/vst_joint.rds",
  RBL = "results/rbl/inputs/vst_joint.rds"
)

for (disease in names(vst_files)) {
  path <- vst_files[[disease]]
  if (!file.exists(path)) {
    cat("  ⚠ WARNING: Missing", disease, "VST:", path, "\n")
    next
  }
  
  expr <- readRDS(path)
  if (is.matrix(expr)) {
    vst_genes <- rownames(expr)
  } else if (is.list(expr) && "expr" %in% names(expr)) {
    vst_genes <- rownames(expr$expr)
  } else {
    cat("  ⚠ WARNING: Unexpected structure for", disease, "\n")
    next
  }
  
  overlap <- sum(genes %in% vst_genes)
  overlap_pct <- round(overlap / length(genes) * 100, 1)
  cat("  ", disease, ": ", overlap, "/", length(genes), 
      " genes (", overlap_pct, "%)\n", sep="")
  
  if (overlap < 200) {
    cat("    ⚠ WARNING: Low overlap! Check gene ID format.\n")
  }
}

# Check 3: Metadata structure
cat("\n[3] Checking metadata structure...\n")
meta_file <- "data/dsmz/DSMZ_metadata.csv"
if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}

meta <- fread(meta_file, nrows=5)  # Just check structure
cat("  ✓ Metadata columns:", paste(colnames(meta), collapse=", "), "\n")

# Check for required columns
required <- c("sample_name", "Cell_Line", "DSMZ_Cell_line_norm", "Disease", "group2")
found <- intersect(required, colnames(meta))
if (length(found) > 0) {
  cat("  ✓ Found key columns:", paste(found, collapse=", "), "\n")
}

# Check 4: HEME samples in metadata
cat("\n[4] Checking for hematologic samples in metadata...\n")
meta_full <- fread(meta_file)
heme_cols <- c("Disease", "group2", "Description")
heme_found <- FALSE

for (col in heme_cols) {
  if (col %in% colnames(meta_full)) {
    heme_count <- sum(grepl("hematologic|leukemia|lymphoma|LL-100|HEME", 
                            meta_full[[col]], ignore.case=TRUE))
    if (heme_count > 0) {
      cat("  ✓ Found", heme_count, "potential HEME samples in column", col, "\n")
      heme_found <- TRUE
    }
  }
}

if (!heme_found) {
  cat("  ⚠ WARNING: No hematologic samples detected in metadata\n")
}

cat("\n========================================\n")
cat("PRE-FLIGHT CHECKS COMPLETE\n")
cat("========================================\n")
