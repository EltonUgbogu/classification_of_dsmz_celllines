#!/usr/bin/env Rscript
# Preflight check: verify all required R packages are available
# This script fails fast if any packages are missing

required_pkgs <- c(
  # CLI and core utilities
  "optparse",
  "dplyr",
  "readr",
  "tidyr",
  "tibble",
  "stringr",
  "magrittr",
  
  # Plotting
  "ggplot2",
  "ggrepel",
  "viridis",
  "RColorBrewer",
  "scales",
  "pheatmap",
  "ggraph",
  
  # Matrix and math
  "Matrix",
  "matrixStats",
  
  # Clustering and dimensionality reduction
  "cluster",
  "factoextra",
  "dbscan",
  "umap",
  "uwot",
  
  # Graph analysis
  "igraph",
  
  # Bioconductor
  "Biobase",
  "BiocGenerics",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "limma",
  "ComplexHeatmap",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

cat("Checking required R packages...\n")
cat("Total packages to check:", length(required_pkgs), "\n\n")

missing <- character(0)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
    cat("  ✗ Missing:", pkg, "\n")
  } else {
    cat("  ✓ Found:", pkg, "\n")
  }
}

cat("\n")

# Optional: Check R version
cat("R version:", as.character(getRversion()), "\n")
cat("R library paths:\n")
for (path in .libPaths()) {
  cat("  -", path, "\n")
}
cat("\n")

if (length(missing) > 0) {
  cat("ERROR: Missing R packages:\n")
  for (pkg in missing) {
    cat("  -", pkg, "\n")
  }
  cat("\n")
  cat("Please install missing packages or update the conda environment.\n")
  quit(status = 1)
} else {
  cat("SUCCESS: All required R packages are available.\n")
  quit(status = 0)
}
