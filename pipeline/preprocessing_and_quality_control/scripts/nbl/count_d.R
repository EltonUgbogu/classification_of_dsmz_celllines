#!/usr/bin/env Rscript
# ==============================================================================
# fastq_2_countdata.R
# Merge STAR counts → KEEP GENE NAMES AS ENSEMBL IDs
# ==============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

# ==============================
# Config
# ==============================
ROOT_DIR <- "/home/chu25/data/nbl"
COHORTS  <- c("GSE100148", "GSE189367", "GSE49711", "SRP409177", "TARGET-NBL")

OUT_MATRIX <- file.path(ROOT_DIR, "nbl_tumour_count.rds")
OUT_META   <- file.path(ROOT_DIR, "nbl_tumour_sample_metadata.csv")
OUT_GENES  <- file.path(ROOT_DIR, "nbl_tumour_gene_list_hgnc.txt")

# ==============================
# Detect strandedness (same as Python)
# ==============================
detect_stranded <- function(dt) {
  dt <- dt[!grepl("^N_", V1)]
  if (ncol(dt) < 4) return(2)
  s2 <- sum(suppressWarnings(as.numeric(dt[[2]])), na.rm = TRUE)
  s3 <- sum(suppressWarnings(as.numeric(dt[[3]])), na.rm = TRUE)
  s4 <- sum(suppressWarnings(as.numeric(dt[[4]])), na.rm = TRUE)
  if (s3 >= 1.5 * s4 && s3 >= 0.6 * s2) return(3)
  if (s4 >= 1.5 * s3 && s4 >= 0.6 * s2) return(4)
  return(2)
}

# ==============================
# Load one file
# ==============================
load_one <- function(path, sample_name) {
  cat("  →", sample_name, "\n")
  dt <- fread(path, header = FALSE, sep = "\t", colClasses = "character")
  col_idx <- detect_stranded(dt)
  dt <- dt[!grepl("^N_", V1), .(gene = V1, count = as.integer(V(col_idx)))]
  dt[, gene := sub("\\.\\d+$", "", gene)]  # strip version
  vec <- dt$count
  names(vec) <- dt$gene
  return(vec)
}

# ==============================
# Main loop
# ==============================
cat("Starting merge (Ensembl IDs only)...\n")
all_counts <- list()
meta <- list()

for (cohort in COHORTS) {
  dir <- file.path(ROOT_DIR, cohort, "results", "valid_count")
  files <- list.files(dir, pattern = "\\.tab$", full.names = TRUE)
  if (length(files) == 0) { cat("  No files in", cohort, "\n"); next }
  
  cat("\n[", cohort, "] →", length(files), "samples\n")
  for (f in files) {
    sname <- paste0(cohort, "_", tools::file_path_sans_ext(basename(f)))
    tryCatch({
      vec <- load_one(f, sname)
      all_counts[[sname]] <- vec
      meta[[sname]] <- data.frame(sample = sname, cohort = cohort, file = f)
    }, error = function(e) cat("  ERROR:", e$message, "\n"))
  }
}

if (length(all_counts) == 0) stop("No samples loaded.")

cat("\nMerging matrix...\n")
# Find all genes
all_genes <- unique(unlist(lapply(all_counts, names)))
mat <- matrix(0L, nrow = length(all_counts), ncol = length(all_genes),
              dimnames = list(names(all_counts), all_genes))

for (s in names(all_counts)) {
  mat[s, names(all_counts[[s]])] <- all_counts[[s]]
}

final_mat <- t(mat)  # genes x samples → samples x genes
cat("FINAL: ", nrow(final_mat), "samples ×", ncol(final_mat), "Ensembl genes\n")

# ==============================
# Save
# ==============================
saveRDS(final_mat, OUT_MATRIX)
cat("SAVED matrix →", OUT_MATRIX, "\n")

meta_df <- do.call(rbind, meta)
write.csv(meta_df, OUT_META, row.names = FALSE)
cat("SAVED metadata →", OUT_META, "\n")

writeLines(colnames(final_mat), OUT_GENES)
cat("SAVED gene list →", OUT_GENES, "\n")
cat("Total Ensembl genes:", ncol(final_mat), "\n")

cat("\nDONE! Load with:\n")
cat("  counts <- readRDS('", OUT_MATRIX, "')   # samples x Ensembl IDs\n")
cat("  meta   <- read.csv('", OUT_META, "')\n")