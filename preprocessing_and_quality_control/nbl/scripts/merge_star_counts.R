#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

# ==============================================================================
# Merge STAR count files into one combined genes x samples count matrix
#
# Sources:
#   1. GEO cohorts (configured by COHORTS; currently GSE189367, SRP409177)
#      {ROOT_DIR}/{cohort}/results/star/valid_count/{sample}.tab
#      4-column STAR ReadsPerGene format; sample ID = SRR accession
#
#   2. TARGET-NBL (GDC)
#      {ROOT_DIR}/target_nbl/merged/target_nbl_count.rds
#      Pre-merged by merge_target_nbl_counts.R; columns = aliquot_submitter_id
#
# Genes are matched on stripped Ensembl ID (no version suffix).
# The final matrix is the intersection of genes across both sources.
# ==============================================================================

# ------------------------------
# Config
# ------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = FALSE)
ROOT_DIR <- Sys.getenv("ROOT_DIR", file.path(repo_root, "data", "nbl"))
COHORTS_ENV <- Sys.getenv("COHORTS", "GSE189367,SRP409177")
COHORTS <- trimws(strsplit(COHORTS_ENV, ",", fixed = TRUE)[[1]])
COHORTS <- COHORTS[nzchar(COHORTS)]
if (length(COHORTS) == 0) {
  stop("No GEO cohorts configured. Set COHORTS to a comma-separated project list.")
}
TARGET_RDS <- file.path(ROOT_DIR, "target_nbl", "merged", "target_nbl_count.rds")
TARGET_META_FILE <- file.path(ROOT_DIR, "target_nbl", "merged", "target_nbl_sample_metadata.csv")

OUT_DIR    <- file.path(ROOT_DIR, "count_data")
OUT_MATRIX <- file.path(OUT_DIR, "nbl_tumour_count.rds")
OUT_META   <- file.path(OUT_DIR, "nbl_tumour_sample_metadata.csv")
OUT_GENES  <- file.path(OUT_DIR, "nbl_tumour_gene_list_hgnc.txt")

# ------------------------------
# Detect strandedness from STAR ReadsPerGene.out.tab format
# col 2 = unstranded
# col 3 = stranded forward
# col 4 = stranded reverse
# ------------------------------
detect_stranded <- function(dt) {
  dt2 <- dt[!grepl("^N_", V1)]
  if (ncol(dt2) < 4) return(2L)

  s2 <- sum(suppressWarnings(as.numeric(dt2[[2]])), na.rm = TRUE)
  s3 <- sum(suppressWarnings(as.numeric(dt2[[3]])), na.rm = TRUE)
  s4 <- sum(suppressWarnings(as.numeric(dt2[[4]])), na.rm = TRUE)

  if (s3 >= 1.5 * s4 && s3 >= 0.6 * s2) return(3L)
  if (s4 >= 1.5 * s3 && s4 >= 0.6 * s2) return(4L)
  return(2L)
}

# ------------------------------
# Read one STAR count file
# Input:  path to {sample}.tab (valid_count/ output from Snakefile)
# Returns:
#   genes   = Ensembl IDs without version suffix
#   counts  = integer vector
#   col_idx = chosen count column (2=unstranded, 3=fwd, 4=rev)
# ------------------------------
read_star_counts <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t", colClasses = "character")

  if (ncol(dt) < 4) {
    stop("Expected at least 4 columns in STAR count file: ", path)
  }

  col_idx <- detect_stranded(dt)

  dt <- dt[!grepl("^N_", V1)]
  genes  <- sub("\\.\\d+$", "", dt[[1]])
  counts <- suppressWarnings(as.integer(dt[[col_idx]]))
  counts[is.na(counts)] <- 0L

  list(genes = genes, counts = counts, col_idx = col_idx)
}

# ------------------------------
# Collect all sample files
# Scans {ROOT_DIR}/{cohort}/results/star/valid_count/*.tab
# ------------------------------
collect_sample_table <- function(root_dir, cohorts) {
  rows <- list()

  for (cohort in cohorts) {
    count_dir <- file.path(root_dir, cohort, "results", "star", "valid_count")

    if (!dir.exists(count_dir)) {
      message("[WARN] Missing directory: ", count_dir)
      next
    }

    files <- sort(list.files(count_dir, pattern = "\\.tab$", full.names = TRUE))

    if (length(files) == 0) {
      message("[WARN] No .tab files found in: ", count_dir)
      next
    }

    sample_ids <- tools::file_path_sans_ext(basename(files))

    rows[[cohort]] <- data.table(
      cohort = cohort,
      sample = sample_ids,
      file   = files
    )
  }

  if (length(rows) == 0) {
    stop("No STAR .tab files found in any cohort.")
  }

  sample_table <- rbindlist(rows, use.names = TRUE, fill = TRUE)

  dup_ids <- sample_table$sample[duplicated(sample_table$sample)]
  if (length(dup_ids) > 0) {
    stop(
      "Duplicate sample IDs detected across cohorts: ",
      paste(unique(dup_ids), collapse = ", "),
      "\nUse cohort-prefixed sample names if these are genuine duplicates."
    )
  }

  sample_table
}

# ------------------------------
# Main
# ------------------------------
cat("Collecting STAR count files...\n")
cat("Input root :", ROOT_DIR, "\n")
cat("Output dir :", OUT_DIR, "\n\n")

sample_table <- collect_sample_table(ROOT_DIR, COHORTS)

cat("Found", nrow(sample_table), "samples across",
    length(unique(sample_table$cohort)), "cohorts.\n\n")
print(sample_table[, .N, by = cohort])

# Read first sample to establish gene universe and order
cat("\nReading first sample as template:\n")
cat("  ", sample_table$sample[1], "\n")
first <- read_star_counts(sample_table$file[1])

gene_ids  <- first$genes
n_genes   <- length(gene_ids)
n_samples <- nrow(sample_table)

count_mat <- matrix(
  0L,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_ids, sample_table$sample)
)

count_mat[, 1]             <- first$counts
sample_table$stranded_col  <- NA_integer_
sample_table$stranded_col[1] <- first$col_idx

# Read remaining samples
if (n_samples > 1) {
  for (i in 2:n_samples) {
    sample_id <- sample_table$sample[i]
    file_path <- sample_table$file[i]

    cat("Reading sample", i, "of", n_samples, ":", sample_id, "\n")
    res <- read_star_counts(file_path)

    idx <- match(gene_ids, res$genes)

    if (anyNA(idx)) {
      stop(
        "Gene mismatch in sample ", sample_id,
        ": ", sum(is.na(idx)), " template genes not found."
      )
    }

    count_mat[, i]            <- res$counts[idx]
    sample_table$stranded_col[i] <- res$col_idx
  }
}

cat("\nGEO matrix dimensions:\n")
cat("  genes   =", nrow(count_mat), "\n")
cat("  samples =", ncol(count_mat), "\n")

# ------------------------------
# Load and merge TARGET-NBL matrix
# ------------------------------
if (!file.exists(TARGET_RDS)) {
  stop("TARGET-NBL matrix not found: ", TARGET_RDS,
       "\nRun merge_target_nbl_counts.R first.")
}

cat("\nLoading TARGET-NBL matrix:", TARGET_RDS, "\n")
target_mat <- readRDS(TARGET_RDS)

cat("TARGET-NBL matrix dimensions:\n")
cat("  genes   =", nrow(target_mat), "\n")
cat("  samples =", ncol(target_mat), "\n")

# Intersect genes — GEO uses GENCODE v44, TARGET uses GENCODE v36
common_genes <- intersect(rownames(count_mat), rownames(target_mat))
cat("\nCommon genes (intersection):", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop("No common genes between GEO and TARGET-NBL matrices.")
}

geo_sub    <- count_mat[common_genes, , drop = FALSE]
target_sub <- target_mat[common_genes, , drop = FALSE]

combined_mat <- cbind(geo_sub, target_sub)
cat("\nCombined matrix dimensions:\n")
cat("  genes   =", nrow(combined_mat), "\n")
cat("  samples =", ncol(combined_mat), "\n")

# Build combined metadata
target_meta <- fread(TARGET_META_FILE)
target_meta[, cohort := "TARGET-NBL"]

if (!"aliquot_id" %in% names(target_meta)) {
  stop("TARGET metadata is missing required column: aliquot_id")
}

setnames(target_meta, "count_col", "stranded_col", skip_absent = TRUE)

geo_meta <- sample_table[, .(cohort, sample, stranded_col = as.character(stranded_col))]
target_slim <- target_meta[, .(cohort, sample = aliquot_id,
                               file_id, file, stranded_col)]
target_slim[, stranded_col := as.character(stranded_col)]

combined_meta <- rbindlist(list(geo_meta, target_slim),
                           use.names = TRUE, fill = TRUE)

# ------------------------------
# Save outputs
# ------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(combined_mat, OUT_MATRIX)
fwrite(combined_meta, OUT_META)
writeLines(rownames(combined_mat), OUT_GENES)

cat("\nSaved:\n")
cat("  Matrix   :", OUT_MATRIX, "\n")
cat("  Metadata :", OUT_META, "\n")
cat("  Gene IDs :", OUT_GENES, "\n")
cat("\nDone.\n")
