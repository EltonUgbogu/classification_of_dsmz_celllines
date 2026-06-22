#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

# ==============================================================================
# Merge STAR count files into one combined genes x samples count matrix for RBL
#
# Sources:
#   {ROOT_DIR}/{cohort}/results/star/valid_count/{sample}.tab
#   4-column STAR ReadsPerGene format; sample ID = SRR accession
#
# Genes are matched on stripped Ensembl ID (no version suffix).
# The final matrix contains the shared gene order from the first sample and all
# subsequent samples are aligned to that same gene set.
# ==============================================================================

# ------------------------------
# Config
# ------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = FALSE)
ROOT_DIR <- Sys.getenv("ROOT_DIR", file.path(repo_root, "data", "rbl"))

default_cohorts <- c(
  "GSE111168_primary_tumours",
  "GSE196420_primary_tumours",
  "GSE268136_primary_tumours"
)
cohort_env <- Sys.getenv("PROJECTS", "")
COHORTS <- if (nzchar(cohort_env)) {
  strsplit(cohort_env, ",", fixed = TRUE)[[1]]
} else {
  default_cohorts
}
COHORTS <- trimws(COHORTS)
COHORTS <- COHORTS[nzchar(COHORTS)]

OUT_DIR    <- file.path(ROOT_DIR, "count_data")
OUT_MATRIX <- file.path(OUT_DIR, "rbl_tumour_count.rds")
OUT_META   <- file.path(OUT_DIR, "rbl_tumour_sample_metadata.csv")
OUT_GENES  <- file.path(OUT_DIR, "rbl_tumour_gene_list_hgnc.txt")

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
cat("Output dir :", OUT_DIR, "\n")
cat("Cohorts    :", paste(COHORTS, collapse = ", "), "\n\n")

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

count_mat[, 1]               <- first$counts
sample_table$stranded_col    <- NA_integer_
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

    count_mat[, i]               <- res$counts[idx]
    sample_table$stranded_col[i] <- res$col_idx
  }
}

cat("\nRBL matrix dimensions:\n")
cat("  genes   =", nrow(count_mat), "\n")
cat("  samples =", ncol(count_mat), "\n")

meta <- sample_table[, .(
  cohort,
  sample,
  stranded_col = as.character(stranded_col)
)]

# ------------------------------
# Save outputs
# ------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(count_mat, OUT_MATRIX)
fwrite(meta, OUT_META)
writeLines(rownames(count_mat), OUT_GENES)

cat("\nSaved:\n")
cat("  Matrix   :", OUT_MATRIX, "\n")
cat("  Metadata :", OUT_META, "\n")
cat("  Gene IDs :", OUT_GENES, "\n")
cat("\nDone.\n")
