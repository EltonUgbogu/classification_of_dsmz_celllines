#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Merge TARGET NBL GDC STAR augmented count TSVs into a genes x samples matrix
#
# Input layout:
#   {BASE}/count_data/{file_uuid}/*.rna_seq.augmented_star_gene_counts.tsv
#
# Sample identity is resolved via:
#   {BASE}/metadata/target_nbl_file_to_aliquot.tsv
#   (maps file_uuid -> aliquot_submitter_id, e.g. TARGET-30-PARACS-01A-01R)
#
# Output:
#   {BASE}/merged/target_nbl_count.rds
#   {BASE}/merged/target_nbl_sample_metadata.csv
#   {BASE}/merged/target_nbl_gene_list.txt
# ==============================================================================

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
BASE     <- Sys.getenv("TARGET_NBL_DIR",
                        "/work/ugbogu/pipeline/data/nbl/target_nbl")
COUNT_DIR <- file.path(BASE, "count_data")
META_FILE <- file.path(BASE, "metadata", "target_nbl_file_to_aliquot.tsv")
OUT_DIR   <- file.path(BASE, "merged")

OUT_MATRIX <- file.path(OUT_DIR, "target_nbl_count.rds")
OUT_META   <- file.path(OUT_DIR, "target_nbl_sample_metadata.csv")
OUT_GENES  <- file.path(OUT_DIR, "target_nbl_gene_list.txt")

# ------------------------------------------------------------------------------
# Load UUID -> aliquot mapping
# ------------------------------------------------------------------------------
if (!file.exists(META_FILE)) {
  stop("Missing mapping file: ", META_FILE)
}
mapping <- fread(META_FILE)
# Expected columns: file_id, file_name, case_submitter_id,
#                   sample_submitter_id, aliquot_submitter_id
required_cols <- c("file_id", "aliquot_submitter_id")
missing_cols  <- setdiff(required_cols, names(mapping))
if (length(missing_cols) > 0) {
  stop("Missing columns in mapping file: ", paste(missing_cols, collapse = ", "))
}
uuid_to_aliquot <- setNames(mapping$aliquot_submitter_id, mapping$file_id)

# ------------------------------------------------------------------------------
# Detect strandedness from augmented STAR count columns
# unstranded    = col 4 (after gene_id, gene_name, gene_type)
# stranded_first  = col 5
# stranded_second = col 6
# ------------------------------------------------------------------------------
detect_stranded_col <- function(dt) {
  data_rows <- dt[!grepl("^N_", gene_id)]

  s_un    <- sum(data_rows$unstranded,      na.rm = TRUE)
  s_first <- sum(data_rows$stranded_first,  na.rm = TRUE)
  s_sec   <- sum(data_rows$stranded_second, na.rm = TRUE)

  if (s_first >= 1.5 * s_sec  && s_first >= 0.6 * s_un) return("stranded_first")
  if (s_sec   >= 1.5 * s_first && s_sec   >= 0.6 * s_un) return("stranded_second")
  return("unstranded")
}

# ------------------------------------------------------------------------------
# Read one augmented STAR count TSV
# Returns: genes (Ensembl ID without version), counts, chosen column name
# ------------------------------------------------------------------------------
read_target_counts <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE,
              colClasses = "character")

  required <- c("gene_id", "unstranded", "stranded_first", "stranded_second")
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  # Convert count columns to numeric
  for (col in c("unstranded", "stranded_first", "stranded_second")) {
    dt[[col]] <- suppressWarnings(as.integer(dt[[col]]))
  }

  count_col <- detect_stranded_col(dt)

  # Drop STAR summary rows and strip gene version suffix
  dt <- dt[!grepl("^N_", gene_id)]
  genes  <- sub("\\.\\d+$", "", dt$gene_id)
  counts <- dt[[count_col]]
  counts[is.na(counts)] <- 0L

  list(genes = genes, counts = counts, count_col = count_col)
}

# ------------------------------------------------------------------------------
# Discover all UUID sample directories
# ------------------------------------------------------------------------------
uuid_dirs <- list.dirs(COUNT_DIR, full.names = TRUE, recursive = FALSE)
uuid_dirs <- uuid_dirs[dir.exists(uuid_dirs)]

if (length(uuid_dirs) == 0) {
  stop("No sample directories found under: ", COUNT_DIR)
}

uuids <- basename(uuid_dirs)

# Map UUIDs to aliquot IDs; warn and skip any without a mapping
aliquot_ids <- uuid_to_aliquot[uuids]
no_mapping  <- is.na(aliquot_ids)
if (any(no_mapping)) {
  message("[WARN] ", sum(no_mapping), " UUID(s) have no mapping and will be skipped:")
  message(paste(" ", uuids[no_mapping], collapse = "\n"))
  uuid_dirs   <- uuid_dirs[!no_mapping]
  uuids       <- uuids[!no_mapping]
  aliquot_ids <- aliquot_ids[!no_mapping]
}

# Find the TSV file within each UUID directory
find_tsv <- function(uuid_dir) {
  f <- list.files(uuid_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(f) == 0) return(NA_character_)
  if (length(f) > 1) {
    message("[WARN] Multiple TSV files in ", uuid_dir, "; using first.")
  }
  f[1]
}

tsv_files <- vapply(uuid_dirs, find_tsv, character(1))
missing_tsv <- is.na(tsv_files)
if (any(missing_tsv)) {
  message("[WARN] ", sum(missing_tsv), " director(ies) have no TSV file and will be skipped:")
  message(paste(" ", uuid_dirs[missing_tsv], collapse = "\n"))
  uuid_dirs   <- uuid_dirs[!missing_tsv]
  uuids       <- uuids[!missing_tsv]
  aliquot_ids <- aliquot_ids[!missing_tsv]
  tsv_files   <- tsv_files[!missing_tsv]
}

n_samples <- length(uuids)
if (n_samples == 0) stop("No valid samples remain after filtering.")

sample_table <- data.table(
  file_id    = uuids,
  aliquot_id = aliquot_ids,
  file       = tsv_files
)

cat("==> TARGET NBL count merge\n")
cat("    Input dir  :", COUNT_DIR, "\n")
cat("    Output dir :", OUT_DIR, "\n")
cat("    Samples    :", n_samples, "\n\n")

# ------------------------------------------------------------------------------
# Build count matrix
# ------------------------------------------------------------------------------
cat("Reading first sample as template:", sample_table$aliquot_id[1], "\n")
first <- read_target_counts(sample_table$file[1])

gene_ids  <- first$genes
n_genes   <- length(gene_ids)

count_mat <- matrix(
  0L,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_ids, sample_table$aliquot_id)
)

count_mat[, 1]             <- first$counts
sample_table$count_col     <- NA_character_
sample_table$count_col[1]  <- first$count_col

if (n_samples > 1) {
  for (i in 2:n_samples) {
    cat("Reading sample", i, "of", n_samples, ":", sample_table$aliquot_id[i], "\n")
    res <- read_target_counts(sample_table$file[i])

    idx <- match(gene_ids, res$genes)
    if (anyNA(idx)) {
      stop("Gene mismatch in sample ", sample_table$aliquot_id[i],
           ": ", sum(is.na(idx)), " template genes not found.")
    }

    count_mat[, i]            <- res$counts[idx]
    sample_table$count_col[i] <- res$count_col
  }
}

cat("\nFinal matrix dimensions:\n")
cat("  genes   =", nrow(count_mat), "\n")
cat("  samples =", ncol(count_mat), "\n")

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(count_mat, OUT_MATRIX)
fwrite(sample_table, OUT_META)
writeLines(rownames(count_mat), OUT_GENES)

cat("\nSaved:\n")
cat("  Matrix   :", OUT_MATRIX, "\n")
cat("  Metadata :", OUT_META, "\n")
cat("  Gene IDs :", OUT_GENES, "\n")
cat("\nDone.\n")
