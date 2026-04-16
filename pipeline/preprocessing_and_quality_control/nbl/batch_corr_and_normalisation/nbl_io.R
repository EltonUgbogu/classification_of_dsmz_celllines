#!/usr/bin/env Rscript


#============================================================
# nbl_io.R: Load and preprocess NBL data
#==========================================================
#
# PURPOSE:
# Load NBL Tumour purity corrected data 
# and filter DSMZ raw count matrix by cancer type
#
#
# USAGE:
# Rscript nbl_io.R --config config/config.yaml --tumour_counts /path/to/tumour_counts.rds --dsmz_counts /path/to/dsmz_counts.rds --dsmz_meta /path/to/dsmz_meta.csv --dsmz_sample_filter "Breast"
#
# ARGUMENTS:
# --config: Path to YAML config file
# --tumour_counts: Path to Tumour count matrix
# --dsmz_counts: Path to DSMZ raw count matrix
# --dsmz_meta: Path to DSMZ metadata file
# --dsmz_sample_filter: Filter DSMZ metadata and count matrix by cancer type

#============================================================

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to YAML config file [default: %default]"),
  make_option("--tumour_counts", type = "character",
             help = "Tumour count matrix for neuroblastoma"),
  make_option("--dsmz_counts", type = "character", default = NULL,
             help = "DSMZ raw count matrix "),
  make_option("--dsmz_meta", type = "character", default = NULL,
             help = "DSMZ metadata file"),
  make_option("--dsmz_sample_filter", type = "character", default = NULL,
            help = "Filter DSMZ metadata and count matrix by cancer type"),
)


# Load DSMZ raw wide table (counts + gene columns) and metadata CSV
#counts_path - path to DSMZ count data RDS found in config
#meta_path - path to DSMZ metadata CSV found in config
load_dsmz_data <- function(counts_path, meta_path) {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) stop(sprintf("[ERROR] DSMZ RDS not found: %s", counts_path))
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  raw  <- readRDS(counts_path)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  # Validate raw data columns
  stopifnot(all(c("Ensembl_ID", "gene_name") %in% colnames(raw)))
  # Validate metadata has sample_name
  stopifnot("sample_name" %in% colnames(meta))
  cat(sprintf("[INFO] DSMZ raw counts: %d rows x %d cols\n", nrow(raw), ncol(raw))) # expect 167 samples + 3 column annotatation metedata  = 170 colnames
  return(list(raw = raw, meta = meta))
}

#-----------NBL Tumour loading function-----------
load_tumour_nbl_data <- function(file_path) {
  cat("[INFO] Loading tumour data...\n")
  if (!file.exists(file_path)) stop(sprintf("[ERROR] NBL Tumour RDS not found: %s", file_path))
  nbl_tumour_se <- readRDS(file_path)
  tumour_counts <- assay(tumour_se)
  
  if (any(duplicated(rownames(tumour_counts)))) {
    dup_genes <- rownames(tumour_counts)[duplicated(rownames(tumour_counts))]
    cat(sprintf("[INFO] Collapsing %d duplicate tumour genes: %s\n",
                length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
    tumour_counts <- rowsum(tumour_counts, group = rownames(tumour_counts), reorder = TRUE)
  }
  
  storage.mode(tumour_counts) <- "double"
  cat(sprintf("[DIMS] tumour dims: %d genes x %d samples\n", nrow(tumour_counts), ncol(tumour_counts)))  
  
  cat("[DIMS] tumour_nbl_counts dimensions: ", dim(tumour_nbl_counts), "\n")
  return(tumour_nbl_counts)
}

#-----------DSMZ loading function-----------
# count - dsmz count data after buld_dsmz_matrix function
load_dsmz_nbl_data <- function(counts, meta_path, dsmz_disease_filter) {
  cat("[INFO] Loading DSMZ neuroblastoma data...\n")
  # Validate inputs
  if (is.null(counts) || !is.matrix(counts) && !is.data.frame(counts)) {
    stop("[ERROR] 'counts' must be a non-null matrix or data frame")
  }
  if (!file.exists(meta_path)) {
    stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  }
  # Load metadata
  dsmz_meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  # Validate required columns
  stopifnot("sample_name" %in% colnames(dsmz_meta))
  stopifnot("organ" %in% colnames(dsmz_meta))
  # Filter by organ if specified
  if (!is.null(dsmz_organ_filter) && length(dsmz_organ_filter) > 0) {
    keep_ids <- dsmz_meta$sample_name[dsmz_meta$organ %in% dsmz_organ_filter]
  } else {
    keep_ids <- dsmz_meta$sample_name  # Keep all if no filter
  }
  # Validate sample overlap
  missing <- setdiff(keep_ids, colnames(counts))
  if (length(missing) > 0) {
    stop(sprintf("[ERROR] %d samples in metadata not found in counts: %s", 
                 length(missing), paste(head(missing), collapse = ", ")))
  }
  # Subset counts
  dsmz_nbl_counts <- counts[, keep_ids, drop = FALSE]
  
  cat(sprintf("[INFO] Selected %d samples for organ(s): %s\n", 
              ncol(dsmz_nbl_counts), paste(dsmz_organ_filter, collapse = ", ")))
  
  return(dsmz_nbl_counts)
}

