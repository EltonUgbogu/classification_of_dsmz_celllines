config <- list(
  tcga_se_rds = "/home/chu25/data/final_nbl_count_matrix_hgnc.parquet",
  dsmz_rds = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",
  tcga_project_filter = c("TCGA-BRCA"),
  dsmz_organ_filter = c("Breast"),
  tcga_nbl_subtype_metadata = "/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/TCGA_BRCA_project_subtype.csv",
  outdir = "/home/chu25/dsmz/old_scripts/nbl_analysis/nbl/unsupervised/results"
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

#-----------TCGA-BRCA loading function-----------
load_tcga_nbl_data <- function(file_path, tcga_project_filter) {
  cat("[INFO] Loading TCGA data...\n")
  if (!file.exists(file_path)) stop(sprintf("[ERROR] TCGA RDS not found: %s", file_path))
  tcga_se <- readRDS(file_path)
  tcga_counts <- assay(tcga_se)
  
  if (any(duplicated(rownames(tcga_counts)))) {
    dup_genes <- rownames(tcga_counts)[duplicated(rownames(tcga_counts))]
    cat(sprintf("[INFO] Collapsing %d duplicate TCGA genes: %s\n",
                length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
    tcga_counts <- rowsum(tcga_counts, group = rownames(tcga_counts), reorder = TRUE)
  }
  
  storage.mode(tcga_counts) <- "double"
  cat(sprintf("[DIMS] TCGA dims: %d genes x %d samples\n", nrow(tcga_counts), ncol(tcga_counts)))
  
  if (is.null(tcga_project_filter) || length(tcga_project_filter) == 0) {
    stop("tcga_project_filter must not be null or empty")
  }
  
  sample_names <- colnames(tcga_counts)
  project_ids <- as.character(colData(tcga_se)[sample_names, "project_id"])
  
  keep_tcga_nbl <- project_ids %in% tcga_project_filter
  if (sum(keep_tcga_nbl) == 0) {
    stop("No TCGA samples match the project filter: ", paste(tcga_project_filter, collapse = ", "))
  }
  
  tcga_nbl_counts <- tcga_counts[, keep_tcga_nbl, drop = FALSE]
  cat("[DIMS] tcga_nbl_counts dimensions: ", dim(tcga_nbl_counts), "\n")
  return(tcga_nbl_counts)
}

#-----------DSMZ loading function-----------
# count - dsmz count data after buld_dsmz_matrix function
load_dsmz_breast_data <- function(counts, meta_path, dsmz_organ_filter) {
  cat("[INFO] Loading DSMZ breast data...\n")
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

#-----------PAM50 subtype filtering-----------
load_tcga_nbl_data_by_pam50 <- function(tcga_nbl_counts, tcga_nbl_subtype_metadata) {
  if (is.character(tcga_nbl_subtype_metadata)) {
    tcga_nbl_subtype_df <- read.csv(tcga_nbl_subtype_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(tcga_nbl_subtype_metadata)) {
    tcga_nbl_subtype_df <- tcga_nbl_subtype_metadata
  } else {
    stop("tcga_nbl_subtype_metadata must be a CSV file path (character) or data.frame")
  }
  
  stopifnot("sample_barcode" %in% colnames(tcga_nbl_subtype_df),
            "PAM50_Subtype" %in% colnames(tcga_nbl_subtype_df))
  
  common_samples <- intersect(colnames(tcga_nbl_counts), tcga_nbl_subtype_df$sample_barcode)
  if (length(common_samples) == 0) {
    stop("No matching samples between nbl_counts and sample_barcode in metadata")
  }
  
  matched_meta <- tcga_nbl_subtype_df[tcga_nbl_subtype_df$sample_barcode %in% common_samples, , drop = FALSE]
  tcga_subtype_counts <- tcga_nbl_counts[, common_samples, drop = FALSE]
  
  labels <- setNames(as.character(matched_meta$PAM50_Subtype), matched_meta$sample_barcode)
  labels <- labels[colnames(tcga_subtype_counts)]
  
  keep <- !is.na(labels)
  if (any(!keep)) {
    n_dropped <- sum(!keep)
    message(sprintf("Dropping %d sample(s) with missing PAM50_Subtype", n_dropped))
    tcga_subtype_counts <- tcga_subtype_counts[, keep, drop = FALSE]
    labels <- labels[keep]
  }
  
  stopifnot(identical(colnames(tcga_subtype_counts), names(labels)))
  cat("[DIMS] tcga_subtype_counts dimensions: ", dim(tcga_subtype_counts), "\n")
  cat("[DIMS] tcga_subtype_labels length: ", length(labels), "\n")
  return(list(tcga_subtype_counts = tcga_subtype_counts, labels = labels))
}