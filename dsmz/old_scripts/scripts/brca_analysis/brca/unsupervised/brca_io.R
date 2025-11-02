config <- list(
  tcga_se_rds = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds",
  dsmz_rds = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",
  tcga_project_filter = c("TCGA-BRCA"),
  dsmz_organ_filter = c("Breast"),
  tcga_brca_subtype_metadata = "/Users/eltonugbogu/classification_of_dsmz_celllines/dsmz/old_scripts/scripts/brca_analysis/brca/unsupervised/TCGA_BRCA_project_subtype.csv",
  outdir = "/home/chu25/dsmz/dsmz_tcga_brca/results"
)

#-----------TCGA-BRCA loading function-----------
load_tcga_brca_data <- function(file_path, tcga_project_filter) {
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
  
  keep_tcga_brca <- project_ids %in% tcga_project_filter
  if (sum(keep_tcga_brca) == 0) {
    stop("No TCGA samples match the project filter: ", paste(tcga_project_filter, collapse = ", "))
  }
  
  tcga_brca_counts <- tcga_counts[, keep_tcga_brca, drop = FALSE]
  cat("[DIMS] tcga_brca_counts dimensions: ", dim(tcga_brca_counts), "\n")
  return(tcga_brca_counts)
}

#-----------DSMZ loading function-----------
load_dsmz_breast_data <- function(counts_path, meta_path, dsmz_organ_filter) {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) stop(sprintf("[ERROR] DSMZ RDS not found: %s", counts_path))
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  
  dsmz_raw <- readRDS(counts_path)
  dsmz_meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  
  stopifnot(all(c("Ensembl_ID", "gene_name") %in% colnames(dsmz_raw)))
  stopifnot("sample_name" %in% names(dsmz_meta))
  
  cat(sprintf("[INFO] DSMZ table: %d rows x %d cols\n", nrow(dsmz_raw), ncol(dsmz_raw)))
  
  dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)
  matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_raw))
  
  if (length(matched) == 0) {
    stop("No overlapping sample names between DSMZ counts and metadata")
  }
  
  dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched)
  dsmz_counts <- dsmz_raw[, dsmz_meta$sample_id, drop = FALSE]
  
  stopifnot("organ" %in% colnames(dsmz_meta))
  if (is.null(dsmz_organ_filter) || length(dsmz_organ_filter) == 0) {
    stop("dsmz_organ_filter must not be null or empty")
  }
  
  keep_dsmz_brca <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
  if (length(keep_dsmz_brca) == 0) {
    stop("No DSMZ samples match the organ filter: ", paste(dsmz_organ_filter, collapse = ", "))
  }
  
  dsmz_brca_counts <- dsmz_counts[, keep_dsmz_brca, drop = FALSE]
  cat("[DIMS] dsmz_brca_counts dimensions: ", dim(dsmz_brca_counts), "\n")
  return(dsmz_brca_counts)
}

#-----------PAM50 subtype filtering-----------
load_tcga_brca_data_by_pam50 <- function(tcga_brca_counts, tcga_brca_subtype_metadata) {
  if (is.character(tcga_brca_subtype_metadata)) {
    tcga_brca_subtype_df <- read.csv(tcga_brca_subtype_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(tcga_brca_subtype_metadata)) {
    tcga_brca_subtype_df <- tcga_brca_subtype_metadata
  } else {
    stop("tcga_brca_subtype_metadata must be a CSV file path (character) or data.frame")
  }
  
  stopifnot("sample_barcode" %in% colnames(tcga_brca_subtype_df),
            "PAM50_Subtype" %in% colnames(tcga_brca_subtype_df))
  
  common_samples <- intersect(colnames(tcga_brca_counts), tcga_brca_subtype_df$sample_barcode)
  if (length(common_samples) == 0) {
    stop("No matching samples between brca_counts and sample_barcode in metadata")
  }
  
  matched_meta <- tcga_brca_subtype_df[tcga_brca_subtype_df$sample_barcode %in% common_samples, , drop = FALSE]
  tcga_subtype_counts <- tcga_brca_counts[, common_samples, drop = FALSE]
  
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