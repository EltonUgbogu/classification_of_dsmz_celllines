# I/O: load/save helpers and dimension reporting

config <- list(
  tcga_project_filter = c("TCGA-BRCA"),
  dsmz_organ_filter = c("Breast"),
  tcga_brca_subtype_metadata = c("/Users/eltonugbogu/classification_of_dsmz_celllines/dsmz/old_scripts/scripts/brca_analysis/brca/unsupervised/TCGA_BRCA_project_subtype.csv")
)

# --- TCGA-BRCA loading function ---
# tcga_counts - TCGA expression matrix for all projects
# tcga_se - TCGA SummarizedExperiment for all projects
# tcga_project_filter - TCGA project filter
# return: tcga_brca_counts - TCGA expression matrix for all projects matching the project filter
load_tcga_brca_data <- function(tcga_counts, tcga_se, tcga_project_filter) {
  if (is.null(tcga_project_filter)) {
    stop("tcga_project_filter is null")
  }
  
  # Get project_id for samples in tcga_counts
  sample_names <- colnames(tcga_counts)
  project_ids <- as.character(colData(tcga_se)[sample_names, "project_id"])
  
  # Keep only samples with matching project_id
  keep_tcga_brca <- project_ids %in% tcga_project_filter
  
  if (sum(keep_tcga_brca) == 0) {
    stop("No TCGA samples match the project filter: ", paste(tcga_project_filter, collapse = ", "))
  }
  
  tcga_brca_counts <- tcga_counts[, keep_tcga_brca, drop = FALSE]
  cat("tcga_brca_counts dimensions: ", dim(tcga_brca_counts), "\n") # print dimensions of tcga_brca_counts
  
  return(tcga_brca_counts)
}



# --- DSMZ-BRCA loading function ---
# dsmz_counts - DSMZ expression matrix for all samples
# dsmz_meta - DSMZ metadata
# dsmz_organ_filter - DSMZ organ filter
# return: dsmz_brca_counts - DSMZ expression matrix for all samples matching the organ filter
load_dsmz_brca_data <- function(dsmz_counts, dsmz_meta, dsmz_organ_filter) {
  if (is.null(dsmz_organ_filter)) {
    stop("dsmz_organ_filter is null")
  }
  
  stopifnot("organ" %in% colnames(dsmz_meta))
  
  keep_dsmz_brca <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
  
  if (length(keep_dsmz_brca) == 0) {
    stop("No DSMZ samples match the organ filter: ", paste(dsmz_organ_filter, collapse = ", "))
  }
  
  dsmz_brca_counts <- dsmz_counts[, keep_dsmz_brca, drop = FALSE]
  log_dims("DSMZ (after organ subset)", dsmz_brca_counts)
  
  return(dsmz_brca_counts)
}

# brca_counts - tcga_brca_counts from load_tcga_brca_data function
# brca_subtype_metadata - brca_subtype_metadata 
# subset TCGA_BRCA based on PAM50 subtype annotation
load_tcga_brca_data_by_pam50 <- function(brca_counts,
                                         brca_subtype_metadata) {
  
  # Load metadata (CSV or data.frame)
  if (is.character(brca_subtype_metadata)) {
    brca_subtype_df <- read.csv(brca_subtype_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(brca_subtype_metadata)) {
    brca_subtype_df <- brca_subtype_metadata
  } else {
    stop("brca_subtype_metadata must be a CSV file path (character) or data.frame")
  }
  
  # Validate required columns
  stopifnot("sample_barcode" %in% colnames(brca_subtype_df),
            "PAM50_Subtype" %in% colnames(brca_subtype_df))
  
  # Find overlapping samples between counts and metadata
  common_samples <- intersect(colnames(brca_counts), brca_subtype_df$sample_barcode)
  
  if (length(common_samples) == 0) {
    stop("No matching samples between brca_counts and sample_barcode in metadata")
  }
  
  # Subset metadata to only matched samples
  matched_meta <- brca_subtype_df[brca_subtype_df$sample_barcode %in% common_samples, , drop = FALSE]
  
  # Subset counts matrix to matched samples (preserve original order)
  tcga_subtype_counts <- brca_counts[, common_samples, drop = FALSE]
  
  # Create labels vector with sample_barcode as names
  labels <- setNames(as.character(matched_meta$PAM50_Subtype), matched_meta$sample_barcode)
  
  # Reorder labels to match column order of tcga_subtype_counts
  labels <- labels[colnames(tcga_subtype_counts)]
  
  # Drop samples with NA in PAM50 subtype
  keep <- !is.na(labels)
  if (any(!keep)) {
    n_dropped <- sum(!keep)
    message(sprintf("Dropping %d sample(s) with missing PAM50_Subtype", n_dropped))
    tcga_subtype_counts <- tcga_subtype_counts[, keep, drop = FALSE]
    labels <- labels[keep]
  }
  
  # Final sanity check: ensure names match and are in sync
  stopifnot(identical(colnames(tcga_subtype_counts), names(labels)))
  
  cat("brca_subtype_counts dimensions: ", dim(tcga_subtype_counts), "\n")
  cat("labels length: ", length(labels), "\n")
  
  return(list(brca_subtype_counts = tcga_subtype_counts, labels = labels))
}

