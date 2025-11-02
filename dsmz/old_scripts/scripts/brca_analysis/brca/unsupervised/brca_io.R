# I/O: load/save helpers and dimension reporting

config <- list(
  tcga_project_filter = c("TCGA-BRCA"),
  dsmz_organ_filter = c("Breast"),
  tcga_brca_subtype_metadata = c("/Users/eltonugbogu/classification_of_dsmz_celllines/dsmz/old_scripts/scripts/brca_analysis/brca/unsupervised/TCGA_BRCA_project_subtype.csv")
)


if (!is.null(dsmz_organ_filter)) {
  keep_ids <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
  dsmz_meta <- dsmz_meta[dsmz_meta$sample_id %in% keep_ids, , drop = FALSE]
  dsmz_counts <- dsmz_counts[, keep_ids, drop = FALSE]
}
# tcga_se - tcga_raw$se from load_tcga_data function
# tcga_counts - tcga_raw$counts from load_tcga_data function
# tcga_project_filter - subset of tcga organ/tissue projects_id
load_tcga_brca_data <- function(tcga_counts, tcga_se, tcga_project_filter) {
  if (!is.null(tcga_project_filter)) {
    keep_tcga_brca <- colnames(tcga_counts) %in% colnames(assay(tcga_se)) &
               as.character(colData(tcga_se)[colnames(tcga_counts), "project_id"]) %in% tcga_project_filter
    tcga_brca_counts <- tcga_counts[, keep_tcga_brca, drop = FALSE]
    log_dims("TCGA (after project subset)", tcga_brca_counts)
    stopifnot(ncol(tcga_brca_counts) > 0L)
    return(tcga_brca_counts)
  }
}

# dsmz_counts - dsmz_raw$counts from load_dsmz_data function
# dsmz_meta - dsmz_raw$meta from load_dsmz_data function
# dsmz_organ_filter - subset of dsmz organ/tissue
load_dsmz_brca_data <- function(dsmz_counts, dsmz_meta, dsmz_organ_filter) {
  stopifnot("organ" %in% colnames(dsmz_meta))
  if (!is.null(dsmz_organ_filter)) {
    keep_dsmz_brca <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter]
    dsmz_brca_counts <- dsmz_counts[, keep_dsmz_brca, drop = FALSE]
    log_dims("DSMZ (after organ subset)", dsmz_brca_counts)
    stopifnot(ncol(dsmz_brca_counts) > 0L)
    return(dsmz_brca_counts)
  }
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




write_dimension_report <- function(outdir, dims_list) {
  dim_report <- file.path(outdir, "dimension_report.txt")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  con <- file(dim_report, open = "wt")
  on.exit(close(con), add = TRUE)
  lapply(names(dims_list), function(nm) write_dims_line(con, nm, dims_list[[nm]]))
}

safe_rds_save <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}

log_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark=","), format(ns, big.mark=",")))
}

write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con)
}

## form clustering.R
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

safe_rds_save <- function(object, file) {
  ensure_dir(dirname(file))
  saveRDS(object, file = file)
}

run_pca_and_save <- function(V_adj, outdir, n_hvg = 3000, max_pc = 30) {
  ensure_dir(outdir)
  pcs <- make_pcs_matrix(V_adj, n_hvg = n_hvg, max_pc = max_pc)
  safe_rds_save(pcs, file.path(outdir, "pca_objects.rds"))
  pcs
}