


# Build a DSMZ count matrix with Ensembl IDs as row names
# and sample names as column names
build_dsmz_matrix <- function(dsmz_raw) {
  cat("[INFO] Building DSMZ count matrix...\n")
  annot <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot)
  # Convert sample columns to numeric
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) {
    as.numeric(as.character(x))
  })  
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE]) 
  # Set row names: Ensembl_ID without version
  rownames(M) <- sub("\\..*$", "", dsmz_raw$Ensembl_ID) 
  # Collapse duplicates by summing
  if (any(duplicated(rownames(M)))) {
    dup_genes <- rownames(M)[duplicated(rownames(M))]
    cat(sprintf("[INFO] Collapsing %d duplicate DSMZ genes: %s\n",
                length(dup_genes), paste(head(dup_genes, 5), collapse = ", ")))
    M <- rowsum(M, group = rownames(M), reorder = TRUE)
  }
  # Identify and drop all-NA columns
  keep <- colSums(!is.na(M)) > 0
  if (!all(keep)) {
    cat(sprintf("[WARN] Dropping %d DSMZ columns that are all NA\n", sum(!keep)))
    M <- M[, keep, drop = FALSE]
  }
  storage.mode(M) <- "double"
  # Print first 5 sample names (colnames) 
  cat("[INFO] First 5 sample names (colnames):\n")
  cat(paste("   ", head(colnames(M), 5), collapse = "\n    "), "\n")
  # Print first 5 gene IDs (rownames) 
  cat("[INFO] First 5 gene IDs (rownames):\n")
  cat(paste("   ", head(rownames(M), 5), collapse = "\n    "), "\n")
  cat(sprintf("[INFO] Final matrix: %d genes x %d samples\n", nrow(M), ncol(M)))
  M
}


#col_name - column name upon which to align the metadata to the counts
# Align DSMZ metadata to count matrix; write out unmatched for audit
align_counts_metadata <- function(meta_data=dsmz_meta, counts_data=dsmz_counts, col_name="sample_name", outdir) {  # Define function to align metadata and counts
  cat("[INFO] Aligning counts to metadata...\n")       # Print message indicating alignment
  meta_data <- meta_data %>% mutate(sample_id = col_name)  # Create sample_id column from sample_name
  matched <- intersect(meta_data$sample_id, colnames(counts_data))  # Find common samples between metadata and counts
  md_only <- setdiff(meta_data$sample_id, colnames(counts_data))  # Find metadata samples not in counts
  ct_only <- setdiff(colnames(counts_data), meta_data$sample_id)  # Find count samples not in metadata
  if (length(md_only))                                 # Check if there are unmatched metadata samples
    write.csv(tibble(sample_name = md_only),            # Save unmatched metadata samples to CSV
              file.path(outdir, "unmatched_metadata_samples.csv"),
              row.names = FALSE)
  if (length(ct_only))                                 # Check if there are unmatched count samples
    write.csv(tibble(sample_name = ct_only),            # Save unmatched count samples to CSV
              file.path(outdir, "unmatched_count_columns.csv"),
              row.names = FALSE)
  meta_data <- meta_data %>% filter(sample_id %in% matched)  # Filter metadata to matched samples
  counts_data <- counts_data[, meta_data$sample_id, drop=FALSE]  # Subset counts to matched samples
  stopifnot(identical(colnames(counts_data), meta_data$sample_id))  # Verify counts and metadata alignment
  cat(sprintf("[INFO] Matched: %d | Unmatched(meta): %d | Unmatched(counts): %d\n",
              length(matched), length(md_only), length(ct_only)))  # Print alignment summary
  list(meta = meta_data, counts = counts_data)         # Return list with aligned metadata and counts
}








