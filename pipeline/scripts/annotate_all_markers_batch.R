#!/usr/bin/env Rscript

# ============================================================
# Batch annotate all DESeq2 marker files with gene symbols
# ============================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
})

# Find all marker files
find_marker_files <- function() {
  cmd <- "find results/unsupervised -type f -name '*.tsv' | grep -E '(deseq2_markers|component_consensus_markers_v2)' | grep -E '(UP|DOWN)\\.(sig|core|extended|all)\\.tsv$'"
  files <- system(cmd, intern = TRUE)
  files
}

# Annotate a single file
annotate_file <- function(file_path) {
  tryCatch({
    # Read file
    df <- read.delim(file_path, sep = "\t", stringsAsFactors = FALSE)
    
    # Skip if already annotated
    if ("symbol" %in% colnames(df)) {
      return(list(success = FALSE, message = "Already annotated"))
    }
    
    # Find gene_id column
    gene_col <- NULL
    for (col in c("gene_id", "GeneID", "gene")) {
      if (col %in% colnames(df)) {
        gene_col <- col
        break
      }
    }
    
    if (is.null(gene_col)) {
      return(list(success = FALSE, message = "No gene_id column found"))
    }
    
    # Strip version from Ensembl IDs
    df$ensembl_id <- sub("\\..*$", "", df[[gene_col]])
    unique_ids <- unique(df$ensembl_id[!is.na(df$ensembl_id)])
    
    if (length(unique_ids) == 0) {
      return(list(success = FALSE, message = "No Ensembl IDs found"))
    }
    
    # Map to symbols
    map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "GENENAME")
    )
    
    # Merge
    df_annotated <- merge(
      df,
      map,
      by.x = "ensembl_id",
      by.y = "ENSEMBL",
      all.x = TRUE
    )
    
    # Rename SYMBOL to symbol
    if ("SYMBOL" %in% colnames(df_annotated)) {
      df_annotated <- df_annotated %>%
        rename(symbol = SYMBOL)
    }
    
    # Reorder columns: put symbol early
    if ("symbol" %in% colnames(df_annotated)) {
      other_cols <- setdiff(colnames(df_annotated), c("ensembl_id", "symbol", "GENENAME"))
      df_annotated <- df_annotated[, c(gene_col, "ensembl_id", "symbol", "GENENAME", other_cols)]
      df_annotated <- df_annotated[, !is.na(colnames(df_annotated))]
    }
    
    # Write output (overwrite original)
    write.table(df_annotated, file_path, sep = "\t", quote = FALSE, row.names = FALSE)
    
    n_mapped <- sum(!is.na(df_annotated$symbol))
    n_total <- nrow(df_annotated)
    
    return(list(
      success = TRUE,
      message = sprintf("Annotated %d/%d genes (%.1f%%)", n_mapped, n_total, 100*n_mapped/n_total)
    ))
  }, error = function(e) {
    return(list(success = FALSE, message = paste0("Error: ", e$message)))
  })
}

# Main execution
main <- function() {
  cat("[INFO] Finding marker files...\n")
  
  # Find all marker files
  marker_files <- find_marker_files()
  
  if (length(marker_files) == 0) {
    cat("[WARN] No marker files found\n")
    return(1)
  }
  
  cat(sprintf("[INFO] Found %d marker files\n", length(marker_files)))
  
  # Process each file
  success_count <- 0
  skip_count <- 0
  error_count <- 0
  
  for (i in seq_along(marker_files)) {
    file_path <- marker_files[i]
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(marker_files), basename(file_path)))
    
    result <- annotate_file(file_path)
    
    if (result$success) {
      cat(sprintf("  [OK] %s\n", result$message))
      success_count <- success_count + 1
    } else if (grepl("Already annotated", result$message)) {
      cat(sprintf("  [SKIP] %s\n", result$message))
      skip_count <- skip_count + 1
    } else {
      cat(sprintf("  [WARN] %s\n", result$message))
      error_count <- error_count + 1
    }
  }
  
  cat("\n[SUMMARY]\n")
  cat(sprintf("  Successfully annotated: %d\n", success_count))
  cat(sprintf("  Skipped (already annotated): %d\n", skip_count))
  cat(sprintf("  Errors/warnings: %d\n", error_count))
  cat(sprintf("  Total processed: %d\n", length(marker_files)))
  
  return(0)
}

if (!interactive()) {
  quit(status = main())
}
