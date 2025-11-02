# tcga_counts - tcga$counts from load_tcga_data function
# dsmz_counts - dsmz$counts from load_dsmz_data function
# min_shared - minimum number of shared genes optional -> required 
# min_shared - it can also be used to select features (genes)

common_genes <- function(tcga_counts, dsmz_counts, min_shared=NULL) {
  cat("[INFO] Finding common genes...\n")
  cat(sprintf("[INFO] TCGA count matrix: %d genes x %d sample\n", nrow(tcga_counts), ncol(tcga_counts)))  # Print dimensions of raw DSMZ data
  cat(sprintf("[INFO] DSMZ count matrix: %d genes x %d sample\n", nrow(dsmz_counts), ncol(dsmz_counts)))  # Print dimensions of raw DSMZ data
  print(head(rownames(tcga_counts)))
  print(head(rownames(dsmz_counts)))
  common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))  # Find common genes
  cat(sprintf("[INFO] Shared genes: %d\n", length(common)))  # Print number of shared genes
  if(!is.null(min_shared)) {
    if (length(common) < min_shared) cat("[WARN] Few shared genes; check Ensembl IDs\n")  # Warn if too few shared genes
  }
  tcga_counts <- tcga_counts[common, , drop=FALSE]    # Subset TCGA counts to common genes
  dsmz_counts <- dsmz_counts[common, , drop=FALSE]    # Subset DSMZ counts to common genes
  list(tcga_counts = tcga_counts, dsmz_counts = dsmz_counts)
}

# tcga_counts - tcga$counts from load_tcga_data function
# dsmz_counts - dsmz$counts from load_dsmz_data function
# dataset_names - vector of dataset names e.g. c("TCGA", "DSMZ") here it is of length 2
# return list with merged counts matrix and batch factor

merge_counts <- function(tcga_counts, dsmz_counts, dataset_names) {
  cat("[INFO] Merging TCGA and DSMZ counts...\n")
  # Combine TCGA and DSMZ counts matrices after common genes have been identified
  Xc_raw <- cbind(tcga_counts, dsmz_counts) 
  # ensure dataset_name is a vector of length n ; n=2
  if(length(dataset_names) != 2) stop("dataset_names must be a vector of length >= 2")
  batch_levels <- dataset_names
  batch <- factor(c(rep(batch_levels[1], ncol(tcga_counts)), rep(batch_levels[2], ncol(dsmz_counts))), levels = batch_levels) # define batch factor with dataset names
  
  if (any(duplicated(rownames(Xc_raw)))) {
    dup_genes <- rownames(Xc_raw)[duplicated(rownames(Xc_raw))]                  # get duplicate genes
    cat(sprintf("[INFO] Collapsing %d duplicate genes in merged matrix: %s\n",    # print number of duplicate genes and first 5 duplicate genes
                length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))    # print duplicate genes
    Xc_raw <- rowsum(Xc_raw, rownames(Xc_raw), reorder = TRUE)                   # collapse duplicate genes by summing counts for each gene
  }
  storage.mode(Xc_raw) <- "double"                                  # convert counts to double precision
  list(Xc_raw = Xc_raw, batch = batch)                              # return list with merged counts matrix and batch factor
}   

