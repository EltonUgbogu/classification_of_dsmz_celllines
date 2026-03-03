# tcga_counts - tcga$counts from load_tcga_data function
# dsmz_counts - dsmz$counts from load_dsmz_data function
# min_shared - minimum number of shared genes optional -> required 
# min_shared - it can also be used to select features (genes)

common_genes <- function(tcga_counts, dsmz_counts, min_shared=NULL) {
  cat("[INFO] Finding common genes...\n")
  cat(sprintf("[INFO] TCGA-BRCA count matrix: %d genes x %d sample\n", nrow(tcga_counts), ncol(tcga_counts)))  # Print dimensions of raw DSMZ data
  cat(sprintf("[INFO] DSMZ-BREAST count matrix: %d genes x %d sample\n", nrow(dsmz_counts), ncol(dsmz_counts)))  # Print dimensions of raw DSMZ data
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


#' @title Merge TCGA and DSMZ count matrices (full gene union)
#'
#' @description
#' Combines full TCGA and DSMZ count matrices (no gene filtering).
#' Duplicate gene names are collapsed by summing.
#' Returns merged matrix, batch factor, and original sample IDs.
#'
#' @param tcga_counts Numeric matrix (genes x TCGA samples)
#' @param dsmz_counts Numeric matrix (genes x DSMZ samples)
#' @param dataset_names Character vector of length 2: e.g., `c("TCGA", "DSMZ")`
#'
#' @return List with:
#'   - `Xc_raw`: merged count matrix (genes x samples)
#'   - `batch`: batch factor
#'   - `sample_id`: original sample names
#' @export
merge_counts <- function(tcga_counts, dsmz_counts, dataset_names) {
  cat("[INFO] Merging TCGA and DSMZ counts (full gene union)...\n")

  # -------------------------------------------------------------------------
  # 1. Input validation
  # -------------------------------------------------------------------------
  if (!is.matrix(tcga_counts) || !is.matrix(dsmz_counts))
    stop("Both tcga_counts and dsmz_counts must be matrices")
  if (!is.numeric(tcga_counts) || !is.numeric(dsmz_counts))
    stop("Count matrices must contain numeric values")
  if (length(dataset_names) < 2)
    stop("dataset_names must be a vector of length >= 2")
  if (ncol(tcga_counts) == 0 || ncol(dsmz_counts) == 0)
    stop("One or both count matrices have zero samples")

  # -------------------------------------------------------------------------
  # 2. Direct merge (full gene sets)
  # -------------------------------------------------------------------------
  Xc_raw <- cbind(tcga_counts, dsmz_counts)

  # -------------------------------------------------------------------------
  # 3. Collapse duplicate gene names (sum across samples)
  # -------------------------------------------------------------------------
  if (any(duplicated(rownames(Xc_raw)))) {
    dup_genes <- rownames(Xc_raw)[duplicated(rownames(Xc_raw))]
    cat(sprintf("[INFO] Collapsing %d duplicate gene(s): %s\n",
                length(dup_genes),
                paste(head(dup_genes, 5), collapse = ", ")))
    Xc_raw <- rowsum(Xc_raw, group = rownames(Xc_raw))
  }

  # -------------------------------------------------------------------------
  # 4. Build batch factor
  # -------------------------------------------------------------------------
  batch <- factor(
    c(rep(dataset_names[1], ncol(tcga_counts)),
      rep(dataset_names[2], ncol(dsmz_counts))),
    levels = dataset_names[1:2]
  )
  names(batch) <- colnames(Xc_raw)

  # -------------------------------------------------------------------------
  # 5. Preserve original sample IDs
  # -------------------------------------------------------------------------
  sample_id <- colnames(Xc_raw)
  names(sample_id) <- sample_id

  # -------------------------------------------------------------------------
  # 6. Final matrix setup
  # -------------------------------------------------------------------------
  storage.mode(Xc_raw) <- "double"
  cat(sprintf("[DIMS] Merged data: %d genes x %d samples (%d TCGA + %d DSMZ)\n",
              nrow(Xc_raw), ncol(Xc_raw), ncol(tcga_counts), ncol(dsmz_counts)))

  # -------------------------------------------------------------------------
  # 7. Return
  # -------------------------------------------------------------------------
  list(
    Xc_raw     = Xc_raw,
    batch      = batch,
    sample_id  = sample_id
  )
}
