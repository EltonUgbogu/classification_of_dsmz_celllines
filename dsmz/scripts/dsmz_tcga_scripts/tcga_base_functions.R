#tcga_counts - tcga$counts from load_tcga_data function

# Remove Ensembl version suffixes and sum duplicates
harmonize_gene_ids <- function(tcga_counts) {           # Define function to harmonize TCGA gene IDs
  cat("[INFO] Harmonizing TCGA gene IDs...\n")         # Print message indicating gene ID harmonization
  rownames(tcga_counts) <- sub("\\..*$","", rownames(tcga_counts))         # Remove version suffixes from Ensembl IDs
  if (any(duplicated(rownames(tcga_counts)))) {
  dup_genes <- rownames(tcga_counts)[duplicated(rownames(tcga_counts))]
  cat(sprintf("[INFO] Collapsing %d duplicate TCGA genes: %s\n", 
              length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))
  tcga_counts <- rowsum(tcga_counts, rownames(tcga_counts), reorder = TRUE)
  }                
  storage.mode(tcga_counts) <- "double"                # Convert count matrix to double precision
  tcga_counts                                          # Return harmonized count matrix
}


# Load tumor purity vector from SummarizedExperiment or external TSV
# tcga_se - tcga$se from load_tcga_data function
# purity_file - path to purity TSV file found in config -> optional
# return purity vector
verify_purity_data <- function(tcga_se, purity_file=NULL) {  
  p <- NULL                                            # Initialize purity vector as NULL
  if ("purity" %in% colnames(colData(tcga_se))) {      # Check if purity data is in SummarizedExperiment
    p <- as.numeric(colData(tcga_se)$purity)           # Extract purity values as numeric
    names(p) <- colnames(assay(tcga_se))               # Assign sample names to purity values
    cat("[INFO] Using purity from SummarizedExperiment\n")  # Print message indicating source
  } else if (!is.null(purity_file) && file.exists(purity_file)) {  # Check if external purity file is provided and exists
    tab <- read.delim(purity_file, stringsAsFactors=FALSE)  # Read purity TSV file
    stopifnot(all(c("sample","purity") %in% colnames(tab)))  # Verify required columns in TSV
    p <- setNames(tab$purity, tab$sample)              # Create named purity vector
    cat("[INFO] Using external purity file\n")         # Print message indicating source
  } else {
    cat("[INFO] No purity data available\n")           # Print message if no purity data
  }
  p                                                    # Return purity vector
}

# Account for genes that are strongly negatively correlated with purity (likely infiltration) -> not expecting significant correlations
# tcga_norm - Normmalized TCGA counts matrix -> VST (Variance Stabilizing Transformation)
# purity - purity vector from load_purity_data function
# return adjusted TCGA counts matrix
purity_adjust_log <- function(tcga_norm, purity=NULL) {   
  if (is.null(purity)) return(tcga_norm)                # Return input if no purity data
  out <- tcga_norm                                      # Initialize output matrix
  cat("[INFO] Applying purity adjustment on VST...\n")  # Print message indicating adjustment
  keep_samp <- intersect(colnames(tcga_norm), names(purity)[!is.na(purity)])  # Find samples with valid purity values
  M <- tcga_norm[, keep_samp, drop=FALSE]               # Subset VST matrix to valid samples
  pu <- purity[keep_samp]                              # Subset purity to valid samples
  if (ncol(M) < 3) {                                  # Check if enough samples for adjustment
    cat("[WARN] Too few samples with purity to adjust; returning input.\n")  # Warn if insufficient samples
    return(out)                                   # Return input matrix
  }
  # Identify & drop genes strongly NEG correlated with purity (infiltration)
  ct <- apply(M, 1, function(v) suppressWarnings(cor.test(as.numeric(v), pu, method="spearman")))  # Compute Spearman correlations
  pv <- vapply(ct, `[[`, numeric(1), "p.value")        # Extract p-values
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))  # Extract correlation coefficients
  padj <- p.adjust(pv, "BH")                           # Adjust p-values for multiple testing (Benjamini-Hochberg)
  drop <- names(which(padj < 0.01 & rh < -0.4))        # Identify genes with significant negative correlations
  if (length(drop)) {                                  # Check if genes need to be dropped
    cat(sprintf("[INFO] Removing %d purity-associated genes\n", length(drop)))  # Print number of dropped genes
    M <- M[setdiff(rownames(M), drop), , drop=FALSE]   # Remove identified genes
  }
  # Regress out infiltration (1 - purity) on the log scale
  infilt <- 1 - pu                                     # Calculate infiltration (1 - purity)
  design <- model.matrix(~ infilt)                      # Create design matrix for regression
  fit <- lmFit(M, design)                              # Fit linear model using limma
  beta <- fit$coefficients[, "infilt", drop=FALSE]     # Extract infiltration coefficients
  Madj <- as.matrix(M) - beta %*% t(infilt)            # Adjust matrix by subtracting infiltration effect
  out <- tcga_norm                                      # Initialize output matrix
  common_g <- intersect(rownames(Madj), rownames(out))  # Find common genes
  out[common_g, colnames(Madj)] <- Madj[common_g, ]    # Update adjusted values
  list(V = out, dropped = drop)                        # Return adjusted matrix and dropped genes
}
