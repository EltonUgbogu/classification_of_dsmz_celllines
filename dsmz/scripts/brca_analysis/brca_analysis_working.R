#' @file brca_vst_dual_analysis.R
#' @brief Performs integrative analysis of TCGA and DSMZ breast cancer datasets using VST normalization, batch correction, and clustering.
#' @details This script processes RNA-seq data from TCGA and DSMZ, applies variance stabilizing transformation (VST), performs batch correction,
#'          subsets to PAM50 genes, and conducts clustering (k-means, hierarchical, HDBSCAN) and differential expression analysis to identify
#'          breast cancer subtypes. Results are saved as CSVs, PDFs, and RDS files.
#' @author [Your Name]
#' @date 2025-10-23

#' @section Setup:
#' Configures R library paths and loads required packages for data processing, visualization, and statistical analysis.
.libPaths(c("/home/chu25/miniconda3/envs/r_sva_fix/lib/R/library", .libPaths()))

#' Suppress startup messages for cleaner console output.
# nolint start
suppressPackageStartupMessages({
  library(tidyverse)         #' @note Loads tidyverse for data manipulation and visualization.
  library(SummarizedExperiment) #' @note Loads SummarizedExperiment for handling genomic data.
  library(DESeq2)           #' @note Loads DESeq2 for differential expression analysis.
  library(matrixStats)      #' @note Loads matrixStats for matrix operations (e.g., rowIQRs).
  library(pheatmap)         #' @note Loads pheatmap for heatmap visualization.
  library(uwot)             #' @note Loads uwot for UMAP dimensionality reduction.
  library(cluster)          #' @note Loads cluster for clustering algorithms (e.g., silhouette).
  library(fpc)              #' @note Loads fpc for clustering validation metrics.
  library(clusterCrit)      #' @note Loads clusterCrit for clustering criteria (e.g., Davies-Bouldin).
  library(limma)            #' @note Loads limma for linear modeling (e.g., purity adjustment).
  library(sva)              #' @note Loads sva for batch effect correction (e.g., ComBat).
  library(dbscan)           #' @note Loads dbscan for HDBSCAN clustering.
  library(PMCMRplus)        #' @note Loads PMCMRplus for post-hoc statistical tests.
  library(AnnotationDbi)    #' @note Loads AnnotationDbi for gene annotation.
  library(org.Hs.eg.db)     #' @note Loads org.Hs.eg.db for human gene annotations.
})

```{r configuration}
#' @section Configuration:
#' Defines file paths and key parameters for the analysis.
tcga_se_rds   <- "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds" #' @note Path to TCGA SummarizedExperiment RDS file.
dsmz_rds      <- "/home/chu25/data/dsmz/DSMZ_count_gene.rds" #' @note Path to DSMZ count gene RDS file.
dsmz_meta_csv <- "/home/chu25/data/dsmz/DSMZ_metadata.csv" #' @note Path to DSMZ metadata CSV file.
outdir <- "/home/chu25/dsmz/results/brca_vst_dual" #' @note Output directory for results.
dir.create(outdir, showWarnings = FALSE, recursive = TRUE) #' @note Creates output directory if it doesn't exist.

tcga_project_filter <- c("TCGA-BRCA") #' @note Filters TCGA data to include only BRCA project.
dsmz_organ_filter   <- c("Breast") #' @note Filters DSMZ data to include only breast samples.
marker_symbols <- c("ESR1","PGR","ERBB2","MKI67","KRT5","KRT14","KRT17","FOXC1","GRB7",
                    "BCL2","EGFR","CCNB1","BIRC5","MYBL2","KIF2C","KRT8","KRT18","GATA3","XBP1") #' @note List of marker genes for breast cancer subtypes.
pam50_symbols <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC6","CDC20",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT5","KRT14","KRT17","KRT23","MDM2","MELK",
  "MIA","MKI67","MLPH","MMP11","MYBL2","NAT1","ORC6L","PGR","PHGDH","PTTG1",
  "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T","WHSC1L1","KRT7","KRT15"
) #' @note List of PAM50 genes for subtype classification.
config <- list(batch_adjust_method = "combat") #' @note Configuration for batch adjustment method (ComBat).
```

```{r harmonize_ensembl_ids}
#' @section Helper Functions:
#' Defines utility functions for data processing and visualization.

#' @brief Strips Ensembl version numbers from gene IDs.
#' @param x Character vector of Ensembl IDs.
#' @return Character vector with version suffixes removed.
strip_ensver <- function(x) sub("\\..*$","", x)

#' @brief Converts input to a data frame with strings as factors disabled.
#' @param x Input object (e.g., matrix or data frame).
#' @return Data frame with stringsAsFactors set to FALSE.
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

#' @brief Creates a mapping of Ensembl IDs to gene symbols from DSMZ data.
#' @param dsmz_raw Data frame containing DSMZ raw data with Ensembl_ID and gene_name columns.
#' @return Data frame with ensembl, symbol, and symbol_unique columns, indexed by ensembl.
```

```{r map_EnsemblIDs_to_gene_symbols}
#' @brief Creates a mapping of Ensembl IDs to gene symbols from DSMZ data.
#' @param dsmz_raw Data frame containing DSMZ raw data with Ensembl_ID and gene_name columns.
#' @return Data frame with ensembl, symbol, and symbol_unique columns, indexed by ensembl.(Data frame mapping Ensembl IDs to gene symbols.)
make_ens2sym <- function(dsmz_raw) {
  df <- as_df(dsmz_raw) #' @note Converts input to data frame.
  ens2sym <- df %>%
    transmute(
      ensembl = strip_ensver(Ensembl_ID), #' @note Strips version from Ensembl IDs.
      hgnc_gene_symbol  = as.character(gene_name) #' @note Ensures gene names are characters.
    ) %>%
    filter(!is.na(ensembl) & nzchar(ensembl)) %>% #' @note Filters out invalid Ensembl IDs.
    arrange(ensembl, dplyr::desc(nchar(hgnc_gene_symbol)), hgnc_gene_symbol) %>% #' @note Sorts by Ensembl ID and symbol length.
    distinct(ensembl, .keep_all = TRUE) %>% #' @note Keeps unique Ensembl IDs.
    mutate(hgnc_gene_symbol_unique = make.unique(ifelse(is.na(hgnc_gene_symbol) | hgnc_gene_symbol == "", ensembl, hgnc_gene_symbol))) #' @note Creates unique symbols.
  ens2sym <- as.data.frame(ens2sym, stringsAsFactors = FALSE) #' @note Converts to data frame.
  rownames(ens2sym) <- ens2sym$ensembl #' @note Sets row names to Ensembl IDs.
  ens2sym #' @return Final mapping data frame.
}
```

```{r top_var_genes}
#' @brief Selects top variable genes based on interquartile range (IQR).
#' @param mat Numeric matrix with genes as rows and samples as columns.
#' @param n Number of top variable genes to select (default: 3000).
#' @return Character vector of top gene IDs.
top_var_genes <- function(mat, n = 3000) {
  iqr <- matrixStats::rowIQRs(mat) #' @note Computes IQR for each gene.
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))] #' @note Returns top n gene names sorted by IQR.
}
```

```{r zscore_by_gene}
#' @brief Applies z-score normalization to a matrix by gene.
#' @param mat Numeric matrix with genes as rows and samples as columns.
#' @return Z-score normalized matrix.
zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE) #' @note Computes row means.
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1 #' @note Computes row SDs, sets zero SDs to 1.
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/") #' @note Applies z-score transformation.
}
```

```{r safe_pdf}
#' @brief Safely creates a PDF file and evaluates an expression.
#' @param path File path for the PDF output.
#' @param expr Expression to evaluate for plotting.
#' @return Invisible result of the evaluated expression.
safe_pdf <- function(path, expr) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE) #' @note Creates output directory.
  pdf(path) #' @note Opens PDF device.
  on.exit({
    if (dev.cur() != 1) dev.off() #' @note Ensures device is closed on exit.
  }, add = TRUE)
  force(expr) #' @note Evaluates the plotting expression.
}
```

```{r build_dsmz_count_matrix}
#' @brief Builds a count matrix from DSMZ raw data.
#' @param dsmz_raw Data frame with DSMZ count data.
#' @return Numeric matrix with genes as rows and samples as columns.
build_dsmz_matrix <- function(dsmz_raw) {
  stopifnot(all(c("Ensembl_ID","gene_name") %in% names(dsmz_raw))) #' @note Checks required columns exist.
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version") #' @note Defines annotation columns.
  sample_cols <- setdiff(colnames(dsmz_raw), annot) #' @note Identifies sample columns.
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(v) as.numeric(as.character(v))) #' @note Converts sample columns to numeric.
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE]) #' @note Creates matrix from sample columns.
  rownames(M) <- strip_ensver(dsmz_raw$Ensembl_ID) #' @note Sets row names to stripped Ensembl IDs.
  if (any(duplicated(rownames(M)))) {
    dup_genes <- rownames(M)[duplicated(rownames(M))] #' @note Identifies duplicate genes.
    cat(sprintf("[INFO] Collapsing %d duplicate DSMZ genes: %s\n", 
                length(dup_genes), paste(head(dup_genes, 5), collapse=", "))) #' @note Logs duplicate genes.
    M <- rowsum(M, rownames(M), reorder = TRUE) #' @note Collapses duplicates by summing.
  }
  M #' @return Final count matrix.
}
```

```{r report_dims}
#' @brief Report dimensions of a matrix to the console.
#' @param tag Descriptive tag for the matrix.
#' @param mat Matrix to log dimensions for.
report_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_ #' @note Gets number of genes.
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_ #' @note Gets number of samples.
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark=","), format(ns, big.mark=","))) #' @note Prints dimensions.
}
```

```{r write_dims_line_to_file}
#' @brief Writes matrix dimensions to a file.
#' @param con File connection for output.
#' @param tag Descriptive tag for the matrix.
#' @param mat Matrix to write dimensions for.
write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_ #' @note Gets number of genes.
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_ #' @note Gets number of samples.
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con) #' @note Writes dimensions to file.
}
```

```{r compute_logCPM}
#' @brief Computes log2 counts per million (CPM).
#' @param x Numeric matrix of counts.
#' @param prior.count Prior count for log transformation (default: 1).
#' @param lib.size Optional library sizes (default: NULL).
#' @return Log2 CPM matrix.
logCPM <- function(x, prior.count = 1, lib.size = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x) #' @note Ensures input is a matrix.
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size) #' @note Computes log2 CPM using edgeR.
}
```

```{r load_purity_data}
#' @brief Loads tumor purity data for TCGA samples.
#' @param tcga_se SummarizedExperiment object with TCGA data.
#' @param purity_file Optional file path for external purity data.
#' @return Named numeric vector of purity values or NULL if unavailable.
load_purity_data <- function(tcga_se, purity_file = NULL) {
  p <- NULL
  if ("purity" %in% colnames(colData(tcga_se))) {
    p <- as.numeric(colData(tcga_se)$purity) #' @note Extracts purity from SummarizedExperiment.
    names(p) <- colnames(assay(tcga_se)) #' @note Sets sample names.
    cat("[INFO] Using purity from SummarizedExperiment\n") #' @note Logs source.
  } else if (!is.null(purity_file) && file.exists(purity_file)) {
    tab <- read.delim(purity_file, stringsAsFactors = FALSE) #' @note Reads external purity file.
    stopifnot(all(c("sample","purity") %in% colnames(tab))) #' @note Checks required columns.
    p <- setNames(tab$purity, tab$sample) #' @note Creates named purity vector.
    cat("[INFO] Using external purity file\n") #' @note Logs source.
  } else {
    cat("[INFO] No purity data available\n") #' @note Logs absence of purity data.
  }
  p #' @return Purity vector or NULL.
}
```

```{r purity_adjustment}
#' @brief Adjusts VST-transformed matrix data for tumor purity.
#' @param tcga_log VST-transformed matrix for TCGA samples.
#' @param purity Named numeric vector of purity values.
#' @return Adjusted VST matrix.
purity_adjust_log <- function(tcga_log, purity) {
  if (is.null(purity)) return(tcga_log) #' @note Returns input if no purity data.
  cat("[INFO] Applying purity adjustment on VST data...\n") #' @note Logs adjustment start.
  keep_samp <- intersect(colnames(tcga_log), names(purity)[!is.na(purity)]) #' @note Identifies samples with valid purity.
  if (length(keep_samp) < 3) {
    cat("[WARN] Too few samples with valid purity to adjust; returning input.\n") #' @note Warns if too few samples.
    return(tcga_log)
  }
  M <- tcga_log[, keep_samp, drop = FALSE] #' @note Subsets matrix to valid samples.
  pu <- purity[keep_samp] #' @note Subsets purity values.
  ct <- apply(M, 1, function(v) {
    if (var(v) == 0 || var(pu) == 0) return(list(p.value = 1, estimate = 0)) #' @note Handles zero-variance cases.
    suppressWarnings(cor.test(as.numeric(v), pu, method = "spearman")) #' @note Computes Spearman correlation.
  })
  pv <- vapply(ct, `[[`, numeric(1), "p.value") #' @note Extracts p-values.
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1)) #' @note Extracts correlation coefficients.
  pv[is.na(pv)] <- 1 #' @note Sets NA p-values to 1.
  padj <- p.adjust(pv, "BH") #' @note Applies BH correction.
  drop <- names(which(padj < 0.01 & rh < -0.4)) #' @note Identifies genes with significant negative correlation.
  if (length(drop)) {
    cat(sprintf("[INFO] Removing %d purity-associated genes\n", length(drop))) #' @note Logs removed genes.
    M <- M[setdiff(rownames(M), drop), , drop = FALSE] #' @note Removes correlated genes.
  }
  infilt <- 1 - pu #' @note Computes infiltration (1 - purity).
  design <- model.matrix(~ infilt) #' @note Creates design matrix for linear model.
  fit <- lmFit(M, design) #' @note Fits linear model using limma.
  beta <- fit$coefficients[, "infilt", drop = FALSE] #' @note Extracts infiltration coefficients.
  Madj <- as.matrix(M) - beta %*% t(infilt) #' @note Adjusts matrix for infiltration.
  out <- tcga_log #' @note Initializes output matrix.
  common_g <- intersect(rownames(Madj), rownames(out)) #' @note Finds common genes.
  out[common_g, colnames(Madj)] <- Madj[common_g, ] #' @note Updates adjusted values.
  out #' @return Adjusted matrix.
}
```



#' @brief Subsets a matrix to PAM50 genes, collapsing multiple Ensembl IDs per symbol.
#' @param mat Input matrix with Ensembl IDs as row names.
#' @param ens2sym Data frame mapping Ensembl IDs to gene symbols.
#' @param pam_symbols Vector of PAM50 gene symbols.
#' @param collapse_fun Method to collapse multiple Ensembl IDs ("sum", "mean", "max").
#' @return Matrix subset to PAM50 genes with symbols as row names.
subset_to_pam50 <- function(mat, ens2sym, pam_symbols, collapse_fun = c("sum","mean","max")) {
  cf <- match.arg(collapse_fun) #' @note Validates collapse function.
  ensembl_ids <- strip_ensver(rownames(mat)) #' @note Strips version from row names.

  idx_map <- which(ensembl_ids %in% rownames(ens2sym)) #' @note Finds mapped Ensembl IDs.
  if (!length(idx_map)) stop("No Ensembl IDs mapped to symbols. Check inputs.") #' @note Errors if no mappings.
  sym <- ens2sym[ensembl_ids[idx_map], "hgnc_gene_symbol"] #' @note Gets symbols for mapped IDs.
  hit <- tolower(sym) %in% tolower(pam_symbols) #' @note Identifies PAM50 matches.

  keep_idx <- idx_map[hit] #' @note Selects rows matching PAM50 symbols.
  if (!length(keep_idx)) stop("No PAM50 symbols matched after mapping.") #' @note Errors if no matches.
  sub <- mat[keep_idx, , drop = FALSE] #' @note Subsets matrix to matched rows.
  sym_keep <- ens2sym[strip_ensver(rownames(sub)), "hgnc_gene_symbol"] #' @note Gets symbols for subset.

  split_idx <- split(seq_len(nrow(sub)), tolower(sym_keep)) #' @note Groups rows by symbol.
  collapsed <- lapply(split_idx, function(ix) {
    if (cf == "sum")      colSums(sub[ix, , drop = FALSE], na.rm = TRUE) #' @note Sums rows for "sum".
    else if (cf == "mean")  colMeans(sub[ix, , drop = FALSE], na.rm = TRUE) #' @note Averages rows for "mean".
    else                    apply(sub[ix, , drop = FALSE], 2, max, na.rm = TRUE) #' @note Takes max for "max".
  })
  M <- do.call(rbind, collapsed) #' @note Combines collapsed rows.
  rownames(M) <- toupper(names(split_idx)) #' @note Sets row names to uppercase symbols.
  M #' @return Subset matrix.
}





#' @section Data Loading and Filtering:
#' Loads and filters TCGA and DSMZ datasets to focus on breast cancer samples and common genes.
cat("[INFO] Starting pipeline...\n") #' @note Logs pipeline start.

tcga_se <- readRDS(tcga_se_rds) #' @note Loads TCGA SummarizedExperiment.
stopifnot("project_id" %in% colnames(colData(tcga_se))) #' @note Ensures project_id column exists.
dsmz_raw  <- readRDS(dsmz_rds) #' @note Loads DSMZ raw count data.
dsmz_meta <- read.csv(dsmz_meta_csv, check.names = FALSE) #' @note Loads DSMZ metadata.

tcga_counts <- assay(tcga_se) #' @note Extracts count matrix from TCGA SummarizedExperiment.
rownames(tcga_counts) <- strip_ensver(rownames(tcga_counts)) #' @note Strips version from Ensembl IDs.
if (any(duplicated(rownames(tcga_counts)))) {
  dup_genes <- rownames(tcga_counts)[duplicated(rownames(tcga_counts))] #' @note Identifies duplicate genes.
  cat(sprintf("[INFO] Collapsing %d duplicate TCGA genes: %s\n", 
              length(dup_genes), paste(head(dup_genes, 5), collapse=", "))) #' @note Logs duplicates.
  tcga_counts <- rowsum(tcga_counts, rownames(tcga_counts), reorder = TRUE) #' @note Collapses duplicates.
}
storage.mode(tcga_counts) <- "double" #' @note Ensures double precision for counts.

dsmz_counts <- build_dsmz_matrix(dsmz_raw) #' @note Builds DSMZ count matrix.
storage.mode(dsmz_counts) <- "double" #' @note Ensures double precision.

tcga_counts_raw <- tcga_counts #' @note Stores raw TCGA counts.
dsmz_counts_raw <- dsmz_counts #' @note Stores raw DSMZ counts.
dsmz_meta_raw <- dsmz_meta #' @note Stores raw DSMZ metadata.

report_dims("TCGA (raw)", tcga_counts_raw) #' @note Logs raw TCGA dimensions.
report_dims("DSMZ (raw)", dsmz_counts_raw) #' @note Logs raw DSMZ dimensions.

dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name) #' @note Standardizes sample ID column.
matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts)) #' @note Finds matching samples.
dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched) #' @note Filters metadata to matched samples.
dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop = FALSE] #' @note Subsets counts to matched samples.
stopifnot("organ" %in% colnames(dsmz_meta)) #' @note Ensures organ column exists.

if (!is.null(dsmz_organ_filter)) {
  keep_ids <- dsmz_meta$sample_id[dsmz_meta$organ %in% dsmz_organ_filter] #' @note Filters to breast samples.
  dsmz_meta <- dsmz_meta[dsmz_meta$sample_id %in% keep_ids, , drop = FALSE] #' @note Updates metadata.
  dsmz_counts <- dsmz_counts[, keep_ids, drop = FALSE] #' @note Updates counts.
}

report_dims("DSMZ (after organ subset)", dsmz_counts) #' @note Logs filtered DSMZ dimensions.

if (!is.null(dsmz_organ_filter) && ncol(dsmz_counts) == 0L) {
  stop(sprintf("[FATAL] DSMZ organ filter %s left 0 samples.", paste(dsmz_organ_filter, collapse=", "))) #' @note Errors if no samples remain.
}

if (!is.null(tcga_project_filter)) {
  keep_tcga <- colnames(tcga_counts) %in% colnames(assay(tcga_se)) &
               as.character(colData(tcga_se)[colnames(tcga_counts), "project_id"]) %in% tcga_project_filter #' @note Filters to TCGA-BRCA.
  tcga_counts <- tcga_counts[, keep_tcga, drop = FALSE] #' @note Updates counts.
}

report_dims("TCGA (after project subset)", tcga_counts) #' @note Logs filtered TCGA dimensions.

if (!is.null(tcga_project_filter) && ncol(tcga_counts) == 0L) {
  stop(sprintf("[FATAL] TCGA project filter %s left 0 samples.", paste(tcga_project_filter, collapse=", "))) #' @note Errors if no samples remain.
}

if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA-BRCA filter left 0 samples.") #' @note Errors if TCGA has no samples.
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ Breast filter left 0 samples.") #' @note Errors if DSMZ has no samples.

common <- intersect(rownames(tcga_counts), rownames(dsmz_counts)) #' @note Finds common genes.
cat(sprintf("[INFO] Shared genes: %d\n", length(common))) #' @note Logs number of shared genes.
if (length(common) < 1000) cat("[INFO] Few shared genes; check Ensembl IDs\n") #' @note Warns if few shared genes.

tcga_counts <- tcga_counts[common, , drop = FALSE] #' @note Subsets TCGA to common genes.
dsmz_counts <- dsmz_counts[common, , drop = FALSE] #' @note Subsets DSMZ to common genes.

report_dims("TCGA (post-intersect)", tcga_counts) #' @note Logs TCGA dimensions after intersection.
report_dims("DSMZ (post-intersect)", dsmz_counts) #' @note Logs DSMZ dimensions after intersection.

Xc_raw <- cbind(tcga_counts, dsmz_counts) #' @note Merges TCGA and DSMZ counts.
if (any(duplicated(rownames(Xc_raw)))) {
  dup_genes <- rownames(Xc_raw)[duplicated(rownames(Xc_raw))] #' @note Identifies duplicates in merged matrix.
  cat(sprintf("[INFO] Collapsing %d duplicate genes in merged matrix: %s\n", 
              length(dup_genes), paste(head(dup_genes, 5), collapse=", "))) #' @note Logs duplicates.
  Xc_raw <- rowsum(Xc_raw, rownames(Xc_raw), reorder = TRUE) #' @note Collapses duplicates.
}

report_dims("Merged (Xc_raw)", Xc_raw) #' @note Logs merged matrix dimensions.

if (ncol(tcga_counts) == 0L) stop("[FATAL] TCGA subset has 0 samples.") #' @note Errors if TCGA subset is empty.
if (ncol(dsmz_counts) == 0L) stop("[FATAL] DSMZ subset has 0 samples.") #' @note Errors if DSMZ subset is empty.
if (nrow(Xc_raw) == 0L) stop("[FATAL] No common genes after intersection.") #' @note Errors if no common genes.
if (ncol(Xc_raw) == 0L) stop("[FATAL] Merged matrix has 0 samples.") #' @note Errors if merged matrix is empty.

dim_report <- file.path(outdir, "dimension_report.txt") #' @note Path for dimension report file.
con <- file(dim_report, open = "wt") #' @note Opens file connection for writing.
on.exit(close(con), add = TRUE) #' @note Ensures file is closed on exit.

write_dims_line(con, "TCGA (raw)", tcga_counts_raw) #' @note Writes raw TCGA dimensions.
write_dims_line(con, "DSMZ (raw)", dsmz_counts_raw) #' @note Writes raw DSMZ dimensions.
write_dims_line(con, "DSMZ (after organ subset)", dsmz_counts) #' @note Writes filtered DSMZ dimensions.
write_dims_line(con, "TCGA (after project subset)", tcga_counts) #' @note Writes filtered TCGA dimensions.
writeLines(sprintf("Shared genes\t%d", length(common)), con) #' @note Writes number of shared genes.
write_dims_line(con, "TCGA (post-intersect)", tcga_counts) #' @note Writes post-intersect TCGA dimensions.
write_dims_line(con, "DSMZ (post-intersect)", dsmz_counts) #' @note Writes post-intersect DSMZ dimensions.
write_dims_line(con, "Merged (Xc_raw)", Xc_raw) #' @note Writes merged matrix dimensions.

dir.create(file.path(outdir, "shared_genes"), showWarnings = FALSE, recursive = TRUE) #' @note Creates directory for shared genes.
dir.create(file.path(outdir, "unshared_genes"), showWarnings = FALSE, recursive = TRUE) #' @note Creates directory for unshared genes.
writeLines(common, file.path(outdir, "shared_genes", "shared_genes.txt")) #' @note Saves shared genes.
writeLines(setdiff(rownames(tcga_counts_raw), common), file.path(outdir, "unshared_genes", "tcga_unshared_genes.txt")) #' @note Saves TCGA unshared genes.
writeLines(setdiff(rownames(dsmz_counts_raw), common), file.path(outdir, "unshared_genes", "dsmz_unshared_genes.txt")) #' @note Saves DSMZ unshared genes.

batch <- factor(c(rep("TCGA", ncol(tcga_counts)), rep("DSMZ", ncol(dsmz_counts)))) #' @note Creates batch factor for TCGA and DSMZ.
all_counts <- Xc_raw #' @note Stores merged counts.

if (nrow(Xc_raw) > 50000 || ncol(Xc_raw) > 5000) {
  cat("[WARN] Large merged matrix (%d genes x %d samples); consider subsampling.\n", 
      nrow(Xc_raw), ncol(Xc_raw)) #' @note Warns if merged matrix is large.
}

coldata <- tibble(
  sample = colnames(all_counts),
  dataset = ifelse(colnames(all_counts) %in% colnames(tcga_counts), "TCGA", "DSMZ")
) %>% column_to_rownames("sample") #' @note Creates colData with sample and dataset information.

if (any(duplicated(rownames(all_counts)))) {
  all_counts <- rowsum(all_counts, rownames(all_counts), reorder = TRUE) #' @note Collapses any duplicate genes.
}
dds <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ 1) #' @note Creates DESeq2 dataset.
vsd <- vst(dds, blind = TRUE) #' @note Applies VST normalization.
V <- assay(vsd) #' @note Extracts VST matrix.
V_tcga <- V[, colnames(tcga_counts), drop = FALSE] #' @note Subsets VST for TCGA samples.
V_dsmz <- V[, colnames(dsmz_counts), drop = FALSE] #' @note Subsets VST for DSMZ samples.

#' @section VST Per Dataset:
#' Applies VST normalization individually to TCGA and DSMZ datasets.
cat("[INFO] VST per-dataset (TCGA only)...\n") #' @note Logs TCGA VST start.
dds_tcga_ind <- DESeqDataSetFromMatrix(round(tcga_counts), data.frame(row.names = colnames(tcga_counts)), design = ~ 1) #' @note Creates TCGA DESeq2 dataset.
vsd_tcga_ind <- vst(dds_tcga_ind, blind = TRUE) #' @note Applies VST to TCGA.
V_tcga_ind   <- assay(vsd_tcga_ind) #' @note Extracts TCGA VST matrix.
saveRDS(V_tcga_ind, file.path(outdir, "VST_tcga_only.rds")) #' @note Saves TCGA VST matrix.

cat("[INFO] VST per-dataset (DSMZ only)...\n") #' @note Logs DSMZ VST start.
dds_dsmz_ind <- DESeqDataSetFromMatrix(round(dsmz_counts), data.frame(row.names = colnames(dsmz_counts)), design = ~ 1) #' @note Creates DSMZ DESeq2 dataset.
vsd_dsmz_ind <- vst(dds_dsmz_ind, blind = TRUE) #' @note Applies VST to DSMZ.
V_dsmz_ind   <- assay(vsd_dsmz_ind) #' @note Extracts DSMZ VST matrix.
saveRDS(V_dsmz_ind, file.path(outdir, "VST_dsmz_only.rds")) #' @note Saves DSMZ VST matrix.

purity <- load_purity_data(tcga_se, NULL) #' @note Loads purity data for TCGA.
if (!is.null(purity)) {
  cat("[INFO] Applying purity adjustment to TCGA VST data...\n") #' @note Logs purity adjustment.
  V_tcga <- purity_adjust_log(V_tcga, purity) #' @note Adjusts TCGA VST for purity.
}

cat("[INFO] Batch adjustment (global)...\n") #' @note Logs batch adjustment start.
method <- tolower(config$batch_adjust_method) #' @note Gets batch adjustment method.
V_adj <- V #' @note Initializes adjusted VST matrix.
if (method == "none") {
  cat("[INFO] Skipping batch correction.\n") #' @note Logs if no batch correction.
} else if (method == "combat_seq") {
  cat("[WARN] ComBat_seq requires raw counts; skipping batch correction for VST data.\n") #' @note Warns if ComBat_seq is selected.
} else if (method == "combat") {
  V_adj[is.na(V_adj)] <- 0 #' @note Replaces NA values with 0.
  set.seed(42) #' @note Sets seed for reproducibility.
  V_adj <- sva::ComBat(dat = as.matrix(V_adj), batch = batch, mod = NULL, 
                       par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE) #' @note Applies ComBat batch correction.
} else {
  stop("[ERROR] Unknown batch_adjust_method. Use 'none', 'combat_seq', or 'combat'.") #' @note Errors for invalid method.
}
V_tcga <- V_adj[, colnames(tcga_counts), drop = FALSE] #' @note Updates TCGA VST matrix.
V_dsmz <- V_adj[, colnames(dsmz_counts), drop = FALSE] #' @note Updates DSMZ VST matrix.

stopifnot(!is.null(colnames(V_tcga)), !is.null(colnames(V_dsmz))) #' @note Ensures column names exist.

#' @section PCA Analysis:
#' Performs PCA before and after batch correction to visualize batch effects.
cat("[INFO] PCA before and after batch correction...\n") #' @note Logs PCA start.

#' @brief Creates PCA plots for a matrix.
#' @param M Input matrix (genes as rows, samples as columns).
#' @param lab_batch Batch labels for samples.
#' @param title Plot title.
#' @param path_pdf Output PDF file path.
#' @param n_top Number of top variable genes to use (default: 5000).
#' @return List with PCA scores and explained variance.
make_pca_plot <- function(M, lab_batch, title, path_pdf, n_top = 5000) {
  sel <- top_var_genes(M, n = min(n_top, nrow(M))) #' @note Selects top variable genes.
  X   <- t(scale(M[sel, , drop = FALSE])) #' @note Scales and transposes matrix.
  pr  <- prcomp(X, center = FALSE, scale. = FALSE) #' @note Performs PCA.

  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100 #' @note Computes explained variance.
  df <- data.frame(PC1 = pr$x[,1], PC2 = pr$x[,2],
                   batch = lab_batch,
                   sample = rownames(pr$x),
                   stringsAsFactors = FALSE) #' @note Creates PCA data frame.
  safe_pdf(path_pdf, {
    layout(matrix(c(1,2), nrow = 1)) #' @note Sets plot layout.
    plot(df$PC1, df$PC2, pch = 19, cex = 0.7, col = as.factor(df$batch),
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title) #' @note Plots PCA scatter.
    legend("topright", legend = levels(as.factor(df$batch)),
           col = seq_along(levels(as.factor(df$batch))), pch = 19, cex = 0.8) #' @note Adds legend.
    barplot(var_expl[1:20], las = 2, main = "Explained variance (top 20 PCs)",
            ylab = "% variance", xlab = "PC") #' @note Plots variance barplot.
    layout(1) #' @note Resets layout.
  })
  invisible(list(scores = df, var_expl = var_expl)) #' @return PCA results.
}

pca_pre  <- make_pca_plot(V,      batch, "PCA (pre-batch correction)",  file.path(outdir, "pca_pre_batch.pdf")) #' @note PCA before batch correction.
pca_post <- make_pca_plot(V_adj,  batch, "PCA (post-batch correction)", file.path(outdir, "pca_post_batch.pdf")) #' @note PCA after batch correction.

saveRDS(V_tcga, file.path(outdir, "VST_tcga_brca.rds")) #' @note Saves adjusted TCGA VST matrix.
saveRDS(V_dsmz, file.path(outdir, "VST_dsmz_brca.rds")) #' @note Saves adjusted DSMZ VST matrix.

#' @section PAM50 Subsetting:
#' Subsets data to PAM50 genes for subtype analysis.
cat("[INFO] Subsetting to PAM50 genes...\n") #' @note Logs PAM50 subsetting start.

ens2sym_union <- bind_rows(
  tryCatch({
    cat("[INFO] Attempting gene symbol mapping via org.Hs.eg.db...\n") #' @note Logs annotation attempt.
    result <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(strip_ensver(rownames(V_adj))),
                         keytype = "ENSEMBL", columns = c("SYMBOL")) %>%
      filter(!is.na(SYMBOL) & nzchar(SYMBOL)) %>%
      distinct(ENSEMBL, .keep_all = TRUE) %>%
      transmute(ensembl = ENSEMBL, hgnc_gene_symbol = SYMBOL) #' @note Maps Ensembl IDs to symbols.
    cat(sprintf("[INFO] org.Hs.eg.db mapping successful: %d mappings\n", nrow(result))) #' @note Logs mapping success.
    result
  }, error = function(e) {
    cat(sprintf("[WARN] org.Hs.eg.db mapping failed: %s\n", e$message)) #' @note Logs mapping failure.
    NULL
  }),
  as_df(dsmz_raw) %>%
    transmute(ensembl = strip_ensver(Ensembl_ID), hgnc_gene_symbol = as.character(gene_name)) %>%
    filter(!is.na(hgnc_gene_symbol) & nzchar(hgnc_gene_symbol)) %>%
    distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE), #' @note Adds DSMZ mappings.
  tryCatch({
    cat("[INFO] Attempting gene symbol mapping via TCGA rowData...\n") #' @note Logs TCGA mapping attempt.
    rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se)) #' @note Gets TCGA rowData.
    symcol <- intersect(tolower(colnames(rd)), c("gene_name","symbol","hgnc_symbol")) #' @note Finds symbol column.
    if (length(symcol)) {
      result <- data.frame(
        ensembl = strip_ensver(rownames(rd)),
        hgnc_gene_symbol = as.character(rd[[colnames(rd)[match(symcol[1], tolower(colnames(rd)))] ]]),
        stringsAsFactors = FALSE
      ) %>% filter(!is.na(hgnc_gene_symbol) & hgnc_gene_symbol != "") %>%
        distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE) #' @note Creates TCGA mappings.
      cat(sprintf("[INFO] TCGA rowData mapping successful: %d mappings using column '%s'\n", 
                  nrow(result), symcol[1])) #' @note Logs TCGA mapping success.
      result
    } else {
      cat("[WARN] TCGA rowData mapping failed: no suitable symbol column found\n") #' @note Logs TCGA mapping failure.
      NULL
    }
  }, error = function(e) {
    cat(sprintf("[WARN] TCGA rowData mapping failed: %s\n", e$message)) #' @note Logs TCGA mapping error.
    NULL
  })
) %>% distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE) #' @note Combines and deduplicates mappings.
rownames(ens2sym_union) <- ens2sym_union$ensembl #' @note Sets row names to Ensembl IDs.

cat(sprintf("[INFO] Final gene symbol mapping summary: %d total mappings\n", nrow(ens2sym_union))) #' @note Logs total mappings.

M_pam50 <- subset_to_pam50(V_adj, ens2sym_union, pam50_symbols, collapse_fun = "mean") #' @note Subsets to PAM50 genes.
report_dims("PAM50 Matrix", M_pam50) #' @note Logs PAM50 matrix dimensions.
saveRDS(M_pam50, file.path(outdir, "PAM50_matrix_VST_adj.rds")) #' @note Saves PAM50 matrix.
missing_syms <- setdiff(toupper(pam50_symbols), rownames(M_pam50)) #' @note Identifies missing PAM50 symbols.
cat(sprintf("[PAM50] Matched %d/%d genes\n", nrow(M_pam50), length(pam50_symbols))) #' @note Logs matched genes.
if (length(missing_syms)) {
  cat("[PAM50] Missing symbols:", paste(missing_syms, collapse=", "), "\n") #' @note Logs missing symbols.
  writeLines(missing_syms, file.path(outdir, "PAM50_missing_symbols.txt")) #' @note Saves missing symbols.
}

marker_ens <- ens2sym_union$ensembl[tolower(ens2sym_union$hgnc_gene_symbol) %in% tolower(marker_symbols)] #' @note Maps marker symbols to Ensembl IDs.
marker_ens <- intersect(unique(marker_ens), rownames(V_adj)) #' @note Intersects with V_adj genes.
if (length(marker_ens) == 0L) {
  stop("[FATAL] No marker genes found in VST matrix. Check gene ID mapping.") #' @note Errors if no markers found.
} else {
  cat(sprintf("[INFO] Found %d/%d marker genes in VST matrix.\n", 
              length(marker_ens), length(marker_symbols))) #' @note Logs found markers.
  missing_markers <- setdiff(tolower(marker_symbols), tolower(ens2sym_union$hgnc_gene_symbol[ens2sym_union$ensembl %in% marker_ens])) #' @note Identifies missing markers.
  if (length(missing_markers) > 0) {
    cat(sprintf("[WARN] Missing %d markers: %s\n", 
                length(missing_markers), paste(missing_markers, collapse=", "))) #' @note Logs missing markers.
  }
}

#' @section Unsupervised Discovery Clustering (K-Means):
#' Performs k-means clustering on highly variable genes.
set.seed(42) #' @note Sets seed for reproducibility.
hvgs <- top_var_genes(V_adj, n = 3000) #' @note Selects top 3000 variable genes.
X_discovery <- t(V_adj[hvgs, , drop = FALSE]) #' @note Transposes subset matrix.
if (nrow(X_discovery) > 5000) {
  cat("[WARN] Large sample set (%d) for k-means; consider subsampling.\n", nrow(X_discovery)) #' @note Warns if large sample set.
}
X_scaled <- scale(X_discovery) #' @note Scales matrix for clustering.
writeLines(hvgs, file.path(outdir, "HVGs_used_for_discovery_3k.txt")) #' @note Saves highly variable genes.

k_grid <- 2:8 #' @note Defines range of k values to test.
res_list <- list() #' @note Stores k-means results.
metrics <- lapply(k_grid, function(k) {
  km <- kmeans(X_scaled, centers = k, nstart = 50, iter.max = 100) #' @note Performs k-means clustering.
  wcss <- sum(km$withinss) #' @note Computes within-cluster sum of squares.
  sil <- silhouette(km$cluster, dist(X_scaled)) #' @note Computes silhouette scores.
  sil_mean <- mean(sil[, "sil_width"]) #' @note Computes mean silhouette score.
  ch <- tryCatch(calinhara(X_scaled, km$cluster), error = function(e) NA_real_) #' @note Computes Calinski-Harabasz index.
  db <- tryCatch({
    m <- intCriteria(as.matrix(X_scaled), as.integer(km$cluster), "davies_bouldin") #' @note Computes Davies-Bouldin index.
    as.numeric(m$davies_bouldin)
  }, error = function(e) NA_real_)
  res_list[[as.character(k)]] <<- list(km = km, sil = sil) #' @note Stores results.
  data.frame(k = k, WCSS = wcss, Silhouette = sil_mean, CH = ch, DB = db) #' @return Metrics data frame.
})
metrics <- do.call(rbind, metrics) #' @note Combines metrics.

write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), row.names = FALSE) #' @note Saves k-means metrics.

safe_pdf(file.path(outdir, "kmeans_elbow_HVG.pdf"), {
  plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
       ylab = "Within-Cluster SS", main = "Elbow (HVG space)") #' @note Plots elbow curve.
})
safe_pdf(file.path(outdir, "kmeans_metrics_panel_HVG.pdf"), {
  par(mfrow = c(2,2)) #' @note Sets 2x2 plot layout.
  plot(metrics$k, metrics$Silhouette, type = "b", xlab = "k", 
       ylab = "Avg silhouette", main = "Silhouette (higher is better)") #' @note Plots silhouette scores.
  plot(metrics$k, metrics$CH, type = "b", xlab = "k", 
       ylab = "Calinski-Harabasz", main = "Calinski-Harabasz (higher)") #' @note Plots CH index.
  plot(metrics$k, metrics$DB, type = "b", xlab = "k", 
       ylab = "Davies-Bouldin", main = "Davies-Bouldin (lower)") #' @note Plots DB index.
  plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
       ylab = "WCSS", main = "Elbow (lower)") #' @note Plots WCSS.
  par(mfrow = c(1,1)) #' @note Resets layout.
})

best_k <- with(metrics, {
  cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)] #' @note Finds k with max silhouette.
  if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand #' @note Resolves ties with CH index.
})
message(sprintf("[kmeans] best k by silhouette→CH = %s", best_k)) #' @note Logs optimal k.

km_final <- res_list[[as.character(best_k)]]$km #' @note Gets final k-means model.
sil_final <- res_list[[as.character(best_k)]]$sil #' @note Gets final silhouette scores.
clusters_kmeans <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k))) #' @note Creates cluster labels.

sample_is_tcga <- colnames(V_adj) %in% colnames(V_tcga) #' @note Identifies TCGA samples.
dataset_lab <- ifelse(sample_is_tcga, "TCGA", "DSMZ") #' @note Labels samples by dataset.

cluster_df <- tibble(sample = rownames(X_scaled),
                     dataset = dataset_lab,
                     cluster = as.character(clusters_kmeans)) #' @note Creates cluster data frame.
write.csv(cluster_df, file.path(outdir, "clusters_kmeans.csv"), row.names = FALSE) #' @note Saves k-means clusters.

safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_HVG.pdf", best_k)), {
  plot(sil_final, main = sprintf("Silhouette (k = %d, HVG space)", best_k), cex.names = 0.8) #' @note Plots silhouette for best k.
})

sil_dir <- file.path(outdir, "silhouettes_per_k") #' @note Directory for silhouette plots.
dir.create(sil_dir, showWarnings = FALSE, recursive = TRUE) #' @note Creates directory.
for (k in k_grid) {
  silk <- res_list[[as.character(k)]]$sil #' @note Gets silhouette for k.
  if (!is.null(silk)) {
    pdf(file.path(sil_dir, sprintf("silhouette_k=%d.pdf", k)), width = 8, height = 6) #' @note Opens PDF for silhouette.
    plot(silk, main = sprintf("Silhouette (k = %d, HVG space)", k), cex.names = 0.7) #' @note Plots silhouette.
    dev.off() #' @note Closes PDF.
  }
}

cat("Optimal k (k-means):", best_k, "with mean silhouette:", round(mean(sil_final[, "sil_width"]), 3), "\n") #' @note Logs optimal k and silhouette.
print(table(clusters_kmeans, dataset_lab)) #' @note Prints cluster distribution.

#' @section UMAP for K-Means:
#' Visualizes k-means clusters using UMAP.
set.seed(42) #' @note Sets seed for reproducibility.
emb <- uwot::umap(X_discovery, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs UMAP.
emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2") #' @note Converts to data frame.
emb$sample <- rownames(X_discovery) #' @note Adds sample names.
emb$dataset <- dataset_lab #' @note Adds dataset labels.
emb$cluster <- clusters_kmeans #' @note Adds cluster labels.

p_umap <- ggplot(emb, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (k-means, HVG space)") #' @note Creates UMAP plot.

ggsave(file.path(outdir, "umap_kmeans_HVG.pdf"), p_umap, width = 7, height = 6) #' @note Saves UMAP plot.
write.csv(emb, file.path(outdir, "umap_kmeans_HVG_coords.csv"), row.names = FALSE) #' @note Saves UMAP coordinates.

#' @section Hierarchical Clustering (HVG):
#' Performs hierarchical clustering on highly variable genes.
cat("[INFO] Performing hierarchical clustering (HVG)...\n") #' @note Logs hierarchical clustering start.
if (nrow(X_scaled) < 5) {
  stop(sprintf("[FATAL] Too few samples (%d) for hierarchical clustering (need at least 5).", nrow(X_scaled))) #' @note Errors if too few samples.
}

cat("[INFO] Using Euclidean distance for hierarchical clustering (standard for gene expression)\n") #' @note Logs distance metric choice.
d <- dist(X_scaled, method = "euclidean") #' @note Computes Euclidean distance.
hc <- hclust(d, method = "average") #' @note Performs hierarchical clustering with average linkage.
clusters_hc <- cutree(hc, k = 5) #' @note Cuts tree into 5 clusters.
subtype_labels <- c("HER2-high", "HER2-low", "LumA", "LumB", "Basal") #' @note Defines subtype labels.

tcga_sub_col <- c("PAM50","Subtype","BRCA_Subtype","molecular_subtype") #' @note Possible subtype columns.
tcga_sub_col <- tcga_sub_col[tcga_sub_col %in% colnames(colData(tcga_se))][1] #' @note Selects first valid column.
if (!is.na(tcga_sub_col)) {
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Gets TCGA subtype labels.
  sub_lvls <- sort(unique(na.omit(lab_tcga))) #' @note Gets unique subtype levels.
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE)) #' @note Computes subtype means.
  cl_means <- sapply(1:5, function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes cluster means.
  g_common <- intersect(rownames(sub_means), hvgs) #' @note Finds common genes.
  corr <- cor(cl_means[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman") #' @note Correlates clusters with subtypes.
  cluster_to_subtype <- apply(corr, 1, function(v) sub_lvls[which.max(v)]) #' @note Maps clusters to subtypes.
  clusters_hc <- factor(clusters_hc, labels = subtype_labels[match(cluster_to_subtype, sub_lvls)]) #' @note Labels clusters.
} else {
  clusters_hc <- factor(clusters_hc, labels = subtype_labels) #' @note Uses default labels if no TCGA subtypes.
}

cluster_hc_df <- tibble(sample = rownames(X_scaled),
                        dataset = dataset_lab,
                        cluster = as.character(clusters_hc)) #' @note Creates hierarchical cluster data frame.
write.csv(cluster_hc_df, file.path(outdir, "clusters_hierarchical_HVG.csv"), row.names = FALSE) #' @note Saves clusters.

safe_pdf(file.path(outdir, "dendrogram_hierarchical_HVG.pdf"), {
  plot(hc, main = "Hierarchical Clustering (HVG, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5) #' @note Plots dendrogram.
  abline(h = hc$height[length(hc$height) - 4], col = "red", lty = 2) #' @note Adds cut line.
})

#' @section UMAP for Hierarchical Clustering (HVG):
#' Visualizes hierarchical clusters using UMAP.
set.seed(42) #' @note Sets seed for reproducibility.
emb_hc <- uwot::umap(X_scaled, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs UMAP.
emb_hc <- as.data.frame(emb_hc); colnames(emb_hc) <- c("UMAP1","UMAP2") #' @note Converts to data frame.
emb_hc$sample <- rownames(X_scaled) #' @note Adds sample names.
emb_hc$dataset <- dataset_lab #' @note Adds dataset labels.
emb_hc$cluster <- clusters_hc #' @note Adds cluster labels.

p_hc_umap <- ggplot(emb_hc, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (Hierarchical Clustering, HVG space)") #' @note Creates UMAP plot.

ggsave(file.path(outdir, "umap_hierarchical_HVG.pdf"), p_hc_umap, width = 7, height = 6) #' @note Saves UMAP plot.
write.csv(emb_hc, file.path(outdir, "umap_hierarchical_HVG_coords.csv"), row.names = FALSE) #' @note Saves UMAP coordinates.

#' @section Hierarchical Clustering (PAM50):
#' Performs hierarchical clustering on PAM50 genes.
cat("[INFO] Performing hierarchical clustering (PAM50)...\n") #' @note Logs PAM50 clustering start.
X_pam50 <- t(scale(M_pam50)) #' @note Scales and transposes PAM50 matrix.
if (nrow(X_pam50) < 5) {
  stop(sprintf("[FATAL] Too few samples (%d) for PAM50 hierarchical clustering (need at least 5).", nrow(X_pam50))) #' @note Errors if too few samples.
}

cat("[INFO] Using Euclidean distance for PAM50 hierarchical clustering (consistent with HVG analysis)\n") #' @note Logs distance metric.
d_pam50 <- dist(X_pam50, method = "euclidean") #' @note Computes Euclidean distance.
hc_pam50 <- hclust(d_pam50, method = "average") #' @note Performs hierarchical clustering.
clusters_hc_pam50 <- cutree(hc_pam50, k = 5) #' @note Cuts tree into 5 clusters.

if (!is.na(tcga_sub_col)) {
  cl_means_pam50 <- sapply(1:5, function(cn) rowMeans(M_pam50[, clusters_hc_pam50 == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes cluster means.
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50)) #' @note Finds common genes.
  corr_pam50 <- cor(cl_means_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman") #' @note Correlates clusters with subtypes.
  cluster_to_subtype_pam50 <- apply(corr_pam50, 1, function(v) sub_lvls[which.max(v)]) #' @note Maps clusters to subtypes.
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = subtype_labels[match(cluster_to_subtype_pam50, sub_lvls)]) #' @note Labels clusters.
} else {
  clusters_hc_pam50 <- factor(clusters_hc_pam50, labels = subtype_labels) #' @note Uses default labels.
}

cluster_hc_pam50_df <- tibble(sample = rownames(X_pam50),
                              dataset = dataset_lab,
                              cluster = as.character(clusters_hc_pam50)) #' @note Creates PAM50 cluster data frame.
write.csv(cluster_hc_pam50_df, file.path(outdir, "clusters_hierarchical_PAM50.csv"), row.names = FALSE) #' @note Saves clusters.

safe_pdf(file.path(outdir, "dendrogram_hierarchical_PAM50.pdf"), {
  plot(hc_pam50, main = "Hierarchical Clustering (PAM50, Euclidean, Average Linkage)", xlab = "", sub = "", cex = 0.5) #' @note Plots dendrogram.
  abline(h = hc_pam50$height[length(hc_pam50$height) - 4], col = "red", lty = 2) #' @note Adds cut line.
})

#' @section UMAP for Hierarchical Clustering (PAM50):
#' Visualizes PAM50 hierarchical clusters using UMAP.
set.seed(42) #' @note Sets seed for reproducibility.
emb_hc_pam50 <- uwot::umap(X_pam50, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs UMAP.
emb_hc_pam50 <- as.data.frame(emb_hc_pam50); colnames(emb_hc_pam50) <- c("UMAP1","UMAP2") #' @note Converts to data frame.
emb_hc_pam50$sample <- rownames(X_pam50) #' @note Adds sample names.
emb_hc_pam50$dataset <- dataset_lab #' @note Adds dataset labels.
emb_hc_pam50$cluster <- clusters_hc_pam50 #' @note Adds cluster labels.

p_hc_pam50_umap <- ggplot(emb_hc_pam50, aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (Hierarchical Clustering, PAM50 space)") #' @note Creates UMAP plot.

ggsave(file.path(outdir, "umap_hierarchical_PAM50.pdf"), p_hc_pam50_umap, width = 7, height = 6) #' @note Saves UMAP plot.
write.csv(emb_hc_pam50, file.path(outdir, "umap_hierarchical_PAM50_coords.csv"), row.names = FALSE) #' @note Saves UMAP coordinates.

#' @section Differential Expression Analysis:
#' Performs differential expression analysis using DESeq2.
cat("[INFO] Performing differential expression analysis...\n") #' @note Logs DE analysis start.

subtype_levels <- c("Basal","HER2-high","HER2-low","LumA","LumB") #' @note Defines subtype levels.
cl_here <- factor(clusters_hc, levels = subtype_levels) #' @note Sets hierarchical clusters as subtypes.
cl_here <- droplevels(cl_here) #' @note Drops unused levels.
coldata$subtype <- cl_here #' @note Adds subtypes to colData.

levels(coldata$subtype) <- gsub("[^A-Za-z0-9_.]", "_", levels(coldata$subtype)) #' @note Sanitizes factor levels.

dds_de <- DESeqDataSetFromMatrix(round(all_counts), coldata, design = ~ subtype) #' @note Creates DESeq2 dataset for DE.
dds_de <- DESeq(dds_de) #' @note Runs DESeq2 analysis.

present_lvls <- levels(coldata$subtype) #' @note Gets present subtype levels.
targets <- setdiff(present_lvls, "Basal") #' @note Defines non-Basal subtypes for comparison.
if (!length(targets)) {
  stop("[DE] No non-Basal subtypes present to contrast against Basal.") #' @note Errors if no non-Basal subtypes.
}

de_results <- list() #' @note Stores DE results.
for (sub in targets) {
  res <- results(dds_de, contrast = c("subtype", sub, "Basal"), alpha = 0.05) #' @note Computes DE results vs. Basal.
  res <- as.data.frame(res) %>%
    mutate(ensembl = rownames(res),
           hgnc_gene_symbol  = ens2sym_union$hgnc_gene_symbol[match(ensembl, ens2sym_union$ensembl)]) %>%
    filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange), abs(log2FoldChange) > 1) %>%
    arrange(padj) #' @note Filters significant genes and adds symbols.
  de_results[[sub]] <- res #' @note Stores results for subtype.
}
dir.create(file.path(outdir, "de_results"), showWarnings = FALSE, recursive = TRUE) #' @note Creates DE results directory.
for (sub in names(de_results)) {
  write.csv(de_results[[sub]], file.path(outdir, "de_results", sprintf("de_%s_vs_Basal.csv", sub)), row.names = FALSE) #' @note Saves DE results.
}

#' @section DE Heatmap:
#' Creates heatmap of top differentially expressed genes.
top_de_genes <- unique(unlist(lapply(de_results, function(res) head(res$ensembl, 10)))) #' @note Selects top 10 DE genes per contrast.
if (length(top_de_genes) > 0) {
  V_de <- V_adj[top_de_genes, , drop = FALSE] #' @note Subsets VST matrix to DE genes.
  if (any(duplicated(rownames(V_de)))) {
    rownames(V_de) <- make.unique(rownames(V_de)) #' @note Ensures unique row names.
  }
  lab <- ens2sym_union$hgnc_gene_symbol[match(rownames(V_de), ens2sym_union$ensembl)] #' @note Maps Ensembl IDs to symbols.
  lab[is.na(lab) | lab == ""] <- rownames(V_de) #' @note Uses Ensembl IDs for missing symbols.
  rownames(V_de) <- lab #' @note Sets row names to symbols.
  ann <- data.frame(dataset = dataset_lab, subtype = clusters_hc, row.names = colnames(V_de)) #' @note Creates annotation data frame.
  safe_pdf(file.path(outdir, "heatmap_de_genes.pdf"), {
    pheatmap(V_de[, order(clusters_hc), drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_col = ann,
             main = "Top DE Genes by Subtype") #' @note Plots DE heatmap.
  })
}

#' @section UMAP + HDBSCAN (HVG):
#' Performs iterative UMAP and HDBSCAN clustering on highly variable genes.
cat("[INFO] Performing iterative UMAP + HDBSCAN (HVG)...\n") #' @note Logs HDBSCAN start.
set.seed(42) #' @note Sets seed for reproducibility.
if (nrow(X_scaled) < 20) {
  stop(sprintf("[FATAL] Too few samples (%d) for UMAP+HDBSCAN (need at least 20).", nrow(X_scaled))) #' @note Errors if too few samples.
}
if (nrow(X_scaled) > 5000) {
  cat("[WARN] Large sample set (%d) for UMAP+HDBSCAN; consider subsampling.\n", nrow(X_scaled)) #' @note Warns if large sample set.
}

emb1 <- uwot::umap(X_scaled, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs initial UMAP.
hdb1 <- hdbscan(emb1, minPts = 10) #' @note Performs initial HDBSCAN.
clusters_hdb <- rep("Unassigned", nrow(X_scaled)) #' @note Initializes cluster assignments.
basal_idx <- which(hdb1$cluster == which.max(table(hdb1$cluster[hdb1$cluster != 0]))) #' @note Identifies Basal cluster.
clusters_hdb[basal_idx] <- "Basal" #' @note Assigns Basal labels.

non_basal <- which(clusters_hdb == "Unassigned") #' @note Identifies non-Basal samples.
if (length(non_basal) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification.\n", length(non_basal)) #' @note Warns if too few samples.
} else {
  emb2 <- uwot::umap(X_scaled[non_basal, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine") #' @note UMAP for non-Basal samples.
  hdb2 <- hdbscan(emb2, minPts = 5) #' @note HDBSCAN for non-Basal samples.
  her2_high_idx <- non_basal[which(hdb2$cluster == which.max(table(hdb2$cluster[hdb2$cluster != 0])))] #' @note Identifies HER2-high cluster.
  clusters_hdb[her2_high_idx] <- "HER2-high" #' @note Assigns HER2-high labels.
}

non_her2_high <- which(clusters_hdb == "Unassigned") #' @note Identifies non-HER2-high samples.
if (length(non_her2_high) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification.\n", length(non_her2_high)) #' @note Warns if too few samples.
} else {
  emb3 <- uwot::umap(X_scaled[non_her2_high, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine") #' @note UMAP for non-HER2-high samples.
  hdb3 <- hdbscan(emb3, minPts = 5) #' @note HDBSCAN for non-HER2-high samples.
  her2_low_idx <- non_her2_high[which(hdb3$cluster == which.max(table(hdb3$cluster[hdb3$cluster != 0])))] #' @note Identifies HER2-low cluster.
  clusters_hdb[her2_low_idx] <- "HER2-low" #' @note Assigns HER2-low labels.
}

luminal_idx <- which(clusters_hdb == "Unassigned") #' @note Identifies remaining samples.
if (length(luminal_idx) < 2) {
  cat("[WARN] Too few Luminal samples (%d); assigning all as LumA.\n", length(luminal_idx)) #' @note Warns if too few samples.
  clusters_hdb[luminal_idx] <- "LumA" #' @note Assigns LumA labels.
} else {
  emb_lum <- uwot::umap(X_scaled[luminal_idx, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "cosine") #' @note UMAP for Luminal samples.
  km_lum <- kmeans(emb_lum, centers = 2, nstart = 50) #' @note K-means for Luminal samples.
  lum_clusters <- ifelse(km_lum$cluster == 1, "LumA", "LumB") #' @note Assigns LumA or LumB labels.
  clusters_hdb[luminal_idx] <- lum_clusters #' @note Updates cluster assignments.
}

clusters_hdb <- factor(clusters_hdb, levels = subtype_labels) #' @note Converts to factor with subtype levels.
hdb_df <- tibble(sample = rownames(X_scaled),
                 dataset = dataset_lab,
                 subtype = clusters_hdb) #' @note Creates HDBSCAN data frame.
write.csv(hdb_df, file.path(outdir, "clusters_hdbscan_HVG.csv"), row.names = FALSE) #' @note Saves HDBSCAN clusters.

emb_hdb <- as.data.frame(emb1); colnames(emb_hdb) <- c("UMAP1","UMAP2") #' @note Converts initial UMAP to data frame.
emb_hdb$sample <- rownames(X_scaled) #' @note Adds sample names.
emb_hdb$dataset <- dataset_lab #' @note Adds dataset labels.
emb_hdb$subtype <- clusters_hdb #' @note Adds subtype labels.

p_hdb_umap <- ggplot(emb_hdb, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + theme_bw() +
  ggtitle("UMAP (HDBSCAN Subtypes, HVG space)") #' @note Creates UMAP plot.

ggsave(file.path(outdir, "umap_hdbscan_HVG.pdf"), p_hdb_umap, width = 7, height = 6) #' @note Saves UMAP plot.
write.csv(emb_hdb, file.path(outdir, "umap_hdbscan_HVG_coords.csv"), row.names = FALSE) #' @note Saves UMAP coordinates.

#' @section UMAP + HDBSCAN (PAM50):
#' Performs iterative UMAP and HDBSCAN clustering on PAM50 genes.
cat("[INFO] Performing iterative UMAP + HDBSCAN (PAM50)...\n") #' @note Logs PAM50 HDBSCAN start.
set.seed(42) #' @note Sets seed for reproducibility.
if (nrow(X_pam50) < 20) {
  stop(sprintf("[FATAL] Too few samples (%d) for PAM50 UMAP+HDBSCAN (need at least 20).", nrow(X_pam50))) #' @note Errors if too few samples.
}
if (nrow(X_pam50) > 5000) {
  cat("[WARN] Large sample set (%d) for PAM50 UMAP+HDBSCAN; consider subsampling.\n", nrow(X_pam50)) #' @note Warns if large sample set.
}

emb1_pam50 <- uwot::umap(X_pam50, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs initial UMAP.
hdb1_pam50 <- hdbscan(emb1_pam50, minPts = 10) #' @note Performs initial HDBSCAN.
clusters_hdb_pam50 <- rep("Unassigned", nrow(X_pam50)) #' @note Initializes cluster assignments.
basal_idx_pam50 <- which(hdb1_pam50$cluster == which.max(table(hdb1_pam50$cluster[hdb1_pam50$cluster != 0]))) #' @note Identifies Basal cluster.
clusters_hdb_pam50[basal_idx_pam50] <- "Basal" #' @note Assigns Basal labels.

non_basal_pam50 <- which(clusters_hdb_pam50 == "Unassigned") #' @note Identifies non-Basal samples.
if (length(non_basal_pam50) < 10) {
  cat("[WARN] Too few non-Basal samples (%d); skipping HER2-high identification (PAM50).\n", length(non_basal_pam50)) #' @note Warns if too few samples.
} else {
  emb2_pam50 <- uwot::umap(X_pam50[non_basal_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine") #' @note UMAP for non-Basal samples.
  hdb2_pam50 <- hdbscan(emb2_pam50, minPts = 5) #' @note HDBSCAN for non-Basal samples.
  her2_high_idx_pam50 <- non_basal_pam50[which(hdb2_pam50$cluster == which.max(table(hdb2_pam50$cluster[hdb2_pam50$cluster != 0])))] #' @note Identifies HER2-high cluster.
  clusters_hdb_pam50[her2_high_idx_pam50] <- "HER2-high" #' @note Assigns HER2-high labels.
}

non_her2_high_pam50 <- which(clusters_hdb_pam50 == "Unassigned") #' @note Identifies non-HER2-high samples.
if (length(non_her2_high_pam50) < 10) {
  cat("[WARN] Too few non-HER2-high samples (%d); skipping HER2-low identification (PAM50).\n", length(non_her2_high_pam50)) #' @note Warns if too few samples.
} else {
  emb3_pam50 <- uwot::umap(X_pam50[non_her2_high_pam50, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "cosine") #' @note UMAP for non-HER2-high samples.
  hdb3_pam50 <- hdbscan(emb3_pam50, minPts = 5) #' @note HDBSCAN for non-HER2-high samples.
  her2_low_idx_pam50 <- non_her2_high_pam50[which(hdb3_pam50$cluster == which.max(table(hdb3_pam50$cluster[hdb3_pam50$cluster != 0])))] #' @note Identifies HER2-low cluster.
  clusters_hdb_pam50[her2_low_idx_pam50] <- "HER2-low" #' @note Assigns HER2-low labels.
}

#' @file brca_vst_dual_analysis_continued.R
#' @brief Continues integrative analysis of TCGA and DSMZ breast cancer datasets with UMAP, HDBSCAN, statistical tests, and marker gene heatmaps.
#' @details This script completes the UMAP and HDBSCAN clustering on PAM50 genes, performs statistical tests (Fisher’s Exact, Kruskal-Wallis, DSCF) to compare clusters and subtypes, and generates heatmaps and UMAP visualizations for marker genes and joint TCGA/DSMZ data. Results are saved as CSVs and PDFs.
#' @author [Your Name]
#' @date 2025-10-23

#' @section UMAP + HDBSCAN (PAM50) Continued:
#' Finalizes iterative UMAP and HDBSCAN clustering on PAM50 genes by assigning Luminal subtypes and visualizing results.
if (length(luminal_idx_pam50) < 2) {
  cat("[WARN] Too few Luminal samples (%d); assigning all as LumA (PAM50).\n", length(luminal_idx_pam50)) #' @note Warns if fewer than 2 Luminal samples are available.
  clusters_hdb_pam50[luminal_idx_pam50] <- "LumA" #' @note Assigns all Luminal samples to LumA subtype.
} else {
  emb_lum_pam50 <- uwot::umap(X_pam50[luminal_idx_pam50, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "cosine") #' @note Performs UMAP dimensionality reduction on Luminal samples with 10 neighbors and cosine metric.
  km_lum_pam50 <- kmeans(emb_lum_pam50, centers = 2, nstart = 50) #' @note Applies k-means clustering to split Luminal samples into two clusters (LumA and LumB) with 50 random starts.
  lum_clusters_pam50 <- ifelse(km_lum_pam50$cluster == 1, "LumA", "LumB") #' @note Assigns LumA or LumB labels based on k-means cluster assignments.
  clusters_hdb_pam50[luminal_idx_pam50] <- lum_clusters_pam50 #' @note Updates cluster assignments for Luminal samples.
}

clusters_hdb_pam50 <- factor(clusters_hdb_pam50, levels = subtype_labels) #' @note Converts HDBSCAN PAM50 clusters to a factor with predefined subtype levels (e.g., Basal, HER2-high, HER2-low, LumA, LumB).
hdb_pam50_df <- tibble(
  sample = rownames(X_pam50), #' @note Extracts sample names from PAM50 matrix.
  dataset = dataset_lab, #' @note Uses dataset labels (TCGA or DSMZ).
  subtype = clusters_hdb_pam50 #' @note Includes HDBSCAN PAM50 subtype assignments.
) #' @note Creates a data frame with sample, dataset, and subtype information.
write.csv(hdb_pam50_df, file.path(outdir, "clusters_hdbscan_PAM50.csv"), row.names = FALSE) #' @note Saves HDBSCAN PAM50 cluster assignments to a CSV file.

emb_hdb_pam50 <- as.data.frame(emb1_pam50); colnames(emb_hdb_pam50) <- c("UMAP1","UMAP2") #' @note Converts initial UMAP coordinates to a data frame and names columns UMAP1 and UMAP2.
emb_hdb_pam50$sample <- rownames(X_pam50) #' @note Adds sample names to the UMAP data frame.
emb_hdb_pam50$dataset <- dataset_lab #' @note Adds dataset labels (TCGA or DSMZ).
emb_hdb_pam50$subtype <- clusters_hdb_pam50 #' @note Adds HDBSCAN PAM50 subtype labels.

p_hdb_pam50_umap <- ggplot(emb_hdb_pam50, aes(UMAP1, UMAP2, color = subtype, shape = dataset)) +
  geom_point(alpha = 0.9, size = 1.8) + #' @note Creates a scatter plot of UMAP coordinates with points colored by subtype and shaped by dataset.
  theme_bw() + #' @note Applies a clean black-and-white theme.
  ggtitle("UMAP (HDBSCAN Subtypes, PAM50 space)") #' @note Sets plot title.

ggsave(file.path(outdir, "umap_hdbscan_PAM50.pdf"), p_hdb_pam50_umap, width = 7, height = 6) #' @note Saves UMAP plot to a PDF file.
write.csv(emb_hdb_pam50, file.path(outdir, "umap_hdbscan_PAM50_coords.csv"), row.names = FALSE) #' @note Saves UMAP coordinates to a CSV file.

#' @section Statistical Analysis:
#' Performs Fisher’s Exact Tests to compare cluster assignments with dataset and TCGA PAM50 subtypes, followed by Kruskal-Wallis and DSCF tests for marker gene expression.
cat("[INFO] Performing statistical analysis...\n") #' @note Logs the start of statistical analysis.
dir.create(file.path(outdir, "stats_results"), showWarnings = FALSE) #' @note Creates a directory for statistical results, ignoring warnings if it exists.

#' @subsection Fisher’s Exact Tests:
#' Tests association between clustering methods (k-means, hierarchical, HDBSCAN, HDBSCAN PAM50) and dataset (TCGA vs. DSMZ) or TCGA PAM50 subtypes.
fisher_results <- list() #' @note Initializes a list to store Fisher’s Exact Test results.
cont_table_kmeans <- table(clusters_kmeans, dataset_lab) #' @note Creates contingency table for k-means clusters vs. dataset.
cont_table_hc <- table(clusters_hc, dataset_lab) #' @note Creates contingency table for hierarchical clusters vs. dataset.
cont_table_hdb <- table(clusters_hdb, dataset_lab) #' @note Creates contingency table for HDBSCAN clusters vs. dataset.
cont_table_hdb_pam50 <- table(clusters_hdb_pam50, dataset_lab) #' @note Creates contingency table for HDBSCAN PAM50 clusters vs. dataset.

if (all(dim(cont_table_kmeans) >= 2) && sum(cont_table_kmeans) >= 5) {
  fisher_kmeans <- fisher.test(cont_table_kmeans, simulate.p.value = nrow(cont_table_kmeans) > 2) #' @note Performs Fisher’s Exact Test for k-means vs. dataset, using simulation if table is larger than 2x2.
  fisher_results[["kmeans_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_kmeans_vs_dataset",
    p_value = fisher_kmeans$p.value,
    statistic = NA_real_
  ) #' @note Stores k-means test results in a data frame.
} else {
  cat("[WARN] Insufficient data for Fisher’s test (k-means vs. dataset).\n") #' @note Warns if contingency table is too small or sparse.
}

if (all(dim(cont_table_hc) >= 2) && sum(cont_table_hc) >= 5) {
  fisher_hc <- fisher.test(cont_table_hc, simulate.p.value = nrow(cont_table_hc) > 2) #' @note Performs Fisher’s Exact Test for hierarchical clusters vs. dataset.
  fisher_results[["hc_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hc_vs_dataset",
    p_value = fisher_hc$p.value,
    statistic = NA_real_
  ) #' @note Stores hierarchical test results.
} else {
  cat("[WARN] Insufficient data for Fisher’s test (hierarchical vs. dataset).\n") #' @note Warns if contingency table is too small or sparse.
}

if (all(dim(cont_table_hdb) >= 2) && sum(cont_table_hdb) >= 5) {
  fisher_hdb <- fisher.test(cont_table_hdb, simulate.p.value = nrow(cont_table_hdb) > 2) #' @note Performs Fisher’s Exact Test for HDBSCAN clusters vs. dataset.
  fisher_results[["hdb_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_vs_dataset",
    p_value = fisher_hdb$p.value,
    statistic = NA_real_
  ) #' @note Stores HDBSCAN test results.
} else {
  cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN vs. dataset).\n") #' @note Warns if contingency table is too small or sparse.
}

if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
  fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2) #' @note Performs Fisher’s Exact Test for HDBSCAN PAM50 clusters vs. dataset.
  fisher_results[["hdb_pam50_vs_dataset"]] <- data.frame(
    test = "Fisher_Exact_hdb_pam50_vs_dataset",
    p_value = fisher_hdb_pam50$p.value,
    statistic = NA_real_
  ) #' @note Stores HDBSCAN PAM50 test results.
} else {
  cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN PAM50 vs. dataset).\n") #' @note Warns if contingency table is too small or sparse.
}

#' @subsection Fisher’s Exact Tests vs. TCGA PAM50 Subtypes:
#' Tests association between clustering methods and TCGA PAM50 subtypes for TCGA samples.
if (!is.na(tcga_sub_col)) {
  tcga_samples <- colnames(V_tcga) #' @note Extracts TCGA sample names.
  lab_tcga <- as.character(colData(tcga_se)[tcga_samples, tcga_sub_col]) #' @note Gets TCGA PAM50 subtype labels.
  cont_table_kmeans_pam50 <- table(clusters_kmeans[tcga_samples], lab_tcga) #' @note Creates contingency table for k-means vs. TCGA subtypes.
  cont_table_hc_pam50 <- table(clusters_hc[tcga_samples], lab_tcga) #' @note Creates contingency table for hierarchical vs. TCGA subtypes.
  cont_table_hdb_pam50 <- table(clusters_hdb[tcga_samples], lab_tcga) #' @note Creates contingency table for HDBSCAN vs. TCGA subtypes.
  cont_table_hdb_pam50_sub <- table(clusters_hdb_pam50[tcga_samples], lab_tcga) #' @note Creates contingency table for HDBSCAN PAM50 vs. TCGA subtypes.

  if (all(dim(cont_table_kmeans_pam50) >= 2) && sum(cont_table_kmeans_pam50) >= 5) {
    fisher_kmeans_pam50 <- fisher.test(cont_table_kmeans_pam50, simulate.p.value = nrow(cont_table_kmeans_pam50) > 2) #' @note Performs Fisher’s Exact Test for k-means vs. TCGA PAM50.
    fisher_results[["kmeans_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_kmeans_vs_PAM50",
      p_value = fisher_kmeans_pam50$p.value,
      statistic = NA_real_
    ) #' @note Stores k-means vs. PAM50 test results.
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (k-means vs. PAM50).\n") #' @note Warns if contingency table is too small or sparse.
  }

  if (all(dim(cont_table_hc_pam50) >= 2) && sum(cont_table_hc_pam50) >= 5) {
    fisher_hc_pam50 <- fisher.test(cont_table_hc_pam50, simulate.p.value = nrow(cont_table_hc_pam50) > 2) #' @note Performs Fisher’s Exact Test for hierarchical vs. TCGA PAM50.
    fisher_results[["hc_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hc_vs_PAM50",
      p_value = fisher_hc_pam50$p.value,
      statistic = NA_real_
    ) #' @note Stores hierarchical vs. PAM50 test results.
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (hierarchical vs. PAM50).\n") #' @note Warns if contingency table is too small or sparse.
  }

  if (all(dim(cont_table_hdb_pam50) >= 2) && sum(cont_table_hdb_pam50) >= 5) {
    fisher_hdb_pam50 <- fisher.test(cont_table_hdb_pam50, simulate.p.value = nrow(cont_table_hdb_pam50) > 2) #' @note Performs Fisher’s Exact Test for HDBSCAN vs. TCGA PAM50.
    fisher_results[["hdb_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_vs_PAM50",
      p_value = fisher_hdb_pam50$p.value,
      statistic = NA_real_
    ) #' @note Stores HDBSCAN vs. PAM50 test results.
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN vs. PAM50).\n") #' @note Warns if contingency table is too small or sparse.
  }

  if (all(dim(cont_table_hdb_pam50_sub) >= 2) && sum(cont_table_hdb_pam50_sub) >= 5) {
    fisher_hdb_pam50_sub <- fisher.test(cont_table_hdb_pam50_sub, simulate.p.value = nrow(cont_table_hdb_pam50_sub) > 2) #' @note Performs Fisher’s Exact Test for HDBSCAN PAM50 vs. TCGA PAM50.
    fisher_results[["hdb_pam50_vs_pam50"]] <- data.frame(
      test = "Fisher_Exact_hdb_pam50_vs_PAM50",
      p_value = fisher_hdb_pam50_sub$p.value,
      statistic = NA_real_
    ) #' @note Stores HDBSCAN PAM50 vs. PAM50 test results.
  } else {
    cat("[WARN] Insufficient data for Fisher’s test (HDBSCAN PAM50 vs. PAM50).\n") #' @note Warns if contingency table is too small or sparse.
  }
}

fisher_df <- do.call(rbind, fisher_results) #' @note Combines Fisher’s Exact Test results into a single data frame.
write.csv(fisher_df, file.path(outdir, "stats_results", "fisher_tests.csv"), row.names = FALSE) #' @note Saves Fisher’s Exact Test results to a CSV file.

#' @subsection Kruskal-Wallis and DSCF Tests:
#' Performs Kruskal-Wallis tests on marker gene expression across clusters, followed by DSCF pairwise tests for significant genes.
kw_results <- list() #' @note Initializes a list for Kruskal-Wallis test results.
dscf_results <- list() #' @note Initializes a list for DSCF test results.
for (method in c("kmeans", "hc", "hdb", "hdb_pam50")) {
  clusters <- switch(method,
                     kmeans = clusters_kmeans,
                     hc = clusters_hc,
                     hdb = clusters_hdb,
                     hdb_pam50 = clusters_hdb_pam50) #' @note Selects cluster assignments based on clustering method.
  mat <- if (method == "hdb_pam50") M_pam50 else V_adj #' @note Uses PAM50 matrix for HDBSCAN PAM50, otherwise VST matrix.
  genes <- if (method == "hdb_pam50") rownames(M_pam50) else marker_ens #' @note Uses PAM50 genes or marker Ensembl IDs.
  for (gene in genes) {
    gene_sym <- if (method == "hdb_pam50") gene else ens2sym_union$hgnc_gene_symbol[match(gene, ens2sym_union$ensembl)] #' @note Gets gene symbol, using Ensembl ID if missing.
    gene_sym <- ifelse(is.na(gene_sym) || gene_sym == "", gene, gene_sym) #' @note Falls back to Ensembl ID if symbol is missing or empty.
    expr <- mat[gene, ] #' @note Extracts expression values for the gene.
    if (sum(!is.na(expr)) < length(expr) * 0.5) {
      cat(sprintf("[WARN] Skipping Kruskal-Wallis for %s (%s, %s): too many NA values.\n", gene_sym, gene, method)) #' @note Warns if more than 50% of expression values are NA.
      next
    }
    kw_test <- kruskal.test(expr ~ clusters) #' @note Performs Kruskal-Wallis test to compare expression across clusters.
    if (is.na(kw_test$p.value)) {
      cat(sprintf("[WARN] Kruskal-Wallis failed for %s (%s, %s): insufficient variation.\n", gene_sym, gene, method)) #' @note Warns if test fails due to insufficient variation.
      next
    }
    kw_results[[paste(method, gene_sym, sep = "_")]] <- data.frame(
      test = sprintf("Kruskal_Wallis_%s_%s", method, gene_sym),
      gene = gene,
      symbol = gene_sym,
      p_value = kw_test$p.value,
      statistic = kw_test$statistic
    ) #' @note Stores Kruskal-Wallis test results.
    if (kw_test$p.value < 0.05) {
      keep <- !is.na(expr) #' @note Identifies samples with non-NA expression.
      expr_ok <- expr[keep] #' @note Subsets expression to non-NA samples.
      grp_ok <- droplevels(clusters[keep]) #' @note Subsets and drops unused cluster levels.

      if (nlevels(grp_ok) >= 2 && length(unique(grp_ok)) >= 2) {
        dscf_df <- NULL #' @note Initializes DSCF results data frame.
        dscf_test <- try(PMCMRplus::dscfAllPairsTest(expr_ok, grp_ok, p.adjust.method = "BH"), silent = TRUE) #' @note Performs DSCF test with BH correction on non-NA samples.
        if (!inherits(dscf_test, "try-error") && !is.null(dscf_test$p.value)) {
          p_mat <- try(as.matrix(dscf_test$p.value), silent = TRUE) #' @note Converts DSCF p-values to a matrix.
          if (inherits(p_mat, "try-error") || length(dim(p_mat)) != 2) next #' @note Skips if p-value matrix is invalid.
          if (nrow(p_mat) < 2 || ncol(p_mat) < 2) next #' @note Skips if matrix is too small.

          lev_present <- levels(grp_ok) #' @note Gets present cluster levels.
          if (is.null(rownames(p_mat)) || any(!nzchar(rownames(p_mat)))) {
            rownames(p_mat) <- lev_present[seq_len(nrow(p_mat))] #' @note Sets row names to present levels if missing.
          }
          if (is.null(colnames(p_mat)) || any(!nzchar(colnames(p_mat)))) {
            colnames(p_mat) <- lev_present[seq_len(ncol(p_mat))] #' @note Sets column names to present levels if missing.
          }

          pairs_idx <- which(upper.tri(p_mat, diag = FALSE), arr.ind = TRUE) #' @note Identifies upper triangle indices for pairwise comparisons.
          swap <- FALSE
          if (nrow(pairs_idx) == 0) {
            pairs_idx <- which(lower.tri(p_mat, diag = FALSE), arr.ind = TRUE) #' @note Falls back to lower triangle if upper is empty.
            swap <- TRUE
          }
          if (nrow(pairs_idx) == 0) next #' @note Skips if no pairwise comparisons are available.

          subtype1 <- rownames(p_mat)[pairs_idx[, 1]] #' @note Extracts first subtype for each pair.
          subtype2 <- colnames(p_mat)[pairs_idx[, 2]] #' @note Extracts second subtype for each pair.
          pvals <- p_mat[cbind(pairs_idx[, 1], pairs_idx[, 2])] #' @note Extracts p-values for pairs.

          if (swap) {
            tmp <- subtype1; subtype1 <- subtype2; subtype2 <- tmp #' @note Swaps labels for consistent direction if using lower triangle.
          }

          keep <- is.finite(pvals) & !is.na(subtype1) & !is.na(subtype2) #' @note Filters valid comparisons.
          if (!any(keep)) next #' @note Skips if no valid comparisons.

          dscf_df <- data.frame(
            test = sprintf("DSCF_%s_%s_%s_vs_%s", method, gene_sym, subtype1[keep], subtype2[keep]),
            gene = gene,
            symbol = gene_sym,
            subtype1 = subtype1[keep],
            subtype2 = subtype2[keep],
            p_value = pvals[keep],
            stringsAsFactors = FALSE
          ) #' @note Creates data frame with DSCF test results.

          if (is.data.frame(dscf_df) && nrow(dscf_df) > 0) {
            dscf_results[[paste(method, gene_sym, sep = "_")]] <- dscf_df #' @note Stores DSCF results.
          }
        }
      }
      p_box <- ggplot(data.frame(expr = mat[gene, ], subtype = clusters),
                      aes(x = subtype, y = expr, fill = subtype)) +
        geom_boxplot() + theme_bw() +
        labs(title = sprintf("Expression of %s (%s)", gene_sym, method), y = "VST Expression") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) #' @note Creates boxplot for significant genes.
      ggsave(file.path(outdir, "stats_results", sprintf("boxplot_%s_%s.pdf", method, gene_sym)), p_box, width = 7, height = 5) #' @note Saves boxplot to PDF.
    }
  }
}

kw_df <- tryCatch(do.call(rbind, kw_results), error = function(e) NULL) #' @note Combines Kruskal-Wallis results into a data frame, handling errors.
if (!is.null(kw_df) && nrow(kw_df) > 0) {
  write.csv(kw_df, file.path(outdir, "stats_results", "kruskal_wallis_tests.csv"), row.names = FALSE) #' @note Saves Kruskal-Wallis results to CSV.
} else {
  message("[KW] No significant Kruskal–Wallis results to write.") #' @note Logs if no significant results.
}

dscf_df <- tryCatch(do.call(rbind, dscf_results), error = function(e) NULL) #' @note Combines DSCF results into a data frame, handling errors.
if (!is.null(dscf_df) && nrow(dscf_df) > 0) {
  write.csv(dscf_df, file.path(outdir, "stats_results", "dscf_tests.csv"), row.names = FALSE) #' @note Saves DSCF results to CSV.
} else {
  message("[DSCF] No pairwise DSCF results to write.") #' @note Logs if no DSCF results.

#' @section Marker Heatmaps:
#' Generates heatmaps for PAM50 and marker gene expression across clusters and subtypes.
ens2sym <- make_ens2sym(dsmz_raw) #' @note Creates Ensembl-to-symbol mapping using DSMZ data.
ens_vec <- strip_ensver(as_df(dsmz_raw)$Ensembl_ID) #' @note Strips version numbers from DSMZ Ensembl IDs.
dup_ids <- sum(duplicated(ens_vec)) #' @note Counts duplicate Ensembl IDs.
dup_unique <- sum(table(ens_vec) > 1, na.rm = TRUE) #' @note Counts unique duplicate IDs.
cat(sprintf("[INFO] DSMZ Ensembl duplicate rows: %d | duplicate IDs: %d (collapsed)\n",
            dup_ids, dup_unique)) #' @note Logs duplicate ID information.

marker_ens <- rownames(ens2sym)[tolower(ens2sym$hgnc_gene_symbol) %in% tolower(marker_symbols)] #' @note Identifies Ensembl IDs for marker symbols.
marker_ens <- intersect(unique(marker_ens), rownames(V_adj)) #' @note Filters to Ensembl IDs present in VST matrix.

if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50 #' @note Uses PAM50 matrix for heatmap.
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark)) #' @note Ensures unique row names for heatmap.
  }
  ann <- data.frame(
    dataset = dataset_lab,
    kmeans_cluster = clusters_kmeans,
    hc_subtype = clusters_hc,
    hdb_subtype = clusters_hdb,
    hdb_pam50_subtype = clusters_hdb_pam50,
    row.names = colnames(V_mark)
  ) #' @note Creates annotation data frame with dataset and cluster assignments.
  safe_pdf(file.path(outdir, "heatmap_pam50_all_clusters.pdf"), {
    pheatmap(V_mark[, order(clusters_hdb_pam50), drop = FALSE],
             show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_col = ann,
             main = "PAM50 Gene Expression (k-means, HC, HDBSCAN, HDBSCAN PAM50)") #' @note Plots heatmap of PAM50 genes sorted by HDBSCAN PAM50 clusters.
  }) #' @note Saves heatmap to PDF using safe_pdf function.
} else {
  message("[PAM50 heatmap] No PAM50 genes found in V_adj — skipping.") #' @note Logs if no PAM50 genes are available.
}

#' @subsection PAM50/TCGA Subtype Correlation:
#' Computes Spearman correlations between cluster means and TCGA PAM50 subtype means.
if (!is.na(tcga_sub_col)) {
  lab_tcga <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Gets TCGA PAM50 subtype labels.
  sub_lvls <- sort(unique(na.omit(lab_tcga))) #' @note Extracts unique subtype levels.
  sub_means <- sapply(sub_lvls, function(s) rowMeans(V_tcga[hvgs, lab_tcga == s, drop = FALSE], na.rm = TRUE)) #' @note Computes mean expression per TCGA subtype for HVGs.
  cl_means_kmeans <- sapply(levels(clusters_kmeans), function(cn) rowMeans(V_adj[hvgs, clusters_kmeans == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes k-means cluster means for HVGs.
  cl_means_hc <- sapply(levels(clusters_hc), function(cn) rowMeans(V_adj[hvgs, clusters_hc == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes hierarchical cluster means for HVGs.
  cl_means_hdb <- sapply(levels(clusters_hdb), function(cn) rowMeans(V_adj[hvgs, clusters_hdb == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes HDBSCAN cluster means for HVGs.
  cl_means_hdb_pam50 <- sapply(levels(clusters_hdb_pam50), function(cn) rowMeans(M_pam50[, clusters_hdb_pam50 == cn, drop = FALSE], na.rm = TRUE)) #' @note Computes HDBSCAN PAM50 cluster means for PAM50 genes.
  g_common <- intersect(rownames(sub_means), hvgs) #' @note Finds common genes between subtype means and HVGs.
  g_common_pam50 <- intersect(rownames(sub_means), rownames(M_pam50)) #' @note Finds common genes between subtype means and PAM50 genes.
  corr_kmeans <- cor(cl_means_kmeans[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman") #' @note Computes Spearman correlation for k-means clusters.
  corr_hc <- cor(cl_means_hc[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman") #' @note Computes Spearman correlation for hierarchical clusters.
  corr_hdb <- cor(cl_means_hdb[g_common, , drop = FALSE], sub_means[g_common, , drop = FALSE], method = "spearman") #' @note Computes Spearman correlation for HDBSCAN clusters.
  corr_hdb_pam50 <- cor(cl_means_hdb_pam50[g_common_pam50, , drop = FALSE], sub_means[g_common_pam50, , drop = FALSE], method = "spearman") #' @note Computes Spearman correlation for HDBSCAN PAM50 clusters.
  write.csv(corr_kmeans, file.path(outdir, "cor_kmeans_vs_TCGA_subtype.csv")) #' @note Saves k-means correlation matrix.
  write.csv(corr_hc, file.path(outdir, "cor_hierarchical_vs_TCGA_subtype.csv")) #' @note Saves hierarchical correlation matrix.
  write.csv(corr_hdb, file.path(outdir, "cor_hdbscan_vs_TCGA_subtype.csv")) #' @note Saves HDBSCAN correlation matrix.
  write.csv(corr_hdb_pam50, file.path(outdir, "cor_hdbscan_pam50_vs_TCGA_subtype.csv")) #' @note Saves HDBSCAN PAM50 correlation matrix.
  safe_pdf(file.path(outdir, "heatmap_cluster_vs_TCGA_subtype_corr.pdf"), {
    par(mfrow = c(2, 2)) #' @note Sets 2x2 plot layout for heatmaps.
    pheatmap::pheatmap(corr_kmeans, main = "k-means vs TCGA Subtypes (ρ)") #' @note Plots k-means correlation heatmap.
    pheatmap::pheatmap(corr_hc, main = "Hierarchical vs TCGA Subtypes (ρ)") #' @note Plots hierarchical correlation heatmap.
    pheatmap::pheatmap(corr_hdb, main = "HDBSCAN vs TCGA Subtypes (ρ)") #' @note Plots HDBSCAN correlation heatmap.
    pheatmap::pheatmap(corr_hdb_pam50, main = "HDBSCAN PAM50 vs TCGA Subtypes (ρ)") #' @note Plots HDBSCAN PAM50 correlation heatmap.
    par(mfrow = c(1, 1)) #' @note Resets plot layout.
  }) #' @note Saves correlation heatmaps to PDF.
}

#' @section DSMZ to TCGA Mapping:
#' Maps DSMZ samples to TCGA subtypes using Spearman correlations.
tcga_mean <- rowMeans(V_tcga, na.rm = TRUE) #' @note Computes mean expression across TCGA samples.
rho_overall <- sapply(colnames(V_dsmz), function(s) {
  suppressWarnings(cor(V_dsmz[, s], tcga_mean, method = "spearman", use = "pairwise.complete.obs")) #' @note Computes Spearman correlation between each DSMZ sample and TCGA mean.
})

rho_by_sub <- NULL #' @note Initializes matrix for subtype-specific correlations.
best_tcga_subtype <- NULL #' @note Initializes vector for best-matching TCGA subtypes.
best_tcga_rho <- NULL #' @note Initializes vector for best correlation values.

if (!is.na(tcga_sub_col)) {
  lab <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Gets TCGA subtype labels.
  sub_levels <- sort(unique(na.omit(lab))) #' @note Extracts unique subtype levels.
  sub_means <- sapply(sub_levels, function(s) {
    rowMeans(V_tcga[, lab == s, drop = FALSE], na.rm = TRUE) #' @note Computes mean expression per TCGA subtype.
  }) #' @note Creates matrix of subtype means.
  rho_by_sub <- sapply(colnames(V_dsmz), function(s) {
    apply(sub_means, 2, function(ref) {
      suppressWarnings(cor(V_dsmz[, s], ref, method = "spearman", use = "pairwise.complete.obs")) #' @note Computes correlation between DSMZ sample and each TCGA subtype mean.
    })
  }) #' @note Creates matrix of correlations for each DSMZ sample.
  colnames(rho_by_sub) <- colnames(V_dsmz) #' @note Sets column names to DSMZ sample IDs.
  best_tcga_subtype <- apply(rho_by_sub, 2, function(v) names(which.max(v))) #' @note Identifies best-matching TCGA subtype for each DSMZ sample.
  best_tcga_rho <- apply(rho_by_sub, 2, max) #' @note Extracts maximum correlation value for each DSMZ sample.
  write.csv(rho_by_sub, file.path(outdir, "rho_dsmz_to_tcga_subtypes.csv")) #' @note Saves correlation matrix to CSV.
}

map_tbl <- tibble(
  sample = colnames(V_dsmz),
  rho_overall = unname(rho_overall),
  kmeans_cluster = clusters_kmeans[colnames(V_dsmz)],
  hc_subtype = clusters_hc[colnames(V_dsmz)],
  hdb_subtype = clusters_hdb[colnames(V_dsmz)],
  hdb_pam50_subtype = clusters_hdb_pam50[colnames(V_dsmz)]
) #' @note Creates data frame mapping DSMZ samples to clusters and overall TCGA correlation.
if (!is.null(best_tcga_subtype)) {
  map_tbl$best_tcga_subtype <- best_tcga_subtype[map_tbl$sample] #' @note Adds best-matching TCGA subtype.
  map_tbl$best_tcga_subtype_rho <- best_tcga_rho[map_tbl$sample] #' @note Adds best correlation values.
}
write.csv(map_tbl, file.path(outdir, "mapping_dsmz_to_tcga.csv"), row.names = FALSE) #' @note Saves mapping table to CSV.
print(head(map_tbl)) #' @note Prints first few rows of mapping table for inspection.

#' @section Marker Expression Heatmap:
#' Creates a heatmap of PAM50 gene expression for TCGA and DSMZ samples with cluster and subtype annotations.
if (nrow(M_pam50) > 0) {
  V_mark <- M_pam50 #' @note Uses PAM50 matrix for heatmap.
  if (any(duplicated(rownames(V_mark)))) {
    rownames(V_mark) <- make.unique(rownames(V_mark)) #' @note Ensures unique row names for heatmap.
  }
  tcga_sub_vec <- rep(NA_character_, ncol(V_tcga)) #' @note Initializes TCGA subtype vector with NA.
  if (!is.na(tcga_sub_col)) {
    tcga_sub_vec <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Populates TCGA subtype vector.
  }
  ann_tcga <- data.frame(
    dataset = "TCGA",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_tcga)]),
    hc_subtype = as.character(clusters_hc[colnames(V_tcga)]),
    hdb_subtype = as.character(clusters_hdb[colnames(V_tcga)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50[colnames(V_tcga)]),
    TCGA_subtype = tcga_sub_vec,
    row.names = colnames(V_tcga),
    check.names = FALSE
  ) #' @note Creates annotation data frame for TCGA samples.
  ann_dsmz <- data.frame(
    dataset = "DSMZ",
    kmeans_cluster = as.character(clusters_kmeans[colnames(V_dsmz)]),
    hc_subtype = as.character(clusters_hc[colnames(V_dsmz)]),
    hdb_subtype = as.character(clusters_hdb[colnames(V_dsmz)]),
    hdb_pam50_subtype = as.character(clusters_hdb_pam50[colnames(V_dsmz)]),
    TCGA_subtype = NA_character_,
    row.names = colnames(V_dsmz),
    check.names = FALSE
  ) #' @note Creates annotation data frame for DSMZ samples.
  cols_order <- c(colnames(V_tcga), colnames(V_dsmz)) #' @note Defines column order (TCGA followed by DSMZ).
  V_mark <- V_mark[, cols_order, drop = FALSE] #' @note Reorders PAM50 matrix columns.
  ann_col <- rbind(ann_tcga, ann_dsmz)[cols_order, , drop = FALSE] #' @note Combines and reorders annotations.
  non_empty_annot <- colSums(!is.na(ann_col)) > 0 #' @note Identifies non-empty annotation columns.
  ann_col <- if (any(non_empty_annot)) ann_col[, non_empty_annot, drop = FALSE] else NULL #' @note Filters out empty annotation columns.
  if (nrow(V_mark) >= 2 && ncol(V_mark) >= 2) {
    pdf(file.path(outdir, "heatmap_pam50_vst.pdf"), width = 12, height = 8) #' @note Opens PDF device for heatmap.
    pheatmap(V_mark, show_rownames = TRUE, show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = TRUE,
             annotation_col = ann_col,
             main = "Marker expression (VST) across TCGA-BRCA & DSMZ-BRCA") #' @note Plots PAM50 expression heatmap with annotations.
    dev.off() #' @note Closes PDF device.
  } else {
    message(sprintf("[marker heatmap] Skipping: matrix dim = %d x %d (need ≥ 2 x 2).",
                    nrow(V_mark), ncol(V_mark))) #' @note Logs if heatmap matrix is too small.
  }
}

#' @section Median Z-Scores Heatmap:
#' Creates a heatmap of median z-scores for marker genes by subtype/cluster.
#' @subsection Marker Diagnostics and Ensembl Mapping:
#' Maps marker gene symbols to Ensembl IDs using DSMZ and TCGA data.
cat("[MARKER] Requested symbols:", paste(marker_symbols, collapse=", "), "\n") #' @note Logs requested marker gene symbols.

map_dsmz <- as_df(dsmz_raw) %>%
  transmute(ensembl = strip_ensver(Ensembl_ID), hgnc_gene_symbol = as.character(gene_name)) %>%
  distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE) #' @note Creates DSMZ Ensembl-to-symbol mapping, removing duplicates.

map_tcga <- tryCatch({
  cat("[INFO] Attempting TCGA gene symbol mapping via rowData...\n") #' @note Logs attempt to map TCGA gene symbols.
  rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se)) #' @note Extracts TCGA rowData.
  symcol <- intersect(tolower(colnames(rd)), c("gene_name","symbol","hgnc_symbol")) #' @note Finds column containing gene symbols.
  if (length(symcol)) {
    result <- data.frame(
      ensembl = strip_ensver(rownames(rd)),
      hgnc_gene_symbol = as.character(rd[[colnames(rd)[match(symcol[1], tolower(colnames(rd)))] ]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(hgnc_gene_symbol) & hgnc_gene_symbol != "") %>%
      distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE) #' @note Creates TCGA Ensembl-to-symbol mapping.
    cat(sprintf("[INFO] TCGA rowData mapping successful: %d mappings using column '%s'\n",
                nrow(result), symcol[1])) #' @note Logs successful TCGA mapping.
    result
  } else {
    cat("[WARN] TCGA rowData mapping failed: no suitable symbol column found\n") #' @note Warns if no symbol column is found.
    NULL
  }
}, error = function(e) {
  cat(sprintf("[WARN] TCGA rowData mapping failed: %s\n", e$message)) #' @note Warns if TCGA mapping fails.
  NULL
}) #' @note Handles errors in TCGA mapping.

ens2sym_union <- bind_rows(filter(map_dsmz, !is.na(hgnc_gene_symbol)),
                           filter(map_tcga, !is.null(map_tcga))) %>%
  distinct(hgnc_gene_symbol, .keep_all = TRUE) %>% distinct(ensembl, .keep_all = TRUE) #' @note Combines DSMZ and TCGA mappings, removing duplicates.

cat(sprintf("[INFO] Marker gene symbol mapping summary: %d total mappings\n", nrow(ens2sym_union))) #' @note Logs total number of mappings.

marker_ens <- ens2sym_union$ensembl[
  tolower(ens2sym_union$hgnc_gene_symbol) %in% tolower(marker_symbols)
] #' @note Identifies Ensembl IDs for marker symbols.
marker_ens <- intersect(unique(marker_ens), rownames(V)) #' @note Filters to Ensembl IDs present in VST matrix.

cat(sprintf("[MARKER] Found %d of %d PAM50/PAM-like markers in V\n",
            length(marker_ens), length(marker_symbols))) #' @note Logs number of matched marker genes.
if (length(marker_ens)) {
  cat("[MARKER] Example matches:", paste(head(marker_ens, 5), collapse=", "), "...\n") #' @note Logs example matched Ensembl IDs.
} else {
  stop("[FATAL] No marker genes from the requested list are present in V. Check symbol mapping.") #' @note Errors if no marker genes are found.
}

V_tcga_z <- zscore_by_gene(V_tcga) #' @note Applies z-score normalization to TCGA VST data by gene.
V_dsmz_z <- zscore_by_gene(V_dsmz) #' @note Applies z-score normalization to DSMZ VST data by gene.

#' @brief Aggregates median expression per group.
#' @param M Numeric matrix with genes as rows and samples as columns.
#' @param groups Factor of group labels (e.g., clusters or subtypes).
#' @return Matrix of median expression values per group.
aggregate_median <- function(M, groups) {
  grp <- droplevels(factor(groups)) #' @note Converts groups to factor and drops unused levels.
  lev <- levels(grp) #' @note Gets group levels.
  out <- sapply(lev, function(g) {
    matrixStats::rowMedians(M[, grp == g, drop = FALSE], na.rm = TRUE) #' @note Computes median expression per group.
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1) #' @note Ensures output is a matrix.
  rownames(out) <- rownames(M) #' @note Sets row names to gene IDs.
  colnames(out) <- lev #' @note Sets column names to group levels.
  out #' @return Median expression matrix.
}

if (!is.na(tcga_sub_col)) {
  tcga_groups <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Gets TCGA subtype labels.
  TCGA_med <- aggregate_median(V_tcga_z, tcga_groups) #' @note Computes median z-scores per TCGA subtype.
} else {
  TCGA_med <- matrixStats::rowMedians(V_tcga_z, na.rm = TRUE) #' @note Computes overall median z-scores for TCGA.
  TCGA_med <- matrix(TCGA_med, ncol = 1, dimnames = list(rownames(V_tcga_z), "TCGA_BRCA")) #' @note Converts to matrix with single column.
}

DSMZ_med <- aggregate_median(V_dsmz_z, clusters_kmeans[colnames(V_dsmz)]) #' @note Computes median z-scores per k-means cluster for DSMZ.

if (is.null(rownames(TCGA_med))) rownames(TCGA_med) <- rownames(V_tcga_z) #' @note Ensures TCGA median matrix has gene names.
if (is.null(rownames(DSMZ_med))) rownames(DSMZ_med) <- rownames(V_dsmz_z) #' @note Ensures DSMZ median matrix has gene names.

#' @subsection Median Z-Scores Heatmap:
#' Creates a heatmap of median z-scores for marker genes across TCGA subtypes and DSMZ clusters.
common_mark <- intersect(marker_ens, intersect(rownames(TCGA_med), rownames(DSMZ_med))) #' @note Finds common marker genes in TCGA and DSMZ median matrices.
if (length(common_mark) == 0L) {
  cat("[DEBUG] marker_ens length =", length(marker_ens), "\n") #' @note Logs debug information for marker genes.
  cat("[DEBUG] nrow(TCGA_med) =", nrow(TCGA_med), " nrow(DSMZ_med) =", nrow(DSMZ_med), "\n") #' @note Logs matrix dimensions.
  cat("[DEBUG] Example TCGA_med genes:", paste(head(rownames(TCGA_med), 5), collapse=", "), "\n") #' @note Logs example TCGA genes.
  cat("[DEBUG] Example DSMZ_med genes:", paste(head(rownames(DSMZ_med), 5), collapse=", "), "\n") #' @note Logs example DSMZ genes.
  stop("[FATAL] No overlap between marker and (TCGA_med ∩ DSMZ_med). Check ID harmonization (Ensembl, version stripping) and markers.") #' @note Errors if no common markers.
} else {
  H <- cbind(TCGA_med[common_mark, , drop = FALSE], DSMZ_med[common_mark, , drop = FALSE]) #' @note Combines TCGA and DSMZ median z-scores for common markers.
  pdf(file.path(outdir, "heatmap_marker_median_z.pdf"), width = 10, height = 7) #' @note Opens PDF device for heatmap.
  pheatmap(H, show_rownames = TRUE,
           main = "Median z-scores of marker genes by subtype/cluster") #' @note Plots heatmap of median z-scores.
  dev.off() #' @note Closes PDF device.
}

#' @section UMAP Embedding:
#' Creates UMAP visualizations of joint TCGA and DSMZ data using highly variable genes.
set.seed(42) #' @note Sets random seed for reproducibility.
genes_umap <- top_var_genes(V, n = 3000) #' @note Selects top 3000 variable genes for UMAP.
Xu <- t(V[genes_umap, , drop = FALSE]) #' @note Transposes VST matrix to samples x genes for UMAP.
emb <- uwot::umap(Xu, n_neighbors = 20, min_dist = 0.3, metric = "cosine") #' @note Performs UMAP with 20 neighbors and cosine metric.
emb <- as.data.frame(emb) #' @note Converts UMAP coordinates to data frame.
colnames(emb) <- c("UMAP1","UMAP2") #' @note Names UMAP dimensions.
emb$sample <- rownames(Xu) #' @note Adds sample names.
emb$dataset <- ifelse(emb$sample %in% colnames(V_tcga), "TCGA", "DSMZ") #' @note Adds dataset labels (TCGA or DSMZ).
emb$discovery_cluster <- clusters_kmeans[emb$sample] #' @note Adds k-means cluster assignments.

emb$TCGA_subtype <- NA #' @note Initializes TCGA subtype column with NA.
if (!is.na(tcga_sub_col)) {
  emb$TCGA_subtype[match(colnames(V_tcga), emb$sample)] <- as.character(colData(tcga_se)[colnames(V_tcga), tcga_sub_col]) #' @note Adds TCGA subtype labels for TCGA samples.
}

p1 <- ggplot(emb, aes(UMAP1, UMAP2, color = dataset)) +
  geom_point(alpha = 0.8, size = 1.6) +
  theme_bw() +
  ggtitle("UMAP (dataset)") #' @note Creates UMAP plot colored by dataset.

p2 <- ggplot(subset(emb, dataset=="TCGA"), aes(UMAP1, UMAP2, color = TCGA_subtype)) +
  geom_point(alpha = 0.8, size = 1.6) +
  theme_bw() +
  ggtitle("UMAP (TCGA BRCA subtype)") #' @note Creates UMAP plot for TCGA samples colored by subtype.

p3 <- ggplot(subset(emb, dataset=="DSMZ"), aes(UMAP1, UMAP2, color = discovery_cluster)) +
  geom_point(alpha = 0.9, size = 1.9) +
  theme_bw() +
  ggtitle("UMAP (discovery clusters)") #' @note Creates UMAP plot for DSMZ samples colored by k-means clusters.

ggsave(file.path(outdir, "umap_dataset.pdf"), p1, width = 7, height = 6) #' @note Saves dataset-colored UMAP plot.
ggsave(file.path(outdir, "umap_tcga_subtype.pdf"), p2, width = 7, height = 6) #' @note Saves TCGA subtype-colored UMAP plot.
ggsave(file.path(outdir, "umap_dsmz_cluster.pdf"), p3, width = 7, height = 6) #' @note Saves DSMZ cluster-colored UMAP plot.

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork) #' @note Loads patchwork for combining plots.
  combined <- (p1 | p2 | p3) + plot_layout(guides = "collect") #' @note Combines UMAP plots with shared legend.
  ggsave(file.path(outdir, "umap_combined.pdf"), combined, width = 15, height = 5) #' @note Saves combined UMAP plot.
} else {
  pdf(file.path(outdir, "umap_combined.pdf"), width = 15, height = 5) #' @note Opens PDF device for combined plot.
  print(p1) #' @note Prints dataset plot.
  print(p2) #' @note Prints TCGA subtype plot.
  print(p3) #' @note Prints DSMZ cluster plot.
  dev.off() #' @note Closes PDF device.
}

write.csv(emb, file.path(outdir, "umap_embeddings.csv"), row.names = FALSE) #' @note Saves UMAP coordinates to CSV.

cat("\n[OK] Done. Results in:", outdir, "\n") #' @note Logs completion and output directory.