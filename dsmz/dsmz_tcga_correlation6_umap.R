suppressPackageStartupMessages({
  library(matrixStats)
  library(uwot)
  library(RANN)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(edgeR)  
  library(limma)  
  library(sva)
  library(RColorBrewer)
  library(SummarizedExperiment)
  library(tibble)
})

set.seed(42)


# Helpers
`%||%` <- function(x, y) if (is.null(x)) y else x      # Define custom operator to return y if x is NULL

# Log2-CPM helper using edgeR::cpm (prior.count avoids -Inf for zeros)
logCPM <- function(x, prior.count = 1, lib.size = NULL) {  # Define function to compute log2-CPM
  if (!is.matrix(x)) x <- as.matrix(x)                  # Convert input to matrix if not already
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size)  # Compute log2-CPM using edgeR
}

config <- list(                                         # Define configuration list for pipeline parameters
  # Inputs
  tcga_se_rds   = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds",  # Path to TCGA SummarizedExperiment RDS
  dsmz_rds      = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",  # Path to DSMZ count data RDS
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",    # Path to DSMZ metadata CSV
  purity_tsv    = NULL,                                # Optional path to purity TSV (sample, purity columns); NULL if not used

  # Outputs
  outdir        = "/home/chu25/dsmz/results/correlation6",  # Output directory for results
  dsmz_cache_rds = "/home/chu25/data/dsmz/DSMZ_aligned_cache.rds",  # Path to cached aligned DSMZ data

  # Analysis
  var_gene_set  = "tcga_5k",                           # Gene set for variable gene selection ("tcga_5k", "tcga_10k", "all")
  min_shared    = 1000,                                # Minimum number of shared genes required
  batch_adjust_method = "combat_seq",                  # Batch correction method ("none", "combat_seq", or "combat")
  force_realign = FALSE                                # Force realignment even if cache exists (set to TRUE to bypass cache)
)

fig_dir <- file.path(config$outdir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Check output directory
if (!dir.exists(config$outdir) && !dir.create(config$outdir, recursive = TRUE)) {
  stop(sprintf("[ERROR] Cannot create output directory: %s", config$outdir))
}

# Validate configuration parameters
valid_batch_methods <- c("none", "combat_seq", "combat")
if (!config$batch_adjust_method %in% valid_batch_methods) {
  stop(sprintf("[ERROR] Invalid batch_adjust_method: %s. Must be one of %s", 
               config$batch_adjust_method, paste(valid_batch_methods, collapse = ", ")))
}

valid_gene_sets <- c("tcga_5k", "tcga_10k", "all")
if (!config$var_gene_set %in% valid_gene_sets) {
  stop(sprintf("[ERROR] Invalid var_gene_set: %s. Must be one of %s", 
               config$var_gene_set, paste(valid_gene_sets, collapse = ", ")))
}


# Load TCGA SummarizedExperiment and extract counts
load_tcga_data <- function(file_path) {                 # Define function to load TCGA data
  cat("[INFO] Loading TCGA data...\n")                 # Print message indicating TCGA data loading
  if (!file.exists(file_path)) stop(sprintf("[ERROR] TCGA RDS not found: %s", file_path))  # Stop if file does not exist
  tryCatch({
    se <- readRDS(file_path)                             # Read TCGA SummarizedExperiment from RDS
    counts <- assay(se)                                  # Extract count matrix from SummarizedExperiment
    cat(sprintf("[INFO] TCGA dims: %d genes x %d samples\n", nrow(counts), ncol(counts)))  # Print dimensions of TCGA count matrix
    list(se = se, counts = counts)                       # Return list with SummarizedExperiment and counts
  }, error = function(e) {
    stop(sprintf("[ERROR] Failed to load TCGA RDS: %s", e$message))
  })
}

# Load DSMZ raw wide table (counts + gene columns) and metadata CSV
load_dsmz_data <- function(counts_path, meta_path) {    # Define function to load DSMZ data
  cat("[INFO] Loading DSMZ data...\n")                 # Print message indicating DSMZ data loading
  if (!file.exists(counts_path)) stop(sprintf("[ERROR] DSMZ RDS not found: %s", counts_path))  # Stop if counts file does not exist
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))  # Stop if metadata file does not exist
  raw <- readRDS(counts_path)                          # Read DSMZ raw count data from RDS
  meta <- read.csv(meta_path)                          # Read DSMZ metadata from CSV
  stopifnot(all(c("Ensembl_ID","gene_name") %in% colnames(raw)))  # Verify required columns in raw data
  stopifnot("sample_name" %in% names(meta))            # Verify sample_name column in metadata
  cat(sprintf("[INFO] DSMZ table: %d rows x %d cols\n", nrow(raw), ncol(raw)))  # Print dimensions of raw DSMZ data
  list(raw = raw, meta = meta)                         # Return list with raw counts and metadata
}

harmonize_gene_ids <- function(tcga_counts) {           # Define function to harmonize TCGA gene IDs
  cat("[INFO] Harmonizing TCGA gene IDs...\n")         # Print message indicating gene ID harmonization
  g <- sub("\\..*$","", rownames(tcga_counts))         # Remove version suffixes from Ensembl IDs
  if (any(duplicated(g))) {                            # Check for duplicate gene IDs
    tcga_counts <- rowsum(tcga_counts, g, reorder = TRUE)  # Sum counts for duplicate genes
    rownames(tcga_counts) <- sort(unique(g))           # Set unique sorted gene IDs as row names
  } else {
    rownames(tcga_counts) <- g                         # Set harmonized gene IDs as row names
  }
  tcga_counts                                          # Return harmonized count matrix
}

# Build a numeric DSMZ count matrix from the raw wide table
build_dsmz_matrix <- function(dsmz_raw) {               # Define function to build DSMZ count matrix
  cat("[INFO] Building DSMZ count matrix...\n")        # Print message indicating matrix construction
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version")  # Define annotation columns to exclude
  sample_cols <- setdiff(colnames(dsmz_raw), annot)    # Identify sample columns
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) as.numeric(as.character(x)))  # Convert sample columns to numeric
  M <- as.matrix(dsmz_raw[, sample_cols, drop=FALSE])  # Create matrix from sample columns
  rownames(M) <- dsmz_raw$Ensembl_ID                   # Set Ensembl IDs as row names
  if (any(duplicated(rownames(M))))                    # Check for duplicate gene IDs
    M <- rowsum(M, rownames(M), reorder = TRUE)        # Sum counts for duplicate genes
  keep <- vapply(as.data.frame(M), function(v) any(!is.na(v)), logical(1))  # Identify columns with non-NA values
  if (!all(keep)) {                                    # Check if any columns are all NA
    cat(sprintf("[WARN] Dropping %d DSMZ columns that are all NA\n", sum(!keep)))  # Warn about dropped columns
    M <- M[, keep, drop=FALSE]                         # Drop all-NA columns
  }
  M                                                    # Return DSMZ count matrix
}

align_metadata <- function(dsmz_meta, dsmz_counts, outdir) {  # Define function to align metadata and counts
  cat("[INFO] Aligning metadata to counts...\n")       # Print message indicating alignment
  dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)  # Create sample_id column from sample_name
  matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))  # Find common samples between metadata and counts
  md_only <- setdiff(dsmz_meta$sample_id, colnames(dsmz_counts))  # Find metadata samples not in counts
  ct_only <- setdiff(colnames(dsmz_counts), dsmz_meta$sample_id)  # Find count samples not in metadata
  if (length(md_only))                                 # Check if there are unmatched metadata samples
    write.csv(tibble(sample_name = md_only),            # Save unmatched metadata samples to CSV
              file.path(outdir, "unmatched_metadata_samples.csv"),
              row.names = FALSE)
  if (length(ct_only))                                 # Check if there are unmatched count samples
    write.csv(tibble(sample_name = ct_only),            # Save unmatched count samples to CSV
              file.path(outdir, "unmatched_count_columns.csv"),
              row.names = FALSE)
  dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched)  # Filter metadata to matched samples
  dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop=FALSE]  # Subset counts to matched samples
  stopifnot(identical(colnames(dsmz_counts), dsmz_meta$sample_id))  # Verify counts and metadata alignment
  cat(sprintf("[INFO] Matched: %d | Unmatched(meta): %d | Unmatched(counts): %d\n",
              length(matched), length(md_only), length(ct_only)))  # Print alignment summary
  list(meta = dsmz_meta, counts = dsmz_counts)         # Return list with aligned metadata and counts
}

# Load or create aligned DSMZ data with caching
load_or_align_dsmz_data <- function(dsmz_raw, dsmz_meta, cache_file, outdir, force_realign = FALSE) {  # Define function to load or create aligned DSMZ data
  if (!force_realign && file.exists(cache_file)) {     # Check if cache file exists and realignment is not forced
    cat("[INFO] Loading cached aligned DSMZ data from:", cache_file, "\n")  # Print message indicating cache loading
    cached_data <- readRDS(cache_file)                 # Load cached data
    # Verify the cached data has the expected structure
    if (all(c("meta", "counts") %in% names(cached_data)) && 
        nrow(cached_data$meta) > 0 && ncol(cached_data$counts) > 0) {
      cat(sprintf("[INFO] Cached data: %d samples, %d genes\n", 
                  nrow(cached_data$meta), nrow(cached_data$counts)))  # Print cached data summary
      return(cached_data)                              # Return cached data
    } else {
      cat("[WARN] Cached data appears corrupted, realigning...\n")  # Warn about corrupted cache
    }
  }
  
  cat("[INFO] Building DSMZ count matrix and aligning metadata...\n")  # Print message indicating alignment
  dsmz_counts <- build_dsmz_matrix(dsmz_raw)          # Build DSMZ count matrix
  storage.mode(dsmz_counts) <- "double"               # Ensure double precision
  aligned_data <- align_metadata(dsmz_meta, dsmz_counts, outdir)  # Align metadata with counts
  
  # Save aligned data to cache
  cat("[INFO] Saving aligned DSMZ data to cache:", cache_file, "\n")  # Print message indicating cache saving
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)  # Create cache directory
  saveRDS(aligned_data, cache_file)                   # Save aligned data to RDS
  aligned_data                                        # Return aligned data
}

# Load tumor purity vector from SummarizedExperiment or external TSV
load_purity_data <- function(tcga_se, purity_file=NULL) {  # Define function to load tumor purity data
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

# Adjust TCGA LOG matrix for purity: (1) drop negative-purity-correlated genes, (2) regress out infiltration
# NOTE: operates on log2-CPM (not raw counts)
purity_adjust_log <- function(tcga_log, purity) {       # Define function to adjust TCGA log2-CPM for purity
  if (is.null(purity)) return(tcga_log)                # Return input if no purity data
  cat("[INFO] Applying purity adjustment on log2-CPM...\n")  # Print message indicating adjustment
  keep_samp <- intersect(colnames(tcga_log), names(purity)[!is.na(purity)])  # Find samples with valid purity values
  if (length(keep_samp) < 3) {                         # Check if enough samples for adjustment
    cat("[WARN] Too few samples with valid purity to adjust; returning input.\n")  # Warn if insufficient samples
    return(tcga_log)                                   # Return input matrix
  }
  M <- tcga_log[, keep_samp, drop=FALSE]               # Subset log matrix to valid samples
  pu <- purity[keep_samp]                              # Subset purity to valid samples
  # Identify & drop genes strongly NEG correlated with purity (likely infiltration)
  ct <- apply(M, 1, function(v) {
    if (var(v) == 0 || var(pu) == 0) return(list(p.value = 1, estimate = 0))  # Handle zero variance cases
    suppressWarnings(cor.test(as.numeric(v), pu, method="spearman"))  # Compute Spearman correlations
  })
  pv <- vapply(ct, `[[`, numeric(1), "p.value")        # Extract p-values
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))  # Extract correlation coefficients
  pv[is.na(pv)] <- 1                                   # Replace NA p-values with 1 (non-significant)
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
  out <- tcga_log                                      # Initialize output matrix
  common_g <- intersect(rownames(Madj), rownames(out))  # Find common genes
  out[common_g, colnames(Madj)] <- Madj[common_g, ]    # Update adjusted values
  out                                                  # Return adjusted matrix
}

tcga <- load_tcga_data(config$tcga_se_rds)         # Load TCGA data
dsmz <- load_dsmz_data(config$dsmz_rds, config$dsmz_meta_csv)  # Load DSMZ data
tcga_counts <- harmonize_gene_ids(tcga$counts); storage.mode(tcga_counts) <- "double"  # Harmonize TCGA gene IDs and ensure double precision

# Check matrix sizes for memory warnings
if (nrow(tcga_counts) > 50000 || ncol(tcga_counts) > 10000) {
  cat("[WARN] Large TCGA matrix detected; consider subsampling for testing.\n")
}
  
  # Load or create aligned DSMZ data with caching
aligned_dsmz <- load_or_align_dsmz_data(dsmz$raw, dsmz$meta, config$dsmz_cache_rds, config$outdir, config$force_realign)  # Load or create aligned DSMZ data
dsmz_meta <- aligned_dsmz$meta                     # Extract aligned metadata
dsmz_counts <- aligned_dsmz$counts                 # Extract aligned counts

# Check DSMZ matrix size for memory warnings
if (nrow(dsmz_counts) > 50000 || ncol(dsmz_counts) > 10000) {
  cat("[WARN] Large DSMZ matrix detected; consider subsampling for testing.\n")
}
  
dsmz_log_counts <- logCPM(dsmz_counts)                           # Compute log2-CPM for uncorrected data
tcga_log_counts <- logCPM(tcga_counts)                           # Compute log2-CPM for uncorrected data

# Align genes before combining to prevent misalignment
common_genes <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
if (length(common_genes) < config$min_shared) {
  stop(sprintf("[ERROR] Too few shared genes: %d (required: %d)", 
               length(common_genes), config$min_shared))
}
tcga_counts <- tcga_counts[common_genes, , drop=FALSE]
dsmz_counts <- dsmz_counts[common_genes, , drop=FALSE]

Xc_raw <- cbind(tcga_counts, dsmz_counts)          # Combine TCGA and DSMZ raw counts
batch <- factor(c(rep("TCGA", ncol(tcga_counts)), rep("DSMZ", ncol(dsmz_counts))))  # Define batch factor

# Check combined matrix size for memory warnings
if (nrow(Xc_raw) > 50000 || ncol(Xc_raw) > 10000) {
  cat("[WARN] Large combined matrix detected; consider subsampling for testing.\n")
}

# --------------------- Batch correction (ONE place) --------------------------
cat("[INFO] Batch adjustment (global)...\n")       # Print message indicating batch correction
Xc <- Xc_raw                                   # Use raw counts
Xc[is.na(Xc)] <- 0                            # Replace NA with 0

# Verify that counts are non-negative integers for ComBat_seq
if (!all(Xc >= 0) || !all(Xc == floor(Xc))) {
  stop("[ERROR] Input to ComBat_seq contains non-integer or negative values")
}

set.seed(42)                                   # Set seed for reproducibility
Xc_adj <- sva::ComBat_seq(as.matrix(Xc), batch = batch)  # Apply ComBat-seq
M1_log <- logCPM(Xc_adj)                       # Log-transform corrected counts
# --------------------- Optional purity adjustment on TCGA (log scale) -------
purity <- load_purity_data(tcga$se, config$purity_tsv)  # Load purity data
if (!is.null(purity)) {                            # Check if purity data exists
  tcga_cols <- colnames(tcga_counts)               # Get TCGA sample columns
  tcga_log_adj <- purity_adjust_log(M1_log[, tcga_cols, drop=FALSE], purity)  # Adjust TCGA data for purity
  M1_log[, tcga_cols] <- tcga_log_adj             # Update adjusted TCGA data
}

dsmz_cols <- colnames(dsmz_counts)                 # Get DSMZ sample columns
  
tcga_log_mat <- M1_log[, tcga_cols, drop=FALSE]    # Subset TCGA log matrix
dsmz_log_mat <- M1_log[, dsmz_cols, drop=FALSE]    # Subset DSMZ log matrix


# -------------------- 0) Prepare matrices & features -------------------------
common_genes <- intersect(rownames(dsmz_log_mat), rownames(tcga_log_mat))
Xd <- dsmz_log_mat[common_genes, , drop=FALSE]
Xt <- tcga_log_mat[common_genes, , drop=FALSE]

# Gene selection based on config$var_gene_set
if (config$var_gene_set == "tcga_5k") {
  iqr <- matrixStats::rowIQRs(Xt)
  sel <- names(sort(setNames(iqr, rownames(Xt)), decreasing = TRUE))[seq_len(min(5000, length(iqr)))]
} else if (config$var_gene_set == "tcga_10k") {
  iqr <- matrixStats::rowIQRs(Xt)
  sel <- names(sort(setNames(iqr, rownames(Xt)), decreasing = TRUE))[seq_len(min(10000, length(iqr)))]
} else if (config$var_gene_set == "all") {
  sel <- common_genes
} else {
  stop(sprintf("[ERROR] Invalid var_gene_set: %s", config$var_gene_set))
}
Xd <- Xd[sel, , drop=FALSE]
Xt <- Xt[sel, , drop=FALSE]

# Filter zero-variance genes before scaling
zero_var_d <- matrixStats::rowVars(Xd) == 0
zero_var_t <- matrixStats::rowVars(Xt) == 0
zero_var <- zero_var_d | zero_var_t
if (any(zero_var)) {
  cat(sprintf("[WARN] Removing %d genes with zero variance\n", sum(zero_var)))
  Xd <- Xd[!zero_var, , drop=FALSE]
  Xt <- Xt[!zero_var, , drop=FALSE]
}

# Center/scale per gene before PCA (optional but helps stability)
Xd_sc <- t(scale(t(Xd), center = TRUE, scale = TRUE))
Xt_sc <- t(scale(t(Xt), center = TRUE, scale = TRUE))

# Check memory usage before PCA
if (nrow(Xd_sc) > 50000 || ncol(Xd_sc) + ncol(Xt_sc) > 10000) {
  cat("[WARN] Large scaled matrix for PCA; consider subsampling genes or samples.\n")
}

# -------------------- 1) PCA (combined, unsupervised) ------------------------
# PCA on samples (transpose to samples × genes)
tryCatch({
  pc <- prcomp(t(cbind(Xd_sc, Xt_sc)), center = FALSE, scale. = FALSE)
  pcs_keep <- 30
  PC_all   <- pc$x[, 1:pcs_keep, drop = FALSE]

  # Split back to DSMZ/TCGA PC scores (rows are samples in same order as columns of cbind)
  n_dsmz <- ncol(Xd_sc)
  PC_dsmz <- PC_all[1:n_dsmz, , drop=FALSE]
  PC_tcga <- PC_all[(n_dsmz + 1):nrow(PC_all), , drop=FALSE]

  # Make sure row order matches sample IDs
  stopifnot(identical(rownames(PC_dsmz), colnames(Xd_sc)))
  stopifnot(identical(rownames(PC_tcga), colnames(Xt_sc)))
  
  cat("[INFO] PCA completed successfully\n")
}, error = function(e) {
  stop(sprintf("[ERROR] PCA failed: %s", e$message))
})

# -------------------- 2) DSMZ-anchored UMAP on PCs ---------------------------
# Check memory usage before UMAP
if (nrow(PC_dsmz) > 10000 || ncol(PC_dsmz) > 30) {
  cat("[WARN] Large PC matrix for UMAP; consider reducing PCs or samples.\n")
}

tryCatch({
  if (nrow(PC_dsmz) < 5) {
    stop("[ERROR] Too few DSMZ samples for UMAP (need at least 5)")
  }
  umap_model <- uwot::umap(
    X = PC_dsmz,              # DSMZ samples × 30 PCs
    metric = "euclidean",     # PCs → euclidean is fine; cosine also OK
    n_neighbors = 15,
    min_dist = 0.2,
    scale = FALSE,            # PCs already scaled by PCA
    ret_model = TRUE,
    verbose = TRUE
  )
  emb_dsmz <- umap_model$embedding
  colnames(emb_dsmz) <- c("UMAP1","UMAP2")

  # Project ALL TCGA PCs into DSMZ UMAP space (unsupervised)
  stopifnot(identical(colnames(PC_dsmz), colnames(PC_tcga)))  # Verify PC columns match
  emb_tcga <- tryCatch({
    uwot::umap_transform(PC_tcga, umap_model)
  }, error = function(e) {
    stop(sprintf("[ERROR] UMAP transform failed: %s", e$message))
  })
  colnames(emb_tcga) <- c("UMAP1","UMAP2")
  
  cat("[INFO] UMAP embedding completed successfully\n")
}, error = function(e) {
  stop(sprintf("[ERROR] UMAP failed: %s", e$message))
})

# -------------------- 3) Build plotting metadata -----------------------------
dsmz_df <- dsmz_meta %>%
  transmute(sample_id = sample_id,
            dataset = "DSMZ",
            organ = ifelse(is.na(organ) | organ == "", "Other", as.character(organ))) %>%
  bind_cols(as.data.frame(emb_dsmz))

tcga_df <- tibble(
  sample_id = rownames(emb_tcga),
  dataset = "TCGA",
  organ = NA_character_
) %>% bind_cols(as.data.frame(emb_tcga))

# -------------------- 4) Unsupervised organ inference for TCGA (KNN on PCs) --
k <- 15
nn_pc <- RANN::nn2(data = PC_dsmz, query = PC_tcga, k = k)
dsmz_org_map <- dsmz_df %>% select(sample_id, organ) %>% tibble::column_to_rownames("sample_id")
nn_ids <- rownames(PC_dsmz)[nn_pc$nn.idx]    # DSMZ neighbor sample_ids per TCGA row
tcga_df$organ_inferred <- vapply(seq_len(nrow(nn_ids)), function(i){
  labs <- dsmz_org_map[nn_ids[i,], , drop=FALSE][,1]
  names(sort(table(labs), decreasing = TRUE))[1]
}, character(1))

# -------------------- 5) Plotting (two panels with clear legends) ------------
# Check number of points before plotting
if (nrow(dsmz_df) + nrow(tcga_df) > 100000) {
  cat("[WARN] Large number of points in UMAP plot; consider subsampling for visualization.\n")
}

# Palette (extendable)
pal <- c("#E69F00","#56B4E9","#009E73","#F0E442",
         "#0072B2","#D55E00","#CC79A7","#000000",
         "#999999","#CCBB44","#66CCEE","#AA3377")

# Dynamically extend palette if needed
n_colors_needed <- max(length(unique(dsmz_df$organ)), length(unique(tcga_df$organ_inferred)))
if (n_colors_needed > length(pal)) {
  pal <- c(pal, brewer.pal(min(n_colors_needed - length(pal), 8), "Set2"))
}

# Panel A: DSMZ colored by organ; TCGA grey
pA <- bind_rows(
    dsmz_df %>% mutate(plot_color = organ,      plot_alpha = 0.95, plot_size = 1.6),
    tcga_df %>%  mutate(plot_color = "TCGA",    plot_alpha = 0.30, plot_size = 0.8)
  ) %>%
  mutate(plot_color = factor(plot_color, levels = c(sort(unique(dsmz_df$organ)), "TCGA")),
         dataset    = factor(dataset, levels = c("DSMZ","TCGA"))) %>%
  ggplot(aes(UMAP1, UMAP2, color = plot_color, shape = dataset, alpha = plot_alpha, size = plot_size)) +
  geom_point(stroke = 0.25) +
  scale_color_manual(values = c(setNames(pal[seq_len(min(length(pal)-1, length(unique(dsmz_df$organ))))],
                                         sort(unique(dsmz_df$organ))),
                                "TCGA" = "grey70"),
                     name = "Organ (DSMZ) / TCGA") +
  scale_shape_manual(values = c(16, 17), name = "Dataset") +
  scale_alpha_identity() + scale_size_identity() +
  coord_equal() + theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  labs(title = "UMAP on Top 30 PCs (DSMZ-anchored)\nPanel A: DSMZ colored by organ; TCGA in grey",
       subtitle = "Unsupervised: UMAP fit on DSMZ PCs; all TCGA projected",
       x = "UMAP1", y = "UMAP2")

# Panel B: DSMZ colored by organ; TCGA colored by inferred organ (unsupervised)
organ_levels <- sort(unique(c(dsmz_df$organ, tcga_df$organ_inferred)))
pB <- bind_rows(
    dsmz_df %>% mutate(plot_color = organ,           plot_alpha = 0.95, plot_size = 1.6),
    tcga_df %>%  mutate(plot_color = organ_inferred, plot_alpha = 0.35, plot_size = 0.9)
  ) %>%
  mutate(plot_color = factor(plot_color, levels = organ_levels),
         dataset    = factor(dataset, levels = c("DSMZ","TCGA"))) %>%
  ggplot(aes(UMAP1, UMAP2, color = plot_color, shape = dataset, alpha = plot_alpha, size = plot_size)) +
  geom_point(stroke = 0.25) +
  scale_color_manual(values = setNames(pal[seq_len(min(length(pal), length(organ_levels)))], organ_levels),
                     name = "Organ / Inferred Organ") +
  scale_shape_manual(values = c(16, 17), name = "Dataset") +
  scale_alpha_identity() + scale_size_identity() +
  coord_equal() + theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  labs(title = "UMAP on Top 30 PCs (DSMZ-anchored)\nPanel B: TCGA colored by inferred organ (KNN on PCs)",
       subtitle = paste0("Unsupervised KNN (k=", k, ") in 30D PC space"),
       x = "UMAP1", y = "UMAP2")

plt <- pA + pB + plot_layout(ncol = 2, widths = c(1,1))

# Check write permissions before saving
if (!file.access(dirname(file.path(fig_dir, "UMAP_DSMZcentered_PCA30_dual.pdf")), 2) == 0) {
  stop(sprintf("[ERROR] No write permission for figure directory: %s", fig_dir))
}

tryCatch({
  ggsave(file.path(fig_dir, "UMAP_DSMZcentered_PCA30_dual.pdf"),  plt, width = 14, height = 6)
  ggsave(file.path(fig_dir, "UMAP_DSMZcentered_PCA30_dual.png"),  plt, width = 14, height = 6, dpi = 300)
  message("[INFO] Saved: ", file.path(fig_dir, "UMAP_DSMZcentered_PCA30_dual.{pdf,png}"))
}, error = function(e) {
  stop(sprintf("[ERROR] Failed to save plots: %s", e$message))
})
