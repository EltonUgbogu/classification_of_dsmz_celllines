#tcga_se_rds - TCGA SummarizedExperiment RDS filtered for Purity
#dsmz_rds - DSMZ count data from gdc-style processing pipeline
#dsmz_meta_csv - DSMZ metadata CSV from CellDrive #https://celldive.dsmz.de
#purity_tsv - Optional path to purity information from TCGA
#outdir - Output directory for results
#dsmz_cache_rds - Path to cached aligned DSMZ data

config <- list(                                         # Define configuration list for pipeline parameters
  # Inputs
  tcga_se_rds   = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds",  # Path to TCGA SummarizedExperiment RDS
  dsmz_rds      = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",  # Path to DSMZ count data RDS
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",    # Path to DSMZ metadata CSV
  purity_tsv    = NULL,                                # Optional path to purity TSV (sample, purity columns); NULL if not used
  # Outputs
  outdir        = "/home/chu25/dsmz/results/dsmz_tcga",  # Output directory for results
  dsmz_cache_rds = "/home/chu25/data/dsmz/DSMZ_aligned_cache.rds" # Path to cached aligned DSMZ data
)

# Create necessary directories
dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)  # Create output directory, allow recursive creation, suppress warnings

# =============================================================================
# DATA LOADING
# =============================================================================

# Load TCGA SummarizedExperiment and extract counts
#file_path - path to TCGA SummarizedExperiment RDS found in config
load_tcga_data <- function(file_path) {                 # Define function to load TCGA data
  cat("[INFO] Loading TCGA data...\n")                 # Print message indicating TCGA data loading
  if (!file.exists(file_path)) stop(sprintf("[ERROR] TCGA RDS not found: %s", file_path))  # Stop if file does not exist
  se <- readRDS(file_path)                             # Read TCGA SummarizedExperiment from RDS
  counts <- assay(se)                                  # Extract count matrix from SummarizedExperiment
  if (any(duplicated(rownames(counts)))) {                                     # check for duplicate genes  
    dup_genes <- rownames(counts)[duplicated(rownames(counts))]                  # get duplicate genes
    cat(sprintf("[INFO] Collapsing %d duplicate TCGA genes: %s\n",                 # print number of duplicate genes and first 5 duplicate genes
                length(dup_genes), paste(head(dup_genes, 5), collapse=", ")))    # print duplicate genes
    counts <- rowsum(counts, rownames(counts), reorder = TRUE)                   # collapse duplicate genes by summing counts for each gene
  }
  storage.mode(counts) <- "double"                                            # convert counts to double precision
  cat(sprintf("[INFO] TCGA dims: %d genes x %d samples\n", nrow(counts), ncol(counts)))  # print dimensions of TCGA count matrix
  list(se = se, counts = counts)                       # return list with SummarizedExperiment and counts
}

# Load DSMZ raw wide table (counts + gene columns) and metadata CSV
#counts_path - path to DSMZ count data RDS found in config
#meta_path - path to DSMZ metadata CSV found in config
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