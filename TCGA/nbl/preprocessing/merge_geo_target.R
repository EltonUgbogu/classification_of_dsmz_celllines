#!/usr/bin/env Rscript
# --------------------------------------------------------------
# merge_geo_target.R
# MERGE GEO + TARGET ON ENSEMBL IDs (MAX OVERLAP)
# KEEP HGNC IN rowData(se)
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(data.table)
  library(SummarizedExperiment)
  library(arrow)
  library(tidyr)
})

## ------------------- USER SETTINGS ---------------------------------
geo_rds      <- "/home/chu25/data/nbl/nbl_raw_count_matrix.rds"
target_rds   <- "/home/chu25/data/nbl/TARGET-NBL/TARGET_NBL_STAR_Counts_data.rds"
map_tsv      <- "/home/chu25/data/nbl/final_nbl_ensembl_to_hgnc.tsv"
out_se_rds   <- "/home/chu25/data/nbl/nbl_combined_count_matrix_ensembl.rds"
## -----------------------------------------------------------------

message("[INFO] Loading GEO matrix ...")
geo_mat <- readRDS(geo_rds)
if (!is.matrix(geo_mat)) stop("GEO matrix is not a base R matrix!")

# --- ENSURE SAMPLES × GENES ---
if (nrow(geo_mat) < ncol(geo_mat)) {
  message("   GEO: samples × genes: ", nrow(geo_mat), " x ", ncol(geo_mat))
} else {
  message("   GEO: genes × samples → transposing...")
  geo_mat <- t(geo_mat)
}

message("[INFO] Loading TARGET object ...")
target_obj <- readRDS(target_rds)

# --- Extract counts from TARGET ---
if (inherits(target_obj, "SummarizedExperiment")) {
  assay_names <- names(assays(target_obj))
  message("   TARGET assays: ", paste(assay_names, collapse = ", "))
  target_mat <- assay(target_obj, "unstranded")
} else if (is.matrix(target_obj)) {
  target_mat <- target_obj
} else {
  stop("Unknown TARGET object type")
}

# --- ENSURE SAMPLES × GENES ---
if (nrow(target_mat) < ncol(target_mat)) {
  message("   TARGET: samples × genes: ", nrow(target_mat), " x ", ncol(target_mat))
} else {
  message("   TARGET: genes × samples → transposing...")
  target_mat <- t(target_mat)
}

message("[INFO] GEO matrix (final):    ", nrow(geo_mat),    " samples x ", ncol(geo_mat),    " genes")
message("[INFO] TARGET matrix (final): ", nrow(target_mat), " samples x ", ncol(target_mat), " genes")

## -----------------------------------------------------------------
## 1. Extract clean Ensembl IDs (strip version)
## -----------------------------------------------------------------
geo_genes    <- sub("\\.\\d+$", "", colnames(geo_mat))
target_genes <- sub("\\.\\d+$", "", colnames(target_mat))

colnames(geo_mat)    <- geo_genes
colnames(target_mat) <- target_genes

## -----------------------------------------------------------------
## 2. Find common Ensembl IDs
## -----------------------------------------------------------------
common_ensembl <- intersect(geo_genes, target_genes)
if (length(common_ensembl) == 0) stop("No common Ensembl IDs!")

geo_mat    <- geo_mat[, common_ensembl,   drop = FALSE]
target_mat <- target_mat[, common_ensembl, drop = FALSE]

message("[INFO] Common Ensembl genes: ", length(common_ensembl))

## -----------------------------------------------------------------
## 3. Load HGNC map → rowData
## -----------------------------------------------------------------
message("[INFO] Loading Ensembl → HGNC map ...")
hgnc_map <- fread(map_tsv, sep = "\t", header = TRUE, data.table = FALSE)
hgnc_map <- hgnc_map[!is.na(hgnc_map$HGNC_Symbol) & hgnc_map$HGNC_Symbol != "", ]
hgnc_dict <- setNames(hgnc_map$HGNC_Symbol, hgnc_map$Ensembl_ID)

row_data <- data.frame(
  Ensembl_ID  = common_ensembl,
  HGNC_Symbol = hgnc_dict[common_ensembl],
  stringsAsFactors = FALSE
)
row_data$HGNC_Symbol[is.na(row_data$HGNC_Symbol)] <- ""

message("[INFO] HGNC mapped: ", sum(row_data$HGNC_Symbol != ""), " / ", nrow(row_data), " genes")

## -----------------------------------------------------------------
## 4. Sample metadata
## -----------------------------------------------------------------
geo_coldata    <- data.frame(sample = rownames(geo_mat),    source = "GEO",    stringsAsFactors = FALSE)
target_coldata <- data.frame(sample = rownames(target_mat), source = "TARGET", stringsAsFactors = FALSE)
coldata        <- rbind(geo_coldata, target_coldata)
rownames(coldata) <- coldata$sample

## -----------------------------------------------------------------
## 5. Combine counts (samples × genes)
## -----------------------------------------------------------------
combined_samples_x_genes <- rbind(geo_mat, target_mat)  # 282 × N

## -----------------------------------------------------------------
## 6. Transpose to genes × samples
## -----------------------------------------------------------------
combined_counts <- t(combined_samples_x_genes)  # N × 282

# --- ALIGN rowData to rownames ---
row_data <- row_data[match(rownames(combined_counts), row_data$Ensembl_ID), , drop = FALSE]
stopifnot(nrow(row_data) == nrow(combined_counts))
stopifnot(identical(rownames(combined_counts), row_data$Ensembl_ID))

## -----------------------------------------------------------------
## 7. Build SummarizedExperiment
## -----------------------------------------------------------------
se <- SummarizedExperiment(
  assays  = list(counts = combined_counts),
  colData = DataFrame(coldata),
  rowData = DataFrame(row_data)
)

# Explicit rownames
rownames(se) <- row_data$Ensembl_ID
rownames(rowData(se)) <- row_data$Ensembl_ID

message("[INFO] Final SE: ", nrow(se), " genes (Ensembl) x ", ncol(se), " samples")
message("[INFO]   → HGNC in rowData: ", sum(rowData(se)$HGNC_Symbol != ""), " symbols")
message("[INFO] Saving to ", out_se_rds)
saveRDS(se, out_se_rds)
message("[SUCCESS] Merged on Ensembl, HGNC in metadata!")