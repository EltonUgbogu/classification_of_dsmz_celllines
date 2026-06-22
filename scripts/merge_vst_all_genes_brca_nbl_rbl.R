#!/usr/bin/env Rscript

# merge_vst_all_genes_brca_nbl_rbl.R
# ==============================================================================
# Combine BRCA, NBL, and RBL joint VST tumour/cell-line matrices into one
# samples x genes matrix on the intersection of Ensembl gene IDs (versions
# stripped), suitable for an unsupervised UMAP that does NOT depend on the
# pan-cancer DEG / feature-set selection.
#
# Filtering rules (none of them use any class label):
#   - orient each input as samples x genes;
#   - strip Ensembl version suffixes from gene IDs;
#   - within each cohort, deduplicate gene IDs by KEEPING THE FIRST OCCURRENCE
#     (deterministic: input column order is preserved);
#   - keep only genes present in ALL three cohorts;
#   - drop genes with any NA in any sample of the merged matrix;
#   - drop genes with zero variance across the merged matrix.
#
# Also writes a matching metadata TSV by subsetting an existing joint
# metadata TSV (no relabelling, no inference).
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option("--brca_rds",  type = "character", help = "BRCA joint VST RDS"),
  make_option("--nbl_rds",   type = "character", help = "NBL joint VST RDS"),
  make_option("--rbl_rds",   type = "character", help = "RBL joint VST RDS"),
  make_option("--joint_meta", type = "character",
              help = "Source joint metadata TSV (sample_id, cancer_type, sample_type[, cohort])"),
  make_option("--out_expr",  type = "character", help = "Output merged RDS path"),
  make_option("--out_meta",  type = "character", help = "Output filtered metadata TSV path"),
  make_option("--out_log",   type = "character", default = NULL,
              help = "Optional human-readable merge report path (defaults to <out_meta>.merge_log.tsv)")
)
opt <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# --- helpers ----------------------------------------------------------------

clean_ensg <- function(ids) sub("\\.\\d+$", "", ids)

read_vst <- function(path, label) {
  if (!file.exists(path)) stop("Missing input: ", path)
  obj <- readRDS(path)
  # Some upstream RDS objects wrap the matrix in a list with an `expr` element.
  if (is.list(obj) && "expr" %in% names(obj)) obj <- obj$expr
  X <- as.matrix(obj)
  rn <- rownames(X) %||% character(0)
  cn <- colnames(X) %||% character(0)
  # Orient as samples x genes (genes are Ensembl IDs starting "ENSG").
  rn_has_ensg <- length(rn) > 0 && mean(grepl("^ENSG", rn)) > 0.5
  cn_has_ensg <- length(cn) > 0 && mean(grepl("^ENSG", cn)) > 0.5
  if (rn_has_ensg && !cn_has_ensg) {
    X <- t(X)
    rn <- rownames(X); cn <- colnames(X)
    orient_decision <- "genes_in_rows -> transposed to samples_x_genes"
  } else if (cn_has_ensg && !rn_has_ensg) {
    orient_decision <- "samples_in_rows (already samples_x_genes)"
  } else if (rn_has_ensg && cn_has_ensg) {
    orient_decision <- "ENSG in BOTH dimensions; defaulting to as-is (samples_x_genes assumed)"
  } else {
    orient_decision <- "no ENSG in either dimension; defaulting to as-is"
  }
  colnames(X) <- clean_ensg(colnames(X))
  # Deduplicate gene columns deterministically: keep first occurrence.
  dup <- duplicated(colnames(X))
  n_dup <- sum(dup)
  if (n_dup > 0) X <- X[, !dup, drop = FALSE]
  cat(sprintf("[merge] %s: %d samples x %d genes  (orientation: %s; deduplicated %d genes, kept first occurrence)\n",
              label, nrow(X), ncol(X), orient_decision, n_dup))
  list(
    X = X,
    n_dup_dropped_keep_first = n_dup,
    orientation = orient_decision
  )
}

# --- load + align -----------------------------------------------------------

brca <- read_vst(opt$brca_rds, "BRCA")
nbl  <- read_vst(opt$nbl_rds,  "NBL")
rbl  <- read_vst(opt$rbl_rds,  "RBL")

shared <- Reduce(intersect, list(
  colnames(brca$X), colnames(nbl$X), colnames(rbl$X)
))
cat(sprintf("[merge] gene intersection: %d genes shared across BRCA(%d) / NBL(%d) / RBL(%d)\n",
            length(shared), ncol(brca$X), ncol(nbl$X), ncol(rbl$X)))
if (length(shared) <= 2) {
  stop("Intersection of gene IDs is empty or near-empty (", length(shared),
       "). Check input orientation and ENSG ID conventions.")
}

# Re-order each cohort's columns to the SAME shared gene order so rbind aligns.
brca_X <- brca$X[, shared, drop = FALSE]
nbl_X  <- nbl$X[,  shared, drop = FALSE]
rbl_X  <- rbl$X[,  shared, drop = FALSE]

merged <- rbind(brca_X, nbl_X, rbl_X)
cat(sprintf("[merge] merged matrix before NA + zero-var filtering: %d samples x %d genes\n",
            nrow(merged), ncol(merged)))

# --- NA + zero-variance filtering on the merged matrix ----------------------

na_cols <- colSums(is.na(merged)) > 0
n_dropped_na <- sum(na_cols)
if (n_dropped_na > 0) {
  merged <- merged[, !na_cols, drop = FALSE]
}
col_vars <- apply(merged, 2, stats::var, na.rm = FALSE)
zv_cols  <- !is.finite(col_vars) | col_vars == 0
n_dropped_zv <- sum(zv_cols)
if (n_dropped_zv > 0) {
  merged <- merged[, !zv_cols, drop = FALSE]
}
cat(sprintf("[merge] dropped %d NA-bearing genes and %d zero-variance genes (post-merge)\n",
            n_dropped_na, n_dropped_zv))
cat(sprintf("[merge] final merged matrix: %d samples x %d genes\n",
            nrow(merged), ncol(merged)))

# --- write merged expression matrix -----------------------------------------

dir.create(dirname(opt$out_expr), recursive = TRUE, showWarnings = FALSE)
saveRDS(merged, opt$out_expr)
cat("[merge] wrote", opt$out_expr, "\n")

# --- subset metadata --------------------------------------------------------

meta <- read_tsv(opt$joint_meta, show_col_types = FALSE)
required <- c("sample_id", "cancer_type", "sample_type")
miss <- setdiff(required, names(meta))
if (length(miss) > 0) stop("joint_meta missing required columns: ",
                           paste(miss, collapse = ", "))

# Normalise sample type spelling (matches the plotting script's convention).
meta <- meta %>%
  mutate(
    sample_type = case_when(
      sample_type %in% c("tumour", "Tumour", "TUMOUR") ~ "Tumour",
      sample_type %in% c("cell_line", "Cell Line", "CELL_LINE", "cellline") ~ "Cell Line",
      TRUE ~ sample_type
    ),
    cancer_type = toupper(trimws(cancer_type))
  )

# Keep only metadata rows whose sample_id is in the merged matrix and only the
# three cohorts of interest. Order metadata to match merged matrix row order.
keep <- meta$sample_id %in% rownames(merged) &
        meta$cancer_type %in% c("BRCA", "NBL", "RBL")
meta_kept <- meta[keep, , drop = FALSE]
meta_kept <- meta_kept[match(rownames(merged), meta_kept$sample_id), , drop = FALSE]

n_in_meta <- sum(!is.na(meta_kept$sample_id))
n_in_X    <- nrow(merged)
cat(sprintf("[merge] metadata rows kept: %d (out of %d samples in merged matrix)\n",
            n_in_meta, n_in_X))
if (n_in_meta != n_in_X) {
  missing_ids <- rownames(merged)[is.na(meta_kept$sample_id)]
  cat("[WARN] Samples in merged matrix but absent from joint metadata (showing up to 20):\n")
  print(head(missing_ids, 20))
}

# Add display lineage column for downstream convenience (not used during fitting).
meta_kept$display_lineage <- dplyr::recode(
  meta_kept$cancer_type,
  BRCA = "Breast cancer",
  NBL  = "Neuroblastoma",
  RBL  = "Retinoblastoma",
  .default = NA_character_
)

dir.create(dirname(opt$out_meta), recursive = TRUE, showWarnings = FALSE)
write_tsv(meta_kept, opt$out_meta)
cat("[merge] wrote", opt$out_meta, "\n")

# --- merge log (machine-readable) -------------------------------------------

log_path <- opt$out_log %||% paste0(opt$out_meta, ".merge_log.tsv")
log_df <- tibble::tibble(
  step = c(
    "brca_n_samples", "brca_n_genes_in", "brca_n_dup_dropped_keep_first", "brca_orientation",
    "nbl_n_samples",  "nbl_n_genes_in",  "nbl_n_dup_dropped_keep_first",  "nbl_orientation",
    "rbl_n_samples",  "rbl_n_genes_in",  "rbl_n_dup_dropped_keep_first",  "rbl_orientation",
    "n_genes_shared_across_all_three",
    "n_samples_merged_pre_filter",
    "n_genes_merged_pre_filter",
    "n_genes_dropped_na_post_merge",
    "n_genes_dropped_zerovar_post_merge",
    "n_samples_final",
    "n_genes_final",
    "out_expr_rds",
    "out_meta_tsv"
  ),
  value = c(
    as.character(nrow(brca$X)), as.character(ncol(brca$X) + brca$n_dup_dropped_keep_first),
      as.character(brca$n_dup_dropped_keep_first), brca$orientation,
    as.character(nrow(nbl$X)),  as.character(ncol(nbl$X) + nbl$n_dup_dropped_keep_first),
      as.character(nbl$n_dup_dropped_keep_first),  nbl$orientation,
    as.character(nrow(rbl$X)),  as.character(ncol(rbl$X) + rbl$n_dup_dropped_keep_first),
      as.character(rbl$n_dup_dropped_keep_first),  rbl$orientation,
    as.character(length(shared)),
    as.character(nrow(brca_X) + nrow(nbl_X) + nrow(rbl_X)),
    as.character(length(shared)),
    as.character(n_dropped_na),
    as.character(n_dropped_zv),
    as.character(nrow(merged)),
    as.character(ncol(merged)),
    opt$out_expr,
    opt$out_meta
  )
)
write_tsv(log_df, log_path)
cat("[merge] wrote", log_path, "\n")
cat("[OK]    merge complete\n")
