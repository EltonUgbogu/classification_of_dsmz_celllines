#!/usr/bin/env Rscript
# =============================================================================
# compute_nn_purity.R
# =============================================================================
#
# Compute nearest-neighbour purity test for hematologic rejection validation
#
# For each sample, finds k nearest neighbours and computes fraction that are
# from the same lineage. This is the KEY METRIC for hematologic rejection.
#
# USAGE:
# Rscript compute_nn_purity.R \
#   --cor-matrix results/unsupervised/pan_cancer/pan_cancer_cor.rds \
#   --expr-data results/unsupervised/pan_cancer/pan_cancer_expr.rds \
#   --k 10 \
#   --output results/unsupervised/pan_cancer/nn_purity.tsv
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("--cor-matrix"), type="character", default=NULL,
              help="Path to correlation matrix RDS file"),
  make_option(c("--expr-data"), type="character", default=NULL,
              help="Path to pan-cancer expression RDS file (for metadata)"),
  make_option(c("--k"), type="integer", default=10,
              help="Number of nearest neighbours (default: 10)"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output path for NN purity table (TSV)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`cor-matrix`) || is.null(opt$`expr-data`) || is.null(opt$output)) {
  stop("--cor-matrix, --expr-data, and --output are required")
}

cat("Loading correlation matrix...\n")
cor_mat <- readRDS(opt$`cor-matrix`)
cat("  Correlation matrix: ", nrow(cor_mat), " x ", ncol(cor_mat), "\n", sep="")

cat("Loading expression data for metadata...\n")
dat <- readRDS(opt$`expr-data`)
meta <- dat$meta
cat("  Metadata: ", nrow(meta), " samples\n", sep="")

# Create lineage lookup
lineage_lookup <- setNames(meta$lineage, meta$sample_id)

cat("Computing nearest-neighbour purity (k=", opt$k, ")...\n", sep="")
purity_list <- lapply(colnames(cor_mat), function(s) {
  # Get correlations for this sample (excluding self)
  cor_vec <- cor_mat[s, ]
  cor_vec <- cor_vec[names(cor_vec) != s]  # Remove self
  
  # Get k nearest neighbours
  nn_sorted <- sort(cor_vec, decreasing=TRUE)
  nn_ids <- names(nn_sorted)[1:min(opt$k, length(nn_sorted))]
  
  # Get lineages
  sample_lineage <- lineage_lookup[s]
  nn_lineages <- lineage_lookup[nn_ids]
  
  # Compute same-lineage fraction
  same_lineage_frac <- mean(nn_lineages == sample_lineage, na.rm=TRUE)
  heme_frac <- mean(nn_lineages == "HEME", na.rm=TRUE)
  
  data.table(
    sample = s,
    lineage = sample_lineage,
    n_nn = length(nn_ids),
    same_lineage_frac = same_lineage_frac,
    heme_frac = heme_frac,
    nn_lineages = paste(nn_lineages, collapse=";")
  )
})

purity <- rbindlist(purity_list)

cat("Saving NN purity results...\n")
fwrite(purity, file=opt$output, sep="\t")

cat("[OK] NN purity table saved to:", opt$output, "\n")

# Print summary for hematologic samples
cat("\n=== HEMATOLOGIC NN PURITY SUMMARY ===\n")
heme_purity <- purity[lineage == "HEME"]
if (nrow(heme_purity) > 0) {
  cat("Hematologic samples: ", nrow(heme_purity), "\n", sep="")
  cat("Mean same-lineage fraction: ", 
      round(mean(heme_purity$same_lineage_frac, na.rm=TRUE), 3), "\n", sep="")
  cat("Mean heme_frac: ", 
      round(mean(heme_purity$heme_frac, na.rm=TRUE), 3), "\n", sep="")
  cat("Samples with heme_frac >= 0.8: ", 
      sum(heme_purity$heme_frac >= 0.8, na.rm=TRUE), " / ", nrow(heme_purity), "\n", sep="")
  cat("\nPASS CRITERION: heme_frac >= 0.8 for most heme samples\n")
} else {
  cat("No hematologic samples found in data\n")
}

cat("\n=== ALL LINEAGES NN PURITY SUMMARY ===\n")
for (lineage_val in unique(purity$lineage)) {
  sub <- purity[lineage == lineage_val]
  cat(lineage_val, ": mean same-lineage frac = ", 
      round(mean(sub$same_lineage_frac, na.rm=TRUE), 3), 
      " (n=", nrow(sub), ")\n", sep="")
}
