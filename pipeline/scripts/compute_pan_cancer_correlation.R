#!/usr/bin/env Rscript
# =============================================================================
# compute_pan_cancer_correlation.R
# =============================================================================
#
# Compute Spearman correlation matrix for pan-cancer expression data
#
# USAGE:
# Rscript compute_pan_cancer_correlation.R \
#   --input results/unsupervised/pan_cancer/pan_cancer_expr.rds \
#   --output results/unsupervised/pan_cancer/pan_cancer_cor.rds
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("--input"), type="character", default=NULL,
              help="Path to pan-cancer expression RDS file"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output path for correlation matrix (RDS)"),
  make_option(c("--method"), type="character", default="spearman",
              help="Correlation method: spearman or pearson (default: spearman)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  stop("--input and --output are required")
}

cat("Loading pan-cancer expression matrix...\n")
dat <- readRDS(opt$input)
expr <- dat$expr

cat("  Expression matrix: ", nrow(expr), " genes x ", ncol(expr), " samples\n", sep="")

cat("Computing", opt$method, "correlation matrix...\n")
cor_mat <- cor(expr, method=opt$method)

cat("  Correlation matrix: ", nrow(cor_mat), " x ", ncol(cor_mat), "\n", sep="")

cat("Saving correlation matrix...\n")
saveRDS(cor_mat, file=opt$output)

cat("[OK] Correlation matrix saved to:", opt$output, "\n")
cat("  Correlation range: [", min(cor_mat, na.rm=TRUE), ", ", 
    max(cor_mat, na.rm=TRUE), "]\n", sep="")
