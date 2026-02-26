#!/usr/bin/env Rscript
# =============================================================================
# quantify_cross_lineage_edges.R
# =============================================================================
#
# Quantify cross-lineage edges (HEME ↔ solid cancers)
# This is the HARD PROOF for hematologic rejection validation
#
# PASS CRITERION: nrow(cross) == 0 (or extremely close to 0)
#
# USAGE:
# Rscript quantify_cross_lineage_edges.R \
#   --edges results/unsupervised/pan_cancer/graph/pan_cancer_graph_edges.tsv \
#   --expr-data results/unsupervised/pan_cancer/pan_cancer_expr.rds \
#   --output results/unsupervised/pan_cancer/cross_lineage_edges.tsv
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("--edges"), type="character", default=NULL,
              help="Path to graph edges TSV file"),
  make_option(c("--expr-data"), type="character", default=NULL,
              help="Path to pan-cancer expression RDS file (for metadata)"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output path for cross-lineage edges (TSV)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$edges) || is.null(opt$`expr-data`) || is.null(opt$output)) {
  stop("--edges, --expr-data, and --output are required")
}

cat("Loading graph edges...\n")
edges <- fread(opt$edges)
cat("  Edges: ", nrow(edges), "\n", sep="")

cat("Loading expression data for metadata...\n")
dat <- readRDS(opt$`expr-data`)
meta <- dat$meta
cat("  Metadata: ", nrow(meta), " samples\n", sep="")

# Create lineage lookup
lineage_lookup <- setNames(meta$lineage, meta$sample_id)

# Add lineage information to edges
cat("Annotating edges with lineage information...\n")
edges[, from_lineage := lineage_lookup[from]]
edges[, to_lineage := lineage_lookup[to]]

# Identify cross-lineage edges (HEME ↔ solid)
cat("Identifying cross-lineage edges (HEME ↔ solid)...\n")
cross <- edges[
  (from_lineage == "HEME" & to_lineage != "HEME") |
  (from_lineage != "HEME" & to_lineage == "HEME")
]

cat("  Cross-lineage edges (HEME ↔ solid): ", nrow(cross), "\n", sep="")

# Also identify within-lineage edges for comparison
within_heme <- edges[from_lineage == "HEME" & to_lineage == "HEME"]
within_solid <- edges[
  from_lineage %in% c("BRCA", "NBL", "RBL") & 
  to_lineage %in% c("BRCA", "NBL", "RBL")
]

cat("  Within-HEME edges: ", nrow(within_heme), "\n", sep="")
cat("  Within-solid edges: ", nrow(within_solid), "\n", sep="")

# Save cross-lineage edges
if (nrow(cross) > 0) {
  fwrite(cross, file=opt$output, sep="\t")
  cat("  Cross-lineage edges saved to:", opt$output, "\n")
} else {
  # Create empty file with header
  fwrite(cross, file=opt$output, sep="\t")
  cat("  No cross-lineage edges found (PASS!)\n")
}

# Print summary
cat("\n=== CROSS-LINEAGE EDGE SUMMARY ===\n")
cat("Total edges: ", nrow(edges), "\n", sep="")
cat("Cross-lineage (HEME ↔ solid): ", nrow(cross), "\n", sep="")
cat("Within-HEME: ", nrow(within_heme), "\n", sep="")
cat("Within-solid: ", nrow(within_solid), "\n", sep="")

if (nrow(cross) == 0) {
  cat("\n✅ PASS: No cross-lineage edges (HEME ↔ solid)\n")
} else {
  cat("\n⚠️  WARNING: ", nrow(cross), " cross-lineage edges found\n", sep="")
  cat("   Cross-lineage fraction: ", 
      round(nrow(cross) / nrow(edges) * 100, 2), "%\n", sep="")
  
  # Show breakdown by lineage pair
  cat("\n   Breakdown by lineage pair:\n")
  cross[, pair := paste(pmin(from_lineage, to_lineage), 
                        pmax(from_lineage, to_lineage), sep=" ↔ ")]
  print(table(cross$pair))
}

# Also check component composition
cat("\n=== COMPONENT COMPOSITION ===\n")
# This would require loading components file, but we can compute it here
# For now, just report the edge statistics

cat("[OK] Cross-lineage edge analysis complete\n")
