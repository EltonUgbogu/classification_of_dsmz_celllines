#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
stop_with <- function(...) stop(paste0(...), call. = FALSE)

option_list <- list(
  make_option("--direction",  type="character", default=NULL),
  make_option("--hc_rds",     type="character", default=NULL),
  make_option("--kmeans_rds", type="character", default=NULL),
  make_option("--meta_tsv",   type="character", default=NULL),  # not strictly needed here but kept for compatibility
  make_option("--view",       type="character", default=NULL),  # cell, tumour, or cell_tumour
  make_option("--out_root",   type="character", default=NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$direction))  stop_with("Missing --direction")
if (is.null(opt$hc_rds))     stop_with("Missing --hc_rds")
if (is.null(opt$kmeans_rds)) stop_with("Missing --kmeans_rds")
if (is.null(opt$view))        stop_with("Missing --view")
if (is.null(opt$out_root))   stop_with("Missing --out_root")

if (!opt$view %in% c("cell", "tumour", "cell_tumour")) {
  stop_with("Invalid --view: ", opt$view, ". Must be: cell, tumour, or cell_tumour")
}

if (!file.exists(opt$hc_rds))     stop_with("HC RDS not found: ", opt$hc_rds)
if (!file.exists(opt$kmeans_rds)) stop_with("KMeans RDS not found: ", opt$kmeans_rds)

# Helper: extract a named cluster vector from various possible formats
extract_clusters <- function(obj, label) {
  # Common patterns seen in your pipeline:
  # - list(clusters=...)
  # - list(clusters_labeled=...)
  # - list(final_clusters=...)
  # - named vector itself
  cand <- NULL

  if (is.atomic(obj) && !is.null(names(obj))) {
    cand <- obj
  } else if (is.list(obj)) {
    for (nm in c("clusters_labeled","clusters","final_clusters","consensusClass","cluster","labels")) {
      if (!is.null(obj[[nm]])) {
        cand <- obj[[nm]]
        break
      }
    }
  }

  if (is.null(cand)) {
    stop_with("Could not extract clusters from ", label, " object. Keys seen: ",
              if (is.list(obj)) paste(names(obj), collapse=", ") else "<non-list>")
  }
  if (is.null(names(cand))) stop_with("Extracted clusters from ", label, " have no names (sample IDs).")
  cand
}

hc_obj <- readRDS(opt$hc_rds)
km_obj <- readRDS(opt$kmeans_rds)

hc_clusters <- extract_clusters(hc_obj, "HC")
km_clusters <- extract_clusters(km_obj, "KMeans")

# Output paths MUST match Snakefile expectations (use --view to construct paths)
out_hc_dir <- file.path(opt$out_root, paste0("ccp_hc_expr_", opt$view))
out_km_dir <- file.path(opt$out_root, paste0("ccp_kmeans_expr_", opt$view))
dir.create(out_hc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_km_dir, showWarnings = FALSE, recursive = TRUE)

out_hc_rds <- file.path(out_hc_dir, paste0("ccp_hc_expr_", opt$view, "_clusters_optimal.rds"))
out_km_rds <- file.path(out_km_dir, paste0("ccp_kmeans_expr_", opt$view, "_clusters_optimal.rds"))

info("direction=%s, view=%s", opt$direction, opt$view)
info("HC input: %s", opt$hc_rds)
info("KM input: %s", opt$kmeans_rds)
info("Writing:")
info("  %s", out_hc_rds)
info("  %s", out_km_rds)

# Save in a consistent structure (future scripts can rely on this)
saveRDS(list(
  direction = opt$direction,
  alg = "hc",
  clusters = hc_clusters,
  source_rds = opt$hc_rds,
  timestamp = Sys.time()
), out_hc_rds)

saveRDS(list(
  direction = opt$direction,
  alg = "kmeans",
  clusters = km_clusters,
  source_rds = opt$kmeans_rds,
  timestamp = Sys.time()
), out_km_rds)

info("DONE.")

