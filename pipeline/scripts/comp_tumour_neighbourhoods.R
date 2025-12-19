#!/usr/bin/env Rscript

# comp_tumour_neighbourhoods.R
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
  library(optparse)
})

cat("=== Starting computation of tumour neighbourhoods ===\n")

# ─────────────────────────────────────────────────────────────
# 0. CLI: get config path (and optional overrides)
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to YAML config file [default: %default]"),
  make_option("--expr-rds", type = "character", default = NULL,
              help = "Optional override: integrated DSMZ+TCGA expression RDS"),
  make_option("--mapping-rds", type = "character", default = NULL,
              help = "Optional override: cell_line → original_sample_id mapping RDS"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr)")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Using config file:", opt$config, "\n")
cfg <- yaml::read_yaml(opt$config)

# Small helper: x %||% y
`%||%` <- function(x, y) if (!is.null(x)) x else y

# -------------------------------------------------------------------
# 0a. Load base_functions helpers from config
# -------------------------------------------------------------------
if (is.null(cfg$paths$tumour_nh_base_functions)) {
  stop("Config must define paths$tumour_nh_base_functions")
}

base_fun_dir <- cfg$paths$tumour_nh_base_functions

if (!dir.exists(base_fun_dir)) {
  stop("Base-functions directory does not exist: ", base_fun_dir)
}

helper_files <- list.files(
  base_fun_dir,
  pattern = "\\.R$",
  full.names = TRUE
)

if (length(helper_files) == 0L) {
  stop("No .R helper files found in base_functions dir: ", base_fun_dir)
}

cat("[INFO] Loading", length(helper_files), "helper files from:", base_fun_dir, "\n")
invisible(lapply(helper_files, source))

# -------------------------------------------------------------------
# 0b. Validate direction (config-driven)
# -------------------------------------------------------------------
if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Directions are config-driven (BRCA may have PAM50; NBL/RBL won't)
directions <- cfg$tumour_neighbourhoods$directions
if (is.null(directions) || length(directions) == 0) {
  directions <- c("hvg_euc", "hvg_corr")  # fallback
}

if (!direction %in% directions) {
  stop("Invalid --direction: ", direction,
       "\nAllowed (from config tumour_neighbourhoods$directions): ",
       paste(directions, collapse = ", "))
}

# -------------------------------------------------------------------
# 0c. Get paths from config.yaml (brca_analysis paths)
# -------------------------------------------------------------------
unsup_root <- cfg$paths$unsup_root

# Extract gene set from direction
gene_set <- if (startsWith(direction, "pam50")) "PAM50" else "HVG500"

# Extract distance type from direction
dist_type <- if (grepl("_corr$", direction)) "correlation" else "euclidean"

cat("[INFO] Direction: ", direction, "\n", sep = "")
cat("[INFO] Gene set: ", gene_set, "\n", sep = "")
cat("[INFO] Distance metric: ", dist_type, "\n", sep = "")

# Base dirs in the *brca_analysis* tree
base_functions_dir <- cfg$paths$tumour_nh_base_functions
tn_results_root    <- file.path(unsup_root, "tumour_neighbourhoods", direction)

# -------------------------------------------------------------------
# 0b. Paths for expression & mapping (from config.yaml)
# -------------------------------------------------------------------
expr_mat_path <- opt$expr_rds %||%
  if (gene_set == "PAM50") {
    cfg$paths$tumour_nh_expr_pam50
  } else if (gene_set == "HVG500") {
    cfg$paths$tumour_nh_expr_hvg
  } else {
    stop("Unsupported gene_set: ", gene_set, ". Must be PAM50 or HVG500.")
  }

# Direction-aware map file selection (PAM50 vs HVG have separate maps)
mapping_path <- opt$mapping_rds %||% {
  # Infer map file from expression path directory + direction
  map_suffix <- if (startsWith(direction, "pam50")) {
    "cell_line_to_original_sample_id_pam50.rds"
  } else {
    "cell_line_to_original_sample_id_hvg.rds"
  }
  file.path(dirname(expr_mat_path), map_suffix)
}

cat("[INFO] Expression matrix path:\n  ", expr_mat_path, "\n", sep = "")
cat("[INFO] Cell line mapping path:\n  ", mapping_path, "\n", sep = "")

if (!file.exists(expr_mat_path)) {
  stop("Expression matrix RDS not found at: ", expr_mat_path)
}
if (!file.exists(mapping_path)) {
  stop("cell_line_to_original_sample_id RDS not found at: ", mapping_path)
}

# ─────────────────────────────────────────────────────────────
# 1. Source helpers (tumour_nh_io.R + tumour_neighbourhood.R)
# ─────────────────────────────────────────────────────────────
cat("[INFO] Sourcing base functions from:\n  ", base_functions_dir, "\n", sep = "")

for (f in c("tumour_nh_io.R", "tumour_neighbourhood.R")) {
  src_path <- file.path(base_functions_dir, f)
  if (!file.exists(src_path)) {
    stop("Required helper not found: ", src_path)
  }
  source(src_path)   # no local=TRUE so nh_paths & functions are global
}

# Safety checks
if (!exists("make_nh_paths") && !exists("nh_paths")) {
  stop("Neither nh_paths nor make_nh_paths() is defined after sourcing tumour_nh_io.R")
}
if (!exists("compute_tumour_neighbourhoods")) {
  stop("compute_tumour_neighbourhoods() not found after sourcing tumour_neighbourhood.R")
}

# Build nh_paths using unsup_root from config
if (exists("make_nh_paths")) {
  nh_paths <- make_nh_paths(unsup_root)
  cat("[INFO] Built nh_paths using unsup_root:", unsup_root, "\n")
}

# ─────────────────────────────────────────────────────────────
# helper: normalize sample IDs consistently
# ─────────────────────────────────────────────────────────────
normalize_id <- function(x) {
  # strip AGN prefixes (CELL: or TUMOUR:)
  x <- sub("^(CELL:|TUMOUR:)", "", x)
  
  # keep existing normalizations
  x <- gsub("-", "_", x)    # CAL-120  → CAL_120
  x <- gsub(" ", "_", x)    # spaces   → underscores
  x
}

# ─────────────────────────────────────────────────────────────
# helper: extract cell line name from technical DSMZ ID
# e.g. "NG_29643_CAL_120_lib581301_8005_3" → "CAL_120"
# Note: input is expected to be normalized (underscores, not dashes)
# ─────────────────────────────────────────────────────────────
extract_cell_line_from_technical_id <- function(x) {
  # Remove NG_XXXXX_ prefix (normalized version, with underscore)
  x <- gsub("^NG_[0-9]+_", "", x)
  # Remove _libXXXXXX_... suffix (library prep info)
  x <- gsub("_lib[0-9_]+$", "", x)
  x
}

# ─────────────────────────────────────────────────────────────
# 2. Load integrated expression matrix
# ─────────────────────────────────────────────────────────────
expr_mat <- readRDS(expr_mat_path)
cat("[INFO] Expression matrix:", nrow(expr_mat), "samples x", ncol(expr_mat), "genes\n")

# ─────────────────────────────────────────────────────────────
# 3. Load original_sample_id → cell_line mapping
#    Mapping structure (as built by build_dsmz_tcga_pam50matrix.R):
#      names = original sample IDs (NG-29643_CAL_120_lib..., ...)
#      values = cell line names (CAL-120, CAL-51, ...)
# ─────────────────────────────────────────────────────────────
orig_to_cellline <- readRDS(mapping_path)
cat("[INFO] Loaded original_sample_id → cell_line mapping (",
    length(orig_to_cellline), " entries)\n", sep = "")

original_ids_raw    <- names(orig_to_cellline)      # NG-29643_CAL_120_lib..., ...
cell_line_names_raw <- unname(orig_to_cellline)     # CAL-120, CAL-51, ...

cat("[DEBUG] Example mapping keys (original IDs):\n")
print(head(original_ids_raw, 5))
cat("[DEBUG] Example mapping values (cell lines):\n")
print(head(cell_line_names_raw, 5))

current_rownames <- rownames(expr_mat)
cat("[DEBUG] Example expr_mat rownames:\n")
print(head(current_rownames, 10))

# ─────────────────────────────────────────────────────────────
# 4. Identify DSMZ rows by original IDs (BEFORE renaming)
#    Build dataset_vec BEFORE renaming so we don't lose track
# ─────────────────────────────────────────────────────────────
dsmz_mask <- current_rownames %in% original_ids_raw
n_dsmz    <- sum(dsmz_mask)

if (n_dsmz == 0) {
  stop("FATAL: No DSMZ samples found in expr_mat rownames using mapping values.\n",
       "  Example rownames: ", paste(head(current_rownames, 5), collapse = ", "),
       "\n  Example mapping values: ", paste(head(original_ids_raw, 5), collapse = ", "))
}

cat("[INFO] Detected ", n_dsmz, " DSMZ samples in expr_mat via mapping\n", sep = "")

# Build dataset_vec BEFORE renaming (this is the key fix!)
dataset_vec <- ifelse(dsmz_mask, "DSMZ", "TCGA")
names(dataset_vec) <- current_rownames

# ─────────────────────────────────────────────────────────────
# 5. Build mapping: original ID -> normalized cell line name
# ─────────────────────────────────────────────────────────────
cell_line_norm <- normalize_id(cell_line_names_raw)          # CAL-120 -> CAL_120
orig_to_cell   <- setNames(cell_line_norm, original_ids_raw) # NG-... -> CAL_120

# Also keep a normalized version for run_single_neighbourhood()
orig_ids_norm <- normalize_id(original_ids_raw)
orig_to_cell_norm <- setNames(cell_line_norm, orig_ids_norm)

# ─────────────────────────────────────────────────────────────
# 6. Rename DSMZ rows from original IDs -> clean cell line names
# ─────────────────────────────────────────────────────────────
ids_to_replace <- current_rownames[dsmz_mask]
current_rownames[dsmz_mask] <- orig_to_cell[ids_to_replace]

rownames(expr_mat) <- current_rownames

# Final normalization + sync dataset_vec names
rownames(expr_mat) <- normalize_id(rownames(expr_mat))
names(dataset_vec) <- rownames(expr_mat)

cat("[INFO] Successfully renamed ", n_dsmz, " DSMZ samples in expr_mat.\n", sep = "")
cat("Verification of DSMZ sample names (first 5):\n")
print(head(rownames(expr_mat)[dataset_vec == "DSMZ"], 5))

cat("[INFO] Total samples: ", nrow(expr_mat),
    " ( ", sum(dataset_vec == "DSMZ"), " DSMZ )\n", sep = "")

# ─────────────────────────────────────────────────────────────
# 5. Method table: pipeline-native discovery
#    Uses get_nh_methods() to enumerate actual pipeline outputs
# ─────────────────────────────────────────────────────────────
# after sourcing tumour_nh_io.R and tumour_neighbourhood.R
methods <- get_nh_methods(unsup_root = unsup_root, direction = direction)

cat("\n[INFO] Methods discovered:\n")
print(methods)

methods_exist <- methods %>% dplyr::filter(exists)
if (nrow(methods_exist) == 0) {
  stop("FATAL: No cluster RDS found for direction=", direction,
       "\nChecked:\n", paste(methods$path, collapse = "\n"))
}

cat(sprintf("\n[INFO] %d/%d methods available (existing files)\n",
            nrow(methods_exist), nrow(methods)))

# ─────────────────────────────────────────────────────────────
# 6. Helper to run a single method
# ─────────────────────────────────────────────────────────────
run_single_neighbourhood <- function(path, method_id, outdir) {
  cat("\n============================================================\n")
  cat("Method:", method_id, "\n")
  cat("Loading clusters from:", path, "\n")

  # Guard: skip if cluster file doesn't exist
  if (!file.exists(path)) {
    message("[WARN] Cluster file not found for method ", method_id, ": ", path,
            " — skipping.")
    return(NULL)
  }

  clust_obj <- readRDS(path)
  if (is.null(clust_obj$clusters)) {
    stop("Object at ", path, " does not have a $clusters element.")
  }

  cluster_vec <- clust_obj$clusters
  if (is.null(names(cluster_vec))) {
    stop("Cluster vector in ", path, " has no names.")
  }

  # Normalize cluster IDs to same convention as expr_mat
  names(cluster_vec) <- normalize_id(names(cluster_vec))

  # Map any original NG ids → cell-line names using *normalized* map
  nm  <- names(cluster_vec)
  hit <- nm %in% names(orig_to_cell_norm)
  if (any(hit)) {
    nm[hit] <- orig_to_cell_norm[nm[hit]]
    names(cluster_vec) <- nm
    cat("[INFO] Mapped", sum(hit),
        "cluster IDs from NG_* to cell-line names for", method_id, "\n")
  }
  
  # Fallback: if still looks like NG_* technical id, extract cell line
  # (handles cases where mapping might be incomplete)
  is_ng <- grepl("^NG_[0-9]+_", nm)
  if (any(is_ng)) {
    nm[is_ng] <- extract_cell_line_from_technical_id(nm[is_ng])
    names(cluster_vec) <- nm
    cat("[INFO] Fallback: extracted cell-line names from", sum(is_ng),
        "NG_* technical IDs for", method_id, "\n")
  }

  # Quick peek at DSMZ names in cluster_vec
  dsmz_names_in_clusters <- names(cluster_vec)[!grepl("^TCGA-", names(cluster_vec))]
  cat("[INFO] First DSMZ names in cluster_vec: ",
      paste(head(dsmz_names_in_clusters, 5), collapse = ", "), "\n")

  # Align to expr_mat
  common_ids <- intersect(rownames(expr_mat), names(cluster_vec))

  if (length(common_ids) == 0) {
    stop("No overlapping sample IDs between expr_mat and cluster_vec for ", method_id,
         "\n  expr_mat sample IDs (first 5): ", paste(head(rownames(expr_mat), 5), collapse = ", "),
         "\n  cluster_vec sample IDs (first 5): ", paste(head(names(cluster_vec), 5), collapse = ", "))
  }

  ## A) Samples in expr_mat but not in cluster_vec
  missing_in_clusters <- setdiff(rownames(expr_mat), names(cluster_vec))
  if (length(missing_in_clusters) > 0) {
    cat("[WARNING] ", length(missing_in_clusters),
        " samples in expr_mat not found in cluster_vec for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_clusters, 5), collapse = ", "), "\n", sep = "")
  }

  ## B) Samples in cluster_vec but not in expr_mat
  missing_in_expr <- setdiff(names(cluster_vec), rownames(expr_mat))
  if (length(missing_in_expr) > 0) {
    dsmz_only <- missing_in_expr[!grepl("^TCGA-", missing_in_expr)]
    cat("[WARNING] ", length(missing_in_expr),
        " samples in cluster_vec not found in expr_mat for ", method_id, "\n", sep = "")
    cat("  Examples: ",
        paste(head(missing_in_expr, 5), collapse = ", "), "\n", sep = "")
    if (length(dsmz_only) > 0) {
      cat("  DSMZ-only in cluster_vec (not in expr_mat): ",
          paste(head(dsmz_only, 10), collapse = ", "), "\n", sep = "")
    }
  }

  ## C) Subset to common_ids
  if (length(common_ids) < nrow(expr_mat)) {
    cat("  Using only ", length(common_ids), " common samples\n", sep = "")
    expr_mat_subset    <- expr_mat[common_ids, , drop = FALSE]
    dataset_vec_subset <- dataset_vec[common_ids]
  } else {
    expr_mat_subset    <- expr_mat
    dataset_vec_subset <- dataset_vec
  }

  cluster_vec <- cluster_vec[rownames(expr_mat_subset)]
  if (any(is.na(cluster_vec))) {
    missing_ids <- rownames(expr_mat_subset)[is.na(cluster_vec)]
    stop("Internal error: Some samples still missing after subsetting for ", method_id,
         ". Examples: ", paste(head(missing_ids, 5), collapse = ", "))
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # === Main call – correct argument name is cluster_m ===
  # Pass distance type based on direction (euclidean vs correlation)
  res <- compute_tumour_neighbourhoods(
    emb_mat   = expr_mat_subset,
    cluster_m = cluster_vec,
    dataset   = dataset_vec_subset,
    top_frac  = 0.10,
    top_n_min = 30,
    top_n_max = 200,
    method_id = method_id,
    distance  = dist_type  # NEW: pass distance metric from direction
  )

  # Save outputs
  nh_rds   <- file.path(outdir, paste0("Top_m_neighbourhoods_", method_id, ".rds"))
  long_rds <- file.path(outdir, paste0("Top_m_long_", method_id, ".rds"))
  long_csv <- file.path(outdir, paste0("Top_m_long_", method_id, ".csv"))

  saveRDS(res$neighbourhoods, nh_rds)
  saveRDS(res$long_df,       long_rds)
  write_csv(res$long_df,     long_csv)

  cat("\nTumour neighbourhoods computed successfully!\n")
  cat("Method      :", res$method_id, "\n")
  cat("Cell lines  :", length(res$neighbourhoods), "\n")
  cat("Total pairs :", nrow(res$long_df), "\n\n")

  print(
    res$long_df %>%
      count(cell_line, in_top) %>%
      arrange(desc(in_top), cell_line)
  )

  cat("\nAll results saved to:", outdir, "\n")
  invisible(res)
}

# ─────────────────────────────────────────────────────────────
# 7. Run all methods
# ─────────────────────────────────────────────────────────────
all_results <- vector("list", nrow(methods_exist))
names(all_results) <- methods_exist$method_id

n_success <- 0
n_skipped <- 0

for (i in seq_len(nrow(methods_exist))) {
  m <- methods_exist[i, ]
  result <- run_single_neighbourhood(
    path      = m$path,
    method_id = m$method_id,
    outdir    = m$outdir
  )
  
  if (is.null(result)) {
    n_skipped <- n_skipped + 1
  } else {
    all_results[[m$method_id]] <- result
    n_success <- n_success + 1
  }
}

# Remove NULL entries
all_results <- Filter(Negate(is.null), all_results)

cat("\n============================================================\n")
cat("SUMMARY:\n")
cat(sprintf("  Methods attempted: %d\n", nrow(methods_exist)))
cat(sprintf("  Successful: %d\n", n_success))
cat(sprintf("  Skipped (missing files): %d\n", n_skipped))

if (n_success == 0) {
  stop("FATAL: No clustering methods completed successfully. Check cluster file paths.")
}

cat("\nAll available methods completed.\n")

