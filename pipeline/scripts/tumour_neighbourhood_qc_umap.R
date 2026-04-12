#!/usr/bin/env Rscript

# tumour_neighbourhood_qc_umap.R
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
})

log_info <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))
log_warn <- function(...) cat(sprintf("[WARN] %s\n", sprintf(...)))
log_stop <- function(...) stop(sprintf("[ERROR] %s", sprintf(...)), call. = FALSE)

# -----------------------------
# CLI
# -----------------------------
opt_list <- list(
  make_option(c("-c", "--config"), type="character", help="Path to config.yaml"),
  make_option("--profile", type="character", default=NULL,
              help="Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  make_option("--direction",    type="character", help="Direction, e.g. HVG_euc, HVG_corr, pam50_euc, pam50_corr"),
  make_option("--consensus_rds",type="character", help="Final_consensus_tumour_neighbourhoods_*.rds"),
  make_option("--expr_rds",     type="character", help="Integrated expression matrix RDS (samples x genes)"),
  make_option("--map_rds",      type="character", help="cell_line_to_original_sample_id_*.rds"),
  make_option("--qc_tsv",       type="character", help="Output QC summary TSV"),
  make_option("--umap_pdf",     type="character", help="Output UMAP PDF"),
  make_option("--umap_tsv",     type="character", help="Output UMAP coordinates TSV")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Source lib_config.R using absolute path (works with --directory)
# Derive pipeline root from config path
if (!is.null(opt$config) && file.exists(opt$config)) {
  config_dir <- dirname(normalizePath(opt$config))
  pipe_root <- normalizePath(file.path(config_dir, ".."))
  lib_config_path <- file.path(pipe_root, "scripts", "lib_config.R")
} else {
  # Fallback: try relative path (for manual runs)
  lib_config_path <- "scripts/lib_config.R"
}
if (file.exists(lib_config_path)) {
  source(lib_config_path)
} else {
  stop("Cannot find lib_config.R. Tried: ", lib_config_path)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

req <- c("config","direction","consensus_rds","expr_rds","map_rds","qc_tsv","umap_pdf","umap_tsv")
missing <- req[!nzchar(sapply(req, function(k) opt[[k]] %||% ""))]
if (length(missing) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse=", ")), call. = FALSE)
}

# -----------------------------
# IO helpers
# -----------------------------
dir_create <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_rds_or_stop <- function(p) {
  if (!file.exists(p)) log_stop("File not found: %s", p)
  readRDS(p)
}

# Try to extract a character vector of "cell sample IDs" from map_rds,
# regardless of whether it is a named vector, list, or data frame.
# In your pipeline, map_rds is typically:
# names(map)  = original technical IDs (NG-29643_...)
# values(map) = cell line names (CAL-120, ...)
# After comp_tumour_neighbourhoods.R, expr rownames become CAL_120 etc.
extract_cell_ids <- function(map_obj) {
  normalize_id <- function(x) {
    x <- sub("^(CELL:|TUMOUR:)", "", x)
    x <- gsub("-", "_", x)
    x <- gsub(" ", "_", x)
    x
  }

  # Named vector: use VALUES as the real cell IDs (not names)
  if (is.atomic(map_obj) && !is.null(names(map_obj))) {
    return(unique(normalize_id(as.character(unname(map_obj)))))
  }

  # Character vector
  if (is.atomic(map_obj) && is.character(map_obj)) {
    return(unique(normalize_id(as.character(map_obj))))
  }

  # List
  if (is.list(map_obj) && !is.data.frame(map_obj)) {
    if (!is.null(names(map_obj))) {
      # if it's a named list, names often still not the right thing; prefer values if possible
      vals <- unlist(map_obj, use.names = FALSE)
      if (length(vals) > 0) return(unique(normalize_id(as.character(vals))))
      return(unique(normalize_id(as.character(names(map_obj)))))
    }
    return(unique(normalize_id(as.character(unlist(map_obj, use.names = FALSE)))))
  }

  # Data frame
  if (is.data.frame(map_obj)) {
    cn <- tolower(colnames(map_obj))
    # prefer cell_line-ish column
    candidates <- c("cell_line","cellline","dsmz","cell","cell_id")
    hit <- intersect(candidates, cn)
    if (length(hit) > 0) {
      col <- colnames(map_obj)[match(hit[1], cn)]
      return(unique(normalize_id(as.character(map_obj[[col]]))))
    }
    return(unique(normalize_id(as.character(map_obj[[1]]))))
  }

  character(0)
}

# Ensure matrix is numeric (samples x genes)
as_numeric_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) log_stop("Expression object is not a matrix/data.frame.")
  storage.mode(x) <- "numeric"
  x
}

# -----------------------------
# Load config (optional usage)
# -----------------------------
# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)

log_info("Direction: %s", opt$direction)
log_info("expr_rds: %s", opt$expr_rds)
log_info("map_rds:  %s", opt$map_rds)
log_info("consensus_rds: %s", opt$consensus_rds)

expr <- read_rds_or_stop(opt$expr_rds)
expr <- as_numeric_matrix(expr)

if (is.null(rownames(expr))) {
  log_stop("Expression matrix has no rownames; expected sample IDs as rownames.")
}
if (nrow(expr) < 5 || ncol(expr) < 5) {
  log_stop("Expression matrix too small (nrow=%d, ncol=%d).", nrow(expr), ncol(expr))
}

map_obj <- read_rds_or_stop(opt$map_rds)
cell_ids <- extract_cell_ids(map_obj)

if (length(cell_ids) == 0) {
  log_warn("Could not extract cell IDs from map_rds. All samples will be labelled 'UNKNOWN'.")
}

sample_ids <- rownames(expr)
sample_type <- if (length(cell_ids) == 0) {
  rep("UNKNOWN", length(sample_ids))
} else {
  ifelse(sample_ids %in% cell_ids, "CELL", "TUMOUR")
}

# -----------------------------
# Load consensus and summarize (NEW: long-table consensus)
# -----------------------------
cons <- read_rds_or_stop(opt$consensus_rds)

# Expect the new format: a data.frame with columns:
# cell_line, tumour_id, p_consensus  (plus n_methods etc.)
if (!is.data.frame(cons)) {
  log_stop("consensus_rds is not a data.frame. Expected long table with p_consensus.")
}

if (!("cell_line" %in% colnames(cons))) {
  if ("cell_line_canonical" %in% colnames(cons)) {
    cons$cell_line <- cons$cell_line_canonical
  } else if ("cell_line_display" %in% colnames(cons)) {
    cons$cell_line <- cons$cell_line_display
  } else if ("cell_tech_id" %in% colnames(cons)) {
    cons$cell_line <- cons$cell_tech_id
  }
}

if (!"tumour_id" %in% colnames(cons)) {
  if ("tumor_id" %in% colnames(cons)) {
    cons <- cons %>% dplyr::rename(tumour_id = tumor_id)
  } else {
    log_stop("consensus_rds missing required column: tumour_id")
  }
}

need_cols <- c("cell_line","tumour_id","p_consensus")
miss_cols <- setdiff(need_cols, colnames(cons))
if (length(miss_cols) > 0) {
  log_stop("consensus_rds missing required columns: %s", paste(miss_cols, collapse=", "))
}

# Clean p_consensus
cons2 <- cons %>%
  dplyr::filter(!is.na(p_consensus), is.finite(p_consensus), p_consensus >= 0, p_consensus <= 1)

# Per-cell-line summaries
cell_sum <- cons2 %>%
  group_by(cell_line) %>%
  summarise(
    n_pairs = n(),
    max_p = max(p_consensus),
    median_p = median(p_consensus),
    frac_ge_0_7 = mean(p_consensus >= 0.7),
    .groups = "drop"
  )

# Per-tumour summaries (optional QC: how "claimed" are tumours)
tum_sum <- cons2 %>%
  group_by(tumour_id) %>%
  summarise(
    best_p = max(p_consensus),
    n_cell_lines_ge_0_7 = sum(p_consensus >= 0.7),
    .groups = "drop"
  )

pcons_summary <- tibble(
  direction = opt$direction,
  n_samples = nrow(expr),
  n_genes   = ncol(expr),
  n_cell_samples   = sum(sample_type == "CELL"),
  n_tumour_samples = sum(sample_type == "TUMOUR"),
  frac_missing_expr = mean(is.na(expr)),
  pcons_n_celllines = nrow(cell_sum),
  pcons_pairs_total = nrow(cons2),
  pcons_cell_max_p_median = median(cell_sum$max_p, na.rm = TRUE),
  pcons_cell_max_p_mean   = mean(cell_sum$max_p, na.rm = TRUE),
  pcons_tum_best_p_median = median(tum_sum$best_p, na.rm = TRUE),
  pcons_tum_best_p_mean   = mean(tum_sum$best_p, na.rm = TRUE)
)

# Keep for plotting later if you want
# cell_sum, tum_sum, cons2 are available

# Gene variance summary (quick)
gene_vars <- apply(expr, 2, var, na.rm = TRUE)
pcons_summary <- pcons_summary %>%
  mutate(
    gene_var_median = median(gene_vars, na.rm = TRUE),
    gene_var_mean   = mean(gene_vars, na.rm = TRUE)
  )

# -----------------------------
# UMAP computation
# -----------------------------
set.seed(42)

is_corr <- grepl("_corr$", opt$direction)

# Work on a clean matrix (remove all-NA genes)
keep_genes <- which(colSums(!is.na(expr)) > 0)
expr2 <- expr[, keep_genes, drop = FALSE]

# Impute remaining NA with gene means (minimal, stable)
if (anyNA(expr2)) {
  gene_means <- colMeans(expr2, na.rm = TRUE)
  for (j in seq_len(ncol(expr2))) {
    idx <- which(is.na(expr2[, j]))
    if (length(idx) > 0) expr2[idx, j] <- gene_means[j]
  }
}

umap_res <- NULL

if (is_corr) {
  log_info("UMAP mode: correlation distance (1 - Pearson correlation between samples)")
  
  cor_mat <- suppressWarnings(cor(t(expr2), use = "pairwise.complete.obs", method = "pearson"))
  cor_mat[!is.finite(cor_mat)] <- 0
  cor_mat <- pmax(pmin(cor_mat, 1), -1)
  diag(cor_mat) <- 1
  
  dist_mat <- 1 - cor_mat
  dist_mat[!is.finite(dist_mat)] <- 1
  dist_mat <- pmax(dist_mat, 0)
  diag(dist_mat) <- 0
  
  dist_obj <- as.dist(dist_mat)
  
  umap_res <- uwot::umap(
    dist_obj,
    n_neighbors = 15,
    min_dist = 0.3,
    n_components = 2,
    verbose = TRUE
  )
} else {
  log_info("UMAP mode: euclidean (PCA -> UMAP)")
  # scale genes, PCA, keep min(50, ncol-1)
  x <- scale(expr2)
  npc <- min(50, ncol(x) - 1, nrow(x) - 1)
  if (npc < 2) log_stop("Not enough dimensions for PCA/UMAP after filtering.")
  pca <- prcomp(x, center = FALSE, scale. = FALSE)
  Xp <- pca$x[, seq_len(npc), drop = FALSE]

  umap_res <- uwot::umap(
    Xp,
    metric = "euclidean",
    n_neighbors = 15,
    min_dist = 0.3,
    n_components = 2,
    verbose = TRUE
  )
}

umap_df <- tibble(
  sample_id = sample_ids,
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2],
  sample_type = sample_type
)

# -----------------------------
# Write outputs
# -----------------------------
dir_create(opt$qc_tsv)
dir_create(opt$umap_tsv)
dir_create(opt$umap_pdf)

readr::write_tsv(pcons_summary, opt$qc_tsv)
readr::write_tsv(umap_df, opt$umap_tsv)

# Plot
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = sample_type, colour = sample_type)) +
  geom_point(size = 1.3, alpha = 0.8) +
  theme_classic(base_size = 11) +
  labs(
    title = sprintf("Tumour-Cell Integration UMAP (%s)", opt$direction),
    subtitle = "Colour/shape indicate sample type inferred from map_rds",
    x = "UMAP1", y = "UMAP2"
  )

ggsave(opt$umap_pdf, p, width = 8, height = 6, units = "in")

log_info("Wrote QC TSV: %s", opt$qc_tsv)
log_info("Wrote UMAP TSV: %s", opt$umap_tsv)
log_info("Wrote UMAP PDF: %s", opt$umap_pdf)
log_info("Done.")

