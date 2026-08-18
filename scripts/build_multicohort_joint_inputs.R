#!/usr/bin/env Rscript

# build_multicohort_joint_inputs.R

# ---- Thread limits (HPC-safe) ----
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
if (is.na(n_threads) || n_threads < 1) n_threads <- 8

Sys.setenv(
  OMP_NUM_THREADS = n_threads,
  OPENBLAS_NUM_THREADS = n_threads,
  MKL_NUM_THREADS = n_threads,
  VECLIB_MAXIMUM_THREADS = n_threads,
  NUMEXPR_NUM_THREADS = n_threads
)

message(sprintf("[INFO] Thread caps set to %d (SLURM_CPUS_PER_TASK=%s)",
                n_threads, Sys.getenv("SLURM_CPUS_PER_TASK", unset = "NA")))

# Additional thread caps for R-specific parallel backends
options(mc.cores = n_threads)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = n_threads)

#
# Joint tumour vs joint cell line multicohort benchmarking
# - Build ONE joint expression matrix from multiple "joint VST" disease matrices
# - Run UMAP (no PCA / no cPCA)
# - Benchmark feature selection methods × UMAP distance metrics
#
# INPUTS (via CLI):
#   --profiles: comma-separated profile names (e.g., "brca,nbl,rbl")
#   --vst_list: space-separated paths to joint VST RDS files (one per profile, in same order)
#
# OUTPUT:
#   results/multicohort_cancer_benchmark/
#     inputs/joint_expr_matrix.rds
#     inputs/joint_metadata.tsv
#     benchmarks/summary.tsv
#     benchmarks/coords_<feature>_<dist>.tsv
#     benchmarks/01_main_overview_<feature>_<dist>.pdf/png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(uwot)
  library(optparse)
  library(scales)
  library(matrixStats)
  library(RColorBrewer)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--pipe_root", type="character", default=".",
              help="Pipeline root [default: %default]"),
  make_option("--outdir", type="character", default="results/multicohort_cancer_benchmark",
              help="Output directory [default: %default]"),

  # Profiles and their VST paths (from Snakemake)
  make_option("--profiles", type="character", default="brca,nbl,rbl,heme",
              help="Comma-separated profile names [default: %default]"),
  make_option("--vst_list", type="character", default="",
              help="Space-separated paths to joint VST RDS files (one per profile, in same order as --profiles)"),

  # Legacy: individual paths (for manual runs)
  make_option("--brca_joint", type="character", default="",
              help="BRCA joint VST RDS [default: auto-detect from --vst_list]"),
  make_option("--nbl_joint", type="character", default="",
              help="NBL joint VST RDS [default: auto-detect from --vst_list]"),
  make_option("--rbl_joint", type="character", default="",
              help="RBL joint VST RDS [default: auto-detect from --vst_list]"),

  # Benchmark grid
  make_option("--feature_methods", type="character",
              default="Variance,MAD,Entropy",
              help="Comma-separated feature methods [default: %default]"),
  make_option("--dist_metrics", type="character",
              default="cosine,euclidean",
              help="Comma-separated UMAP metrics (uwot) [default: %default]"),
  make_option("--n_genes", type="integer", default=500,
              help="Number of selected genes [default: %default]"),

  # UMAP params
  make_option("--n_neighbors", type="integer", default=30,
              help="UMAP n_neighbors [default: %default]"),
  make_option("--min_dist", type="double", default=0.3,
              help="UMAP min_dist [default: %default]"),
  make_option("--seed", type="integer", default=42,
              help="Random seed [default: %default]"),

  # Evaluation params
  make_option("--eval_k", type="integer", default=30,
              help="k for kNN evaluation in UMAP space [default: %default]"),

  # Optional: downsample for speed in very large cohorts (0 = no downsample)
  make_option("--max_samples", type="integer", default=0,
              help="If >0, downsample to this many total samples (stratified by cancer_type×sample_type). [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

pipe_root <- normalizePath(opt$pipe_root, mustWork = TRUE)
abs_from_root <- function(p) {
  if (is.null(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(pipe_root, p)
}

outdir <- abs_from_root(opt$outdir)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "inputs"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "benchmarks"), recursive=TRUE, showWarnings=FALSE)

feature_methods <- trimws(strsplit(opt$feature_methods, ",")[[1]])
dist_metrics    <- trimws(strsplit(opt$dist_metrics, ",")[[1]])

# Parse profiles and VST paths
profiles <- strsplit(opt$profiles, ",")[[1]] %>% trimws()

# Map profile names to cancer type labels
profile_to_cancer <- function(prof) {
  prof_upper <- toupper(prof)
  case_when(
    prof_upper == "BRCA" ~ "BRCA",
    prof_upper == "NBL"  ~ "NBL",
    prof_upper == "RBL"  ~ "RBL",
    prof_upper %in% c("HEME", "HEM", "HEMATOLOGIC", "HEMATOLOGY",
                      "LEU", "LYM", "LEU_LYM", "LEU/LYM") ~ "HEME",
    TRUE ~ prof_upper
  )
}

cancer_label_map <- c(
  BRCA = "Breast Cancer",
  NBL  = "Neuroblastoma",
  RBL  = "Retinoblastoma",
  HEME = "Hematologic Malignancy"
)

# Parse VST paths
vst_paths <- character(0)
if (nzchar(opt$vst_list)) {
  # From Snakemake: space-separated list
  vst_paths <- strsplit(opt$vst_list, " +")[[1]] %>% trimws()
  vst_paths <- vst_paths[nzchar(vst_paths)]
} else {
  # Fallback: try legacy individual options
  if (nzchar(opt$brca_joint)) vst_paths <- c(vst_paths, opt$brca_joint)
  if (nzchar(opt$nbl_joint))  vst_paths <- c(vst_paths, opt$nbl_joint)
  if (nzchar(opt$rbl_joint))  vst_paths <- c(vst_paths, opt$rbl_joint)
}

# Validate: profiles and vst_paths must match in length
if (length(vst_paths) == 0) {
  stop("No VST paths provided. Use --vst_list or legacy --brca_joint/--nbl_joint/--rbl_joint")
}
if (length(profiles) != length(vst_paths)) {
  stop(sprintf("Mismatch: %d profiles but %d VST paths. Profiles: %s", 
               length(profiles), length(vst_paths), paste(profiles, collapse=", ")))
}

cat("[INFO] pipe_root:", pipe_root, "\n")
cat("[INFO] outdir   :", outdir, "\n")
cat("[INFO] profiles :", paste(profiles, collapse=", "), "\n")
cat("[INFO] feature_methods:", paste(feature_methods, collapse=", "), "\n")
cat("[INFO] dist_metrics   :", paste(dist_metrics, collapse=", "), "\n")

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

clean_ensg <- function(ids) sub("\\.\\d+$", "", ids)

as_samples_x_genes <- function(x) {
  # Accept matrix / data.frame; return samples x genes
  X <- as.matrix(x)

  rn <- rownames(X) %||% character(0)
  cn <- colnames(X) %||% character(0)

  # Heuristic:
  # - sample IDs often include TCGA-/TARGET/GSE or DSMZ-like names (letters/numbers/underscore)
  # - gene IDs often start with ENSG
  rn_has_ensg <- length(rn) > 0 && mean(grepl("^ENSG", rn)) > 0.5
  cn_has_ensg <- length(cn) > 0 && mean(grepl("^ENSG", cn)) > 0.5

  if (rn_has_ensg && !cn_has_ensg) {
    # genes x samples -> transpose
    X <- t(X)
  }
  # Now we *expect* columns are genes
  if (ncol(X) < 50) {
    warning("Matrix has <50 columns after orientation; check input RDS (ncol=", ncol(X), ").")
  }
  if (!is.null(colnames(X))) colnames(X) <- clean_ensg(colnames(X))
  X
}

infer_sample_type <- function(sample_id) {
  ifelse(grepl("^TCGA-|^TARGET|^GSE|^SRP[0-9]+|^SRR[0-9]+|^GDC-", sample_id), "Tumour", "Cell Line")
}

theme_umap_main <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(linewidth = 0.4, colour = "grey20"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text         = element_text(size = 9, colour = "grey20"),
      axis.title        = element_text(size = 10, colour = "grey10"),
      legend.position   = "right",
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9),
      legend.key.size   = grid::unit(0.8, "lines"),
      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
      plot.title        = element_text(size = 13, face = "bold", colour = "grey10"),
      plot.subtitle     = element_text(size = 10, colour = "grey40"),
      plot.caption      = element_text(size = 7, colour = "grey50", hjust = 0),
      plot.margin       = margin(12, 12, 8, 8),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

prepare_plot_data <- function(df) {
  disease <- cancer_label_map[df$cancer_type]
  disease[is.na(disease)] <- df$cancer_type[is.na(disease)]
  disease_levels <- unique(disease)
  df %>%
    mutate(
      type = factor(ifelse(sample_type == "Cell Line", "Cell Line", "Tumour"),
                    levels = c("Tumour", "Cell Line")),
      disease = factor(disease, levels = disease_levels)
    )
}

build_disease_palette <- function(labels) {
  base_palette <- c(
    "Breast Cancer"         = "#E41A1C",
    "Neuroblastoma"         = "#377EB8",
    "Retinoblastoma"        = "#FF7F00",
    "Hematologic Malignancy" = "#4DAF4A"
  )
  palette <- setNames(rep(NA_character_, length(labels)), labels)
  overlap <- intersect(labels, names(base_palette))
  palette[overlap] <- base_palette[overlap]
  missing <- labels[is.na(palette)]
  if (length(missing) > 0) {
    extra_cols <- RColorBrewer::brewer.pal(max(3, length(missing)), "Dark2")
    palette[missing] <- extra_cols[seq_along(missing)]
  }
  palette
}

plot_main_overview <- function(df, palette, subtitle_text) {
  ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(
      data   = df %>% filter(type == "Tumour"),
      aes(colour = disease, shape = type),
      size   = 2.5,
      alpha  = 0.95,
      stroke = 0.5
    ) +
    geom_point(
      data   = df %>% filter(type == "Cell Line"),
      aes(colour = disease, shape = type),
      size   = 2.0,
      alpha  = 0.95,
      stroke = 0.5
    ) +
    scale_colour_manual(values = palette, name = "Cancer type") +
    scale_shape_manual(values = c("Tumour" = 16, "Cell Line" = 4), name = "Sample type") +
    coord_fixed(ratio = 1) +
    theme_umap_main() +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1)),
      shape  = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1))
    ) +
    labs(
      title    = "Multicohort Tumour–Cell Line Alignment",
      subtitle = subtitle_text,
      caption  = "Colour: cancer type; shape: circle = tumour, cross = cell line",
      x = "UMAP 1", y = "UMAP 2"
    )
}

save_plot <- function(plot, path_base, width = 18, height = 10) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, units = "in")
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
}

# Helper: Select candidate genes using cheap methods (for Spearman/MX efficiency)
# Returns union of top N genes from Variance, MAD, MeanAbsDev, and Entropy
select_candidate_genes <- function(X, top_n_per_method = 3000) {
  # X is samples x genes
  n_genes <- ncol(X)
  if (n_genes <= top_n_per_method) {
    # If we have fewer genes than candidates needed, return all
    return(colnames(X))
  }
  
  # METHOD 1: Variance
  var_scores <- matrixStats::colVars(X, na.rm = TRUE)
  names(var_scores) <- colnames(X)
  top_var <- names(sort(var_scores, decreasing = TRUE))[1:min(top_n_per_method, length(var_scores))]
  
  # METHOD 2: MAD (using matrixStats for speed)
  medians <- matrixStats::colMedians(X, na.rm = TRUE)
  mad_scores <- matrixStats::colMads(X, center = medians, na.rm = TRUE)
  names(mad_scores) <- colnames(X)
  top_mad <- names(sort(mad_scores, decreasing = TRUE))[1:min(top_n_per_method, length(mad_scores))]
  
  # METHOD 3: Mean Absolute Deviation
  mean_vals <- colMeans(X, na.rm = TRUE)
  absdev <- colMeans(abs(sweep(X, 2, mean_vals, "-")), na.rm = TRUE)
  names(absdev) <- colnames(X)
  top_mean_absdev <- names(sort(absdev, decreasing = TRUE))[1:min(top_n_per_method, length(absdev))]
  
  # METHOD 4: Entropy
  entropy_scores <- apply(X, 2, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < 2 || var(v) == 0) return(0)
    h <- hist(v, plot=FALSE, breaks=30)$counts
    p <- h / sum(h)
    p <- p[p > 0]
    -sum(p * log(p))
  })
  names(entropy_scores) <- colnames(X)
  top_entropy <- names(sort(entropy_scores, decreasing = TRUE))[1:min(top_n_per_method, length(entropy_scores))]
  
  # Union of all candidates
  candidates <- unique(c(top_var, top_mad, top_mean_absdev, top_entropy))
  return(candidates)
}

# Feature selection (works on X: samples x genes)
compute_feature_scores <- function(X, method) {
  if (method == "Variance") {
    matrixStats::colVars(X, na.rm = TRUE)
  } else if (method == "MAD") {
    med <- matrixStats::colMedians(X, na.rm = TRUE)
    matrixStats::colMads(X, center = med, na.rm = TRUE)
  } else if (method == "MeanAbsDev") {
    # Mean Absolute Deviation from mean (vectorized)
    mu <- colMeans(X, na.rm = TRUE)
    colMeans(abs(sweep(X, 2, mu, "-")), na.rm = TRUE)
  } else if (method == "Entropy") {
    apply(X, 2, function(v) {
      h <- hist(v, plot=FALSE, breaks=30)$counts
      p <- h / sum(h)
      p <- p[p > 0]
      -sum(p * log(p))
    })
  } else if (method == "Spearman") {
    # EFFICIENT: Compute Spearman connectivity only on candidate genes
    # This reduces computation from O(n²) on all genes to O(m²) on candidates (m << n)
    
    # Step 1: Select candidate genes using cheap methods
    candidates <- select_candidate_genes(X, top_n_per_method = 3000)
    cat(sprintf("[INFO] Spearman: computing on %d candidate genes (from %d total)\n", 
                length(candidates), ncol(X)))
    
    # Step 2: Subset to candidates (X is samples x genes, so subset columns)
    X_cand <- X[, candidates, drop = FALSE]
    
    # Step 3: Transpose to genes x samples for ranking
    X_cand_t <- t(X_cand)
    
    # Step 4: Pre-rank all genes once (Spearman = Pearson on ranks)
    # matrixStats::rowRanks is much faster than apply(..., rank)
    Xr <- matrixStats::rowRanks(X_cand_t, ties.method = "average", preserveShape = TRUE)
    rownames(Xr) <- rownames(X_cand_t)
    colnames(Xr) <- colnames(X_cand_t)
    
    # Step 5: Block-wise correlation computation (memory-efficient)
    g <- nrow(Xr)  # Number of candidate genes
    block <- 500   # Block size
    avg_corr <- numeric(g)
    names(avg_corr) <- rownames(Xr)
    
    cat(sprintf("[INFO] Spearman: computing block-wise correlations (block size: %d)\n", block))
    for (i in seq(1, g, by = block)) {
      i2 <- min(g, i + block - 1)
      Xi <- Xr[i:i2, , drop = FALSE]
      
      # Compute Pearson correlation on ranks (equivalent to Spearman)
      C <- cor(t(Xi), t(Xr), method = "pearson", use = "pairwise.complete.obs")
      
      # Remove self-correlation (diagonal = 1.0) for proper network connectivity
      idx <- match(rownames(Xi), rownames(Xr))
      ok <- !is.na(idx)
      if (any(ok)) {
        C[cbind(which(ok), idx[ok])] <- NA_real_
      }
      
      # Store mean absolute correlation for each gene in this block
      avg_corr[i:i2] <- rowMeans(abs(C), na.rm = TRUE)
      
      # Log every 10 blocks (corrected condition)
      block_id <- ((i - 1) %/% block) + 1
      if (block_id %% 10 == 1 || i2 == g) {
        cat(sprintf("[INFO] Spearman blocks: %d-%d / %d\n", i, i2, g))
      }
      
      rm(C)
    }
    
    rm(Xr)
    gc()  # Single gc() call after loop
    
    # Step 6: Return scores for all genes (candidates have real scores, others are NA)
    all_scores <- rep(NA_real_, ncol(X))
    names(all_scores) <- colnames(X)
    all_scores[names(avg_corr)] <- avg_corr
    
    return(all_scores)
  } else if (method == "MX") {
    # EFFICIENT: MX = Spearman connectivity × variance
    # Reuse the same candidate-based Spearman computation
    
    # Step 1: Select candidate genes
    candidates <- select_candidate_genes(X, top_n_per_method = 3000)
    cat(sprintf("[INFO] MX: computing on %d candidate genes (from %d total)\n", 
                length(candidates), ncol(X)))
    
    # Step 2: Subset to candidates
    X_cand <- X[, candidates, drop = FALSE]
    X_cand_t <- t(X_cand)
    
    # Step 3: Compute Spearman connectivity (same as Spearman method above)
    Xr <- matrixStats::rowRanks(X_cand_t, ties.method = "average", preserveShape = TRUE)
    rownames(Xr) <- rownames(X_cand_t)
    colnames(Xr) <- colnames(X_cand_t)
    
    g <- nrow(Xr)
    block <- 500
    avg_corr <- numeric(g)
    names(avg_corr) <- rownames(Xr)
    
    cat(sprintf("[INFO] MX: computing Spearman connectivity (block size: %d)\n", block))
    for (i in seq(1, g, by = block)) {
      i2 <- min(g, i + block - 1)
      Xi <- Xr[i:i2, , drop = FALSE]
      C <- cor(t(Xi), t(Xr), method = "pearson", use = "pairwise.complete.obs")
      
      # Remove self-correlation (diagonal = 1.0) for proper network connectivity
      idx <- match(rownames(Xi), rownames(Xr))
      ok <- !is.na(idx)
      if (any(ok)) {
        C[cbind(which(ok), idx[ok])] <- NA_real_
      }
      
      avg_corr[i:i2] <- rowMeans(abs(C), na.rm = TRUE)
      
      # Log every 10 blocks (corrected condition)
      block_id <- ((i - 1) %/% block) + 1
      if (block_id %% 10 == 1 || i2 == g) {
        cat(sprintf("[INFO] MX Spearman blocks: %d-%d / %d\n", i, i2, g))
      }
      
      rm(C)
    }
    
    rm(Xr)
    gc()  # Single gc() call after loop
    
    # Step 4: Compute variance for candidate genes
    gene_var <- matrixStats::colVars(X_cand, na.rm = TRUE)
    names(gene_var) <- colnames(X_cand)
    
    # Align gene_var to match avg_corr ordering (critical for element-wise multiplication)
    gene_var <- gene_var[names(avg_corr)]
    stopifnot(!any(is.na(gene_var)))
    
    # Step 5: Scale variance to [0, 1]
    max_var <- max(gene_var, na.rm = TRUE)
    gene_var_scaled <- gene_var / max_var
    
    # Step 6: Compute MX score = connectivity × scaled variance
    mx_score <- avg_corr * gene_var_scaled
    names(mx_score) <- names(avg_corr)
    
    # Step 7: Return scores for all genes
    all_scores <- rep(NA_real_, ncol(X))
    names(all_scores) <- colnames(X)
    all_scores[names(mx_score)] <- mx_score
    
    return(all_scores)
  } else if (method == "kTotal") {
    # EFFICIENT: kTotal connectivity computed only on candidate genes
    # kTotal = sum of absolute correlations (total connectivity)
    # This is similar to Spearman but uses Pearson correlation directly
    
    # Step 1: Select candidate genes using cheap methods
    candidates <- select_candidate_genes(X, top_n_per_method = 3000)
    cat(sprintf("[INFO] kTotal: computing on %d candidate genes (from %d total)\n", 
                length(candidates), ncol(X)))
    
    # Step 2: Subset to candidates (X is samples x genes, so subset columns)
    X_cand <- X[, candidates, drop = FALSE]
    
    # Step 3: Standardize columns (Z-score) for correlation stability
    X_cand_scaled <- scale(X_cand, center = TRUE, scale = TRUE)
    
    # Remove constant/NA columns after scaling (genes with zero variance)
    keep <- colSums(is.na(X_cand_scaled)) < nrow(X_cand_scaled)
    if (sum(!keep) > 0) {
      cat(sprintf("[INFO] kTotal: removing %d constant/NA genes after scaling\n", sum(!keep)))
      X_cand_scaled <- X_cand_scaled[, keep, drop = FALSE]
    }
    
    # Step 4: Block-wise correlation computation (memory-efficient)
    p <- ncol(X_cand_scaled)  # Number of candidate genes
    block <- 500
    kt <- numeric(p)
    names(kt) <- colnames(X_cand_scaled)
    
    cat(sprintf("[INFO] kTotal: computing block-wise correlations (block size: %d)\n", block))
    for (start in seq(1, p, by = block)) {
      end <- min(start + block - 1, p)
      Xi <- X_cand_scaled[, start:end, drop = FALSE]
      
      # Compute Pearson correlation between block and all candidates
      C <- cor(Xi, X_cand_scaled, method = "pearson", use = "pairwise.complete.obs")
      
      # Remove self-correlation (diagonal = 1.0) for proper network connectivity
      # For each row in C (corresponding to a gene in Xi), find its column index in X_cand_scaled
      gene_names_block <- colnames(Xi)
      gene_names_all <- colnames(X_cand_scaled)
      idx <- match(gene_names_block, gene_names_all)
      ok <- !is.na(idx)
      if (any(ok)) {
        C[cbind(which(ok), idx[ok])] <- NA_real_
      }
      
      # kTotal = sum of absolute correlations per gene (total connectivity)
      kt[start:end] <- rowSums(abs(C), na.rm = TRUE)
      # Ensure non-finite values are set to NA
      kt[start:end][!is.finite(kt[start:end])] <- NA_real_
      
      # Log every 10 blocks (corrected condition)
      block_id <- ((start - 1) %/% block) + 1
      if (block_id %% 10 == 1 || end == p) {
        cat(sprintf("[INFO] kTotal blocks: %d-%d / %d\n", start, end, p))
      }
      
      rm(C)
    }
    gc()  # Single gc() call after loop instead of every iteration
    
    # Step 5: Return scores for all genes (candidates have real scores, others are NA)
    all_scores <- rep(NA_real_, ncol(X))
    names(all_scores) <- colnames(X)
    all_scores[names(kt)] <- kt
    
    return(all_scores)
  } else {
    stop("Unknown feature method: ", method)
  }
}

# Fast-ish kNN in 2D (UMAP space)
# Tries RANN if available, else falls back to dist() (OK for smaller N)
knn_indices_2d <- function(M, k) {
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data=M, query=M, k=k+1)$nn.idx
    nn[, -1, drop=FALSE]
  } else {
    # fallback: O(N^2) memory/time
    D <- as.matrix(dist(M))
    apply(D, 1, function(d) order(d)[2:(k+1)])
  }
}

evaluate_embedding <- function(coords_df, k=30) {
  # coords_df must have: sample_id, UMAP1, UMAP2, cancer_type, sample_type
  M <- as.matrix(coords_df[, c("UMAP1","UMAP2")])
  rownames(M) <- coords_df$sample_id

  nn <- knn_indices_2d(M, k=k)  # returns matrix: k x N (or N x k depending on fallback)
  if (nrow(nn) == k) nn <- t(nn) # ensure N x k

  is_tumour <- coords_df$sample_type == "Tumour"
  cancer    <- coords_df$cancer_type

  frac_tumour_neighbors <- sapply(seq_len(nrow(coords_df)), function(i) {
    mean(is_tumour[nn[i, ]])
  })

  same_cancer_purity <- sapply(seq_len(nrow(coords_df)), function(i) {
    mean(cancer[nn[i, ]] == cancer[i])
  })

  cell_idx <- which(!is_tumour)

  tibble(
    mean_frac_tumour_neighbors_celllines = mean(frac_tumour_neighbors[cell_idx], na.rm=TRUE),
    median_frac_tumour_neighbors_celllines = median(frac_tumour_neighbors[cell_idx], na.rm=TRUE),
    mean_same_cancer_purity_all = mean(same_cancer_purity, na.rm=TRUE),
    median_same_cancer_purity_all = median(same_cancer_purity, na.rm=TRUE)
  )
}

# Stratified downsample if requested
stratified_downsample <- function(X, meta, max_n) {
  if (max_n <= 0 || nrow(meta) <= max_n) return(list(X=X, meta=meta))

  meta <- meta %>%
    mutate(stratum = paste(cancer_type, sample_type, sep="__"))

  n_total <- nrow(meta)
  frac <- max_n / n_total

  set.seed(opt$seed)
  keep_ids <- meta %>%
    group_by(stratum) %>%
    slice_sample(prop = frac) %>%
    ungroup() %>%
    pull(sample_id)

  keep_ids <- unique(keep_ids)
  X2 <- X[keep_ids, , drop=FALSE]
  meta2 <- meta %>% filter(sample_id %in% keep_ids) %>% slice(match(keep_ids, sample_id))

  cat(sprintf("[INFO] Downsampled from %d -> %d samples (stratified)\n", n_total, nrow(meta2)))
  list(X=X2, meta=meta2)
}

# ─────────────────────────────────────────────────────────────
# Load disease-level joint matrices
# ─────────────────────────────────────────────────────────────
inputs <- tibble(
  profile = profiles,
  path = vst_paths,
  cancer_type = profile_to_cancer(profiles)
) %>%
  mutate(path = vapply(path, abs_from_root, character(1)))

for (i in seq_len(nrow(inputs))) {
  if (!file.exists(inputs$path[i])) {
    stop("Missing input RDS: ", inputs$path[i], " (profile: ", inputs$profile[i], ")")
  }
}

source(file.path(opt$pipe_root, "R", "multicohort_gene_space.R"))

cat("[INFO] Loading joint matrices:\n")
print(inputs)

X_list <- map2(inputs$path, inputs$cancer_type, function(p, ct) {
  cat("  -", ct, ":", p, "\n")
  x <- readRDS(p)
  X <- as_samples_x_genes(x)
  # ensure rownames exist
  if (is.null(rownames(X))) stop("No rownames (sample IDs) in: ", p)
  # make sample IDs unique across cancers (they should already be, but be safe)
  if (anyDuplicated(rownames(X)) > 0) {
    warning("Duplicated sample IDs inside ", ct, " matrix; making unique via make.unique()")
    rownames(X) <- make.unique(rownames(X))
  }
  X
})
names(X_list) <- inputs$cancer_type

# Harmonise genes across all included cancers.
# G_common = intersection of the cohort gene sets, matched on exact gene
# identifiers. harmonise_gene_space() rejects duplicate or blank identifiers,
# subsets every cohort to exactly G_common, and asserts identical gene order
# across cohorts so that concatenation can never align genes by position.
.gene_space <- harmonise_gene_space(X_list)
X_list      <- .gene_space$X_list
genes_shared <- .gene_space$genes_common

cat("[INFO] Shared genes across all cancers:", length(genes_shared), "\n")
cat("[INFO] Common-gene-space audit:\n")
print(.gene_space$audit, row.names = FALSE)

dir.create(file.path(outdir, "inputs"), recursive = TRUE, showWarnings = FALSE)
.gene_space_audit_path <- file.path(outdir, "inputs", "common_gene_space_audit.tsv")
write.table(.gene_space$audit, .gene_space_audit_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[INFO] Wrote common-gene-space audit:", .gene_space_audit_path, "\n")

# Combine samples
X_joint <- do.call(rbind, X_list)
cat("[INFO] Joint matrix:", nrow(X_joint), "samples x", ncol(X_joint), "genes\n")

# Build metadata
meta <- imap_dfr(X_list, function(X, ct) {
  tibble(
    sample_id   = rownames(X),
    cancer_type = ct,
    sample_type = infer_sample_type(rownames(X))
  )
}) %>%
  mutate(
    cohort = ifelse(sample_type == "Tumour",
                    ifelse(grepl("^TARGET", sample_id), "TARGET",
                           ifelse(grepl("^GSE", sample_id), "GSE", "TCGA")),
                    "DSMZ")
  )

# Ensure ordering matches X_joint
meta <- meta %>% slice(match(rownames(X_joint), sample_id))
stopifnot(all(meta$sample_id == rownames(X_joint)))

cat("[INFO] Sample summary:\n")
print(table(meta$cancer_type, meta$sample_type))

# Optional downsampling
if (opt$max_samples > 0) {
  ds <- stratified_downsample(X_joint, meta, opt$max_samples)
  X_joint <- ds$X
  meta    <- ds$meta
}

# Save joint inputs in the pipeline-wide genes x samples convention.
X_joint_out <- t(X_joint)
stopifnot(all(colnames(X_joint_out) == meta$sample_id))
cat("[INFO] Saving joint matrix:", nrow(X_joint_out), "genes x",
    ncol(X_joint_out), "samples\n")
saveRDS(X_joint_out, file.path(outdir, "inputs", "joint_expr_matrix.rds"))
write_tsv(meta, file.path(outdir, "inputs", "joint_metadata.tsv"))
cat("[SAVED] inputs/joint_expr_matrix.rds\n")
cat("[SAVED] inputs/joint_metadata.tsv\n")

# ─────────────────────────────────────────────────────────────
# Extra outputs for downstream ensemble reuse
# ─────────────────────────────────────────────────────────────
featset_dir  <- file.path(outdir, sprintf("feature_sets_top%d", opt$n_genes))
submat_root  <- file.path(outdir, "featuresets")   # mirrors per-disease layout
dir.create(featset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(submat_root, recursive = TRUE, showWarnings = FALSE)

all_summaries <- list()

for (fm in feature_methods) {
  cat("\n[INFO] Feature method:", fm, "\n")
  
  # Check if outputs already exist (resume/skip logic)
  gene_list_path <- file.path(featset_dir, sprintf("genes_top%d_%s.txt", opt$n_genes, fm))
  fm_outdir <- file.path(submat_root, fm)
  sub_rds_path <- file.path(fm_outdir, "expr_submatrix.rds")
  
  # Check if all UMAP outputs exist for this method
  all_coords_exist <- TRUE
  for (dm in dist_metrics) {
    coords_path <- file.path(outdir, "benchmarks", paste0("coords_", fm, "_", dm, ".tsv"))
    if (!file.exists(coords_path)) {
      all_coords_exist <- FALSE
      break
    }
  }
  
  if (file.exists(gene_list_path) && file.exists(sub_rds_path) && all_coords_exist) {
    # A cached feature list is only reusable when every gene it names is still
    # present in the current common gene space. Cohort matrices regenerated with
    # different gene retention shrink G_common, which leaves earlier lists
    # referencing absent identifiers; subsetting on those would fail.
    genes_keep <- readLines(gene_list_path)
    .fl <- validate_feature_list(genes_keep, colnames(X_joint), basename(gene_list_path))
    .cache_usable <- .fl$ok
    if (!.cache_usable) {
      cat("  [STALE] ", .fl$label, ": ", .fl$n_missing,
          " gene(s) absent from the current common gene space; recomputing.\n", sep = "")
      cat("          examples: ",
          paste(utils::head(.fl$missing, 5), collapse = ", "), "\n", sep = "")
    }
  } else {
    .cache_usable <- FALSE
  }

  if (.cache_usable) {
    cat("  [SKIP] Outputs exist for ", fm, ", skipping feature selection and UMAP.\n", sep="")
    X_sub <- X_joint[, genes_keep, drop=FALSE]
  } else {
    # Compute feature scores
    scores <- compute_feature_scores(X_joint, fm)
    scores <- scores[!is.na(scores)]
    ord <- order(scores, decreasing=TRUE)
    genes_keep <- names(scores)[ord][seq_len(min(opt$n_genes, length(ord)))]
    if (length(genes_keep) < 50) stop("Too few genes selected for ", fm, ": ", length(genes_keep))

    # Save gene list (downstream reproducibility)
    writeLines(genes_keep, con = gene_list_path)

    # Save submatrix RDS in a per-method folder (mirrors per-disease pipeline)
    dir.create(fm_outdir, recursive = TRUE, showWarnings = FALSE)

    X_sub <- X_joint[, genes_keep, drop=FALSE]
    saveRDS(X_sub, sub_rds_path)
    
    cat("  [SAVED] ", file.path(sprintf("feature_sets_top%d", opt$n_genes), basename(gene_list_path)), "\n", sep="")
    cat("  [SAVED] ", file.path("featuresets", fm, "expr_submatrix.rds"), "\n", sep="")
  }

  for (dm in dist_metrics) {
    tag <- paste0(fm, "_", dm)
    coords_path <- file.path(outdir, "benchmarks", paste0("coords_", tag, ".tsv"))
    
    # Skip UMAP if coords already exist
    if (file.exists(coords_path)) {
      cat("  [SKIP] ", tag, " (coords exist)\n", sep="")
      coords <- read_tsv(coords_path, show_col_types = FALSE)
      eval_tbl <- evaluate_embedding(coords, k=opt$eval_k)
    } else {
      cat("  [RUN] ", tag, " ...\n", sep="")

      set.seed(opt$seed)
      emb <- uwot::umap(
        X           = X_sub,
        n_neighbors = opt$n_neighbors,
        min_dist    = opt$min_dist,
        metric      = dm,
        n_threads   = n_threads,
        scale       = FALSE,
        verbose     = TRUE
      )

      coords <- as.data.frame(emb)
      colnames(coords) <- c("UMAP1","UMAP2")
      coords$sample_id <- rownames(X_sub)

      coords <- coords %>%
        left_join(meta, by="sample_id")

      # Evaluate in UMAP space (kNN)
      eval_tbl <- evaluate_embedding(coords, k=opt$eval_k)

      # Save coords
      write_tsv(coords, coords_path)
    }  # closes the else{}

    # Plot (main overview style)
    plot_df <- prepare_plot_data(coords)
    subtitle <- sprintf(
      "Joint tumours vs joint cell lines; genes=%d (%s top-%d); metric=%s; nn=%d; min_dist=%.2f",
      length(genes_keep), fm, length(genes_keep), dm, opt$n_neighbors, opt$min_dist
    )

    palette <- build_disease_palette(levels(droplevels(plot_df$disease)))
    p <- plot_main_overview(plot_df, palette, subtitle)

    plot_base <- file.path(outdir, "benchmarks", paste0("01_main_overview_", tag))
    save_plot(p, plot_base, width=18, height=10)

    # Summary row
    srow <- tibble(
      feature_method = fm,
      dist_metric    = dm,
      n_genes        = length(genes_keep),
      n_samples      = nrow(X_joint),
      n_tumours      = sum(meta$sample_type == "Tumour"),
      n_cell_lines   = sum(meta$sample_type == "Cell Line"),
      eval_k         = opt$eval_k,
      mean_frac_tumour_neighbors_celllines = eval_tbl$mean_frac_tumour_neighbors_celllines,
      median_frac_tumour_neighbors_celllines = eval_tbl$median_frac_tumour_neighbors_celllines,
      mean_same_cancer_purity_all = eval_tbl$mean_same_cancer_purity_all,
      median_same_cancer_purity_all = eval_tbl$median_same_cancer_purity_all,
      coords_tsv     = basename(coords_path),
      plot_pdf       = basename(paste0(plot_base, ".pdf")),
      plot_png       = basename(paste0(plot_base, ".png"))
    )

    all_summaries[[tag]] <- srow
    cat("  [SAVED] ", basename(coords_path), " + plot\n", sep="")
  }
}

summary_tbl <- bind_rows(all_summaries) %>%
  # Ranking idea:
  # - Prefer cell lines surrounded by tumours (higher mean_frac_tumour_neighbors_celllines)
  # - But also preserve cancer structure (higher mean_same_cancer_purity_all)
  arrange(
    desc(mean_frac_tumour_neighbors_celllines),
    desc(mean_same_cancer_purity_all)
  )

summary_path <- file.path(outdir, "benchmarks", "summary.tsv")
write_tsv(summary_tbl, summary_path)

# Write a simple manifest of feature lists produced
feat_manifest <- tibble(
  feature_method = feature_methods,
  gene_list = file.path(sprintf("feature_sets_top%d", opt$n_genes), sprintf("genes_top%d_%s.txt", opt$n_genes, feature_methods)),
  submatrix_rds = file.path("featuresets", feature_methods, "expr_submatrix.rds")
)
write_tsv(feat_manifest, file.path(outdir, "feature_sets_manifest.tsv"))

cat("\n[SAVED] benchmarks/summary.tsv\n")
cat("[SAVED] feature_sets_manifest.tsv\n")
cat("\n=== TOP 10 (by ranking) ===\n")
print(head(summary_tbl, 10))

cat("\n[SUCCESS] Joint benchmark complete.\n")
