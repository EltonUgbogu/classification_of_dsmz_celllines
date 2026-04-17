#!/usr/bin/env Rscript
# tcga_consensus_clustering_methods.R
# Build TCGA-only consensus clustering from multiple methods
# K is chosen data-driven by silhouette on the co-association matrix

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(readr)
  library(cluster)
})

source("scripts/lib_config.R")

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------------
# Load config.yaml
# ------------------------------------------------------------------
# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
unsup_root <- cfg$paths$unsup_root

# ------------------------------------------------------------------
# Function: Build consensus from multiple clustering methods
# ------------------------------------------------------------------
run_tcga_inductive_consensus <- function(
    base_dir = file.path(unsup_root, "inductive_tcga"),
    k_range = 2:8,
    out_prefix = "BRCA_TCGA_PAM50_CONSENSUS"
) {
    cat("\n=== TCGA Inductive Consensus Clustering ===\n")
    cat("Base directory:", base_dir, "\n")
    
    # Define input files (canonical names)
    cluster_files <- list(
        km_tcga = file.path(
            base_dir,
            "consensusclustering",
            "km_tcga",
            "BRCA_TCGA_PAM50_km_clusters_bestk.csv"
        ),
        km_pca_cc = file.path(
            base_dir,
            "consensusclustering",
            "kmeans_pca_tcga",
            "consensus",
            "consensus_clusters_K3.csv"
        ),
        hc_pca_opt = file.path(
            base_dir,
            "pca_hc_tcga",
            "PCA_TCGA_BRCA_PAM50_samples_clusters_optimal_k=2.csv"
        ),
        km_pca_tcga = file.path(
            base_dir,
            "pca_kmeans_tcga",
            "patient_clustering_kmeans_tcga",
            "final_clusters.csv"
        )
    )
    
    method_names <- names(cluster_files)
    
    # Check all files exist
    missing <- cluster_files[!sapply(cluster_files, file.exists)]
    if (length(missing) > 0) {
        stop("Missing cluster files:\n  ", paste(missing, collapse = "\n  "))
    }
    
    cat("[INFO] Loading cluster assignments from", length(cluster_files), "methods\n")
    
    # Load all cluster assignments
    clusters_list <- list()
    for (method in method_names) {
        df <- readr::read_csv(cluster_files[[method]], show_col_types = FALSE)
        
        # Handle different column names
        if ("sample_id" %in% colnames(df)) {
            sample_col <- "sample_id"
        } else if ("sample" %in% colnames(df)) {
            sample_col <- "sample"
        } else {
            sample_col <- colnames(df)[1]
        }
        
        cluster_col <- if ("cluster" %in% colnames(df)) "cluster" else colnames(df)[2]
        
        clusters_list[[method]] <- setNames(df[[cluster_col]], df[[sample_col]])
        cat("  ", method, ":", length(clusters_list[[method]]), "samples\n")
    }
    
    # Find common samples
    common_samples <- Reduce(intersect, lapply(clusters_list, names))
    cat("[INFO] Number of common tumours across all methods:", length(common_samples), "\n")
    
    if (length(common_samples) < 10) {
        stop("Too few common samples (", length(common_samples), ") for consensus clustering")
    }
    
    # Build co-association matrix
    n_samples <- length(common_samples)
    coassoc <- matrix(0, nrow = n_samples, ncol = n_samples)
    rownames(coassoc) <- colnames(coassoc) <- common_samples
    
    for (method in method_names) {
        cl <- clusters_list[[method]][common_samples]
        for (i in seq_len(n_samples)) {
            for (j in seq_len(n_samples)) {
                if (cl[i] == cl[j]) {
                    coassoc[i, j] <- coassoc[i, j] + 1
                }
            }
        }
    }
    
    coassoc <- coassoc / length(method_names)
    cat("[INFO] Co-association matrix built:", paste(dim(coassoc), collapse = " x "), "\n")
    
    # Convert to distance matrix
    dist_mat <- as.dist(1 - coassoc)
    
    # Hierarchical clustering on co-association
    hc <- hclust(dist_mat, method = "ward.D2")
    
    # Find optimal K via silhouette
    sil_scores <- numeric(length(k_range))
    names(sil_scores) <- k_range
    
    for (k in k_range) {
        cut_k <- cutree(hc, k = k)
        sil <- silhouette(cut_k, dist_mat)
        sil_scores[as.character(k)] <- mean(sil[, "sil_width"])
    }
    
    best_k <- k_range[which.max(sil_scores)]
    best_sil <- max(sil_scores)
    
    cat("[INFO] Chosen K =", best_k, "(max silhouette =", round(best_sil, 4), ")\n")
    cat("[INFO] Silhouette scores (K", paste(k_range, collapse = ", "), "):",
        paste(round(sil_scores, 4), collapse = ", "), "\n")
    
    # Get final clusters
    final_clusters <- cutree(hc, k = best_k)
    final_clusters <- paste0("C", final_clusters)
    names(final_clusters) <- common_samples
    
    # Create output directory
    out_dir <- base_dir
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write canonical consensus clusters file
    clusters_df <- data.frame(
        sample_id = names(final_clusters),
        cluster = final_clusters,
        stringsAsFactors = FALSE
    )
    
    clusters_file <- file.path(out_dir, paste0(out_prefix, "_clusters_bestk.csv"))
    readr::write_csv(clusters_df, clusters_file)
    cat("[OUTPUT] Canonical consensus clusters saved to:", clusters_file, "\n")
    
    # Write best-K summary metadata
    summary_df <- data.frame(
        method = "consensus_from_methods",
        best_k = best_k,
        silhouette_score = round(best_sil, 4),
        n_samples = length(common_samples),
        n_methods = length(method_names),
        methods = paste(method_names, collapse = ";"),
        stringsAsFactors = FALSE
    )
    
    summary_file <- file.path(out_dir, paste0(out_prefix, "_best_k_summary.tsv"))
    readr::write_tsv(summary_df, summary_file)
    cat("[OUTPUT] Canonical best-K summary saved to:", summary_file, "\n")
    
    # Also write K-specific file for reference (if best_k != 3, still write k3 for compatibility)
    if (best_k == 3) {
        k3_file <- file.path(out_dir, paste0(out_prefix, "_k3_SAMPLES.csv"))
        readr::write_csv(clusters_df, k3_file)
    } else {
        # Write k3 version for backward compatibility
        k3_clusters <- cutree(hc, k = 3)
        k3_clusters <- paste0("C", k3_clusters)
        names(k3_clusters) <- common_samples
        k3_df <- data.frame(
            sample_id = names(k3_clusters),
            cluster = k3_clusters,
            stringsAsFactors = FALSE
        )
        k3_file <- file.path(out_dir, paste0(out_prefix, "_k3_SAMPLES.csv"))
        readr::write_csv(k3_df, k3_file)
        cat("[OUTPUT] K=3 clusters (for compatibility) saved to:", k3_file, "\n")
    }
    
    invisible(list(
        clusters = final_clusters,
        best_k = best_k,
        silhouette = best_sil,
        coassoc = coassoc
    ))
}

# ------------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------------
run_tcga_inductive_consensus(
    base_dir = file.path(unsup_root, "inductive_tcga"),
    k_range = 2:8,
    out_prefix = "BRCA_TCGA_PAM50_CONSENSUS"
)

