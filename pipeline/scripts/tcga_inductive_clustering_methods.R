#!/usr/bin/env Rscript

# tcga_inductive_clustering_methods.R
#
# Build a TCGA-only inductive tumour clustering using the
# best clustering per clustering algorithm, then derive a
# consensus clustering via a co-association matrix.
#
# NOW CONFIG-DRIVEN:
#   - base_dir = cfg$paths$unsup_root from config.yaml
#   - CLI: --config, --k_min, --k_max
#   - Outputs under: <unsup_root>/inductive_tcga/inductive_clustering_methods

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(cluster)
  library(optparse)
  library(yaml)
})

# -------------------------------------------------------------------
# 1. Helper: robust loader for cluster label CSVs
# -------------------------------------------------------------------
load_cluster_labels <- function(path, out_dir = NULL, method_name = NULL) {
  message("[INFO] Loading clustering labels from: ", path)
  if (!file.exists(path)) {
    stop("[ERROR] File not found: ", path)
  }

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L) {
    stop("[ERROR] File has zero rows: ", path)
  }

  # Try to detect sample IDs and cluster column
  col_lower <- tolower(colnames(df))

  # 1) Sample column
  sample_candidates <- c(
    "sample", "samples", "id", "patient", "tumour", "tumor"
  )
  sample_col_idx <- which(col_lower %in% sample_candidates)

  if (length(sample_col_idx) > 0) {
    samples <- df[[sample_col_idx[1]]]
  } else if (!is.null(rownames(df)) && any(grepl("^TCGA", rownames(df)))) {
    samples <- rownames(df)
  } else {
    message(
      "[WARN] Could not detect explicit sample column in ", path,
      " - using first column as sample ID."
    )
    samples <- df[[1]]
  }

  # 2) Cluster column
  cluster_candidates <- c(
    "cluster", "clusters", "class",
    "consensusclass", "consensus_class",
    "clusterid", "cluster_id"
  )
  cluster_col_idx <- which(col_lower %in% cluster_candidates)

  if (length(cluster_col_idx) == 0) {
    if (ncol(df) == 2L) {
      cluster_col_idx <- setdiff(seq_len(2L), sample_col_idx[1])
    } else {
      message(
        "[WARN] Could not detect explicit cluster column in ", path,
        " - using last column as cluster label."
      )
      cluster_col_idx <- ncol(df)
    }
  }

  clusters <- df[[cluster_col_idx[1]]]

  # Clean up: drop NAs, coerce to factor then integer codes 1..K
  keep <- !is.na(samples) & !is.na(clusters)
  samples <- samples[keep]
  clusters <- clusters[keep]

  clusters_factor <- as.integer(as.factor(clusters))
  names(clusters_factor) <- samples

  n_samples <- length(clusters_factor)
  n_clusters <- length(unique(clusters_factor))

  message(
    "[INFO] Loaded ", n_samples, " labelled tumours from ", path,
    " (K = ", n_clusters, ")"
  )

  # Save parsed output if out_dir provided
  if (!is.null(out_dir) && !is.null(method_name)) {
    parsed_df <- data.frame(
      sample = names(clusters_factor),
      cluster = clusters_factor,
      stringsAsFactors = FALSE
    )
    parsed_out <- file.path(
      out_dir,
      paste0("parsed_", method_name, "_labels.csv")
    )
    write.csv(parsed_df, parsed_out, row.names = FALSE)
    message("[OUTPUT] Parsed labels saved to: ", parsed_out)
  }

  clusters_factor
}

# -------------------------------------------------------------------
# 2. Build co-association matrix from a list of clusterings
# -------------------------------------------------------------------
build_coassoc_matrix <- function(label_list, out_dir = NULL) {
  if (length(label_list) < 2L) {
    stop("[ERROR] Need at least 2 clustering runs to build a consensus.")
  }

  # Common tumour set across all runs
  common_ids <- Reduce(intersect, lapply(label_list, names))
  common_ids <- sort(common_ids)

  if (length(common_ids) < 5L) {
    stop("[ERROR] Fewer than 5 common tumours across clusterings.")
  }

  message(
    "[INFO] Number of common tumours across all methods: ",
    length(common_ids)
  )

  # Subset each label vector to common IDs in the same order
  label_list <- lapply(label_list, function(x) x[common_ids])

  n <- length(common_ids)
  m <- length(label_list)
  co <- matrix(
    0, nrow = n, ncol = n,
    dimnames = list(common_ids, common_ids)
  )

  # For each clustering run, add 1 to pairs co-clustered
  for (r in seq_len(m)) {
    labs <- label_list[[r]]
    for (cl in unique(labs)) {
      idx <- which(labs == cl)
      co[idx, idx] <- co[idx, idx] + 1
    }
  }

  # Normalise to [0,1] by number of runs
  co <- co / m

  message(
    "[INFO] Co-association matrix built: ", n, " x ", n,
    " from ", m, " clustering runs."
  )

  # Save outputs
  if (!is.null(out_dir)) {
    # Save common tumour IDs
    common_ids_out <- file.path(out_dir, "common_tumour_ids.txt")
    writeLines(common_ids, common_ids_out)
    message("[OUTPUT] Common tumour IDs saved to: ", common_ids_out)

    # Save aligned label matrix (tumours x methods)
    aligned_mat <- do.call(cbind, label_list)
    aligned_df <- as.data.frame(aligned_mat)
    aligned_df <- cbind(sample = rownames(aligned_df), aligned_df)
    aligned_out <- file.path(out_dir, "aligned_labels_all_methods.csv")
    write.csv(aligned_df, aligned_out, row.names = FALSE)
    message("[OUTPUT] Aligned labels matrix saved to: ", aligned_out)

    # Save co-association matrix as RDS
    co_rds <- file.path(out_dir, "coassociation_matrix.rds")
    saveRDS(co, co_rds)
    message("[OUTPUT] Co-association matrix (RDS) saved to: ", co_rds)

    # Save co-association matrix as CSV (for inspection)
    co_df <- as.data.frame(co)
    co_df <- cbind(sample = rownames(co_df), co_df)
    co_csv <- file.path(out_dir, "coassociation_matrix.csv")
    write.csv(co_df, co_csv, row.names = FALSE)
    message("[OUTPUT] Co-association matrix (CSV) saved to: ", co_csv)

    # Save summary statistics
    co_summary <- data.frame(
      metric = c(
        "n_tumours",
        "n_methods",
        "min_coassoc",
        "max_coassoc",
        "mean_coassoc",
        "median_coassoc"
      ),
      value = c(
        n,
        m,
        round(min(co), 4),
        round(max(co), 4),
        round(mean(co), 4),
        round(median(co), 4)
      )
    )
    summary_out <- file.path(out_dir, "coassociation_summary.csv")
    write.csv(co_summary, summary_out, row.names = FALSE)
    message("[OUTPUT] Co-association summary saved to: ", summary_out)
  }

  co
}

# -------------------------------------------------------------------
# 3. Choose best K by silhouette on 1 - coassoc
# -------------------------------------------------------------------
choose_best_k_by_silhouette <- function(co, k_range = 2:8, out_dir = NULL) {
  d <- as.dist(1 - co)
  sil_scores <- numeric(length(k_range))
  sil_details <- list()

  for (i in seq_along(k_range)) {
    k <- k_range[i]
    hc <- hclust(d, method = "average")
    cl <- cutree(hc, k = k)
    sil <- cluster::silhouette(cl, d)
    sil_scores[i] <- mean(sil[, "sil_width"])

    # Store per-sample silhouette for this K
    sil_details[[as.character(k)]] <- data.frame(
      sample = names(cl),
      cluster = cl,
      sil_width = sil[, "sil_width"],
      stringsAsFactors = FALSE
    )
  }

  best_idx <- which.max(sil_scores)
  best_k <- k_range[best_idx]

  # Save outputs
  if (!is.null(out_dir)) {
    # K selection summary
    k_summary <- data.frame(
      K = k_range,
      mean_silhouette = round(sil_scores, 4),
      is_best = k_range == best_k
    )
    k_summary_out <- file.path(out_dir, "k_selection_silhouette.csv")
    write.csv(k_summary, k_summary_out, row.names = FALSE)
    message("[OUTPUT] K selection summary saved to: ", k_summary_out)

    # Per-sample silhouette for each K
    for (k in k_range) {
      sil_k_out <- file.path(
        out_dir,
        sprintf("silhouette_per_sample_K%d.csv", k)
      )
      write.csv(sil_details[[as.character(k)]], sil_k_out, row.names = FALSE)
    }
    message(
      "[OUTPUT] Per-sample silhouette files saved for K = ",
      paste(k_range, collapse = ", ")
    )
  }

  list(
    best_k = best_k,
    k_range = k_range,
    sil_scores = sil_scores,
    sil_details = sil_details
  )
}

# -------------------------------------------------------------------
# 4. Final clustering and comprehensive outputs
# -------------------------------------------------------------------
perform_final_clustering <- function(co, best_k, out_dir = NULL) {
  d <- as.dist(1 - co)
  hc <- hclust(d, method = "average")
  clusters <- cutree(hc, k = best_k)

  # Silhouette at best K
  sil <- cluster::silhouette(clusters, d)
  mean_sil <- mean(sil[, "sil_width"])

  message("[INFO] Final consensus clusters: K = ", best_k)
  print(table(clusters))
  message(
    sprintf("[INFO] Mean silhouette width at K = %d: %.4f", best_k, mean_sil)
  )

  # Save outputs
  if (!is.null(out_dir)) {
    # hclust object
    hc_out <- file.path(out_dir, "consensus_hclust.rds")
    saveRDS(hc, hc_out)
    message("[OUTPUT] hclust object saved to: ", hc_out)

    # Dendrogram as PDF
    dend_out <- file.path(out_dir, "consensus_dendrogram.pdf")
    pdf(dend_out, width = 12, height = 6)
    plot(
      hc,
      main = sprintf(
        "Consensus Clustering Dendrogram (K = %d, mean sil = %.3f)",
        best_k, mean_sil
      ),
      xlab = "Tumour samples",
      sub = "",
      cex = 0.3
    )
    rect.hclust(hc, k = best_k, border = 2:10)
    dev.off()
    message("[OUTPUT] Dendrogram saved to: ", dend_out)

    # Final cluster assignments
    consensus_df <- data.frame(
      sample = names(clusters),
      cluster = clusters,
      stringsAsFactors = FALSE
    )
    clusters_out <- file.path(
      out_dir,
      sprintf("final_consensus_clusters_K%d.csv", best_k)
    )
    write.csv(consensus_df, clusters_out, row.names = FALSE)
    message("[OUTPUT] Final cluster assignments saved to: ", clusters_out)

    # Cluster size summary
    cluster_sizes <- as.data.frame(table(clusters))
    colnames(cluster_sizes) <- c("cluster", "n_tumours")
    cluster_sizes$proportion <- round(
      cluster_sizes$n_tumours / sum(cluster_sizes$n_tumours), 4
    )
    sizes_out <- file.path(
      out_dir,
      sprintf("cluster_sizes_K%d.csv", best_k)
    )
    write.csv(cluster_sizes, sizes_out, row.names = FALSE)
    message("[OUTPUT] Cluster sizes saved to: ", sizes_out)

    # Per-sample silhouette at best K
    sil_df <- data.frame(
      sample = names(clusters),
      cluster = clusters,
      sil_width = sil[, "sil_width"],
      stringsAsFactors = FALSE
    )
    sil_out <- file.path(
      out_dir,
      sprintf("silhouette_per_sample_final_K%d.csv", best_k)
    )
    write.csv(sil_df, sil_out, row.names = FALSE)
    message("[OUTPUT] Final silhouette per sample saved to: ", sil_out)

    # Silhouette summary
    sil_summary <- sil_df %>%
      group_by(cluster) %>%
      summarise(
        n_samples = n(),
        mean_sil = round(mean(sil_width), 4),
        min_sil = round(min(sil_width), 4),
        max_sil = round(max(sil_width), 4),
        n_negative = sum(sil_width < 0),
        .groups = "drop"
      )
    sil_summary <- rbind(
      sil_summary,
      data.frame(
        cluster = "OVERALL",
        n_samples = nrow(sil_df),
        mean_sil = round(mean_sil, 4),
        min_sil = round(min(sil_df$sil_width), 4),
        max_sil = round(max(sil_df$sil_width), 4),
        n_negative = sum(sil_df$sil_width < 0)
      )
    )
    sil_summary_out <- file.path(
      out_dir,
      sprintf("silhouette_summary_K%d.csv", best_k)
    )
    write.csv(sil_summary, sil_summary_out, row.names = FALSE)
    message("[OUTPUT] Silhouette summary saved to:", sil_summary_out)

    # Silhouette plot
    sil_plot_out <- file.path(
      out_dir,
      sprintf("silhouette_plot_K%d.pdf", best_k)
    )
    pdf(sil_plot_out, width = 8, height = 6)
    plot(
      sil,
      main = sprintf(
        "Silhouette Plot (K = %d, mean = %.3f)",
        best_k, mean_sil
      ),
      col = 1:best_k,
      border = NA
    )
    dev.off()
    message("[OUTPUT] Silhouette plot saved to:", sil_plot_out)
  }

  list(
    hclust = hc,
    clusters = clusters,
    silhouette = sil,
    mean_sil = mean_sil
  )
}

# -------------------------------------------------------------------
# 5. Main driver (now config + CLI)
# -------------------------------------------------------------------
run_tcga_inductive_consensus <- function(
  base_dir,
  k_range = 2:8,
  out_prefix = "BRCA_TCGA_PAM50_CONSENSUS"
) {
  message("[INFO] Base directory: ", base_dir)
  base_dir <- normalizePath(base_dir, mustWork = TRUE)

  # Create output directory (this is where EVERYTHING goes)
  out_dir <- file.path(base_dir, "inductive_tcga", "inductive_clustering_methods")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("[INFO] Output directory: ", out_dir)

  # ---- 5.1. Auto-detect PCA+HC files (no hard-coded k) ----
  pca_hc_dir <- file.path(base_dir, "pca_hc_tcga")

  # Optimal K file
  opt_candidates <- list.files(
    pca_hc_dir,
    pattern = "PCA_TCGA_BRCA_PAM50_samples_clusters_optimal_k=.*[.]csv",
    full.names = TRUE
  )
  if (length(opt_candidates) == 0L) {
    stop("[ERROR] No optimal PCA+HC clusters file found in: ", pca_hc_dir)
  }
  if (length(opt_candidates) > 1L) {
    opt_candidates <- sort(opt_candidates)
    message("[WARN] Multiple optimal PCA+HC files found. Using: ", opt_candidates[1])
  }
  pca_hc_opt_file <- opt_candidates[1]

  # Dynamic K file
  dyn_candidates <- list.files(
    pca_hc_dir,
    pattern = "PCA_TCGA_BRCA_PAM50_samples_clusters_dynamic_k=.*[.]csv",
    full.names = TRUE
  )
  if (length(dyn_candidates) == 0L) {
    stop("[ERROR] No dynamic PCA+HC clusters file found in: ", pca_hc_dir)
  }
  if (length(dyn_candidates) > 1L) {
    dyn_candidates <- sort(dyn_candidates)
    message("[WARN] Multiple dynamic PCA+HC files found. Using: ", dyn_candidates[1])
  }
  pca_hc_dyn_file <- dyn_candidates[1]

  # ---- 5.1b. Auto-detect raw HC TCGA files (optimal / dynamic) ----
  hc_dir <- file.path(base_dir, "hc_tcga")

  hc_opt_candidates <- list.files(
    hc_dir,
    pattern = "TCGA-BRCA_PAM50_patients_clusters_optimal_k=.*[.]csv",
    full.names = TRUE
  )
  if (length(hc_opt_candidates) == 0L) {
    stop("[ERROR] No optimal raw HC clusters file found in: ", hc_dir)
  }
  if (length(hc_opt_candidates) > 1L) {
    hc_opt_candidates <- sort(hc_opt_candidates)
    message("[WARN] Multiple optimal raw HC files found. Using: ", hc_opt_candidates[1])
  }
  hc_tcga_opt_file <- hc_opt_candidates[1]

  hc_dyn_candidates <- list.files(
    hc_dir,
    pattern = "TCGA-BRCA_PAM50_patients_clusters_dynamic_k=.*[.]csv",
    full.names = TRUE
  )
  if (length(hc_dyn_candidates) == 0L) {
    stop("[ERROR] No dynamic raw HC clusters file found in: ", hc_dir)
  }
  if (length(hc_dyn_candidates) > 1L) {
    hc_dyn_candidates <- sort(hc_dyn_candidates)
    message("[WARN] Multiple dynamic raw HC files found. Using: ", hc_dyn_candidates[1])
  }
  hc_tcga_dyn_file <- hc_dyn_candidates[1]

  # ---- 5.2. Define the "best per method" cluster files ----
  # All methods correspond to TCGA-only clustering outputs.
  cluster_files <- c(
    file.path(
      base_dir, "kmeans_tcga", "patient_clustering", "final_clusters.csv"
    ),                                                       # 1) kmeans_tcga
    file.path(
      base_dir, "pca_km_tcga",
      "patient_clustering_kmeans_tcga", "final_clusters.csv"
    ),                                                       # 2) pca_km_tcga
    pca_hc_opt_file,                                         # 3) pca_hc_tcga_opt
    pca_hc_dyn_file,                                         # 4) pca_hc_tcga_dyn
    hc_tcga_opt_file,                                        # 5) hc_tcga_opt
    hc_tcga_dyn_file,                                        # 6) hc_tcga_dyn
    file.path(
      base_dir, "consensusclustering", "km_tcga", "tcga_km_clusters_K3.csv"
    ),                                                       # 7) ccp_km_tcga
    file.path(
      base_dir, "consensusclustering", "hc_tcga", "tcga_hc_clusters_K3.csv"
    ),                                                       # 8) ccp_hc_tcga
    file.path(
      base_dir, "consensusclustering", "kmeans_pca_tcga",
      "consensus", "consensus_clusters_K3.csv"
    ),                                                       # 9) ccp_kmeans_pca_tcga
    file.path(
      base_dir, "consensusclustering", "hc_pca_tcga",
      "consensus", "consensus_clusters_K3.csv"
    )                                                        # 10) ccp_hc_pca_tcga
  )

  method_names <- c(
    "kmeans_tcga",
    "pca_km_tcga",
    "pca_hc_tcga_opt",
    "pca_hc_tcga_dyn",
    "hc_tcga_opt",
    "hc_tcga_dyn",
    "ccp_km_tcga",
    "ccp_hc_tcga",
    "ccp_kmeans_pca_tcga",
    "ccp_hc_pca_tcga"
  )

  # Sanity check: all files exist
  missing <- !file.exists(cluster_files)
  if (any(missing)) {
    stop(
      "[ERROR] The following cluster files are missing:\n",
      paste0("  - ", cluster_files[missing], collapse = "\n")
    )
  }

  # Save method-to-file mapping
  method_mapping <- data.frame(
    method = method_names,
    source_file = cluster_files,
    stringsAsFactors = FALSE
  )
  mapping_out <- file.path(out_dir, "input_method_mapping.csv")
  write.csv(method_mapping, mapping_out, row.names = FALSE)
  message("[OUTPUT] Method mapping saved to:", mapping_out)

  # ---- 5.3. Load all label vectors ----
  label_list <- list()

  for (i in seq_along(cluster_files)) {
    path <- cluster_files[i]
    labs <- load_cluster_labels(
      path,
      out_dir = out_dir,
      method_name = method_names[i]
    )
    label_list[[method_names[i]]] <- labs
  }

  # Save loading summary
  loading_summary <- data.frame(
    method = method_names,
    n_samples = sapply(label_list, length),
    n_clusters = sapply(label_list, function(x) length(unique(x))),
    stringsAsFactors = FALSE
  )
  loading_out <- file.path(out_dir, "loading_summary.csv")
  write.csv(loading_summary, loading_out, row.names = FALSE)
  message("[OUTPUT] Loading summary saved to:", loading_out)

  # ---- 5.4. Build co-association matrix ----
  co <- build_coassoc_matrix(label_list, out_dir = out_dir)

  # ---- 5.5. Choose best K by silhouette ----
  best_result <- choose_best_k_by_silhouette(
    co,
    k_range = k_range,
    out_dir = out_dir
  )

  best_k <- best_result$best_k
  message(
    sprintf(
      "[INFO] Chosen K = %d (max silhouette = %.4f)",
      best_k, max(best_result$sil_scores)
    )
  )

  # ---- 5.6. Final hierarchical clustering at best K ----
  final_result <- perform_final_clustering(co, best_k, out_dir = out_dir)

  # ---- 5.7. Save final summary report ----
  files_vec <- list.files(
    out_dir,
    pattern = "\\.(csv|rds|pdf|txt)$",
    full.names = FALSE
  )
  files_lines <- paste0("  - ", files_vec)
  sep_line <- paste(rep("=", 60), collapse = "")

  summary_report <- c(
    sep_line,
    "TCGA INDUCTIVE CONSENSUS CLUSTERING - SUMMARY REPORT",
    sep_line,
    "",
    sprintf("Output directory: %s", out_dir),
    sprintf("Number of methods combined: %d", length(method_names)),
    sprintf("Methods: %s", paste(method_names, collapse = ", ")),
    "",
    sprintf("Common tumours across all methods: %d", nrow(co)),
    sprintf("K range tested: %d - %d", min(k_range), max(k_range)),
    sprintf("Best K (by silhouette): %d", best_k),
    sprintf("Mean silhouette at best K: %.4f", final_result$mean_sil),
    "",
    "Cluster sizes:",
    capture.output(print(table(final_result$clusters))),
    "",
    "Output files generated:",
    files_lines,
    "",
    sep_line
  )

  report_out <- file.path(out_dir, "summary_report.txt")
  writeLines(summary_report, report_out)
  message("[OUTPUT] Summary report saved to:", report_out)

  # ---- 5.8. Canonical outputs (USED BY OTHER SCRIPTS) ----
  consensus_df <- data.frame(
    sample = names(final_result$clusters),
    cluster = final_result$clusters,
    stringsAsFactors = FALSE
  )

  canonical_clusters <- file.path(
    out_dir,
    sprintf("%s_clusters_bestk.csv", out_prefix)
  )
  write.csv(consensus_df, canonical_clusters, row.names = FALSE)
  message("[OUTPUT] Canonical consensus clusters saved to:", canonical_clusters)

  bestk_meta <- data.frame(
    best_k = best_k,
    mean_silhouette = round(final_result$mean_sil, 4)
  )
  canonical_bestk <- file.path(
    out_dir,
    sprintf("%s_best_k_summary.tsv", out_prefix)
  )
  readr::write_tsv(bestk_meta, canonical_bestk)
  message("[OUTPUT] Canonical best-K summary saved to:", canonical_bestk)

  message("\n[SUCCESS] Inductive TCGA consensus clustering complete.")
  message("[SUCCESS] All outputs saved to:", out_dir)

  invisible(list(
    coassoc = co,
    hclust = final_result$hclust,
    clusters = final_result$clusters,
    best_k = best_k,
    mean_sil = final_result$mean_sil,
    out_dir = out_dir
  ))
}

# -------------------------------------------------------------------
# 6. CLI entry point (config + K range)
# -------------------------------------------------------------------
if (sys.nframe() == 0) {
  option_list <- list(
    make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
                help = "Path to config.yaml [default: %default]"),
    make_option("--profile", type = "character", default = NULL,
                help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
    make_option("--k_min", type = "integer", default = 2,
                help = "Minimum K for silhouette search [default: %default]"),
    make_option("--k_max", type = "integer", default = 8,
                help = "Maximum K for silhouette search [default: %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  source("scripts/lib_config.R")
  # Load config using profiled config loader (merges defaults + profile)
  profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
  cfg <- read_profiled_config(opt$config, profile)
  base_dir <- cfg$paths$unsup_root

  k_range <- seq.int(opt$k_min, opt$k_max)

  run_tcga_inductive_consensus(
    base_dir  = base_dir,
    k_range   = k_range,
    out_prefix = "BRCA_TCGA_PAM50_CONSENSUS"
  )
}

