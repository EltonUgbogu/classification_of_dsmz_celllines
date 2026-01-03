#!/usr/bin/env Rscript

# inductive_tcga_tumour_clustering.R
# Evaluate the final TCGA consensus clustering vs PAM50,
# generate UMAPs, confusion matrices, and clustering metrics.

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(uwot)
  library(cluster)
  library(pheatmap)
  library(grid)
  library(aricode)
  library(scales)
  library(patchwork)
  library(yaml)
  library(optparse)
})

cat("=== TCGA inductive consensus: clustering evaluation ===\n")

# ----------------------------------------------------------------------
# 0) CLI: get method parameter
# ----------------------------------------------------------------------
option_list <- list(
  make_option("--method", type = "character", default = "pam50_euc",
              help = "Method ID: pam50_corr, pam50_euc, hvg_corr, or hvg_euc [default: %default]"),
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
method_id <- opt$method

# Load method-aware helper
source("R/brca_unsup_methods.R")

# ----------------------------------------------------------------------
# 0a) Load config.yaml
# ----------------------------------------------------------------------
cfg <- yaml::read_yaml(opt$config)
unsup_root <- cfg$paths$unsup_root
data_root <- cfg$paths$data_root

# Get method context
m_ctx <- get_method_context(cfg, method_id)

cat("[INFO] Method:", method_id,
    "| feature:", m_ctx$feature,
    "| distance:", m_ctx$distance, "\n")
cat("[INFO] Cluster prefix:", m_ctx$cluster_prefix, "\n")

# ----------------------------------------------------------------------
# 0a) Paths and parameters (from config.yaml)
# ----------------------------------------------------------------------

# PAM50 subtype annotation
subtype_csv <- cfg$paths$tcga_subtype_csv

# Final consensus clusters (best K chosen data-driven) - same for all methods
# Now produced by tcga_inductive_clustering_methods.R under inductive_clustering_methods/
consensus_csv <- file.path(
  unsup_root,
  "inductive_tcga",
  "inductive_clustering_methods",
  "BRCA_TCGA_PAM50_CONSENSUS_clusters_bestk.csv"
)

# Load expression matrix based on method feature
if (m_ctx$feature == "PAM50") {
  vst_rds <- cfg$paths$tcga_pam50_expr
  use_joint_hvg <- FALSE
} else if (m_ctx$feature == "HVG") {
  # Use the joint DSMZ+TCGA HVG500 matrix that already exists
  # We'll subset to TCGA samples after loading
  vst_rds <- file.path(
    unsup_root,
    "tumour_neighbourhoods_input",
    "BRCA_TCGA-DSMZ_HVG500_samples_x_genes.rds"
  )
  use_joint_hvg <- TRUE
} else {
  stop("Unknown feature type: ", m_ctx$feature)
}

cat("[INFO] Using expression matrix:", vst_rds, "\n")
if (use_joint_hvg) {
  cat("[INFO] Will subset to TCGA samples from joint HVG matrix\n")
}

# Output directory for evaluation
out_dir <- file.path(
  unsup_root,
  "inductive_tcga",
  "inductive_tcga_clusters"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] VST RDS         : ", vst_rds, "\n", sep = "")
cat("[INFO] PAM50 subtype   : ", subtype_csv, "\n", sep = "")
cat("[INFO] Consensus CSV   : ", consensus_csv, "\n", sep = "")
cat("[INFO] Output dir      : ", out_dir, "\n\n", sep = "")

# UMAP / silhouette outputs (for TCGA-only consensus)
out_umap_pdf <- file.path(
  unsup_root,
  "inductive_tcga",
  "TCGA_inductive_consensus_UMAP_clusters.pdf"
)
out_umap_subtype <- file.path(
  unsup_root,
  "inductive_tcga",
  "TCGA_inductive_consensus_UMAP_PAM50.pdf"
)
out_umap_combined <- file.path(
  out_dir,
  "TCGA_inductive_consensus_UMAP_combined.pdf"
)
out_sil_pdf <- file.path(
  unsup_root,
  "inductive_tcga",
  "TCGA_inductive_consensus_silhouette.pdf"
)

# Confusion & mapping outputs
conf_counts_csv <- file.path(
  out_dir,
  "cluster_PAM50_confusion_counts.csv"
)
conf_props_csv <- file.path(
  out_dir,
  "cluster_PAM50_confusion_props_rowwise.csv"
)
conf_heat_counts <- file.path(
  out_dir,
  "Fig_cluster_PAM50_confusion_heatmap_counts.pdf"
)
conf_heat_props <- file.path(
  out_dir,
  "Fig_cluster_PAM50_confusion_heatmap_props.pdf"
)
cluster_mapping_csv <- file.path(
  out_dir,
  "cluster_majority_PAM50_mapping_accuracy.csv"
)

# Canonical metrics file Snakemake depends on (method-aware filename)
metrics_tsv <- file.path(
  out_dir,
  sprintf("%s_CONSENSUS_CLUSTERING_METRICS.tsv", m_ctx$cluster_prefix)
)

# Optional: neighbourhood tumour matrix (if tumour-neighbourhood pipeline ran)
nh_tumours_path <- file.path(
  data_root,
  "bcc_preprocessing",
  "results",
  "TCGA-BRCA_induced_PAM50_expr_matrix_neighbourhood_tumours.rds"
)

# ----------------------------------------------------------------------
# 1. LOAD VST + CONSENSUS CLUSTERS
# ----------------------------------------------------------------------
if (!file.exists(vst_rds)) {
  stop("VST RDS not found: ", vst_rds)
}
if (!file.exists(consensus_csv)) {
  stop("Consensus CSV not found: ", consensus_csv)
}
if (!file.exists(subtype_csv)) {
  stop("Subtype CSV not found: ", subtype_csv)
}

vst <- readRDS(vst_rds)

# For HVG: the joint matrix has samples as ROWS (transposed relative to PAM50)
# Subset to TCGA samples only
if (use_joint_hvg) {
  cat("[INFO] Loaded joint HVG matrix: ", nrow(vst), " samples × ",
      ncol(vst), " genes\n", sep = "")
  
  # Matrix is samples × genes; identify TCGA samples by barcode pattern
  tcga_idx <- grepl("^TCGA-", rownames(vst))
  if (!any(tcga_idx)) {
    stop("No TCGA samples detected in HVG matrix rownames.")
  }
  vst <- vst[tcga_idx, , drop = FALSE]
  cat("[INFO] Subset to TCGA: ", nrow(vst), " samples\n", sep = "")
  
  # Transpose to genes × samples format (like PAM50 matrix)
  vst <- t(vst)
  cat("[INFO] Transposed to: ", nrow(vst), " genes × ",
      ncol(vst), " tumours\n", sep = "")
} else {
  cat("[INFO] Loaded VST matrix: ", nrow(vst), " genes × ",
      ncol(vst), " tumours\n", sep = "")
}

clusters <- read.csv(consensus_csv, header = TRUE, stringsAsFactors = FALSE)
if (ncol(clusters) == 2) {
  colnames(clusters) <- c("sample", "cluster")
} else {
  stop("Consensus CSV must contain exactly two columns: sample, cluster")
}
clusters$sample <- as.character(clusters$sample)

# Align VST and clusters by sample
common <- intersect(colnames(vst), clusters$sample)
if (length(common) == 0) {
  stop("No overlap between VST colnames and consensus cluster samples.")
}

vst <- vst[, common, drop = FALSE]
clusters <- clusters[match(common, clusters$sample), , drop = FALSE]
stopifnot(identical(colnames(vst), clusters$sample))

# Derive K from data (no hard-coded k)
k_value <- length(unique(clusters$cluster))
cat("[INFO] Final consensus K = ", k_value, "\n", sep = "")

# ----------------------------------------------------------------------
# 1b. LOAD PAM50 SUBTYPES AND ALIGN
# ----------------------------------------------------------------------
pam50_df <- read.csv(subtype_csv, stringsAsFactors = FALSE) %>%
  transmute(
    sample = sample_barcode,
    PAM50_Subtype = PAM50_Subtype
  )

pam50_df <- pam50_df[pam50_df$sample %in% colnames(vst), , drop = FALSE]

cat("[INFO] PAM50 annotations available for ",
    nrow(pam50_df), " tumours\n", sep = "")

# ----------------------------------------------------------------------
# 1c. OPTIONAL: PAM50 SUBTYPES AMONG NEIGHBOURHOOD TUMOURS
# ----------------------------------------------------------------------
if (file.exists(nh_tumours_path)) {
  cat("\n[INFO] Counting PAM50 subtypes among neighbourhood tumours...\n")
  nh_mat <- readRDS(nh_tumours_path)
  nh_tumour_ids <- colnames(nh_mat)
  n_tumours <- length(nh_tumour_ids)

  cat("  Number of neighbourhood tumours: ", n_tumours, "\n", sep = "")

  pam50_df_harmonized <- pam50_df %>%
    mutate(sample = gsub("_", "-", sample, fixed = TRUE))

  nh_subtypes <- pam50_df_harmonized %>%
    filter(sample %in% nh_tumour_ids)

  cat("  Tumours with PAM50 annotation: ", nrow(nh_subtypes), "\n", sep = "")

  subtype_counts_df <- nh_subtypes %>%
    count(PAM50_Subtype, name = "n") %>%
    mutate(
      proportion = n / sum(n),
      percentage = round(100 * proportion, 1)
    ) %>%
    arrange(desc(n))

  cat("\nPAM50 subtype distribution among neighbourhood tumours:\n")
  print(subtype_counts_df)

  subtype_counts_csv <- file.path(
    out_dir,
    "PAM50_subtype_counts_neighbourhood_tumours.csv"
  )
  write.csv(subtype_counts_df, subtype_counts_csv, row.names = FALSE)
  cat("[OUTPUT] Subtype counts saved to: ", subtype_counts_csv, "\n", sep = "")
} else {
  cat("[INFO] Neighbourhood tumour matrix not found at:\n  ",
      nh_tumours_path, "\n",
      "→ Skipping subtype summary for neighbourhood tumours.\n\n", sep = "")
}

# ----------------------------------------------------------------------
# 2. UMAP
# ----------------------------------------------------------------------
set.seed(42)
cat("[INFO] Computing UMAP on VST (TCGA only)...\n")
um <- umap(
  t(vst),
  n_neighbors = 20,
  min_dist = 0.3,
  metric = "cosine"
)

um_df <- as.data.frame(um)
colnames(um_df) <- c("UMAP1", "UMAP2")
um_df$sample <- rownames(um)
um_df$cluster <- factor(clusters$cluster)

# Attach PAM50 subtype (may be NA for some samples)
um_df <- um_df %>%
  left_join(pam50_df, by = "sample")

um_df_sub <- um_df %>% filter(!is.na(PAM50_Subtype))

# ----------------------------------------------------------------------
# 2a. COLOURS & THEMES (now in R/brca_unsup_methods.R)
# ----------------------------------------------------------------------
# pam50_colours, cluster_colours, and theme_dimred are loaded from shared helpers
theme_umap <- theme_dimred  # alias for backward compatibility

# ----------------------------------------------------------------------
# 2b. UMAP plotting helpers
# ----------------------------------------------------------------------
plot_umap_pam50 <- function(um_df,
                            out_pdf = NULL,
                            title = "TCGA BRCA inductive consensus clusters",
                            subtitle = "UMAP coloured by PAM50 subtype") {

  required_cols <- c("UMAP1", "UMAP2", "PAM50_Subtype")
  if (!all(required_cols %in% colnames(um_df))) {
    stop("um_df must contain columns: UMAP1, UMAP2, PAM50_Subtype")
  }

  plot_df <- um_df %>%
    filter(!is.na(PAM50_Subtype)) %>%
    mutate(subtype = factor(PAM50_Subtype, levels = names(pam50_colours)))

  if (nrow(plot_df) == 0) {
    stop("No samples with valid PAM50_Subtype found")
  }

  present_subtypes <- levels(droplevels(plot_df$subtype))
  plot_colours <- pam50_colours[present_subtypes]

  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, colour = subtype)) +
    geom_point(size = 1.8, alpha = 0.75, shape = 16) +
    scale_colour_manual(
      values = plot_colours,
      name = "PAM50 subtype",
      na.value = "grey70"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 11) +
    guides(
      colour = guide_legend(
        title = "PAM50 subtype",
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2")

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(out_png, p, width = 7, height = 5.5, units = "in", dpi = 300)
    cat("[SUCCESS] Saved PAM50 UMAP to:\n  ", out_pdf, "\n  ", out_png, "\n")
  }

  p
}

plot_umap_cluster <- function(um_df,
                              out_pdf = NULL,
                              title = "TCGA BRCA inductive consensus clusters",
                              subtitle = NULL,
                              k = NULL) {

  required_cols <- c("UMAP1", "UMAP2", "cluster")
  if (!all(required_cols %in% colnames(um_df))) {
    stop("um_df must contain columns: UMAP1, UMAP2, cluster")
  }

  plot_df <- um_df %>%
    filter(!is.na(cluster)) %>%
    mutate(cluster = factor(cluster))

  if (nrow(plot_df) == 0) {
    stop("No samples with valid cluster assignment found")
  }

  if (is.null(k)) {
    k <- length(unique(plot_df$cluster))
  }
  if (is.null(subtitle)) {
    subtitle <- sprintf("UMAP coloured by inductive cluster (K = %d)", k)
  }

  present_clusters <- as.character(
    sort(unique(as.integer(as.character(plot_df$cluster))))
  )
  plot_colours <- cluster_colours[present_clusters]

  if (any(is.na(plot_colours))) {
    extra_colours <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra_colours
  }

  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 1.8, alpha = 0.75, shape = 16) +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 11) +
    guides(
      colour = guide_legend(
        title = "Cluster",
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2")

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(out_png, p, width = 7, height = 5.5, units = "in", dpi = 300)
    cat("[SUCCESS] Saved cluster UMAP to:\n  ", out_pdf, "\n  ", out_png, "\n")
  }

  p
}

plot_umap_combined <- function(um_df,
                               out_pdf = NULL,
                               k = NULL) {

  required_cols <- c("UMAP1", "UMAP2", "cluster", "PAM50_Subtype")
  if (!all(required_cols %in% colnames(um_df))) {
    stop(
      "um_df must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  plot_df <- um_df %>%
    filter(!is.na(PAM50_Subtype), !is.na(cluster)) %>%
    mutate(
      subtype = factor(PAM50_Subtype, levels = names(pam50_colours)),
      cluster = factor(cluster)
    )

  if (nrow(plot_df) == 0) {
    stop("No samples with both valid PAM50_Subtype and cluster found")
  }

  present_subtypes <- levels(droplevels(plot_df$subtype))
  pam50_plot_colours <- pam50_colours[present_subtypes]

  present_clusters <- as.character(
    sort(unique(as.integer(as.character(plot_df$cluster))))
  )
  cluster_plot_colours <- cluster_colours[present_clusters]

  p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 1.5, alpha = 0.75, shape = 16) +
    scale_colour_manual(values = cluster_plot_colours, name = "Cluster") +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 10) +
    guides(
      colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))
    ) +
    labs(title = "Inductive cluster", x = "UMAP 1", y = "UMAP 2")

  p2 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, colour = subtype)) +
    geom_point(size = 1.5, alpha = 0.75, shape = 16) +
    scale_colour_manual(
      values = pam50_plot_colours,
      name = "PAM50 subtype"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 10) +
    guides(
      colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))
    ) +
    labs(title = "PAM50 subtype", x = "UMAP 1", y = "UMAP 2")

  combined <- p1 + p2 +
    plot_annotation(
      title = "TCGA BRCA inductive consensus clusters",
      subtitle = "Comparison of inductive clustering vs PAM50 molecular subtype",
      theme = theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
          colour = "grey10"
        ),
        plot.subtitle = element_text(size = 11, colour = "grey40")
      )
    )

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, combined, width = 12, height = 5.5, units = "in", dpi = 300)
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(out_png, combined, width = 12, height = 5.5, units = "in", dpi = 300)
    cat(
      "[SUCCESS] Saved combined UMAP to:\n  ",
      out_pdf, "\n  ", out_png, "\n"
    )
  }

  combined
}

# ----------------------------------------------------------------------
# 2c. Generate UMAP plots
# ----------------------------------------------------------------------
plot_umap_cluster(um_df, out_pdf = out_umap_pdf, k = k_value)
plot_umap_pam50(um_df, out_pdf = out_umap_subtype)
plot_umap_combined(um_df, out_pdf = out_umap_combined, k = k_value)

# ----------------------------------------------------------------------
# 3. SILHOUETTE
# ----------------------------------------------------------------------
cat("[INFO] Computing silhouette on consensus clusters...\n")
d <- dist(t(vst), method = "euclidean")
sil <- silhouette(as.integer(as.factor(clusters$cluster)), d)

pdf(out_sil_pdf, width = 6, height = 5)
plot(
  sil,
  col = as.integer(as.factor(clusters$cluster)),
  main = sprintf(
    "Silhouette Plot (K = %d, mean = %.3f)",
    k_value, mean(sil[, "sil_width"])
  )
)
dev.off()

cat("[SUCCESS] Silhouette plot saved to: ", out_sil_pdf, "\n", sep = "")

# ----------------------------------------------------------------------
# 4. CLUSTER × PAM50 CONFUSION TABLES + HEATMAPS
# ----------------------------------------------------------------------
cat("[INFO] Building cluster × PAM50 confusion tables...\n")

um_df_sub <- um_df %>% filter(!is.na(PAM50_Subtype))

conf_long <- um_df_sub %>%
  count(cluster, PAM50_Subtype, name = "n") %>%
  arrange(cluster, PAM50_Subtype)

write.csv(conf_long, conf_counts_csv, row.names = FALSE)
cat("[OUTPUT] Confusion counts saved to: ", conf_counts_csv, "\n", sep = "")

conf_counts_mat <- conf_long %>%
  pivot_wider(
    id_cols = cluster,
    names_from = PAM50_Subtype,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

conf_props_mat <- conf_counts_mat / rowSums(conf_counts_mat)

conf_props_df <- conf_props_mat %>%
  as.data.frame() %>%
  rownames_to_column("cluster")

write.csv(conf_props_df, conf_props_csv, row.names = FALSE)
cat("[OUTPUT] Row-wise confusion proportions saved to: ",
    conf_props_csv, "\n", sep = "")

pdf(conf_heat_counts, width = 7, height = 5)
pheatmap::pheatmap(
  conf_counts_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.0f",
  main = "Cluster × PAM50 subtype (counts)"
)
dev.off()
cat("[OUTPUT] Confusion heatmap (counts) saved to: ",
    conf_heat_counts, "\n", sep = "")

pdf(conf_heat_props, width = 7, height = 5)
pheatmap::pheatmap(
  conf_props_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Cluster × PAM50 subtype (row-wise proportions)"
)
dev.off()
cat("[OUTPUT] Confusion heatmap (proportions) saved to: ",
    conf_heat_props, "\n", sep = "")

# ----------------------------------------------------------------------
# 5. MAJORITY-VOTE MAPPING + ACCURACY
# ----------------------------------------------------------------------
cat("[INFO] Computing majority-vote mapping: clusters → PAM50 subtype...\n")

majority_map <- conf_long %>%
  group_by(cluster) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    cluster,
    majority_PAM50 = PAM50_Subtype,
    majority_count = n
  )

um_df_map <- um_df_sub %>%
  left_join(majority_map, by = "cluster") %>%
  mutate(
    predicted_PAM50 = majority_PAM50,
    correct = (predicted_PAM50 == PAM50_Subtype)
  )

overall_accuracy <- mean(um_df_map$correct)

cat(sprintf(
  "[RESULT] Overall majority-vote accuracy: %.3f\n",
  overall_accuracy
))

per_cluster_acc <- um_df_map %>%
  group_by(cluster) %>%
  summarise(
    n_tumours = n(),
    majority_PAM50 = first(majority_PAM50),
    majority_count = first(majority_count),
    accuracy_cluster = mean(correct),
    .groups = "drop"
  )

mapping_out <- per_cluster_acc %>%
  mutate(overall_accuracy = overall_accuracy)

write.csv(mapping_out, cluster_mapping_csv, row.names = FALSE)
cat("[OUTPUT] Majority-vote mapping + accuracies saved to: ",
    cluster_mapping_csv, "\n", sep = "")

# ----------------------------------------------------------------------
# 6. COMPREHENSIVE CLUSTERING EVALUATION METRICS
# ----------------------------------------------------------------------

calculate_ARI <- function(true_labels, pred_labels) {
  aricode::ARI(true_labels, pred_labels)
}

calculate_NMI <- function(true_labels, pred_labels) {
  aricode::NMI(true_labels, pred_labels)
}

calculate_AMI <- function(true_labels, pred_labels) {
  aricode::AMI(true_labels, pred_labels)
}

calculate_v_measure <- function(true_labels, pred_labels, beta = 1.0) {
  cont <- table(true_labels, pred_labels)
  n <- sum(cont)

  p_c <- rowSums(cont) / n
  p_k <- colSums(cont) / n

  H_C <- -sum(p_c * log(p_c + 1e-10))
  H_K <- -sum(p_k * log(p_k + 1e-10))

  H_C_K <- 0
  for (k in seq_len(ncol(cont))) {
    n_k <- sum(cont[, k])
    if (n_k > 0) {
      p_ck <- cont[, k] / n_k
      p_ck <- p_ck[p_ck > 0]
      H_C_K <- H_C_K + (n_k / n) * (-sum(p_ck * log(p_ck)))
    }
  }

  H_K_C <- 0
  for (c in seq_len(nrow(cont))) {
    n_c <- sum(cont[c, ])
    if (n_c > 0) {
      p_kc <- cont[c, ] / n_c
      p_kc <- p_kc[p_kc > 0]
      H_K_C <- H_K_C + (n_c / n) * (-sum(p_kc * log(p_kc)))
    }
  }

  homogeneity <- ifelse(H_C == 0, 1, 1 - H_C_K / H_C)
  completeness <- ifelse(H_K == 0, 1, 1 - H_K_C / H_K)

  if (homogeneity + completeness == 0) {
    v_measure <- 0
  } else {
    v_measure <- ((1 + beta) * homogeneity * completeness) /
      (beta * homogeneity + completeness)
  }

  list(
    homogeneity = homogeneity,
    completeness = completeness,
    v_measure = v_measure
  )
}

calculate_FMI <- function(true_labels, pred_labels) {
  cont <- table(true_labels, pred_labels)

  TP <- sum(choose(cont, 2))
  same_true <- sum(choose(table(true_labels), 2))
  same_pred <- sum(choose(table(pred_labels), 2))

  precision <- ifelse(same_pred == 0, 0, TP / same_pred)
  recall <- ifelse(same_true == 0, 0, TP / same_true)

  sqrt(precision * recall)
}

calculate_purity <- function(true_labels, pred_labels) {
  cont <- table(pred_labels, true_labels)
  sum(apply(cont, 1, max)) / length(true_labels)
}

calculate_f_measure <- function(true_labels, pred_labels) {
  cont <- table(true_labels, pred_labels)
  n <- sum(cont)

  f_scores <- numeric(nrow(cont))

  for (i in seq_len(nrow(cont))) {
    n_i <- sum(cont[i, ])
    best_f <- 0

    for (j in seq_len(ncol(cont))) {
      n_j <- sum(cont[, j])
      n_ij <- cont[i, j]

      if (n_ij > 0) {
        precision <- n_ij / n_j
        recall <- n_ij / n_i
        f <- 2 * precision * recall / (precision + recall)
        best_f <- max(best_f, f)
      }
    }

    f_scores[i] <- best_f * (n_i / n)
  }

  sum(f_scores)
}

evaluate_clustering <- function(true_labels, pred_labels, expr_matrix = NULL) {

  valid_idx <- !is.na(true_labels) & !is.na(pred_labels)
  true_labels <- true_labels[valid_idx]
  pred_labels <- pred_labels[valid_idx]

  true_labels_chr <- as.character(true_labels)
  pred_labels_int <- as.integer(as.factor(pred_labels))

  v_res <- calculate_v_measure(true_labels_chr, pred_labels_int)

  results <- list(
    ARI = calculate_ARI(true_labels_chr, pred_labels_int),
    AMI = calculate_AMI(true_labels_chr, pred_labels_int),
    NMI = calculate_NMI(true_labels_chr, pred_labels_int),
    FMI = calculate_FMI(true_labels_chr, pred_labels_int),
    Purity = calculate_purity(true_labels_chr, pred_labels_int),
    F_measure = calculate_f_measure(true_labels_chr, pred_labels_int),
    Homogeneity = v_res$homogeneity,
    Completeness = v_res$completeness,
    V_measure = v_res$v_measure,
    n_samples = length(true_labels_chr),
    n_true_classes = length(unique(true_labels_chr)),
    n_pred_clusters = length(unique(pred_labels_int))
  )

  if (!is.null(expr_matrix)) {
    if (ncol(expr_matrix) >= length(pred_labels_int)) {
      expr_sub <- expr_matrix[, valid_idx, drop = FALSE]
      d <- dist(t(expr_sub), method = "euclidean")
      sil <- silhouette(pred_labels_int, d)
      results$Silhouette_mean <- mean(sil[, "sil_width"])
    }
  }

  results
}

run_comprehensive_evaluation <- function(um_df_sub, vst = NULL) {

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("COMPREHENSIVE CLUSTERING EVALUATION\n")
  cat(strrep("=", 60), "\n\n")

  true_labels <- um_df_sub$PAM50_Subtype
  pred_labels <- um_df_sub$cluster

  metrics <- evaluate_clustering(true_labels, pred_labels, expr_matrix = vst)

  cat("EXTERNAL VALIDATION METRICS (vs PAM50 ground truth):\n")
  cat(strrep("-", 50), "\n")

  cat(sprintf("  Adjusted Rand Index (ARI):       %.4f\n", metrics$ARI))
  cat(sprintf("  Adjusted Mutual Info (AMI):      %.4f\n", metrics$AMI))
  cat(sprintf("  Normalized Mutual Info (NMI):    %.4f\n", metrics$NMI))
  cat(sprintf("  Fowlkes-Mallows Index (FMI):     %.4f\n", metrics$FMI))
  cat(sprintf("  V-measure:                       %.4f\n", metrics$V_measure))
  cat(sprintf("    - Homogeneity:                 %.4f\n", metrics$Homogeneity))
  cat(sprintf(
    "    - Completeness:                %.4f\n",
    metrics$Completeness
  ))
  cat(sprintf("  F-measure:                       %.4f\n", metrics$F_measure))
  cat(sprintf("  Purity (majority-vote):          %.4f\n", metrics$Purity))

  if (!is.null(metrics$Silhouette_mean)) {
    cat("\nINTERNAL VALIDATION:\n")
    cat(sprintf(
      "  Mean Silhouette Width:           %.4f\n",
      metrics$Silhouette_mean
    ))
  }

  cat("\nSAMPLE INFO:\n")
  cat(sprintf("  Total samples:                   %d\n", metrics$n_samples))
  cat(sprintf(
    "  True classes (PAM50):            %d\n",
    metrics$n_true_classes
  ))
  cat(sprintf(
    "  Predicted clusters:              %d\n",
    metrics$n_pred_clusters
  ))
  cat("\n", strrep("=", 60), "\n\n", sep = "")

  metrics
}

compute_cluster_score <- function(metrics) {
  (metrics$ARI + metrics$NMI + metrics$F_measure) / 3
}

cat("[INFO] Running comprehensive clustering evaluation...\n")
vst_sub <- vst[, um_df_sub$sample, drop = FALSE]

metrics <- run_comprehensive_evaluation(um_df_sub, vst_sub)

# Save metrics as TSV (canonical Snakemake target)
metrics_df <- as.data.frame(metrics)
readr::write_tsv(metrics_df, metrics_tsv)
cat("[OUTPUT] Comprehensive evaluation metrics (TSV) saved to:\n  ",
    metrics_tsv, "\n", sep = "")

cluster_score <- compute_cluster_score(metrics)
cat(sprintf("[RESULT] Overall Cluster Score: %.4f\n", cluster_score))

cat("\n[SUCCESS] Inductive TCGA tumour clustering, subtype alignment and ",
    "evaluation finished.\n", sep = "")
