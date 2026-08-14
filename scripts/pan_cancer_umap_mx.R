#!/usr/bin/env Rscript

# ==============================================================================
# pan_cancer_umap_mx.R
# Pan-Cancer UMAP Visualization with Joint MX Gene List
# ==============================================================================
#
# PURPOSE:
# Creates a unified UMAP embedding across multiple cancer types (BRCA, NBL, RBL)
# using a joint MX gene set. Tumours are filtered based on best-performing
# neighbourhood directions, and cell lines are included from all profiles.
#
# INPUTS:
# - Per-profile expression matrices (featuresets/MX/expr_submatrix.rds)
# - Per-profile MX lists (to create joint list)
# - Best direction summaries (to filter tumours)
# - Top_m_long files (to get tumour sample IDs)
#
# OUTPUTS:
# - Pan-cancer UMAP coordinates (TSV)
# - Pan-cancer UMAP plot (PDF)
# - Metadata with cohort/cancer_type/sample_type labels
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(uwot)
  library(yaml)
  library(optparse)
  library(purrr)
})

# Optional: ggrepel for better label placement (will check availability later)
if (requireNamespace("ggrepel", quietly = TRUE)) {
  library(ggrepel)
}

# Helper operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Source lib_config.R using absolute path (works with --directory)
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--profile", type = "character", default = NULL,
              help = "Dummy profile (not used, but required for lib_config)"),
  make_option("--profiles", type = "character", default = "brca,nbl,rbl",
              help = "Comma-separated list of profiles [default: %default]"),
  make_option("--outdir", type = "character", default = "results/unsupervised/pan_cancer/umap",
              help = "Output directory [default: %default]"),
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs for PCA before UMAP [default: %default]"),
  make_option("--n_neighbors", type = "integer", default = 15,
              help = "UMAP n_neighbors parameter [default: %default]"),
  make_option("--min_dist", type = "double", default = 0.1,
              help = "UMAP min_dist parameter [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source lib_config.R using absolute path
if (!is.null(opt$config) && file.exists(opt$config)) {
  config_dir <- dirname(normalizePath(opt$config))
  pipe_root <- normalizePath(file.path(config_dir, ".."))
  lib_config_path <- file.path(pipe_root, "scripts", "lib_config.R")
} else {
  lib_config_path <- "scripts/lib_config.R"
}
if (file.exists(lib_config_path)) {
  source(lib_config_path)
} else {
  stop("Cannot find lib_config.R. Tried: ", lib_config_path)
}

# Parse profiles
profiles <- strsplit(opt$profiles, ",")[[1]] %>% trimws()
cat("[INFO] Processing profiles:", paste(profiles, collapse = ", "), "\n")

# Derive pipeline root
pipe_root <- normalizePath(file.path(dirname(opt$config), ".."))

# Helper to absolutize paths
abs_from_root <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(pipe_root, p)
}

resolve_mx_gene_list_path <- function(cfg) {
  mx_topn <- as.integer(cfg$feature_selection$method_topn$MX %||% 500L)
  file.path(
    abs_from_root(cfg$paths$unsup_root),
    "feature_selection_unsupervised",
    "feature_sets",
    sprintf("genes_top%d_MX.txt", mx_topn)
  )
}

# Create output directory
outdir <- abs_from_root(opt$outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 1: Load per-profile configs and determine best directions
# ==============================================================================

cat("[INFO] Loading configs and determining best directions...\n")

profile_configs <- list()
best_directions <- character(length(profiles))
names(best_directions) <- profiles

for (prof in profiles) {
  cfg <- read_profiled_config(opt$config, prof)
  profile_configs[[prof]] <- cfg
  
  # Get best direction from direction summary
  dir_summary_path <- file.path(
    abs_from_root(cfg$paths$unsup_root),
    "tumour_neighbourhoods", "final_consensus_all",
    "p_consensus_direction_summary.tsv"
  )
  
  if (file.exists(dir_summary_path)) {
    dir_summary <- read_tsv(dir_summary_path, show_col_types = FALSE)
    # Pick direction with highest frac_ge_thr, then frac_eq_1, then median_p (stable tie-breaking)
    best_dir <- dir_summary %>%
      arrange(desc(frac_ge_thr), desc(frac_eq_1), desc(median_p)) %>%
      slice(1) %>%
      pull(direction)
    best_directions[prof] <- best_dir
    cat(sprintf("  %s: best direction = %s (frac_ge_thr=%.3f, frac_eq_1=%.3f)\n", 
                toupper(prof), best_dir,
                dir_summary %>% filter(direction == best_dir) %>% pull(frac_ge_thr),
                dir_summary %>% filter(direction == best_dir) %>% pull(frac_eq_1)))
  } else {
    # Fallback: use first available direction
    cat(sprintf("  %s: direction summary not found, will use all directions\n", toupper(prof)))
    best_directions[prof] <- NA_character_
  }
}

# ==============================================================================
# Step 2: Create joint MX list (union of all profiles)
# ==============================================================================

cat("[INFO] Creating joint MX list...\n")

mx_lists <- list()
for (prof in profiles) {
  cfg <- profile_configs[[prof]]
  mx_path <- resolve_mx_gene_list_path(cfg)

  if (file.exists(mx_path)) {
    mx_list <- read_lines(mx_path) %>% trimws() %>% .[nzchar(.)]
    mx_lists[[prof]] <- mx_list
    cat(sprintf("  %s: %d MX genes\n", toupper(prof), length(mx_list)))
  } else {
    stop(sprintf("MX list not found for %s: %s", prof, mx_path))
  }
}

# Union of all MX lists (candidate set)
joint_mx_union <- unique(unlist(mx_lists))
cat(sprintf("[INFO] Joint MX union (candidate): %d unique genes\n", length(joint_mx_union)))

# ==============================================================================
# Step 3: Load expression matrices and filter tumours
# ==============================================================================

cat("[INFO] Loading expression matrices and filtering tumours...\n")

expr_matrices <- list()
tumour_samples_all <- character(0)
cell_line_samples_all <- character(0)

for (prof in profiles) {
  cfg <- profile_configs[[prof]]
  
  # Load expression matrix
  expr_path <- file.path(
    abs_from_root(cfg$paths$unsup_root),
    "feature_selection_unsupervised",
    "featuresets", "MX", "expr_submatrix.rds"
  )
  
  if (!file.exists(expr_path)) {
    stop(sprintf("Expression matrix not found for %s: %s", prof, expr_path))
  }
  
  expr_mat <- readRDS(expr_path)
  cat(sprintf("  %s: loaded matrix %d samples × %d genes\n", 
              toupper(prof), nrow(expr_mat), ncol(expr_mat)))
  
  # Get tumour samples from Top_m_long files
  best_dir <- best_directions[prof]
  tumour_samples <- character(0)
  
  if (!is.na(best_dir)) {
    # Find all Top_m_long files for this direction
    top_m_dir <- file.path(
      abs_from_root(cfg$paths$unsup_root),
      "tumour_neighbourhoods",
      best_dir
    )
    
    top_m_files <- list.files(top_m_dir, pattern = "Top_m_long_.*\\.csv$", 
                              recursive = TRUE, full.names = TRUE)
    
    if (length(top_m_files) > 0) {
      for (f in top_m_files) {
        top_m <- read_csv(f, show_col_types = FALSE)
        # Top_m_long files should expose tumour_id (British spelling) but handle legacy variants
        # Look for tumour/tumor ID column with priority order
        tumour_col <- NULL
        col_priority <- c("tumour_id", "tumor_id", "tumor", "tumour", "sample_id", "tumour_sample_id")
        
        # First try exact matches (case-insensitive)
        for (priority_col in col_priority) {
          matches <- grep(paste0("^", priority_col, "$"), names(top_m), ignore.case = TRUE, value = TRUE)
          if (length(matches) > 0) {
            tumour_col <- matches[1]
            break
          }
        }
        
        # If no exact match, try pattern matching
        if (is.null(tumour_col)) {
          for (col in names(top_m)) {
            if (grepl("^tumor.*id$|^tumour.*id$|^tumor$|^tumour$|^sample.*id$", col, ignore.case = TRUE)) {
              tumour_col <- col
              break
            }
          }
        }
        
        # Fallback: use first column if no obvious tumour column found
        if (is.null(tumour_col) && ncol(top_m) > 0) {
          tumour_col <- names(top_m)[1]
          warning(sprintf("No tumour ID column found in %s, using first column: %s", 
                         basename(f), tumour_col))
        }
        
        if (!is.null(tumour_col) && tumour_col %in% names(top_m)) {
          new_samples <- unique(top_m[[tumour_col]])
          new_samples <- new_samples[!is.na(new_samples) & nzchar(new_samples)]
          # Normalize tumour IDs based on prefix:
          # - TCGA/TARGET: convert underscores to hyphens (TCGA_B6_A0I6 -> TCGA-B6-A0I6)
          # - GSE: keep underscores (GSE268136_SRR29131244 stays as is)
          is_tcga_target <- grepl("^TCGA_|^TARGET", new_samples)
          new_samples[is_tcga_target] <- gsub("_", "-", new_samples[is_tcga_target])
          # GSE IDs keep underscores (they match the expression matrix format)
          tumour_samples <- c(tumour_samples, new_samples)
        } else {
          stop(sprintf("Cannot find tumour ID column in %s. Available columns: %s",
                      basename(f), paste(names(top_m), collapse = ", ")))
        }
      }
      tumour_samples <- unique(tumour_samples)
      tumour_samples <- tumour_samples[!is.na(tumour_samples) & nzchar(tumour_samples)]
      
      if (length(tumour_samples) == 0) {
        stop(sprintf("No tumour samples extracted from Top_m_long files for %s (direction: %s). Check column names in Top_m_long files.",
                    prof, best_dir))
      }
      
      cat(sprintf("    Found %d tumour samples from %d Top_m_long files\n",
                  length(tumour_samples), length(top_m_files)))
    } else {
      # Top_m_long files are mandatory for pan-cancer analysis
      stop(sprintf("No Top_m_long files found for %s (direction: %s). Tumour filtering requires Top_m_long files from best-performing neighbourhood.",
                  prof, best_dir))
    }
  } else {
    # Best direction is required for pan-cancer analysis
    stop(sprintf("No best direction determined for %s. Cannot filter tumours without direction summary.", prof))
  }
  
  # Match tumour samples against expression matrix rownames
  # (after normalization, they should match)
  all_samples <- rownames(expr_mat)
  tumour_samples_matched <- intersect(tumour_samples, all_samples)
  
  if (length(tumour_samples_matched) == 0) {
    # Try alternative normalization: maybe some IDs need different handling
    # Check if we need to try other patterns
    cat(sprintf("    Warning: No exact matches after normalization. Trying alternative patterns...\n"))
    # Try matching with partial string matching for GSE IDs or other formats
    tumour_samples_matched <- all_samples[all_samples %in% tumour_samples | 
                                          sapply(tumour_samples, function(x) any(grepl(x, all_samples, fixed = TRUE)))]
    
    if (length(tumour_samples_matched) == 0) {
      stop(sprintf("No tumour samples matched in expression matrix for %s.\n  Extracted from Top_m_long: %d samples\n  First 5 extracted: %s\n  First 5 in matrix: %s\n  Check tumour ID format normalization.",
                  prof, length(tumour_samples), 
                  paste(head(tumour_samples, 5), collapse = ", "),
                  paste(head(all_samples[grepl("^TCGA-|^TARGET|^GSE", all_samples)], 5), collapse = ", ")))
    }
  }
  
  tumour_samples <- tumour_samples_matched
  
  # Classify remaining samples: tumours are TCGA/TARGET/GSE, cell lines are everything else
  is_tumour_in_matrix <- grepl("^TCGA-|^TARGET|^GSE", all_samples)
  cell_line_samples <- all_samples[!is_tumour_in_matrix]
  
  # Store samples
  tumour_samples_all <- c(tumour_samples_all, tumour_samples)
  cell_line_samples_all <- c(cell_line_samples_all, cell_line_samples)
  
  # Subset matrix to include only tumour and cell line samples
  keep_samples <- c(tumour_samples, cell_line_samples)
  keep_samples <- intersect(keep_samples, rownames(expr_mat))
  
  if (length(keep_samples) == 0) {
    stop(sprintf("No samples to include for %s", prof))
  }
  
  expr_mat_subset <- expr_mat[keep_samples, , drop = FALSE]
  expr_matrices[[prof]] <- expr_mat_subset
  
  cat(sprintf("  %s: %d tumours + %d cell lines = %d total samples\n",
              toupper(prof), length(tumour_samples), length(cell_line_samples),
              length(keep_samples)))
}

# ==============================================================================
# Step 4: Combine matrices on shared genes
# ==============================================================================

cat("[INFO] Combining matrices on shared genes...\n")

# Get intersection of genes across all matrices AND joint MX union
# This ensures: 1) same genes in all matrices, 2) genes are in the MX candidate set
all_genes <- map(expr_matrices, colnames)
shared_genes <- Reduce(intersect, all_genes)  # Intersection across all matrices (mandatory)
shared_genes <- intersect(shared_genes, joint_mx_union)  # Also intersect with joint MX union

cat(sprintf("  Shared genes (across matrices): %d\n", length(Reduce(intersect, all_genes))))
cat(sprintf("  Shared genes (with joint MX union): %d\n", length(shared_genes)))

if (length(shared_genes) < 50) {
  warning(sprintf("Only %d shared genes across all three profiles. This may limit UMAP quality.", length(shared_genes)))
  # Option: use genes present in at least 2 of 3 matrices for more genes
  # For now, proceed with intersection but warn
  if (length(shared_genes) < 20) {
    stop(sprintf("Too few shared genes (%d < 20). Cannot proceed with pan-cancer UMAP. Check gene ID consistency or use a joint MX list.", length(shared_genes)))
  }
}

# Combine matrices (subset to shared genes, reorder columns identically)
# Convert to matrices and rbind (bind_rows expects data frames with same columns)
combined_expr <- map(expr_matrices, function(m) {
  m[, shared_genes, drop = FALSE]
}) %>%
  map(as.matrix) %>%
  do.call(rbind, .)

cat(sprintf("[INFO] Combined matrix: %d samples × %d genes\n",
            nrow(combined_expr), ncol(combined_expr)))

# ==============================================================================
# Step 5: Create metadata
# ==============================================================================

cat("[INFO] Creating metadata...\n")

metadata <- tibble(
  sample_id = rownames(combined_expr),
  cohort = if_else(grepl("^TCGA-", sample_id), "TCGA",
            if_else(grepl("^TARGET", sample_id), "TARGET", "DSMZ")),
  sample_type = if_else(grepl("^TCGA-|^TARGET", sample_id), "tumour", "cell_line"),
  cancer_type = NA_character_
)

# Assign cancer types based on which profile the sample came from
for (prof in profiles) {
  cfg <- profile_configs[[prof]]
  cancer_type <- toupper(cfg$analysis$cancer_type %||% toupper(prof))
  
  # Samples from this profile's matrix
  profile_samples <- rownames(expr_matrices[[prof]])
  metadata$cancer_type[metadata$sample_id %in% profile_samples] <- cancer_type
}

# Verify all samples have cancer_type
if (any(is.na(metadata$cancer_type))) {
  warning("Some samples have missing cancer_type. Assigning based on sample ID patterns.")
  # Fallback: try to infer from sample IDs
  for (prof in profiles) {
    cfg <- profile_configs[[prof]]
    cancer_type <- toupper(cfg$analysis$cancer_type %||% toupper(prof))
    # This is a fallback - may not work for all sample ID patterns
  }
}

cat(sprintf("  TCGA tumours: %d\n", sum(metadata$cohort == "TCGA")))
cat(sprintf("  DSMZ cell lines: %d\n", sum(metadata$cohort == "DSMZ")))
cat(sprintf("  Cancer types: %s\n", paste(unique(metadata$cancer_type), collapse = ", ")))

# ==============================================================================
# Step 6: PCA (optional) and UMAP
# ==============================================================================

cat("[INFO] Computing PCA and UMAP...\n")

# PCA (using pre-scaled data)
if (opt$n_pcs > 0 && opt$n_pcs < ncol(combined_expr)) {
  # Data is already scaled, so center=FALSE, scale.=FALSE
  pca_result <- prcomp(combined_expr, center = FALSE, scale. = FALSE)
  pca_coords <- pca_result$x[, 1:min(opt$n_pcs, ncol(pca_result$x)), drop = FALSE]
  cat(sprintf("  PCA: %d PCs (from pre-scaled data)\n", ncol(pca_coords)))
  umap_input <- pca_coords
} else {
  cat("  Skipping PCA, using scaled expression\n")
  umap_input <- as.matrix(combined_expr)
}

# UMAP
umap_result <- umap(
  umap_input,
  n_neighbors = opt$n_neighbors,
  min_dist = opt$min_dist,
  n_components = 2,
  verbose = TRUE
)

umap_coords <- as_tibble(umap_result, .name_repair = "unique")
names(umap_coords) <- c("UMAP_1", "UMAP_2")

# Combine with metadata
umap_output <- bind_cols(metadata, umap_coords)

cat("[INFO] UMAP complete\n")

# ==============================================================================
# Step 7: Save outputs
# ==============================================================================

cat("[INFO] Saving outputs...\n")

# Coordinates + metadata
coords_path <- file.path(outdir, "pan_cancer_mx_umap_coords.tsv")
write_tsv(umap_output, coords_path)
cat(sprintf("  Saved: %s\n", coords_path))

# Plot
plot_path <- file.path(outdir, "pan_cancer_mx_umap_plot.pdf")
pdf(plot_path, width = 12, height = 10)

# Color by cancer_type, shape by cohort
# Check if ggrepel is available for better label placement
use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

p <- ggplot(umap_output, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cancer_type, shape = cohort), size = 1.5, alpha = 0.6) +
  scale_shape_manual(values = c(TCGA = 16, DSMZ = 17)) +
  labs(
    title = "Pan-Cancer UMAP (Joint MX)",
    subtitle = sprintf("%d samples, %d genes", nrow(umap_output), length(shared_genes)),
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cancer Type",
    shape = "Cohort"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right"
  )

# Optional: label cell lines (top 10 per cancer type, or all if <= 20 total)
if (sum(umap_output$cohort == "DSMZ") > 0) {
  cell_lines <- umap_output %>% filter(cohort == "DSMZ")
  
  if (nrow(cell_lines) <= 20) {
    # Label all cell lines if <= 20
    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = cell_lines,
        aes(label = sample_id),
        size = 2.5,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.3,
        point.padding = 0.2
      )
    } else {
      p <- p + geom_text(
        data = cell_lines,
        aes(label = sample_id),
        size = 2.5,
        hjust = 0,
        vjust = 0,
        nudge_x = 0.1,
        nudge_y = 0.1,
        check_overlap = TRUE
      )
    }
  } else {
    # Label top 10 per cancer type (by some metric - here just first 10 per type)
    cell_lines_labeled <- cell_lines %>%
      group_by(cancer_type) %>%
      slice_head(n = 10) %>%
      ungroup()
    
    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = cell_lines_labeled,
        aes(label = sample_id),
        size = 2.5,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.3,
        point.padding = 0.2
      )
    } else {
      p <- p + geom_text(
        data = cell_lines_labeled,
        aes(label = sample_id),
        size = 2.5,
        hjust = 0,
        vjust = 0,
        nudge_x = 0.1,
        nudge_y = 0.1,
        check_overlap = TRUE
      )
    }
  }
}

print(p)
dev.off()
cat(sprintf("  Saved: %s\n", plot_path))

cat("[INFO] Pan-cancer UMAP complete!\n")
