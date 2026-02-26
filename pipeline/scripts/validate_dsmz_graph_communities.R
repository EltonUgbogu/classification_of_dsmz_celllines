#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(optparse)
})

# Helper: compute Shannon entropy manually (no external package needed)
shannon_entropy <- function(p) {
  p <- p[p > 0]  # remove zeros
  if (length(p) == 0) return(0)
  -sum(p * log2(p))
}

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--graph_dir", type = "character", default = NULL,
              help = "Directory containing DSMZ_DSMZ_graph_node_annotations_*.tsv"),
  make_option("--meta_tsv", type = "character", 
              default = "results/unsupervised/pan_cancer/inputs/joint_metadata.tsv",
              help = "Path to joint_metadata.tsv [default: %default]"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (e.g., MAD_euclidean) - auto-detected if not provided")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$graph_dir)) {
  stop("Please provide --graph_dir")
}

if (!dir.exists(opt$graph_dir)) {
  stop("Graph directory not found: ", opt$graph_dir)
}

# Auto-detect direction from graph_dir if not provided
if (is.null(opt$direction)) {
  # Extract direction from path (e.g., .../MAD_euclidean/...)
  path_parts <- strsplit(opt$graph_dir, "/")[[1]]
  direction_idx <- which(path_parts %in% c("global", "within"))
  if (length(direction_idx) > 0 && direction_idx < length(path_parts)) {
    opt$direction <- path_parts[direction_idx + 1]
  } else {
    stop("Cannot auto-detect direction from path. Please provide --direction")
  }
}

cat("[INFO] Graph directory:", opt$graph_dir, "\n")
cat("[INFO] Direction:", opt$direction, "\n")
cat("[INFO] Metadata:", opt$meta_tsv, "\n")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
nodes_file <- file.path(opt$graph_dir, sprintf("DSMZ_DSMZ_graph_node_annotations_%s.tsv", opt$direction))

if (!file.exists(nodes_file)) {
  stop("Node annotations file not found: ", nodes_file)
}

nodes <- read_tsv(nodes_file, show_col_types = FALSE)
cat("[INFO] Loaded", nrow(nodes), "nodes\n")

if (!file.exists(opt$meta_tsv)) {
  stop("Metadata file not found: ", opt$meta_tsv)
}

meta <- read_tsv(opt$meta_tsv, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

# Normalize IDs (strip CELL:/TUMOUR: prefixes)
strip_prefix <- function(x) sub("^(CELL:|TUMOUR:|TUMOR:)", "", x)
meta <- meta %>%
  mutate(sample_id = strip_prefix(sample_id))

# Merge disease labels
df <- nodes %>%
  mutate(cell_line_clean = strip_prefix(as.character(cell_line))) %>%
  left_join(
    meta %>% 
      filter(sample_type == "Cell Line", cohort == "DSMZ") %>%
      select(cell_line = sample_id, cancer_type),
    by = c("cell_line_clean" = "cell_line")
  ) %>%
  select(-cell_line_clean)

if (any(is.na(df$cancer_type))) {
  missing <- sum(is.na(df$cancer_type))
  cat("[WARNING]", missing, "cell lines missing cancer_type labels\n")
  # Try to infer from cell_line name patterns if needed
  df <- df %>%
    mutate(
      cancer_type = ifelse(
        is.na(cancer_type),
        # Fallback: try to extract from sample_id patterns
        case_when(
          grepl("BRCA|BT_474|CAL_|MCF|T47D|SKBR|HCC|MDA", cell_line, ignore.case = TRUE) ~ "BRCA",
          grepl("NBL|SK|BE", cell_line, ignore.case = TRUE) ~ "NBL",
          grepl("RBL|RB|Y79", cell_line, ignore.case = TRUE) ~ "RBL",
          TRUE ~ "UNKNOWN"
        ),
        cancer_type
      )
    )
}

cat("[INFO] Merged disease labels:", length(unique(df$cancer_type[!is.na(df$cancer_type)])), "diseases\n")

# ------------------------------------------------------------------------------
# Helper: compute purity + entropy per community
# ------------------------------------------------------------------------------
summarise_community <- function(data, comm_col) {
  data %>%
    group_by(.data[[comm_col]], cancer_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[comm_col]]) %>%
    mutate(
      frac = n / sum(n)
    ) %>%
    summarise(
      community = first(.data[[comm_col]]),
      size = sum(n),
      top_label = cancer_type[which.max(frac)],
      top_frac = max(frac),
      entropy = shannon_entropy(frac),
      n_diseases = n_distinct(cancer_type),
      .groups = "drop"
    ) %>%
    arrange(desc(size))
}

# ------------------------------------------------------------------------------
# Compute summaries for both methods
# ------------------------------------------------------------------------------
cat("\n=== Computing Louvain community summary ===\n")
louv_summary <- summarise_community(df, "community_louv") %>%
  mutate(method = "Louvain")

cat("\n=== Computing Leiden community summary ===\n")
leid_summary <- summarise_community(df, "community_leid") %>%
  mutate(method = "Leiden")

final_summary <- bind_rows(louv_summary, leid_summary) %>%
  arrange(method, desc(size))

# ------------------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------------------
out_file <- file.path(opt$graph_dir, sprintf("DSMZ_DSMZ_graph_community_validation_%s.tsv", opt$direction))
write_tsv(final_summary, out_file)
cat("\n[SUCCESS] Community validation table written:\n  ", out_file, "\n")

# ------------------------------------------------------------------------------
# Print summary statistics
# ------------------------------------------------------------------------------
cat("\n=== Summary Statistics ===\n")
cat("\nLouvain communities:\n")
print(louv_summary %>% 
  summarise(
    n_communities = n(),
    mean_size = mean(size),
    mean_purity = mean(top_frac),
    mean_entropy = mean(entropy),
    high_purity = sum(top_frac >= 0.7),
    low_entropy = sum(entropy <= 0.6)
  ))

cat("\nLeiden communities:\n")
print(leid_summary %>% 
  summarise(
    n_communities = n(),
    mean_size = mean(size),
    mean_purity = mean(top_frac),
    mean_entropy = mean(entropy),
    high_purity = sum(top_frac >= 0.7),
    low_entropy = sum(entropy <= 0.6)
  ))

cat("\n=== Top communities by size (Louvain) ===\n")
print(head(louv_summary, 10))

cat("\n=== Top communities by size (Leiden) ===\n")
print(head(leid_summary, 10))

cat("\n[INFO] Validation complete.\n")
cat("[INFO] Trust criteria: size >= 5, top_frac >= 0.7, entropy <= 0.6\n")

