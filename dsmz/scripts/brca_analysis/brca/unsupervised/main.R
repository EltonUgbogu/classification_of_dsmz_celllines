# main.R
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dynamicTreeCut)
  library(dbscan)
  library(uwot)
})

# --- Source your functions ---
src <- function(...) lapply(file.path("R", c(...)), source, chdir = TRUE)
src("helpers.R",
    "clustering_fixed.R",
    "clustering_dynamic.R",
    "clustering_hdbscan.R",
    "labeling_tcga.R")

# --- Load inputs (fail early if missing) ---
req <- c("data/PC.rds", "data/V_adj.rds", "data/V_tcga_brca.rds", "data/PAM50_genes.rds", "data/tcga_brca_metadata.csv")
miss <- req[!file.exists(req)]
if (length(miss)) stop(sprintf("Missing required file(s):\n- %s", paste(miss, collapse = "\n- ")))

PC       <- readRDS("data/PC.rds")          # matrix n_samples x n_PCs
V_adj    <- readRDS("data/V_adj.rds")       # matrix n_genes x n_samples
V_tcga   <- readRDS("data/V_tcga_brca.rds") # matrix n_genes x n_tcga
features <- readRDS("data/PAM50_genes.rds") # character vector length ~50
meta_ann <- "data/tcga_brca_metadata.csv"   # has sample_barcode + PAM50_Subtype

# --- Ensure names & alignment ---
if (is.null(rownames(PC))) rownames(PC) <- paste0("sample_", seq_len(nrow(PC)))
stopifnot(identical(colnames(V_adj), rownames(PC)))  # PC rows correspond to V_adj cols

# --- Output folders ---
out_fixed   <- "results/fixed_hc_k4"
out_dynamic <- "results/dynamic_cut"
out_hdb     <- "results/hdbscan_iter"
out_label   <- "results/tcga_labeling"

dir.create("results", showWarnings = FALSE)
dir.create(out_fixed,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dynamic, recursive = TRUE, showWarnings = FALSE)
dir.create(out_hdb,     recursive = TRUE, showWarnings = FALSE)
dir.create(out_label,   recursive = TRUE, showWarnings = FALSE)

# ========== A) Fixed-k hierarchical clustering ==========
res_fixed <- run_fixed_hc(
  PC,
  k = 4,
  dist_method = "euclidean",
  linkage = "average",
  seed = 42,
  outdir = out_fixed
)
fixed_clusters <- res_fixed$clusters

# ========== B) Dynamic Tree Cut ==========
res_dyn <- run_dynamic_cut(
  PC,
  min_cluster_size = NULL,   # uses choose_minCluster() if present (else heuristic)
  deep_split = 2,
  pam_stage = TRUE,
  dist_method = "euclidean",
  linkage = "average",
  outdir = out_dynamic,
  save_unassigned = TRUE
)
dynamic_clusters <- res_dyn$clusters  # factor of integer labels

# ========== C) Iterative HDBSCAN ==========
labs <- paste0("C", 1:4)
hdb_clusters <- run_iterative_hdbscan(
  PC,
  cluster_labels = labs,
  seed = 42,
  metric = "cosine",
  umap_neighbors_first = 20,
  umap_neighbors_next  = 15,
  umap_min_dist = 0.3,
  outdir = out_hdb,
  metadata = NULL,
  extra_cols = NULL
)

# ========== D) Supervised labeling (TCGA PAM50) ==========
labeled_subtypes <- label_by_tcga(
  clusters     = hdb_clusters,  # or fixed_clusters / dynamic_clusters
  V_adj        = V_adj,
  features     = features,
  metadata_ann = meta_ann,
  V_tcga       = V_tcga,
  subtype      = "PAM50_Subtype",
  outdir       = out_label,
  write_by_cluster = TRUE,
  write_by_subtype = TRUE
)

cat("\n[SUMMARY]\n",
    sprintf("Fixed HC clusters:      %d\n", length(unique(as.integer(fixed_clusters)))),
    sprintf("Dynamic cut clusters:    %d\n", length(unique(as.integer(dynamic_clusters)))),
    sprintf("HDBSCAN labels (k=%d):   %s\n", length(labs), paste(levels(hdb_clusters), collapse = ", ")),
    sprintf("Mapped subtypes:         %s\n", paste(levels(labeled_subtypes), collapse = ", ")),
    sep = "")

message("[DONE] All steps completed.")
