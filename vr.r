Got it — here’s your **fully corrected, simplified** script. I added a tiny `logCPM()` helper and replaced the batch-adjustment logic with a minimal, robust version:

* Supports `config$batch_adjust_method`: `"none"`, `"combat_seq"`, or `"combat"`.
* `ComBat_seq` runs **without** `group` (avoids confounding entirely).
* If a batch has fewer than 2 samples, it **auto-falls back** to `ComBat` on logCPM.
* No extra tables/CSV exports or verbose checks.

Paste this over your current Rmd:

````r
---
title: "TCGA–DSMZ Correlation Analysis"
author: "chu25"
date: "`r format(Sys.Date())`"
output:
  pdf_document:
    toc: true
    number_sections: true
geometry: margin=1in
fontsize: 11pt
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
set.seed(42)
# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(matrixStats)
  library(limma)
  library(edgeR)
  library(pheatmap)
  library(uwot)
  library(ggrepel)
  library(sva)         # ComBat / ComBat_seq
  library(patchwork)   # for combining plots
  library(viridis)     # for color palettes
  library(ggpubr)      # for statistical tests
  library(rstatix)     # for statistical tests
  library(grid)        # for pheatmap drawing
  library(RColorBrewer)
})

# Helpers
`%||%` <- function(x, y) if (is.null(x)) y else x

# Log2-CPM helper using edgeR::cpm (prior.count avoids -Inf for zeros)
logCPM <- function(x, prior.count = 1, lib.size = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x)
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size)
}

set.seed(42)
````

```{r tcga_proj_map}
# Human-readable TCGA project names for labeling 
tcga_proj_map <- tibble::tribble(
  ~tcga_code, ~tcga_name,
  "LAML","Acute Myeloid Leukemia","ACC","Adrenocortical carcinoma",
  "BLCA","Bladder Urothelial Carcinoma","LGG","Brain Lower Grade Glioma",
  "BRCA","Breast invasive carcinoma","CESC","Cervical squamous cell carcinoma and endocervical adenocarcinoma",
  "CHOL","Cholangiocarcinoma","LCML","Chronic Myelogenous Leukemia",
  "COAD","Colon adenocarcinoma","CNTL","Controls","ESCA","Esophageal carcinoma",
  "FPPP","FFPE Pilot Phase II","GBM","Glioblastoma multiforme","HNSC","Head and Neck squamous cell carcinoma",
  "KICH","Kidney Chromophobe","KIRC","Kidney renal clear cell carcinoma","KIRP","Kidney renal papillary cell carcinoma",
  "LIHC","Liver hepatocellular carcinoma","LUAD","Lung adenocarcinoma","LUSC","Lung squamous cell carcinoma",
  "DLBC","Lymphoid Neoplasm Diffuse Large B-cell Lymphoma","MESO","Mesothelioma","MISC","Miscellaneous",
  "OV","Ovarian serous cystadenocarcinoma","PAAD","Pancreatic adenocarcinoma","PCPG","Pheochromocytoma and Paraganglioma",
  "PRAD","Prostate adenocarcinoma","READ","Rectum adenocarcinoma","SARC","Sarcoma",
  "SKCM","Skin Cutaneous Melanoma","STAD","Stomach adenocarcinoma","TGCT","Testicular Germ Cell Tumors",
  "THYM","Thymoma","THCA","Thyroid carcinoma","UCS","Uterine Carcinosarcoma",
  "UCEC","Uterine Corpus Endometrial Carcinoma","UVM","Uveal Melanoma"
)

config <- list(
  # Inputs
  tcga_se_rds   = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds",
  dsmz_rds      = "/home/chu25/data/dsmz/DSMZ_count_gene.rds",
  dsmz_meta_csv = "/home/chu25/data/dsmz/DSMZ_metadata.csv",
  purity_tsv    = NULL,  # set TSV with columns: sample, purity if SE lacks purity

  # Outputs
  outdir       = "/home/chu25/dsmz/results/correlation2",

  # Analysis
  var_gene_set = "tcga_5k",        # "tcga_5k", "tcga_10k", "all"
  min_shared   = 1000,
  batch_adjust_method = "combat_seq",  # "none", "combat_seq", or "combat"
  overlay_tcga_centroids = TRUE
)

dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)
knitr::opts_chunk$set(fig.path = paste0(config$outdir, "/figs-"))
```

```{r io_functions}
# =============================================================================
# IO HELPERS
# =============================================================================
load_tcga_data <- function(file_path) {
  cat("[INFO] Loading TCGA data...\n")
  if (!file.exists(file_path)) stop(sprintf("[ERROR] TCGA RDS not found: %s", file_path))
  se <- readRDS(file_path)
  counts <- assay(se)
  cat(sprintf("[INFO] TCGA dims: %d genes x %d samples\n", nrow(counts), ncol(counts)))
  list(se = se, counts = counts)
}

load_dsmz_data <- function(counts_path, meta_path) {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) stop(sprintf("[ERROR] DSMZ RDS not found: %s", counts_path))
  if (!file.exists(meta_path))    stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  raw  <- readRDS(counts_path)
  meta <- read.csv(meta_path)

  stopifnot(all(c("Ensembl_ID","gene_name") %in% colnames(raw)))
  stopifnot("sample_name" %in% names(meta))
  cat(sprintf("[INFO] DSMZ table: %d rows x %d cols\n", nrow(raw), ncol(raw)))
  list(raw = raw, meta = meta)
}

# =============================================================================
# PROCESS: counts + metadata
# =============================================================================
harmonize_gene_ids <- function(tcga_counts){
  cat("[INFO] Harmonizing TCGA gene IDs...\n")
  g <- sub("\\..*$","", rownames(tcga_counts))
  if (any(duplicated(g))) {
    tcga_counts <- rowsum(tcga_counts, g, reorder = TRUE)
    rownames(tcga_counts) <- sort(unique(g))
  } else {
    rownames(tcga_counts) <- g
  }
  tcga_counts
}

# Build DSMZ count matrix from the raw table
build_dsmz_matrix <- function(dsmz_raw){
  cat("[INFO] Building DSMZ count matrix...\n")
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot)
  dsmz_raw[sample_cols] <- lapply(
    dsmz_raw[sample_cols],
    function(x) as.numeric(as.character(x))
  )
  M <- as.matrix(dsmz_raw[, sample_cols, drop=FALSE])
  rownames(M) <- dsmz_raw$Ensembl_ID
  if (any(duplicated(rownames(M))))
    M <- rowsum(M, rownames(M), reorder = TRUE)
  keep <- vapply(as.data.frame(M), function(v) any(!is.na(v)), logical(1))
  if (!all(keep)) {
    cat(sprintf("[WARN] Dropping %d DSMZ columns that are all NA\n", sum(!keep)))
    M <- M[, keep, drop=FALSE]
  }
  M
}

# Align DSMZ metadata with DSMZ counts
align_metadata <- function(dsmz_meta, dsmz_counts, outdir){
  cat("[INFO] Aligning metadata to counts...\n")
  dsmz_meta <- dsmz_meta %>% mutate(sample_id = sample_name)
  matched <- intersect(dsmz_meta$sample_id, colnames(dsmz_counts))
  md_only <- setdiff(dsmz_meta$sample_id, colnames(dsmz_counts))
  ct_only <- setdiff(colnames(dsmz_counts), dsmz_meta$sample_id)
  if (length(md_only))
    write.csv(tibble(sample_name = md_only),
              file.path(outdir, "unmatched_metadata_samples.csv"),
              row.names = FALSE)
  if (length(ct_only))
    write.csv(tibble(sample_name = ct_only),
              file.path(outdir, "unmatched_count_columns.csv"),
              row.names = FALSE)
  dsmz_meta <- dsmz_meta %>% filter(sample_id %in% matched)
  dsmz_counts <- dsmz_counts[, dsmz_meta$sample_id, drop=FALSE]
  stopifnot(identical(colnames(dsmz_counts), dsmz_meta$sample_id))
  cat(sprintf("[INFO] Matched: %d | Unmatched(meta): %d | Unmatched(counts): %d\n",
              length(matched), length(md_only), length(ct_only)))
  list(meta = dsmz_meta, counts = dsmz_counts)
}
```

```{r purity_functions}
# =============================================================================
# PURITY 
# =============================================================================
load_purity_data <- function(tcga_se, purity_file=NULL){
  p <- NULL
  if ("purity" %in% colnames(colData(tcga_se))) {
    p <- as.numeric(colData(tcga_se)$purity)
    names(p) <- colnames(assay(tcga_se))
    cat("[INFO] Using purity from SummarizedExperiment\n")
  } else if (!is.null(purity_file) && file.exists(purity_file)) {
    tab <- read.delim(purity_file, stringsAsFactors=FALSE)
    stopifnot(all(c("sample","purity") %in% colnames(tab)))
    p <- setNames(tab$purity, tab$sample)
    cat("[INFO] Using external purity file\n")
  } else {
    cat("[INFO] No purity data available\n")
  }
  p
}
```

```{r adjust_purity}
apply_purity_adjustment <- function(tcga_counts, purity){
  if (is.null(purity)) return(tcga_counts)
  cat("[INFO] Applying purity adjustment...\n")
  keep <- colnames(tcga_counts)[!is.na(purity)]
  purity <- purity[keep]
  M <- tcga_counts[, keep, drop=FALSE]
  ct <- apply(M, 1, function(v) suppressWarnings(
    cor.test(as.numeric(v), purity, method="spearman")
  ))
  pv <- vapply(ct, `[[`, numeric(1), "p.value")
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))
  padj <- p.adjust(pv, "BH")
  drop <- names(which(padj < 0.01 & rh < -0.4))
  if (length(drop)) {
    cat(sprintf("[INFO] Removing %d purity-associated genes\n", length(drop)))
    M <- M[setdiff(rownames(M), drop), , drop=FALSE]
  }
  infilt <- 1 - purity
  design <- model.matrix(~ infilt)
  fit <- lmFit(M, design)
  beta <- fit$coefficients[,2, drop=FALSE]
  Madj <- as.matrix(M) - beta %*% t(infilt)
  out <- tcga_counts
  common <- intersect(rownames(Madj), rownames(out))
  out[common, colnames(Madj)] <- Madj[common, ]
  out
}
```

```{r compute_tcga_means}
# =============================================================================
# TCGA labels & means (for correlation step)
# =============================================================================
compute_tcga_means <- function(tcga_se, tcga_adj){
  labs <- NULL
  for (nm in c("project_id","study","disease_type")) 
    if (nm %in% colnames(colData(tcga_se))) { 
      labs <- colData(tcga_se)[[nm]] 
      break 
    }
  stopifnot(!is.null(labs))
  labs <- sub("^TCGA-","", as.character(labs))
  names(labs) <- colnames(tcga_adj)
  lev <- sort(unique(labs))
  means <- sapply(
    lev, 
    function(ct) rowMeans(tcga_adj[, labs==ct, drop=FALSE])
  )
  colnames(means) <- lev
  list(means=means, labels=labs)
}
```

```{r select_variable_genes}
select_variable_genes <- function(mat, var_gene_set){
  sel <- rownames(mat)
  if (tolower(var_gene_set) %in% c("tcga_5k","tcga_10k")){
    n <- if (tolower(var_gene_set)=="tcga_5k") 5000 else 10000
    iq <- matrixStats::rowIQRs(mat)
    n <- min(n, length(iq))
    sel <- names(sort(setNames(iq, rownames(mat)), decreasing=TRUE))[seq_len(n)]
  }
  cat(sprintf("[INFO] Variable-gene setting: %s (%d genes used)\n", var_gene_set, length(sel)))
  sel
}
```

```{r compute_correlations}
compute_correlations <- function(dsmz_counts, tcga_means, dsmz_meta, sel_genes){
  cat("[INFO] Computing Spearman correlations...\n")
  spearman <- outer(
    dsmz_meta$sample_id,
    colnames(tcga_means),
    Vectorize(function(cell, cohort) {
      suppressWarnings(
        cor(
          dsmz_counts[sel_genes, cell],
          tcga_means[sel_genes, cohort],
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    })
  )
  dimnames(spearman) <- list(dsmz_meta$sample_id, colnames(tcga_means))
  spearman
}
```

```{r get_organ_mapping}
# =============================================================================
# ORGAN ↔ TCGA mapping (uses your 'organ' column)
# =============================================================================
get_organ_tcga_mapping <- function(){
  list(
    "Adrenal"       = c("ACC","PCPG"),
    "Bladder"       = "BLCA",
    "Brain/CNS"     = c("GBM","LGG"),
    "Breast"        = "BRCA",
    "Cervix"        = "CESC",
    "Liver/Biliary" = c("LIHC","CHOL"),
    "Colon/Rectum"  = c("COAD","READ"),
    "Esophagus"     = "ESCA",
    "Head/Neck"     = "HNSC",
    "Kidney"        = c("KICH","KIRC","KIRP"),
    "Lung"          = c("LUAD","LUSC"),
    "Mesothelium"   = "MESO",
    "Ovary"         = "OV",
    "Pancreas"      = "PAAD",
    "Prostate"      = "PRAD",
    "Sarcoma"       = "SARC",
    "Skin"          = "SKCM",
    "Stomach"       = "STAD",
    "Testis"        = "TGCT",
    "Thyroid"       = "THCA",
    "Thymus"        = "THYM",
    "Uterus"        = c("UCEC","UCS"),
    "Uveal/Eye"     = "UVM",
    "Lymphoid"      = "DLBC",
    "Leukemia"      = "LAML",
    "Myeloid"       = "LAML",
    "Bone Marrow"   = "LAML",
    "Body Cavity"   = "MESO",
    "Hematologic"   = c("DLBC","LAML"),
    "Adrenal/SNS"   = c("ACC","PCPG"),
    "Unknown"       = character(0),
    "Other"         = character(0)
  )
}
```

```{r perform_organ_mapping}
perform_organ_mapping <- function(dsmz_meta, spearman_scores, tcga_means) {
  cat("[INFO] Performing organ-constrained mapping (robust)…\n")
  organ2tcga <- get_organ_tcga_mapping()
  corr_rows <- rownames(spearman_scores)
  corr_cols <- colnames(spearman_scores)
  tcga_codes_available <- intersect(colnames(tcga_means), corr_cols)
  dsmz_meta %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      organ     = ifelse(is.na(organ) | organ == "", "Other", trimws(as.character(organ)))
    ) %>%
    rowwise() %>%
    mutate(
      tcga_candidates = list({
        cands0 <- organ2tcga[[organ]]
        if (is.null(cands0)) character(0) else intersect(cands0, tcga_codes_available)
      }),
      tcga_code = {
        sid <- sample_id
        cand <- tcga_candidates
        if (!is.character(cand)) cand <- as.character(cand)
        if (is.na(sid) || !(sid %in% corr_rows) || length(cand) == 0L) {
          NA_character_
        } else {
          ridx <- which(corr_rows == sid)[1]
          cand <- intersect(cand, corr_cols)
          if (length(cand) == 0L) {
            NA_character_
          } else {
            sc <- spearman_scores[ridx, cand, drop = TRUE]
            sc[!is.finite(sc)] <- -Inf
            cand[which.max(sc)]
          }
        }
      }
    ) %>%
    ungroup()
}
```

```{r tcga_project}
tcga_project_to_organ_map <- function(){
  get_organ_tcga_mapping() %>%               
    enframe(name = "organ", value = "codes") %>%
    unnest_longer(codes, values_to = "tcga_project") %>%
    filter(!is.na(tcga_project)) %>%
    distinct(tcga_project, .keep_all = TRUE) %>%
    select(tcga_project, organ)
}
```

```{r plotting}
# =============================================================================
# - COLOR PALETTES AND PLOTTING
# =============================================================================
get_palettes <- function() {
  list(
    okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                  "#0072B2", "#D55E00", "#CC79A7", "#000000"),
    tol_bright = c("#4477AA", "#EE6677", "#228833", "#CCBB44", 
                   "#66CCEE", "#AA3377", "#BBBBBB"),
    tol_vibrant = c("#EE7733", "#0077BB", "#33BBEE", "#EE3377", 
                    "#CC3311", "#009988", "#BBBBBB"),
    two_group = c("#E69F00", "#0072B2"),
    sequential = c("#FFFFFF", "#FFF2CC", "#FFE699", "#FFD966", 
                   "#FFCC33", "#FFB300", "#E69F00", "#CC8800")
  )
}

create_plots<- function(summary_tbl, spearman_scores, dsmz_map, outdir) {
  cat("[INFO] Creating correlation plots...\n")
  palettes <- get_palettes()

  mean_ci95 <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    m <- mean(x)
    se <- stats::sd(x) / sqrt(n)
    data.frame(y = m, ymin = m - 1.96*se, ymax = m + 1.96*se)
  }

  fig_dir <- file.path(outdir, "figs")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  spearman_long <- spearman_scores %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    tidyr::pivot_longer(-sample_id, names_to = "tcga_project", values_to = "rho") %>%
    dplyr::left_join(dsmz_map %>% dplyr::select(sample_id, organ), by = "sample_id") %>%
    dplyr::mutate(organ = dplyr::coalesce(organ, "Other"))

  best_corr_data <- spearman_long %>%
    dplyr::group_by(sample_id) %>%
    dplyr::slice_max(order_by = rho, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  organ_summary <- best_corr_data %>%
    dplyr::group_by(organ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(rho, na.rm = TRUE),
      sd = sd(rho, na.rm = TRUE),
      se = sd/sqrt(n),
      ci95_low = mean - 1.96*se,
      ci95_high = mean + 1.96*se,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean))

  readr::write_csv(organ_summary, file.path(fig_dir, "violin_rho_by_organ_summary.csv"))

  kw <- rstatix::kruskal_test(best_corr_data, rho ~ organ)
  dunn <- rstatix::dunn_test(best_corr_data, rho ~ organ, p.adjust.method = "BH") %>%
    rstatix::add_significance("p.adj")
  dunn_sig <- dunn %>% dplyr::filter(p.adj < 0.05)
  if (nrow(dunn_sig) > 0) {
    dunn_sig <- rstatix::add_xy_position(dunn_sig, x = "organ", data = best_corr_data, step.increase = 0.06)
  }

  organs <- sort(unique(best_corr_data$organ))
  n_organs <- length(organs)
  if (n_organs <= 8) {
    organ_colors <- palettes$okabe_ito[1:n_organs]
  } else {
    organ_colors <- c(palettes$okabe_ito, palettes$tol_bright, palettes$tol_vibrant)[1:n_organs]
  }
  names(organ_colors) <- organs

  subtitle_txt <- paste0("Kruskal–Wallis p = ", formatC(kw$p, format = "e", digits = 2))

  p1 <- ggplot(best_corr_data, aes(x = organ, y = rho)) +
    geom_violin(aes(fill = organ), trim = FALSE, alpha = 0.6, color = "black", size = 0.5) +
    geom_jitter(aes(color = organ, shape = organ), width = 0.15, alpha = 0.7, size = 1.5) +
    stat_summary(fun.data = mean_ci95, geom = "pointrange", color = "black", size = 0.8, fatten = 3) +
    scale_fill_manual(values = organ_colors, guide = "none") +
    scale_color_manual(values = organ_colors, guide = "none") +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 8, 4, 3, 7), length.out = n_organs), guide = "none") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Organ/Tissue Type", 
      y = "Best Spearman ρ",
      title = "Distribution of Best DSMZ→TCGA Correlations by Organ/Tissue",
      subtitle = subtitle_txt
    )
  if (nrow(dunn_sig) > 0) {
    p1 <- p1 + ggpubr::stat_pvalue_manual(
      dunn_sig, label = "p.adj.signif", hide.ns = TRUE,
      tip.length = 0.01, size = 3.5
    )
  }
  print(p1)
  ggsave(file.path(fig_dir, "violin_rho_by_organ.pdf"), p1, width = 12, height = 7)

  # Heatmap of top 25
  topN <- head(dplyr::arrange(summary_tbl, dplyr::desc(best_rho_overall))$sample_id, 25)
  hm <- spearman_scores[topN, , drop = FALSE]
  colors_div <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  hp <- pheatmap::pheatmap(
    hm, cluster_rows = TRUE, cluster_cols = TRUE,
    color = colorRampPalette(colors_div)(100),
    main = "Top 25 DSMZ Samples: Spearman Correlation with TCGA Cohorts",
    fontsize = 9, fontsize_row = 7, fontsize_col = 8,
    border_color = "grey60", silent = TRUE
  )
  grid::grid.newpage(); grid::grid.draw(hp$gtable)
  pdf(file.path(fig_dir, "heatmap_top25.pdf"), width = 12, height = 10)
  grid::grid.newpage(); grid::grid.draw(hp$gtable)
  dev.off()

  # Stacked bar: DSMZ organ vs assigned TCGA cohort
  assignment_data <- dsmz_map %>%
    dplyr::mutate(tcga_assignment = ifelse(is.na(tcga_code) | tcga_code == "", "Unassigned", tcga_code)) %>%
    dplyr::count(organ, tcga_assignment, name = "count")
  totals <- assignment_data %>% dplyr::group_by(organ) %>% dplyr::summarise(total = sum(count), .groups = "drop") %>% dplyr::arrange(dplyr::desc(total))
  organ_order <- totals$organ
  assignment_data$organ <- factor(assignment_data$organ, levels = organ_order)
  totals$organ          <- factor(totals$organ,          levels = organ_order)
  assignment_labels <- assignment_data %>%
    dplyr::group_by(organ) %>%
    dplyr::arrange(organ, tcga_assignment, .by_group = TRUE) %>%
    dplyr::mutate(pos = cumsum(count) - count/2) %>%
    dplyr::ungroup()
  tcga_assignments <- sort(unique(assignment_data$tcga_assignment))
  n_tcga <- length(tcga_assignments)
  tcga_colors <- c(palettes$okabe_ito, palettes$tol_bright, palettes$tol_vibrant)[1:n_tcga]
  names(tcga_colors) <- tcga_assignments
  label_df <- dplyr::filter(assignment_labels, count >= 2)
  assignment_data <- assignment_data %>% dplyr::mutate(tcga_assignment = factor(tcga_assignment, levels = tcga_assignments))
  label_df <- label_df %>% dplyr::mutate(tcga_assignment = factor(tcga_assignment, levels = tcga_assignments))

  p3 <- ggplot(assignment_data, aes(x = organ, y = count, fill = tcga_assignment)) +
    geom_col(color = "grey15", linewidth = 0.2) +
    geom_text(data = label_df, inherit.aes = FALSE, aes(x = organ, y = pos, label = count), size = 3, color = "black") +
    geom_text(data = totals, inherit.aes = FALSE, aes(x = organ, y = total, label = total), vjust = -0.4, fontface = "bold", size = 3.2) +
    scale_fill_manual(values = tcga_colors, name = "TCGA Assignment") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    theme_bw(base_size = 11) +
    labs(x = "Organ/Tissue Type", y = "Number of Cell Lines", title = "DSMZ Cell Line Assignments to TCGA Projects by Organ/Tissue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  print(p3)
  ggsave(file.path(fig_dir, "organ_tcga_assignment.pdf"), p3, width = 14, height = 9)

  readr::write_csv(kw,   file.path(fig_dir, "violin_kw_results.csv"))
  readr::write_csv(dunn, file.path(fig_dir, "violin_dunn_results.csv"))
  cat("[INFO] Plots saved to: ", fig_dir, "\n")
}
```

```{r pca}
# =============================================================================
#  PCA 
# =============================================================================
pc_dataset_s_ <- function(M_log, dataset_factor, top_n = 10, outdir = NULL, title_tag = "dataset") {
  if (!is.matrix(M_log) && !is.data.frame(M_log)) stop("M_log must be a matrix or data frame")
  if (ncol(M_log) != length(dataset_factor)) stop("Length mismatch between M_log columns and dataset_factor")
  palettes <- get_palettes()
  pc <- prcomp(t(M_log), scale. = TRUE)
  var_explained <- (pc$sdev^2) / sum(pc$sdev^2)
  if (is.null(top_n) || top_n > ncol(pc$x)) {
    cum_var <- cumsum(var_explained)
    top_n <- which.max(cum_var >= 0.9)
    top_n <- min(top_n, ncol(pc$x))
  }
  get_r2 <- function(y) unname(summary(lm(y ~ dataset_factor))$r.squared)
  pcs_to_check <- seq_len(min(top_n, ncol(pc$x)))
  r2_values <- sapply(pcs_to_check, function(i) get_r2(pc$x[, i]))
  df_r2 <- data.frame(
    PC = paste0("PC", pcs_to_check),
    R2 = r2_values,
    VariancePct = var_explained[pcs_to_check] * 100
  )
  df_r2$PC <- factor(df_r2$PC, levels = df_r2$PC)
  p1 <- ggplot(df_r2, aes(x = PC, y = R2, fill = VariancePct)) +
    geom_col(color = "black", size = 0.3) +
    scale_fill_viridis_c(option = "plasma", name = "Variance explained (%)",
                         limits = c(0, max(var_explained) * 100)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Variance in PCs Explained by Dataset Factor",
      subtitle = sprintf("Top %d PCs (%.1f%% cumulative variance)", 
                         top_n, sum(var_explained[pcs_to_check]) * 100),
      x = "Principal Component",
      y = expression(R^2 ~ "(dataset factor)")
    )

  pc_df <- as.data.frame(pc$x[, 1:2])
  pc_df$Group <- dataset_factor
  pc1_var <- round(var_explained[1] * 100, 1)
  pc2_var <- round(var_explained[2] * 100, 1)
  group_colors <- palettes$two_group[1:length(unique(dataset_factor))]
  names(group_colors) <- unique(dataset_factor)

  p2 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = 2.5, alpha = 0.8, stroke = 0.5) +
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.15, 
                 linetype = "dashed", size = 0.8) +
    scale_color_manual(values = group_colors, name = title_tag) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_shape_manual(values = c(16, 17)[1:length(unique(dataset_factor))], name = title_tag) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(
      title = paste("PCA: PC1 vs PC2 –", title_tag, "Comparison"),
      x = paste0("PC1 (", pc1_var, "% variance)"),
      y = paste0("PC2 (", pc2_var, "% variance)")
    ) +
    coord_fixed(ratio = 1)

  combined_plot <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1.2))
  print(combined_plot)
  if (!is.null(outdir)) {
    fig_dir <- file.path(outdir, "figs")
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(fig_dir, paste0("PCA_s_", title_tag, ".pdf")),
           combined_plot, width = 10, height = 12, device = "pdf")
    write.csv(df_r2, file.path(fig_dir, paste0("PCA_s_", title_tag, "_R2_table.csv")), row.names = FALSE)
  }
  return(list(r2_table = df_r2, pc = pc, plot = combined_plot,
              summary = list(total_variance = sum(var_explained[pcs_to_check]),
                             significant_pcs = sum(r2_values > 0.1))))
}

# =============================================================================
#  SCREE PLOT
# =============================================================================
make_scree_plot_ <- function(pc, outdir, title_tag, k = 20) {
  palettes <- get_palettes()
  var_expl <- (pc$sdev^2) / sum(pc$sdev^2)
  k <- min(k, length(var_expl))
  df <- data.frame(
    PC = factor(paste0("PC", seq_len(k)), levels = paste0("PC", seq_len(k))),
    Var = var_expl[seq_len(k)] * 100,
    Cum = cumsum(var_expl[seq_len(k)]) * 100
  )
  p <- ggplot(df, aes(PC, Var)) +
    geom_col(fill = palettes$okabe_ito[5], color = "black", size = 0.4, alpha = 0.85) +
    geom_line(aes(y = Cum, group = 1), color = palettes$okabe_ito[6], size = 1.5) +
    geom_point(aes(y = Cum), color = palettes$okabe_ito[6], fill = "white", shape = 21, size = 3, stroke = 1.5) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, face = "italic"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Scree Plot –", title_tag),
      x = "Principal Component",
      y = "Variance Explained (%)",
      caption = "Blue bars: per-PC variance; Orange line/points: cumulative variance"
    )
  if (!is.null(outdir)) {
    dir.create(file.path(outdir, "figs"), showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(outdir, "figs", paste0("scree_", title_tag, ".pdf")), p, width = 10, height = 6)
  }
  return(p)
}

# =============================================================================
# PC ← covariate R²
# =============================================================================
pc_covariate_r2_matrix_ <- function(pc, sample_df, covariates, n_pcs = 10, 
                                    outdir = NULL, tag = "before", show_plot = TRUE) {
  stopifnot(is.list(pc), !is.null(pc$x)); stopifnot("sample" %in% colnames(sample_df))
  rn <- rownames(pc$x); if (is.null(rn)) stop("pc$x must have rownames (sample IDs).")
  md <- sample_df[match(rn, sample_df$sample), , drop = FALSE]
  k <- min(n_pcs, ncol(pc$x)); pcs <- paste0("PC", seq_len(k))
  Y <- pc$x[, seq_len(k), drop = FALSE]
  covariates <- covariates[covariates %in% colnames(md)]
  if (!length(covariates)) stop("No valid covariates found in metadata.")
  r2_mat <- matrix(NA_real_, nrow = k, ncol = length(covariates), dimnames = list(pcs, covariates))
  for (j in seq_along(covariates)) {
    v <- md[[covariates[j]]]; if (is.character(v) || is.logical(v)) v <- factor(v)
    for (i in seq_len(k)) {
      y <- Y[, i]; ok <- is.finite(y) & !is.na(v)
      if (!any(ok)) next
      if (is.factor(v) && nlevels(droplevels(v[ok])) < 2) next
      fit <- try(lm(y[ok] ~ v[ok]), silent = TRUE)
      if (!inherits(fit, "try-error")) r2_mat[i, j] <- summary(fit)$r.squared
    }
  }
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    write.csv(as.data.frame(r2_mat), file.path(outdir, paste0("PC_covariate_R2_matrix_", tag, ".csv")))
  }
  r2_df <- as.data.frame(r2_mat) %>%
    tibble::rownames_to_column("PC") %>%
    tidyr::pivot_longer(-PC, names_to = "Covariate", values_to = "R2")
  p <- ggplot(r2_df, aes(Covariate, PC, fill = R2)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey95", limits = c(0, 1), name = "R²") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 12),
          panel.grid = element_blank(),
          legend.position = "right") +
    labs(title = paste0("PC ← covariate R² (", tag, ")"), x = "Covariate", y = "Principal Component")
  if (isTRUE(show_plot)) print(p)
  if (!is.null(outdir)) ggsave(file.path(outdir, paste0("PC_covariate_R2_matrix_", tag, ".pdf")), p, width = 9, height = 7)
  invisible(list(r2 = r2_mat, plot = p))
}
```

```{r output_helpers}
# =============================================================================
# OUTPUT HELPERS
# =============================================================================
save_results <- function(spearman_scores, summary_tbl, dsmz_map, outdir){
  cat("[INFO] Saving correlation tables...\n")
  all_wide <- as.data.frame(spearman_scores) %>% tibble::rownames_to_column("sample_id")
  all_long <- all_wide %>% pivot_longer(-sample_id, names_to = "tcga_cohort", values_to = "rho")
  write.csv(all_wide, file.path(outdir, "spearman_all_wide.csv"), row.names = FALSE)
  write.csv(all_long, file.path(outdir, "spearman_all_long.csv"), row.names = FALSE)
  saveRDS(spearman_scores, file.path(outdir, "spearman_all.rds"))
  write.csv(summary_tbl, file.path(outdir, "spearman_best_by_sample.csv"), row.names = FALSE)
  write.csv(dsmz_map %>% select(sample_id, organ, tcga_code),
            file.path(outdir, "dsmz_tcga_mapping.csv"), row.names = FALSE)
}

save_debug_info <- function(tcga_counts, dsmz_raw, common_genes, outdir){
  writeLines(utils::capture.output(head(rownames(tcga_counts))), file.path(outdir, "head_tcga_genes.txt"))
  writeLines(utils::capture.output(head(dsmz_raw$Ensembl_ID)),     file.path(outdir, "head_dsmz_genes.txt"))
  writeLines(common_genes[1:min(100, length(common_genes))],       file.path(outdir, "common_genes_head.txt"))
  sink(file.path(outdir, "sessionInfo.txt")); print(sessionInfo()); sink()
}
```

```{r main_function}
# =============================================================================
# MAIN FUNCTION
# =============================================================================
main <- function(){
  cat("[INFO] Starting pipeline...\n")

  # -------------------------- Load input datasets ----------------------------
  tcga <- load_tcga_data(config$tcga_se_rds)
  dsmz <- load_dsmz_data(config$dsmz_rds, config$dsmz_meta_csv)

  # -------------------------- Preprocess gene expression ---------------------
  tcga_counts <- harmonize_gene_ids(tcga$counts); storage.mode(tcga_counts) <- "double"
  dsmz_counts <- build_dsmz_matrix(dsmz$raw);     storage.mode(dsmz_counts) <- "double"

  al <- align_metadata(dsmz$meta, dsmz_counts, config$outdir)
  dsmz_meta <- al$meta; stopifnot("organ" %in% colnames(dsmz_meta))
  dsmz_counts <- al$counts

  # -------------------------- Find common genes ------------------------------
  common <- intersect(rownames(tcga_counts), rownames(dsmz_counts))
  cat(sprintf("[INFO] Shared genes: %d\n", length(common)))
  if (length(common) < config$min_shared) cat("[WARN] Few shared genes; check Ensembl IDs\n")
  tcga_counts <- tcga_counts[common, , drop=FALSE]
  dsmz_counts <- dsmz_counts[common, , drop=FALSE]

  # -------------------------- Purity adjustment ------------------------------
  purity   <- load_purity_data(tcga$se, config$purity_tsv)
  tcga_adj <- apply_purity_adjustment(tcga_counts, purity)

  # -------------------------- Compute correlations ---------------------------
  tcga_res   <- compute_tcga_means(tcga$se, tcga_adj)
  tcga_means <- tcga_res$means
  stopifnot(identical(rownames(tcga_means), rownames(dsmz_counts)))

  sel_genes_corr <- select_variable_genes(tcga_means, config$var_gene_set)
  spearman <- compute_correlations(dsmz_counts, tcga_means, dsmz_meta, sel_genes_corr)
  dsmz_map <- perform_organ_mapping(dsmz_meta, spearman, tcga_means)

  # -------------------------- Summary statistics -----------------------------
  best_idx <- max.col(spearman, ties.method="first")
  best_overall <- tibble(
    sample_id = rownames(spearman),
    best_tcga_overall = colnames(spearman)[best_idx],
    best_rho_overall  = spearman[cbind(seq_len(nrow(spearman)), best_idx)]
  )
  rho_allowed <- rep(NA_real_, nrow(dsmz_map))
  ok <- !is.na(dsmz_map$tcga_code)
  if (any(ok)) {
    for (i in which(ok)) {
      sample_id <- dsmz_map$sample_id[i]; tcga_code <- dsmz_map$tcga_code[i]
      rho_allowed[i] <- spearman[sample_id, tcga_code]
    }
  }
  summary_tbl <- dsmz_map %>% select(sample_id, organ, tcga_code) %>%
    left_join(best_overall, by="sample_id") %>% mutate(rho_allowed = rho_allowed)

  # -------------------------- Plots + save results ---------------------------
  create_plots(summary_tbl, spearman, dsmz_map, config$outdir)
  save_results(spearman, summary_tbl, dsmz_map, config$outdir)
  save_debug_info(tcga_counts, dsmz$raw, common, config$outdir)

  # ========================= PCA ANALYSIS ONLY (UMAP REMOVED) =================
  tcga_sample_counts <- tcga_counts
  dsmz_sample_counts <- dsmz_counts

  # TCGA metadata
  tcga_lab <- sub("^TCGA-","", as.character(tcga_res$labels))
  tcga_lab <- tcga_lab[colnames(tcga_sample_counts)]
  tcga_df  <- tibble(sample = colnames(tcga_sample_counts),
                     dataset = "TCGA",
                     tcga_project = tcga_lab)
  proj_map <- tcga_project_to_organ_map()
  tcga_df <- tcga_df %>% left_join(proj_map, by = "tcga_project") %>% mutate(organ = coalesce(organ, "Other"))

  # DSMZ metadata
  dsmz_df <- tibble(sample = colnames(dsmz_sample_counts),
                    dataset = "DSMZ",
                    tcga_project = NA_character_) %>%
             left_join(dsmz_meta %>% select(sample_id, organ), by = c("sample"="sample_id")) %>%
             mutate(organ = coalesce(organ, "Other"))

  # Combine metadata
  sample_df <- bind_rows(tcga_df, dsmz_df)

  # Feature selection for PCA
  M0_log <- logCPM(cbind(tcga_sample_counts, dsmz_sample_counts))
  sel_genes_pca <- select_variable_genes(M0_log, config$var_gene_set)
  M0_log <- M0_log[sel_genes_pca, , drop=FALSE]

  # -------------------------- BEFORE batch adjustment ------------------------
  cat("[INFO] PCA BEFORE batch adjustment...\n")
  r2_before <- pc_dataset_s_(M0_log, factor(sample_df$dataset), top_n = 10, outdir = config$outdir, title_tag = "before")
  scree_before <- make_scree_plot_(r2_before$pc, config$outdir, "before")
  r2_cov_before <- pc_covariate_r2_matrix_(
    pc         = r2_before$pc,
    sample_df  = sample_df,
    covariates = c("dataset", "organ", "tcga_project"),
    n_pcs      = 10,
    outdir     = config$outdir,
    tag        = "before",
    show_plot  = TRUE
  )

  # -------------------------- Batch adjustment (SIMPLIFIED) ------------------
  cat("[INFO] Batch adjustment...\n")
  method <- tolower(config$batch_adjust_method)
  batch  <- factor(sample_df$dataset)   # e.g., TCGA vs DSMZ

  if (method == "none") {
    cat("[INFO] Skipping batch correction (config$batch_adjust_method = 'none').\n")
    M1_log <- M0_log

  } else if (method == "combat_seq") {
    # Minimal safeguard: ComBat_seq requires >=2 samples per batch
    if (min(table(batch)) < 2) {
      cat("[WARN] Some batches have <2 samples; falling back to ComBat on logCPM.\n")
      M1_log_all <- sva::ComBat(dat = as.matrix(M0_log), batch = batch,
                                mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                                mean.only = FALSE)
      M1_log <- M1_log_all
    } else {
      Xc <- cbind(tcga_sample_counts, dsmz_sample_counts)
      Xc[is.na(Xc)] <- 0
      Xc <- round(Xc)
      set.seed(42)
      Xc_adj <- sva::ComBat_seq(as.matrix(Xc), batch = batch)  # no group -> no confounding
      M1_log <- logCPM(Xc_adj)
    }

  } else if (method == "combat") {
    M1_log_all <- sva::ComBat(dat = as.matrix(M0_log), batch = batch,
                              mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                              mean.only = FALSE)
    M1_log <- M1_log_all

  } else {
    stop("[ERROR] Unknown batch_adjust_method. Use 'none', 'combat_seq', or 'combat'.")
  }

  # Keep selected genes for downstream plots
  M1_log <- M1_log[sel_genes_pca, , drop = FALSE]

  # -------------------------- AFTER batch adjustment -------------------------
  cat("[INFO] PCA AFTER batch adjustment...\n")
  r2_after <- pc_dataset_s_(M1_log, factor(sample_df$dataset), top_n = 10, outdir = config$outdir, title_tag = "after")
  scree_after <- make_scree_plot_(r2_after$pc, config$outdir, "after")
  r2_cov_after <- pc_covariate_r2_matrix_(
    pc         = r2_after$pc,
    sample_df  = sample_df,
    covariates = c("dataset", "organ", "tcga_project"),
    n_pcs      = 10,
    outdir     = config$outdir,
    tag        = "after",
    show_plot  = TRUE
  )

  # -------------------------- Summarize R² comparison ------------------------
  cat("\n**PC variance explained by dataset (R^2):**\n\n")
  r2_comparison <- tibble(
    Stage = c("Before Batch Correction", "After Batch Correction"),
    PC1_R2 = c(r2_before$r2_table$R2[r2_before$r2_table$PC == "PC1"],
               r2_after$r2_table$R2[r2_after$r2_table$PC == "PC1"]),
    PC2_R2 = c(r2_before$r2_table$R2[r2_before$r2_table$PC == "PC2"],
               r2_after$r2_table$R2[r2_after$r2_table$PC == "PC2"])
  )
  print(r2_comparison)

  cat(sprintf("[INFO] Pipeline completed successfully. Results saved to: %s\n", config$outdir))
  list(spearman_scores = spearman,
       summary_table   = summary_tbl,
       dsmz_mapping    = dsmz_map,
       r2_comparison   = r2_comparison)
}

# =============================================================================
# EXECUTION
# =============================================================================
if (!interactive()) {
  results <- main()
} else {
  cat("[INFO] Running in interactive mode. Call main() to execute pipeline.\n")
}
```

````

If you ever want to completely skip correction for a quick run, just set:
```r
config$batch_adjust_method <- "none"
````

Ping me with the next error (if any) and I’ll squash it.
