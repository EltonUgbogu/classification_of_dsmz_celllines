#!/usr/bin/env Rscript

# characterize_clusters.R
# Post-hoc biological characterisation of consensus clusters
# - supports consensus *_clusters_optimal.rds (list with $clusters named vector)
# - supports legacy CSV: sample, cluster
# - optional post-hoc label enrichment/validation from metadata TSV

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggrepel)

  # CLI + optional external validation
  library(optparse)
  library(mclust)   # ARI
  library(aricode)  # NMI
})

# ============================================================================
# CLI
# ============================================================================

option_list <- list(
  make_option(c("--vst_rds"), type = "character", default = NULL,
              help = "RDS containing expression matrix (genes x samples). REQUIRED."),
  make_option(c("--clusters_rds"), type = "character", default = NULL,
              help = "RDS containing cluster assignments (e.g. *_clusters_optimal.rds). Preferred."),
  make_option(c("--clusters_csv"), type = "character", default = NULL,
              help = "Legacy CSV with at least two columns: sample, cluster."),
  make_option(c("--outdir"), type = "character", default = NULL,
              help = "Output directory. REQUIRED."),

  # Optional metadata / label agreement (post-hoc)
  make_option(c("--meta_tsv"), type = "character", default = NULL,
              help = "TSV metadata file (e.g., joint_metadata.tsv). Optional."),
  make_option(c("--meta_id_col"), type = "character", default = "sample",
              help = "Column in meta_tsv corresponding to sample IDs."),
  make_option(c("--label_col"), type = "character", default = "cancer_type",
              help = "Metadata column to enrich/validate against."),
  make_option(c("--cross_threshold"), type = "double", default = 0.20,
              help = "Flag clusters as cross-label if top label fraction < this threshold."),

  # GSEA
  make_option(c("--gsea_collection"), type = "character", default = "H",
              help = "MSigDB collection to use via msigdbr. Default: H (Hallmark)."),
  make_option(c("--min_genes_gsea"), type = "integer", default = 20,
              help = "Minimum ranked genes required to run GSEA."),

  # Cell-line characterisation (tumour-informed)
  make_option(c("--cell_pconsensus_tsv"), type = "character", default = NULL,
              help = "TSV with cell_line vs tumour_cluster p-consensus scores (winner map). Optional.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vst_rds) || !nzchar(opt$vst_rds)) stop("[ERROR] --vst_rds is required.")
if (is.null(opt$outdir)  || !nzchar(opt$outdir))  stop("[ERROR] --outdir is required.")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Helpers
# ============================================================================

harmonise_ids <- function(x) {
  gsub("_", "-", x, fixed = TRUE)
}

as_named_cluster_factor <- function(x) {
  # Ensure: named factor
  if (is.null(names(x))) stop("[ERROR] cluster vector must be named (names = sample IDs).")
  x <- as.character(x)
  x <- factor(x)
  x
}

load_cluster_vec <- function(clusters_rds = NULL, clusters_csv = NULL) {
  cluster_vec <- NULL

  if (!is.null(clusters_rds) && nzchar(clusters_rds)) {
    if (!file.exists(clusters_rds)) stop("[ERROR] clusters_rds not found: ", clusters_rds)
    cat("[INFO] Loading clusters from RDS:", clusters_rds, "\n")
    obj <- readRDS(clusters_rds)

    if (is.list(obj) && "clusters" %in% names(obj)) {
      cluster_vec <- obj$clusters
    } else if (is.vector(obj)) {
      cluster_vec <- obj
    } else if (is.data.frame(obj) && ncol(obj) >= 2) {
      cluster_vec <- obj[[2]]
      names(cluster_vec) <- obj[[1]]
    } else {
      stop("[ERROR] Unsupported clusters_rds format.")
    }

    names(cluster_vec) <- harmonise_ids(names(cluster_vec))
    cluster_vec <- as_named_cluster_factor(cluster_vec)
    return(cluster_vec)
  }

  if (!is.null(clusters_csv) && nzchar(clusters_csv)) {
    if (!file.exists(clusters_csv)) stop("[ERROR] clusters_csv not found: ", clusters_csv)
    cat("[INFO] Loading clusters from CSV:", clusters_csv, "\n")
    df <- read.csv(clusters_csv, stringsAsFactors = FALSE)
    if (ncol(df) < 2) stop("[ERROR] clusters_csv must have at least two columns: sample, cluster.")
    colnames(df)[1:2] <- c("sample", "cluster")
    df$sample <- harmonise_ids(df$sample)
    cluster_vec <- df$cluster
    names(cluster_vec) <- df$sample
    cluster_vec <- as_named_cluster_factor(cluster_vec)
    return(cluster_vec)
  }

  stop("[ERROR] Provide either --clusters_rds or --clusters_csv.")
}

map_ensembl_to_symbol <- function(vst) {
  gene_ids <- rownames(vst)
  use_symbols <- TRUE

  if (!any(grepl("^ENSG", gene_ids))) {
    cat("[INFO] Row names do not look like Ensembl IDs; assuming symbols.\n")
    rownames(vst) <- toupper(rownames(vst))
    return(list(vst = vst, use_symbols = TRUE))
  }

  cat("[INFO] Detected Ensembl-like IDs; mapping to HGNC symbols.\n")
  ensg_stripped <- sub("\\.\\d+$", "", gene_ids)

  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensg_stripped,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )

  valid <- !is.na(mapped) & mapped != ""
  cat(sprintf("[INFO] Mapped %d / %d Ensembl IDs to symbols (%.1f%%).\n",
              sum(valid), length(valid), 100 * mean(valid)))

  if (sum(valid) == 0L) {
    cat("[WARN] No Ensembl IDs mapped; keeping Ensembl IDs.\n")
    return(list(vst = vst, use_symbols = FALSE))
  }

  vst2 <- vst[valid, , drop = FALSE]
  rownames(vst2) <- toupper(mapped[valid])

  # Collapse duplicates (same symbol) by mean
  if (any(duplicated(rownames(vst2)))) {
    cat("[INFO] Collapsing duplicated symbols by mean expression.\n")
    vst2 <- limma::avereps(vst2, ID = rownames(vst2))
  }

  list(vst = vst2, use_symbols = TRUE)
}

# ============================================================================
# 1. Load expression + clusters
# ============================================================================

cat("[INFO] Loading expression matrix...\n")
if (!file.exists(opt$vst_rds)) stop("[ERROR] VST RDS not found: ", opt$vst_rds)
vst <- readRDS(opt$vst_rds)  # genes x samples
cat("[INFO] VST dims:", paste(dim(vst), collapse = " x "), "\n")

colnames(vst) <- harmonise_ids(colnames(vst))

cluster_vec <- load_cluster_vec(opt$clusters_rds, opt$clusters_csv)

# Align samples
common <- intersect(colnames(vst), names(cluster_vec))
if (length(common) == 0L) stop("[ERROR] No overlap between VST samples and cluster labels.")
common <- sort(common)
vst <- vst[, common, drop = FALSE]
cluster_vec <- cluster_vec[common]
stopifnot(identical(colnames(vst), names(cluster_vec)))

cat(sprintf("[INFO] Using %d aligned samples.\n", length(common)))
print(table(cluster_vec))

# ============================================================================
# 2. Gene IDs: Ensembl -> SYMBOL (recommended for GSEA)
# ============================================================================

mapped_out <- map_ensembl_to_symbol(vst)
vst <- mapped_out$vst
genes <- rownames(vst)

# ============================================================================
# 3. Cluster centroids
# ============================================================================

cat("[INFO] Computing cluster centroids...\n")
centroids <- sapply(levels(cluster_vec), function(cl) {
  rowMeans(vst[, cluster_vec == cl, drop = FALSE])
})

saveRDS(centroids, file.path(opt$outdir, "cluster_centroids.rds"))
write.csv(
  as.data.frame(centroids) %>% tibble::rownames_to_column("gene"),
  file.path(opt$outdir, "cluster_centroids.csv"),
  row.names = FALSE
)

# ============================================================================
# 4. limma DE: each cluster vs rest
# ============================================================================

cat("[INFO] Running limma DE per cluster vs rest...\n")
marker_list <- list()

for (cl in levels(cluster_vec)) {
  cat(sprintf("  - Cluster %s\n", cl))
  group <- factor(cluster_vec == cl, levels = c(FALSE, TRUE))
  design <- model.matrix(~ group)
  fit <- lmFit(vst, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  marker_list[[cl]] <- tt
  write.csv(head(tt, 50),
            file.path(opt$outdir, sprintf("cluster_%s_top50_markers.csv", cl)),
            row.names = FALSE)
}

saveRDS(marker_list, file.path(opt$outdir, "cluster_DE_markers.rds"))

top10_list <- lapply(marker_list, function(tt) {
  tt2 <- tt[order(-tt$t), ]
  head(tt2, 10)
})
saveRDS(top10_list, file.path(opt$outdir, "cluster_top10_markers_by_t.rds"))
write.csv(
  purrr::imap_dfr(top10_list, ~ dplyr::mutate(.x, cluster = .y)),
  file.path(opt$outdir, "cluster_top10_markers_by_t.csv"),
  row.names = FALSE
)

# ============================================================================
# 5. GSEA (MSigDB via msigdbr + fgseaMultilevel)
# ============================================================================

cat("[INFO] Preparing gene sets via msigdbr...\n")
msig <- msigdbr(species = "Homo sapiens", collection = opt$gsea_collection)
pathways_list <- split(msig$gene_symbol, msig$gs_name)

gsea_results <- list()

cat("[INFO] Running GSEA per cluster...\n")
for (cl in levels(cluster_vec)) {
  df <- marker_list[[cl]]
  ranks <- df$t
  names(ranks) <- df$gene
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < opt$min_genes_gsea) {
    cat("[WARN] Cluster ", cl, ": too few ranked genes; skipping GSEA.\n")
    gsea_results[[cl]] <- tibble()
    next
  }

  res <- fgseaMultilevel(pathways = pathways_list, stats = ranks)
  res <- as.data.frame(res)
  gsea_results[[cl]] <- res

  res_clean <- res
  if (nrow(res_clean) > 0) {
    res_clean <- res_clean %>%
      mutate(across(where(is.list),
                    ~ vapply(.x, function(z) paste(z, collapse = ";"), character(1))))
  }

  write.csv(
    res_clean %>% arrange(padj, desc(NES)),
    file.path(opt$outdir, sprintf("cluster_%s_GSEA_msig_%s.csv", cl, opt$gsea_collection)),
    row.names = FALSE
  )
}

saveRDS(gsea_results, file.path(opt$outdir, sprintf("cluster_GSEA_msig_%s_results.rds", opt$gsea_collection)))

# ============================================================================
# 6. Marker heatmaps (top 30 upregulated per cluster)
# ============================================================================

cat("[INFO] Building per-cluster marker heatmaps...\n")
for (cl in levels(cluster_vec)) {
  tt <- marker_list[[cl]]
  tt <- tt[order(-tt$t), ]
  genes_cl <- head(tt$gene, 30)
  genes_cl <- intersect(genes_cl, rownames(vst))

  if (length(genes_cl) < 2) next

  samples_cl <- names(cluster_vec)[cluster_vec == cl]
  ann_col <- data.frame(cluster = cluster_vec[samples_cl], stringsAsFactors = FALSE)
  rownames(ann_col) <- samples_cl
  mat_cl <- vst[genes_cl, samples_cl, drop = FALSE]

  pdf(file.path(opt$outdir, sprintf("Cluster_%s_marker_heatmap_top30.pdf", cl)),
      width = 7, height = 8)
  pheatmap(mat_cl,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = ann_col,
           show_rownames = TRUE,
           fontsize_row = 6,
           fontsize_col = 6,
           main = sprintf("Cluster %s – top 30 up-regulated markers", cl))
  dev.off()
}

# ============================================================================
# 7. Volcano plots
# ============================================================================

cat("[INFO] Generating volcano plots...\n")
for (cl in levels(cluster_vec)) {
  tt <- marker_list[[cl]] %>%
    mutate(
      log10padj = -log10(adj.P.Val + 1e-300),
      signif = adj.P.Val < 0.05 & abs(logFC) > 1
    ) %>%
    arrange(adj.P.Val) %>%
    mutate(label = ifelse(signif & row_number() <= 10, gene, ""))

  p <- ggplot(tt, aes(x = logFC, y = log10padj)) +
    geom_point(aes(color = signif), alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 3, max.overlaps = 10,
                             min.segment.length = 0) +
    labs(
      title = paste("Cluster", cl, "vs rest"),
      x = "log2 fold-change",
      y = "-log10(adj. P-value)",
      color = "Significant"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path(opt$outdir, sprintf("cluster_%s_volcano.pdf", cl)),
         p, width = 6, height = 5)
}

# ============================================================================
# 8. NES heatmap across clusters (filtered)
# ============================================================================

cat("[INFO] Building NES heatmap across clusters...\n")
nes_long <- purrr::map2_dfr(gsea_results, names(gsea_results), function(df, cl) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% transmute(pathway, cluster = cl, NES, padj)
})

if (!is.null(nes_long) && nrow(nes_long) > 0) {
  sel_paths <- nes_long %>% filter(NES > 0, padj < 0.05) %>% pull(pathway) %>% unique()

  if (length(sel_paths) >= 2) {
    nes_mat <- nes_long %>%
      filter(pathway %in% sel_paths) %>%
      pivot_wider(id_cols = pathway, names_from = cluster, values_from = NES, values_fill = 0) %>%
      column_to_rownames("pathway") %>%
      as.matrix()

    pdf(file.path(opt$outdir, sprintf("MSig_%s_NES_heatmap_filtered.pdf", opt$gsea_collection)),
        width = 7, height = 6)
    pheatmap(nes_mat,
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, fontsize_row = 6,
             main = sprintf("MSig %s NES (NES>0, padj<0.05)", opt$gsea_collection))
    dev.off()
  } else {
    cat("[INFO] NES heatmap skipped: <2 significant pathways.\n")
  }
} else {
  cat("[INFO] No GSEA results; skipping NES heatmap.\n")
}

# ============================================================================
# 9. Post-hoc label enrichment / agreement (optional)
# ============================================================================

if (!is.null(opt$meta_tsv) && nzchar(opt$meta_tsv) && file.exists(opt$meta_tsv)) {

  cat("[INFO] Loading metadata TSV for post-hoc evaluation...\n")
  meta <- readr::read_tsv(opt$meta_tsv, show_col_types = FALSE)

  if (!(opt$meta_id_col %in% colnames(meta))) stop("[ERROR] meta_id_col not in metadata: ", opt$meta_id_col)
  if (!(opt$label_col %in% colnames(meta)))   stop("[ERROR] label_col not in metadata: ", opt$label_col)

  meta[[opt$meta_id_col]] <- harmonise_ids(as.character(meta[[opt$meta_id_col]]))

  meta_sub <- meta %>%
    select(!!opt$meta_id_col, !!opt$label_col) %>%
    rename(sample = !!opt$meta_id_col, label = !!opt$label_col) %>%
    filter(sample %in% names(cluster_vec)) %>%
    distinct(sample, .keep_all = TRUE)

  meta_sub <- meta_sub[match(names(cluster_vec), meta_sub$sample), , drop = FALSE]
  label_vec <- meta_sub$label
  names(label_vec) <- meta_sub$sample

  df_eval <- tibble(sample = names(cluster_vec),
                    cluster = as.character(cluster_vec),
                    label = as.character(label_vec)) %>%
    filter(!is.na(label), label != "")

  if (nrow(df_eval) > 0) {

    # Purity
    purity_tbl <- df_eval %>%
      count(cluster, label, name = "n") %>%
      group_by(cluster) %>%
      mutate(frac = n / sum(n)) %>%
      arrange(cluster, desc(frac)) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      rename(top_label = label, top_label_n = n, top_label_frac = frac) %>%
      mutate(flag_cross_label = top_label_frac < opt$cross_threshold)

    write.csv(purity_tbl,
              file.path(opt$outdir, sprintf("cluster_purity_%s.csv", opt$label_col)),
              row.names = FALSE)

    # Agreement
    ari <- mclust::adjustedRandIndex(df_eval$cluster, df_eval$label)
    nmi <- aricode::NMI(df_eval$cluster, df_eval$label)

    eval_tbl <- tibble(
      label_col = opt$label_col,
      n_samples = nrow(df_eval),
      ARI = ari,
      NMI = nmi,
      cross_threshold = opt$cross_threshold
    )
    write.csv(eval_tbl,
              file.path(opt$outdir, sprintf("cluster_label_agreement_%s.csv", opt$label_col)),
              row.names = FALSE)

    # Fisher enrichment
    enrich <- expand.grid(
      cluster = sort(unique(df_eval$cluster)),
      label = sort(unique(df_eval$label)),
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(
        a = sum(df_eval$cluster == cluster & df_eval$label == label),
        b = sum(df_eval$cluster == cluster & df_eval$label != label),
        c = sum(df_eval$cluster != cluster & df_eval$label == label),
        d = sum(df_eval$cluster != cluster & df_eval$label != label),
        p = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value,
        odds_ratio = (a * d + 1e-9) / (b * c + 1e-9)
      ) %>%
      ungroup() %>%
      group_by(cluster) %>%
      mutate(p_adj = p.adjust(p, method = "BH")) %>%
      ungroup() %>%
      mutate(log10_fdr = -log10(p_adj + 1e-300))

    write.csv(enrich,
              file.path(opt$outdir, sprintf("cluster_enrichment_fisher_%s.csv", opt$label_col)),
              row.names = FALSE)

  } else {
    cat("[INFO] No non-missing labels to evaluate; skipping enrichment.\n")
  }
} else {
  cat("[INFO] No metadata TSV provided; skipping label enrichment/validation.\n")
}

# ============================================================================
# 10. Summary table (markers + top pathways)
# ============================================================================

cat("[INFO] Writing cluster summary table...\n")

key_markers_vec <- purrr::map_chr(levels(cluster_vec), function(cl) {
  tt <- marker_list[[cl]]
  up <- tt[order(-tt$t), ]
  paste(head(up$gene, 10), collapse = ", ")
})

top_pathways_vec <- purrr::map_chr(levels(cluster_vec), function(cl) {
  res <- gsea_results[[cl]]
  if (is.null(res) || nrow(res) == 0) return("None")
  res_filt <- res %>%
    as.data.frame() %>%
    filter(NES > 0, padj < 0.05) %>%
    arrange(padj, desc(NES))
  if (nrow(res_filt) == 0) return("None")
  head_names <- head(res_filt$pathway, 5) %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
  paste(head_names, collapse = "; ")
})

n_per_cluster <- as.integer(table(cluster_vec)[levels(cluster_vec)])

summary_tbl <- tibble(
  Cluster = levels(cluster_vec),
  N_samples = n_per_cluster,
  Key_markers = key_markers_vec,
  Top_pathways = top_pathways_vec,
  Notes = ""
)

write.csv(summary_tbl, file.path(opt$outdir, "cluster_summary.csv"), row.names = FALSE)

cat("\n================ CLUSTER CHARACTERISATION SUMMARY ================\n")
print(summary_tbl %>% select(Cluster, N_samples, Top_pathways))
cat("==================================================================\n")

# ============================================================================
# 11. Tumour-informed cell-line characterisation (derived)
# ============================================================================

if (!is.null(opt$cell_pconsensus_tsv) &&
    nzchar(opt$cell_pconsensus_tsv) &&
    file.exists(opt$cell_pconsensus_tsv)) {

  cat("[INFO] Building tumour-informed cell-line characterisation...\n")

  # Expected columns:
  # cell_line | cluster | p_consensus
  pc <- readr::read_tsv(opt$cell_pconsensus_tsv, show_col_types = FALSE)

  required_cols <- c("cell_line", "cluster", "p_consensus")
  if (!all(required_cols %in% colnames(pc))) {
    stop("[ERROR] cell_pconsensus_tsv must contain: cell_line, cluster, p_consensus")
  }

  # Harmonise IDs
  pc$cell_line <- harmonise_ids(pc$cell_line)

  # Dominant cluster per cell line
  cell_dom <- pc %>%
    group_by(cell_line) %>%
    arrange(desc(p_consensus)) %>%
    summarise(
      dominant_cluster = first(cluster),
      p_consensus = first(p_consensus),
      specificity = first(p_consensus) / sum(p_consensus + 1e-9),
      .groups = "drop"
    )

  # Attach tumour program annotations
  cluster_prog <- read.csv(file.path(opt$outdir, "cluster_summary.csv"),
                            stringsAsFactors = FALSE) %>%
    select(Cluster, Top_pathways) %>%
    rename(dominant_cluster = Cluster,
           dominant_program = Top_pathways)

  cell_char <- cell_dom %>%
    left_join(cluster_prog, by = "dominant_cluster") %>%
    arrange(desc(p_consensus))

  write.table(
    cell_char,
    file.path(opt$outdir, "cell_line_characterisation.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  cat("[OUTPUT] Tumour-informed cell-line characterisation written to:\n",
      file.path(opt$outdir, "cell_line_characterisation.tsv"), "\n")

} else {
  cat("[INFO] No cell_pconsensus_tsv provided; skipping cell-line characterisation.\n")
}

cat("[SUCCESS] Cluster characterization finished.\n")

