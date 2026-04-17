#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(tibble)
})

opt_list <- list(
  make_option("--features_txt", type="character", help="pan_cancer_features_clean.txt (one gene per line)"),
  make_option("--pan_expr_rds", type="character", help="pan_cancer_expr.rds (VST/log expr matrix or list containing it)"),
  make_option("--meta_tsv", type="character",
  help="joint_metadata.tsv with sample_id + (lineage or cancer_type) + (type or sample_type). role optional. cohort/component optional."),
  make_option("--brca_tables", type="character", help="BRCA deseq2_markers/tables dir"),
  make_option("--nbl_tables",  type="character", help="NBL deseq2_markers/tables dir"),
  make_option("--rbl_tables",  type="character", help="RBL deseq2_markers/tables dir"),
  make_option("--brca_unique", type="character", help="BRCA unique_feature_set_recurrence_ge_2.txt"),
  make_option("--nbl_unique",  type="character", help="NBL unique_feature_set_recurrence_ge_2.txt"),
  make_option("--rbl_unique",  type="character", help="RBL unique_feature_set_recurrence_ge_2.txt"),
  make_option("--out_dir", type="character", default="deg_eval"),
  make_option("--padj_anchor", type="double", default=0.05),
  make_option("--padj_isolate", type="double", default=0.01),
  make_option("--label_top", type="integer", default=12),
  make_option("--heatmap_top", type="integer", default=40)
)
opt <- parse_args(OptionParser(option_list=opt_list))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
read_gene_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- x[!grepl("^\\s*$", x)]
  x <- x[!grepl("^\\s*#", x)]
  unique(trimws(x))
}

infer_role <- function(path) {
  b <- tolower(basename(path))
  if (startsWith(b, "anchor_")) return("anchor")
  if (startsWith(b, "isolate_")) return("isolate")
  return(NA_character_)
}

read_deseq <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  if (!("gene" %in% names(df))) {
    gene_col <- intersect(names(df), c("gene_id","Gene","geneid","ensg","ENSG"))
    if (length(gene_col) == 0) stop("No gene column found in: ", path)
    df <- df |> rename(gene = !!gene_col[1])
  }
  req <- c("gene","log2FoldChange","padj","baseMean")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing columns in ", path, ": ", paste(miss, collapse=", "))
  df |> mutate(source_file = basename(path))
}

# robust extraction of expression matrix from pan_cancer_expr.rds
extract_expr_matrix <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) return(as.matrix(obj))

  if (is.list(obj)) {
    for (k in c("expr", "expr_mat", "mat", "matrix", "vst", "vst_mat", "assay", "x")) {
      if (!is.null(obj[[k]])) {
        m <- obj[[k]]
        if (is.data.frame(m)) m <- as.matrix(m)
        if (is.matrix(m)) return(m)
      }
    }
  }
  stop("Could not extract expression matrix from RDS. Provide a matrix or a list containing one of keys: expr/expr_mat/mat/matrix/vst/vst_mat/assay/x")
}

# ----------------------------
# Inputs
# ----------------------------
features_retained <- read_gene_list(opt$features_txt)

unique_brca <- read_gene_list(opt$brca_unique)
unique_nbl  <- read_gene_list(opt$nbl_unique)
unique_rbl  <- read_gene_list(opt$rbl_unique)

meta <- read_tsv(opt$meta_tsv, show_col_types = FALSE)

# ----------------------------
# Standardize ID + disease column
# ----------------------------
if ("sample_id" %in% names(meta) && !("sample" %in% names(meta))) {
  meta <- meta |> rename(sample = sample_id)
}

# disease can come from lineage OR cancer_type
if (!("disease" %in% names(meta))) {
  if ("lineage" %in% names(meta)) {
    meta <- meta |> rename(disease = lineage)
  } else if ("cancer_type" %in% names(meta)) {
    meta <- meta |> rename(disease = cancer_type)
  } else {
    stop("meta_tsv must contain either 'lineage' or 'cancer_type' to define disease.")
  }
}

# ----------------------------
# Standardize sample_type
# ----------------------------
# pan_cancer file uses `type`; benchmark file uses `sample_type`
if (!("sample_type" %in% names(meta))) {
  if ("type" %in% names(meta)) {
    meta <- meta |> rename(sample_type = type)
  } else {
    stop("meta_tsv must contain 'type' (pan_cancer) or 'sample_type' (benchmark).")
  }
}

meta <- meta |>
  mutate(
    sample_type = tolower(trimws(as.character(sample_type))),
    sample_type = case_when(
      sample_type %in% c("cell line", "cellline", "cell_line") ~ "cell_line",
      sample_type %in% c("tumor") ~ "tumour",
      TRUE ~ sample_type
    ),
    disease = toupper(trimws(as.character(disease)))
  )

# ----------------------------
# Optional columns (keep if present)
# ----------------------------
# DO NOT infer anchor/isolate here
if (!("role" %in% names(meta))) meta <- meta |> mutate(role = "other")

# cohort exists in benchmark meta; ignore if missing
# component may exist in unsupervised meta; ignore if missing

# Final checks
stopifnot(all(c("sample", "disease", "sample_type", "role") %in% names(meta)))

expr_obj <- readRDS(opt$pan_expr_rds)
expr_mat <- extract_expr_matrix(expr_obj)

if (is.null(rownames(expr_mat))) stop("Expression matrix has no rownames (genes).")
if (is.null(colnames(expr_mat))) stop("Expression matrix has no colnames (samples).")

common_samples <- intersect(colnames(expr_mat), meta$sample)
if (length(common_samples) < 5) stop("Too few overlapping samples between expr matrix and metadata.")
expr_mat <- expr_mat[, common_samples, drop = FALSE]
meta <- meta |> filter(sample %in% common_samples) |> arrange(match(sample, colnames(expr_mat)))

message("[INFO] sample_type counts in expr/meta overlap:")
print(table(meta$sample_type))
message("[INFO] disease counts in expr/meta overlap:")
print(table(meta$disease))

# ----------------------------
# Load DE tables (background for volcano/MA)
# ----------------------------
load_tables_from_dir <- function(dir, disease) {
  files <- list.files(dir, pattern="\\.tsv(\\.gz)?$", full.names=TRUE)
  if (length(files) == 0) return(tibble())
  info <- tibble(path = files) |>
    mutate(role = vapply(path, infer_role, character(1))) |>
    filter(!is.na(role)) |>
    mutate(disease = disease)

  if (nrow(info) == 0) return(tibble())

  out <- info |>
    rowwise() |>
    do({
      m <- .
      read_deseq(m$path) |>
        mutate(disease = m$disease, role = m$role)
    }) |>
    ungroup()

  out
}

de_brca <- load_tables_from_dir(opt$brca_tables, "BRCA")
de_nbl  <- load_tables_from_dir(opt$nbl_tables,  "NBL")
de_rbl  <- load_tables_from_dir(opt$rbl_tables,  "RBL")
all_de <- bind_rows(de_brca, de_nbl, de_rbl)

if (nrow(all_de) == 0) stop("No DE tables loaded. Check *tables dirs and filename prefixes (anchor_/isolate_).")

all_de <- all_de |>
  filter(!is.na(log2FoldChange), !is.na(padj), padj > 0, !is.na(baseMean)) |>
  mutate(in_retained = gene %in% features_retained,
         class = case_when(
           !in_retained ~ "other",
           log2FoldChange >= 0 ~ "retained_up",
           TRUE ~ "retained_down"
         ),
         neglog10_padj = -log10(padj),
         padj_thr = if_else(role == "isolate", opt$padj_isolate, opt$padj_anchor))

# ----------------------------
# Volcano plots
# ----------------------------
plot_volcano_disease <- function(disease) {
  df <- all_de |> filter(disease == !!disease)

  thr_df <- df |> distinct(role, padj_thr)

  lab <- df |>
    filter(in_retained) |>
    group_by(role) |>
    mutate(abs_lfc = abs(log2FoldChange)) |>
    arrange(desc(abs_lfc)) |>
    slice_head(n = opt$label_top) |>
    ungroup()

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = class, shape = role), alpha = 0.55, size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(
      data = thr_df,
      aes(yintercept = -log10(padj_thr)),
      linetype = "dashed"
    ) +
    facet_wrap(~role, nrow = 1) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3
    ) +
    scale_color_manual(values = c(
      other = "grey75",
      retained_up = "red3",
      retained_down = "blue3"
    )) +
    labs(
      title = paste0("Volcano — ", disease, " (highlight: pan-cancer retained genes)"),
      x = "log2 fold change",
      y = expression(-log[10]("padj"))
    ) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(file.path(opt$out_dir, paste0("Volcano_", disease, "_facet_role.png")),
         p, width = 10, height = 4, dpi = 300)
}

for (d in c("BRCA","NBL","RBL")) plot_volcano_disease(d)

thr_df_comb <- all_de |> distinct(disease, role, padj_thr)

p_comb <- ggplot(all_de, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(aes(color = class), alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(
    data = thr_df_comb,
    aes(yintercept = -log10(padj_thr)),
    linetype = "dashed"
  ) +
  facet_grid(disease ~ role) +
  scale_color_manual(values = c(
    other = "grey80",
    retained_up = "red3",
    retained_down = "blue3"
  )) +
  labs(
    title = "Combined Volcano — by disease × role (highlight: pan-cancer retained genes)",
    x = "log2 fold change",
    y = expression(-log[10]("padj"))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(file.path(opt$out_dir, "Volcano_combined_disease_by_role.png"),
       p_comb, width = 10, height = 7, dpi = 300)

# ----------------------------
# Per-contrast volcano (one contrast = one plot; interpretable)
# ----------------------------
plot_volcano_contrast <- function(disease, role, source_file) {
  df <- all_de |> filter(disease == !!disease, role == !!role, source_file == !!source_file)
  if (nrow(df) == 0) return(invisible(NULL))

  padj_thr <- unique(df$padj_thr)[1]

  lab <- df |> filter(in_retained) |>
    mutate(abs_lfc = abs(log2FoldChange)) |>
    arrange(desc(abs_lfc)) |>
    slice_head(n = opt$label_top)

  p <- ggplot(df, aes(log2FoldChange, neglog10_padj)) +
    geom_point(aes(color = class), alpha = 0.55, size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = Inf) +
    scale_color_manual(values = c(other = "grey75", retained_up = "red3", retained_down = "blue3")) +
    labs(
      title = paste0("Volcano — ", disease, " ", role, " (", source_file, ")"),
      x = "log2 fold change",
      y = expression(-log[10]("padj"))
    ) +
    theme_classic() +
    theme(legend.position = "none")

  base_name <- sub("\\.tsv(\\.gz)?$", "", source_file)
  out <- file.path(opt$out_dir, paste0("Volcano_", disease, "_", role, "_", base_name, ".png"))
  ggsave(out, p, width = 7, height = 5, dpi = 300)
}

uniq_contrasts <- all_de |> distinct(disease, role, source_file)
for (i in seq_len(nrow(uniq_contrasts))) {
  row <- uniq_contrasts[i, ]
  plot_volcano_contrast(row$disease, row$role, row$source_file)
}

# ----------------------------
# MA plots
# ----------------------------
plot_ma_grid <- function() {
  df <- all_de |>
    mutate(ma_sig = in_retained & (padj <= padj_thr))

  p <- ggplot(df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = ma_sig), alpha = 0.55, size = 0.9) +
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(disease ~ role) +
    scale_color_manual(values = c(`FALSE`="grey80", `TRUE`="red3")) +
    labs(
      title = "MA plots — by disease × role (red = in retained set AND padj below role threshold)",
      x = "baseMean (log scale)",
      y = "log2 fold change"
    ) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(file.path(opt$out_dir, "MA_combined_disease_by_role.png"),
         p, width = 10, height = 7, dpi = 300)
}
plot_ma_grid()

for (d in c("BRCA","NBL","RBL")) {
  for (r in c("anchor","isolate")) {
    df <- all_de |> filter(disease == d, role == r) |>
      mutate(ma_sig = in_retained & (padj <= padj_thr))
    if (nrow(df) == 0) next

    p <- ggplot(df, aes(x = baseMean, y = log2FoldChange)) +
      geom_point(aes(color = ma_sig), alpha = 0.55, size = 1.0) +
      scale_x_log10() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c(`FALSE`="grey80", `TRUE`="red3")) +
      labs(
        title = paste0("MA — ", d, " (", r, ")"),
        x = "baseMean (log scale)",
        y = "log2 fold change"
      ) +
      theme_classic() +
      theme(legend.position = "none")

    ggsave(file.path(opt$out_dir, paste0("MA_", d, "_", r, ".png")),
           p, width = 7, height = 5, dpi = 300)
  }
}

# ----------------------------
# UpSet (recurrence≥2 unique sets)
# ----------------------------
sets <- tibble(
  gene = c(unique_brca, unique_nbl, unique_rbl),
  set  = c(rep("BRCA", length(unique_brca)),
           rep("NBL",  length(unique_nbl)),
           rep("RBL",  length(unique_rbl)))
) |> distinct()

wide <- sets |>
  mutate(value = 1) |>
  pivot_wider(names_from = set, values_from = value, values_fill = 0)

wide$pattern <- apply(wide[, c("BRCA","NBL","RBL")], 1, paste0, collapse = "")
intersections <- wide |>
  count(pattern, sort = TRUE) |>
  mutate(
    BRCA = substr(pattern,1,1),
    NBL  = substr(pattern,2,2),
    RBL  = substr(pattern,3,3)
  )

write_tsv(intersections, file.path(opt$out_dir, "UpSet_unique_recurrence_ge2_intersections.tsv"))

p_up <- intersections |>
  ggplot(aes(x = reorder(pattern, -n), y = n)) +
  geom_col() +
  labs(
    title = "Overlap of recurrence≥2 unique marker sets (BRCA/NBL/RBL)",
    x = "Intersection pattern (1=in set)",
    y = "Number of genes"
  ) +
  theme_classic()

ggsave(file.path(opt$out_dir, "UpSet_unique_recurrence_ge2_barplot.png"),
       p_up, width = 8, height = 4, dpi = 300)

# ----------------------------
# Heatmaps (per disease): genes from disease-specific unique recurrence≥2 sets
# - Samples: disease only
# - Genes: top N from unique_feature_set_recurrence_ge_2.txt
#   ranked by max |log2FC| in that disease (fallback: input order)
# - Values: pan_cancer_expr.rds matrix; z-score per gene
# - Drop genes with zero variance across selected samples
# ----------------------------
ann_col <- meta |>
  select(sample, disease, sample_type, any_of(c("cohort", "component"))) |>
  column_to_rownames("sample")

unique_sets <- list(
  BRCA = unique_brca,
  NBL  = unique_nbl,
  RBL  = unique_rbl
)

for (d in c("BRCA","NBL","RBL")) {

  cand <- intersect(unique_sets[[d]], rownames(expr_mat))
  if (length(cand) < 5) next

  rank_tbl <- all_de |>
    filter(disease == d, gene %in% cand) |>
    group_by(gene) |>
    summarise(max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE), .groups="drop") |>
    arrange(desc(max_abs_lfc))

  ranked_genes <- rank_tbl |> pull(gene)
  if (length(ranked_genes) < 5) ranked_genes <- cand

  top_genes <- head(ranked_genes, opt$heatmap_top)
  if (length(top_genes) < 5) next

  samp <- meta |>
    filter(disease == d, sample_type == "cell_line") |>
    pull(sample)
  samp <- intersect(samp, colnames(expr_mat))
  if (length(samp) < 5) next

  sub <- expr_mat[top_genes, samp, drop = FALSE]

  keep <- apply(sub, 1, sd, na.rm = TRUE) > 0
  sub <- sub[keep, , drop = FALSE]
  if (nrow(sub) < 5) next

  sub_z <- t(scale(t(sub)))
  sub_z[is.na(sub_z)] <- 0

  base_name <- paste0("Heatmap_top", opt$heatmap_top, "_", d, "_uniqueRecGE2")
  out_png <- file.path(opt$out_dir, paste0(base_name, ".png"))
  out_pdf <- file.path(opt$out_dir, paste0(base_name, ".pdf"))

  png(out_png, width = 1600, height = 1200, res = 200)
  pheatmap(
    sub_z,
    annotation_col = ann_col[samp, , drop = FALSE],
    show_colnames = FALSE,
    fontsize_row = 8,
    main = paste0("Top ", opt$heatmap_top, " genes from unique recurrence≥2 set — ", d, " (z-scored)")
  )
  dev.off()

  pdf(out_pdf, width = 8, height = 6)
  pheatmap(
    sub_z,
    annotation_col = ann_col[samp, , drop = FALSE],
    show_colnames = FALSE,
    fontsize_row = 8,
    main = paste0("Top ", opt$heatmap_top, " genes from unique recurrence≥2 set — ", d, " (z-scored)")
  )
  dev.off()
}

writeLines(as.character(Sys.time()), file.path(opt$out_dir, "deg_eval.done"))
message("[DONE] Wrote DEG-eval outputs to: ", opt$out_dir)
