#!/usr/bin/env Rscript

# =============================================================================
# cellline_centred_precision_at_k.R
# =============================================================================
#
# Cell line centred tumour class similarity analysis using Precision@k.
#
# Scientific question
#   For each cell line, which clinical patient tumour samples does it most
#   closely resemble, and how specifically and broadly does it rank same
#   lineage patient samples at the top of the similarity order?
#
# Computation
#   Arithmetic-mean aggregation across replicate profiles, grouped by cell_line_group,
#   cell_lineage, tumour, tum_lineage.  Precision@k = n_same_topk / n_topk.
#   Bootstrap uncertainty: B = 2000, sampling with replacement over biological
#   cell line groups; bootstrap statistic = mean Precision@k; lower = 2.5th
#   percentile, upper = 97.5th percentile.  Fixed seed = 20260603.
#
# Inputs (read from the canonical pipeline outputs):
#   - cellline_tumour_scores_long.tsv.gz    long table of cell line by tumour
#                                          Spearman correlation scores
#   - tumour_cellline_group_scores_long.tsv.gz    optional, used for reciprocal
#                                          consistency
#   - tumour_components.tsv                tumour graph component assignments
#
# Outputs (written to cellline_similarity_precision_bootstrap/):
#   - cellline_centred_rank_summary.tsv
#   - cellline_topk_metrics.tsv
#   - cellline_topk_lineage_composition.tsv
#   - cellline_component_mapping_summary.tsv
#   - replicate_collapse_mapping.tsv
#   - reciprocal_mapping_summary.tsv (if reciprocal scores are available)
#   - Fig_cellline_to_tumour_top1_lineage_agreement.{pdf,png}
#   - Fig_cellline_to_tumour_precision_at_k.{pdf,png}
#   - Fig_cellline_to_tumour_same_lineage_rank_percentile.{pdf,png}
#   - Fig_cellline_to_tumour_top50_lineage_composition.{pdf,png}
#   - Fig_cellline_to_tumour_top50_component_composition.{pdf,png}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#   - cellline_centred_QC_log.txt


# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(svglite)
})

# -----------------------------------------------------------------------------
# Command-line argument parsing and path configuration
# -----------------------------------------------------------------------------
#
# Usage:
#   Rscript cellline_centred_precision_at_k.R \
#     --mapping-dir <path>           parent of cellline_to_tumour_similarity/
#                                    and tumour_to_cellline_similarity/
#                                    (default: results/unsupervised/pan_cancer/
#                                               tumour_mapping)
#     [--outdir <path>]              output directory
#                                    (default: <mapping-dir>/
#                                     cellline_to_tumour_similarity/
#                                     cellline_similarity_precision_bootstrap)
#     [--tumour-components <path>]   optional TSV with columns: tumour, component
#                                    component-dependent annotations skipped
#                                    if absent; Precision@k bootstrap unaffected
#
get_arg <- function(flag, default = NA_character_) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (is.na(i) || i >= length(a)) default else a[[i + 1L]]
}

mapping_dir <- get_arg("--mapping-dir",
                        "results/unsupervised/pan_cancer/tumour_mapping")
outdir_arg  <- get_arg("--outdir")
tc_comp_arg <- get_arg("--tumour-components")

C2T_DIR <- file.path(mapping_dir, "cellline_to_tumour_similarity")
T2C_DIR <- file.path(mapping_dir, "tumour_to_cellline_similarity")

# Output directory: named after the statistical method, not the manuscript use.
# Default: cellline_similarity_precision_bootstrap/ inside the c2t directory.
OUT_DIR <- if (!is.na(outdir_arg)) outdir_arg else
             file.path(C2T_DIR, "cellline_similarity_precision_bootstrap")

C2T_LONG_FILE <- file.path(C2T_DIR, "cellline_tumour_scores_long.tsv.gz")
T2C_LONG_FILE <- file.path(T2C_DIR, "tumour_cellline_group_scores_long.tsv.gz")

# tumour_components.tsv is optional.  When absent, component-dependent
# annotations (Panel D labels, top1_tumour_component) are omitted;
# Precision@k bootstrap analysis (Panel B) is unaffected.
if (!is.na(tc_comp_arg)) {
  TUMOUR_COMPONENTS_FILE <- tc_comp_arg
} else {
  TUMOUR_COMPONENTS_FILE <- file.path(dirname(mapping_dir), "tumour_components.tsv")
}
HAS_COMPONENTS <- file.exists(TUMOUR_COMPONENTS_FILE)
BOOTSTRAP_SEED <- 20260603L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

QC_LOG <- file.path(OUT_DIR, "cellline_centred_QC_log.txt")
qc_con <- file(QC_LOG, open = "w")
qc <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = qc_con)
}
on.exit(close(qc_con), add = TRUE)

die <- function(...) stop(paste0(...), call. = FALSE)

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) die("Missing required input (", label, "): ", path)
}

fread_tsv <- function(path, ...) {
  if (grepl("\\.gz$", path)) {
    return(fread(cmd = paste("gzip -dc", shQuote(path)), ...))
  }
  fread(path, ...)
}

assert_cols <- function(dt, cols, label = "table") {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0) {
    die("Missing required columns in ", label, ": ", paste(missing, collapse = ", "))
  }
}

# -----------------------------------------------------------------------------
# Visual style
# -----------------------------------------------------------------------------
# Okabe-Ito colour blind safe palette mapped to lineages.
LINEAGE_PALETTE <- c(
  BRCA = "#D55E00",   # vermilion
  NBL  = "#0072B2",   # blue
  RBL  = "#009E73",   # bluish green
  HEME = "#CC79A7",   # reddish purple
  Other = "#999999"
)
EVAL_LINEAGES <- c("BRCA", "NBL", "RBL")

cellline_to_tumour_ranking_plot_config <- list(
  png_dpi = 300L,
  base_font_size = 15,
  axis_title_font_size = 17,
  axis_text_font_size = 14,
  legend_title_font_size = 14,
  legend_text_font_size = 14,
  facet_label_font_size = 15,
  title_font_size = 18,
  tile_label_size = 5.2,
  line_width = 1.05,
  point_size = 2.8,
  jitter_width = 0.18,
  small_point_size = 0.65,
  bar_width = 0.86,
  top1_width_in = 7.2,
  top1_height_in = 6.2,
  precision_width_in = 8.2,
  precision_height_in = 5.8,
  rank_percentile_width_in = 9.5,
  rank_percentile_height_in = 14.5,
  top50_lineage_width_in = 9.5,
  top50_lineage_height_in = 14.5
)

cellline_to_tumour_component_plot_config <- list(
  output_pdf = "Fig_cellline_to_tumour_top50_component_composition.pdf",
  output_png = "Fig_cellline_to_tumour_top50_component_composition.png",
  png_dpi = 300L,
  base_font_size = 15,
  axis_title_font_size = 17,
  axis_text_font_size = 13,
  legend_title_font_size = 14,
  legend_text_font_size = 13,
  facet_label_font_size = 15,
  title_font_size = 18,
  bar_width = 0.86,
  plot_width_in = 12.0,
  plot_height_in = 14.5
)

theme_plos <- function(base_size = cellline_to_tumour_ranking_plot_config$base_font_size) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(
        face = "bold",
        size = cellline_to_tumour_ranking_plot_config$facet_label_font_size
      ),
      legend.key = element_blank(),
      axis.text = element_text(size = cellline_to_tumour_ranking_plot_config$axis_text_font_size),
      axis.title = element_text(
        face = "bold",
        size = cellline_to_tumour_ranking_plot_config$axis_title_font_size
      ),
      legend.text = element_text(size = cellline_to_tumour_ranking_plot_config$legend_text_font_size),
      legend.title = element_text(
        face = "bold",
        size = cellline_to_tumour_ranking_plot_config$legend_title_font_size
      ),
      plot.title = element_text(
        face = "bold",
        size = cellline_to_tumour_ranking_plot_config$title_font_size
      ),
      plot.subtitle = element_text(size = cellline_to_tumour_ranking_plot_config$axis_text_font_size),
      plot.tag = element_blank(),
      plot.margin = margin(t = 12, r = 18, b = 12, l = 16)
    )
}

save_figure <- function(p, base_path, width, height,
                        dpi = cellline_to_tumour_ranking_plot_config$png_dpi) {
  ggsave(paste0(base_path, ".pdf"), p, width = width, height = height,
         device = grDevices::pdf, units = "in")
  ggsave(paste0(base_path, ".png"), p, width = width, height = height,
         dpi = dpi, units = "in")
  qc("  wrote ", base_path, ".{pdf,png}")
}

# =============================================================================
# Step 1: load and validate inputs
# =============================================================================

qc("=========================================================================")
qc("Cell line centred Precision@k similarity analysis")
qc("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
qc("Working directory: ", normalizePath("."))
qc("Bootstrap seed: ", BOOTSTRAP_SEED)
qc("=========================================================================")

require_file(C2T_LONG_FILE, "cellline by tumour long scores")
if (!HAS_COMPONENTS) {
  qc("WARNING: tumour_components.tsv not found at: ", TUMOUR_COMPONENTS_FILE)
  qc("  Component annotations and Panel D component labels will be omitted.")
  qc("  Precision@k bootstrap analysis (Panel B) is computed regardless.")
}

scores_long <- fread_tsv(C2T_LONG_FILE)
assert_cols(scores_long,
            c("cell_line", "cell_lineage", "tumour", "tum_lineage", "score"),
            "cellline_tumour_scores_long")

if (anyNA(scores_long$score)) die("NA scores present; cannot rank.")
if (anyNA(scores_long$cell_lineage)) die("Missing cell_lineage for some cell lines.")
if (anyNA(scores_long$tum_lineage)) die("Missing tum_lineage for some tumours.")

if (HAS_COMPONENTS) {
  components <- fread(TUMOUR_COMPONENTS_FILE)
  if (!"tumour" %in% names(components) && "sample" %in% names(components)) {
    setnames(components, "sample", "tumour")
  }
  assert_cols(components, c("tumour", "component"), "tumour_components")
  components <- components[, .(tumour, component)]
} else {
  components <- data.table(tumour = character(0), component = integer(0))
}

# -----------------------------------------------------------------------------
# Replicate collapse mapping
#   The raw cell line identifiers carry a sequencing prefix (NG-<run>_) and a
#   profile-level suffix (_lib<id>_<run>_<lane>). Stripping these yields a
#   biological cell line group label. Repeated profile-level observations
#   collapse naturally.
# -----------------------------------------------------------------------------

derive_cell_line_group <- function(x) {
  g <- sub("^NG-[A-Za-z0-9]+_", "", x)
  g <- sub("_lib.*$", "", g)
  g
}

cellline_meta <- unique(scores_long[, .(cell_line_raw = cell_line, cell_lineage)])
cellline_meta[, cell_line_group := derive_cell_line_group(cell_line_raw)]

# Sanity: groups must be lineage consistent
group_lineage_check <- cellline_meta[, .(n_lineages = uniqueN(cell_lineage)),
                                     by = cell_line_group][n_lineages > 1]
if (nrow(group_lineage_check) > 0) {
  die("Cell line group spans multiple lineages: ",
      paste(group_lineage_check$cell_line_group, collapse = ", "))
}

replicate_counts <- cellline_meta[, .N, by = .(cell_line_group, cell_lineage)]
multi_replicate <- replicate_counts[N > 1]

cellline_meta[, replicate_status := ifelse(
  cell_line_group %in% multi_replicate$cell_line_group,
  "multi_profile", "single_profile"
)]
cellline_meta[, collapse_rule := "mean_score_across_replicate_profiles"]
cellline_meta[, retained_for_primary_analysis := TRUE]

setcolorder(cellline_meta,
            c("cell_line_raw", "cell_line_group", "cell_lineage",
              "replicate_status", "collapse_rule",
              "retained_for_primary_analysis"))

fwrite(cellline_meta, file.path(OUT_DIR, "replicate_collapse_mapping.tsv"),
       sep = "\t")
qc("Replicate mapping written to replicate_collapse_mapping.tsv")

n_raw <- uniqueN(cellline_meta$cell_line_raw)
n_groups <- uniqueN(cellline_meta$cell_line_group)
n_tumours <- uniqueN(scores_long$tumour)

qc("")
qc("Cohort sizes:")
qc("  raw cell line entries        : ", n_raw)
qc("  biological cell line groups  : ", n_groups)
qc("  tumour samples (denominator) : ", n_tumours)

per_lineage_cl <- unique(cellline_meta[, .(cell_line_group, cell_lineage)])[, .N,
                       by = cell_lineage][order(cell_lineage)]
per_lineage_t <- unique(scores_long[, .(tumour, tum_lineage)])[, .N,
                       by = tum_lineage][order(tum_lineage)]
qc("Cell line groups per lineage:")
for (i in seq_len(nrow(per_lineage_cl))) {
  qc("  ", per_lineage_cl$cell_lineage[i], " : ", per_lineage_cl$N[i])
}
qc("Tumour samples per lineage:")
for (i in seq_len(nrow(per_lineage_t))) {
  qc("  ", per_lineage_t$tum_lineage[i], " : ", per_lineage_t$N[i])
}

# -----------------------------------------------------------------------------
# Aggregate replicate profiles and recompute ranks
# -----------------------------------------------------------------------------

if (!"cell_line_group" %in% names(scores_long)) {
  scores_long <- merge(
    scores_long,
    cellline_meta[, .(cell_line_raw, cell_line_group)],
    by.x = "cell_line", by.y = "cell_line_raw", all.x = TRUE
  )
}

agg <- scores_long[, .(score = mean(score)),
                   by = .(cell_line_group, cell_lineage, tumour, tum_lineage)]

# Rank tumours per cell line group with ties.method = "average"
agg[, rank_all := frank(-score, ties.method = "average"),
    by = cell_line_group]

# Tie diagnostics per cell line group
tie_summary <- agg[, .(
  n_tied_groups = sum(duplicated(score) | duplicated(score, fromLast = TRUE)) /
                    pmax(1, sum(duplicated(score) | duplicated(score, fromLast = TRUE)))
), by = cell_line_group]

tie_groups <- agg[, {
  ord <- order(-score)
  s <- score[ord]
  rl <- rle(s)
  n_tied_groups <- sum(rl$lengths > 1)
  n_affected_ranks <- sum(rl$lengths[rl$lengths > 1])
  list(n_tied_groups = n_tied_groups, n_affected_ranks = n_affected_ranks)
}, by = cell_line_group]

qc("")
qc("Tie statistics across cell line groups:")
qc("  cell line groups with at least one tied score group: ",
   sum(tie_groups$n_tied_groups > 0))
qc("  total tied score groups across all cell lines      : ",
   sum(tie_groups$n_tied_groups))
qc("  total ranks affected by ties                        : ",
   sum(tie_groups$n_affected_ranks))

# Evaluation denominator
# In the present pipeline only BRCA, NBL, RBL tumours are scored.
# rank_all and rank_eval are therefore identical and reported as such.
eval_lineages_present <- intersect(EVAL_LINEAGES, unique(agg$tum_lineage))
agg[, rank_eval := rank_all]

n_total_all <- uniqueN(agg$tumour)
n_total_eval <- agg[tum_lineage %in% eval_lineages_present, uniqueN(tumour)]

qc("Ranking denominator:")
qc("  n_total_tumours_all  : ", n_total_all)
qc("  n_total_tumours_eval : ", n_total_eval, "  (eval lineages: ",
   paste(eval_lineages_present, collapse = ", "), ")")
if (n_total_all == n_total_eval) {
  qc("  note: rank_all == rank_eval because every scored tumour belongs to ",
     "an evaluation lineage.")
}

agg[, rank_percentile := rank_all / n_total_all]
agg[, same_lineage := cell_lineage == tum_lineage]

# Attach tumour components
agg <- merge(agg, components, by = "tumour", all.x = TRUE)

# Component annotation for top-50 tumour-neighbourhood composition plots
component_anno <- agg[, .(
  n_tumours = uniqueN(tumour),
  dominant_lineage = names(sort(table(tum_lineage), decreasing = TRUE))[1],
  lineage_purity = max(prop.table(table(tum_lineage)))
), by = component][order(-n_tumours)]
component_anno[, short_label := paste0("Component ", component)]
agg <- merge(agg, component_anno[, .(component, short_label)],
             by = "component", all.x = TRUE)
setnames(agg, "short_label", "component_label")

qc("")
qc("Tumour components: ", nrow(component_anno),
   " components covering ", sum(component_anno$n_tumours), " tumours")

# =============================================================================
# Step 2: per cell line summary metrics
# =============================================================================

iqr_lo <- function(x) as.numeric(quantile(x, 0.25, na.rm = TRUE))
iqr_hi <- function(x) as.numeric(quantile(x, 0.75, na.rm = TRUE))

per_cl <- agg[, {
  setorder(.SD, rank_all)
  top1 <- .SD[1]
  top2 <- if (.N >= 2) .SD[2]$score else NA_real_
  top10 <- if (.N >= 10) .SD[10]$score else NA_real_
  same_ranks <- .SD[same_lineage == TRUE]$rank_all
  same_pct <- same_ranks / n_total_all
  n_same <- length(same_ranks)
  list(
    cell_lineage = top1$cell_lineage,
    n_total_tumours_all = .N,
    n_total_tumours_eval = sum(.SD$tum_lineage %in% eval_lineages_present),
    n_same_lineage_tumours = n_same,
    top1_tumour = top1$tumour,
    top1_tumour_lineage = top1$tum_lineage,
    top1_score = top1$score,
    top1_same_lineage = top1$same_lineage,
    top1_rank = top1$rank_all,
    top1_rank_percentile = top1$rank_percentile,
    top1_tumour_component = top1$component,
    top1_component_label = top1$component_label,
    best_same_lineage_rank = if (n_same > 0) min(same_ranks) else NA_real_,
    median_same_lineage_rank = if (n_same > 0) median(same_ranks) else NA_real_,
    same_lineage_rank_IQR_lower = if (n_same > 0) iqr_lo(same_ranks) else NA_real_,
    same_lineage_rank_IQR_upper = if (n_same > 0) iqr_hi(same_ranks) else NA_real_,
    median_same_lineage_rank_percentile = if (n_same > 0) median(same_pct) else NA_real_,
    same_lineage_rank_percentile_IQR_lower = if (n_same > 0) iqr_lo(same_pct) else NA_real_,
    same_lineage_rank_percentile_IQR_upper = if (n_same > 0) iqr_hi(same_pct) else NA_real_,
    delta_top1_top2 = top1$score - top2,
    delta_top1_top10 = top1$score - top10,
    score_rank2 = top2,
    score_rank10 = top10,
    score_IQR_top10 = if (.N >= 10) IQR(.SD[1:10]$score) else NA_real_,
    fraction_same_lineage_in_top_1pct = if (n_same > 0) mean(same_pct <= 0.01) else NA_real_,
    fraction_same_lineage_in_top_5pct = if (n_same > 0) mean(same_pct <= 0.05) else NA_real_,
    fraction_same_lineage_in_top_10pct = if (n_same > 0) mean(same_pct <= 0.10) else NA_real_,
    fraction_same_lineage_in_top_20pct = if (n_same > 0) mean(same_pct <= 0.20) else NA_real_
  )
}, by = cell_line_group]

# Add a representative raw label (first replicate)
raw_examples <- cellline_meta[, .(cell_line_raw_examples = paste(sort(unique(cell_line_raw)),
                                                                 collapse = ";")),
                              by = cell_line_group]
per_cl <- merge(raw_examples, per_cl, by = "cell_line_group")

setcolorder(per_cl,
            c("cell_line_group", "cell_line_raw_examples", "cell_lineage",
              "n_total_tumours_all", "n_total_tumours_eval",
              "n_same_lineage_tumours",
              "top1_tumour", "top1_tumour_lineage", "top1_score",
              "top1_same_lineage", "top1_rank", "top1_rank_percentile",
              "best_same_lineage_rank", "median_same_lineage_rank",
              "same_lineage_rank_IQR_lower", "same_lineage_rank_IQR_upper",
              "median_same_lineage_rank_percentile",
              "same_lineage_rank_percentile_IQR_lower",
              "same_lineage_rank_percentile_IQR_upper",
              "delta_top1_top2", "delta_top1_top10",
              "top1_tumour_component", "top1_component_label"))

fwrite(per_cl, file.path(OUT_DIR, "cellline_centred_rank_summary.tsv"),
       sep = "\t")
qc("Wrote cellline_centred_rank_summary.tsv (", nrow(per_cl), " rows)")

# Top 1 accuracy with Wilson score CI, overall and per lineage
wilson_ci <- function(success, total, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  phat <- success / total
  denom <- 1 + z^2 / total
  centre <- (phat + z^2 / (2 * total)) / denom
  half_width <- z * sqrt(phat * (1 - phat) / total + z^2 / (4 * total^2)) / denom
  c(max(0, centre - half_width), min(1, centre + half_width))
}

binom_row <- function(success, total, label) {
  if (total == 0) return(data.table(group = label, n_correct = 0, n_total = 0,
                                    accuracy = NA_real_, ci_lower = NA_real_,
                                    ci_upper = NA_real_))
  ci <- wilson_ci(success, total)
  data.table(group = label, n_correct = success, n_total = total,
             accuracy = success / total, ci_lower = ci[1], ci_upper = ci[2])
}

top1_overall <- binom_row(sum(per_cl$top1_same_lineage),
                          nrow(per_cl), "overall")
top1_per_lin <- per_cl[, binom_row(sum(top1_same_lineage), .N, cell_lineage),
                       by = cell_lineage]
top1_per_lin[, cell_lineage := NULL]
top1_table <- rbind(top1_overall, top1_per_lin)
fwrite(top1_table, file.path(OUT_DIR, "cellline_top1_accuracy.tsv"), sep = "\t")

# Balanced accuracy = mean of per-lineage accuracies
bal_acc <- mean(top1_per_lin$accuracy)
qc("")
qc("Top 1 lineage accuracy:")
for (i in seq_len(nrow(top1_table))) {
  r <- top1_table[i]
  qc(sprintf("  %-7s : %d / %d = %.3f  [95%% CI %.3f, %.3f]",
             r$group, r$n_correct, r$n_total, r$accuracy,
             r$ci_lower, r$ci_upper))
}
qc(sprintf("  balanced accuracy across lineages : %.3f", bal_acc))

# =============================================================================
# Step 3: top k metrics, enrichment, FDR
# =============================================================================

abs_k <- c(1, 5, 10, 25, 50, 100)
abs_k <- abs_k[abs_k <= n_total_all]
pct_k <- c(0.01, 0.05, 0.10)
k_table <- rbind(
  data.table(k = abs_k, k_type = "absolute"),
  data.table(k = ceiling(pct_k * n_total_all), k_type = "percentile")
)
k_table <- unique(k_table[order(k)])

topk_metrics <- vector("list", 0)
topk_compositions <- vector("list", 0)
topk_components <- vector("list", 0)

for (i in seq_len(nrow(k_table))) {
  ki <- k_table$k[i]
  ktype <- k_table$k_type[i]

  # Per cell line
  for (cl in per_cl$cell_line_group) {
    sub <- agg[cell_line_group == cl][order(rank_all)]
    top <- sub[rank_all <= ki]
    n_topk <- nrow(top)
    n_same_topk <- sum(top$same_lineage)
    n_same_total <- per_cl[cell_line_group == cl]$n_same_lineage_tumours
    expected_frac <- n_same_total / n_total_all
    precision <- n_same_topk / n_topk
    recall <- if (n_same_total > 0) n_same_topk / n_same_total else NA_real_
    enrichment <- if (expected_frac > 0) precision / expected_frac else NA_real_
    # Hypergeometric: q = same in top, m = same total, n = other total, k = topk
    p_hg <- phyper(n_same_topk - 1, n_same_total,
                   n_total_all - n_same_total, n_topk, lower.tail = FALSE)
    null_expected_count <- n_topk * expected_frac
    topk_metrics[[length(topk_metrics) + 1]] <- data.table(
      cell_line_group = cl,
      cell_lineage = per_cl[cell_line_group == cl]$cell_lineage,
      k = ki, k_type = ktype,
      n_topk = n_topk,
      n_same_lineage_topk = n_same_topk,
      precision_at_k = precision,
      recall_at_k = recall,
      expected_same_lineage_fraction = expected_frac,
      null_expected_count = null_expected_count,
      observed_count = n_same_topk,
      enrichment_at_k = enrichment,
      p_value_hypergeometric = p_hg
    )
    # Lineage composition of the top k
    comp_lin <- top[, .(n_tumours = .N), by = tum_lineage]
    comp_lin[, fraction_topk := n_tumours / n_topk]
    comp_lin[, k := ki][, k_type := ktype]
    comp_lin[, cell_line_group := cl]
    comp_lin[, cell_lineage := per_cl[cell_line_group == cl]$cell_lineage]
    setnames(comp_lin, "tum_lineage", "tumour_lineage")
    topk_compositions[[length(topk_compositions) + 1]] <- comp_lin

    # Component composition of the top k
    comp_top <- top[, .(n_topk_tumours = .N), by = .(component, component_label)]
    comp_top[, fraction_topk := n_topk_tumours / n_topk]
    comp_top[, k := ki][, k_type := ktype]
    comp_top[, cell_line_group := cl]
    comp_top[, cell_lineage := per_cl[cell_line_group == cl]$cell_lineage]
    setnames(comp_top, c("component", "component_label"),
             c("tumour_component", "component_label"))
    topk_components[[length(topk_components) + 1]] <- comp_top
  }
}

topk_metrics <- rbindlist(topk_metrics)
topk_metrics[, q_value_BH := p.adjust(p_value_hypergeometric, method = "BH")]
setcolorder(topk_metrics,
            c("cell_line_group", "cell_lineage", "k", "k_type",
              "n_topk", "n_same_lineage_topk",
              "precision_at_k", "recall_at_k",
              "expected_same_lineage_fraction",
              "null_expected_count", "observed_count",
              "enrichment_at_k",
              "p_value_hypergeometric", "q_value_BH"))
fwrite(topk_metrics, file.path(OUT_DIR, "cellline_topk_metrics.tsv"),
       sep = "\t")
qc("Wrote cellline_topk_metrics.tsv (", nrow(topk_metrics), " rows)")

topk_compositions <- rbindlist(topk_compositions, use.names = TRUE, fill = TRUE)
setcolorder(topk_compositions,
            c("cell_line_group", "cell_lineage", "k", "k_type",
              "tumour_lineage", "n_tumours", "fraction_topk"))
fwrite(topk_compositions,
       file.path(OUT_DIR, "cellline_topk_lineage_composition.tsv"), sep = "\t")
qc("Wrote cellline_topk_lineage_composition.tsv (", nrow(topk_compositions), " rows)")

topk_components <- rbindlist(topk_components, use.names = TRUE, fill = TRUE)
setcolorder(topk_components,
            c("cell_line_group", "cell_lineage", "k", "k_type",
              "tumour_component", "component_label",
              "n_topk_tumours", "fraction_topk"))
fwrite(topk_components,
       file.path(OUT_DIR, "cellline_component_mapping_summary.tsv"), sep = "\t")
qc("Wrote cellline_component_mapping_summary.tsv (", nrow(topk_components), " rows)")

# =============================================================================
# Step 4: reciprocal consistency
# =============================================================================

reciprocal_summary_path <- file.path(OUT_DIR, "reciprocal_mapping_summary.tsv")
have_reciprocal <- file.exists(T2C_LONG_FILE)

if (have_reciprocal) {
  qc("")
  qc("Computing reciprocal cell line tumour matching ...")
  t2c_long <- fread_tsv(T2C_LONG_FILE)
  assert_cols(t2c_long,
              c("tumour", "tum_lineage", "cell_line", "cell_lineage", "score"),
              "tumour_cellline_scores_long")

  if (!"cell_line_group" %in% names(t2c_long)) {
    t2c_long <- merge(
      t2c_long,
      cellline_meta[, .(cell_line_raw, cell_line_group)],
      by.x = "cell_line", by.y = "cell_line_raw", all.x = TRUE
    )
  }
  t2c_agg <- t2c_long[, .(score = mean(score)),
                      by = .(tumour, tum_lineage, cell_line_group, cell_lineage)]
  t2c_agg[, rank_t2c := frank(-score, ties.method = "average"), by = tumour]
  n_groups_total <- uniqueN(t2c_agg$cell_line_group)
  t2c_agg[, rank_percentile_t2c := rank_t2c / n_groups_total]

  # Join on (cell_line_group, tumour)
  recip <- merge(
    agg[, .(cell_line_group, cell_lineage, tumour, tum_lineage,
            cellline_to_tumour_rank = rank_all,
            cellline_to_tumour_rank_percentile = rank_percentile,
            score_c2t = score)],
    t2c_agg[, .(cell_line_group, tumour,
                tumour_to_cellline_rank = rank_t2c,
                tumour_to_cellline_rank_percentile = rank_percentile_t2c,
                score_t2c = score)],
    by = c("cell_line_group", "tumour")
  )
  recip[, reciprocal_top10 := cellline_to_tumour_rank <= 10 &
          tumour_to_cellline_rank <= 10]
  recip[, reciprocal_top10_percent := cellline_to_tumour_rank_percentile <= 0.10 &
          tumour_to_cellline_rank_percentile <= 0.10]
  recip[, reciprocal_rank_score := 1 / cellline_to_tumour_rank +
          1 / tumour_to_cellline_rank]
  setcolorder(recip,
              c("cell_line_group", "tumour", "cell_lineage", "tum_lineage",
                "cellline_to_tumour_rank", "tumour_to_cellline_rank",
                "cellline_to_tumour_rank_percentile",
                "tumour_to_cellline_rank_percentile",
                "reciprocal_top10", "reciprocal_top10_percent",
                "reciprocal_rank_score"))
  fwrite(recip, reciprocal_summary_path, sep = "\t")
  qc("Wrote reciprocal_mapping_summary.tsv (", nrow(recip), " rows)")

  # Per cell line summary
  recip_per_cl <- recip[, .(
    n_reciprocal_top10_pairs = sum(reciprocal_top10),
    fraction_same_lineage_reciprocal_pairs = if (sum(reciprocal_top10) > 0)
      mean(cell_lineage[reciprocal_top10] == tum_lineage[reciprocal_top10])
    else NA_real_,
    median_reciprocal_rank_score = median(reciprocal_rank_score)
  ), by = cell_line_group]
  fwrite(recip_per_cl, file.path(OUT_DIR, "reciprocal_per_cell_line.tsv"),
         sep = "\t")
} else {
  qc("Tumour to cell line long score table not found at ", T2C_LONG_FILE,
     " - skipping reciprocal analysis.")
}

# =============================================================================
# Step 5: retained standalone figures
# =============================================================================

# Order of cell line groups for figures: by lineage then by median same-lineage
# rank percentile.
order_cl <- per_cl[order(cell_lineage, median_same_lineage_rank_percentile)]
cl_levels <- order_cl$cell_line_group

agg[, cell_line_group := factor(cell_line_group, levels = cl_levels)]
per_cl[, cell_line_group := factor(cell_line_group, levels = cl_levels)]
per_cl[, cell_lineage := factor(cell_lineage, levels = EVAL_LINEAGES)]

# Top-ranked tumour-lineage agreement by biological cell-line group lineage
conf <- per_cl[, .N, by = .(cell_lineage, top1_tumour_lineage)]
conf[, row_total := sum(N), by = cell_lineage]
conf[, row_pct := N / row_total]

all_lin <- unique(c(as.character(per_cl$cell_lineage),
                    unique(per_cl$top1_tumour_lineage)))
conf_full <- CJ(cell_lineage = all_lin, top1_tumour_lineage = all_lin)
conf <- merge(conf_full, conf,
              by = c("cell_lineage", "top1_tumour_lineage"),
              all.x = TRUE)
conf[is.na(N), N := 0]
conf[, row_total := NULL]
conf[, row_total := sum(N), by = cell_lineage]
conf[, row_pct := ifelse(row_total > 0, N / row_total, 0)]

p_top1 <- ggplot(conf, aes(x = top1_tumour_lineage, y = cell_lineage,
                           fill = row_pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d\n%.0f%%", N, 100 * row_pct)),
    size = cellline_to_tumour_ranking_plot_config$tile_label_size
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#08306B",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    name = "Row fraction"
  ) +
  scale_x_discrete(limits = all_lin) +
  scale_y_discrete(limits = rev(all_lin)) +
  labs(
    title = "Top-ranked tumour lineage by cell-line lineage",
    x = "Top-ranked tumour lineage",
    y = "Cell-line lineage"
  ) +
  theme_plos() +
  theme(legend.position = "right", panel.grid = element_blank())

save_figure(
  p_top1,
  file.path(OUT_DIR, "Fig_cellline_to_tumour_top1_lineage_agreement"),
  width = cellline_to_tumour_ranking_plot_config$top1_width_in,
  height = cellline_to_tumour_ranking_plot_config$top1_height_in
)

# Precision at k by cell-line lineage
panB_data <- topk_metrics[k_type == "absolute"]

set.seed(BOOTSTRAP_SEED)
boot_mean_ci <- function(x, B = 2000) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(mean = mean(x), lo = NA_real_, hi = NA_real_))
  bs <- replicate(B, mean(sample(x, replace = TRUE)))
  c(mean = mean(x),
    lo = as.numeric(quantile(bs, 0.025)),
    hi = as.numeric(quantile(bs, 0.975)))
}

panB_overall <- panB_data[, {
  ci <- boot_mean_ci(precision_at_k)
  list(group = "Overall", precision_mean = ci["mean"], lo = ci["lo"], hi = ci["hi"])
}, by = .(k)]

panB_per_lin <- panB_data[, {
  ci <- boot_mean_ci(precision_at_k)
  list(group = as.character(cell_lineage[1]),
       precision_mean = ci["mean"], lo = ci["lo"], hi = ci["hi"])
}, by = .(k, cell_lineage)]
panB_per_lin[, cell_lineage := NULL]

panB_combined <- rbind(panB_overall, panB_per_lin)
panB_combined[, group := factor(group, levels = c("Overall", EVAL_LINEAGES))]
curve_palette <- c(Overall = "#000000", LINEAGE_PALETTE[EVAL_LINEAGES])

p_precision <- ggplot(panB_combined, aes(x = k, y = precision_mean,
                                         colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = cellline_to_tumour_ranking_plot_config$line_width) +
  geom_point(size = cellline_to_tumour_ranking_plot_config$point_size) +
  scale_x_continuous(
    breaks = sort(unique(panB_data$k)),
    trans = "log10",
    labels = function(x) formatC(x, format = "d")
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = curve_palette, name = "Cell-line lineage") +
  scale_fill_manual(values = curve_palette, name = "Cell-line lineage") +
  labs(
    title = "Precision at k by cell-line lineage",
    x = "Top k tumours",
    y = "Mean Precision at k"
  ) +
  theme_plos() +
  theme(legend.position = "right")

save_figure(
  p_precision,
  file.path(OUT_DIR, "Fig_cellline_to_tumour_precision_at_k"),
  width = cellline_to_tumour_ranking_plot_config$precision_width_in,
  height = cellline_to_tumour_ranking_plot_config$precision_height_in
)

# Same-lineage rank percentile distribution per biological cell-line group
panC_data <- agg[same_lineage == TRUE]
panC_data[, cell_line_group := factor(cell_line_group, levels = cl_levels)]

p_rank_percentile <- ggplot(panC_data, aes(x = cell_line_group, y = rank_percentile,
                                           fill = cell_lineage)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.45, alpha = 0.60) +
  geom_jitter(
    width = cellline_to_tumour_ranking_plot_config$jitter_width,
    alpha = 0.35,
    size = cellline_to_tumour_ranking_plot_config$small_point_size,
    colour = "grey20"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = LINEAGE_PALETTE, name = "Cell-line lineage") +
  labs(
    title = "Same-lineage rank percentile by cell-line group",
    x = NULL,
    y = "Same-lineage rank percentile"
  ) +
  coord_flip() +
  theme_plos() +
  theme(axis.text.y = element_text(size = cellline_to_tumour_ranking_plot_config$axis_text_font_size))

save_figure(
  p_rank_percentile,
  file.path(OUT_DIR, "Fig_cellline_to_tumour_same_lineage_rank_percentile"),
  width = cellline_to_tumour_ranking_plot_config$rank_percentile_width_in,
  height = cellline_to_tumour_ranking_plot_config$rank_percentile_height_in
)

# Top-50 tumour-lineage composition per biological cell-line group
k_for_D <- 50L
if (!(k_for_D %in% topk_compositions$k)) k_for_D <- max(topk_compositions$k)

panD_data <- topk_compositions[k == k_for_D & k_type == "absolute"]
panD_data[, tumour_lineage := factor(tumour_lineage,
                                      levels = c(EVAL_LINEAGES, "HEME", "Other"))]
prec_at_kD <- topk_metrics[k == k_for_D & k_type == "absolute",
                           .(cell_line_group, precision_at_k)]
order_d <- prec_at_kD[order(-precision_at_k)]$cell_line_group
panD_data[, cell_line_group := factor(cell_line_group, levels = order_d)]

p_top50_lineage <- ggplot(panD_data, aes(x = cell_line_group, y = fraction_topk,
                                         fill = tumour_lineage)) +
  geom_col(width = cellline_to_tumour_ranking_plot_config$bar_width) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.0001), expand = c(0, 0)) +
  scale_fill_manual(values = LINEAGE_PALETTE, name = "Tumour lineage") +
  labs(
    title = "Top-50 tumour-lineage composition by cell-line group",
    x = NULL,
    y = "Fraction of top-50 tumours"
  ) +
  coord_flip() +
  theme_plos() +
  theme(axis.text.y = element_text(size = cellline_to_tumour_ranking_plot_config$axis_text_font_size))

save_figure(
  p_top50_lineage,
  file.path(OUT_DIR, "Fig_cellline_to_tumour_top50_lineage_composition"),
  width = cellline_to_tumour_ranking_plot_config$top50_lineage_width_in,
  height = cellline_to_tumour_ranking_plot_config$top50_lineage_height_in
)

# Top-50 consensus-network component composition per biological cell-line group
k_comp <- 50L
if (!(k_comp %in% topk_components$k)) k_comp <- max(topk_components$k)

comp_data <- topk_components[k == k_comp & k_type == "absolute"]
comp_data[, cell_line_group := factor(cell_line_group, levels = order_d)]
comp_data[, component_label := factor(component_label,
                                      levels = component_anno[order(-n_tumours)]$short_label)]

component_palette <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(12, "Set3")
)(uniqueN(comp_data$component_label))

p_component <- ggplot(comp_data, aes(x = cell_line_group, y = fraction_topk,
                                     fill = component_label)) +
  geom_col(width = cellline_to_tumour_component_plot_config$bar_width) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0),
                     limits = c(0, 1.0001)) +
  scale_fill_manual(values = component_palette, name = "Consensus-network component") +
  labs(
    title = "Top-50 tumour composition by consensus-network component",
    x = NULL,
    y = "Fraction of top-50 tumours"
  ) +
  coord_flip() +
  theme_plos(base_size = cellline_to_tumour_component_plot_config$base_font_size) +
  theme(
    plot.title = element_text(
      size = cellline_to_tumour_component_plot_config$title_font_size,
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold",
      size = cellline_to_tumour_component_plot_config$axis_title_font_size
    ),
    axis.text.x = element_text(size = cellline_to_tumour_component_plot_config$axis_text_font_size),
    axis.text.y = element_text(size = cellline_to_tumour_component_plot_config$axis_text_font_size),
    legend.title = element_text(
      face = "bold",
      size = cellline_to_tumour_component_plot_config$legend_title_font_size
    ),
    legend.text = element_text(size = cellline_to_tumour_component_plot_config$legend_text_font_size),
    legend.key.size = unit(0.42, "cm")
  )

save_figure(
  p_component,
  file.path(OUT_DIR, sub("\\.pdf$", "", cellline_to_tumour_component_plot_config$output_pdf)),
  width = cellline_to_tumour_component_plot_config$plot_width_in,
  height = cellline_to_tumour_component_plot_config$plot_height_in,
  dpi = cellline_to_tumour_component_plot_config$png_dpi
)

# Optional confidence-margin and top-10 score diagnostics remain available in
# cellline_centred_rank_summary.tsv; no retained ranking-plot output is written.

# =============================================================================
# QC summary
# =============================================================================
qc("")
qc("Output directory: ", OUT_DIR)
qc("Files written:")
for (f in list.files(OUT_DIR, full.names = FALSE)) {
  qc("  ", f)
}
qc("")
qc("Done.")
