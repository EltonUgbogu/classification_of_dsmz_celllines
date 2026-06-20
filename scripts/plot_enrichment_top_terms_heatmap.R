#!/usr/bin/env Rscript
# plot_enrichment_top_terms_heatmap.R
# Produces a ggplot2 tile-heatmap of selected enriched terms across
# disease-specific and pan-cancer marker-gene sets.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
input_file <- "/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/enrichment_summary_top_terms.tsv"
outdir     <- "/Users/eltonugbogu/Library/CloudStorage/GoogleDrive-eltontobi8@gmail.com/My Drive/dsmz/pipeline/pipeline/results/unsupervised/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load and validate ───────────────────────────────────────────────────────
dt <- fread(input_file)

required_cols <- c("query_id", "term_name", "source", "p_value")
missing_cols  <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Input table is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

cat("Loaded", nrow(dt), "rows from input table.\n")

# ── 2. Query-ID → display group label ─────────────────────────────────────────
# Internal labels (used for data processing)
map_query <- function(q) {
  if (grepl("^brca__per_contrast__", q))                                        return("BRCA per-contrast markers")
  if (grepl("^nbl__per_contrast__", q))                                         return("NBL per-contrast markers")
  if (grepl("^rbl__per_contrast__", q))                                         return("RBL per-contrast markers")
  if (grepl("^brca__disease_recurrent_ge__2__", q))                             return("BRCA recurrent")
  if (grepl("^nbl__disease_recurrent_ge__2__", q))                              return("NBL recurrent")
  if (grepl("^rbl__disease_recurrent_ge__2__", q))                              return("RBL recurrent")
  if (q == "strict_union__down")                                                 return("Strict union down")
  if (q == "strict_union__up")                                                   return("Strict union up")
  if (q == "strict_union__all")                                                  return("Strict union all")
  if (q == "operative_feature_set__all")                                         return("Operative feature set")
  if (q == "comparison_strict_vs_operative__operative_minus_strict")             return("Operative minus strict")
  # Excluded groups:
  # ordered_strict_union__down_ordered  → "Ordered strict union down"
  # ordered_strict_union__up_ordered    → "Ordered strict union up"
  return(NA_character_)
}

dt[, group_label := vapply(query_id, map_query, character(1))]
# Drop excluded groups (ordered strict union and any unmapped)
dt_mapped <- dt[!is.na(group_label)]
cat("Rows after group mapping:", nrow(dt_mapped), "\n")

# ── 3. Generic / root term removal ────────────────────────────────────────────
generic_terms <- c(
  "biological_process",
  "cellular_component",
  "molecular_function",
  "binding",
  "cellular process",
  "response to stimulus",
  "Transfac [TF]",
  "WIKIPATHWAYS [WP]",
  # root terms found in data
  "multicellular organismal process",
  "cellular anatomical structure",
  "cell periphery",
  "HPA root",
  "KEGG root term",
  "REACTOME root term",
  "Transfac",
  "WIKIPATHWAYS"
)

n_before_generic <- nrow(dt_mapped)
dt_mapped <- dt_mapped[!term_name %in% generic_terms]
n_removed_generic <- n_before_generic - nrow(dt_mapped)
cat("Removed", n_removed_generic, "rows with generic/root term_name values.\n")

# ── 4. Score: −log10(p), capped at 30 ─────────────────────────────────────────
dt_mapped[, score := pmin(-log10(pmax(p_value, 1e-300)), 30)]

# ── 5. Collapse TF motif duplicate rows (same Factor+motif, differ by "match class") ──
# Pattern: keep lowest p_value when rows share Factor name and motif but differ only by
# the "match class: 1" suffix.

# Helper: strip "match class: 1" from TF term_name for grouping
strip_match_class <- function(x) {
  gsub(";?\\s*match class:\\s*[0-9]+", "", x, perl = TRUE)
}

dt_mapped[source == "TF", tf_base := strip_match_class(term_name)]

n_before_dedup <- nrow(dt_mapped)
# For TF rows: keep the row with lowest p_value per (query_id, tf_base); keep original term_name of winner
dt_tf   <- dt_mapped[source == "TF"][order(p_value)][, .SD[1], by = .(query_id, tf_base)]
dt_notf <- dt_mapped[source != "TF"]
dt_mapped <- rbind(dt_notf, dt_tf, fill = TRUE)
n_removed_dups <- n_before_dedup - nrow(dt_mapped)
cat("Removed", n_removed_dups, "duplicate TF motif rows (kept lowest p_value per Factor+motif).\n")

# ── 6. Simplify TF labels ──────────────────────────────────────────────────────
# Input examples:
#   "Factor: Spi-B; motif: TTCYBC; match class: 1"  → "Spi-B motif (TTCYBC) [TF]"
#   "Factor: E2F-2; motif: GCGCGCGCNCS"              → "E2F-2 motif (GCGCGCGCNCS) [TF]"
simplify_tf_label <- function(x) {
  # Extract Factor name
  fac <- str_match(x, "Factor:\\s*([^;]+)")[, 2]
  # Extract motif
  mot <- str_match(x, "motif:\\s*([^;]+)")[, 2]
  ifelse(
    !is.na(fac) & !is.na(mot),
    paste0(trimws(fac), " motif (", trimws(mot), ") [TF]"),
    paste0(x, " [TF]")  # fallback
  )
}

# Source-to-bracket mapping
source_bracket <- c(
  "GO:BP"  = "[GO:BP]",
  "GO:CC"  = "[GO:CC]",
  "GO:MF"  = "[GO:MF]",
  "REAC"   = "[REAC]",
  "KEGG"   = "[KEGG]",
  "WP"     = "[WP]",
  "TF"     = "[TF]",
  "CORUM"  = "[CORUM]",
  "HPA"    = "[HPA]"
)

dt_mapped[, term_label := {
  br <- source_bracket[source]
  br[is.na(br)] <- paste0("[", source[is.na(br)], "]")
  ifelse(
    source == "TF",
    simplify_tf_label(term_name),
    paste0(term_name, " ", br)
  )
}]

# ── 7. Prioritised term selection for main figure (≤30 rows) ──────────────────
# Priority list: terms mentioned in the report (or closely matching)
priority_terms <- c(
  # NBL: IL-10, inflammatory/defense, cytokine, neuron
  "Interleukin-10 signaling",
  "Interleukin-10 signalling",
  "inflammatory response",
  "defense response",
  "immune response",
  "Cytokine Signaling in Immune system",
  "cytokine activity",
  "Cytokine-cytokine receptor interaction",
  "neuron development",
  "neuron differentiation",
  "modulation of chemical synaptic transmission",
  # RBL: signalling/cell communication, basement membrane/collagen, photoreceptor, CORUM
  "signaling",
  "cell communication",
  "signal transduction",
  "network-forming collagen trimer",
  "collagen type IV trimer",
  "collagenous component of basement membrane",
  "glomerular basement membrane development",
  "Anchoring fibril formation",
  "Crosslinking of collagen fibrils",
  "Retina; photoreceptor cells[≥Low]",
  "Retina",
  "CIN85-SH3GL3 complex",
  "PXN-ITGB5-PTK2 complex",
  # BRCA: membrane/cell periphery/extracellular, ETS-related TF
  "plasma membrane",
  "extracellular region",
  "extracellular space",
  "tube development",
  "2q37 copy number variation syndrome",
  # Pan-cancer: developmental terms, tube development, cell migration, E2F
  "multicellular organism development",
  "cell migration",
  "anatomical structure morphogenesis",
  # Operative: focal adhesion/EMT WikiPathways, TF motifs
  "Epithelial to mesenchymal transition in colorectal cancer",
  "T cell modulation and desmoplasia in pancreatic cancer",
  "Focal adhesion",
  "Signal Transduction",
  "Developmental Biology"
)

# Score-based selection per group: top 10 per group by score, keeping only non-generic
per_group_top <- dt_mapped[order(-score), head(.SD, 10), by = group_label]

# Build the priority-selected set: any row whose term_name is in priority list
priority_rows <- dt_mapped[term_name %in% priority_terms]

# Combine priority rows + per-group top rows
candidate_rows <- unique(rbind(priority_rows, per_group_top), by = c("group_label", "term_label"))

# For each unique term_label, pick the best (lowest p) per group
candidate_rows <- candidate_rows[order(p_value), .SD[1], by = .(group_label, term_label)]

# Now select which term_labels appear in the final figure:
# 1. Compute per-term maximum score across groups
term_max_score <- candidate_rows[, .(max_score = max(score)), by = term_label]
term_max_score <- term_max_score[order(-max_score)]

# 2. Flag priority terms
term_max_score[, is_priority := term_label %in% {
  priority_rows[, unique(term_label)]
}]

# 3. Select up to 30 terms: priority first, then highest-scoring non-priority
priority_selected    <- term_max_score[is_priority == TRUE]
nonpriority_selected <- term_max_score[is_priority == FALSE]

n_priority  <- min(nrow(priority_selected), 30L)
n_fill      <- 30L - n_priority
fill_terms  <- head(nonpriority_selected$term_label, n_fill)
final_terms <- c(priority_selected$term_label[seq_len(n_priority)], fill_terms)
final_terms <- unique(final_terms)

# 4. Keep selected terms only
plot_data <- candidate_rows[term_label %in% final_terms]

# Write excluded rows (before further filtering below)
excluded_rows <- dt_mapped[!term_label %in% final_terms]

# ── 8. Quality checks ─────────────────────────────────────────────────────────
dup_labels <- plot_data[, .N, by = term_label][N > 1 | duplicated(term_label)]
# (allow same label in different groups — duplicates within groups are the issue)
per_group_dup <- plot_data[, .N, by = .(group_label, term_label)][N > 1]
if (nrow(per_group_dup) > 0) {
  stop("Duplicate term_label within the same group_label detected:\n",
       paste(per_group_dup$term_label, collapse = "\n"))
}

n_selected_terms <- length(unique(plot_data$term_label))
cat("Number of selected unique terms before figure-level filtering:", n_selected_terms, "\n")

cat("Selected terms per group (before figure-level filtering):\n")
print(plot_data[, .N, by = group_label])

# ── 9. Figure-level term filtering (2e, 2f) ───────────────────────────────────
# 2e: Remove overly broad terms (case-insensitive substring match on term_label)
broad_patterns <- c(
  "Signal Transduction",
  "^signaling",
  "^cell communication",
  "^cellular process",
  "^response to stimulus",
  "^biological regulation",
  "^metabolic process"
)

is_broad <- function(label) {
  lc <- tolower(label)
  any(sapply(broad_patterns, function(pat) {
    if (startsWith(pat, "^")) {
      startsWith(lc, tolower(substring(pat, 2)))
    } else {
      grepl(tolower(pat), lc, fixed = TRUE)
    }
  }))
}

broad_labels <- unique(plot_data$term_label)[sapply(unique(plot_data$term_label), is_broad)]
if (length(broad_labels) > 0) {
  cat("Removing overly broad terms (2e):\n")
  cat(paste0("  - ", broad_labels, "\n"))
  # Move to excluded before removing from plot_data
  excluded_broad <- plot_data[term_label %in% broad_labels]
  excluded_rows  <- rbind(
    excluded_rows[, .(query_id, group_label, source, term_name, term_label, p_value, score)],
    excluded_broad[, .(query_id, group_label, source, term_name, term_label, p_value, score)],
    fill = TRUE
  )
  plot_data <- plot_data[!term_label %in% broad_labels]
}

# 2f: Deduplicate Focal adhesion — keep [WP], drop [KEGG]
kegg_fa_label <- "Focal adhesion [KEGG]"
wp_fa_label   <- "Focal adhesion [WP]"
if (kegg_fa_label %in% unique(plot_data$term_label) &&
    wp_fa_label   %in% unique(plot_data$term_label)) {
  cat("Both 'Focal adhesion [KEGG]' and 'Focal adhesion [WP]' present; keeping [WP], removing [KEGG].\n")
  excluded_kegg_fa <- plot_data[term_label == kegg_fa_label]
  excluded_rows <- rbind(
    excluded_rows,
    excluded_kegg_fa[, .(query_id, group_label, source, term_name, term_label, p_value, score)],
    fill = TRUE
  )
  plot_data <- plot_data[term_label != kegg_fa_label]
}

n_fig_terms <- length(unique(plot_data$term_label))
cat("Number of unique terms for figure after 2e/2f filtering:", n_fig_terms, "\n")

# Fill back from excluded candidate pool if below 22 terms
# (using the original candidate_rows pool, excluding already-selected/excluded terms)
if (n_fig_terms < 22) {
  current_labels <- unique(plot_data$term_label)
  excluded_labels <- unique(excluded_rows$term_label)
  all_used <- union(current_labels, excluded_labels)

  fillback_pool <- candidate_rows[!term_label %in% all_used]
  fillback_pool_max <- fillback_pool[, .(max_score = max(score)), by = term_label][order(-max_score)]

  n_needed <- 22L - n_fig_terms
  fill_labels <- head(fillback_pool_max$term_label, n_needed)

  if (length(fill_labels) > 0) {
    cat("Filling back", length(fill_labels), "terms to reach 22 minimum:\n")
    cat(paste0("  + ", fill_labels, "\n"))
    plot_data <- rbind(plot_data, candidate_rows[term_label %in% fill_labels], fill = TRUE)
  }
  n_fig_terms <- length(unique(plot_data$term_label))
}

cat("Final figure term count:", n_fig_terms, "\n")
if (n_fig_terms > 28) {
  warning("More than 28 rows in figure (", n_fig_terms, "). Figure may be crowded.")
}

# ── 10. Fix display labels (2g) ───────────────────────────────────────────────
fix_display_label <- function(lbl) {
  if (grepl("Retina photoreceptor", lbl, ignore.case = TRUE) ||
      grepl("Retina; photoreceptor", lbl, ignore.case = TRUE)) {
    return("Retina photoreceptor cells (>=low) [HPA]")
  }
  if (grepl("Interleukin-10", lbl, ignore.case = TRUE) &&
      grepl("signali", lbl, ignore.case = TRUE)) {
    return("Interleukin-10 signalling [REAC]")
  }
  if (grepl("Cytokine [Ss]ignal", lbl) && grepl("[Ii]mmune", lbl)) {
    return("Cytokine signalling in immune system [REAC]")
  }
  return(lbl)
}

plot_data[, term_label := vapply(term_label, fix_display_label, character(1))]

# ── 11. Build wide matrix for output TSV ──────────────────────────────────────
# (uses internal group_label, before x-axis renaming)
wide_mat <- dcast(
  plot_data,
  term_label ~ group_label,
  value.var = "score",
  fun.aggregate = max,
  fill = NA_real_
)

# ── 12. X-axis group labels with line breaks ──────────────────────────────────
# All 11 groups in prescribed order with \n line breaks for compactness
group_order_display <- c(
  "BRCA\nper-contrast",
  "NBL\nper-contrast",
  "RBL\nper-contrast",
  "BRCA\nrecurrent",
  "NBL\nrecurrent",
  "RBL\nrecurrent",
  "Recurrent core\ndown",
  "Recurrent core\nup",
  "Recurrent core\nall",
  "Final 257-gene\nfeature set",
  "Isolate-rescued\nsubset"
)

# Internal → display label mapping
# Internal group_label values (from map_query) → manuscript-facing x-axis labels
internal_to_display <- c(
  "BRCA per-contrast markers" = "BRCA\nper-contrast",
  "NBL per-contrast markers"  = "NBL\nper-contrast",
  "RBL per-contrast markers"  = "RBL\nper-contrast",
  "BRCA recurrent"            = "BRCA\nrecurrent",
  "NBL recurrent"             = "NBL\nrecurrent",
  "RBL recurrent"             = "RBL\nrecurrent",
  "Strict union down"         = "Recurrent core\ndown",
  "Strict union up"           = "Recurrent core\nup",
  "Strict union all"          = "Recurrent core\nall",
  "Operative feature set"     = "Final 257-gene\nfeature set",
  "Operative minus strict"    = "Isolate-rescued\nsubset"
)

plot_data[, group_display := internal_to_display[group_label]]

# ── 13. Prepare ggplot data ────────────────────────────────────────────────────
# Apply str_wrap AFTER display label fixes (2g)
plot_data[, term_label_wrap := str_wrap(term_label, width = 50)]

# ── 14. Row ordering into biological blocks (2h) ──────────────────────────────
# Compute per-term maximum score and which group holds it
term_group_max <- plot_data[, {
  idx <- which.max(score)
  .(max_score = max(score, na.rm = TRUE), best_group = group_label[idx])
}, by = term_label]

# Block assignment function
assign_block <- function(best_group) {
  brca_groups <- c("BRCA per-contrast markers", "BRCA recurrent")
  nbl_groups  <- c("NBL per-contrast markers", "NBL recurrent")
  rbl_groups  <- c("RBL per-contrast markers", "RBL recurrent")
  strict_groups <- c("Strict union down", "Strict union up", "Strict union all")
  operative_groups <- c("Operative feature set", "Operative minus strict")

  dplyr_case <- function(g) {
    if (g %in% brca_groups)     return(1L)
    if (g %in% nbl_groups)      return(2L)
    if (g %in% rbl_groups)      return(3L)
    if (g %in% strict_groups)   return(4L)
    if (g %in% operative_groups) return(5L)
    return(6L)  # cross-cutting
  }
  vapply(best_group, dplyr_case, integer(1))
}

term_group_max[, block := assign_block(best_group)]
# Sort: block ascending, then max_score descending within block
term_group_max <- term_group_max[order(block, -max_score)]

ordered_terms      <- term_group_max$term_label
ordered_terms_wrap <- str_wrap(ordered_terms, width = 50)

# Ensure term_label_wrap is wrapped consistently
plot_data[, term_label_wrap := str_wrap(term_label, width = 50)]

# Build ordered factor: bottom of heatmap = first in ordered_terms (low block / high score)
# ggplot y-axis reads bottom-to-top, so reverse for top-to-bottom display
plot_data[, term_label_wrap := factor(term_label_wrap, levels = rev(ordered_terms_wrap))]

# ── 15. Build complete grid with all 11 groups (2c) ───────────────────────────
# Use display labels for x-axis, including empty groups
plot_data[, group_display := factor(group_display, levels = group_order_display)]

all_combos <- CJ(
  term_label_wrap = factor(rev(ordered_terms_wrap), levels = rev(ordered_terms_wrap)),
  group_display   = factor(group_order_display, levels = group_order_display)
)
plot_full <- merge(
  all_combos,
  plot_data[, .(term_label_wrap, group_display, score)],
  by = c("term_label_wrap", "group_display"),
  all.x = TRUE
)

# ── 16. Build figure ──────────────────────────────────────────────────────────
n_rows <- length(ordered_terms_wrap)
n_cols <- length(group_order_display)  # always 11

p <- ggplot(plot_full, aes(x = group_display, y = term_label_wrap, fill = score)) +
  geom_tile(colour = "grey60", linewidth = 0.25) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    na.value  = "#F5F5F5",
    limits    = c(0, 30),
    name      = expression(atop(-log[10](p)*",", "capped at 30"))
  ) +
  scale_x_discrete(position = "bottom") +
  labs(
    title = "Selected enriched terms across disease-specific\nand pan-cancer marker-gene sets",
    x     = NULL,
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 11, face = "bold"),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    legend.key.width = unit(0.6, "cm")
  )

# ── 17. Save outputs ──────────────────────────────────────────────────────────
pdf_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.pdf")
png_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.png")
mat_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_matrix.tsv")
sel_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_selected_terms.tsv")
exc_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv")

ggsave(pdf_path, plot = p, width = 16, height = 11, device = "pdf")
ggsave(png_path, plot = p, width = 16, height = 11, dpi = 300, device = "png")
fwrite(wide_mat,   mat_path, sep = "\t")
fwrite(plot_data,  sel_path, sep = "\t")
fwrite(excluded_rows, exc_path, sep = "\t")

cat("PDF  :", pdf_path, "\n")
cat("PNG  :", png_path, "\n")
cat("Matrix TSV :", mat_path, "\n")
cat("Selected terms TSV:", sel_path, "\n")
cat("Excluded terms TSV:", exc_path, "\n")

# ── 18. Report groups with data vs empty ──────────────────────────────────────
cat("\n--- Group data status ---\n")
groups_with_data <- unique(plot_data$group_display)
for (g in group_order_display) {
  has_data <- g %in% as.character(groups_with_data)
  n_terms  <- if (has_data) sum(!is.na(plot_full[group_display == g, score])) else 0L
  cat(sprintf("  %-30s %s  (%d non-NA tiles)\n",
              gsub("\n", " ", g), if (has_data) "DATA" else "EMPTY", n_terms))
}
cat("Final row count:", n_rows, "\n")
cat("Done.\n")
