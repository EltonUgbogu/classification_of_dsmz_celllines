#!/usr/bin/env Rscript
# plot_enrichment_top_terms_heatmap.R
# Build an enrichment heatmap from the ranked marker-source-panel top-term summary,
# writing selected-term, excluded-term, matrix, and figure-provenance sidecars.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
pipeline_root <- Sys.getenv("PIPELINE_ROOT", unset = repo_root)

parse_cli <- function(args) {
  defaults <- list(
    input = "results/unsupervised/enrichment_summary_top_terms.tsv",
    outdir = "results/unsupervised/figures",
    feature_table = "results/unsupervised/pan_cancer/feature_space/pan_cancer_features.tsv"
  )
  if (length(args) == 0L) return(defaults)
  if (any(args %in% c("-h", "--help"))) {
    cat(
      "Usage: Rscript plot_enrichment_top_terms_heatmap.R --input enrichment_summary_top_terms.tsv --outdir figures [--feature-table pan_cancer_features.tsv]\n",
      "\nOptions:\n",
      "  --input PATH          Enrichment summary TSV with query_id, term_name, source, p_value.\n",
      "  --outdir PATH         Output directory for PDF/PNG/TSV sidecars.\n",
      "  --feature-table PATH  Current pan-cancer feature table used only to report feature count.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--input", "--outdir", "--feature-table")) stop("Unknown argument: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    arg_name <- gsub("-", "_", sub("^--", "", key))
    defaults[[arg_name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  defaults
}

feature_count_from_table <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  x <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(x) || !"gene_id" %in% names(x)) return(NA_integer_)
  length(unique(sub("\\.[0-9]+$", "", trimws(as.character(x$gene_id)))))
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
input_file <- cli$input
outdir <- cli$outdir
feature_table <- cli$feature_table
observed_feature_count <- feature_count_from_table(feature_table)
observed_feature_count_label <- ifelse(is.na(observed_feature_count), "unavailable", as.character(observed_feature_count))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file does not exist: ", input_file)

dt <- fread(input_file)
required_cols <- c("query_id", "term_name", "source", "p_value")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Input table is missing required columns: ", paste(missing_cols, collapse = ", "))
}
cat("Loaded", nrow(dt), "rows from input table.\n")
cat("Observed current feature count:", observed_feature_count_label, "\n")

map_query <- function(q, category, cohort, marker_source_class, evidence_class, direction) {
  ql <- tolower(trimws(as.character(q)))
  category <- tolower(trimws(as.character(category)))
  cohort <- toupper(trimws(as.character(cohort)))
  marker_source_class <- tolower(trimws(as.character(marker_source_class)))
  evidence_class <- tolower(trimws(as.character(evidence_class)))
  direction <- toupper(trimws(as.character(direction)))
  ql[is.na(ql)] <- ""
  category[is.na(category)] <- ""
  cohort[is.na(cohort)] <- ""
  marker_source_class[is.na(marker_source_class)] <- ""
  evidence_class[is.na(evidence_class)] <- ""
  direction[is.na(direction)] <- ""

  if (category == "final_panel" || grepl("__final_panel_all$", ql)) {
    return("Final pan-cancer feature set")
  }
  if (category == "cohort_owner" && cohort %in% c("BRCA", "NBL", "RBL")) {
    return(paste(cohort, "cohort-derived markers"))
  }
  if (category == "per_contrast") {
    return("Per-contrast marker sets")
  }
  if (category %in% c("evidence_class", "cohort_evidence_class")) {
    if (grepl("isolate_source", evidence_class)) return("Isolate-source evidence")
    if (grepl("source_recurrent", evidence_class)) return("Recurrent marker-source evidence")
    if (grepl("anchor_singleton", evidence_class)) return("Anchor singleton evidence")
    if (grepl("nonrecurrent", evidence_class)) return("Ranked nonrecurrent markers")
  }
  if (category %in% c("marker_source_class", "cohort_marker_source_class",
                      "marker_source_class_direction", "cohort_marker_source_class_direction")) {
    if (marker_source_class == "isolate") return("Isolate-source evidence")
    if (marker_source_class == "anchor") return("Anchor marker-source evidence")
  }
  if (category %in% c("direction", "cohort_direction")) {
    if (direction == "UP") return("UP marker sets")
    if (direction == "DOWN") return("DOWN marker sets")
    if (direction == "MIXED") return("Mixed-direction marker sets")
    return("Direction-stratified marker sets")
  }
  return(NA_character_)
}

for (col in c("category", "cohort", "marker_source_class", "evidence_class", "direction")) {
  if (!col %in% names(dt)) dt[, (col) := NA_character_]
}
dt[, group_label := mapply(
  map_query,
  dt[["query_id"]], dt[["category"]], dt[["cohort"]],
  dt[["marker_source_class"]], dt[["evidence_class"]], dt[["direction"]],
  USE.NAMES = FALSE
)]
dt_mapped <- dt[!is.na(group_label)]
cat("Rows after group mapping:", nrow(dt_mapped), "\n")
if (nrow(dt_mapped) == 0L) stop("No enrichment rows mapped to current marker-source-panel heatmap groups")

generic_terms <- c(
  "biological_process", "cellular_component", "molecular_function", "binding",
  "cellular process", "response to stimulus", "biological regulation", "metabolic process",
  "multicellular organismal process", "cellular anatomical structure", "cell periphery",
  "HPA root", "KEGG root term", "REACTOME root term", "Transfac", "Transfac [TF]",
  "WIKIPATHWAYS", "WIKIPATHWAYS [WP]"
)
n_before_generic <- nrow(dt_mapped)
dt_mapped <- dt_mapped[!term_name %in% generic_terms]
cat("Removed", n_before_generic - nrow(dt_mapped), "generic/root term rows.\n")
if (nrow(dt_mapped) == 0L) stop("No enrichment rows remain after generic/root term filtering")

dt_mapped[, p_value_numeric := suppressWarnings(as.numeric(p_value))]
dt_mapped <- dt_mapped[!is.na(p_value_numeric) & p_value_numeric > 0]
dt_mapped[, score := pmin(-log10(pmax(p_value_numeric, 1e-300)), 30)]

strip_match_class <- function(x) gsub(";?\\s*match class:\\s*[0-9]+", "", x, perl = TRUE)
dt_mapped[source == "TF", tf_base := strip_match_class(term_name)]
dt_tf <- dt_mapped[source == "TF"][order(p_value_numeric)][, .SD[1], by = .(query_id, tf_base)]
dt_notf <- dt_mapped[source != "TF"]
dt_mapped <- rbind(dt_notf, dt_tf, fill = TRUE)

simplify_tf_label <- function(x) {
  fac <- str_match(x, "Factor:\\s*([^;]+)")[, 2]
  mot <- str_match(x, "motif:\\s*([^;]+)")[, 2]
  ifelse(!is.na(fac) & !is.na(mot), paste0(trimws(fac), " motif (", trimws(mot), ") [TF]"), paste0(x, " [TF]"))
}
source_bracket <- c(
  "GO:BP" = "[GO:BP]", "GO:CC" = "[GO:CC]", "GO:MF" = "[GO:MF]",
  "REAC" = "[REAC]", "KEGG" = "[KEGG]", "WP" = "[WP]", "TF" = "[TF]",
  "CORUM" = "[CORUM]", "HPA" = "[HPA]"
)
dt_mapped[, term_label := {
  br <- source_bracket[source]
  br[is.na(br)] <- paste0("[", source[is.na(br)], "]")
  ifelse(source == "TF", simplify_tf_label(term_name), paste0(term_name, " ", br))
}]

# Keep the strongest entry per group/term, then select a compact cross-group panel.
candidate_rows <- dt_mapped[order(p_value_numeric), .SD[1], by = .(group_label, term_label)]
per_group_top <- candidate_rows[order(-score), head(.SD, 10), by = group_label]
term_max <- per_group_top[, .(max_score = max(score), n_groups = uniqueN(group_label)), by = term_label]
term_max <- term_max[order(-n_groups, -max_score)]
selected_terms <- head(term_max$term_label, 30)
plot_data <- per_group_top[term_label %in% selected_terms]
if (nrow(plot_data) == 0L) stop("No selected enrichment terms available for plotting")
excluded_rows <- dt_mapped[!term_label %in% selected_terms]

fix_display_label <- function(lbl) {
  if (grepl("Retina photoreceptor", lbl, ignore.case = TRUE) || grepl("Retina; photoreceptor", lbl, ignore.case = TRUE)) {
    return("Retina photoreceptor cells (>=low) [HPA]")
  }
  if (grepl("Interleukin-10", lbl, ignore.case = TRUE) && grepl("signali", lbl, ignore.case = TRUE)) {
    return("Interleukin-10 signalling [REAC]")
  }
  if (grepl("Cytokine [Ss]ignal", lbl) && grepl("[Ii]mmune", lbl)) {
    return("Cytokine signalling in immune system [REAC]")
  }
  lbl
}
plot_data[, term_label := vapply(term_label, fix_display_label, character(1))]

wide_mat <- dcast(plot_data, term_label ~ group_label, value.var = "score", fun.aggregate = max, fill = NA_real_)

group_order_display <- c(
  "Final pan-cancer\nfeature set",
  "BRCA cohort-derived\nmarkers",
  "NBL cohort-derived\nmarkers",
  "RBL cohort-derived\nmarkers",
  "Per-contrast\nmarker sets",
  "Recurrent marker-source\nevidence",
  "Isolate-source\nevidence",
  "Anchor marker-source\nevidence",
  "Anchor singleton\nevidence",
  "Ranked nonrecurrent\nmarkers",
  "UP marker\nsets",
  "DOWN marker\nsets",
  "Mixed-direction\nmarker sets"
)
internal_to_display <- c(
  "Final pan-cancer feature set" = "Final pan-cancer\nfeature set",
  "BRCA cohort-derived markers" = "BRCA cohort-derived\nmarkers",
  "NBL cohort-derived markers" = "NBL cohort-derived\nmarkers",
  "RBL cohort-derived markers" = "RBL cohort-derived\nmarkers",
  "Per-contrast marker sets" = "Per-contrast\nmarker sets",
  "Recurrent marker-source evidence" = "Recurrent marker-source\nevidence",
  "Isolate-source evidence" = "Isolate-source\nevidence",
  "Anchor marker-source evidence" = "Anchor marker-source\nevidence",
  "Anchor singleton evidence" = "Anchor singleton\nevidence",
  "Ranked nonrecurrent markers" = "Ranked nonrecurrent\nmarkers",
  "UP marker sets" = "UP marker\nsets",
  "DOWN marker sets" = "DOWN marker\nsets",
  "Mixed-direction marker sets" = "Mixed-direction\nmarker sets"
)
plot_data[, group_display := internal_to_display[group_label]]
plot_data <- plot_data[!is.na(group_display)]

term_group_max <- plot_data[, {
  idx <- which.max(score)
  .(max_score = max(score, na.rm = TRUE), best_group = group_label[idx])
}, by = term_label]
assign_block <- function(best_group) {
  vapply(best_group, function(g) {
    if (g %in% c("Final pan-cancer feature set")) return(1L)
    if (g %in% c("BRCA cohort-derived markers", "NBL cohort-derived markers", "RBL cohort-derived markers")) return(2L)
    if (g %in% c("Per-contrast marker sets")) return(3L)
    if (g %in% c("Recurrent marker-source evidence", "Isolate-source evidence",
                 "Anchor marker-source evidence", "Anchor singleton evidence",
                 "Ranked nonrecurrent markers")) return(4L)
    if (g %in% c("UP marker sets", "DOWN marker sets", "Mixed-direction marker sets")) return(5L)
    6L
  }, integer(1))
}
term_group_max[, block := assign_block(best_group)]
term_group_max <- term_group_max[order(block, -max_score)]
ordered_terms <- term_group_max$term_label
ordered_terms_wrap <- str_wrap(ordered_terms, width = 50)
plot_data[, term_label_wrap := factor(str_wrap(term_label, width = 50), levels = rev(ordered_terms_wrap))]
plot_data[, group_display := factor(group_display, levels = group_order_display)]

all_combos <- CJ(
  term_label_wrap = factor(rev(ordered_terms_wrap), levels = rev(ordered_terms_wrap)),
  group_display = factor(group_order_display, levels = group_order_display)
)
plot_full <- merge(all_combos, plot_data[, .(term_label_wrap, group_display, score)], by = c("term_label_wrap", "group_display"), all.x = TRUE)

n_rows <- length(ordered_terms_wrap)
plot_height <- max(7, 3 + 0.34 * n_rows)
p <- ggplot(plot_full, aes(x = group_display, y = term_label_wrap, fill = score)) +
  geom_tile(colour = "grey60", linewidth = 0.25) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "#F5F5F5", limits = c(0, 30), name = "-log10(p),\ncapped at 30") +
  scale_x_discrete(position = "bottom") +
  labs(
    title = "Selected enriched terms across current marker-source-panel query sets",
    subtitle = paste("Observed current feature count:", observed_feature_count_label),
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.6, "cm")
  )

pdf_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.pdf")
png_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.png")
mat_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_matrix.tsv")
sel_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_selected_terms.tsv")
exc_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv")
prov_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_provenance.tsv")

ggsave(pdf_path, plot = p, width = 16, height = plot_height, device = "pdf")
ggsave(png_path, plot = p, width = 16, height = plot_height, dpi = 300, device = "png")
fwrite(wide_mat, mat_path, sep = "\t")
fwrite(plot_data, sel_path, sep = "\t")
fwrite(excluded_rows, exc_path, sep = "\t")

# Record enrichment-heatmap provenance. Endpoint and g:Profiler run provenance
# remain in the upstream g:Profiler output directory.
prov <- data.table(
  figure_name = "Fig_enrichment_top_terms_heatmap.pdf",
  script = "scripts/plot_enrichment_top_terms_heatmap.R",
  command = paste(c("Rscript", "scripts/plot_enrichment_top_terms_heatmap.R", "--input", input_file, "--outdir", outdir, "--feature-table", feature_table), collapse = " "),
  git_commit = Sys.getenv("GIT_COMMIT", unset = "unavailable_not_git_worktree"),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_files = normalizePath(c(input_file, feature_table), mustWork = FALSE),
  output_files = paste(normalizePath(c(pdf_path, png_path, mat_path, sel_path, exc_path), mustWork = FALSE), collapse = ";"),
  upstream_tables = normalizePath(input_file, mustWork = FALSE),
  observed_feature_count = observed_feature_count_label,
  key_parameters = "group_mapping=current_marker_source_panel;generic_terms_removed=TRUE;top_terms_heatmap=TRUE;feature_count_dynamic=TRUE",
  software_versions = paste0("R=", getRversion(), ";data.table=", as.character(packageVersion("data.table")), ";ggplot2=", as.character(packageVersion("ggplot2"))),
  figure_type = "enrichment",
  source_pipeline_root = pipeline_root,
  copied_to_figure_export_path = "",
  legacy_source_path = "",
  notes = paste("Neutral labels; observed current feature count", observed_feature_count_label)
)
fwrite(prov, prov_path, sep = "\t")

cat("PDF:", pdf_path, "\n")
cat("PNG:", png_path, "\n")
cat("Matrix TSV:", mat_path, "\n")
cat("Selected terms TSV:", sel_path, "\n")
cat("Excluded terms TSV:", exc_path, "\n")
cat("Provenance TSV:", prov_path, "\n")
cat("Done.\n")
