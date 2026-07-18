#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--maximum-terms", type = "integer"),
  make_option("--term-wrap-width", type = "integer"),
  make_option("--figure-width", type = "numeric"),
  make_option("--minimum-figure-height", type = "numeric"),
  make_option("--row-height-increment", type = "numeric")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (nm in c("input", "outdir", "maximum-terms", "term-wrap-width", "figure-width",
             "minimum-figure-height", "row-height-increment")) {
  if (is.null(opt[[nm]]) || (is.character(opt[[nm]]) && !nzchar(opt[[nm]]))) {
    stop("--", nm, " is required", call. = FALSE)
  }
}

dt <- fread(normalizePath(opt$input, mustWork = TRUE), sep = "\t")
required <- c("cancer_type", "marker_evidence_stratum", "source", "term_id", "term_name",
              "significant_contrast_count", "runnable_contrast_count", "term_support_fraction")
missing <- setdiff(required, names(dt))
if (length(missing)) stop("Contrast support table missing columns: ", paste(missing, collapse = ","), call. = FALSE)
if (!nrow(dt)) stop("No contrast-support terms available for plotting", call. = FALSE)

dt[, stratum_id := paste(cancer_type, marker_evidence_stratum, sep = " | ")]
term_rank <- dt[
  ,
  .(
    maximum_term_support_fraction = max(term_support_fraction, na.rm = TRUE),
    total_significant_contrast_count = sum(significant_contrast_count, na.rm = TRUE),
    representative_term_name = term_name[which.max(term_support_fraction)]
  ),
  by = .(source, term_id)
]
setorder(term_rank, -maximum_term_support_fraction, -total_significant_contrast_count, source, term_id)
selected_terms <- term_rank[seq_len(min(.N, opt$`maximum-terms`))]
selected_terms[, term_key := paste(source, term_id, sep = "::")]
dt[, term_key := paste(source, term_id, sep = "::")]
plot_data <- merge(
  CJ(stratum_id = sort(unique(dt$stratum_id)), term_key = selected_terms$term_key, unique = TRUE),
  dt,
  by = c("stratum_id", "term_key"),
  all.x = TRUE
)
plot_data <- merge(
  plot_data,
  selected_terms[, .(term_key, representative_term_name, term_order = seq_len(.N))],
  by = "term_key",
  all.x = TRUE
)
plot_data[is.na(term_support_fraction), term_support_fraction := 0]
plot_data[, term_label := vapply(strwrap(representative_term_name, width = opt$`term-wrap-width`, simplify = FALSE), paste, character(1), collapse = "\n")]
plot_data[, term_label := factor(term_label, levels = rev(unique(plot_data[order(term_order)]$term_label)))]
plot_data[, stratum_id := factor(stratum_id, levels = sort(unique(stratum_id)))]

p <- ggplot(plot_data, aes(x = stratum_id, y = term_label, fill = term_support_fraction)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#eef1f4", high = "#31688e", limits = c(0, 1), name = "Contrast term-support fraction") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
height <- max(opt$`minimum-figure-height`, length(unique(plot_data$term_key)) * opt$`row-height-increment` + 2)
pdf_path <- file.path(outdir, "Fig_enrichment_contrast_term_support_heatmap.pdf")
png_path <- file.path(outdir, "Fig_enrichment_contrast_term_support_heatmap.png")
matrix_path <- file.path(outdir, "Fig_enrichment_contrast_term_support_heatmap_matrix.tsv")
selected_path <- file.path(outdir, "Fig_enrichment_contrast_term_support_heatmap_selected_terms.tsv")
provenance_path <- file.path(outdir, "Fig_enrichment_contrast_term_support_heatmap_provenance.tsv")
ggsave(pdf_path, p, width = opt$`figure-width`, height = height, units = "in")
ggsave(png_path, p, width = opt$`figure-width`, height = height, units = "in", dpi = 300)
fwrite(plot_data, matrix_path, sep = "\t", na = "NA")
fwrite(selected_terms, selected_path, sep = "\t", na = "NA")
fwrite(
  data.table(
    figure_name = "Fig_enrichment_contrast_term_support_heatmap.pdf",
    figure_type = "contrast_term_support_heatmap",
    input_files = normalizePath(opt$input, mustWork = FALSE),
    output_files = paste(normalizePath(c(pdf_path, png_path, matrix_path, selected_path), mustWork = FALSE), collapse = ";"),
    statistic = "term_support_fraction",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ),
  provenance_path,
  sep = "\t"
)
message("[Contrast interpretation] Contrast term-support heatmap written.")
