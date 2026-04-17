#!/usr/bin/env Rscript

# =============================================================================
# plot_pan_cancer_graph.R
# =============================================================================
# Pan-cancer (pan-disease) cell line → cell line consensus graph plots
#
# HARD CONSTRAINTS (project-specific):
#   1) Working on 258 genes (RDS must contain an expr matrix with 258 genes)
#   2) Only concerned with cell lines
#   3) In the .rds, cell line IDs must start with "NG-"
#
# Inputs:
#  --edges       graph/pan_cancer_graph_edges.tsv
#  --components  graph/pan_cancer_components.tsv
#  --meta        pan_cancer_expr.rds   (for lineage labels; cell lines only)
#
# Outputs (PDF):
#  Fig_pan_cancer_graph.pdf
#  Fig_pan_cancer_graph_components.pdf
#  Fig_pan_cancer_component_size.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("--edges"), type="character", default=NULL,
              help="Path to pan_cancer_graph_edges.tsv"),
  make_option(c("--components"), type="character", default=NULL,
              help="Path to pan_cancer_components.tsv"),
  make_option(c("--meta"), type="character", default=NULL,
              help="Path to pan_cancer_expr.rds (contains $meta with sample_id, lineage, type; and optionally expr)"),
  make_option(c("--outdir"), type="character", default=".",
              help="Output directory for plots"),
  make_option(c("--layout"), type="character", default="fr",
              help="Layout: fr (Fruchterman-Reingold), kk (Kamada-Kawai), lgl"),
  make_option(c("--seed"), type="integer", default=1,
              help="Random seed for layout reproducibility"),
  make_option(c("--label-top-hubs"), type="integer", default=0,
              help="Label top N hubs by degree (0 = none)"),
  make_option(c("--require_ng_prefix"), type="logical", default=TRUE,
              help="If TRUE, enforce node IDs start with NG- (default TRUE)"),
  make_option(c("--require_258_genes"), type="logical", default=TRUE,
              help="If TRUE, enforce expr matrix contains exactly 258 genes (default TRUE)")
)

opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opt$edges), file.exists(opt$edges))
stopifnot(!is.null(opt$components), file.exists(opt$components))
stopifnot(!is.null(opt$meta), file.exists(opt$meta))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

is_ng <- function(x) {
  !is.na(x) & grepl("^NG-", x)
}

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
edges <- fread(opt$edges)
comps <- fread(opt$components)

# Accept flexible column names for edges/components
edge_cols <- names(edges)
from_col <- intersect(edge_cols, c("from","source","cell_line_1","cell1","u","A"))[1]
to_col   <- intersect(edge_cols, c("to","target","cell_line_2","cell2","v","B"))[1]
if (is.na(from_col) || is.na(to_col)) {
  stop("Could not infer edge columns in edges TSV. Expected columns like from/to or cell_line_1/cell_line_2.")
}
setnames(edges, c(from_col, to_col), c("from", "to"))

comp_cols <- names(comps)
sample_col <- intersect(comp_cols, c("sample","sample_id","cell_line","node"))[1]
comp_col   <- intersect(comp_cols, c("component","comp","community"))[1]
if (is.na(sample_col) || is.na(comp_col)) {
  stop("Could not infer component columns in components TSV. Expected sample/cell_line and component.")
}
setnames(comps, c(sample_col, comp_col), c("sample_id", "component"))

# meta from RDS
cell_dat <- readRDS(opt$meta)
if (!("meta" %in% names(cell_dat))) {
  stop("RDS object must contain $meta.")
}
meta <- as.data.table(cell_dat$meta)

# -----------------------------------------------------------------------------
# Enforce: cell lines only + NG- prefix
# -----------------------------------------------------------------------------
req_meta <- c("sample_id","lineage")
missing <- setdiff(req_meta, names(meta))
if (length(missing) > 0) stop("meta missing required columns: ", paste(missing, collapse=", "))

# If a type column exists, keep only cell lines
if ("type" %in% names(meta)) {
  meta <- meta[type == "cell_line"]
}

# Enforce NG- prefix in meta
if (isTRUE(opt$require_ng_prefix)) {
  meta <- meta[is_ng(sample_id)]
}

meta <- unique(meta[, .(sample_id, lineage)])

# -----------------------------------------------------------------------------
# Enforce: 258 genes in the RDS expr (sanity check)
# -----------------------------------------------------------------------------
if (isTRUE(opt$require_258_genes)) {
  expr_obj_name <- intersect(names(cell_dat), c("expr", "expression", "mat", "matrix", "counts", "tpm"))[1]
  if (is.na(expr_obj_name)) {
    warning("No obvious expr matrix found in RDS (expected one of: expr/expression/mat/matrix/counts/tpm). Skipping 258-gene validation.")
  } else {
    expr <- cell_dat[[expr_obj_name]]
    if (!is.matrix(expr) && !is.data.frame(expr)) {
      stop("RDS $", expr_obj_name, " exists but is not matrix/data.frame; cannot validate gene count.")
    }
    nr <- nrow(expr); nc <- ncol(expr)

    # genes might be rows or columns depending on how it was saved
    if (!(nr == 258 || nc == 258)) {
      stop(
        "Expected expr matrix to have 258 genes in either rows or columns, but got: nrow=",
        nr, ", ncol=", nc, "."
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Filter edges/components to NG- only (and to those present in meta when possible)
# -----------------------------------------------------------------------------
if (isTRUE(opt$require_ng_prefix)) {
  edges <- edges[is_ng(from) & is_ng(to)]
  comps <- comps[is_ng(sample_id)]
}

# Optionally restrict to nodes we actually have in meta (prevents stray nodes)
keep_ids <- meta$sample_id
edges <- edges[from %in% keep_ids & to %in% keep_ids]
comps <- comps[sample_id %in% keep_ids]

if (nrow(edges) == 0) stop("After filtering to NG- cell lines and meta IDs, edges is empty. Check inputs.")
if (nrow(comps) == 0) warning("After filtering, components table is empty; component coloring will be NA.")

# -----------------------------------------------------------------------------
# Build graph
# -----------------------------------------------------------------------------
g <- graph_from_data_frame(edges[, .(from, to)], directed=FALSE)

# Attach vertex annotations
vtab <- data.table(sample_id = V(g)$name)
vtab <- merge(vtab, meta, by="sample_id", all.x=TRUE)
vtab <- merge(vtab, comps, by="sample_id", all.x=TRUE)

V(g)$lineage <- vtab$lineage
V(g)$component <- vtab$component

# Degree for optional hub labelling
deg <- degree(g)
V(g)$degree <- deg

# -----------------------------------------------------------------------------
# Layout
# -----------------------------------------------------------------------------
set.seed(opt$seed)
lay <- switch(
  opt$layout,
  "kk"  = layout_with_kk(g),
  "lgl" = layout_with_lgl(g),
  "fr"  = layout_with_fr(g),
  layout_with_fr(g)
)

# -----------------------------------------------------------------------------
# Plot helper (base plotting to avoid extra deps)
# -----------------------------------------------------------------------------
pdf_graph <- function(path, color_by=c("lineage","component"), main=NULL, label_top=0,
                      show_legend=TRUE, show_component_labels=FALSE) {
  color_by <- match.arg(color_by)

  vals <- if (color_by == "lineage") V(g)$lineage else as.character(V(g)$component)
  vals[is.na(vals)] <- "NA"

  f <- factor(vals)
  pal <- grDevices::hcl.colors(nlevels(f), palette = "Spectral")
  col_map <- setNames(pal, levels(f))
  vcol <- unname(col_map[as.character(f)])

  vs <- pmax(3, pmin(10, 3 + sqrt(V(g)$degree)))

  lab <- rep("", vcount(g))
  if (label_top > 0) {
    top_idx <- order(V(g)$degree, decreasing=TRUE)[seq_len(min(label_top, vcount(g)))]
    lab[top_idx] <- V(g)$name[top_idx]
  }

  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  n_comp <- length(unique(V(g)$component[!is.na(V(g)$component)]))
  subtitle <- sprintf("Nodes = %d NG- cell lines; Edges = %d; Components = %d; Layout = %s",
                      n_nodes, n_edges, n_comp, toupper(opt$layout))

  grDevices::pdf(path, width=11, height=8.5)
  par(mar=c(0,0,3,0))
  plot(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vs,
    vertex.label = lab,
    vertex.label.cex = 0.6,
    vertex.label.color = "black",
    vertex.frame.color = "white",
    edge.color = grDevices::adjustcolor("black", alpha.f = 0.15),
    edge.width = 0.6,
    main = ifelse(is.null(main),
                  paste0("Pan-cancer cell line graph (colored by ", color_by, ")"),
                  main)
  )

  mtext(subtitle, side=3, line=0.5, cex=0.9, col="gray40")

  if (show_component_labels && color_by == "component") {
    comp_ids <- unique(V(g)$component[!is.na(V(g)$component)])
    for (comp_id in comp_ids) {
      comp_nodes <- which(V(g)$component == comp_id)
      if (length(comp_nodes) > 0) {
        comp_coords <- lay[comp_nodes, , drop=FALSE]
        centroid <- colMeans(comp_coords)
        text(centroid[1], centroid[2], labels = paste("Community", comp_id), cex = 1.2, font = 2,
             col = "black", bg = "white", border = "black", lwd = 1.5)
      }
    }
  }

  if (show_legend) {
    if (color_by == "component" && nlevels(f) > 10) {
      # skip legend
    } else {
      legend(
        "topleft",
        legend = if (color_by == "component") paste("Community", names(col_map)) else names(col_map),
        col = unname(col_map),
        pch = 16,
        pt.cex = 1.2,
        cex = 0.8,
        bty = "n",
        title = color_by
      )
    }
  }

  dev.off()
}

# -----------------------------------------------------------------------------
# Write plots
# -----------------------------------------------------------------------------
n_components <- length(unique(V(g)$component[!is.na(V(g)$component)]))

pdf_graph(
  file.path(opt$outdir, "Fig_pan_cancer_graph.pdf"),
  color_by = "lineage",
  main = "Pan-cancer NG- cell line similarity graph (colored by lineage)",
  label_top = opt$`label-top-hubs`,
  show_legend = TRUE,
  show_component_labels = FALSE
)

pdf_graph(
  file.path(opt$outdir, "Fig_pan_cancer_graph_components.pdf"),
  color_by = "component",
  main = "Pan-cancer NG- cell line similarity graph (colored by component)",
  label_top = 0,
  show_legend = (n_components <= 10),
  show_component_labels = (n_components > 10)
)

comp_sizes <- as.data.table(components(g)$membership)
setnames(comp_sizes, c("V1"), c("component"))
comp_sizes[, component := as.character(component)]
comp_sizes <- comp_sizes[, .N, by=component][order(-N)]

p_comp <- ggplot(comp_sizes, aes(x=reorder(component, N), y=N)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Pan-cancer graph component sizes (NG- cell lines only)",
    x = "Component",
    y = "Number of cell lines"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(
  filename = file.path(opt$outdir, "Fig_pan_cancer_component_size.pdf"),
  plot = p_comp,
  device = "pdf",
  width = 8, height = 6, units = "in"
)

cat("Done.\nSaved:\n",
    " - Fig_pan_cancer_graph.pdf\n",
    " - Fig_pan_cancer_graph_components.pdf\n",
    " - Fig_pan_cancer_component_size.pdf\n",
    "in: ", opt$outdir, "\n", sep="")