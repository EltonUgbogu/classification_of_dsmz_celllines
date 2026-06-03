#!/usr/bin/env Rscript
# =============================================================================
# plot_pan_cancer_lineage_community_figure.R
# =============================================================================
#
# Reproducible wrapper for the thesis two-panel pan-cancer network figure:
#
#   figures/Fig_pan_cancer_cell_line_similarity_network_lineage_community.pdf
#   figures/Fig_pan_cancer_cell_line_similarity_network_lineage_community.png
#
# Panel A – nodes coloured by cancer lineage (Okabe-Ito palette).
# Panel B – same fixed layout, nodes coloured by Louvain community (Dark2
#           palette) with convex-hull overlays and C1-C5 centroid labels.
# Right column – shared legends for lineage, community, and degree centrality.
#
# Node size is proportional to degree centrality (scaled linearly across the
# observed degree range).
#
# The script uses the FIXED layout coordinates from --layout so the geometry
# is deterministic. No stochastic layout algorithm is re-run.
#
# THESIS-REPORTED VALIDATION TARGETS
#   Nodes              = 167
#   Edges              = 2,130
#   Louvain communities = 5
#   Modularity Q       = 0.629  (read from --modularity if supplied)
#   Lineage assort. r  = 0.801  (read from --assortativity if supplied)
#
# REQUIRED INPUTS
#   --edges          path to pan_cancer_graph_edges.tsv
#                    (columns: from, to, weight)
#   --communities    path to pan_cancer_communities.tsv
#                    (columns: sample, component OR community, lineage)
#   --layout         path to recovered_svg_layout.tsv
#                    (columns: sample, x, y)
#   --out-dir        output directory (default: figures)
#
# OPTIONAL INPUTS
#   --modularity     validation_modularity.tsv for Q check (optional)
#   --assortativity  validation_assortativity.tsv for r check (optional)
#
# UPSTREAM PIPELINE ORDER (if required inputs are missing)
#   1. build_pan_cancer_expression_matrix.R   → pan_cancer_expr.rds
#   2. compute_pan_cancer_correlation.R       → pan_cancer_cor.rds
#   3. build_pan_cancer_graph.R               → pan_cancer_graph_edges.tsv
#   4. compute_pan_cancer_communities.R       → pan_cancer_communities.tsv
#                                               + validation_modularity.tsv
#                                               + validation_assortativity.tsv
#   5. (recover or compute layout coords)     → recovered_svg_layout.tsv
#   6. plot_pan_cancer_lineage_community_figure.R (this script)
#
# SOURCE DATA (archived on Google Drive, dsmz/results/pan_cancer 11.02.47):
#   graph/pan_cancer_graph_edges.tsv
#   graph/pan_cancer_communities.tsv
#   graph/community_validation/recovered_svg_layout.tsv
#   graph/community_validation/validation_modularity.tsv
#
# =============================================================================

# Load only minimal packages first so the missing-file check can run cleanly
# without requiring the full visualisation stack.
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# =============================================================================
# COMMAND-LINE ARGUMENTS
# =============================================================================

option_list <- list(
  make_option("--edges",        type = "character", default = NULL,
              help = "[REQUIRED] Edge list TSV (columns: from, to, weight)"),
  make_option("--communities",  type = "character", default = NULL,
              help = "[REQUIRED] Community/lineage TSV (columns: sample, component/community, lineage)"),
  make_option("--layout",       type = "character", default = NULL,
              help = "[REQUIRED] Fixed layout TSV (columns: sample, x, y)"),
  make_option("--out-dir",      type = "character", default = "figures",
              help = "Output directory for PDF and PNG (default: figures)"),
  make_option("--modularity",   type = "character", default = NULL,
              help = "Optional validation_modularity.tsv for Q verification"),
  make_option("--assortativity",type = "character", default = NULL,
              help = "Optional validation_assortativity.tsv for r verification")
)

opt <- parse_args(OptionParser(option_list = option_list))

# =============================================================================
# INPUT VALIDATION
# =============================================================================

cat("=================================================\n")
cat("PAN-CANCER LINEAGE / COMMUNITY NETWORK FIGURE\n")
cat("=================================================\n\n")

MISSING_CRITICAL <- FALSE

check_file <- function(path, label, note = NULL) {
  if (is.null(path) || !nzchar(path)) {
    cat(sprintf("[MISSING] %s -- argument not supplied\n", label))
    if (!is.null(note)) cat(sprintf("          %s\n", note))
    MISSING_CRITICAL <<- TRUE
    return(FALSE)
  }
  if (!file.exists(path)) {
    cat(sprintf("[MISSING] %s\n          %s\n", label, path))
    if (!is.null(note)) cat(sprintf("          %s\n", note))
    MISSING_CRITICAL <<- TRUE
    return(FALSE)
  }
  cat(sprintf("[OK]      %s\n", path))
  return(TRUE)
}

invisible(check_file(opt$edges,       "Edge list (pan_cancer_graph_edges.tsv)",
           "Run: build_pan_cancer_graph.R"))
invisible(check_file(opt$communities, "Community/lineage table (pan_cancer_communities.tsv)",
           "Run: compute_pan_cancer_communities.R"))
invisible(check_file(opt$layout,      "Fixed layout coordinates (recovered_svg_layout.tsv)",
           "Use saved layout or run layout_with_fr() and save coordinates"))

cat("\n")

if (MISSING_CRITICAL) {
  cat("=================================================\n")
  cat("REQUIRED INPUTS MISSING -- cannot generate figure.\n\n")
  cat("Typical command once inputs exist:\n\n")
  cat("  Rscript scripts/plot_pan_cancer_lineage_community_figure.R \\\n")
  cat("    --edges       results/unsupervised/pan_cancer/graph/pan_cancer_graph_edges.tsv \\\n")
  cat("    --communities results/unsupervised/pan_cancer/graph/pan_cancer_communities.tsv \\\n")
  cat("    --layout      results/unsupervised/pan_cancer/graph/community_validation/recovered_svg_layout.tsv \\\n")
  cat("    --modularity  results/unsupervised/pan_cancer/graph/community_validation/validation_modularity.tsv \\\n")
  cat("    --out-dir     figures\n")
  cat("=================================================\n")
  quit(save = "no", status = 1)
}

# All required inputs present — now load visualisation packages.
suppressPackageStartupMessages({
  library(igraph)
  library(RColorBrewer)
})

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("[1] Loading data...\n")

edges      <- fread(opt$edges,       quote = "")
comm_tbl   <- fread(opt$communities, quote = "")
layout_tbl <- fread(opt$layout,      quote = "")

# Validate required columns
req_edge <- c("from", "to", "weight")
miss_edge <- setdiff(req_edge, names(edges))
if (length(miss_edge) > 0) stop("Edge file missing columns: ", paste(miss_edge, collapse = ", "))

# communities: accept both 'component' (historical) and 'community'
if ("sample_id" %in% names(comm_tbl) && !"sample" %in% names(comm_tbl)) {
  setnames(comm_tbl, "sample_id", "sample")
}
if ("component" %in% names(comm_tbl) && !"community" %in% names(comm_tbl)) {
  setnames(comm_tbl, "component", "community")
}
req_comm <- c("sample", "community", "lineage")
miss_comm <- setdiff(req_comm, names(comm_tbl))
if (length(miss_comm) > 0) stop("Community file missing columns: ", paste(miss_comm, collapse = ", "))

# Clean lineage: strip stray quotes, set blank -> UNKNOWN
comm_tbl[, lineage := gsub('^""|""$', '', lineage)]
comm_tbl[lineage == "" | is.na(lineage), lineage := "UNKNOWN"]

# layout: first three cols are sample, x, y (column names may vary)
names(layout_tbl)[1:3] <- c("sample", "x", "y")

cat(sprintf("  Edges: %d rows\n", nrow(edges)))
cat(sprintf("  Nodes (community table): %d rows\n", nrow(comm_tbl)))
cat(sprintf("  Layout entries: %d rows\n", nrow(layout_tbl)))

# =============================================================================
# STEP 2: BUILD IGRAPH OBJECT
# =============================================================================

cat("[2] Building igraph object...\n")

g <- graph_from_data_frame(
  d        = edges[, .(from, to, weight)],
  directed = FALSE,
  vertices = comm_tbl[, .(name = sample)]
)
E(g)$weight <- edges[match(paste(ends(g, E(g))[,1], ends(g, E(g))[,2]),
                           paste(edges$from, edges$to)), weight]

# Attach vertex attributes
v_idx <- match(V(g)$name, comm_tbl$sample)
V(g)$community <- comm_tbl$community[v_idx]
V(g)$lineage   <- comm_tbl$lineage[v_idx]

# =============================================================================
# STEP 3: VALIDATE AGAINST THESIS TARGETS
# =============================================================================

cat("[3] Validating against thesis targets...\n")

n_nodes <- vcount(g)
n_edges <- ecount(g)
n_comm  <- length(unique(V(g)$community))

cat(sprintf("  Nodes:       %d  (target: 167)%s\n",
            n_nodes, if (n_nodes == 167) "" else "  [MISMATCH]"))
cat(sprintf("  Edges:       %d  (target: 2,130)%s\n",
            n_edges, if (n_edges == 2130) "" else "  [MISMATCH]"))
cat(sprintf("  Communities: %d  (target: 5)%s\n",
            n_comm, if (n_comm == 5) "" else "  [MISMATCH]"))

q_recovered  <- NA_real_;  q_source <- "not supplied"
r_recovered  <- NA_real_;  r_source <- "not supplied"

if (!is.null(opt$modularity) && nzchar(opt$modularity) && file.exists(opt$modularity)) {
  mod_dt   <- fread(opt$modularity)
  louv_row <- mod_dt[method == "louvain"]
  if (nrow(louv_row) > 0) {
    q_recovered <- louv_row$modularity[1]
    q_source    <- opt$modularity
    cat(sprintf("  Modularity Q: %.6f  (target ≈ 0.629)%s\n",
                q_recovered, if (abs(q_recovered - 0.629) < 0.02) "" else "  [MISMATCH]"))
  }
} else {
  cat("  Modularity Q: (--modularity not supplied)\n")
}

if (!is.null(opt$assortativity) && nzchar(opt$assortativity) && file.exists(opt$assortativity)) {
  assort_dt  <- fread(opt$assortativity)
  assort_row <- assort_dt[grepl("assortativity", metric, ignore.case = TRUE)]
  if (nrow(assort_row) > 0) {
    r_recovered <- assort_row$value[1]
    r_source    <- opt$assortativity
    cat(sprintf("  Assortativity r: %.6f  (target ≈ 0.801)%s\n",
                r_recovered, if (abs(r_recovered - 0.801) < 0.02) "" else "  [see note]"))
  }
} else {
  cat("  Assortativity r: (--assortativity not supplied)\n")
}

cat("\n")

# =============================================================================
# STEP 4: NODE AESTHETICS
# =============================================================================

cat("[4] Setting up aesthetics...\n")

deg   <- degree(g)
vsize <- 2.5 + 7 * (deg - min(deg)) / (max(deg) - min(deg))

# Okabe-Ito palette (colourblind-safe) assigned to lineages in sorted order
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442")
lineages  <- setdiff(sort(unique(V(g)$lineage)), "UNKNOWN")
lin_col   <- setNames(okabe_ito[seq_along(lineages)], lineages)
lin_col["UNKNOWN"] <- "#999999"

# Dark2 palette scaled to number of communities
comm_ids <- as.character(sort(unique(V(g)$community)))
n_c      <- length(comm_ids)
comm_col <- setNames(
  colorRampPalette(brewer.pal(min(n_c, 8), "Dark2"))(n_c),
  comm_ids
)

node_lin  <- lin_col[V(g)$lineage]
node_comm <- comm_col[as.character(V(g)$community)]

# Degree legend values at four quantile positions
leg_deg <- round(quantile(deg, probs = c(0.1, 0.4, 0.7, 1.0)))
leg_vsz <- (2.5 + 7 * (leg_deg - min(deg)) / (max(deg) - min(deg))) * 0.52

# =============================================================================
# STEP 5: FIXED LAYOUT (deterministic geometry)
# =============================================================================

layout_m <- as.matrix(layout_tbl[match(V(g)$name, layout_tbl$sample), .(x, y)])

# Normalise to [-1, 1] range (needed for draw_hulls which uses igraph rescale)
lx_n <- 2 * (layout_m[, 1] - min(layout_m[, 1])) /
            (max(layout_m[, 1]) - min(layout_m[, 1])) - 1
ly_n <- 2 * (layout_m[, 2] - min(layout_m[, 2])) /
            (max(layout_m[, 2]) - min(layout_m[, 2])) - 1
layout_n <- cbind(lx_n, ly_n)

# =============================================================================
# STEP 6: DRAWING HELPERS
# =============================================================================

# Convex hulls expanded outward from each community centroid
draw_hulls <- function(lsc, mem, cols, expand = 1.18) {
  for (cid in sort(unique(mem))) {
    pts <- lsc[mem == cid, , drop = FALSE]
    col <- cols[as.character(cid)]
    if (nrow(pts) >= 3) {
      hi <- chull(pts); hp <- pts[hi, ]
      cx <- mean(hp[, 1]); cy <- mean(hp[, 2])
      hx <- cx + expand * (hp[, 1] - cx)
      hy <- cy + expand * (hp[, 2] - cy)
      polygon(c(hx, hx[1]), c(hy, hy[1]),
              col    = adjustcolor(col, 0.12),
              border = adjustcolor(col, 0.50),
              lwd    = 1.1)
    }
  }
}

# Centroid labels with filled background box
draw_labels <- function(lsc, mem, cols) {
  for (cid in sort(unique(mem))) {
    pts <- lsc[mem == cid, , drop = FALSE]
    cx  <- median(pts[, 1]); cy <- median(pts[, 2])
    col <- cols[as.character(cid)]
    lbl <- paste0("C", cid)
    tw  <- strwidth(lbl,  cex = 0.73) * 1.7
    th  <- strheight(lbl, cex = 0.73) * 2.1
    rect(cx - tw/2, cy - th/2, cx + tw/2, cy + th/2,
         col = adjustcolor(col, 0.88), border = NA)
    text(cx, cy, lbl, col = "white", font = 2, cex = 0.73)
  }
}

# Single panel plot using the fixed layout
plot_panel <- function(vcol, hull = FALSE) {
  plot(g,
       layout             = layout_m,
       rescale            = TRUE,
       asp                = 0,
       vertex.color       = vcol,
       vertex.size        = vsize,
       vertex.label       = NA,
       vertex.frame.color = "white",
       vertex.frame.width = 1.0,
       edge.width         = 0.35,
       edge.color         = adjustcolor("#C0C0C0", 0.30),
       xlim               = c(-1.15, 1.15),
       ylim               = c(-1.15, 1.15))
  if (hull) {
    par(xpd = TRUE)
    draw_hulls(layout_n, V(g)$community, comm_col)
    draw_labels(layout_n, V(g)$community, comm_col)
    par(xpd = FALSE)
  }
  usr <- par("usr")
  text(usr[1] + 0.02 * (usr[2] - usr[1]),
       usr[4] - 0.01 * (usr[4] - usr[3]),
       if (hull) "B   Louvain community" else "A   Lineage",
       font = 2, cex = 1.05, adj = c(0, 1))
}

# Bold legend section header
legend_header <- function(x, y, txt) {
  text(x, y, txt, font = 2, cex = 0.82, adj = c(0, 1))
}

# =============================================================================
# STEP 7: RENDER PDF AND PNG
# =============================================================================

OUT   <- opt[["out-dir"]]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fname <- "Fig_pan_cancer_cell_line_similarity_network_lineage_community"

cat(sprintf("[5] Writing figures to: %s\n", OUT))

for (fmt in c("pdf", "png")) {
  out_path <- file.path(OUT, paste0(fname, ".", fmt))

  if (fmt == "pdf") {
    pdf(out_path, width = 10.8, height = 5.0, useDingbats = FALSE)
  } else {
    png(out_path, width = 10.8, height = 5.0, units = "in", res = 600, bg = "white")
  }

  # Three-column layout: [Panel A] [Panel B] [Legends]
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4.2, 4.2, 2.4))

  # Panel A — lineage
  par(mar = c(0.4, 0.4, 1.6, 0.3), bg = "white")
  plot_panel(node_lin, hull = FALSE)

  # Panel B — community with convex hulls and centroid labels
  par(mar = c(0.4, 0.4, 1.6, 0.3), bg = "white")
  plot_panel(node_comm, hull = TRUE)

  # Legend column
  par(mar = c(0.4, 0.2, 1.6, 0.4), bg = "white")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  legend_header(0.04, 1.00, "Lineage")
  legend(0.04, 0.95,
         legend = names(lin_col), fill = unname(lin_col), border = "white",
         bty = "n", cex = 0.80, pt.cex = 1.1, y.intersp = 1.05, x.intersp = 0.6)

  legend_header(0.04, 0.58, "Louvain community")
  legend(0.04, 0.53,
         legend = paste0("C", comm_ids), fill = unname(comm_col), border = "white",
         bty = "n", cex = 0.80, pt.cex = 1.1, y.intersp = 1.05, x.intersp = 0.6)

  legend_header(0.04, 0.21, "Degree centrality")
  legend(0.04, 0.16,
         legend  = as.character(unname(leg_deg)),
         pch = 21, col = "grey35", pt.bg = "grey65", pt.cex = leg_vsz,
         bty = "n", cex = 0.80, y.intersp = 1.45, x.intersp = 0.6)

  dev.off()
  cat(sprintf("       %s\n", out_path))
}

# =============================================================================
# VALIDATION REPORT
# =============================================================================

cat("\n=================================================\n")
cat("VALIDATION REPORT\n")
cat("=================================================\n")
cat(sprintf("Figure generated:           YES\n"))
cat(sprintf("PDF: %s\n", file.path(OUT, paste0(fname, ".pdf"))))
cat(sprintf("PNG: %s\n", file.path(OUT, paste0(fname, ".png"))))
cat(sprintf("Edge input:                 %s\n", opt$edges))
cat(sprintf("Communities input:          %s\n", opt$communities))
cat(sprintf("Layout input:               %s\n", opt$layout))
cat(sprintf("Node count:                 %d  (thesis: 167)%s\n",
            n_nodes, if (n_nodes == 167) " [OK]" else " [MISMATCH]"))
cat(sprintf("Edge count:                 %d  (thesis: 2,130)%s\n",
            n_edges, if (n_edges == 2130) " [OK]" else " [MISMATCH]"))
cat(sprintf("Community count:            %d  (thesis: 5)%s\n",
            n_comm,  if (n_comm  == 5)    " [OK]" else " [MISMATCH]"))
if (!is.na(q_recovered)) {
  cat(sprintf("Modularity Q (Louvain):     %.6f  (thesis ≈ 0.629)  [from: %s]\n",
              q_recovered, q_source))
} else {
  cat(sprintf("Modularity Q:               not recovered (%s)\n", q_source))
}
if (!is.na(r_recovered)) {
  cat(sprintf("Lineage assortativity r:    %.6f  (thesis ≈ 0.801)  [from: %s]\n",
              r_recovered, r_source))
} else {
  cat(sprintf("Lineage assortativity r:    not recovered (%s)\n", r_source))
}
cat(sprintf("Existing scripts modified:  NO\n"))
cat("=================================================\n")
