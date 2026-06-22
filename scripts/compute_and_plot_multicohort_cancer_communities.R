#!/usr/bin/env Rscript

# =============================================================================
# compute_and_plot_multicohort_cancer_communities.R
#
# Unweighted Leiden community detection on the final resolved MULTICOHORT_CANCER
# DSMZ cell-line graph.
#
# Expected graph: 56 nodes, 155 edges, 9 connected components.
# Isolates: DU_4475, SKNBE2, SK_BR_3.
#
# Panel A: lineage colours (Okabe-Ito)
# Panel B: Leiden community colours (Dark2)
#
# IMPORTANT: Uses igraph::cluster_leiden() for unweighted Leiden community
# detection. Edge support columns in the input are blank/fallback values from
# resolved neighbours and are intentionally NOT used as weights.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
  library(RColorBrewer)
})

# =============================================================================
# CLI arguments
# =============================================================================
opt_list <- list(
  make_option("--edges",                  type = "character", default = NULL),
  make_option("--node-stats",             type = "character", default = NULL),
  make_option("--anchor-audit",           type = "character", default = NULL),
  make_option("--shortnames",             type = "character", default = NULL),
  make_option("--metadata",               type = "character", default = NULL),
  make_option("--helper-node-annotations",type = "character", default = NULL),
  make_option("--out-communities",        type = "character", default = NULL),
  make_option("--out-community-summary",  type = "character", default = NULL),
  make_option("--out-modularity",         type = "character", default = NULL),
  make_option("--out-layout",             type = "character", default = NULL),
  make_option("--out-pdf",                type = "character", default = NULL),
  make_option("--out-png",                type = "character", default = NULL),
  make_option("--seed",                   type = "integer",   default = 42L)
)
opt <- parse_args(OptionParser(option_list = opt_list))

required_args <- c("edges", "node-stats", "metadata",
                   "out-communities", "out-community-summary",
                   "out-modularity", "out-layout", "out-pdf", "out-png")
missing_args <- required_args[sapply(required_args, function(a) is.null(opt[[a]]))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

for (f in c(opt[["edges"]], opt[["node-stats"]], opt[["metadata"]])) {
  if (!is.null(f) && !file.exists(f)) stop("Input not found: ", f)
}

# =============================================================================
# Check Leiden is available before proceeding
# =============================================================================
if (!exists("cluster_leiden", where = asNamespace("igraph"), inherits = FALSE)) {
  stop(
    "cluster_leiden() is not available in the installed igraph version (",
    as.character(packageVersion("igraph")), "). ",
    "Leiden is required. Do not substitute Louvain. ",
    "Upgrade r-igraph to >= 1.3.0 in envs/tcga-r-env.yaml."
  )
}
message("[INFO] igraph version: ", packageVersion("igraph"), " — cluster_leiden() available")

set.seed(opt$seed)

# =============================================================================
# Helper: lineage from cancer_type string
# =============================================================================
lineage_from_string <- function(x) {
  out <- rep(NA_character_, length(x))
  out[grepl("BRCA|Breast",          x, ignore.case = TRUE)] <- "BRCA"
  out[grepl("NBL|Neuroblast",       x, ignore.case = TRUE)] <- "NBL"
  out[grepl("RBL|Retinoblast",      x, ignore.case = TRUE)] <- "RBL"
  out[grepl("HEME|Hema|Leukemia|Lymphoma|LL-100", x, ignore.case = TRUE)] <- "HEME"
  out
}

# =============================================================================
# 1. Load edges
# =============================================================================
message("[INFO] Loading edges: ", opt[["edges"]])
edges_dt <- fread(opt[["edges"]], quote = "")
stopifnot(all(c("node1", "node2") %in% names(edges_dt)))
edges_dt <- unique(edges_dt[, .(node1, node2)])
message("[INFO] Edge rows (after dedup): ", nrow(edges_dt))

# =============================================================================
# 2. Load node stats
# =============================================================================
message("[INFO] Loading node stats: ", opt[["node-stats"]])
ns_dt <- fread(opt[["node-stats"]], quote = "")
stopifnot("cell_line" %in% names(ns_dt))
setkey(ns_dt, cell_line)
message("[INFO] Node stats rows: ", nrow(ns_dt))

# =============================================================================
# 3. Build igraph (unweighted, undirected)
# =============================================================================
all_nodes <- unique(c(ns_dt$cell_line,
                      edges_dt$node1, edges_dt$node2))
g <- graph_from_data_frame(
  d        = edges_dt[, .(from = node1, to = node2)],
  directed = FALSE,
  vertices = data.frame(name = all_nodes, stringsAsFactors = FALSE)
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

n_nodes <- vcount(g)
n_edges <- ecount(g)
message("[INFO] Graph: ", n_nodes, " nodes, ", n_edges, " edges")

# Validate expected graph properties
if (n_nodes != 56) warning("[WARN] Expected 56 nodes, found ", n_nodes)
if (n_edges != 155) warning("[WARN] Expected 155 edges, found ", n_edges)

# =============================================================================
# 4. Attach node metadata from node_stats
# =============================================================================
ns_match <- ns_dt[match(V(g)$name, ns_dt$cell_line), ]
V(g)$degree_ns     <- ifelse(is.na(ns_match$degree), 0L, as.integer(ns_match$degree))
V(g)$component     <- ifelse(is.na(ns_match$component), NA_character_, as.character(ns_match$component))
V(g)$is_isolate    <- ifelse(is.na(ns_match$is_isolate), FALSE, ns_match$is_isolate == "True")
V(g)$is_central    <- ifelse(is.na(ns_match$is_central), FALSE, ns_match$is_central == "True")

bridge_col <- if ("canonical_bridge_selected" %in% names(ns_match)) {
  ns_match$canonical_bridge_selected == "True"
} else {
  rep(FALSE, vcount(g))
}
bridge_col[is.na(bridge_col)] <- FALSE
V(g)$canonical_bridge_selected <- bridge_col

# =============================================================================
# 5. Attach anchor/bridge audit (optional, graceful if absent)
# =============================================================================
if (!is.null(opt[["anchor-audit"]]) && file.exists(opt[["anchor-audit"]])) {
  message("[INFO] Loading anchor audit: ", opt[["anchor-audit"]])
  anc_dt <- fread(opt[["anchor-audit"]], quote = "")
  if ("node_id" %in% names(anc_dt) && "canonical_bridge_selected" %in% names(anc_dt)) {
    anc_match <- anc_dt[match(V(g)$name, anc_dt$node_id), ]
    bridgev <- anc_match$canonical_bridge_selected == "True"
    bridgev[is.na(bridgev)] <- FALSE
    V(g)$canonical_bridge_selected <- bridgev
    message("[INFO] Bridge/anchor annotations attached from anchor audit")
  }
} else {
  message("[INFO] Anchor audit not found or not provided — using node_stats bridge column")
}

# =============================================================================
# 6. Map lineage from joint_metadata
# =============================================================================
message("[INFO] Loading metadata: ", opt[["metadata"]])
meta_dt <- fread(opt[["metadata"]], quote = "")
# Extract cell_line name from sample_id (format: NG-XXXXX_CellLineName_libYYYYY...)
# Filter to DSMZ cell lines only
if ("sample_type" %in% names(meta_dt)) {
  meta_cl <- meta_dt[grepl("Cell Line", sample_type, ignore.case = TRUE)]
} else {
  meta_cl <- meta_dt
}
if ("sample_id" %in% names(meta_cl) && "cancer_type" %in% names(meta_cl)) {
  meta_cl[, cell_line_extracted := gsub("_lib.*$", "",
                                         gsub("^NG-[A-Za-z0-9]+_", "", sample_id))]
  # Deduplicate: take first cancer_type per extracted cell line
  meta_uniq <- unique(meta_cl[, .(cell_line_extracted, cancer_type)])[
    , .(cancer_type = cancer_type[1]), by = cell_line_extracted]
  lineage_map <- setNames(lineage_from_string(meta_uniq$cancer_type),
                           meta_uniq$cell_line_extracted)
  V(g)$lineage <- ifelse(
    V(g)$name %in% names(lineage_map),
    lineage_map[V(g)$name],
    NA_character_
  )
  message("[INFO] Lineage mapped from metadata: ",
          sum(!is.na(V(g)$lineage)), " / ", vcount(g), " nodes")
} else {
  V(g)$lineage <- NA_character_
  message("[WARN] metadata missing sample_id or cancer_type columns — lineage will be UNKNOWN")
}

# Fallback: helper node annotations
n_missing_lineage <- sum(is.na(V(g)$lineage))
if (n_missing_lineage > 0 && !is.null(opt[["helper-node-annotations"]]) &&
    file.exists(opt[["helper-node-annotations"]])) {
  message("[INFO] Filling ", n_missing_lineage,
          " missing lineages from helper node annotations: ",
          opt[["helper-node-annotations"]])
  help_dt <- fread(opt[["helper-node-annotations"]], quote = "")
  if ("cell_line" %in% names(help_dt)) {
    for (col in c("lineage", "cancer_type", "Disease")) {
      if (col %in% names(help_dt)) {
        help_map <- setNames(lineage_from_string(as.character(help_dt[[col]])),
                             help_dt$cell_line)
        missing_idx <- is.na(V(g)$lineage)
        V(g)$lineage[missing_idx] <- help_map[V(g)$name[missing_idx]]
        break
      }
    }
  }
}

# Final fallback
V(g)$lineage[is.na(V(g)$lineage) | V(g)$lineage == ""] <- "UNKNOWN"
message("[INFO] Lineage distribution: ",
        paste(names(table(V(g)$lineage)), table(V(g)$lineage), sep = "=", collapse = ", "))

# =============================================================================
# 7. Short names
# =============================================================================
if (!is.null(opt[["shortnames"]]) && file.exists(opt[["shortnames"]])) {
  sn_dt <- fread(opt[["shortnames"]], quote = "")
  short_col <- if ("short_id" %in% names(sn_dt)) "short_id" else
               if ("cell_line_short" %in% names(sn_dt)) "cell_line_short" else
               if ("display_label" %in% names(sn_dt)) "display_label" else NULL
  long_col  <- if ("long_id" %in% names(sn_dt)) "long_id" else
               if ("cell_line" %in% names(sn_dt)) "cell_line" else NULL
  if (!is.null(short_col) && !is.null(long_col)) {
    sn_map <- setNames(sn_dt[[short_col]], sn_dt[[long_col]])
    V(g)$cell_line_short <- ifelse(V(g)$name %in% names(sn_map),
                                   sn_map[V(g)$name], V(g)$name)
  } else {
    V(g)$cell_line_short <- V(g)$name
  }
} else {
  V(g)$cell_line_short <- V(g)$name
}

# =============================================================================
# 8. Graph metrics (degree, betweenness, components)
# =============================================================================
V(g)$degree_igraph    <- degree(g)
V(g)$betweenness      <- betweenness(g, normalized = TRUE)
V(g)$component_igraph <- as.character(components(g)$membership)

n_components <- components(g)$no
message("[INFO] Connected components: ", n_components)
if (n_components != 9) warning("[WARN] Expected 9 connected components, found ", n_components)

# Validate isolates
isolate_expected <- c("DU_4475", "SKNBE2", "SK_BR_3")
isolate_found    <- V(g)$name[V(g)$degree_igraph == 0]
message("[INFO] Isolates (degree=0): ", paste(sort(isolate_found), collapse = ", "))
missing_isolates <- setdiff(isolate_expected, isolate_found)
extra_isolates   <- setdiff(isolate_found, isolate_expected)
if (length(missing_isolates) > 0)
  warning("[WARN] Expected isolates not found as degree-0: ",
          paste(missing_isolates, collapse = ", "))
if (length(extra_isolates) > 0)
  warning("[WARN] Extra degree-0 nodes not expected: ",
          paste(extra_isolates, collapse = ", "))

# =============================================================================
# 9. Unweighted Leiden community detection
# =============================================================================
message("[INFO] Running unweighted Leiden community detection (seed=", opt$seed, ")")
set.seed(opt$seed)
leiden_comm <- cluster_leiden(
  g,
  objective_function  = "modularity",
  weights             = NULL,    # unweighted — support columns are blank/fallback
  resolution           = 1.0,
  n_iterations        = 10L
)

V(g)$community_leiden <- as.character(membership(leiden_comm))
leiden_modularity     <- modularity(g, membership(leiden_comm))
n_leiden_communities  <- length(unique(V(g)$community_leiden))

message("[INFO] Leiden communities: ", n_leiden_communities)
message("[INFO] Leiden modularity:  ", round(leiden_modularity, 4))

# =============================================================================
# 10. Layout (seed 42, computed once, saved)
# =============================================================================
message("[INFO] Computing graph layout (seed=", opt$seed, ")")
set.seed(opt$seed)
layout_mat <- layout_with_fr(g)
rownames(layout_mat) <- V(g)$name

# =============================================================================
# 11. Build node-level community table
# =============================================================================
comm_tbl <- data.table(
  cell_line                = V(g)$name,
  cell_line_short          = V(g)$cell_line_short,
  lineage                  = V(g)$lineage,
  component                = V(g)$component_igraph,
  community_leiden         = V(g)$community_leiden,
  degree                   = V(g)$degree_igraph,
  betweenness              = round(V(g)$betweenness, 6),
  is_isolate               = V(g)$is_isolate,
  is_central               = V(g)$is_central,
  canonical_bridge_selected = V(g)$canonical_bridge_selected
)

# =============================================================================
# 12. Community summary table
# =============================================================================
comm_summary <- comm_tbl[, {
  members     <- .SD$cell_line
  lin_tab     <- table(.SD$lineage)
  maj_lin     <- names(which.max(lin_tab))
  maj_lin_n   <- as.integer(max(lin_tab))
  # within-community edges
  node_set    <- .SD$cell_line
  subg        <- induced_subgraph(g, vids = which(V(g)$name %in% node_set))
  within_edges <- ecount(subg)
  list(
    community_size        = .N,
    within_community_edges = within_edges,
    majority_lineage      = maj_lin,
    majority_lineage_n    = maj_lin_n,
    lineage_purity        = round(maj_lin_n / .N, 3),
    members               = paste(sort(members), collapse = ";"),
    lineages_present      = paste(names(lin_tab), lin_tab, sep = "=", collapse = ";")
  )
}, by = community_leiden][order(as.integer(community_leiden))]

# =============================================================================
# 13. Modularity table
# =============================================================================
mod_tbl <- data.table(
  graph_id              = "MULTICOHORT_CANCER_resolved_DSMZ_cellline_graph",
  community_method      = "Leiden",
  weighted              = FALSE,
  n_nodes               = n_nodes,
  n_edges               = n_edges,
  n_connected_components = n_components,
  n_communities         = n_leiden_communities,
  modularity            = round(leiden_modularity, 6),
  seed                  = opt$seed,
  edge_source           = "multicohort_cancer_resolved_cell_line_neighbourhood_graph_edges.tsv"
)

# =============================================================================
# 14. Layout table
# =============================================================================
layout_tbl <- data.table(
  cell_line                 = V(g)$name,
  x                         = layout_mat[, 1],
  y                         = layout_mat[, 2],
  lineage                   = V(g)$lineage,
  community_leiden          = V(g)$community_leiden,
  degree                    = V(g)$degree_igraph,
  component                 = V(g)$component_igraph,
  is_isolate                = V(g)$is_isolate,
  is_central                = V(g)$is_central,
  canonical_bridge_selected = V(g)$canonical_bridge_selected
)

# =============================================================================
# 15. Write TSV outputs
# =============================================================================
for (path in c(opt[["out-communities"]], opt[["out-community-summary"]],
               opt[["out-modularity"]], opt[["out-layout"]])) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
fwrite(comm_tbl,     opt[["out-communities"]],       sep = "\t", quote = FALSE)
fwrite(comm_summary, opt[["out-community-summary"]], sep = "\t", quote = FALSE)
fwrite(mod_tbl,      opt[["out-modularity"]],        sep = "\t", quote = FALSE)
fwrite(layout_tbl,   opt[["out-layout"]],            sep = "\t", quote = FALSE)
message("[INFO] Community TSVs written")

# =============================================================================
# 16. Figure: colour palettes
# =============================================================================
okabe_ito <- c(BRCA = "#E69F00", NBL = "#56B4E9", RBL = "#009E73",
               HEME = "#CC79A7", UNKNOWN = "#999999")
lineages_present <- sort(unique(V(g)$lineage))
lin_col <- okabe_ito[lineages_present]
names(lin_col) <- lineages_present
# Fill any lineage not in Okabe-Ito set (shouldn't occur, but be safe)
extra_lin <- setdiff(lineages_present, names(okabe_ito))
if (length(extra_lin) > 0) {
  extra_cols <- colorRampPalette(c("#8dd3c7", "#fb8072", "#80b1d3"))(length(extra_lin))
  lin_col[extra_lin] <- extra_cols
}

comm_ids <- as.character(sort(unique(as.integer(V(g)$community_leiden))))
n_comm   <- length(comm_ids)
set.seed(opt$seed)
comm_col <- setNames(
  colorRampPalette(brewer.pal(min(n_comm, 8), "Dark2"))(n_comm),
  comm_ids
)

node_lin_col  <- lin_col[V(g)$lineage]
node_comm_col <- comm_col[V(g)$community_leiden]

# Node sizes: degree-scaled, isolates get minimum size
deg      <- V(g)$degree_igraph
deg_range <- max(deg) - min(deg)
vsize    <- if (deg_range > 0) {
  3.0 + 8.0 * (deg - min(deg)) / deg_range
} else {
  rep(4.0, vcount(g))
}

# Bridge/anchor outline thickness
frame_width <- ifelse(V(g)$canonical_bridge_selected, 2.5, 1.0)
frame_col   <- ifelse(V(g)$canonical_bridge_selected, "#222222", "#FFFFFF")

# Normalise layout to [-1,1]
lx <- layout_mat[, 1]; ly <- layout_mat[, 2]
norm01 <- function(v) {
  r <- range(v); if (diff(r) == 0) return(rep(0, length(v)))
  2 * (v - r[1]) / diff(r) - 1
}
layout_norm <- cbind(norm01(lx), norm01(ly))

# Convex hull helper (only for communities with >= 3 members)
draw_hulls <- function(layout_n, mem, cols, expand = 1.18) {
  for (cid in sort(unique(mem))) {
    idx <- which(mem == cid)
    pts <- layout_n[idx, , drop = FALSE]
    col <- cols[as.character(cid)]
    if (nrow(pts) >= 3) {
      hi <- chull(pts); hp <- pts[hi, ]
      cx <- mean(hp[, 1]); cy <- mean(hp[, 2])
      hx <- cx + expand * (hp[, 1] - cx)
      hy <- cy + expand * (hp[, 2] - cy)
      polygon(c(hx, hx[1]), c(hy, hy[1]),
              col    = adjustcolor(col, 0.10),
              border = adjustcolor(col, 0.50),
              lwd    = 1.0)
    }
  }
}

draw_comm_labels <- function(layout_n, mem, comm_ids_ordered) {
  for (i in seq_along(comm_ids_ordered)) {
    cid <- comm_ids_ordered[i]
    idx <- which(mem == cid)
    pts <- layout_n[idx, , drop = FALSE]
    cx  <- median(pts[, 1]); cy <- median(pts[, 2])
    col <- comm_col[as.character(cid)]
    lbl <- paste0("C", i)
    tw  <- strwidth(lbl, cex = 0.70) * 1.8
    th  <- strheight(lbl, cex = 0.70) * 2.2
    rect(cx - tw/2, cy - th/2, cx + tw/2, cy + th/2,
         col = adjustcolor(col, 0.88), border = NA)
    text(cx, cy, lbl, col = "white", font = 2, cex = 0.70)
  }
}

# Label isolates explicitly: small offset text
label_isolates <- function() {
  iso_idx <- which(V(g)$is_isolate | V(g)$degree_igraph == 0)
  if (length(iso_idx) == 0) return(invisible(NULL))
  for (i in iso_idx) {
    x <- layout_mat[i, 1]; y <- layout_mat[i, 2]
    # ggraph/igraph rescale: we use raw layout_mat here since plot() rescales
    # Label will be drawn by igraph's vertex.label mechanism; this is a fallback
  }
}

# Panel plot function
plot_panel <- function(vcol, panel_label, hull = FALSE) {
  plot(g,
       layout              = layout_mat,
       rescale             = TRUE,
       asp                 = 0,
       vertex.color        = vcol,
       vertex.size         = vsize,
       vertex.label        = V(g)$cell_line_short,
       vertex.label.cex    = 0.40,
       vertex.label.color  = "#111111",
       vertex.label.family = "sans",
       vertex.frame.color  = frame_col,
       vertex.frame.width  = frame_width,
       edge.width          = 0.30,
       edge.color          = adjustcolor("#C0C0C0", 0.35),
       xlim                = c(-1.18, 1.18),
       ylim                = c(-1.18, 1.18))
  if (hull) {
    par(xpd = TRUE)
    draw_hulls(layout_norm, V(g)$community_leiden, comm_col)
    draw_comm_labels(layout_norm, V(g)$community_leiden, comm_ids)
    par(xpd = FALSE)
  }
  usr <- par("usr")
  text(usr[1] + 0.02 * (usr[2] - usr[1]),
       usr[4] - 0.01 * (usr[4] - usr[3]),
       panel_label,
       font = 2, cex = 1.05, adj = c(0, 1))
}

# =============================================================================
# 17. Render figures
# =============================================================================
for (fmt in c("pdf", "png")) {
  out_path <- if (fmt == "pdf") opt[["out-pdf"]] else opt[["out-png"]]
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  if (fmt == "pdf") {
    pdf(out_path, width = 12.0, height = 5.5, useDingbats = FALSE)
  } else {
    png(out_path, width = 12.0, height = 5.5,
        units = "in", res = 300, bg = "white")
  }

  # Three-column layout: Panel A | Panel B | Legends
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4.5, 4.5, 3.0))

  # Panel A — lineage
  par(mar = c(0.4, 0.4, 1.8, 0.3), bg = "white")
  plot_panel(node_lin_col, "A   Lineage", hull = FALSE)

  # Panel B — Leiden community
  par(mar = c(0.4, 0.4, 1.8, 0.3), bg = "white")
  plot_panel(node_comm_col, "B   Leiden community", hull = TRUE)

  # Legend column
  par(mar = c(0.4, 0.2, 1.8, 0.4), bg = "white")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # Lineage legend
  text(0.04, 1.00, "Lineage", font = 2, cex = 0.85, adj = c(0, 1))
  legend(0.04, 0.94,
         legend  = names(lin_col),
         fill    = unname(lin_col),
         border  = "white",
         bty     = "n",
         cex     = 0.80,
         pt.cex  = 1.1,
         y.intersp = 1.05,
         x.intersp = 0.6)

  # Community legend
  comm_y_start <- 0.94 - length(lin_col) * 0.098 - 0.06
  text(0.04, comm_y_start, "Leiden community", font = 2, cex = 0.85, adj = c(0, 1))
  comm_labels <- paste0("C", seq_along(comm_ids),
                        " (n=", comm_summary$community_size[
                          match(comm_ids, comm_summary$community_leiden)], ")")
  legend(0.04, comm_y_start - 0.02,
         legend  = comm_labels,
         fill    = unname(comm_col),
         border  = "white",
         bty     = "n",
         cex     = 0.78,
         pt.cex  = 1.0,
         y.intersp = 1.05,
         x.intersp = 0.6)

  # Degree size legend
  deg_y_start <- comm_y_start - 0.02 - length(comm_ids) * 0.086 - 0.06
  if (deg_y_start > 0.05) {
    text(0.04, deg_y_start, "Node size = degree", font = 2, cex = 0.78, adj = c(0, 1))
    leg_deg <- round(quantile(deg, probs = c(0.1, 0.5, 1.0)))
    leg_vsz <- if (deg_range > 0) {
      (3.0 + 8.0 * (leg_deg - min(deg)) / deg_range) * 0.52
    } else {
      rep(2.0, 3)
    }
    legend(0.04, deg_y_start - 0.02,
           legend  = paste0("deg=", leg_deg),
           pt.cex  = leg_vsz,
           pch     = 21,
           col     = "#555555",
           pt.bg   = "#AAAAAA",
           bty     = "n",
           cex     = 0.75,
           y.intersp = 1.3)
  }

  # Bridge node note
  text(0.04, 0.06,
       "Bold border = bridge/anchor node",
       cex = 0.68, adj = c(0, 0.5), col = "#444444")

  dev.off()
  message("[INFO] Written: ", out_path)
}

# =============================================================================
# 18. Validation summary
# =============================================================================
message("\n=== VALIDATION SUMMARY ===")
message("Graph:           ", n_nodes, " nodes, ", n_edges, " edges")
message("Components:      ", n_components)
message("Isolates (deg=0):", paste(sort(isolate_found), collapse = ", "))
message("Leiden comms:    ", n_leiden_communities)
message("Modularity:      ", round(leiden_modularity, 4))
message("Lineage dist:    ", paste(names(table(V(g)$lineage)),
                                    table(V(g)$lineage), sep = "=", collapse = ", "))
message("Bridge nodes:    ", sum(V(g)$canonical_bridge_selected))
message("Layout saved:    ", opt[["out-layout"]])
message("PDF:             ", opt[["out-pdf"]])
message("PNG:             ", opt[["out-png"]])
message("==========================\n")

message("[INFO] Done.")
