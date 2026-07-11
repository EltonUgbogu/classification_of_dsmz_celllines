#!/usr/bin/env Rscript
# =============================================================================
# plot_pan_cancer_two_panel.R
#
# Two-panel publication figure for the 167-node DSMZ cell-line similarity
# network.
#
# Panel A: nodes coloured by annotated cancer type.
# Panel B: nodes coloured by Louvain community with convex hulls and centroid
#          labels C1-C5.
#
# Outputs (written to OUT):
#   Fig_pan_cancer_cell_line_similarity_network_cancer_type_community.pdf
#   Fig_pan_cancer_cell_line_similarity_network_cancer_type_community.png
# =============================================================================

suppressPackageStartupMessages({
  library(igraph); library(data.table); library(optparse)
})
set.seed(1)

option_list <- list(
  make_option("--edges",       type="character", default=NULL,
    help="[REQUIRED] Edge list TSV with columns: from, to, weight"),
  make_option("--communities", type="character", default=NULL,
    help="[REQUIRED] Community/metadata TSV with columns: sample, component, cancer_type"),
  make_option("--layout",      type="character", default=NULL,
    help="[REQUIRED] Layout TSV with columns: sample, x, y"),
  make_option("--out-dir",     type="character", default=NULL,
    help="[REQUIRED] Output directory for PDF and PNG"),
  make_option("--width", type="double", default=11.2,
    help="Output figure width in inches [default: %default]"),
  make_option("--height", type="double", default=5.0,
    help="Output figure height in inches [default: %default]"),
  make_option("--panel-width", type="double", default=4.75,
    help="Relative width for each network panel [default: %default]"),
  make_option("--legend-width", type="double", default=1.70,
    help="Relative width for the separate legend column [default: %default]"),
  make_option("--legend-mode", type="character", default="separate",
    help="Legend placement mode. Only 'separate' is supported for this figure [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (any(sapply(c("edges","communities","layout","out-dir"), function(k) is.null(opt[[k]])))) {
  stop("--edges, --communities, --layout, and --out-dir are all required", call.=FALSE)
}
OUT <- opt[["out-dir"]]
FIG_WIDTH <- opt[["width"]]
FIG_HEIGHT <- opt[["height"]]
PANEL_WIDTH <- opt[["panel-width"]]
LEGEND_WIDTH <- opt[["legend-width"]]
LEGEND_MODE <- opt[["legend-mode"]]
if (!identical(LEGEND_MODE, "separate")) {
  stop("--legend-mode must be 'separate' for this fixed two-panel layout", call.=FALSE)
}
if (any(!is.finite(c(FIG_WIDTH, FIG_HEIGHT, PANEL_WIDTH, LEGEND_WIDTH))) ||
    any(c(FIG_WIDTH, FIG_HEIGHT, PANEL_WIDTH, LEGEND_WIDTH) <= 0)) {
  stop("Plot width, height, panel width, and legend width must be positive finite values", call.=FALSE)
}

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
edges      <- fread(opt[["edges"]],       quote="")
comm_tbl   <- fread(opt[["communities"]], quote="")
layout_tbl <- fread(opt[["layout"]],      quote="")

if (!"cancer_type" %in% names(comm_tbl)) {
  if ("lineage" %in% names(comm_tbl)) {
    warning("Using legacy community column 'lineage' to populate 'cancer_type'")
    comm_tbl[, cancer_type := lineage]
  } else {
    stop("Community table must contain cancer_type", call.=FALSE)
  }
}
comm_tbl[, cancer_type := gsub('^""|""$', '', cancer_type)]
comm_tbl[cancer_type == "" | is.na(cancer_type), cancer_type := "UNKNOWN"]
setnames(comm_tbl, "component", "community")

# -----------------------------------------------------------------------------
# 2. Build graph
# -----------------------------------------------------------------------------
g <- graph_from_data_frame(edges[, .(from, to, weight)], directed=FALSE,
                           vertices=comm_tbl[, .(name=sample)])
V(g)$community <- comm_tbl[match(V(g)$name, comm_tbl$sample), community]
V(g)$cancer_type <- comm_tbl[match(V(g)$name, comm_tbl$sample), cancer_type]

# -----------------------------------------------------------------------------
# 3. Layout: fixed coordinates, pre-normalised to [-1,1] to match igraph
# -----------------------------------------------------------------------------
layout_m <- as.matrix(layout_tbl[match(V(g)$name, layout_tbl$sample), .(x, y)])
lx_n <- 2*(layout_m[,1]-min(layout_m[,1]))/(max(layout_m[,1])-min(layout_m[,1]))-1
ly_n <- 2*(layout_m[,2]-min(layout_m[,2]))/(max(layout_m[,2])-min(layout_m[,2]))-1
layout_n <- cbind(lx_n, ly_n)

# -----------------------------------------------------------------------------
# 4. Node aesthetics
# -----------------------------------------------------------------------------
deg   <- degree(g)
vsize <- 2.5 + 7*(deg - min(deg)) / (max(deg) - min(deg))

# Colour-blind-friendly qualitative palettes with green excluded.
cancer_type_palette <- c(
  BRCA = "#0072B2",
  HEME = "#E69F00",
  NBL  = "#D55E00",
  RBL  = "#CC79A7",
  UNKNOWN = "#666666"
)
fallback_palette <- c("#0072B2", "#E69F00", "#D55E00", "#CC79A7",
                      "#56B4E9", "#AA4499", "#882255", "#999999",
                      "#000000")
cancer_types <- setdiff(sort(unique(V(g)$cancer_type)), "UNKNOWN")
cancer_type_col <- cancer_type_palette[cancer_types]
if (any(is.na(cancer_type_col))) {
  missing_cancer_types <- cancer_types[is.na(cancer_type_col)]
  extra_cols <- setdiff(fallback_palette, unname(cancer_type_palette))
  if (length(missing_cancer_types) > length(extra_cols)) {
    stop("Not enough no-green fallback colours for cancer-type levels", call.=FALSE)
  }
  cancer_type_col[is.na(cancer_type_col)] <- extra_cols[seq_along(missing_cancer_types)]
}
if (any(V(g)$cancer_type == "UNKNOWN")) {
  cancer_type_col <- c(cancer_type_col, UNKNOWN = cancer_type_palette[["UNKNOWN"]])
}

# Louvain community palette, also green-free.
n_comm   <- length(unique(V(g)$community))
comm_ids <- as.character(sort(unique(V(g)$community)))
if (n_comm > length(fallback_palette)) {
  stop("Not enough no-green colours for Louvain communities", call.=FALSE)
}
comm_col <- setNames(fallback_palette[seq_len(n_comm)], comm_ids)

node_cancer_type <- cancer_type_col[V(g)$cancer_type]
node_comm <- comm_col[as.character(V(g)$community)]

# Degree legend: three quantile-spaced values across the observed range.
leg_deg <- unique(as.integer(round(quantile(deg, probs=c(0.1, 0.55, 1.0)))))
leg_vsz <- (2.5 + 7*(leg_deg - min(deg)) / (max(deg) - min(deg))) * 0.46

# -----------------------------------------------------------------------------
# 5. Helper functions
# -----------------------------------------------------------------------------

# Convex hulls expanded outward from each community centroid
draw_hulls <- function(lsc, mem, cols, expand=1.18) {
  for (cid in sort(unique(mem))) {
    pts <- lsc[mem == cid, , drop=FALSE]
    col <- cols[as.character(cid)]
    if (nrow(pts) >= 3) {
      hi <- chull(pts); hp <- pts[hi, ]
      cx <- mean(hp[,1]); cy <- mean(hp[,2])
      hx <- cx + expand*(hp[,1]-cx)
      hy <- cy + expand*(hp[,2]-cy)
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
    pts <- lsc[mem == cid, , drop=FALSE]
    cx  <- median(pts[,1]); cy <- median(pts[,2])
    col <- cols[as.character(cid)]
    lbl <- paste0("C", cid)
    tw  <- strwidth(lbl,  cex=0.73) * 1.7
    th  <- strheight(lbl, cex=0.73) * 2.1
    rect(cx-tw/2, cy-th/2, cx+tw/2, cy+th/2,
         col=adjustcolor(col, 0.88), border=NA)
    text(cx, cy, lbl, col="white", font=2, cex=0.73)
  }
}

# Single panel plot
plot_panel <- function(vcol, hull=FALSE) {
  plot(g,
       layout             = layout_m,
       rescale            = TRUE,
       asp                = 0,
       vertex.color       = vcol,
       vertex.size        = vsize,
       vertex.label       = NA,
       vertex.frame.color = "white",
       vertex.frame.width = 0.8,
       edge.width         = 0.26,
       edge.color         = adjustcolor("#777777", 0.20),
       xlim               = c(-1.04, 1.04),
       ylim               = c(-1.04, 1.04))
  if (hull) {
    par(xpd=TRUE)
    draw_hulls(layout_n, V(g)$community, comm_col)
    draw_labels(layout_n, V(g)$community, comm_col)
    par(xpd=FALSE)
  }
  usr <- par("usr")
  text(usr[1] + 0.02*(usr[2]-usr[1]),
       usr[4] - 0.01*(usr[4]-usr[3]),
       if (hull) "B   Louvain community" else "A   Annotated cancer type",
       font=2, cex=0.98, adj=c(0,1))
}

# Legend section header (avoids \n in legend title which causes font warnings)
legend_header <- function(x, y, txt) {
  text(x, y, txt, font=2, cex=0.82, adj=c(0,1))
}

# -----------------------------------------------------------------------------
# 6. Render PDF and PNG
# -----------------------------------------------------------------------------
fname <- "Fig_pan_cancer_cell_line_similarity_network_cancer_type_community"

for (fmt in c("pdf", "png")) {
  if (fmt == "pdf") {
    pdf(file.path(OUT, paste0(fname, ".pdf")),
        width=FIG_WIDTH, height=FIG_HEIGHT, useDingbats=FALSE)
  } else {
    png_path <- file.path(OUT, paste0(fname, ".png"))
    if (requireNamespace("ragg", quietly=TRUE)) {
      ragg::agg_png(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, units="in",
                    res=600, background="white")
    } else {
      png(png_path, width=FIG_WIDTH, height=FIG_HEIGHT, units="in", res=600, bg="white")
    }
  }

  # Three-column layout: Panel A | Panel B | Legends
  layout(matrix(c(1, 2, 3), nrow=1), widths=c(PANEL_WIDTH, PANEL_WIDTH, LEGEND_WIDTH))

  # Panel A: annotated cancer type.
  par(mar=c(0.15, 0.15, 1.15, 0.10), bg="white")
  plot_panel(node_cancer_type, hull=FALSE)

  # Panel B — community with hulls and labels
  par(mar=c(0.15, 0.15, 1.15, 0.10), bg="white")
  plot_panel(node_comm, hull=TRUE)

  # Legend column
  par(mar=c(0.15, 0.05, 1.15, 0.15), bg="white")
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))

  legend_header(0.04, 1.00, "Cancer type")
  legend(0.04, 0.95,
         legend=names(cancer_type_col), fill=unname(cancer_type_col), border="white",
         bty="n", cex=0.68, pt.cex=0.95, y.intersp=0.92, x.intersp=0.48)

  legend_header(0.04, 0.66, "Louvain community")
  legend(0.04, 0.61,
         legend=paste0("C", comm_ids), fill=unname(comm_col), border="white",
         bty="n", cex=0.68, pt.cex=0.95, y.intersp=0.92, x.intersp=0.48)

  legend_header(0.04, 0.34, "Node size: degree")
  degree_y <- seq(0.255, 0.055, length.out=length(leg_deg))
  old_xpd <- par("xpd")
  par(xpd=NA)
  points(rep(0.16, length(leg_deg)), degree_y,
         pch=21, col="grey35", bg="grey65", cex=leg_vsz)
  text(rep(0.34, length(leg_deg)), degree_y,
       labels=as.character(unname(leg_deg)), cex=0.68, adj=c(0, 0.5))
  par(xpd=old_xpd)

  dev.off()
  if (fmt == "png" && requireNamespace("png", quietly=TRUE)) {
    tryCatch(
      png::readPNG(file.path(OUT, paste0(fname, ".png"))),
      error=function(e) {
        stop("PNG validation failed for ",
             file.path(OUT, paste0(fname, ".png")), ": ",
             conditionMessage(e), call.=FALSE)
      }
    )
  }
  cat("Saved:", file.path(OUT, paste0(fname, ".", fmt)), "\n")
}
