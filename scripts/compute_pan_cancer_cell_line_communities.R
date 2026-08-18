#!/usr/bin/env Rscript
# =============================================================================
# compute_pan_cancer_cell_line_communities.R
# =============================================================================
#
# Runs weighted Louvain community detection and weighted Leiden community
# detection on the cell-line-only pan-cancer similarity graph. Louvain is the
# primary partition used by the figure; Leiden is written as an algorithmic
# comparison. A Leiden resolution sweep is also written as a sensitivity
# analysis and does not replace the primary Louvain partition.
#
# Outputs:
#   pan_cancer_cell_line_communities.tsv
#   pan_cancer_cell_line_leiden_communities.tsv
#   pan_cancer_cell_line_community_metrics.tsv
#   pan_cancer_cell_line_lineage_discordant_profiles.tsv
#   pan_cancer_cell_line_leiden_resolution_sweep_assignments.tsv
#   pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv
#   pan_cancer_cell_line_leiden_resolution_sweep.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option("--edges", type="character", default=NULL,
    help="[REQUIRED] Edge list TSV (from, to, weight)"),
  make_option("--meta", type="character", default=NULL,
    help="[REQUIRED] Node metadata TSV (sample_id, lineage, ...)"),
  make_option("--out", type="character", default=NULL,
    help="[REQUIRED] Louvain community TSV path"),
  make_option("--leiden-out", type="character", default=NULL,
    help="[REQUIRED] Leiden community TSV path"),
  make_option("--community-metrics-out", type="character", default=NULL,
    help="[REQUIRED] Community metrics TSV path"),
  make_option("--lineage-discordant-out", type="character", default=NULL,
    help="[REQUIRED] Lineage-discordant profile TSV path"),
  make_option("--seed", type="integer", default=1,
    help="Random seed for Louvain and Leiden (default: 1)"),
  make_option("--leiden-resolution", type="double", default=1.0,
    help="Legacy single Leiden resolution parameter (default: 1.0)"),
  make_option("--leiden-resolution-sweep", type="character", default="",
    help="Comma-separated Leiden resolution sweep values"),
  make_option("--leiden-sweep-assignments-out", type="character", default="",
    help="Long-format Leiden resolution sweep assignments TSV path"),
  make_option("--leiden-sweep-summary-out", type="character", default="",
    help="Resolution-level Leiden sweep summary TSV path"),
  make_option("--leiden-sweep-plot-out", type="character", default="",
    help="Optional Leiden sweep PDF path"),
  make_option("--expected-communities", type="integer", default=0,
    help="Expected Louvain community count; 0 = skip check (default: 0)")
)
opt <- parse_args(OptionParser(option_list=option_list))

required <- c(
  "edges", "meta", "out", "leiden-out", "community-metrics-out",
  "lineage-discordant-out"
)
if (any(vapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]), logical(1)))) {
  stop("Missing required argument(s): ", paste(required, collapse=", "))
}
if (!file.exists(opt$edges)) stop("edges file not found: ", opt$edges)
if (!file.exists(opt$meta))  stop("meta file not found: ", opt$meta)

for (path in c(opt$out, opt[["leiden-out"]], opt[["community-metrics-out"]],
               opt[["lineage-discordant-out"]],
               opt[["leiden-sweep-assignments-out"]],
               opt[["leiden-sweep-summary-out"]],
               opt[["leiden-sweep-plot-out"]])) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  }
}

parse_resolution_sweep <- function(raw, legacy_resolution) {
  if (is.null(raw) || !nzchar(trimws(raw))) {
    return(legacy_resolution)
  }
  parts <- trimws(unlist(strsplit(raw, "[,;[:space:]]+")))
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  if (length(vals) == 0L || any(!is.finite(vals))) {
    stop("Invalid --leiden-resolution-sweep values: ", raw)
  }
  vals <- unique(c(vals, legacy_resolution))
  vals[order(vals)]
}

derive_cell_line <- function(profile_id) {
  out <- sub("^NG-[^_]+_", "", profile_id)
  out <- sub("_lib.*$", "", out)
  out
}

format_numeric <- function(x) {
  format(x, scientific=FALSE, trim=TRUE)
}

resolution_sweep <- parse_resolution_sweep(
  opt[["leiden-resolution-sweep"]],
  opt[["leiden-resolution"]]
)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading edges:", opt$edges, "\n")
edges <- fread(opt$edges)
required_edge_cols <- c("from", "to", "weight")
if (!all(required_edge_cols %in% names(edges))) {
  stop("edges must contain columns: ", paste(required_edge_cols, collapse=", "))
}
E_und <- edges[, .(weight=max(as.numeric(weight), na.rm=TRUE)),
               by=.(from=pmin(from, to), to=pmax(from, to))]
setorder(E_und, from, to)
cat("  Edges (undirected):", nrow(E_und), "\n")

cat("[2] Loading metadata:", opt$meta, "\n")
meta <- fread(opt$meta)
required_meta_cols <- c("sample_id", "lineage")
if (!all(required_meta_cols %in% names(meta))) {
  stop("metadata must contain columns: ", paste(required_meta_cols, collapse=", "))
}
meta <- unique(meta, by="sample_id")

# ---------------------------------------------------------------------------
# 2. Build igraph
# ---------------------------------------------------------------------------
cat("[3] Building igraph...\n")
vertices <- data.table(name=meta$sample_id)
g <- graph_from_data_frame(E_und, directed=FALSE, vertices=vertices)
E(g)$weight <- E_und$weight
cat("  Nodes:", vcount(g), "  Edges:", ecount(g), "\n")

meta_map <- setNames(meta$lineage, meta$sample_id)
V(g)$lineage <- meta_map[V(g)$name]
V(g)$lineage[is.na(V(g)$lineage) | V(g)$lineage == ""] <- "UNKNOWN"

cancer_type_assortativity <- tryCatch(
  assortativity_nominal(g, as.integer(factor(V(g)$lineage)), directed=FALSE),
  error=function(e) {
    warning("Could not compute cancer-type assortativity: ", conditionMessage(e))
    NA_real_
  }
)

# ---------------------------------------------------------------------------
# 3. Community detection
# ---------------------------------------------------------------------------
run_louvain <- function(graph, seed) {
  set.seed(seed)
  cluster_louvain(graph, weights=E(graph)$weight)
}

run_leiden <- function(graph, seed, resolution) {
  if (!"cluster_leiden" %in% getNamespaceExports("igraph")) {
    stop("igraph::cluster_leiden() is required for the Methods comparison")
  }
  set.seed(seed)
  tryCatch(
    cluster_leiden(graph, weights=E(graph)$weight, resolution=resolution),
    error=function(e_resolution) {
      tryCatch(
        cluster_leiden(graph, weights=E(graph)$weight,
                       resolution_parameter=resolution),
        error=function(e_resolution_parameter) {
          stop(
            "cluster_leiden() failed with both resolution and ",
            "resolution_parameter arguments. First error: ",
            conditionMessage(e_resolution), "; second error: ",
            conditionMessage(e_resolution_parameter)
          )
        }
      )
    }
  )
}

cat("[4] Running weighted Louvain (seed=", opt$seed, ")...\n", sep="")
cl_louvain <- run_louvain(g, opt$seed)
cat("  Louvain communities:", length(unique(membership(cl_louvain))), "\n")
cat("  Louvain weighted modularity:", round(modularity(cl_louvain), 4), "\n")

if (opt[["expected-communities"]] > 0 &&
    length(unique(membership(cl_louvain))) != opt[["expected-communities"]]) {
  warning("Expected ", opt[["expected-communities"]],
          " Louvain communities but found ",
          length(unique(membership(cl_louvain))))
}

cat("[5] Running weighted Leiden (seed=", opt$seed,
    ", resolution=", opt[["leiden-resolution"]], ")...\n", sep="")
cl_leiden <- run_leiden(g, opt$seed, opt[["leiden-resolution"]])
cat("  Leiden communities:", length(unique(membership(cl_leiden))), "\n")

# ---------------------------------------------------------------------------
# 4. Derived node and community summaries
# ---------------------------------------------------------------------------
make_node_table <- function(algorithm, membership_vector) {
  mem <- as.integer(membership_vector[V(g)$name])
  dt <- data.table(
    sample = V(g)$name,
    component = mem,
    lineage = V(g)$lineage,
    degree = as.integer(degree(g)[V(g)$name]),
    weighted_degree = as.numeric(strength(g, weights=E(g)$weight)[V(g)$name]),
    algorithm = algorithm
  )
  dt[, community_size := as.integer(.N), by=component]
  setcolorder(dt, c(
    "sample", "component", "community_size", "lineage", "degree",
    "weighted_degree", "algorithm"
  ))
  dt[order(component, sample)]
}

make_metrics <- function(node_dt) {
  alg <- unique(node_dt$algorithm)
  if (length(alg) != 1L) {
    stop("Expected one algorithm label per node table")
  }
  comp_map <- setNames(node_dt$component, node_dt$sample)
  edge_dt <- copy(E_und)
  edge_dt[, `:=`(
    component_from = as.integer(comp_map[from]),
    component_to = as.integer(comp_map[to])
  )]
  within <- edge_dt[component_from == component_to & !is.na(component_from)]
  edge_metrics <- within[, .(
    within_community_edge_count = .N,
    mean_within_community_edge_weight = mean(weight),
    total_within_community_edge_weight = sum(weight)
  ), by=.(component=component_from)]

  lineage_counts <- node_dt[, .N, by=.(component, lineage)]
  setorder(lineage_counts, component, -N, lineage)
  dominant <- lineage_counts[, .SD[1], by=component]
  setnames(dominant, c("lineage", "N"), c("dominant_lineage", "dominant_lineage_count"))
  lineage_strings <- lineage_counts[
    , .(lineage_counts=paste(paste0(lineage, "=", N), collapse=";")),
    by=component
  ]

  metrics <- node_dt[, .(
    algorithm = alg,
    community_size = .N,
    mean_degree = as.numeric(mean(degree)),
    median_degree = as.numeric(median(degree)),
    mean_weighted_degree = as.numeric(mean(weighted_degree)),
    median_weighted_degree = as.numeric(median(weighted_degree))
  ), by=component]
  metrics <- merge(metrics, edge_metrics, by="component", all.x=TRUE, sort=FALSE)
  metrics <- merge(metrics, dominant, by="component", all.x=TRUE, sort=FALSE)
  metrics <- merge(metrics, lineage_strings, by="component", all.x=TRUE, sort=FALSE)
  metrics[is.na(within_community_edge_count), within_community_edge_count := 0L]
  metrics[, community_purity := dominant_lineage_count / community_size]
  setcolorder(metrics, c(
    "algorithm", "component", "community_size",
    "within_community_edge_count", "mean_within_community_edge_weight",
    "total_within_community_edge_weight", "mean_degree", "median_degree",
    "mean_weighted_degree", "median_weighted_degree", "dominant_lineage",
    "dominant_lineage_count", "community_purity", "lineage_counts"
  ))
  metrics[order(algorithm, component)]
}

make_discordant <- function(node_dt, metrics_dt) {
  dom <- metrics_dt[, .(
    algorithm, component, dominant_lineage, community_purity
  )]
  out <- merge(node_dt, dom, by=c("algorithm", "component"), all.x=TRUE, sort=FALSE)
  out <- out[lineage != dominant_lineage]
  setcolorder(out, c(
    "algorithm", "sample", "component", "lineage", "dominant_lineage",
    "degree", "weighted_degree", "community_size", "community_purity"
  ))
  out[order(algorithm, component, lineage, sample)]
}

make_sweep_assignment <- function(resolution, membership_vector) {
  node_dt <- make_node_table("Leiden", membership_vector)
  setnames(node_dt, c("sample", "component", "lineage"),
           c("profile_id", "leiden_community", "cancer_type"))
  node_dt[, cell_line := derive_cell_line(profile_id)]

  counts <- node_dt[, .N, by=.(leiden_community, cancer_type)]
  setorder(counts, leiden_community, -N, cancer_type)
  majority <- counts[, .SD[1], by=leiden_community]
  setnames(
    majority,
    c("cancer_type", "N"),
    c("community_majority_cancer_type", "community_majority_cancer_type_count")
  )
  count_strings <- counts[
    , .(community_cancer_type_counts=paste(paste0(cancer_type, "=", N), collapse=";")),
    by=leiden_community
  ]

  out <- merge(node_dt, majority, by="leiden_community", all.x=TRUE, sort=FALSE)
  out <- merge(out, count_strings, by="leiden_community", all.x=TRUE, sort=FALSE)
  out[, community_purity := community_majority_cancer_type_count / community_size]
  out[, is_singleton := community_size == 1L]
  out[, resolution := resolution]
  setcolorder(out, c(
    "resolution", "profile_id", "cell_line", "cancer_type",
    "leiden_community", "community_size", "community_cancer_type_counts",
    "community_majority_cancer_type", "community_purity", "is_singleton",
    "degree", "weighted_degree"
  ))
  out[order(resolution, leiden_community, cancer_type, profile_id)]
}

classify_nbl_rbl <- function(mixed_non_singleton_any,
                             nbl_majority_non_singleton_count,
                             rbl_majority_non_singleton_count,
                             singleton_fraction_all,
                             singleton_fraction_nbl_rbl) {
  if (mixed_non_singleton_any) {
    return("mixed")
  }
  if (singleton_fraction_all >= 0.5 || singleton_fraction_nbl_rbl >= 0.5) {
    return("singleton_degenerate")
  }
  if (nbl_majority_non_singleton_count > 0L &&
      rbl_majority_non_singleton_count > 0L) {
    return("separated_non_degenerate")
  }
  "ambiguous"
}

make_sweep_summary <- function(resolution, assignment_dt, membership_vector) {
  community_dt <- unique(assignment_dt[, .(
    leiden_community, community_size, community_majority_cancer_type,
    community_purity, is_singleton
  )])

  ct_counts <- assignment_dt[
    cancer_type %in% c("NBL", "RBL"),
    .N,
    by=.(leiden_community, cancer_type)
  ]
  if (nrow(ct_counts) == 0L) {
    ct_wide <- data.table(leiden_community=integer(), NBL=integer(), RBL=integer())
  } else {
    ct_wide <- dcast(
      ct_counts,
      leiden_community ~ cancer_type,
      value.var="N",
      fill=0
    )
    if (!"NBL" %in% names(ct_wide)) ct_wide[, NBL := 0L]
    if (!"RBL" %in% names(ct_wide)) ct_wide[, RBL := 0L]
  }
  ct_wide <- merge(
    community_dt[, .(leiden_community, community_size, is_singleton)],
    ct_wide,
    by="leiden_community",
    all.x=TRUE,
    sort=FALSE
  )
  ct_wide[is.na(NBL), NBL := 0L]
  ct_wide[is.na(RBL), RBL := 0L]

  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  n_communities <- nrow(community_dt)
  n_singletons <- community_dt[is_singleton == TRUE, .N]
  nbl_total <- assignment_dt[cancer_type == "NBL", .N]
  rbl_total <- assignment_dt[cancer_type == "RBL", .N]
  nbl_singletons <- assignment_dt[cancer_type == "NBL" & is_singleton == TRUE, .N]
  rbl_singletons <- assignment_dt[cancer_type == "RBL" & is_singleton == TRUE, .N]

  n_mixed <- ct_wide[NBL > 0L & RBL > 0L & community_size > 1L, .N]
  same_any <- ct_wide[NBL > 0L & RBL > 0L, .N] > 0L
  same_non_singleton <- n_mixed > 0L

  nbl_majority_count <- community_dt[community_majority_cancer_type == "NBL", .N]
  rbl_majority_count <- community_dt[community_majority_cancer_type == "RBL", .N]
  nbl_majority_non_singleton <- community_dt[
    community_majority_cancer_type == "NBL" & is_singleton == FALSE,
    .N
  ]
  rbl_majority_non_singleton <- community_dt[
    community_majority_cancer_type == "RBL" & is_singleton == FALSE,
    .N
  ]

  singleton_fraction_all <- n_singletons / n_nodes
  target_total <- nbl_total + rbl_total
  target_singletons <- nbl_singletons + rbl_singletons
  singleton_fraction_nbl_rbl <- if (target_total > 0L) target_singletons / target_total else NA_real_

  status <- classify_nbl_rbl(
    same_non_singleton,
    nbl_majority_non_singleton,
    rbl_majority_non_singleton,
    singleton_fraction_all,
    singleton_fraction_nbl_rbl
  )

  weighted_modularity <- tryCatch(
    modularity(g, as.integer(membership_vector[V(g)$name]), weights=E(g)$weight),
    error=function(e) {
      warning("Could not compute weighted modularity at resolution ",
              resolution, ": ", conditionMessage(e))
      NA_real_
    }
  )

  data.table(
    resolution = resolution,
    n_nodes = n_nodes,
    n_edges = n_edges,
    n_communities = n_communities,
    n_singletons = n_singletons,
    singleton_fraction = singleton_fraction_all,
    weighted_modularity = as.numeric(weighted_modularity),
    cancer_type_assortativity = as.numeric(cancer_type_assortativity),
    n_mixed_NBL_RBL_communities = n_mixed,
    NBL_RBL_same_community_any = same_any,
    NBL_RBL_same_non_singleton_community_any = same_non_singleton,
    NBL_profiles_total = nbl_total,
    RBL_profiles_total = rbl_total,
    NBL_profiles_in_singletons = nbl_singletons,
    RBL_profiles_in_singletons = rbl_singletons,
    NBL_majority_community_count = nbl_majority_count,
    RBL_majority_community_count = rbl_majority_count,
    NBL_majority_non_singleton_community_count = nbl_majority_non_singleton,
    RBL_majority_non_singleton_community_count = rbl_majority_non_singleton,
    NBL_RBL_separation_status = status
  )
}

plot_sweep_summary <- function(summary_dt, path) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(NULL))
  }
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  status_levels <- c(
    "mixed", "separated_non_degenerate", "singleton_degenerate", "ambiguous"
  )
  status_cols <- c(
    mixed = "#b2182b",
    separated_non_degenerate = "#2166ac",
    singleton_degenerate = "#666666",
    ambiguous = "#f4a582"
  )
  status_idx <- match(summary_dt$NBL_RBL_separation_status, status_levels)
  xlog <- all(summary_dt$resolution > 0)

  pdf(path, width=8, height=6)
  op <- par(no.readonly=TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add=TRUE)
  par(mfrow=c(2, 1), mar=c(4, 4.5, 3, 1))
  plot(
    summary_dt$resolution,
    summary_dt$n_communities,
    type="b",
    pch=19,
    log=if (xlog) "x" else "",
    xlab="Leiden resolution",
    ylab="Count",
    main="Leiden resolution sweep",
    col="#2166ac",
    ylim=range(c(summary_dt$n_communities, summary_dt$n_singletons), na.rm=TRUE)
  )
  lines(summary_dt$resolution, summary_dt$n_singletons, type="b", pch=17, col="#b2182b")
  legend(
    "topleft",
    legend=c("Communities", "Singletons"),
    col=c("#2166ac", "#b2182b"),
    pch=c(19, 17),
    bty="n"
  )

  plot(
    summary_dt$resolution,
    status_idx,
    pch=19,
    log=if (xlog) "x" else "",
    xlab="Leiden resolution",
    ylab="NBL/RBL status",
    yaxt="n",
    ylim=c(0.5, length(status_levels) + 0.5),
    col=status_cols[summary_dt$NBL_RBL_separation_status]
  )
  axis(2, at=seq_along(status_levels), labels=status_levels, las=2, cex.axis=0.75)
  grid()
  invisible(NULL)
}

write_table_if_missing <- function(dt, path, label) {
  if (file.exists(path) && file.info(path)$size > 0) {
    cat("[KEEP]", label, "already exists; not overwriting:", path, "\n")
    return(invisible(FALSE))
  }
  fwrite(dt, path, sep="\t")
  cat("[OK]", label, "written to:", path, "\n")
  invisible(TRUE)
}

louvain_dt <- make_node_table("Louvain", membership(cl_louvain))
leiden_dt <- make_node_table("Leiden", membership(cl_leiden))
metrics <- rbindlist(list(make_metrics(louvain_dt), make_metrics(leiden_dt)),
                     use.names=TRUE)
discordant <- rbindlist(
  list(make_discordant(louvain_dt, metrics),
       make_discordant(leiden_dt, metrics)),
  use.names=TRUE
)

cat("[6] Running Leiden resolution sweep:",
    paste(format_numeric(resolution_sweep), collapse=", "), "\n")
sweep_assignments <- vector("list", length(resolution_sweep))
sweep_summary <- vector("list", length(resolution_sweep))
for (i in seq_along(resolution_sweep)) {
  res <- resolution_sweep[[i]]
  cl_res <- run_leiden(g, opt$seed, res)
  mem <- membership(cl_res)
  assign_dt <- make_sweep_assignment(res, mem)
  summary_dt <- make_sweep_summary(res, assign_dt, mem)
  sweep_assignments[[i]] <- assign_dt
  sweep_summary[[i]] <- summary_dt

  cat(
    "  resolution=", format_numeric(res),
    " communities=", summary_dt$n_communities,
    " singletons=", summary_dt$n_singletons,
    " NBL_RBL_same=", summary_dt$NBL_RBL_same_community_any,
    " same_non_singleton=", summary_dt$NBL_RBL_same_non_singleton_community_any,
    " status=", summary_dt$NBL_RBL_separation_status,
    "\n",
    sep=""
  )
  if (summary_dt$singleton_fraction >= 0.5) {
    warning(
      "Leiden resolution ", format_numeric(res),
      " is mostly singleton communities (",
      round(100 * summary_dt$singleton_fraction, 1), "%)"
    )
  }
}
sweep_assignments <- rbindlist(sweep_assignments, use.names=TRUE)
sweep_summary <- rbindlist(sweep_summary, use.names=TRUE)

# ---------------------------------------------------------------------------
# 5. Write outputs
# ---------------------------------------------------------------------------
write_table_if_missing(louvain_dt, opt$out, "Louvain communities")
write_table_if_missing(leiden_dt, opt[["leiden-out"]], "Leiden communities")
write_table_if_missing(metrics, opt[["community-metrics-out"]], "Community metrics")
write_table_if_missing(
  discordant,
  opt[["lineage-discordant-out"]],
  "Lineage-discordant profiles"
)

if (nzchar(opt[["leiden-sweep-assignments-out"]])) {
  fwrite(sweep_assignments, opt[["leiden-sweep-assignments-out"]], sep="\t")
  cat("[OK] Leiden sweep assignments written to:",
      opt[["leiden-sweep-assignments-out"]], "\n")
}

if (nzchar(opt[["leiden-sweep-summary-out"]])) {
  fwrite(sweep_summary, opt[["leiden-sweep-summary-out"]], sep="\t")
  cat("[OK] Leiden sweep summary written to:",
      opt[["leiden-sweep-summary-out"]], "\n")
}

if (nzchar(opt[["leiden-sweep-plot-out"]])) {
  plot_sweep_summary(sweep_summary, opt[["leiden-sweep-plot-out"]])
  cat("[OK] Leiden sweep plot written to:",
      opt[["leiden-sweep-plot-out"]], "\n")
}

cat("  Louvain summary:\n")
print(louvain_dt[, .N, by=.(component, lineage)][order(component, lineage)])
cat("  Leiden summary:\n")
print(leiden_dt[, .N, by=.(component, lineage)][order(component, lineage)])
cat("  Leiden resolution sweep summary:\n")
print(sweep_summary[order(resolution)])
