#!/usr/bin/env Rscript
# Run Louvain community detection across resolution values for sensitivity analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option("--edges", type="character", default=NULL,
              help="[REQUIRED] Edge list TSV with from, to, weight"),
  make_option("--meta", type="character", default=NULL,
              help="[REQUIRED] Node metadata TSV with sample_id and cancer_type"),
  make_option("--feature-source", type="character", default=NULL,
              help="[REQUIRED] Feature source label"),
  make_option("--resolution-sweep", type="character", default=NULL,
              help="[REQUIRED] Comma-separated Louvain resolution values"),
  make_option("--seed", type="integer", default=1,
              help="Random seed for Louvain (default: 1)"),
  make_option("--assignments-out", type="character", default=NULL,
              help="[REQUIRED] Long-format assignment TSV"),
  make_option("--summary-out", type="character", default=NULL,
              help="[REQUIRED] Resolution-level summary TSV"),
  make_option("--plot-out", type="character", default=NULL,
              help="[REQUIRED] Diagnostic PDF")
)
opt <- parse_args(OptionParser(option_list=option_list))

required <- c(
  "edges", "meta", "feature-source", "resolution-sweep",
  "assignments-out", "summary-out", "plot-out"
)
missing <- required[vapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required argument(s): ", paste(missing, collapse=", "))
}
for (path in c(opt$edges, opt$meta)) {
  if (!file.exists(path)) stop("Input file not found: ", path)
}
for (path in c(opt[["assignments-out"]], opt[["summary-out"]], opt[["plot-out"]])) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
}

if (!"resolution" %in% names(formals(cluster_louvain))) {
  stop("Installed igraph::cluster_louvain() does not support a resolution argument")
}

parse_resolution_sweep <- function(raw) {
  parts <- trimws(unlist(strsplit(raw, "[,;[:space:]]+")))
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  if (length(vals) == 0L || any(!is.finite(vals))) {
    stop("Invalid --resolution-sweep values: ", raw)
  }
  vals <- unique(vals)
  vals[order(vals)]
}

derive_cell_line <- function(profile_id) {
  out <- sub("^NG-[^_]+_", "", profile_id)
  sub("_lib.*$", "", out)
}

format_numeric <- function(x) format(x, scientific=FALSE, trim=TRUE)

resolution_sweep <- parse_resolution_sweep(opt[["resolution-sweep"]])

edges <- fread(opt$edges)
if (!all(c("from", "to", "weight") %in% names(edges))) {
  stop("edges must contain columns: from, to, weight")
}
edge_dt <- edges[, .(weight=max(as.numeric(weight), na.rm=TRUE)),
                 by=.(from=pmin(from, to), to=pmax(from, to))]
setorder(edge_dt, from, to)

meta <- fread(opt$meta)
if (!"sample_id" %in% names(meta)) stop("metadata must contain sample_id")
if (!"cancer_type" %in% names(meta)) {
  if ("lineage" %in% names(meta)) {
    warning("Using legacy metadata column 'lineage' to populate 'cancer_type'")
    meta[, cancer_type := lineage]
  } else {
    stop("metadata must contain cancer_type")
  }
}
meta <- unique(meta, by="sample_id")

vertices <- data.table(name=meta$sample_id)
g <- graph_from_data_frame(edge_dt, directed=FALSE, vertices=vertices)
E(g)$weight <- edge_dt$weight
cancer_type_map <- setNames(as.character(meta$cancer_type), meta$sample_id)
V(g)$cancer_type <- cancer_type_map[V(g)$name]
V(g)$cancer_type[is.na(V(g)$cancer_type) | V(g)$cancer_type == ""] <- "UNKNOWN"

cancer_type_assortativity <- tryCatch(
  assortativity_nominal(g, as.integer(factor(V(g)$cancer_type)), directed=FALSE),
  error=function(e) {
    warning("Could not compute cancer-type assortativity: ", conditionMessage(e))
    NA_real_
  }
)

run_louvain <- function(graph, seed, resolution) {
  set.seed(seed)
  cluster_louvain(graph, weights=E(graph)$weight, resolution=resolution)
}

make_assignment <- function(resolution, membership_vector) {
  mem <- as.integer(membership_vector[V(g)$name])
  dt <- data.table(
    feature_source = opt[["feature-source"]],
    resolution = resolution,
    profile_id = V(g)$name,
    cell_line = derive_cell_line(V(g)$name),
    cancer_type = V(g)$cancer_type,
    louvain_community = mem,
    degree = as.integer(degree(g)[V(g)$name]),
    weighted_degree = as.numeric(strength(g, weights=E(g)$weight)[V(g)$name])
  )
  dt[, community_size := as.integer(.N), by=louvain_community]

  counts <- dt[, .N, by=.(louvain_community, cancer_type)]
  setorder(counts, louvain_community, -N, cancer_type)
  majority <- counts[, .SD[1], by=louvain_community]
  setnames(
    majority,
    c("cancer_type", "N"),
    c("community_majority_cancer_type", "community_majority_cancer_type_count")
  )
  count_strings <- counts[
    , .(community_cancer_type_counts=paste(paste0(cancer_type, "=", N), collapse=";")),
    by=louvain_community
  ]
  dt <- merge(dt, majority, by="louvain_community", all.x=TRUE, sort=FALSE)
  dt <- merge(dt, count_strings, by="louvain_community", all.x=TRUE, sort=FALSE)
  dt[, community_purity := community_majority_cancer_type_count / community_size]
  dt[, is_singleton := community_size == 1L]
  setcolorder(dt, c(
    "feature_source", "resolution", "profile_id", "cell_line", "cancer_type",
    "louvain_community", "community_size", "community_cancer_type_counts",
    "community_majority_cancer_type", "community_purity", "is_singleton",
    "degree", "weighted_degree"
  ))
  dt[order(resolution, louvain_community, cancer_type, profile_id)]
}

largest_majority <- function(community_dt, label) {
  out <- community_dt[community_majority_cancer_type == label]
  if (nrow(out) == 0L) {
    return(list(size=0L, composition=NA_character_))
  }
  setorder(out, -community_size, louvain_community)
  list(size=out$community_size[1], composition=out$community_cancer_type_counts[1])
}

classify_nbl_rbl <- function(mixed_non_singleton_any,
                             nbl_non_singleton_majority_count,
                             rbl_non_singleton_majority_count,
                             singleton_fraction_all,
                             singleton_fraction_nbl_rbl) {
  if (mixed_non_singleton_any) {
    return("mixed")
  }
  if (singleton_fraction_all >= 0.5 || singleton_fraction_nbl_rbl >= 0.5) {
    return("singleton_degenerate")
  }
  if (nbl_non_singleton_majority_count > 0L &&
      rbl_non_singleton_majority_count > 0L) {
    return("separated_non_degenerate")
  }
  "ambiguous"
}

compute_modularity <- function(membership_vector, resolution) {
  tryCatch({
    if ("resolution" %in% names(formals(modularity))) {
      modularity(g, as.integer(membership_vector[V(g)$name]),
                 weights=E(g)$weight, resolution=resolution)
    } else {
      modularity(g, as.integer(membership_vector[V(g)$name]),
                 weights=E(g)$weight)
    }
  }, error=function(e) {
    warning("Could not compute weighted modularity at resolution ",
            resolution, ": ", conditionMessage(e))
    NA_real_
  })
}

make_summary <- function(resolution, assignment_dt, membership_vector) {
  community_dt <- unique(assignment_dt[, .(
    louvain_community, community_size, community_cancer_type_counts,
    community_majority_cancer_type, community_purity, is_singleton
  )])

  ct_counts <- assignment_dt[
    cancer_type %in% c("NBL", "RBL"),
    .N,
    by=.(louvain_community, cancer_type)
  ]
  if (nrow(ct_counts) == 0L) {
    ct_wide <- data.table(louvain_community=integer(), NBL=integer(), RBL=integer())
  } else {
    ct_wide <- dcast(ct_counts, louvain_community ~ cancer_type,
                     value.var="N", fill=0)
    if (!"NBL" %in% names(ct_wide)) ct_wide[, NBL := 0L]
    if (!"RBL" %in% names(ct_wide)) ct_wide[, RBL := 0L]
  }
  ct_wide <- merge(
    community_dt[, .(louvain_community, community_size, is_singleton)],
    ct_wide,
    by="louvain_community",
    all.x=TRUE,
    sort=FALSE
  )
  ct_wide[is.na(NBL), NBL := 0L]
  ct_wide[is.na(RBL), RBL := 0L]

  n_nodes <- vcount(g)
  n_singletons <- community_dt[is_singleton == TRUE, .N]
  nbl_total <- assignment_dt[cancer_type == "NBL", .N]
  rbl_total <- assignment_dt[cancer_type == "RBL", .N]
  nbl_singletons <- assignment_dt[cancer_type == "NBL" & is_singleton == TRUE, .N]
  rbl_singletons <- assignment_dt[cancer_type == "RBL" & is_singleton == TRUE, .N]
  n_mixed <- ct_wide[NBL > 0L & RBL > 0L & community_size > 1L, .N]
  same_any <- ct_wide[NBL > 0L & RBL > 0L, .N] > 0L
  same_non_singleton <- n_mixed > 0L
  nbl_non_singleton_majority <- community_dt[
    community_majority_cancer_type == "NBL" & is_singleton == FALSE, .N
  ]
  rbl_non_singleton_majority <- community_dt[
    community_majority_cancer_type == "RBL" & is_singleton == FALSE, .N
  ]
  singleton_fraction <- n_singletons / n_nodes
  nbl_rbl_total <- nbl_total + rbl_total
  nbl_rbl_singletons <- nbl_singletons + rbl_singletons
  target_singleton_fraction <- if (nbl_rbl_total > 0L) nbl_rbl_singletons / nbl_rbl_total else NA_real_
  status <- classify_nbl_rbl(
    same_non_singleton,
    nbl_non_singleton_majority,
    rbl_non_singleton_majority,
    singleton_fraction,
    target_singleton_fraction
  )
  largest_nbl <- largest_majority(copy(community_dt), "NBL")
  largest_rbl <- largest_majority(copy(community_dt), "RBL")

  data.table(
    feature_source = opt[["feature-source"]],
    resolution = resolution,
    n_nodes = n_nodes,
    n_edges = ecount(g),
    n_communities = nrow(community_dt),
    n_singletons = n_singletons,
    singleton_fraction = singleton_fraction,
    weighted_modularity = as.numeric(compute_modularity(membership_vector, resolution)),
    cancer_type_assortativity = as.numeric(cancer_type_assortativity),
    n_mixed_NBL_RBL_communities = n_mixed,
    NBL_RBL_same_community_any = same_any,
    NBL_RBL_same_non_singleton_community_any = same_non_singleton,
    NBL_profiles_total = nbl_total,
    RBL_profiles_total = rbl_total,
    NBL_profiles_in_singletons = nbl_singletons,
    RBL_profiles_in_singletons = rbl_singletons,
    largest_NBL_majority_community_size = largest_nbl$size,
    largest_NBL_majority_community_composition = largest_nbl$composition,
    largest_RBL_majority_community_size = largest_rbl$size,
    largest_RBL_majority_community_composition = largest_rbl$composition,
    NBL_RBL_separation_status = status
  )
}

plot_summary <- function(summary_dt, path) {
  status_levels <- c("mixed", "separated_non_degenerate", "singleton_degenerate", "ambiguous")
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
    xlab="Louvain resolution",
    ylab="Count",
    main=paste("Louvain resolution sweep:", opt[["feature-source"]]),
    col="#2166ac",
    ylim=range(c(summary_dt$n_communities, summary_dt$n_singletons), na.rm=TRUE)
  )
  lines(summary_dt$resolution, summary_dt$n_singletons, type="b", pch=17, col="#b2182b")
  legend("topleft", legend=c("Communities", "Singletons"),
         col=c("#2166ac", "#b2182b"), pch=c(19, 17), bty="n")

  plot(
    summary_dt$resolution,
    status_idx,
    pch=19,
    log=if (xlog) "x" else "",
    xlab="Louvain resolution",
    ylab="NBL/RBL status",
    yaxt="n",
    ylim=c(0.5, length(status_levels) + 0.5),
    col=status_cols[summary_dt$NBL_RBL_separation_status]
  )
  axis(2, at=seq_along(status_levels), labels=status_levels, las=2, cex.axis=0.75)
  grid()
  invisible(NULL)
}

cat("[1] Loading graph for feature source:", opt[["feature-source"]], "\n")
cat("  Nodes:", vcount(g), " Edges:", ecount(g), "\n")
cat("[2] Running Louvain resolution sweep:",
    paste(format_numeric(resolution_sweep), collapse=", "), "\n")

assignments <- vector("list", length(resolution_sweep))
summaries <- vector("list", length(resolution_sweep))
for (i in seq_along(resolution_sweep)) {
  res <- resolution_sweep[[i]]
  cl <- run_louvain(g, opt$seed, res)
  mem <- membership(cl)
  assign_dt <- make_assignment(res, mem)
  summary_dt <- make_summary(res, assign_dt, mem)
  assignments[[i]] <- assign_dt
  summaries[[i]] <- summary_dt
  cat(
    "  feature_source=", opt[["feature-source"]],
    " resolution=", format_numeric(res),
    " communities=", summary_dt$n_communities,
    " singletons=", summary_dt$n_singletons,
    " singleton_fraction=", round(summary_dt$singleton_fraction, 4),
    " mixed_non_singleton=", summary_dt$NBL_RBL_same_non_singleton_community_any,
    " status=", summary_dt$NBL_RBL_separation_status,
    "\n",
    sep=""
  )
  if (summary_dt$singleton_fraction >= 0.5 ||
      ((summary_dt$NBL_profiles_in_singletons + summary_dt$RBL_profiles_in_singletons) /
       (summary_dt$NBL_profiles_total + summary_dt$RBL_profiles_total)) >= 0.5) {
    warning("Louvain resolution ", format_numeric(res),
            " is singleton-heavy or degenerate for feature source ",
            opt[["feature-source"]])
  }
}

assignments <- rbindlist(assignments, use.names=TRUE)
summaries <- rbindlist(summaries, use.names=TRUE)

fwrite(assignments, opt[["assignments-out"]], sep="\t")
fwrite(summaries, opt[["summary-out"]], sep="\t")
plot_summary(summaries, opt[["plot-out"]])

cat("[OK] Louvain sweep assignments written to:", opt[["assignments-out"]], "\n")
cat("[OK] Louvain sweep summary written to:", opt[["summary-out"]], "\n")
cat("[OK] Louvain sweep diagnostic plot written to:", opt[["plot-out"]], "\n")
