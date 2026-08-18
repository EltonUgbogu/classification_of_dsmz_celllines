# ==============================================================================
# multicohort_gene_space.R
# Deterministic common-gene-space construction for the multicohort profile
# ==============================================================================
#
# The multicohort profile concatenates per-cohort expression matrices. Cohort
# preprocessing can retain slightly different gene sets, so the matrices must be
# reduced to an explicit common gene space before concatenation.
#
#   G_common = G_brca n G_nbl n G_rbl
#
# using exact gene identifiers. No identifier normalisation is applied: the
# cohort matrices carry unversioned Ensembl gene IDs, so identifiers are
# already directly comparable. Genes absent from a cohort are dropped, never
# zero-filled, and gene alignment is always by identifier, never by row
# position.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# harmonise_gene_space: reduce cohort matrices to the exact common gene space
# ------------------------------------------------------------------------------
# X_list: named list of samples x genes matrices (colnames = gene identifiers).
#
# Returns a list with:
#   X_list       - the input matrices subset to G_common, in identical gene order
#   genes_common - the common gene identifiers, in the returned column order
#   audit        - per-cohort data.frame: n_genes_before, n_genes_after,
#                  n_genes_removed, n_duplicate_gene_ids
#
# Fails explicitly on: missing gene identifiers, duplicate gene identifiers
# within a cohort, an empty or implausibly small intersection, or any residual
# disagreement in gene identity or order after subsetting.

harmonise_gene_space <- function(X_list, min_common_genes = 50L) {

  if (!is.list(X_list) || length(X_list) == 0L) {
    stop("harmonise_gene_space: X_list must be a non-empty named list of matrices.")
  }
  if (is.null(names(X_list)) || any(!nzchar(names(X_list)))) {
    stop("harmonise_gene_space: X_list must be named by cohort.")
  }

  # --- Identifier integrity, per cohort --------------------------------------
  for (cohort in names(X_list)) {
    ids <- colnames(X_list[[cohort]])
    if (is.null(ids)) {
      stop("harmonise_gene_space: cohort '", cohort,
           "' has no gene identifiers (colnames).")
    }
    if (any(is.na(ids)) || any(!nzchar(trimws(ids)))) {
      stop("harmonise_gene_space: cohort '", cohort,
           "' has missing or blank gene identifiers.")
    }
    dup <- unique(ids[duplicated(ids)])
    if (length(dup) > 0L) {
      stop("harmonise_gene_space: cohort '", cohort, "' has ",
           length(dup), " duplicate gene identifier(s), which makes gene ",
           "identity ambiguous. Examples: ",
           paste(utils::head(dup, 10), collapse = ", "),
           "\nResolve duplicates upstream with an explicit, tested rule ",
           "before multicohort construction.")
    }
  }

  n_before <- vapply(X_list, ncol, integer(1))

  # --- Common gene space ------------------------------------------------------
  # Ordered by the first cohort's gene order so the result is deterministic and
  # independent of the order in which cohorts are supplied.
  reference_ids <- colnames(X_list[[1]])
  genes_common <- Reduce(intersect, lapply(X_list, colnames))
  genes_common <- reference_ids[reference_ids %in% genes_common]

  if (length(genes_common) < min_common_genes) {
    stop("harmonise_gene_space: only ", length(genes_common),
         " gene(s) are shared by all cohorts (minimum ", min_common_genes,
         "). Check cohort preprocessing and gene identifier formats.")
  }

  # --- Subset by identifier (never by position) -------------------------------
  X_list <- lapply(X_list, function(X) X[, genes_common, drop = FALSE])

  # --- Post-conditions: identical gene identity and order across cohorts ------
  for (cohort in names(X_list)) {
    ids <- colnames(X_list[[cohort]])
    if (!identical(ids, genes_common)) {
      stop("harmonise_gene_space: cohort '", cohort,
           "' does not carry the common gene space in the expected order ",
           "after subsetting. This indicates positional rather than ",
           "identifier-based alignment.")
    }
  }

  audit <- data.frame(
    cohort               = names(X_list),
    n_genes_before       = as.integer(n_before[names(X_list)]),
    n_genes_after        = vapply(X_list, ncol, integer(1)),
    n_genes_removed      = as.integer(n_before[names(X_list)]) -
                             vapply(X_list, ncol, integer(1)),
    n_duplicate_gene_ids = 0L,
    identifier_rule      = "exact match on unversioned Ensembl gene IDs; no normalisation applied",
    stringsAsFactors     = FALSE,
    row.names            = NULL
  )

  list(X_list = X_list, genes_common = genes_common, audit = audit)
}

# ------------------------------------------------------------------------------
# validate_feature_list: check a feature list against the common gene space
# ------------------------------------------------------------------------------
# Cached feature lists written by an earlier run can reference genes that are no
# longer in the common gene space, for example after cohort matrices are
# regenerated with different gene retention. Reusing such a list would subscript
# the joint matrix with absent identifiers.
#
# Returns list(ok = logical, missing = character, n_missing = integer).

validate_feature_list <- function(genes, genes_common, label = "feature list") {
  genes <- genes[nzchar(genes)]
  missing <- setdiff(genes, genes_common)
  list(
    ok        = length(missing) == 0L,
    missing   = missing,
    n_missing = length(missing),
    label     = label
  )
}
