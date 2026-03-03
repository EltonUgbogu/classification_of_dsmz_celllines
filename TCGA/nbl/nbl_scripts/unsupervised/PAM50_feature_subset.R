#PAM50_feature_subset.R
# ------------------------------------------------------------------------------
# Marker Sets: Common breast cancer subtype and proliferation-related genes
# ------------------------------------------------------------------------------
marker_symbols <- c(
  "ESR1", "PGR", "ERBB2", "MKI67", "KRT5", "KRT14", "KRT17", "FOXC1", "GRB7",
  "BCL2", "EGFR", "CCNB1", "BIRC5", "MYBL2", "KIF2C", "KRT8", "KRT18", "GATA3", "XBP1"
)

# PAM50 classifier gene symbols (used in intrinsic subtype assignment)
pam50_symbols <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC6","CDC20",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT5","KRT14","KRT17","KRT23","MDM2","MELK",
  "MIA","MKI67","MLPH","MMP11","MYBL2","NAT1","ORC6L","PGR","PHGDH","PTTG1",
  "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T","WHSC1L1","KRT7","KRT15"
)



#' @title Map Ensembl Gene IDs to HGNC Symbols (Robust, Multi-Source)
#'
#' @description
#' Robustly converts Ensembl gene IDs (with or without version) to HGNC symbols
#' using up to three sources in order of preference:
#' \enumerate{
#'   \item **Bioconductor** (`org.Hs.eg.db`)
#'   \item **User-supplied mapping** (e.g., DSMZ annotation table)
#'   \item **`rowData()` of a `SummarizedExperiment`**
#' }
#' Returns a named character vector: `Ensembl ID → Symbol`.  
#' Designed for use with `rownames(V) <- map[strip_ensver(rownames(V))]`  
#'
#' @param expr_mat Matrix or `SummarizedExperiment` with Ensembl IDs in `rownames`.
#' @param user_map Optional `data.frame` with columns `ensembl` (or `Ensembl_ID`) and `symbol` (or `gene_name`).
#'                 If `NULL`, skipped.
#' @param se_obj Optional `SummarizedExperiment`. If provided, `rowData()` is mined for symbol columns.
#' @param strip_version Logical. Remove `.1`, `.2` version suffix from Ensembl IDs? Default `TRUE`.
#' @param quiet Logical. Suppress progress messages? Default `FALSE`.
#'
#' @return Named character vector: `Ensembl (no version) → HGNC symbol`.
#'         `NA` if no mapping found for a given ID.
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   symbol_vec <- map_ensembl_to_symbol(V_adj)
#'   rownames(V_adj) <- symbol_vec[strip_ensver(rownames(V_adj))]
#'
#'   # With DSMZ annotation
#'   symbol_vec <- map_ensembl_to_symbol(V_adj, user_map = dsmz_raw)
#' }
#' @export
map_ensembl_to_symbol <- function(
  expr_mat,
  user_map = NULL,
  se_obj = NULL,
  strip_version = TRUE,
  quiet = FALSE
) {
  # ------------------------------------------------------------------ #
  # 0. Helper: strip Ensembl version
  # ------------------------------------------------------------------ #
  strip_ensver <- function(x) {
    if (strip_version) sub("\\.\\d+$", "", x) else x
  }
  .msg <- function(...) if (!quiet) message(...)

  # ------------------------------------------------------------------ #
  # 1. Extract Ensembl IDs from expression matrix
  # ------------------------------------------------------------------ #
  if (inherits(expr_mat, "SummarizedExperiment")) {
    se_obj <- expr_mat
    expr_mat <- SummarizedExperiment::assay(expr_mat)
  }
  stopifnot(is.matrix(expr_mat) || is(expr_mat, "dgCMatrix"))

  ensembl_ids <- rownames(expr_mat)
  if (is.null(ensembl_ids)) {
    stop("Expression matrix must have rownames (Ensembl IDs).")
  }
  ensembl_clean <- strip_ensver(ensembl_ids)
  unique_ensembl <- unique(ensembl_clean)

  .msg("[INFO] Found ", length(unique_ensembl), " unique Ensembl IDs in expression matrix.")

  # ------------------------------------------------------------------ #
  # 2. Source 1: Bioconductor org.Hs.eg.db (without attaching)
  # ------------------------------------------------------------------ #
  map1 <- NULL
  if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    .msg("[INFO] Attempting mapping via org.Hs.eg.db...")
    orgdb <- tryCatch(
      getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
      error = function(e) NULL
    )
    if (!is.null(orgdb)) {
      map1 <- tryCatch({
        AnnotationDbi::select(
          orgdb,
          keys = unique_ensembl,
          keytype = "ENSEMBL",
          columns = "SYMBOL"
        ) %>%
          dplyr::filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
          dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
          dplyr::transmute(ensembl = ENSEMBL, symbol = SYMBOL)
      }, error = function(e) {
        .msg("[WARN] org.Hs.eg.db failed: ", e$message); NULL
      })
      if (!is.null(map1)) .msg("[INFO] org.Hs.eg.db: ", nrow(map1), " mappings")
    } else {
      .msg("[WARN] Could not retrieve org.Hs.eg.db object.")
    }
  } else {
    .msg("[WARN] org.Hs.eg.db not available.")
  }

  # ------------------------------------------------------------------ #
  # 3. Source 2: User-provided mapping
  # ------------------------------------------------------------------ #
  map2 <- NULL
  if (!is.null(user_map)) {
    .msg("[INFO] Using user-supplied mapping...")
    if (!is.data.frame(user_map)) stop("`user_map` must be a data.frame.")
    # Detect column names
    ens_col <- intersect(names(user_map), c("ensembl", "Ensembl_ID", "ensembl_id"))
    sym_col <- intersect(names(user_map), c("symbol", "gene_name", "hgnc_symbol", "SYMBOL"))
    if (length(ens_col) == 0 || length(sym_col) == 0) {
      stop("`user_map` must contain Ensembl and symbol columns. Found: ", paste(names(user_map), collapse = ", "))
    }
    map2 <- user_map %>%
      dplyr::transmute(
        ensembl = strip_ensver(.data[[ens_col[1]]]),
        symbol  = as.character(.data[[sym_col[1]]])
      ) %>%
      dplyr::filter(!is.na(symbol), nzchar(symbol)) %>%
      dplyr::distinct(symbol, .keep_all = TRUE) %>%
      dplyr::distinct(ensembl, .keep_all = TRUE)
    .msg("[INFO] User map: ", nrow(map2), " mappings")
  }

  # ------------------------------------------------------------------ #
  # 4. Source 3: SummarizedExperiment rowData
  # ------------------------------------------------------------------ #
  map3 <- NULL
  if (!is.null(se_obj) && inherits(se_obj, "SummarizedExperiment")) {
    .msg("[INFO] Extracting mapping from rowData()...")
    rd <- as.data.frame(SummarizedExperiment::rowData(se_obj))
    sym_candidates <- intersect(tolower(names(rd)), c("gene_name", "symbol", "hgnc_symbol"))
    if (length(sym_candidates) > 0) {
      sym_col <- names(rd)[match(sym_candidates[1], tolower(names(rd)))]
      map3 <- data.frame(
        ensembl = strip_ensver(rownames(rd)),
        symbol  = as.character(rd[[sym_col]]),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(!is.na(symbol), symbol != "") %>%
        dplyr::distinct(symbol, .keep_all = TRUE) %>%
        dplyr::distinct(ensembl, .keep_all = TRUE)
      .msg("[INFO] rowData: ", nrow(map3), " mappings (col: '", sym_col, "')")
    } else {
      .msg("[WARN] No symbol column in rowData.")
    }
  }

  # ------------------------------------------------------------------ #
  # 5. Combine all sources + deduplicate
  # ------------------------------------------------------------------ #
  combined <- dplyr::bind_rows(map1, map2, map3)

  if (is.null(combined) || nrow(combined) == 0) {
    warning("No gene symbol mappings found from any source.")
    mapping_vec <- setNames(rep(NA_character_, length(ensembl_ids)), ensembl_ids)
    return(mapping_vec)
  }

  # Final deduplication: one symbol per Ensembl ID (prefer earlier source)
  final_map <- combined %>%
    dplyr::distinct(ensembl, .keep_all = TRUE)

  # Expand to full input set (preserve original versioned names)
  idx <- match(ensembl_clean, final_map$ensembl)
  mapping_vec <- final_map$symbol[idx]
  names(mapping_vec) <- ensembl_ids

  if (!quiet) {
    mapped_ok <- sum(!is.na(mapping_vec))
    message("[INFO] Final mapping: ", mapped_ok, "/", length(mapping_vec),
            " genes mapped (", round(100 * mapped_ok / length(mapping_vec), 1), "%)")
  }

  return(mapping_vec)
}

# 1. Basic: only Bioconductor
#symbol_vec <- map_ensembl_to_symbol(V_adj)

# 2. With DSMZ annotation
#symbol_vec <- map_ensembl_to_symbol(V_adj, user_map = dsmz_raw)

# 3. With TCGA SummarizedExperiment
#symbol_vec <- map_ensembl_to_symbol(V_adj, se_obj = tcga_se)

# 4. All three
#symbol_vec <- map_ensembl_to_symbol(V_adj, user_map = dsmz_raw, se_obj = tcga_se)

# Apply
#rownames(V_adj) <- symbol_vec

#' @title Subset Expression Matrix to PAM50 Genes with Collapsing
#'
#' @description
#' Subsets a gene expression matrix to the **PAM50 gene set**, using a pre-computed
#' Ensembl → HGNC symbol mapping. Handles:
#' - Versioned Ensembl IDs
#' - Multiple Ensembl IDs per symbol
#' - Case-insensitive matching
#' - Collapsing duplicates via `sum`, `mean`, or `max`
#'
#' @param mat          Expression matrix (genes × samples). Rownames = Ensembl IDs.
#' @param ens2sym      Data frame with columns `ensembl` and `symbol`, or named vector.
#' @param pam50_genes  Character vector of official PAM50 gene symbols (case-insensitive).
#' @param collapse     How to combine duplicate genes: `"sum"`, `"mean"`, or `"max"` (default `"sum"`).
#' @param quiet        Suppress progress messages? Default `FALSE`.
#'
#' @return Matrix with **PAM50 genes as rows** (uppercase symbols), samples as columns.
#' @export
#'
#' @examples
#' \dontrun{
#'   pam50_mat <- subset_to_pam50(
#'     mat = vst_counts,
#'     ens2sym = symbol_map,
#'     pam50_genes = pam50_official
#'   )
#' }
subset_to_pam50 <- function(
  mat,
  ens2sym,
  pam50_genes,
  collapse = c("sum", "mean", "max"),
  quiet = FALSE
) {
  collapse <- match.arg(collapse)

  # ------------------------------------------------------------------ #
  # 1. Input validation
  # ------------------------------------------------------------------ #
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (is.null(rownames(mat))) stop("Matrix must have rownames (Ensembl IDs).")

  # Accept named vector or data frame for ens2sym
  if (is.data.frame(ens2sym)) {
    if (!all(c("ensembl", "symbol") %in% colnames(ens2sym))) {
      stop("ens2sym must have columns 'ensembl' and 'symbol'.")
    }
    symbol_lookup <- setNames(ens2sym$symbol, ens2sym$ensembl)
  } else if (is.vector(ens2sym) && !is.null(names(ens2sym))) {
    symbol_lookup <- ens2sym
  } else {
    stop("ens2sym must be a data frame or named vector.")
  }

  # ------------------------------------------------------------------ #
  # 2. Strip version & match Ensembl IDs
  # ------------------------------------------------------------------ #
  ensembl_clean <- sub("\\.\\d+$", "", rownames(mat))
  mapped_idx <- which(ensembl_clean %in% names(symbol_lookup))
  if (length(mapped_idx) == 0) {
    stop("No Ensembl IDs in matrix found in ens2sym mapping.")
  }
  if (!quiet) message("[INFO] ", length(mapped_idx), " / ", nrow(mat), " genes mapped to symbols.")

  # ------------------------------------------------------------------ #
  # 3. Match to PAM50 (case-insensitive) with alias normalization
  # ------------------------------------------------------------------ #
  symbols <- symbol_lookup[ensembl_clean[mapped_idx]]

  # Handle legacy aliases so we reach 50/50 consistently
  alias_map <- c("ORC6L" = "ORC6", "WHSC1L1" = "NSD3")
  sym_up <- toupper(symbols)
  sym_up <- ifelse(sym_up %in% names(alias_map), alias_map[sym_up], sym_up)

  pam50_up <- toupper(pam50_genes)
  pam50_up <- ifelse(pam50_up %in% names(alias_map), alias_map[pam50_up], pam50_up)

  is_pam50 <- !is.na(sym_up) & sym_up %in% pam50_up

  pam50_idx <- mapped_idx[is_pam50]
  if (length(pam50_idx) == 0) {
    stop("None of the PAM50 genes were found in the data.")
  }

  if (!quiet) message("[INFO] ", length(pam50_idx), " rows match PAM50 genes (before collapsing).")

  # ------------------------------------------------------------------ #
  # 4. Subset matrix
  # ------------------------------------------------------------------ #
  sub_mat <- mat[pam50_idx, , drop = FALSE]
  pam50_symbols <- sym_up[is_pam50]

  # ------------------------------------------------------------------ #
  # 5. Collapse duplicates (one row per symbol)
  # ------------------------------------------------------------------ #
  split_by_gene <- split(seq_len(nrow(sub_mat)), pam50_symbols)

  collapse_row <- function(rows) {
    row_data <- sub_mat[rows, , drop = FALSE]
    switch(collapse,
           sum  = colSums(row_data,  na.rm = TRUE),
           mean = colMeans(row_data, na.rm = TRUE),
           max  = apply(row_data, 2, max, na.rm = TRUE)
    )
  }

  collapsed_rows <- lapply(split_by_gene, collapse_row)
  result <- do.call(rbind, collapsed_rows)
  rownames(result) <- toupper(names(split_by_gene))

  # ------------------------------------------------------------------ #
  # 6. Final report
  # ------------------------------------------------------------------ #
  missing_pam50 <- setdiff(toupper(pam50_genes), rownames(result))
  if (length(missing_pam50) > 0 && !quiet) {
    message("[WARN] ", length(missing_pam50), " PAM50 genes missing: ",
            paste(head(missing_pam50, 5), collapse = ", "),
            if (length(missing_pam50) > 5) "..." else "")
  }

  if (!quiet) {
    message("[SUCCESS] Final PAM50 matrix: ", nrow(result), " genes × ", ncol(result), " samples")
  }

  return(result)
}


# Assume you have:
# - V_adj: your expression matrix (Ensembl rownames)
# - ens2sym: from map_ensembl_to_symbol(V_adj)
# - pam50_genes: official list (e.g., from pamr package)

#pam50_matrix <- subset_to_pam50(
#  mat = V_adj,
#  ens2sym = ens2sym,
#  pam50_genes = pam50_official,
#  collapse = "sum"
#)

#' @title Check for Missing Subtype Marker Genes in the Expression Matrix
#'
#' @description
#' This function verifies whether all expected marker genes are present
#' in the expression matrix after Ensembl-to-symbol mapping. It is commonly
#' used prior to supervised classification or subtype annotation pipelines.
#'
#' @param feature_count_matrix Matrix with genes as rows and samples as columns.
#' @param ens2sym_union Data frame with columns `ensembl` and `symbol` mapping Ensembl IDs to gene symbols.
#' @param marker_symbols Character vector of expected gene symbols (e.g., for subtype classifiers).
#'
#' @return None; prints diagnostic messages and missing symbols if any.
#'
missing_subtype_marker <- function(feature_count_matrix, ens2sym_union, marker_symbols) {
  stopifnot(all(c("ensembl", "symbol") %in% colnames(ens2sym_union)))

  # Identify marker Ensembl IDs that map to the expected symbols
  marker_ens <- ens2sym_union$ensembl[
    tolower(ens2sym_union$symbol) %in% tolower(marker_symbols)
  ]
  marker_ens <- intersect(marker_ens, rownames(feature_count_matrix))

  if (length(marker_ens) == 0L) {
    stop("[FATAL] No expected marker genes found in expression matrix. Check ID mapping or input.")
  }

  # Report findings
  cat(sprintf("[INFO] Found %d/%d expected marker genes.\n",
              length(marker_ens), length(marker_symbols)))

  # Identify missing markers (case-insensitive)
  found_symbols <- ens2sym_union$symbol[ens2sym_union$ensembl %in% marker_ens]
  missing_syms <- setdiff(tolower(marker_symbols), tolower(found_symbols))

  if (length(missing_syms) > 0) {
    cat(sprintf("[WARN] %d expected marker genes are missing:\n  %s\n",
                length(missing_syms), paste(missing_syms, collapse = ", ")))
  }
}


#' @title Validate PAM50 Gene Coverage
#'
#' @description
#' Given a subsetted PAM50 expression matrix (e.g., from `subset_to_pam50`),
#' this function reports the number and identity of missing genes relative to the official list.
#'
#' @param pam50_matrix Matrix of gene expression for PAM50 genes (rows = gene symbols).
#' @param pam50_symbols Vector of expected PAM50 gene symbols.
#' @param outdir Optional directory to save list of missing symbols.
#'
#' @return Vector of missing PAM50 gene symbols (uppercase).
#'
missing_pam50_genes <- function(pam50_matrix, pam50_symbols, outdir = NULL) {
  expected <- toupper(pam50_symbols)
  observed <- rownames(pam50_matrix)
  missing <- setdiff(expected, observed)

  cat(sprintf("[PAM50] Matched %d/%d genes.\n", length(observed), length(expected)))
  
  if (length(missing) > 0) {
    cat("[PAM50] Missing symbols:\n  ", paste(missing, collapse = ", "), "\n")
    if (!is.null(outdir)) {
      writeLines(missing, file.path(outdir, "PAM50_missing_symbols.txt"))
    }
  }

  return(missing)
}
