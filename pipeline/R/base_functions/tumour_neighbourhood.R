# ==============================================================================
# tumour_neighbourhood.R
# Tumour Neighbourhood Computation Algorithm
# ==============================================================================
#
# PURPOSE:
# This module provides the core algorithm for computing local tumour
# neighbourhoods for each DSMZ cell line within an integrated expression space.
# For each cell line, the algorithm identifies the most similar tumour
# samples (non-DSMZ cohort) based on expression distance, enabling biological
# interpretation of cell line–tumour relationships.
#
# The algorithm operates within the context of consensus clustering results,
# restricting neighbour candidates to tumours in the same cluster as the query
# cell line. This constraint ensures that neighbourhoods reflect biologically
# coherent sample groupings.
#
# DISTANCE METRICS:
# Two distance metrics are supported:
#   - Euclidean: Standard L2 norm, suitable for comparing absolute expression
#   - Correlation: 1 - Pearson r, captures similarity in expression patterns
#
# ==============================================================================

# ------------------------------------------------------------------------------
# FUNCTION: compute_tumour_neighbourhoods
# ------------------------------------------------------------------------------
# Computes local tumour neighbourhoods for each DSMZ cell line within a given
# clustering solution. The function identifies the nearest tumour samples
# (non-DSMZ cohort) for each cell line based on expression distance.
#
# PARAMETERS:
#   emb_mat:    Numeric matrix of expression values (samples × features).
#               Rownames must correspond to sample identifiers.
#   cluster_m:  Named vector of cluster assignments. Names must match sample
#               identifiers in emb_mat.
#   dataset:    Named vector indicating sample origin ("DSMZ" for cell lines,
#               any other value for tumour samples). Names must match sample
#               identifiers.
#   top_frac:   Fraction of candidate tumours to include in neighbourhood.
#               Default: 0.10 (top 10%).
#   top_n_min:  Minimum neighbourhood size. Ensures robustness for small
#               clusters. Default: 30.
#   top_n_max:  Maximum neighbourhood size. Prevents oversized neighbourhoods.
#               Default: 200.
#   method_id:  String identifier for the clustering method. Used for output
#               labelling and downstream aggregation.
#   distance:   Distance metric to use. One of "euclidean" or "correlation".
#
# RETURNS:
#   A list containing:
#     - neighbourhoods: List mapping each cell line to its top tumour neighbours
#     - long_df:        Tidy data frame with all pairwise distances and ranks
#     - method_id:      The clustering method identifier
#     - distance:       The distance metric used
#     - params:         List of neighbourhood size parameters
#     - n_dsmz:         Number of DSMZ cell lines processed
#     - n_tumour:       Number of tumour (non-DSMZ) samples in the dataset
# ------------------------------------------------------------------------------

compute_tumour_neighbourhoods <- function(
  emb_mat,
  cluster_m,
  dataset,
  top_frac  = 0.10,
  top_n_min = 30,
  top_n_max = 200,
  method_id = "unspecified_method",
  distance  = c("euclidean", "correlation")
) {
  
  # ----------------------------------------------------------------------------
  # SECTION 1: DISTANCE METRIC ARGUMENT MATCHING
  # ----------------------------------------------------------------------------
  # The match.arg() function validates and normalises the distance argument
  # against the allowed values specified in the function signature. This
  # provides both input validation and partial matching (e.g., "euc" matches
  # "euclidean").
  
  distance <- match.arg(distance)
  
  # ----------------------------------------------------------------------------
  # SECTION 2: INPUT VALIDATION
  # ----------------------------------------------------------------------------
  # Comprehensive input validation prevents cryptic downstream errors by
  # catching invalid inputs early. Each validation check provides a descriptive
  # error message to guide the user towards the correct input format.
  
  # Validate expression matrix structure.
  # The matrix must be numeric with samples as rows and features as columns.
  if (!is.matrix(emb_mat) || !is.numeric(emb_mat)) {
    stop("emb_mat must be a numeric matrix (samples × features).")
  }
  
  # Verify that sample identifiers are present as rownames.
  # These identifiers are essential for matching across input objects.
  if (is.null(rownames(emb_mat))) {
    stop("emb_mat must have rownames corresponding to sample IDs.")
  }
  
  # Validate cluster assignment vector.
  if (missing(cluster_m) || is.null(cluster_m)) {
    stop("cluster_m is required. Supply a named vector of cluster labels.")
  }
  
  # Verify that cluster vector has named elements for sample matching.
  if (is.null(names(cluster_m))) {
    stop("cluster_m must be named with sample IDs matching rownames(emb_mat).")
  }
  
  # Validate dataset origin vector.
  if (missing(dataset) || is.null(dataset)) {
    stop("dataset vector is required (use 'DSMZ' for cell lines; any other value for tumour samples).")
  }
  
  # Verify that dataset vector has named elements for sample matching.
  if (is.null(names(dataset))) {
    stop("dataset must be named with sample IDs.")
  }

  # ----------------------------------------------------------------------------
  # SECTION 3: INPUT ALIGNMENT
  # ----------------------------------------------------------------------------
  # The three input objects (expression matrix, cluster vector, dataset vector)
  # may contain different sample sets. This section identifies the common
  # samples and subsets all inputs to ensure perfect alignment.
  #
  # The nested intersect() calls compute the three-way intersection efficiently:
  #   1. First intersect: rownames(emb_mat) ∩ names(cluster_m)
  #   2. Second intersect: result ∩ names(dataset)

  common_ids <- intersect(rownames(emb_mat), intersect(names(cluster_m), names(dataset)))
  
  if (length(common_ids) == 0) {
    stop("No overlapping sample IDs between emb_mat, cluster_m, and dataset.")
  }

  # Subset all inputs to common samples, maintaining consistent ordering.
  # The drop = FALSE argument prevents dimension reduction when subsetting
  # to a single row.
  emb_mat   <- emb_mat[common_ids, , drop = FALSE]
  cluster_m <- cluster_m[common_ids]
  dataset   <- dataset[common_ids]

  # ----------------------------------------------------------------------------
  # SECTION 4: SAMPLE PARTITIONING
  # ----------------------------------------------------------------------------
  # Samples are partitioned into DSMZ cell lines (query set) and tumour samples
  # (target set, defined as all non-DSMZ samples). The neighbourhood algorithm
  # computes distances from each cell line to all candidate tumours.

  # Extract sample identifiers by dataset origin.
  # Cohort-agnostic: tumours are defined as all non-DSMZ samples.
  dsmz_ids   <- names(dataset)[dataset == "DSMZ"]
  tumour_ids <- names(dataset)[dataset != "DSMZ"]   # cohort-agnostic tumour set

  # Validate that both sample types are present.
  if (length(dsmz_ids) == 0) stop("No DSMZ samples found.")
  if (length(tumour_ids) == 0) stop("No tumour (non-DSMZ) samples found.")
  
  # Safety checks: ensure cluster vector and embedding matrix cover all IDs
  # These catch ID mismatches early.
  missing_clust <- setdiff(c(dsmz_ids, tumour_ids), names(cluster_m))
  if (length(missing_clust) > 0) {
    stop("cluster_m missing assignments for ", length(missing_clust), " samples. Examples: ",
         paste(head(missing_clust, 10), collapse = ", "))
  }
  
  missing_emb <- setdiff(c(dsmz_ids, tumour_ids), rownames(emb_mat))
  if (length(missing_emb) > 0) {
    stop("emb_mat missing rows for ", length(missing_emb), " samples. Examples: ",
         paste(head(missing_emb, 10), collapse = ", "))
  }

  # ----------------------------------------------------------------------------
  # SECTION 5: DISTANCE FUNCTION DEFINITIONS
  # ----------------------------------------------------------------------------
  # Two distance functions are defined as closures within the main function.
  # This encapsulation keeps the distance logic modular while avoiding
  # global namespace pollution.
  
  # --------------------------------------------------------------------------
  # euclidean_dist: Computes Euclidean (L2) distance
  # --------------------------------------------------------------------------
  # For a query vector q and target matrix T (rows = samples, cols = features),
  # the Euclidean distance to each target sample is:
  #   d(q, t_i) = sqrt(sum((q - t_i)^2))
  #
  # PARAMETERS:
  #   query_vec:  Numeric vector representing the query sample (cell line)
  #   target_mat: Numeric matrix of target samples (tumours × features)
  #
  # RETURNS:
  #   Named numeric vector of distances (one per target sample)
  
  euclidean_dist <- function(query_vec, target_mat) {
    # Compute the difference between each target row and the query vector.
    # The sweep() function subtracts query_vec from each row of target_mat.
    # MARGIN = 2 indicates column-wise operation (subtract from each column).
    diffs <- sweep(target_mat, 2, query_vec, "-")
    
    # Compute the L2 norm (Euclidean length) of each difference vector.
    # rowSums() efficiently computes the sum of squared differences per row.
    sqrt(rowSums(diffs^2))
  }
  
  # --------------------------------------------------------------------------
  # correlation_dist: Computes correlation-based distance
  # --------------------------------------------------------------------------
  # Correlation distance is defined as 1 - Pearson correlation coefficient.
  # This metric captures similarity in expression patterns rather than
  # absolute expression levels, making it robust to global scaling differences.
  #
  # Distance interpretation:
  #   - d = 0: Perfect positive correlation (identical patterns)
  #   - d = 1: No correlation (orthogonal patterns)
  #   - d = 2: Perfect negative correlation (inverted patterns)
  #
  # PARAMETERS:
  #   query_vec:  Numeric vector representing the query sample (cell line)
  #   target_mat: Numeric matrix of target samples (tumours × features)
  #
  # RETURNS:
  #   Named numeric vector of distances (one per target sample)
  
  correlation_dist <- function(query_vec, target_mat) {
    # Compute Pearson correlation between each target row and the query.
    # The cor() function requires samples as columns, so target_mat is
    # transposed. The result is a column vector of correlations.
    # 
    # use = "pairwise.complete.obs" handles missing values by using all
    # complete pairs of observations for each correlation.
    cors <- cor(t(target_mat), query_vec,
                use = "pairwise.complete.obs",
                method = "pearson")
    
    # Convert correlation to distance and return as a numeric vector.
    # The as.numeric() call drops the matrix structure from cor() output.
    as.numeric(1 - cors)
  }
  
  # Select the appropriate distance function based on the distance argument.
  # This conditional assignment enables clean dispatch without repeated checks.
  dist_fun <- if (distance == "euclidean") euclidean_dist else correlation_dist
  
  cat("[INFO] Using", distance, "distance for tumour neighbourhood computation\n")

  # ----------------------------------------------------------------------------
  # SECTION 6: NEIGHBOURHOOD COMPUTATION LOOP
  # ----------------------------------------------------------------------------
  # The main algorithm iterates over each DSMZ cell line, computing distances
  # to candidate tumours and identifying the nearest neighbours. Results are
  # stored in two complementary formats:
  #   1. neighbourhoods: Compact list of top neighbours per cell line
  #   2. long_records:   Detailed records for consensus analysis
  
  # Pre-allocate result storage for efficiency.
  # Using vector("list", n) is more efficient than growing lists incrementally.
  neighbourhoods <- vector("list", length(dsmz_ids))
  names(neighbourhoods) <- dsmz_ids
  long_records <- vector("list", length(dsmz_ids))

  # Iterate over each DSMZ cell line as the query sample.
  # Note: cell_tech_id uses collapsed DSMZ canonical IDs.
  for (cell_tech_id in dsmz_ids) {
    
    # ------------------------------------------------------------------------
    # STEP 1: RETRIEVE CLUSTER ASSIGNMENT
    # ------------------------------------------------------------------------
    # The cell line's cluster assignment determines which tumours are
    # considered as neighbourhood candidates.
    
    cluster_of_cell <- cluster_m[[cell_tech_id]]
    
    # Handle missing cluster assignments (stop in pipelines to surface data issues).
    if (is.na(cluster_of_cell) || is.null(cluster_of_cell)) {
      stop("Missing cluster assignment for DSMZ sample: ", cell_tech_id)
    }

    # ------------------------------------------------------------------------
    # STEP 2: IDENTIFY CANDIDATE TUMOURS
    # ------------------------------------------------------------------------
    # Candidate tumours are restricted to those in the same cluster as the
    # query cell line. This constraint ensures biological coherence of the
    # neighbourhood.
    
    candidate_tumours <- tumour_ids[cluster_m[tumour_ids] == cluster_of_cell]

    # Fallback: if the cluster contains too few tumours, expand to all tumours.
    # This ensures meaningful neighbourhood computation even for small or
    # singleton clusters.
    if (length(candidate_tumours) < 10) {
      message("[", method_id, "] ", cell_tech_id,
              ": only ", length(candidate_tumours),
              " tumours in cluster → using all tumours.")
      candidate_tumours <- tumour_ids
    }

    # ------------------------------------------------------------------------
    # STEP 3: COMPUTE PAIRWISE DISTANCES
    # ------------------------------------------------------------------------
    # Extract the query expression profile and target expression matrix,
    # then compute distances using the selected metric.
    
    # Extract the cell line's expression profile as a vector.
    # The drop = TRUE argument converts the single-row matrix to a vector.
    query_vec <- emb_mat[cell_tech_id, , drop = TRUE]
    
    # Extract the expression profiles of all candidate tumours.
    target_mat <- emb_mat[candidate_tumours, , drop = FALSE]
    
    # Compute distances using the selected distance function.
    distances <- dist_fun(query_vec, target_mat)

    # ------------------------------------------------------------------------
    # STEP 4: RANK TUMOURS BY PROXIMITY
    # ------------------------------------------------------------------------
    # Tumours are sorted by ascending distance (nearest first).
    
    # The order() function returns indices that would sort the vector.
    ord <- order(distances)
    
    # Apply the sorting order to obtain ranked tumour lists.
    ranked_tumours <- candidate_tumours[ord]
    ranked_dist   <- distances[ord]

    # ------------------------------------------------------------------------
    # STEP 5: DETERMINE ADAPTIVE NEIGHBOURHOOD SIZE
    # ------------------------------------------------------------------------
    # The neighbourhood size is adaptive, balancing three constraints:
    #   1. top_frac:  Include at least this fraction of candidates
    #   2. top_n_min: Ensure a minimum neighbourhood size for robustness
    #   3. top_n_max: Prevent oversized neighbourhoods
    #
    # The final size is: min(max(top_n_min, ceil(top_frac * n)), top_n_max, n)
    
    n_total <- length(ranked_tumours)
    
    # Compute the fraction-based size (ceiling ensures at least 1 if n > 0).
    n_frac  <- ceiling(top_frac * n_total)
    
    # Apply minimum constraint.
    n_keep  <- max(top_n_min, n_frac)
    
    # Apply maximum constraint and ensure we don't exceed available tumours.
    n_keep  <- min(n_keep, top_n_max, n_total)

    # Extract the top-ranked tumours as the neighbourhood.
    top_tumours <- ranked_tumours[seq_len(n_keep)]

    # ------------------------------------------------------------------------
    # STEP 6: STORE RESULTS
    # ------------------------------------------------------------------------
    
    # Store the neighbourhood as a simple character vector.
    neighbourhoods[[cell_tech_id]] <- top_tumours

    # Construct a detailed data frame for consensus analysis.
    # This long format facilitates aggregation across multiple methods.
    # Use cell_tech_id directly (collapsed DSMZ canonical IDs).
    long_records[[cell_tech_id]] <- data.frame(
      method       = method_id,
      cell_tech_id = cell_tech_id,
      tumor_id     = ranked_tumours,
      rank         = seq_along(ranked_tumours),
      distance     = ranked_dist,
      in_top       = ranked_tumours %in% top_tumours,
      cluster      = cluster_of_cell,
      stringsAsFactors = FALSE
    )
  }

  # ----------------------------------------------------------------------------
  # SECTION 7: RESULT AGGREGATION
  # ----------------------------------------------------------------------------
  # Individual cell line records are combined into a single tidy data frame.
  # Note: bind_rows() requires dplyr; this module assumes it is loaded upstream.
  long_df <- bind_rows(long_records)

  # ----------------------------------------------------------------------------
  # SECTION 8: RETURN VALUE CONSTRUCTION
  # ----------------------------------------------------------------------------

  list(
    neighbourhoods = neighbourhoods,
    long_df        = long_df,
    method_id      = method_id,
    distance       = distance,
    params         = list(
      top_frac  = top_frac,
      top_n_min = top_n_min,
      top_n_max = top_n_max
    ),
    n_dsmz         = length(dsmz_ids),
    n_tumour       = length(tumour_ids)
  )
}