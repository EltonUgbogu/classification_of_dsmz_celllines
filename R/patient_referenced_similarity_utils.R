threshold_restrict_mean_profiles <- function(mean_mat, threshold) {
  if (!is.matrix(mean_mat) || !is.numeric(mean_mat)) {
    stop("mean_mat must be a numeric matrix.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("threshold must be a single numeric value.")
  }
  if (threshold <= 0 || threshold > 1) {
    stop("threshold must satisfy 0 < threshold <= 1.")
  }
  ifelse(mean_mat >= threshold, mean_mat, 0)
}


pairwise_active_tumour_union_mask <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  x > 0 | y > 0
}


pairwise_active_tumour_union_counts <- function(x, y) {
  active <- pairwise_active_tumour_union_mask(x, y)
  list(
    active = active,
    n_pairwise_active_tumours = sum(active, na.rm = TRUE),
    n_shared_selected_tumours = sum(x[active] > 0 & y[active] > 0, na.rm = TRUE)
  )
}


compute_pairwise_active_union_similarity <- function(x, y, metric = c("pearson", "jaccard")) {
  metric <- match.arg(metric)
  counts <- pairwise_active_tumour_union_counts(x, y)
  active <- counts$active
  n_pairwise_active_tumours <- counts$n_pairwise_active_tumours
  n_shared_selected_tumours <- counts$n_shared_selected_tumours

  if (n_pairwise_active_tumours == 0L) {
    return(list(
      similarity = NA_real_,
      n_pairwise_active_tumours = 0L,
      n_shared_selected_tumours = 0L,
      undefined_similarity_reason = "empty_pairwise_active_tumour_union"
    ))
  }

  x_active <- x[active]
  y_active <- y[active]

  if (metric == "pearson") {
    if (n_pairwise_active_tumours < 2L) {
      return(list(
        similarity = NA_real_,
        n_pairwise_active_tumours = n_pairwise_active_tumours,
        n_shared_selected_tumours = n_shared_selected_tumours,
        undefined_similarity_reason = "fewer_than_two_active_tumours"
      ))
    }
    if (isTRUE(all.equal(stats::var(x_active), 0)) ||
        isTRUE(all.equal(stats::var(y_active), 0))) {
      return(list(
        similarity = NA_real_,
        n_pairwise_active_tumours = n_pairwise_active_tumours,
        n_shared_selected_tumours = n_shared_selected_tumours,
        undefined_similarity_reason = "zero_variance_active_profile"
      ))
    }
    similarity <- suppressWarnings(stats::cor(x_active, y_active, method = "pearson"))
    reason <- if (is.na(similarity)) "pearson_undefined" else NA_character_
    return(list(
      similarity = similarity,
      n_pairwise_active_tumours = n_pairwise_active_tumours,
      n_shared_selected_tumours = n_shared_selected_tumours,
      undefined_similarity_reason = reason
    ))
  }

  x_binary <- as.integer(x_active > 0)
  y_binary <- as.integer(y_active > 0)
  denom <- sum(x_binary == 1L | y_binary == 1L, na.rm = TRUE)
  if (denom == 0L) {
    return(list(
      similarity = NA_real_,
      n_pairwise_active_tumours = n_pairwise_active_tumours,
      n_shared_selected_tumours = n_shared_selected_tumours,
      undefined_similarity_reason = "empty_pairwise_active_tumour_union"
    ))
  }
  similarity <- sum(x_binary == 1L & y_binary == 1L, na.rm = TRUE) / denom
  list(
    similarity = similarity,
    n_pairwise_active_tumours = n_pairwise_active_tumours,
    n_shared_selected_tumours = n_shared_selected_tumours,
    undefined_similarity_reason = NA_character_
  )
}
