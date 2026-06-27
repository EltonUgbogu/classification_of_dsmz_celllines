load_required_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
}

script_directory <- function() {
  command_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- command_arguments[grepl("^--file=", command_arguments)]
  if (length(file_argument) == 0) {
    return(file.path("research_framework", "optimisation", "scripts"))
  }
  dirname(normalizePath(sub("^--file=", "", file_argument[[1]]), mustWork = TRUE))
}

parse_key_value_arguments <- function(default_config = file.path("research_framework", "optimisation", "config", "optimisation.yaml")) {
  trailing_arguments <- commandArgs(trailingOnly = TRUE)
  parsed_arguments <- list(config = default_config)
  if (length(trailing_arguments) == 0) {
    return(parsed_arguments)
  }
  argument_index <- 1
  while (argument_index <= length(trailing_arguments)) {
    argument_name <- trailing_arguments[[argument_index]]
    if (startsWith(argument_name, "--")) {
      key <- sub("^--", "", argument_name)
      value_index <- argument_index + 1
      if (value_index <= length(trailing_arguments) && !startsWith(trailing_arguments[[value_index]], "--")) {
        parsed_arguments[[key]] <- trailing_arguments[[value_index]]
        argument_index <- argument_index + 2
      } else {
        parsed_arguments[[key]] <- TRUE
        argument_index <- argument_index + 1
      }
    } else {
      argument_index <- argument_index + 1
    }
  }
  parsed_arguments
}

load_optimisation_config <- function(config_path) {
  load_required_packages(c("yaml"))
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path, call. = FALSE)
  }
  yaml::read_yaml(config_path)
}

config_value <- function(config, path, default = NULL) {
  current_value <- config
  for (path_element in path) {
    if (is.null(current_value[[path_element]])) {
      return(default)
    }
    current_value <- current_value[[path_element]]
  }
  current_value
}

repository_path <- function(...) {
  file.path(...)
}

ensure_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

ensure_parent_directory <- function(path) {
  ensure_directory(dirname(path))
}

read_table_if_exists <- function(path) {
  load_required_packages(c("data.table"))
  if (!file.exists(path)) {
    return(NULL)
  }
  data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)
}

write_tsv_table <- function(table, path) {
  load_required_packages(c("data.table"))
  ensure_parent_directory(path)
  data.table::fwrite(table, path, sep = "\t", quote = FALSE, na = "NA")
}

write_tsv_gzip <- function(table, path) {
  ensure_parent_directory(path)
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.table(table, file = connection, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

write_markdown <- function(lines, path, append = FALSE) {
  ensure_parent_directory(path)
  writeLines(lines, con = path, sep = "\n", useBytes = TRUE)
}

append_markdown <- function(lines, path) {
  ensure_parent_directory(path)
  cat(paste0(lines, collapse = "\n"), "\n", file = path, append = TRUE)
}

relative_path_or_na <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  path
}

collapse_examples <- function(paths, n = 5) {
  if (length(paths) == 0) {
    return("")
  }
  paste(utils::head(paths, n), collapse = ";")
}

discover_files <- function(root, pattern) {
  if (!dir.exists(root)) {
    return(character(0))
  }
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

normalise_file_path <- function(path) {
  gsub("\\\\", "/", path)
}

parse_representation_from_path <- function(path) {
  normalised_path <- normalise_file_path(path)
  path_parts <- strsplit(normalised_path, "/", fixed = TRUE)[[1]]
  unsupervised_index <- match("unsupervised", path_parts)
  cohort_id <- if (!is.na(unsupervised_index) && length(path_parts) >= unsupervised_index + 1) {
    path_parts[[unsupervised_index + 1]]
  } else {
    NA_character_
  }
  neighbourhood_index <- match("tumour_neighbourhoods", path_parts)
  representation_id <- if (!is.na(neighbourhood_index) && length(path_parts) >= neighbourhood_index + 1) {
    path_parts[[neighbourhood_index + 1]]
  } else {
    NA_character_
  }
  distance_metric <- if (grepl("_corr$", representation_id)) {
    "correlation"
  } else if (grepl("_euc$", representation_id)) {
    "euclidean"
  } else if (grepl("_cosine$", representation_id)) {
    "cosine"
  } else if (grepl("_euclidean$", representation_id)) {
    "euclidean"
  } else {
    "unspecified"
  }
  feature_set_id <- sub("_(corr|euc|cosine|euclidean)$", "", representation_id)
  transformation_id <- if (grepl("PCA", feature_set_id, ignore.case = TRUE)) {
    "pca_or_loading"
  } else {
    "selected_expression"
  }
  list(
    cohort_id = cohort_id,
    representation_id = representation_id,
    feature_set_id = feature_set_id,
    distance_metric = distance_metric,
    transformation_id = transformation_id
  )
}

final_consensus_paths <- function(config) {
  thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
  discover_files(thesis_results_root, "^Final_consensus_tumour_neighbourhoods_.*\\.tsv$")
}

graph_edge_paths <- function(config) {
  thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
  discover_files(thesis_results_root, "^cell_line_similarity_graph_edges_.*\\.tsv$")
}

node_annotation_paths <- function(config) {
  thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
  discover_files(thesis_results_root, "^cell_line_similarity_graph_node_annotations_.*\\.tsv$")
}

p_consensus_summary_paths <- function(config, filename) {
  thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
  discover_files(thesis_results_root, paste0("^", filename, "$"))
}

metadata_paths <- function(config) {
  thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
  discover_files(thesis_results_root, "^(metadata|metadata_with_components|joint_metadata|pan_cancer_feature_expr_metadata).*\\.tsv$")
}

standard_cell_line_column <- function(table) {
  candidate_columns <- c("cell_line", "cell_line_id", "cell_tech_id", "cell_line1", "from")
  candidate_columns[candidate_columns %in% names(table)][1]
}

standard_tumour_column <- function(table) {
  candidate_columns <- c("tumour_id", "tumour_sample_id", "tumour", "sample_id")
  candidate_columns[candidate_columns %in% names(table)][1]
}

safe_numeric <- function(value, default = NA_real_) {
  converted_value <- suppressWarnings(as.numeric(value))
  converted_value[is.na(converted_value)] <- default
  converted_value
}

classify_quantile <- function(value, high = 0.8, moderate = 0.5, high_label = "high", moderate_label = "moderate", low_label = "low") {
  ifelse(
    is.na(value), "not_estimated",
    ifelse(value >= high, high_label, ifelse(value >= moderate, moderate_label, low_label))
  )
}

bh_q_value <- function(p_values) {
  stats::p.adjust(p_values, method = "BH")
}

