# lib_config.R
# Shared config reading utilities for R scripts.
#
# Profile merging must produce exactly the same effective configuration as the
# Snakefile's Python deep_merge(), because both read config/config.yaml and both
# feed the same analysis. The two implementations are kept deliberately aligned:
#
#   * a mapping (YAML block/flow mapping) present in both defaults and the
#     profile is merged key by key, recursively;
#   * anything else -- scalar, sequence, or a type change -- is replaced
#     wholesale by the profile value.
#
# The sequence case is the one that matters scientifically. A YAML sequence
# parses to an unnamed R list, so a naive `is.list()` test treats it as a
# mapping and recurses over `names()`, which is NULL; the loop body then never
# runs and the profile's sequence is silently discarded in favour of the
# defaults. Feature-method lists, distance lists, k grids and the eligible
# clustering formulations are all sequences, so that failure mode would let
# Snakemake and the R scripts run with different scientific universes without
# raising. is_config_mapping() below draws the distinction explicitly.

suppressPackageStartupMessages({
  library(yaml)
})

# Null coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

# TRUE only for a YAML mapping: a named list whose names are all non-empty.
# A YAML sequence parses to an unnamed list and is therefore not a mapping, so
# it is replaced rather than merged -- matching Python, where only dict/dict
# pairs recurse.
is_config_mapping <- function(x) {
  is.list(x) &&
    length(x) > 0L &&
    !is.null(names(x)) &&
    all(nzchar(names(x)))
}

# Deep merge two configuration trees.
# Mirrors the Snakefile's deep_merge(): recurse only when both sides are
# mappings; otherwise the override wins outright.
deep_merge <- function(base, override) {
  if (!is.list(base)) {
    return(override)
  }
  out <- base
  for (nm in names(override)) {
    if (is_config_mapping(override[[nm]]) && is_config_mapping(out[[nm]])) {
      out[[nm]] <- deep_merge(out[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

# Read config with profile merging.
# Produces the same effective configuration as the Snakefile for a given profile.
read_profiled_config <- function(config_file, profile_override = NULL) {
  cfg0 <- yaml::read_yaml(config_file)

  # If config has a profiles structure, merge defaults + the selected profile.
  # `cfg0[["profile"]]` is an exact lookup: `cfg0$profile` would partially match
  # the `profiles` key and return the whole profile table instead of a name.
  if (!is.null(cfg0[["profiles"]])) {
    prof <- profile_override %||% cfg0[["profile"]] %||% names(cfg0[["profiles"]])[[1]]

    if (is.null(prof) || !is.character(prof) || length(prof) != 1L ||
        !(prof %in% names(cfg0[["profiles"]]))) {
      stop("Profile not found in config: ", paste(as.character(prof), collapse = ", "),
           ". Available profiles: ", paste(names(cfg0[["profiles"]]), collapse = ", "))
    }

    defaults <- cfg0[["defaults"]] %||% list()
    profcfg  <- cfg0[["profiles"]][[prof]]
    cfg <- deep_merge(defaults, profcfg)

    # Record the resolved profile name for reference by callers.
    cfg$profile <- prof

    return(cfg)
  }

  # If no profiles structure, return config as-is.
  cfg0$profile <- cfg0[["profile"]] %||% "default"
  cfg0
}
