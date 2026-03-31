# lib_config.R
# Shared config reading utilities for R scripts
# Handles profile merging (defaults + profile-specific overrides) like Snakefile does

suppressPackageStartupMessages({
  library(yaml)
})

# Null coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

# Deep merge two lists (recursive)
deep_merge <- function(base, override) {
  out <- base
  for (nm in names(override)) {
    if (is.list(override[[nm]]) && is.list(out[[nm]])) {
      out[[nm]] <- deep_merge(out[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

# Read config with profile merging
# Merges defaults + profile-specific config, matching Snakefile behavior
read_profiled_config <- function(config_file, profile_override = NULL) {
  cfg0 <- yaml::read_yaml(config_file)

  # If config has profiles structure, merge defaults + profile
  if (!is.null(cfg0$profiles)) {
    # Determine which profile to use
    prof <- profile_override %||% cfg0$profile %||% names(cfg0$profiles)[[1]]
    
    if (is.null(prof) || !(prof %in% names(cfg0$profiles))) {
      stop("Profile not found in config: ", prof, 
           ". Available profiles: ", paste(names(cfg0$profiles), collapse = ", "))
    }
    
    # Merge defaults + profile-specific config
    defaults <- cfg0$defaults %||% list()
    profcfg  <- cfg0$profiles[[prof]]
    cfg <- deep_merge(defaults, profcfg)
    
    # Add profile name to config for reference
    cfg$profile <- prof
    
    return(cfg)
  }

  # If no profiles structure, return config as-is
  cfg0$profile <- cfg0$profile %||% "default"
  cfg0
}

