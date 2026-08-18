#!/usr/bin/env snakemake

# Classification of DSMZ Cell Lines to Cancer Types Based on Transcriptomic Analysis

# Unsupervised Clustering and Tumour Neighbourhood Analysis Pipeline

import json
import os
import re
import shlex
import sys
from pathlib import Path
from snakemake.shell import shell
from snakemake.io import expand
from snakemake.utils import min_version

min_version("7.32.4")

# =============================================================================
# IMPORTS / CONFIG / ENVIRONMENT PATHS
# =============================================================================
# Stage role: resolves pipeline roots, config, and conda environment paths before any rules are declared.
# Anchors all path resolution to Snakemake's workflow directory, preventing
# inadvertent resolution relative to the caller's current working directory.
REPO_ROOT = Path(workflow.basedir).resolve()

def repo_path(*parts):
    """Return an absolute path to a repository-owned file or directory."""
    return str(REPO_ROOT.joinpath(*parts))

def cfg_path(value):
    """Resolve a configured path against the repository root when relative."""
    path = Path(str(value)).expanduser()
    return str(path if path.is_absolute() else REPO_ROOT / path)

BASE = str(REPO_ROOT)
SCRIPTS_DIR = repo_path("scripts")
PIPE_ROOT = BASE

# Centralises conda environment paths to avoid hardcoded absolute strings in rules.
CONDA_ENV_R = repo_path("envs", "tcga-r-env.yaml")
CONDA_ENV_R_BASE = repo_path("envs", "r-base.yaml")
CONDA_ENV_PY = repo_path("envs", "python-graph-env.yaml")
CONDA_ENV_SHELL = repo_path("envs", "pipeline-shell-tools.yaml")

# Resolves the config file path once at pipeline load time.
CFGFILE_ABS = repo_path("config", "config.yaml")

# Instructs Snakemake to load configuration from the resolved absolute path.
configfile: CFGFILE_ABS

# Canonical seed used by stochastic rules instead of embedding per-rule constants.
PIPELINE_SEED = int(config.get("validation", {}).get("seed", 42))
if PIPELINE_SEED < 0:
    raise ValueError("validation.seed must be a non-negative integer")

# Verbose mode: enables debug-level prints when --config verbose=1 is passed.
VERBOSE = bool(int(config.get("verbose", 0)))

def vprint(*args):
    """Prints debug messages only when verbose mode is active."""
    if VERBOSE:
        print(*args)

# Reads the pipeline profile from --config; falls back to config file, then 'default'.
# This temporary value is overwritten after profile validation below.
profile_name = config.get("pipeline_profile") or config.get("profile") or "default"

# Initialises log root under a per-profile subdirectory; overwritten after validation.
LOGROOT = os.path.join(PIPE_ROOT, "logs", profile_name)
IS_DRYRUN = any(
    arg in ("-n", "--dry-run", "--dryrun")
    or (arg.startswith("-") and not arg.startswith("--") and "n" in arg[1:])
    for arg in sys.argv[1:]
)

# Defines a safe default for PIPELINE_TARGET so DAG parsing never fails before
# profile-specific targets are resolved later in the Snakefile.
PIPELINE_TARGET = []

# -----------------------------------------------------------------------------
# DEFAULT TARGET (FIRST RULE)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# RULE: all
#
# Snakemake builds this profile-aware target when no explicit rule is requested.
# -----------------------------------------------------------------------------
rule all:
    input:
        # Lazily evaluated so the target list is not resolved until after profile
        # setup and PIPELINE_TARGET is overwritten with profile-specific outputs.
        lambda wc: PIPELINE_TARGET,
        # Pan-cancer UMAP (optional, independent target reference)
        # Inactive target path retained for optional default-target assembly:
        # os.path.join(BASE, "results", "unsupervised", "pan_cancer", "umap", "pan_cancer_mx_umap_coords.tsv")


# =============================================================================
# GENERIC HELPERS
# =============================================================================
# Stage role: defines reusable path and configuration helpers used by later workflow sections.
# Utility function: converts a config-relative path to an absolute path.
def abspath(p):
    return cfg_path(p)


STUDY_DESIGN_REL = config.get("study_design", {}).get("file", "config/study_design.yaml")
STUDY_DESIGN_FILE = abspath(STUDY_DESIGN_REL)

# =============================================================================
# PROFILE DETECTION / CONFIG MERGE
# =============================================================================
# Stage role: selects the active profile and merges profile-specific settings into the default configuration.

# Method helper: recursively merges default and profile-specific config without mutating either source dictionary.
def deep_merge(base, override):
    """Recursively merges two configuration dictionaries without mutating either input."""
    merged = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(merged.get(k), dict):
            merged[k] = deep_merge(merged[k], v)
        else:
            merged[k] = v
    return merged


# Method helper: resolves the active profile and returns the merged configuration used by all dependent rules.
def get_profile_cfg(profile):
    """Returns the merged configuration for a named profile, used for cross-profile pan-cancer access."""
    default_cfg = config.get("defaults", {})
    profiles = config.get("profiles", {})
    if profile not in profiles:
        raise ValueError(f"Profile '{profile}' not found in config.profiles")
    return deep_merge(default_cfg, profiles[profile])


# Method helper: resolves per-profile VST inputs used when constructing the multicohort expression matrix.
def validate_authoritative_vst_path(profile, path):
    """Fails closed when an active tumour profile points at a superseded matrix."""
    tumour_profiles = {"brca", "nbl", "rbl"}
    required_name = "joint_vst_purity_filtered_post_bc.rds"
    basename = os.path.basename(str(path))
    if profile in tumour_profiles and basename != required_name:
        raise ValueError(
            f"Profile '{profile}' must use {required_name}; configured path is {path}. "
            "Only the authoritative purity-filtered, post-batch-correction joint VST path is accepted; "
            "legacy downstream matrices are superseded."
        )
    if str(path).endswith("vst_RBL_joint_batch_corrected.rds"):
        raise ValueError("The removed RBL placeholder is not a valid active matrix")
    return path


def profile_vst_joint_abs(profile):
    """Resolves the absolute path to vst_joint_rds for a given profile."""
    cfgp = get_profile_cfg(profile)
    p = cfgp.get("paths", {}).get("vst_joint_rds")
    if not p:
        raise KeyError(f"paths.vst_joint_rds missing for profile '{profile}'")
    return abspath(validate_authoritative_vst_path(profile, p))


if "profiles" in config:
    profiles = config.get("profiles", {})
    # Reads the mandatory pipeline_profile flag; raises an informative error if absent.
    profile_name = config.get("pipeline_profile")

    # Enforces explicit profile selection: pipeline_profile must be supplied via --config.
    if not profile_name:
        available = ", ".join(profiles.keys())
        raise ValueError(
            "Missing --config pipeline_profile=<name>. "
            f"Available profiles: {available}. "
            "Example: snakemake --config pipeline_profile=brca"
        )

    if profile_name not in profiles:
        available = ", ".join(profiles.keys())
        raise ValueError(
            f"Profile '{profile_name}' not found in config. "
            f"Available profiles: {available}"
        )

    vprint(f"[Snakefile] pipeline_profile={profile_name}")

    # Merges default configuration with profile-specific overrides.
    default_cfg = config.get("defaults", {})
    cfg = deep_merge(default_cfg, profiles[profile_name])
    vprint(f"[Snakefile] Using profile: {profile_name}")
    vprint(f"[Snakefile] profile_name={profile_name}")
    vprint(f"[Snakefile] cfg.paths.unsup_root = {cfg.get('paths', {}).get('unsup_root', 'NOT SET')}")
    vprint(f"[Snakefile] cfg.analysis.cancer_type = {cfg.get('analysis', {}).get('cancer_type', 'NOT SET')}")

    # Updates LOGROOT to the validated, profile-specific log directory.
    LOGROOT = os.path.join(PIPE_ROOT, "logs", profile_name)
    if not IS_DRYRUN:
        os.makedirs(LOGROOT, exist_ok=True)
    vprint(f"[Snakefile] LOGROOT = {LOGROOT}")
else:
    # Falls back to flat config when no profiles block is present.
    cfg = config
    profile_name = config.get("profile", "default")


# Method helper: retrieves profile-aware graph layout settings for reproducible plotting parameters.
def graph_layout_section(plot_key):
    layout_cfg = config.get("graph_plot_layout")
    if layout_cfg is None:
        raise ValueError("Missing required config section: graph_plot_layout")

    if plot_key not in layout_cfg:
        raise ValueError(f"Missing graph_plot_layout section for plot: {plot_key}")

    plot_cfg = layout_cfg[plot_key]
    if "default" not in plot_cfg:
        raise ValueError(f"Missing graph_plot_layout.{plot_key}.default")

    profile_key = profile_name.lower()
    default_cfg = plot_cfg["default"]
    profile_cfg = plot_cfg.get(profile_key, {})

    merged = dict(default_cfg)
    merged.update(profile_cfg)
    return merged


# Method helper: reads one plotting parameter from the resolved graph layout section.
def graph_layout_value(plot_key, param_key):
    section = graph_layout_section(plot_key)
    if param_key not in section:
        raise ValueError(
            f"Missing graph layout parameter: graph_plot_layout.{plot_key}.{param_key} "
            f"for profile {profile_name.lower()}"
        )
    return section[param_key]


# =============================================================================
# PATH HELPERS AND GLOBAL CONSTANTS
# =============================================================================
# Stage role: materialises output namespaces and shared constants used across method stages.
# These constants translate profile configuration into the paths and method labels used by later rules.

MULTICOHORT_CFG = config.get("multicohort_cancer", {})
IS_MULTICOHORT_PROFILE = profile_name == "multicohort_cancer"
IS_PAN_CANCER_PROFILE = profile_name == "pan_cancer"
DECLARE_PAN_CANCER_RULES = IS_PAN_CANCER_PROFILE
MULTICOHORT_OUTDIR = os.path.normpath(
    MULTICOHORT_CFG.get(
        "outdir",
        os.path.join("results", "unsupervised", "multicohort_cancer"),
    )
)
MULTICOHORT_CONSENSUS_DIR = os.path.join(MULTICOHORT_OUTDIR, "consensus")
MULTICOHORT_FEATURE_SELECTION_DIR = os.path.join(
    MULTICOHORT_OUTDIR,
    "feature_selection_unsupervised",
)
MULTICOHORT_STUDY_DESIGN_DIR = os.path.join(MULTICOHORT_OUTDIR, "study_design")
MULTICOHORT_TUMOUR_NEIGHBOURHOODS_DIR = os.path.join(
    MULTICOHORT_OUTDIR,
    "tumour_neighbourhoods",
)
MULTICOHORT_VALIDATION_DIR = os.path.join(MULTICOHORT_OUTDIR, "validation")
DISABLE_PCA_EVERYWHERE = bool(MULTICOHORT_CFG.get("disable_pca_everywhere", False))
vprint(f"[Snakefile] DISABLE_PCA_EVERYWHERE = {DISABLE_PCA_EVERYWHERE}")

# =============================================================================
# FEATURE / DIRECTION CONSTANTS
# =============================================================================
# Stage role: defines feature-distance representations and validates direction naming before DAG construction.
#
# The feature-method universe, the distance universe and the per-method top-N
# values are scientific declarations owned by configuration. They are required
# here rather than defaulted: a fallback list embedded in orchestration is a
# second source of truth that can drift from config and silently substitute a
# different analysis. Missing keys fail at parse time, before any job runs.
def require_cfg_section(section, *, description):
    """Return a required mapping from the merged profile configuration."""
    value = cfg.get(section)
    if not isinstance(value, dict) or not value:
        raise ValueError(
            f"Missing required config section: {section} ({description}). "
            f"Declare it under defaults or profiles.{profile_name} in {CFGFILE_ABS}."
        )
    return value


def require_cfg_value(section_mapping, key, *, section, kind=None):
    """Return a required key from an already-resolved config section."""
    if key not in section_mapping:
        raise ValueError(
            f"Missing required config key: {section}.{key}. "
            f"Declare it under defaults or profiles.{profile_name} in {CFGFILE_ABS}."
        )
    value = section_mapping[key]
    if kind is not None and not isinstance(value, kind):
        raise ValueError(f"config {section}.{key} must be of type {kind.__name__}")
    return value


_FEATURE_SETS_CFG = require_cfg_section(
    "feature_sets", description="feature-method and distance universe"
)
FEATURE_METHODS = require_cfg_value(
    _FEATURE_SETS_CFG, "methods", section="feature_sets", kind=list
)
DISTANCES = require_cfg_value(
    _FEATURE_SETS_CFG, "distances", section="feature_sets", kind=list
)
if not FEATURE_METHODS:
    raise ValueError("config feature_sets.methods must be a non-empty list")
if not DISTANCES:
    raise ValueError("config feature_sets.distances must be a non-empty list")

# Per-method top-N gene counts. Every configured feature method must declare
# one; there is no built-in table to fall back to.
_FEATURE_SELECTION_CFG = require_cfg_section(
    "feature_selection", description="per-method top-N gene counts"
)
_topn_cfg = require_cfg_value(
    _FEATURE_SELECTION_CFG, "method_topn", section="feature_selection", kind=dict
)
_missing_topn = [m for m in FEATURE_METHODS if m not in _topn_cfg]
if _missing_topn:
    raise ValueError(
        "Missing feature_selection.method_topn entries for configured feature "
        f"method(s): {', '.join(_missing_topn)}"
    )
METHOD_TOPN = {m: int(_topn_cfg[m]) for m in FEATURE_METHODS}

def method_topn(m):
    """Return the configured top-N gene count for feature method m."""
    if m not in METHOD_TOPN:
        raise ValueError(
            f"Feature method '{m}' has no configured feature_selection.method_topn entry"
        )
    return METHOD_TOPN[m]

# Method helper: constructs feature-distance representation names from configured features and distance metrics.
def build_directions():
    """
    Constructs the full list of analysis directions as method × distance combinations.
    Config-explicit direction lists are authoritative for the active profile.
    If no explicit list is present, fall back to method × distance combinations.
    """
    feature_dirs = [f"{m}_{d}" for m in FEATURE_METHODS for d in DISTANCES]
    explicit = cfg.get("tumour_neighbourhoods", {}).get("directions")
    if explicit is not None:
        return explicit
    return feature_dirs


FEATURE_DISTANCE_DIRECTIONS = build_directions()

# Adds precomputed representation outputs to cross-direction p-consensus
# summaries without admitting them to the generic clustering, neighbourhood,
# similarity-graph, or support-threshold direction universes.
_P_CONS_ADDITIONAL_SUMMARY_DIRECTIONS = cfg.get(
    "tumour_neighbourhoods", {}
).get("additional_summary_directions", [])
if not isinstance(_P_CONS_ADDITIONAL_SUMMARY_DIRECTIONS, list):
    raise ValueError(
        "tumour_neighbourhoods.additional_summary_directions must be a list"
    )
P_CONS_SUMMARY_DIRECTIONS = list(FEATURE_DISTANCE_DIRECTIONS)
for _direction in _P_CONS_ADDITIONAL_SUMMARY_DIRECTIONS:
    _direction = str(_direction).strip()
    if not _direction:
        raise ValueError(
            "tumour_neighbourhoods.additional_summary_directions contains an empty value"
        )
    if _direction in P_CONS_SUMMARY_DIRECTIONS:
        raise ValueError(
            "Duplicate p-consensus summary direction: "
            f"{_direction}"
        )
    P_CONS_SUMMARY_DIRECTIONS.append(_direction)


def cfgrel(*keys):
    """
    Traverses the merged config by key path and returns the value as-is (relative).
    Used for rule output declarations to keep paths portable across environments.
    """
    d = cfg
    for k in keys:
        d = d[k]
    return d


def cfgget_path_rel(default_rel, *keys):
    """
    Returns the config value at the given key path as a relative path, falling back
    to default_rel if any key is absent. Used in rule output declarations.
    """
    d = cfg
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default_rel
        d = d[k]
    return d


def cfgget_path_abs(default_rel, *keys):
    """
    Returns the config value at the given key path as an absolute path, falling back
    to the absolute form of default_rel if any key is absent. Used in params/shell contexts.
    """
    return abspath(cfgget_path_rel(default_rel, *keys))


# Sets fail-fast shell execution flags and pins thread counts for linear algebra libraries
# to prevent unintended CPU over-subscription on shared HPC nodes.
shell.prefix(
    f"set -euo pipefail; export PIPELINE_ROOT={shlex.quote(PIPE_ROOT)}; "
    "export LC_ALL=C; "
    "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; "
)

# Resolves the unsupervised output root: relative for rule output declarations
# (Snakemake DAG portability), absolute for params/shell arguments (script compatibility).
UNSUP_BASE = os.path.normpath(cfgget_path_rel("results/unsupervised", "paths", "unsup_root"))
vprint(f"[Snakefile] UNSUP_BASE (from config) = {UNSUP_BASE}")

# Converts absolute config-sourced paths back to PIPE_ROOT-relative for output declarations.
if os.path.isabs(UNSUP_BASE):
    UNSUP_BASE = os.path.relpath(UNSUP_BASE, PIPE_ROOT)

if IS_MULTICOHORT_PROFILE:
    legacy_multicohort_unsup = os.path.join(MULTICOHORT_OUTDIR, "unsupervised")
    if os.path.normpath(UNSUP_BASE) == os.path.normpath(legacy_multicohort_unsup):
        UNSUP_BASE = MULTICOHORT_OUTDIR

# Ensures the profile name is embedded as the final directory component of the output root,
# so results from different profiles never collide under the same parent directory.
# These paths keep profile-specific outputs separated while preserving portable relative rule outputs.

# Append the profile only when it is not already present in the configured
# path. This avoids doubled roots such as results/unsupervised/brca/brca.
if profile_name not in os.path.normpath(UNSUP_BASE).split(os.sep):
    UNSUP_REL = os.path.join(UNSUP_BASE, profile_name)
else:
    UNSUP_REL = UNSUP_BASE

# Absolute counterpart used only in params, inputs, and script arguments.
UNSUP = os.path.join(PIPE_ROOT, UNSUP_REL)
vprint(f"[Snakefile] UNSUP_REL (final, relative) = {UNSUP_REL}")
vprint(f"[Snakefile] UNSUP (final, absolute) = {UNSUP}")

# =============================================================================
# DATA PREPARATION CONFIGURATION
# =============================================================================
# Stage role: resolves configured expression inputs and multicohort staging paths before data-preparation rules.
# Resolves the joint VST-normalised expression matrix (absolute; external input).
VST_JOINT = profile_vst_joint_abs(profile_name)
ACTIVE_VST_VALIDATION = os.path.join(
    UNSUP_REL, "inputs", "active_joint_vst.validated.ok"
)


# Rule: validate_profile_joint_vst
# Rejects placeholders, truncated objects, duplicated identifiers, and invalid
# values before any active cohort matrix can enter downstream analyses.
rule validate_profile_joint_vst:
    input:
        joint = VST_JOINT
    output:
        marker = ACTIVE_VST_VALIDATION
    params:
        script = os.path.join(SCRIPTS_DIR, "validate_active_expression_matrix.R")
    log:
        os.path.join(LOGROOT, "validate_profile_joint_vst.log")
    conda: CONDA_ENV_R_BASE
    shell:
        r'''
        mkdir -p "$(dirname "{output.marker}")" "$(dirname "{log}")"
        Rscript "{params.script}" \
          --input "{input.joint}" \
          --output "{output.marker}" \
          --min-genes 1000 \
          --min-samples 3 \
          > "{log}" 2>&1
        '''

# Multicohort input construction (only active when running the multicohort_cancer profile).
MC_PROFILES    = [str(p) for p in MULTICOHORT_CFG.get("profiles", [])]
MC_OUTDIR_REL  = MULTICOHORT_OUTDIR
MC_OUTDIR      = abspath(MC_OUTDIR_REL)
MC_INPUT_DIR   = os.path.join(MC_OUTDIR, "inputs")
MC_EXPR        = os.path.join(MC_INPUT_DIR, "joint_expr_matrix.rds")
MC_META        = os.path.join(MC_INPUT_DIR, "joint_metadata.tsv")

if IS_MULTICOHORT_PROFILE and MC_PROFILES:
    MC_VST_INPUTS    = [profile_vst_joint_abs(p) for p in MC_PROFILES]
    MC_PROFILES_ARG  = ",".join(MC_PROFILES)
    MC_VST_LIST_STR  = " ".join(shlex.quote(p) for p in MC_VST_INPUTS)
    if os.path.abspath(VST_JOINT) != os.path.abspath(MC_EXPR):
        raise ValueError(
            "multicohort profile expects paths.vst_joint_rds to match the multicohort inputs output"
        )

    # Rule: build_multicohort_joint_inputs
    # Method role: data-preparation rule that merges cohort-level VST matrices into the multicohort expression space.
    # Flow: configured profile VST matrices -> joint expression matrix and sample metadata.
    # Provides the multicohort_cancer input matrix for feature selection and pan-cancer aggregation.
    rule build_multicohort_joint_inputs:
        """Merge cohort-level VST matrices into the joint expression matrix and metadata."""
        input:
            vst_list = MC_VST_INPUTS
        output:
            joint_matrix = VST_JOINT,
            joint_meta   = MC_META,
        params:
            script       = os.path.join(SCRIPTS_DIR, "build_multicohort_joint_inputs.R"),
            pipe_root    = PIPE_ROOT,
            outdir       = MC_OUTDIR,
            profiles     = MC_PROFILES_ARG,
            vst_list_str = MC_VST_LIST_STR
        log:
            os.path.join(LOGROOT, "multicohort_joint_inputs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            mkdir -p "{params.outdir}/inputs"
            Rscript "{params.script}" \
              --pipe_root "{params.pipe_root}" \
              --outdir    "{params.outdir}" \
              --profiles  "{params.profiles}" \
              --vst_list  "{params.vst_list_str}" \
              > "{log}" 2>&1
            '''


# =============================================================================
# DATA PREPARATION
# =============================================================================
# Stage role: harmonises configured cohort expression matrices and metadata before feature selection.
# Multicohort joint-input construction and per-profile cell/tumour splits live
# here so later workflow sections consume materialised expression matrices.

# Resolves cell line and tumour VST expression matrices (absolute; external inputs).
# Falls back to vst_joint_rds for profiles (e.g. brca, nbl) that do not define
# separate cell_vst_rds / tumour_vst_rds keys.
CELL_VST   = cfgget_path_abs(cfgrel("paths", "vst_joint_rds"), "paths", "cell_vst_rds")
TUMOUR_VST = cfgget_path_abs(cfgrel("paths", "vst_joint_rds"), "paths", "tumour_vst_rds")

if profile_name in ("brca", "nbl", "rbl"):

    # Rule: split_profile_joint_vst
    # Method role: data-preparation rule that separates a profile joint VST object by sample type.
    # Flow: joint profile VST matrix -> cell-line and tumour VST matrices.
    # Provides input for HC/k-means clustering and tumour-neighbourhood input construction.
    rule split_profile_joint_vst:
        input:
            joint = VST_JOINT,
            validation = ACTIVE_VST_VALIDATION
        output:
            cell = CELL_VST,
            tumour = TUMOUR_VST
        params:
            script = os.path.join(SCRIPTS_DIR, "split_joint_vst_by_sample_type.R")
        log: os.path.join(LOGROOT, "split_joint_vst_by_sample_type.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "$(dirname "{output.cell}")" "$(dirname "{output.tumour}")"
            Rscript "{params.script}" \
              --joint_rds "{input.joint}" \
              --out_cell "{output.cell}" \
              --out_tumour "{output.tumour}" \
              > "{log}" 2>&1
            '''

# Feature list paths.
# _REL variant → used in rule output: declarations (relative, for portable DAG tracking).
# Non-suffixed absolute variant → used in rule input:, params:, and shell: references.
# Directory no longer encodes a single top-N value since each method may use a different value.
# File naming: genes_top{N}_{method}.txt  (N is method-specific).
FEATURESETS_DIR_REL = os.path.join(UNSUP_REL, "feature_selection_unsupervised", "feature_sets")
FEATURESETS_DIR     = os.path.join(UNSUP,     "feature_selection_unsupervised", "feature_sets")
FEATURELIST_FILES   = [
    os.path.join(FEATURESETS_DIR_REL, f"genes_top{method_topn(m)}_{m}.txt")
    for m in FEATURE_METHODS
    if m != "MX"
]

# Defines relative and absolute roots for HC/k-means clustering outputs.
FEATURE_METHOD_OUTDIR_REL = os.path.join(UNSUP_REL, "feature_selection_unsupervised", "featuresets")
FEATURE_METHOD_OUTDIR     = os.path.join(UNSUP,     "feature_selection_unsupervised", "featuresets")

# PIPELINE_TARGET remains initialised near rule all; final assembly is handled
# by build_pipeline_targets() after all optional section variables are defined.


# =============================================================================
# FEATURE SELECTION
# =============================================================================
# Stage role: derives unsupervised gene sets that define feature-distance representations for later stages.

# Rule: feature_selection_unsupervised
# Method role: analysis rule that derives unsupervised gene sets for each configured feature method.
# Flow: profile VST expression matrix -> method-specific ranked gene lists.
# Analysis role: defines the feature component of every feature-distance representation.
rule feature_selection_unsupervised:
    """
    Selects the top N genes from the joint VST expression matrix using multiple
    statistical methods (Variance, MAD, Entropy, PCA loadings, MX, HVG, etc.).
    Each method produces an independent ranked gene list that defines one
    feature direction propagated through all clustering rules that consume feature directions.
    The MX gene list is declared explicitly as mx_list so later rules
    can reference it as the canonical MX top-500 set without duplicating the path.
    """
    input:
        vst = VST_JOINT,
        validation = ACTIVE_VST_VALIDATION,
        cfg_file = CFGFILE_ABS
    output:
        # Relative paths for portable DAG tracking; R script writes via absolute --outdir.
        feature_lists = FEATURELIST_FILES,
        # Canonical MX gene list (MX method, top-500) — named output for rules that consume this output.
        mx_list = os.path.join(FEATURESETS_DIR_REL, f"genes_top{method_topn('MX')}_MX.txt")
    params:
        outdir         = lambda wc: os.path.join(UNSUP, "feature_selection_unsupervised"),
        # Comma-separated method:N pairs so the R script knows each method's top-N.
        method_topn    = ",".join(f"{m}:{n}" for m, n in METHOD_TOPN.items())
    # WGCNA::allowWGCNAThreads() requires nThreads >= 2, so this rule must
    # always run with at least 2 threads (regardless of executor defaults).
    threads: 2
    log: os.path.join(LOGROOT, "feature_selection_unsupervised.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p "{params.outdir}" "{FEATURESETS_DIR}"
        Rscript "{SCRIPTS_DIR}/feature_selection_unsupervised.R" \
          --config="{input.cfg_file}" \
          --profile="{profile_name}" \
          --outdir="{params.outdir}" \
          --vst_rds="{input.vst}" \
          --scripts_dir="{SCRIPTS_DIR}" \
          --method_topn="{params.method_topn}" \
          > "{log}" 2>&1
        '''


# =============================================================================
# HC/K-MEANS CLUSTERING INPUT CONSTRUCTION / CONFIGURATION
# =============================================================================
# Stage role: maps each feature-distance representation to filtered expression inputs and clustering parameters.

# Resolves the HC/k-means clustering output root.
# Config override (hclust_kmeans_root) must be a relative path; if absent,
# the default is derived from the profile-scoped UNSUP_REL directory.
_hclust_kmeans_cfg_raw = cfg["paths"].get("hclust_kmeans_root")
if _hclust_kmeans_cfg_raw and not os.path.isabs(_hclust_kmeans_cfg_raw):
    HCLUST_KMEANS_ROOT_REL = _hclust_kmeans_cfg_raw
else:
    HCLUST_KMEANS_ROOT_REL = os.path.join(UNSUP_REL, "hclust_kmeans")
HCLUST_KMEANS_ROOT     = HCLUST_KMEANS_ROOT_REL                  # Backwards-compatibility alias.
HCLUST_KMEANS_ROOT_ABS = os.path.join(PIPE_ROOT, HCLUST_KMEANS_ROOT_REL) if not os.path.isabs(HCLUST_KMEANS_ROOT_REL) else HCLUST_KMEANS_ROOT_REL
HCLUST_KMEANS_BASE_FUN = abspath(cfg["paths"].get("hclust_kmeans_base_functions_dir", cfg["paths"].get("base_functions_dir", os.path.join(UNSUP, "base_functions"))))

# HC/k-means clustering hyperparameters are scientific parameters owned by
# configuration, under the hclust_kmeans section. Every key is
# required, so a missing declaration fails at parse time instead of silently
# substituting a value embedded here. This k range is deliberately distinct from
# clustering.k_grid, which belongs to the ConsensusClusterPlus formulations.
#
# The HCLUST_KMEANS_ variable prefix below names the HC/k-means clustering
# stage and its hyperparameters.
HCLUST_KMEANS_SECTION = "hclust_kmeans"
HCLUST_KMEANS_CFG = require_cfg_section(
    HCLUST_KMEANS_SECTION,
    description="hierarchical-clustering and k-means hyperparameters",
)
if not HCLUST_KMEANS_CFG:
    raise ValueError(f"config {HCLUST_KMEANS_SECTION} must not be empty")

HCLUST_KMEANS_N_PCS    = int(require_cfg_value(HCLUST_KMEANS_CFG, "n_pcs", section=HCLUST_KMEANS_SECTION, kind=int))
HCLUST_KMEANS_MAX_K_HC = int(require_cfg_value(HCLUST_KMEANS_CFG, "max_k_hc", section=HCLUST_KMEANS_SECTION, kind=int))
HCLUST_KMEANS_K_MIN     = int(require_cfg_value(HCLUST_KMEANS_CFG, "k_min", section=HCLUST_KMEANS_SECTION, kind=int))
HCLUST_KMEANS_K_MAX     = int(require_cfg_value(HCLUST_KMEANS_CFG, "k_max", section=HCLUST_KMEANS_SECTION, kind=int))
HCLUST_KMEANS_SEED     = PIPELINE_SEED

if HCLUST_KMEANS_N_PCS < 1:
    raise ValueError(f"config {HCLUST_KMEANS_SECTION}.n_pcs must be >= 1")
if HCLUST_KMEANS_MAX_K_HC < 2:
    raise ValueError(f"config {HCLUST_KMEANS_SECTION}.max_k_hc must be >= 2")
if HCLUST_KMEANS_K_MIN < 2:
    raise ValueError(f"config {HCLUST_KMEANS_SECTION}.k_min must be >= 2")
if HCLUST_KMEANS_K_MAX < HCLUST_KMEANS_K_MIN:
    raise ValueError(
        f"config {HCLUST_KMEANS_SECTION}.k_max must be >= {HCLUST_KMEANS_SECTION}.k_min"
    )

# Reads and filters the consensus direction list independently of FEATURE_DISTANCE_DIRECTIONS.
CONS_DIRECTIONS = cfg.get("tumour_neighbourhoods", {}).get("directions", FEATURE_DISTANCE_DIRECTIONS)

# Builds wildcard constraint patterns from configured directions.
FEATURE_DISTANCE_DIRECTION_PATTERN = "|".join(map(re.escape, FEATURE_DISTANCE_DIRECTIONS))
FEATURE_DISTANCE_EUC_DIRECTIONS    = [d for d in FEATURE_DISTANCE_DIRECTIONS if d.endswith("_euc")]
# Falls back to an impossible pattern ('a^') when no Euclidean directions exist,
# preventing k-means rules from matching any wildcard.
FEATURE_DISTANCE_EUC_PATTERN       = "|".join(map(re.escape, FEATURE_DISTANCE_EUC_DIRECTIONS)) if FEATURE_DISTANCE_EUC_DIRECTIONS else r"a^"

CONS_DIRECTION_PATTERN = "|".join(map(re.escape, CONS_DIRECTIONS))
CONS_EUC_DIRECTIONS    = [d for d in CONS_DIRECTIONS if d.endswith("_euc")]
CONS_EUC_PATTERN       = "|".join(map(re.escape, CONS_EUC_DIRECTIONS)) if CONS_EUC_DIRECTIONS else r"a^"
PAN_CANCER_FEATURE_SET_NAME = "PanCancerFeatureSet"
PAN_CANCER_FEATURE_SET_GENE_LIST = (
    config.get("defaults", {})
    .get("marker_postprocessing", {})
    .get("pan_cancer", {})
    .get("final_features_clean")
)
EXTERNAL_GENE_LIST_FEATURES = {PAN_CANCER_FEATURE_SET_NAME} if PAN_CANCER_FEATURE_SET_GENE_LIST else set()


# Method helper: derives the majority edge-support threshold from the
# active feature-distance representation graphs for this profile.
def derive_majority_support_threshold(active_representations):
    """
    Return m = max(2, floor(|R| / 2) + 1), where |R| is the number of active
    feature-distance representations contributing to the current graph.
    """
    representations = [str(rep).strip() for rep in active_representations if str(rep).strip()]
    if len(representations) != len(set(representations)):
        duplicates = sorted({rep for rep in representations if representations.count(rep) > 1})
        raise ValueError(
            "[Support threshold] Active representation identifiers must be unique; "
            f"duplicates={duplicates}"
        )
    if len(representations) < 2:
        raise ValueError(
            "[Support threshold] At least two active representations are required "
            "for majority-threshold graph construction."
        )
    threshold = max(2, len(representations) // 2 + 1)
    print(
        "[Support threshold] "
        f"profile={profile_name} active_representations={len(representations)} "
        f"majority_support_threshold={threshold}"
    )
    return threshold


SIMILARITY_CONSENSUS_MIN_SUPPORT = derive_majority_support_threshold(CONS_DIRECTIONS)

# Defines the six HC and six k-means clustering kind identifiers.
# PCA-prefixed kinds are optionally removed for pan-cancer runs where
# cross-cancer PCA spaces are not comparable.
HC_KINDS = ["pca_hc_cell", "pca_hc_tumour", "pca_hc_cell_tumour", "hc_cell", "hc_tumour", "hc_cell_tumour"]
KM_KINDS = ["pca_kmeans_cell", "pca_kmeans_tumour", "pca_kmeans_cell_tumour", "kmeans_cell", "kmeans_tumour", "kmeans_cell_tumour"]

if DISABLE_PCA_EVERYWHERE:
    HC_KINDS = [k for k in HC_KINDS if not k.startswith("pca_")]
    KM_KINDS = [k for k in KM_KINDS if not k.startswith("pca_")]
    vprint(f"[Snakefile] Filtered PCA kinds: HC_KINDS={HC_KINDS}, KM_KINDS={KM_KINDS}")

# JOINT cell-line + tumour HC/k-means clustering kinds. These are the only
# HC/k-means outputs eligible for tumour-neighbourhood construction and the
# p_consensus recurrence fraction; cell-only and tumour-only kinds never
# contribute. Derived from HC_KINDS/KM_KINDS so the PCA filter above applies.
HCLUST_JOINT_KINDS = [k for k in HC_KINDS if k.endswith("_cell_tumour")]
KMEANS_JOINT_KINDS = [k for k in KM_KINDS if k.endswith("_cell_tumour")]


# Method helper: extracts the feature component from a feature-distance representation name.
def dir_to_feature(direction):
    """Extracts the feature method name from a direction string (e.g., 'Variance_euc' → 'Variance')."""
    if direction.endswith("_euc") or direction.endswith("_corr"):
        return direction.rsplit("_", 1)[0]
    return direction


# Method helper: extracts the distance metric from a feature-distance representation name.
def dir_to_dist(direction):
    """Maps a direction suffix to its R-compatible distance method string."""
    return "correlation" if direction.endswith("corr") else "euclidean"


# Method helper: maps a representation to the gene-list file consumed by clustering input construction.
def dir_to_gene_list(direction):
    """
    Returns the gene-list path appropriate for a given direction. Feature-selection
    outputs use relative paths that match their producer rule; configured external
    feature sets use their configured absolute paths.

    Raises ValueError if the direction contains an unknown feature method.
    """
    feature = dir_to_feature(direction)

    if feature in EXTERNAL_GENE_LIST_FEATURES:
        # External feature-set directions reuse the configured marker-derived
        # gene list rather than a feature_selection_unsupervised output. The
        # configured value is returned verbatim so that it matches the producer
        # rule (construct_pan_cancer_feature_panel) output declaration exactly;
        # rewriting a relative configured path to an absolute one would make
        # Snakemake treat it as a distinct, unproducible file.
        return PAN_CANCER_FEATURE_SET_GENE_LIST
    elif feature in FEATURE_METHODS:
        # Feature-set directions: filename encodes the method-specific top-N.
        # Uses FEATURESETS_DIR_REL to match the producer rule output.
        return os.path.join(FEATURESETS_DIR_REL, f"genes_top{method_topn(feature)}_{feature}.txt")
    else:
        # Unknown feature: fail loudly to prevent silent reuse of a default gene list.
        raise ValueError(
            f"Unknown feature '{feature}' extracted from direction '{direction}'. "
            f"Known features: {sorted(set(FEATURE_METHODS) | EXTERNAL_GENE_LIST_FEATURES)}"
        )


# Method helper: maps a representation to the script-level gene-set argument used by dependent rules.
def dir_to_geneset(direction):
    """Returns a human-readable gene set label for annotation and reporting purposes."""
    feature = dir_to_feature(direction)
    if feature in EXTERNAL_GENE_LIST_FEATURES:
        return feature
    return f"{feature}_top{method_topn(feature)}"


# Method helper: validates representation naming against the expected feature-list path at parse time.
def validate_direction_to_gene_list(direction):
    """
    Validates that a direction string correctly resolves to its gene-list file.
    Three checks:
      1. The feature method parsed from the direction is recognized.
      2. The resolved filename encodes the correct method name.
      3. Feature-selection output paths are relative, matching the producer rule
         output style and preventing absolute-vs-relative DAG mismatches.

    Runs at Snakefile parse time to catch bugs before any job executes.
    """
    feature = dir_to_feature(direction)

    # 1. Check that the extracted feature is recognized.
    known_features = set(FEATURE_METHODS) | EXTERNAL_GENE_LIST_FEATURES
    assert feature in known_features, \
        f"Unknown feature '{feature}' extracted from direction '{direction}'. " \
        f"Known features: {sorted(known_features)}"

    # 2. Verify the gene-list file path is well-formed and method-specific.
    gene_list_path = dir_to_gene_list(direction)
    filename = os.path.basename(gene_list_path)

    # For feature-selection directions, verify filename format.
    if feature not in EXTERNAL_GENE_LIST_FEATURES:
        if not filename.startswith("genes_top"):
            raise AssertionError(
                f"Gene list file '{gene_list_path}' has unexpected name format. "
                f"Expected 'genes_top<N>_<method>.txt' pattern."
            )

    # Verify that the filename method matches the direction-derived feature.
    if feature not in EXTERNAL_GENE_LIST_FEATURES:
        if not filename.endswith(f"_{feature}.txt"):
            raise AssertionError(
                f"Gene list for direction '{direction}' (method '{feature}') "
                f"resolved to file '{filename}', which does not end with '_{feature}.txt'. "
                f"This indicates a misalignment in the direction-to-gene-list mapping."
            )

    # 3. Feature-selection paths must be relative to match the producer rule
    #    (feature_selection_unsupervised) output declarations.
    #    Absolute paths would cause Snakemake MissingInputException.
    if feature not in EXTERNAL_GENE_LIST_FEATURES and os.path.isabs(gene_list_path):
        raise AssertionError(
            f"Gene list path for direction '{direction}' is absolute: {gene_list_path}\n"
            f"The producer rule (feature_selection_unsupervised) declares outputs as "
            f"relative paths under FEATURESETS_DIR_REL={FEATURESETS_DIR_REL}.\n"
            f"Using an absolute path here will cause a MissingInputException."
        )

    return gene_list_path


# ---------------------------------------------------------------------------
# DIRECTION → GENE-LIST VALIDATION (runs at Snakefile parse time)
# Catches mapping bugs early: every declared direction must resolve to a
# gene-list file whose filename encodes the correct feature method.
# ---------------------------------------------------------------------------
for _d in FEATURE_DISTANCE_DIRECTIONS:
    validate_direction_to_gene_list(_d)
vprint(f"[Snakefile] Direction→gene-list validation passed for {len(FEATURE_DISTANCE_DIRECTIONS)} directions")


# ---------------------------------------------------------------------------
# PAN-CANCER FEATURE-PANEL REPRESENTATION OWNERSHIP
# ---------------------------------------------------------------------------
# PanCancerFeatureSet_* representations are not produced by
# feature_selection_unsupervised: their gene list is the marker-derived
# pan-cancer feature panel, a product of the pan_cancer marker-postprocessing
# stage. Any profile whose configured representation universe declares such a
# representation must therefore build the panel first, so the panel becomes a
# DAG dependency of the representation's clustering, neighbourhood, graph and
# resolution stages rather than a file that merely happens to be on disk.
#
# The representation universe is read from the effective configuration only;
# it is never inferred from directory contents.
PAN_CANCER_FEATURE_SET_DIRECTIONS = sorted(
    {
        d
        for d in list(FEATURE_DISTANCE_DIRECTIONS) + list(CONS_DIRECTIONS)
        if dir_to_feature(d) == PAN_CANCER_FEATURE_SET_NAME
    }
)
REQUIRES_PAN_CANCER_FEATURE_PANEL = bool(PAN_CANCER_FEATURE_SET_DIRECTIONS)

if REQUIRES_PAN_CANCER_FEATURE_PANEL:
    if not PAN_CANCER_FEATURE_SET_GENE_LIST:
        raise ValueError(
            f"Profile '{profile_name}' declares representation(s) "
            f"{PAN_CANCER_FEATURE_SET_DIRECTIONS}, but "
            "defaults.marker_postprocessing.pan_cancer.final_features_clean is "
            "not configured, so the required pan-cancer feature panel cannot be "
            "built."
        )
    if os.path.isabs(PAN_CANCER_FEATURE_SET_GENE_LIST):
        raise ValueError(
            "defaults.marker_postprocessing.pan_cancer.final_features_clean must "
            "be a repository-relative path so that the PanCancerFeatureSet "
            "representations depend on the rule that produces it; got "
            f"{PAN_CANCER_FEATURE_SET_GENE_LIST}"
        )
    vprint(
        "[Snakefile] Pan-cancer feature panel required by "
        f"{len(PAN_CANCER_FEATURE_SET_DIRECTIONS)} representation(s): "
        f"{PAN_CANCER_FEATURE_SET_DIRECTIONS}"
    )


def hclust_kmeans_outdir_dir(direction, kind):
    """Returns the absolute output directory path for a given direction–kind combination."""
    return os.path.join(HCLUST_KMEANS_ROOT_ABS, direction, kind)


# =============================================================================
# HC/K-MEANS CLUSTERING
# =============================================================================
# Stage role: applies hierarchical clustering and k-means to each configured
# feature-distance representation in expression space and PCA-reduced space.
# JOINT (cell-line + tumour) outputs feed tumour-neighbourhood construction and
# the p_consensus recurrence fraction alongside the parallel ConsensusClusterPlus
# formulations; cell-only and tumour-only outputs never contribute to p_consensus.
# Configuration for this stage lives under hclust_kmeans, its
# orchestration variables use the HCLUST_KMEANS_ prefix, and its
# outputs are written under the hclust_kmeans/ namespace.

# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_build_inputs
#
# Method role: constructs representation-specific cell-line and tumour expression matrices.
# Flow: VST matrices plus the configured gene list -> filtered clustering inputs.
# -----------------------------------------------------------------------------
rule hclust_kmeans_build_inputs:
    """
    Subsets the cell line and tumour VST expression matrices to the gene list
    defined by the given direction, producing filtered RDS objects used as
    input to all HC/k-means clustering rules for that direction.
    """
    input:
        cell   = CELL_VST,
        tumour = TUMOUR_VST,
        genes  = lambda wc: dir_to_gene_list(wc.direction)
    output:
        cell_filt   = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour_filt = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    log: os.path.join(LOGROOT, "hclust_kmeans_build_inputs_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "$(dirname "{output.cell_filt}")" "$(dirname "{output.tumour_filt}")"
        Rscript "{SCRIPTS_DIR}/build_hclust_kmeans_matrices.R" \
          --cell_rds   "{input.cell}" \
          --tumour_rds "{input.tumour}" \
          --genes      "{input.genes}" \
          --out_cell   "{output.cell_filt}" \
          --out_tumour "{output.tumour_filt}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_hc_cell
#
# Method role: clusters cell-line samples after PCA using hierarchical clustering.
# Flow: filtered cell-line matrix -> PCA hierarchical-clustering assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_hc_cell:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the cell line expression matrix. Optimal k is selected by silhouette
    or gap statistic within [2, max_k]. Produces a cell-only cluster assignment
    RDS; cell-only outputs never feed tumour-neighbourhood construction or
    p_consensus, which consume JOINT cell-line + tumour outputs only.
    """
    input:
        cell = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_hc_cell", "pca_hc_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_hc_cell"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_hc_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_hc_cell.R" \
          --cell_rds    "{input.cell}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       CELL \
          --n_pcs       {params.n_pcs} \
          --max_k       {params.max_k} \
          --dist_method {params.dist} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_hc_tumour
#
# Method role: clusters tumour samples after PCA within the configured feature-distance representation.
# Flow: filtered tumour expression matrix -> tumour-only PCA hierarchical-clustering assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_hc_tumour:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the tumour expression matrix. Mirrors hclust_kmeans_pca_hc_cell but
    operates on patient tumour samples independently to capture tumour-intrinsic
    structure.
    """
    input:
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_hc_tumour", "pca_hc_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_hc_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_hc_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_hc_tumour.R" \
          --tumour_rds  "{input.tumour}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       TUMOUR \
          --n_pcs       {params.n_pcs} \
          --max_k       {params.max_k} \
          --dist_method {params.dist} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_hc_joint
#
# Method role: clusters combined cell-line and tumour samples after PCA using hierarchical clustering.
# Flow: filtered joint expression matrix -> PCA joint-clustering assignment.
# Joint clustering evaluates whether cell-line and tumour samples occupy neighbouring transcriptomic regions.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_hc_joint:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the concatenated cell line and tumour expression matrix. Joint embedding
    captures shared transcriptional structure across both sample types,
    enabling direct cross-space cluster comparison.
    """
    input:
        cell   = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_hc_cell_tumour", "pca_hc_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_hc_cell_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_hc_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_hc_cell_tumour.R" \
          --cell_rds     "{input.cell}" \
          --tumour_rds   "{input.tumour}" \
          --outdir       "{params.outdir}" \
          --base_dir     "{params.base_dir}" \
          --cell_label   CELL \
          --tumour_label TUMOUR \
          --n_pcs        {params.n_pcs} \
          --max_k        {params.max_k} \
          --dist_method  {params.dist} \
          --seed         {params.seed} \
          --cluster_rds  "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_kmeans_cell
#
# Method role: clusters PCA-reduced cell-line samples using k-means for one representation.
# Flow: filtered cell-line matrix -> PCA k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_kmeans_cell:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the cell line expression matrix. Constrained to Euclidean directions only,
    as k-means minimises Euclidean inertia and is undefined under correlation distance.
    """
    input:
        cell = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_kmeans_cell", "pca_kmeans_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_kmeans_cell"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_kmeans_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_kmeans_cell.R" \
          --cell_rds    "{input.cell}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       CELL \
          --n_pcs       {params.n_pcs} \
          --k_min       {params.kmin} \
          --k_max       {params.kmax} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_kmeans_tumour
#
# Method role: clusters PCA-reduced tumour samples using k-means for one representation.
# Flow: filtered tumour matrix -> PCA k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_kmeans_tumour:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the tumour expression matrix. Euclidean directions only; mirrors
    hclust_kmeans_pca_kmeans_cell but operates on patient samples independently.
    """
    input:
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_kmeans_tumour", "pca_kmeans_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_kmeans_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_kmeans_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_kmeans_tumour.R" \
          --tumour_rds  "{input.tumour}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       TUMOUR \
          --n_pcs       {params.n_pcs} \
          --k_min       {params.kmin} \
          --k_max       {params.kmax} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_pca_kmeans_joint
#
# Method role: clusters PCA-reduced joint cell-line and tumour samples using k-means.
# Flow: filtered joint expression matrix -> PCA joint k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_pca_kmeans_joint:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the joint cell line and tumour expression matrix. Euclidean directions only.
    Joint embedding allows k-means to assign both sample types to shared clusters.
    """
    input:
        cell   = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "pca_kmeans_cell_tumour", "pca_kmeans_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "pca_kmeans_cell_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        n_pcs    = HCLUST_KMEANS_N_PCS,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_pca_kmeans_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/pca_kmeans_cell_tumour.R" \
          --cell_rds     "{input.cell}" \
          --tumour_rds   "{input.tumour}" \
          --outdir       "{params.outdir}" \
          --base_dir     "{params.base_dir}" \
          --cell_label   CELL \
          --tumour_label TUMOUR \
          --n_pcs        {params.n_pcs} \
          --k_min        {params.kmin} \
          --k_max        {params.kmax} \
          --seed         {params.seed} \
          --cluster_rds  "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_hc_cell
#
# Method role: clusters filtered cell-line samples directly using hierarchical clustering.
# Flow: filtered cell-line matrix -> raw-space hierarchical-clustering assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_hc_cell:
    """
    Applies hierarchical clustering directly to the raw gene-subsetted cell line
    expression matrix without PCA pre-processing. Distance metric is determined
    by the direction suffix (Euclidean or correlation). Provides a non-PCA
    baseline for comparison with hclust_kmeans_pca_hc_cell under the same
    feature set.
    """
    input:
        cell = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "hc_cell", "hc_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "hc_cell"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_hc_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/hc_cell.R" \
          --cell_rds    "{input.cell}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       CELL \
          --max_k       {params.max_k} \
          --dist_method {params.dist} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_hc_tumour
#
# Method role: clusters filtered tumour samples directly using hierarchical clustering.
# Flow: filtered tumour matrix -> raw-space hierarchical-clustering assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_hc_tumour:
    """
    Applies hierarchical clustering directly to the raw gene-subsetted tumour
    expression matrix without PCA pre-processing. Provides a non-PCA baseline
    for tumour-intrinsic cluster structure under each feature direction.
    """
    input:
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "hc_tumour", "hc_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "hc_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_hc_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/hc_tumour.R" \
          --tumour_rds  "{input.tumour}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       TUMOUR \
          --max_k       {params.max_k} \
          --dist_method {params.dist} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_hc_joint
#
# Method role: clusters filtered joint cell-line and tumour samples using hierarchical clustering.
# Flow: filtered joint expression matrix -> raw-space joint hierarchical-clustering assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_hc_joint:
    """
    Applies hierarchical clustering to the joint cell line and tumour expression
    matrix without PCA pre-processing. Supports both Euclidean and correlation
    distances, enabling direct assessment of co-clustering fidelity across
    sample types in the raw feature space.
    """
    input:
        cell   = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "hc_cell_tumour", "hc_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "hc_cell_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        max_k    = HCLUST_KMEANS_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_hc_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/hc_cell_tumour.R" \
          --cell_rds     "{input.cell}" \
          --tumour_rds   "{input.tumour}" \
          --outdir       "{params.outdir}" \
          --base_dir     "{params.base_dir}" \
          --cell_label   CELL \
          --tumour_label TUMOUR \
          --max_k        {params.max_k} \
          --dist_method  {params.dist} \
          --seed         {params.seed} \
          --cluster_rds  "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_kmeans_cell
#
# Method role: clusters filtered cell-line samples directly using k-means.
# Flow: filtered cell-line matrix -> raw-space k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_kmeans_cell:
    """
    Applies k-means clustering directly to the raw cell line expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    Serves as a non-PCA counterpart to hclust_kmeans_pca_kmeans_cell for
    sensitivity assessment.
    """
    input:
        cell = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "kmeans_cell", "kmeans_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "kmeans_cell"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_kmeans_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/kmeans_cell.R" \
          --cell_rds    "{input.cell}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       CELL \
          --k_min       {params.kmin} \
          --k_max       {params.kmax} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_kmeans_tumour
#
# Method role: clusters filtered tumour samples directly using k-means.
# Flow: filtered tumour matrix -> raw-space k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_kmeans_tumour:
    """
    Applies k-means clustering directly to the raw tumour expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    """
    input:
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "kmeans_tumour", "kmeans_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "kmeans_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_kmeans_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/kmeans_tumour.R" \
          --tumour_rds  "{input.tumour}" \
          --outdir      "{params.outdir}" \
          --base_dir    "{params.base_dir}" \
          --label       TUMOUR \
          --k_min       {params.kmin} \
          --k_max       {params.kmax} \
          --seed        {params.seed} \
          --cluster_rds "{output.clusters}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: hclust_kmeans_kmeans_joint
#
# Method role: clusters filtered joint cell-line and tumour samples directly using k-means.
# Flow: filtered joint expression matrix -> raw-space joint k-means assignment.
# -----------------------------------------------------------------------------
rule hclust_kmeans_kmeans_joint:
    """
    Applies k-means clustering to the joint cell line and tumour expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    Joint clustering provides a raw-space baseline for cross-type co-assignment.
    """
    input:
        cell   = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "kmeans_cell_tumour", "kmeans_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: hclust_kmeans_outdir_dir(wc.direction, "kmeans_cell_tumour"),
        base_dir = HCLUST_KMEANS_BASE_FUN,
        kmin     = HCLUST_KMEANS_K_MIN,
        kmax     = HCLUST_KMEANS_K_MAX,
        seed     = HCLUST_KMEANS_SEED
    log: os.path.join(LOGROOT, "hclust_kmeans_{direction}_kmeans_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_EUC_PATTERN
    shell:
        r'''
        mkdir -p "{params.outdir}"
        Rscript "{SCRIPTS_DIR}/kmeans_cell_tumour.R" \
          --cell_rds     "{input.cell}" \
          --tumour_rds   "{input.tumour}" \
          --outdir       "{params.outdir}" \
          --base_dir     "{params.base_dir}" \
          --cell_label   CELL \
          --tumour_label TUMOUR \
          --k_min        {params.kmin} \
          --k_max        {params.kmax} \
          --seed         {params.seed} \
          --cluster_rds  "{output.clusters}" \
          > "{log}" 2>&1
        '''

# Rule: hclust_kmeans_all
# Method role: collector target that requires all HC/k-means clustering outputs for configured directions.
# Flow: expected clustering output paths -> Snakemake dependency aggregation only.
# Stage output: provides a target marker without running a new analysis.
rule hclust_kmeans_all:
    """
    Convenience aggregation rule that collects all HC/k-means cluster outputs
    across every direction. Allows the full HC/k-means clustering suite to be
    run independently via 'snakemake hclust_kmeans_all'. The JOINT
    cell-line + tumour outputs are also required inputs of the
    tumour-neighbourhood rules, which pull them on demand.
    """
    input:
        expand(os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=FEATURE_DISTANCE_DIRECTIONS, kind=HC_KINDS),
        expand(os.path.join(HCLUST_KMEANS_ROOT_REL, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=FEATURE_DISTANCE_EUC_DIRECTIONS, kind=KM_KINDS)


# =============================================================================
# CONSENSUS CLUSTERING
# =============================================================================
# Stage role: aggregates clustering evidence within each representation before tumour-neighbourhood analysis.

# Defines consensus output root and clustering kind identifiers.
# PCA-containing kinds are stripped when DISABLE_PCA_EVERYWHERE is set.
CONS_ROOT_REL = os.path.join(UNSUP_REL, "consensus")
CONS_ROOT     = os.path.join(UNSUP,     "consensus")
CONS_KM_KINDS = ["ccp_kmeans_expr_cell_tumour", "ccp_kmeans_pca_cell_tumour", "ccp_kmeans_expr_cell", "ccp_kmeans_pca_cell", "ccp_kmeans_expr_tumour", "ccp_kmeans_pca_tumour"]
CONS_HC_KINDS = ["ccp_hc_expr_cell_tumour", "ccp_hc_pca_cell_tumour", "ccp_hc_expr_cell", "ccp_hc_pca_cell", "ccp_hc_expr_tumour", "ccp_hc_pca_tumour"]

if DISABLE_PCA_EVERYWHERE:
    CONS_KM_KINDS = [k for k in CONS_KM_KINDS if "_pca_" not in k]
    CONS_HC_KINDS = [k for k in CONS_HC_KINDS if "_pca_" not in k]
    vprint(f"[Snakefile] Filtered PCA consensus kinds: CONS_KM_KINDS={CONS_KM_KINDS}, CONS_HC_KINDS={CONS_HC_KINDS}")

# Builds a single regex pattern covering all consensus kind identifiers for wildcard validation.
CONS_KIND_PATTERN = "|".join(map(re.escape, CONS_KM_KINDS + CONS_HC_KINDS))


# ---------------------------------------------------------------------------
# CONFIGURED CLUSTERING-FORMULATION OWNERSHIP GATE (parse time)
# ---------------------------------------------------------------------------
# patient_referenced_graph.clustering_methods_by_distance owns which clustering
# formulations contribute to tumour-neighbourhood construction and to the
# p_consensus denominator. The lists above own only which clustering products
# the workflow builds (implementation dispatch, including the cell-only and
# tumour-only products that never contribute).
#
# This gate ties the two together at parse time: every configured formulation
# identifier must map onto a clustering product this DAG actually declares, for
# the distance component it is declared under. It also rejects the retired
# identifier vocabulary (AGN_*, CCP_HC_*, CCP_KM_*) explicitly, so a stale
# configuration fails loudly instead of resolving to a directory that no rule
# produces.
_RETIRED_FORMULATION_PREFIXES = ("AGN_", "CCP_HC_", "CCP_KM_")


def formulation_id_to_kind(method_id):
    """
    Map a configured clustering-formulation identifier to the (root, kind)
    clustering product the workflow declares for it.

    Mirrors R/base_functions/tumour_nh_io.R::nh_method_directory so that
    orchestration and implementation agree on one vocabulary translation.
    """
    for prefix in _RETIRED_FORMULATION_PREFIXES:
        if method_id.startswith(prefix):
            raise ValueError(
                f"Retired clustering formulation identifier '{method_id}' in "
                "patient_referenced_graph.clustering_methods_by_distance. Use the "
                "HCLUST_*/KMEANS_*/CCP_HCLUST_*/CCP_KMEANS_* vocabulary."
            )
    is_consensus = method_id.startswith("CCP_")
    remainder = method_id[len("CCP_"):] if is_consensus else method_id
    parts = remainder.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unrecognised clustering formulation identifier: {method_id}")
    algorithm, space = parts[0], parts[1]
    scope = "_".join(parts[2:])
    if scope != "cell_tumour":
        raise ValueError(
            "Only JOINT cell-line + tumour formulations are eligible; got: "
            f"{method_id}"
        )
    if algorithm not in ("HCLUST", "KMEANS"):
        raise ValueError(
            f"Unknown clustering algorithm in identifier '{method_id}'; expected "
            "HCLUST or KMEANS."
        )
    if space not in ("expr", "pca"):
        raise ValueError(
            f"Unknown feature space in identifier '{method_id}'; expected expr or pca."
        )
    stem = "hc" if algorithm == "HCLUST" else "kmeans"
    if is_consensus:
        return ("consensus", f"ccp_{stem}_{space}_cell_tumour")
    prefix = "pca_" if space == "pca" else ""
    return ("hclust_kmeans", f"{prefix}{stem}_cell_tumour")


_DECLARED_CLUSTERING_PRODUCTS = {
    "consensus": set(CONS_HC_KINDS) | set(CONS_KM_KINDS),
    "hclust_kmeans": set(HCLUST_JOINT_KINDS) | set(KMEANS_JOINT_KINDS),
}
_KMEANS_ONLY_PRODUCTS = {
    ("consensus", k) for k in CONS_KM_KINDS
} | {
    ("hclust_kmeans", k) for k in KMEANS_JOINT_KINDS
}

_CLUSTERING_METHODS_BY_DISTANCE = (
    cfg.get("patient_referenced_graph", {}).get("clustering_methods_by_distance")
)
if not isinstance(_CLUSTERING_METHODS_BY_DISTANCE, dict) or not _CLUSTERING_METHODS_BY_DISTANCE:
    raise ValueError(
        "Missing required config section: "
        "patient_referenced_graph.clustering_methods_by_distance"
    )

for _distance in sorted({d.rsplit("_", 1)[1] for d in CONS_DIRECTIONS}):
    _declared = _CLUSTERING_METHODS_BY_DISTANCE.get(_distance)
    if not _declared:
        raise ValueError(
            "patient_referenced_graph.clustering_methods_by_distance declares no "
            f"formulations for distance component '{_distance}', which is used by "
            f"profile '{profile_name}'."
        )
    if len(set(_declared)) != len(_declared):
        raise ValueError(
            "Duplicate clustering formulation identifier(s) for distance "
            f"'{_distance}'"
        )
    for _method_id in _declared:
        _root, _kind = formulation_id_to_kind(str(_method_id))
        if _kind not in _DECLARED_CLUSTERING_PRODUCTS[_root]:
            raise ValueError(
                f"Configured clustering formulation '{_method_id}' (distance "
                f"'{_distance}') maps to clustering product {_root}/{_kind}, which "
                "this workflow does not declare. Configuration and DAG products "
                "must describe the same formulation universe."
            )
        if _distance != "euc" and (_root, _kind) in _KMEANS_ONLY_PRODUCTS:
            raise ValueError(
                f"k-means formulation '{_method_id}' is declared for distance "
                f"'{_distance}'; k-means minimises Euclidean inertia and is only "
                "declared for Euclidean representations."
            )

vprint(
    "[Snakefile] Clustering-formulation ownership gate passed for distances: "
    f"{sorted(_CLUSTERING_METHODS_BY_DISTANCE)}"
)


def cons_outdir(direction, kind):
    """Returns the absolute output directory for a consensus direction–kind combination."""
    return os.path.join(CONS_ROOT, direction, kind)


def kind_to_mode(kind):
    """Extracts the feature mode (pca or expr) from a consensus kind identifier."""
    return "pca" if "_pca_" in kind else "expr"


def kind_to_alg(kind):
    """Extracts the clustering algorithm (hc or km) from a consensus kind identifier."""
    return "hc" if kind.startswith("ccp_hc_") else "km"


# -----------------------------------------------------------------------------
# RULE: consensus_cluster_ccp
#
# Method role: computes resampling-based consensus clustering for one representation.
# Flow: profile VST matrices plus the feature list -> consensus-clustering assignment.
# Consensus clustering evaluates grouping stability across resampled clustering runs.
# It is a parallel clustering formulation: it does NOT consume HC/k-means
# clustering outputs; the two branches meet at the p_consensus stage via their
# method-specific tumour neighbourhoods.
# -----------------------------------------------------------------------------
rule consensus_cluster_ccp:
    """
    Runs ConsensusClusterPlus on the joint cell line and tumour expression matrix
    for a given direction and kind. Resampling-based consensus stabilises cluster
    assignments that are sensitive to initialisation in single-run methods.
    Mode (pca/expr) and algorithm (hc/km) are inferred from the kind wildcard,
    allowing a single rule to cover all twelve consensus variants per direction.
    """
    input:
        config   = CFGFILE_ABS,
        # Resolves direction-specific feature-list file using dir_to_gene_list helper.
        # Each direction (e.g. Variance_euc, PCA_corr, MX_euc) maps to its corresponding
        # method's gene list, ensuring clustering uses the correct feature set.
        gene_list  = lambda wc: dir_to_gene_list(wc.direction),
        # Declared so Snakemake tracks split_profile_joint_vst as upstream.
        # The R script resolves these paths from config; inputs wire the DAG.
        cell_vst   = CELL_VST,
        tumour_vst = TUMOUR_VST
    output:
        cluster_rds = os.path.join(CONS_ROOT, "{direction}", "{kind}", "{kind}_clusters_optimal.rds")
    params:
        outdir = lambda wc: cons_outdir(wc.direction, wc.kind),
        mode   = lambda wc: kind_to_mode(wc.kind),
        alg    = lambda wc: kind_to_alg(wc.kind),
        n_pcs  = HCLUST_KMEANS_N_PCS,
        seed   = PIPELINE_SEED
    log: os.path.join(LOGROOT, "consensus_cluster_{direction}_{kind}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=CONS_DIRECTION_PATTERN,
        kind=CONS_KIND_PATTERN
    shell:
        r'''
        mkdir -p {params.outdir}
        Rscript {SCRIPTS_DIR}/consensus_ccp_cell_tumour.R \
          --config {input.config} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --feature_list {input.gene_list} \
          --kind {wildcards.kind} \
          --mode {params.mode} \
          --alg {params.alg} \
          --seed {params.seed} \
          --n_pcs {params.n_pcs} \
          --outdir {params.outdir} \
          --cluster_rds {output.cluster_rds} > {log} 2>&1
        '''

# Rule: consensus_cluster_all
# Method role: collector target that requires all consensus clustering outputs for configured directions.
# Flow: expected consensus output paths -> Snakemake dependency aggregation only.
# Stage output: provides a target marker without running a new analysis.
rule consensus_cluster_all:
    """
    Convenience aggregation rule that collects all consensus cluster outputs
    across every direction and kind. HC kinds run on all directions; k-means
    kinds are constrained to Euclidean directions only.
    """
    input:
        expand(os.path.join(CONS_ROOT, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=CONS_DIRECTIONS, kind=CONS_HC_KINDS),
        expand(os.path.join(CONS_ROOT, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=CONS_EUC_DIRECTIONS, kind=CONS_KM_KINDS)


# =============================================================================
# TUMOUR-NEIGHBOURHOOD ANALYSIS PATHS
# =============================================================================
# Stage role: defines tumour-neighbourhood output paths for configured feature-distance representations.

TUMOUR_NH_ROOT     = os.path.join(UNSUP_REL, "tumour_neighbourhoods")
# Absolute counterpart used in shell find commands to ensure correctness on compute
# nodes where the job working directory may differ from PIPE_ROOT.
TUMOUR_NH_ROOT_ABS = os.path.join(UNSUP,     "tumour_neighbourhoods")

# Path helper functions for tumour neighbourhood outputs (single source of truth).
def nh_final_consensus_rds(direction):
    return os.path.join(
        TUMOUR_NH_ROOT, direction, "final_consensus",
        f"Final_consensus_tumour_neighbourhoods_{direction}.rds"
    )

def nh_qc_dir(direction):
    return os.path.join(TUMOUR_NH_ROOT, direction, "qc")

def nh_qc_tsv(direction):
    return os.path.join(nh_qc_dir(direction), "nh_qc_summary.tsv")

def nh_umap_tsv(direction):
    return os.path.join(nh_qc_dir(direction), "nh_umap.tsv")

def nh_umap_pdf(direction):
    return os.path.join(nh_qc_dir(direction), "nh_umap.pdf")


def nh_consensus_dependency(direction):
    sentinel = ".tumour_neighbourhoods_km_done" if direction.endswith("_euc") else ".tumour_neighbourhoods_done"
    return os.path.join(TUMOUR_NH_ROOT, direction, sentinel)


# =============================================================================
# TUMOUR-NEIGHBOURHOOD ANALYSIS
# =============================================================================
# Stage role: ranks patient tumours around each cell-line within each clustering representation.

# Rule: tumour_nh_hc
# Method role: analysis rule that runs tumour-neighbourhood scoring from hierarchical-clustering evidence.
# Flow: representation-specific JOINT HC/k-means and consensus-clustering inputs -> neighbourhood QC, UMAP, and p-consensus tables.
# Feeds within-representation p-consensus and patient-referenced graph construction.
rule tumour_nh_hc:
    """
    Computes per-cell-line tumour neighbourhoods for all configured directions
    using hierarchical-clustering outputs from BOTH clustering formulations —
    the HC/k-means clustering branch (JOINT HCLUST_*/KMEANS_* outputs) and
    ConsensusClusterPlus (JOINT CCP_* outputs) — with profile-configured
    expression inputs. Only JOINT
    cell-line + tumour clustering outputs feed neighbourhood construction.
    """
    input:
        cfg_file    = CFGFILE_ABS,
        cluster_rds = lambda wc: expand(
            os.path.join(CONS_ROOT, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=CONS_HC_KINDS
        ),
        # JOINT HC/k-means clustering outputs (HC family). Declaring them here
        # connects the HC/k-means clustering formulations to
        # tumour-neighbourhood generation so they contribute to p_consensus.
        hclust_kmeans_cluster_rds = lambda wc: expand(
            os.path.join(HCLUST_KMEANS_ROOT_REL, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=HCLUST_JOINT_KINDS
        )
    output:
        done = touch(os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_done"))
    params:
        script = os.path.join(BASE, "scripts", "compute_tumour_neighbourhoods.R")
    log: os.path.join(LOGROOT, "tumour_nh_hc_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=CONS_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.done})
        Rscript {params.script} \
          --config "{input.cfg_file}" \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --unsup-root {UNSUP} \
          --cluster-family hc \
          > {log} 2>&1
        test "$(find "{TUMOUR_NH_ROOT_ABS}/{wildcards.direction}" -maxdepth 2 -type f -name "Top_m_long_*.csv" | wc -l)" -gt 0 || (echo "ERROR: No Top_m_long_*.csv" >&2; exit 1)
        touch {output.done}
        '''


# Rule: tumour_nh_km
# Method role: analysis rule that runs tumour-neighbourhood scoring from k-means clustering evidence.
# Flow: Euclidean representation inputs -> neighbourhood QC, UMAP, and p-consensus tables.
# Analysis role: complements hierarchical-clustering evidence for p-consensus aggregation.
rule tumour_nh_km:
    """
    k-means neighbourhood pass for all configured Euclidean feature-distance
    representations using profile-configured expression inputs. Consumes BOTH
    clustering formulations of the k-means family: the HC/k-means clustering
    branch (JOINT KMEANS_* outputs) and ConsensusClusterPlus (JOINT CCP_KMEANS_*
    outputs).
    """
    input:
        hc_done     = os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_done"),
        cfg_file    = CFGFILE_ABS,
        cluster_rds = lambda wc: expand(
            os.path.join(CONS_ROOT, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=CONS_KM_KINDS
        ),
        # JOINT HC/k-means clustering outputs (k-means family, Euclidean only).
        hclust_kmeans_cluster_rds = lambda wc: expand(
            os.path.join(HCLUST_KMEANS_ROOT_REL, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=KMEANS_JOINT_KINDS
        )
    output:
        done = touch(os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_km_done"))
    params:
        script = os.path.join(BASE, "scripts", "compute_tumour_neighbourhoods.R")
    log: os.path.join(LOGROOT, "tumour_nh_km_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=CONS_EUC_PATTERN
    shell:
        r'''
        Rscript {params.script} \
          --config "{input.cfg_file}" \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --unsup-root {UNSUP} \
          --cluster-family km \
          > {log} 2>&1
        test "$(find "{TUMOUR_NH_ROOT_ABS}/{wildcards.direction}" -maxdepth 2 -type f -name "Top_m_long_*_KM_*.csv" | wc -l)" -gt 0 || (echo "ERROR: No KM-derived Top_m_long_*_KM_*.csv for {wildcards.direction}" >&2; exit 1)
        touch {output.done}
        '''


# =============================================================================
# WITHIN-REPRESENTATION P-CONSENSUS
# =============================================================================
# Stage role: computes tumour-neighbourhood p-consensus values within each representation.

# Rule: tumour_nh_consensus
# Method role: analysis rule that combines tumour-neighbourhood runs within a representation.
# Flow: HC/KM neighbourhood outputs -> final p-consensus RDS for the direction.
# Provides input for cross-representation p-consensus aggregation.
rule tumour_nh_consensus:
    """
    Within-representation p_consensus for every configured direction:
    the tumour-wise recurrence fraction of each cell line-tumour pair across
    the exact configured set of eligible clustering formulations
    (patient_referenced_graph.clustering_methods_by_distance). Eligible
    formulations are the JOINT HC/k-means and ConsensusClusterPlus outputs:
    8 for Euclidean directions, 4 for correlation directions. Euclidean
    directions depend on the KM sentinel so HC and KM neighbourhood files are
    present; correlation directions depend on HC only.
    """
    input:
        done = lambda wc: nh_consensus_dependency(wc.direction),
        cfg = CFGFILE_ABS
    output:
        consensus_rds = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "Final_consensus_tumour_neighbourhoods_{direction}.rds"),
        consensus_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "Final_consensus_tumour_neighbourhoods_{direction}.tsv")
    params:
        script = os.path.join(BASE, "scripts", "tumour_neighbourhood_p_consensus.R"),
        # The configured p-consensus fraction threshold is passed explicitly so
        # the stage's summary statistics and figures report the value the graph
        # stage actually applies.
        threshold = lambda wc: PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD
    log: os.path.join(LOGROOT, "tumour_nh_consensus_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=CONS_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.consensus_rds})
        Rscript {params.script} \
          --config {input.cfg} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --threshold {params.threshold} \
          --out_rds {output.consensus_rds} \
          --out_tsv {output.consensus_tsv} \
          > {log} 2>&1
        '''


# =============================================================================
# PATIENT-REFERENCED GRAPH SCIENTIFIC PARAMETERS (configuration-owned)
# =============================================================================
# Resolved before any graph rule is declared, so both the primary graph branch
# and the metric-specific branch transmit the same configured values.

PATIENT_REFERENCED_GRAPH_CFG = cfg.get("patient_referenced_graph")
if not isinstance(PATIENT_REFERENCED_GRAPH_CFG, dict):
    raise ValueError("Missing required config section: patient_referenced_graph")

def require_patient_referenced_scalar(key, *, lower=None, upper=None, inclusive_upper=True):
    if key not in PATIENT_REFERENCED_GRAPH_CFG:
        raise ValueError(f"Missing required patient_referenced_graph key: {key}")
    value = PATIENT_REFERENCED_GRAPH_CFG[key]
    if not isinstance(value, (int, float)):
        raise ValueError(f"patient_referenced_graph.{key} must be numeric")
    value = float(value)
    if lower is not None and not (value > lower):
        raise ValueError(f"patient_referenced_graph.{key} must be > {lower}")
    if upper is not None:
        if inclusive_upper and not (value <= upper):
            raise ValueError(f"patient_referenced_graph.{key} must be <= {upper}")
        if not inclusive_upper and not (value < upper):
            raise ValueError(f"patient_referenced_graph.{key} must be < {upper}")
    return value

PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD = require_patient_referenced_scalar(
    "p_consensus_threshold", lower=0, upper=1, inclusive_upper=True
)
PATIENT_REFERENCED_SIMILARITY_QUANTILE = require_patient_referenced_scalar(
    "similarity_quantile", lower=0, upper=1, inclusive_upper=False
)
if "similarity_metrics" not in PATIENT_REFERENCED_GRAPH_CFG:
    raise ValueError("Missing required patient_referenced_graph key: similarity_metrics")
PATIENT_REFERENCED_SIMILARITY_METRICS = PATIENT_REFERENCED_GRAPH_CFG["similarity_metrics"]
if not isinstance(PATIENT_REFERENCED_SIMILARITY_METRICS, list) or not PATIENT_REFERENCED_SIMILARITY_METRICS:
    raise ValueError("patient_referenced_graph.similarity_metrics must be a non-empty list")
if "primary_similarity_metric" not in PATIENT_REFERENCED_GRAPH_CFG:
    raise ValueError("Missing required patient_referenced_graph key: primary_similarity_metric")
PATIENT_REFERENCED_PRIMARY_SIMILARITY_METRIC = PATIENT_REFERENCED_GRAPH_CFG["primary_similarity_metric"]
if PATIENT_REFERENCED_PRIMARY_SIMILARITY_METRIC not in PATIENT_REFERENCED_SIMILARITY_METRICS:
    raise ValueError(
        "patient_referenced_graph.primary_similarity_metric must be present in "
        "patient_referenced_graph.similarity_metrics"
    )
if len(set(PATIENT_REFERENCED_SIMILARITY_METRICS)) != len(PATIENT_REFERENCED_SIMILARITY_METRICS):
    raise ValueError("patient_referenced_graph.similarity_metrics must not contain duplicates")
PATIENT_REFERENCED_SIMILARITY_METRIC_PATTERN = "|".join(map(re.escape, PATIENT_REFERENCED_SIMILARITY_METRICS))

# Implementation capability, not an active-universe declaration: the pairwise
# metric-comparison script implements the Pearson-versus-Jaccard contrast and
# names its output columns after it. Configuration still owns which metrics are
# active; this constant only records which configured metric sets the pairwise
# comparison stage can express, and the comparison rule is declared only when
# the configured set matches it exactly.
PATIENT_REFERENCED_METRIC_COMPARISON_IMPLEMENTED_METRICS = ("pearson", "jaccard")
PATIENT_REFERENCED_METRIC_COMPARISON_ENABLED = (
    set(PATIENT_REFERENCED_SIMILARITY_METRICS)
    == set(PATIENT_REFERENCED_METRIC_COMPARISON_IMPLEMENTED_METRICS)
)
PATIENT_REFERENCED_METRIC_A, PATIENT_REFERENCED_METRIC_B = (
    PATIENT_REFERENCED_METRIC_COMPARISON_IMPLEMENTED_METRICS
)

# =============================================================================
# PATIENT-REFERENCED CELL-LINE SIMILARITY GRAPHS
# =============================================================================
# Stage role: converts per-representation tumour-neighbourhood evidence into patient-referenced cell-line graph sidecars.

# Rule: cell_line_similarity_graph
# Method role: analysis rule that converts p-consensus neighbourhood evidence into a cell-line graph per representation.
# Flow: representation p-consensus output -> graph edges and node annotations.
# Feeds resolved-neighbour and multi-representation support-network stages.
rule cell_line_similarity_graph:
    input:
        consensus_rds = lambda wc: nh_final_consensus_rds(wc.direction),
        cfg = CFGFILE_ABS
    output:
        sim_mat_rds = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_matrix_{direction}.rds"),
        sim_pairs_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_pairs_{direction}.tsv"),
        edges_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_graph_edges_{direction}.tsv"),
        nodes_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_graph_node_summary_{direction}.tsv"),
        nodes_annot = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_graph_node_annotations_{direction}.tsv"),
        comm_summary = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_graph_community_summary_{direction}.tsv"),
        comm_table = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_louvain_vs_leiden_community_table_{direction}.tsv"),
        hist_pdf = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_histogram_{direction}.pdf"),
        scatter_pdf = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_DSMZ_p_consensus_cell_scatter_{direction}.pdf"),
        comm_heatmap = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_Louvain_vs_Leiden_heatmap_{direction}.pdf"),
        graph_leiden = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_graph_Leiden_{direction}.pdf"),
        graph_louvain = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_graph_Louvain_{direction}.pdf"),
        graph_minimal = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_graph_minimal_{direction}.pdf"),
        provenance_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "cell_line_similarity_graph_provenance_{direction}.tsv")
    params:
        script = os.path.join(BASE, "scripts", "compute_cell_line_similarity.R"),
        seed = PIPELINE_SEED,
        # The similarity metric of the primary graph branch is a configured
        # scientific choice (patient_referenced_graph.primary_similarity_metric),
        # validated above to be a member of
        # patient_referenced_graph.similarity_metrics. It is transmitted
        # explicitly so the script never selects a metric of its own.
        similarity_metric = PATIENT_REFERENCED_PRIMARY_SIMILARITY_METRIC,
        # Graph edge threshold quantile and p-consensus fraction threshold are
        # configuration-owned and passed explicitly; the script has no defaults.
        similarity_quantile = PATIENT_REFERENCED_SIMILARITY_QUANTILE,
        consensus_threshold = PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD
    log: os.path.join(LOGROOT, "cell_line_similarity_graph_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.sim_mat_rds})
        Rscript {params.script} \
          --config {input.cfg} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --seed {params.seed} \
          --similarity_metric {params.similarity_metric} \
          --similarity_quantile {params.similarity_quantile} \
          --consensus_threshold {params.consensus_threshold} \
          --provenance_tsv {output.provenance_tsv} \
          > {log} 2>&1
        test -s {output.provenance_tsv} || (echo "ERROR: missing {output.provenance_tsv}" >&2; exit 1)
        '''

# =============================================================================
# P-CONSENSUS AGGREGATION
# =============================================================================
# Stage role: compares p-consensus results across representations and prepares graph-resolution inputs.
P_CONS_ALL_DIR = os.path.join(TUMOUR_NH_ROOT, "final_consensus_all")
P_CONS_PLOTS_DIR = os.path.join(P_CONS_ALL_DIR, "plots")
P_CONS_PATIENT_REFERENCED_COHORTS = ["brca", "nbl", "rbl"]
P_CONS_PATIENT_REFERENCED_COHORT_LABELS = {
    "brca": "BRCA",
    "nbl": "NBL",
    "rbl": "RBL",
}
P_CONS_FEATURE_DISTANCE_REPRESENTATION_DIR_ABS = os.path.join(
    BASE,
    "results",
    "unsupervised",
    "p_consensus",
    "feature_distance_representations",
)
P_CONS_FEATURE_DISTANCE_REPRESENTATION_DIR = os.path.relpath(
    P_CONS_FEATURE_DISTANCE_REPRESENTATION_DIR_ABS,
    BASE,
)

P_CONS_FEATURE_DISTANCE_REPRESENTATION_PREFIX = os.path.join(
    P_CONS_FEATURE_DISTANCE_REPRESENTATION_DIR,
    "Fig_p_consensus_feature_distance_representations_dumbbell",
)

P_CONS_DIRECTION_DUMBBELL_PREFIX = P_CONS_FEATURE_DISTANCE_REPRESENTATION_PREFIX

def p_consensus_direction_summary_for_profile(profile):
    cfgp = get_profile_cfg(profile)
    unsup_root = cfgp.get("paths", {}).get("unsup_root")
    if not unsup_root:
        raise KeyError(f"paths.unsup_root missing for profile '{profile}'")
    return os.path.join(
        abspath(unsup_root),
        "tumour_neighbourhoods",
        "final_consensus_all",
        "p_consensus_direction_summary.tsv",
    )

def disease_display_label(cancer_type):
    """Return a visible full disease name from the active profile cancer type."""
    label = str(cancer_type or "").strip()
    lookup = {
        "BRCA": "Breast Cancer",
        "BREAST CANCER": "Breast Cancer",
        "NBL": "Neuroblastoma",
        "NEUROBLASTOMA": "Neuroblastoma",
        "RBL": "Retinoblastoma",
        "RETINOBLASTOMA": "Retinoblastoma",
    }
    return lookup.get(label.upper(), label or str(profile_name))

def slugify_label(label):
    """Convert a visible label into a stable lowercase filename slug."""
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(label).strip().lower())
    return re.sub(r"_+", "_", slug).strip("_")

P_CONS_DISEASE_LABEL = disease_display_label(
    cfg.get("analysis", {}).get("cancer_type", profile_name)
)
P_CONS_DISEASE_SLUG = slugify_label(P_CONS_DISEASE_LABEL)
P_CONS_PER_CELLLINE_FEATURE_DISTANCE_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"per_cellline_feature_distance_{P_CONS_DISEASE_SLUG}_cleveland",
)

_COHORT_UPPER = profile_name.upper()
# Naming-prefix selector for patient-referenced graph outputs.  The
# multicohort_cancer profile is the canonical multicohort result tree, so graph
# sidecars in that tree use a multicohort_cancer_ prefix rather than the
# pan_cancer_ namespace reserved for feature-space/ranking outputs.  Figure
# names already carry Fig_MULTICOHORT_CANCER and therefore omit the repeated
# multicohort_cancer_ token from the figure stem.
_PR_PREFIX = "multicohort_cancer_" if IS_MULTICOHORT_PROFILE else "patient_referenced_"
_PR_FIGURE_PREFIX = "" if IS_MULTICOHORT_PROFILE else _PR_PREFIX
P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_{_PR_FIGURE_PREFIX}resolved_cell_line_neighbourhood_graph",
)
P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_{_PR_FIGURE_PREFIX}support_threshold_consensus_cell_line_similarity_network",
)
P_CONS_GRAPH_NODE_STATS_DIR = os.path.join(
    P_CONS_ALL_DIR,
    "analysis",
    "resolved_graph_node_statistics",
)
P_CONS_SHORTNAMES_TSV = os.path.join(P_CONS_GRAPH_NODE_STATS_DIR, f"{_PR_PREFIX}cell_line_display_names.tsv")
P_CONS_RESOLVED_EDGES_TSV = os.path.join(P_CONS_GRAPH_NODE_STATS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_edges.tsv")
P_CONS_RESOLVED_NODE_STATS_TSV = os.path.join(P_CONS_GRAPH_NODE_STATS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_node_stats.tsv")
P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}support_threshold_consensus_cell_line_similarity_edges.tsv")
P_CONS_MULTI_REP_EDGE_SUPPORT_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_edge_support_all_pairs.tsv",
)
P_CONS_MULTI_REP_MAJORITY_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_multi_representation_majority_threshold_consensus_network",
)
P_CONS_MULTI_REP_MAJORITY_EDGES_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_majority_threshold_consensus_network_edges.tsv",
)
P_CONS_MULTI_REP_MAJORITY_NODE_STATS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_majority_threshold_consensus_network_node_stats.tsv",
)
P_CONS_MULTI_REP_MAJORITY_COMPONENTS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_majority_threshold_consensus_network_components.tsv",
)
P_CONS_MULTI_REP_UNION_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_multi_representation_union_supported_edges_network",
)
P_CONS_MULTI_REP_UNION_FULL_LABEL_PREFIX = (
    P_CONS_MULTI_REP_UNION_PREFIX + "_full_node_labels"
)
P_CONS_MULTI_REP_UNION_EDGES_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_edges.tsv",
)
P_CONS_MULTI_REP_UNION_NODE_STATS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_node_stats.tsv",
)
P_CONS_MULTI_REP_UNION_COMPONENTS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_components.tsv",
)
P_CONS_MULTI_REP_UNION_SELECTED_NODE_LABELS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_selected_node_labels.tsv",
)
P_CONS_MULTI_REP_UNION_FULL_NODE_LABELS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_full_node_labels.tsv",
)
P_CONS_MULTI_REP_UNION_FULL_NODE_STATS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_full_node_labels_node_stats.tsv",
)
P_CONS_MULTI_REP_UNION_FULL_COMPONENTS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{profile_name}_multi_representation_union_supported_edges_network_full_node_labels_components.tsv",
)
POST_RESOLUTION_COHORTS = ["NBL", "BRCA", "RBL"]
POST_RESOLUTION_TABLE_DIR = os.path.join("results", "tables")
POST_RESOLUTION_STRATIFICATION_TABLES = expand(
    os.path.join(
        POST_RESOLUTION_TABLE_DIR,
        "{cohort}_post_resolution_edge_support_stratification.tsv",
    ),
    cohort=POST_RESOLUTION_COHORTS,
)
POST_RESOLUTION_COMBINED_STRATIFICATION_TSV = os.path.join(
    POST_RESOLUTION_TABLE_DIR,
    "combined_post_resolution_edge_support_stratification.tsv",
)

def _post_resolution_plots_dir(cohort):
    cohort_slug = str(cohort).lower()
    return os.path.join(
        "results",
        "unsupervised",
        cohort_slug,
        "tumour_neighbourhoods",
        "final_consensus_all",
        "plots",
    )


def _post_resolution_union_edges(wildcards):
    cohort_slug = str(wildcards.cohort).lower()
    return os.path.join(
        _post_resolution_plots_dir(wildcards.cohort),
        f"{cohort_slug}_multi_representation_union_supported_edges_network_edges.tsv",
    )


def _post_resolution_consensus_edges(wildcards):
    cohort_slug = str(wildcards.cohort).lower()
    return os.path.join(
        _post_resolution_plots_dir(wildcards.cohort),
        f"{cohort_slug}_multi_representation_majority_threshold_consensus_network_edges.tsv",
    )


def _post_resolution_resolved_edges(wildcards):
    return os.path.join(
        _post_resolution_plots_dir(wildcards.cohort),
        "patient_referenced_resolved_cell_line_neighbourhood_graph_edges.tsv",
    )


def _post_resolution_nodes(wildcards):
    cohort_slug = str(wildcards.cohort).lower()
    return os.path.join(
        _post_resolution_plots_dir(wildcards.cohort),
        f"{cohort_slug}_multi_representation_union_supported_edges_network_node_stats.tsv",
    )

P_CONS_ANCHOR_AUDIT_TSV            = os.path.join(P_CONS_GRAPH_NODE_STATS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv")
P_CONS_DSMZ_META_CSV = cfgget_path_abs(
    cfgget_path_rel("data/dsmz/DSMZ_metadata.csv", "paths", "meta_tsv"),
    "paths", "dsmz_meta_csv",
)
P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_component_panels"
P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_component_panels"
P_CONS_RESOLVED_INTERACTIVE_HTML = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_interactive.html"
P_CONS_SUPPORT_INTERACTIVE_HTML = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_interactive.html"
P_CONS_RESOLVED_COMPONENT_SUMMARY_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_component_summary.tsv",
)
P_CONS_SUPPORT_COMPONENT_SUMMARY_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{_PR_PREFIX}support_threshold_consensus_cell_line_similarity_component_summary.tsv",
)
P_CONS_RESOLVED_NODE_LABELS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_full_node_labels.tsv",
)
P_CONS_SUPPORT_NODE_LABELS_TSV = os.path.join(
    P_CONS_PLOTS_DIR,
    f"{_PR_PREFIX}support_threshold_consensus_cell_line_similarity_full_node_labels.tsv",
)
PAN_CANCER_GRAPH_INSPECTION_TARGETS = (
    [
        P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".pdf",
        P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".png",
        P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".svg",
        P_CONS_RESOLVED_INTERACTIVE_HTML,
        P_CONS_RESOLVED_COMPONENT_SUMMARY_TSV,
        P_CONS_RESOLVED_NODE_LABELS_TSV,
    ]
    if IS_MULTICOHORT_PROFILE else []
)
SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT = config.get("support_threshold_consensus_plot_layout", {})
SUPPORT_THRESHOLD_COMPONENT_LABEL_GAP = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_label_gap",
        2.0 if _COHORT_UPPER in ("BRCA", "MULTICOHORT_CANCER") else 0.40,
    )
)
SUPPORT_THRESHOLD_NBL_COMPONENT_GAP_FRACTION = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("nbl_component_gap_fraction", 0.57)
)
SUPPORT_THRESHOLD_NBL_COMPONENT_MIN_GAP = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("nbl_component_min_gap", 9.5)
)
SUPPORT_THRESHOLD_LEGEND_ANCHOR_X = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "legend_anchor_x",
        0.81 if _COHORT_UPPER == "NBL" else (0.72 if _COHORT_UPPER == "RBL" else 1.01),
    )
)
SUPPORT_THRESHOLD_LEGEND_ANCHOR_Y = SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
    "legend_anchor_y",
    1.30 if _COHORT_UPPER == "BRCA" else None,
)
SUPPORT_THRESHOLD_LEGEND_ANCHOR_Y_ARG = (
    "" if SUPPORT_THRESHOLD_LEGEND_ANCHOR_Y is None
    else f"--legend-anchor-y {float(SUPPORT_THRESHOLD_LEGEND_ANCHOR_Y)}"
)
SUPPORT_THRESHOLD_LEGEND_ALIGN_MODE = SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
    "legend_align_mode",
    "default" if _COHORT_UPPER in ("BRCA", "RBL") else "component-band",
)
SUPPORT_THRESHOLD_LEGEND_GAP_FRAC = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("legend_gap_frac", 0.025)
)
SUPPORT_THRESHOLD_ISOLATE_LABEL_Y_FRAC = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_label_y_frac",
        0.02 if _COHORT_UPPER == "NBL" else 0.010,
    )
)
SUPPORT_THRESHOLD_ISOLATE_BOX_HEIGHT_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_box_height_factor",
        2.00,
    )
)
SUPPORT_THRESHOLD_ISOLATE_GAP = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_gap",
        8.0 if _COHORT_UPPER in ("BRCA", "MULTICOHORT_CANCER") else 4.0,
    )
)
SUPPORT_THRESHOLD_ISOLATE_SPACING = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_spacing",
        18.0 if _COHORT_UPPER in ("BRCA", "MULTICOHORT_CANCER") else 7.0,
    )
)
SUPPORT_THRESHOLD_DENSE_COMPONENT_MIN_SIZE = int(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("dense_component_min_size", 8)
)
SUPPORT_THRESHOLD_DENSE_COMPONENT_MIN_DENSITY = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("dense_component_min_density", 0.25)
)
SUPPORT_THRESHOLD_DENSE_COMPONENT_EXPAND_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("dense_component_expand_factor", 1.85)
)
SUPPORT_THRESHOLD_COMPONENT_GRID_NCOLS = int(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_grid_ncols",
        4 if _COHORT_UPPER == "MULTICOHORT_CANCER" else 3,
    )
)
SUPPORT_THRESHOLD_COMPONENT_GRID_NROWS = int(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_grid_nrows",
        3 if _COHORT_UPPER == "MULTICOHORT_CANCER" else 2,
    )
)
SUPPORT_THRESHOLD_COMPONENT_CELL_WIDTH = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_cell_width",
        84.0 if _COHORT_UPPER == "MULTICOHORT_CANCER" else (62.0 if _COHORT_UPPER == "BRCA" else 54.0),
    )
)
SUPPORT_THRESHOLD_COMPONENT_CELL_HEIGHT = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_cell_height",
        42.0 if _COHORT_UPPER == "MULTICOHORT_CANCER" else (54.0 if _COHORT_UPPER == "BRCA" else 27.0),
    )
)
SUPPORT_THRESHOLD_COMPONENT_ROW_GAP = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_row_gap",
        14.0 if _COHORT_UPPER in ("BRCA", "MULTICOHORT_CANCER") else 5.0,
    )
)
SUPPORT_THRESHOLD_COMPONENT_COL_GAP = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "component_col_gap",
        14.0 if _COHORT_UPPER == "MULTICOHORT_CANCER" else (12.0 if _COHORT_UPPER == "BRCA" else 8.0),
    )
)
SUPPORT_THRESHOLD_LEGEND_WIDTH = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "legend_width",
        64.0 if _COHORT_UPPER == "MULTICOHORT_CANCER" else (62.0 if _COHORT_UPPER == "BRCA" else 44.0),
    )
)
SUPPORT_THRESHOLD_ISOLATE_REGION_HEIGHT = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_region_height",
        20.0 if _COHORT_UPPER == "MULTICOHORT_CANCER" else (18.0 if _COHORT_UPPER == "BRCA" else 8.8),
    )
)
SUPPORT_THRESHOLD_ISOLATE_PANEL_PADDING_X = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("isolate_panel_padding_x", 12.0)
)
SUPPORT_THRESHOLD_ISOLATE_PANEL_PADDING_Y = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("isolate_panel_padding_y", 1.4)
)
SUPPORT_THRESHOLD_ISOLATE_LABEL_BAND_FRAC = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("isolate_label_band_frac", 0.34)
)
SUPPORT_THRESHOLD_ISOLATE_MAX_PER_ROW = int(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_max_per_row",
        7 if _COHORT_UPPER == "MULTICOHORT_CANCER" else 6,
    )
)
SUPPORT_THRESHOLD_ISOLATE_PANEL_WIDTH_FRAC = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "isolate_panel_width_frac",
        0.80 if _COHORT_UPPER == "MULTICOHORT_CANCER" else 0.62,
    )
)
SUPPORT_THRESHOLD_FOOTNOTE_REGION_HEIGHT = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("footnote_region_height", 7.0)
)
SUPPORT_THRESHOLD_COMPONENT_CELL_FILL_X = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("component_cell_fill_x", 0.78)
)
SUPPORT_THRESHOLD_COMPONENT_CELL_FILL_Y = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("component_cell_fill_y", 0.72)
)
SUPPORT_THRESHOLD_DENSE_COMPONENT_CELL_FILL_X = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("dense_component_cell_fill_x", 0.92)
)
SUPPORT_THRESHOLD_DENSE_COMPONENT_CELL_FILL_Y = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("dense_component_cell_fill_y", 0.84)
)
SUPPORT_THRESHOLD_C2_COMPONENT_CELL_FILL_X = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("c2_component_cell_fill_x", 0.88)
)
SUPPORT_THRESHOLD_C2_COMPONENT_CELL_FILL_Y = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get("c2_component_cell_fill_y", 0.80)
)
SUPPORT_THRESHOLD_C1_X_EXPAND_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "c1_x_expand_factor",
        4.00 if _COHORT_UPPER == "BRCA" else 2.60,
    )
)
SUPPORT_THRESHOLD_C1_Y_EXPAND_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "c1_y_expand_factor",
        4.00 if _COHORT_UPPER == "BRCA" else 1.30,
    )
)
SUPPORT_THRESHOLD_C2_X_EXPAND_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "c2_x_expand_factor",
        2.00 if _COHORT_UPPER == "BRCA" else 1.60,
    )
)
SUPPORT_THRESHOLD_C2_Y_EXPAND_FACTOR = float(
    SUPPORT_THRESHOLD_CONSENSUS_PLOT_LAYOUT.get(
        "c2_y_expand_factor",
        2.00 if _COHORT_UPPER == "BRCA" else 1.15,
    )
)
P_CONS_COMMUNITY_DETECTION_DIR     = os.path.join(P_CONS_ALL_DIR, "community_detection")
P_CONS_LEIDEN_COMMUNITIES_TSV      = os.path.join(P_CONS_COMMUNITY_DETECTION_DIR, "multicohort_cancer_communities.tsv")
P_CONS_LEIDEN_SUMMARY_TSV          = os.path.join(P_CONS_COMMUNITY_DETECTION_DIR, "multicohort_cancer_community_summary.tsv")
P_CONS_LEIDEN_MODULARITY_TSV       = os.path.join(P_CONS_COMMUNITY_DETECTION_DIR, "multicohort_cancer_modularity.tsv")
P_CONS_LEIDEN_LAYOUT_TSV           = os.path.join(P_CONS_COMMUNITY_DETECTION_DIR, "multicohort_cancer_layout.tsv")
P_CONS_LEIDEN_FIG_PREFIX           = os.path.join(P_CONS_PLOTS_DIR, "Fig_MULTICOHORT_CANCER_cancer_type_community")
COMMUNITY_STABILITY_DIR = os.path.join(P_CONS_ALL_DIR, "community_stability")
COMMUNITY_STABILITY_ASSIGNMENTS = os.path.join(COMMUNITY_STABILITY_DIR, "community_assignments_long.tsv")
COMMUNITY_STABILITY_COMPONENT_SUMMARY = os.path.join(COMMUNITY_STABILITY_DIR, "community_stability_component_summary.tsv")
COMMUNITY_STABILITY_GLOBAL_SUMMARY = os.path.join(COMMUNITY_STABILITY_DIR, "community_stability_global_summary.tsv")
COMMUNITY_STABILITY_VALIDATION = os.path.join(COMMUNITY_STABILITY_DIR, "community_stability_validation_report.txt")
COMMUNITY_STABILITY_OUTPUTS = [
    COMMUNITY_STABILITY_ASSIGNMENTS,
    os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_matrix_Leiden.tsv"),
    os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_matrix_Louvain.tsv"),
    os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_n_valid_Leiden.tsv"),
    os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_n_valid_Louvain.tsv"),
    COMMUNITY_STABILITY_COMPONENT_SUMMARY,
    COMMUNITY_STABILITY_GLOBAL_SUMMARY,
    os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_coassignment_heatmap_Leiden.pdf"),
    os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_coassignment_heatmap_Louvain.pdf"),
    os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_stability_within_between_Leiden.pdf"),
    os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_stability_within_between_Louvain.pdf"),
    COMMUNITY_STABILITY_VALIDATION,
]
COMMUNITY_STABILITY_CFG = cfg.get("community_stability", {})
COMMUNITY_STABILITY_ENABLED_PROFILES = {
    str(p) for p in COMMUNITY_STABILITY_CFG.get("enabled_profiles", ["brca", "nbl", "rbl"])
}
COMMUNITY_STABILITY_ENABLED = bool(COMMUNITY_STABILITY_CFG.get("enabled", True)) and (
    profile_name in COMMUNITY_STABILITY_ENABLED_PROFILES
)
COMMUNITY_STABILITY_TARGETS = COMMUNITY_STABILITY_OUTPUTS if COMMUNITY_STABILITY_ENABLED else []
STUDY_OUTPUT_DIR_REL = os.path.join(UNSUP_REL, config.get("study_design", {}).get("outputs_dirname", "study_design"))
VALIDATION_OUTPUT_DIR_REL = os.path.join(UNSUP_REL, "validation")

# =============================================================================
# PARALLEL PEARSON–JACCARD GRAPH COMPARISON
# =============================================================================
# Stage role: builds metric-specific patient-referenced graphs in parallel from a shared thresholded pairwise active-tumour-union preprocessing step.

PATIENT_REFERENCED_SIMILARITY_ROOT = os.path.join(UNSUP_REL, "patient_referenced_similarity_metrics")
PATIENT_REFERENCED_SIMILARITY_ROOT_ABS = os.path.join(UNSUP, "patient_referenced_similarity_metrics")
PATIENT_REFERENCED_SIMILARITY_ALL_DIR = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "final_consensus_all")
PATIENT_REFERENCED_SIMILARITY_ALL_DIR_ABS = os.path.join(UNSUP, "patient_referenced_similarity_metrics", "final_consensus_all")
PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "comparisons")
PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR_ABS = os.path.join(UNSUP, "patient_referenced_similarity_metrics", "comparisons")

def patient_referenced_metric_final_consensus_dir(metric, direction):
    return os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, metric, direction, "final_consensus")

def patient_referenced_metric_final_consensus_dir_abs(metric, direction):
    return os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT_ABS, metric, direction, "final_consensus")

def patient_referenced_metric_final_dir(metric):
    return os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, metric)

def patient_referenced_metric_final_dir_abs(metric):
    return os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR_ABS, metric)

rule cell_line_similarity_graph_metric:
    """
    Builds a metric-specific patient-referenced graph branch for each active
    feature-distance representation. Pearson correlates threshold-restricted
    p-consensus fractions over the pairwise active-tumour union. Jaccard
    measures binary overlap over the same pairwise active-tumour union.
    """
    input:
        consensus_rds = lambda wc: nh_final_consensus_rds(wc.direction),
        cfg = CFGFILE_ABS
    output:
        sim_pairs_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "{metric}", "{direction}", "final_consensus", "cell_line_similarity_pairs_{direction}.tsv"),
        edges_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "{metric}", "{direction}", "final_consensus", "cell_line_similarity_graph_edges_{direction}.tsv"),
        nodes_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "{metric}", "{direction}", "final_consensus", "cell_line_similarity_graph_node_summary_{direction}.tsv"),
        provenance_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "{metric}", "{direction}", "final_consensus", "cell_line_similarity_graph_provenance_{direction}.tsv"),
        selected_tumours_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, "{metric}", "{direction}", "final_consensus", "cell_line_similarity_selected_tumours_{direction}.tsv")
    params:
        script = os.path.join(BASE, "scripts", "compute_cell_line_similarity_metric.R"),
        out_base = lambda wc: patient_referenced_metric_final_consensus_dir_abs(wc.metric, wc.direction),
        similarity_metric = lambda wc: wc.metric,
        consensus_threshold = lambda wc: PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD,
        similarity_quantile = lambda wc: PATIENT_REFERENCED_SIMILARITY_QUANTILE
    log: os.path.join(LOGROOT, "cell_line_similarity_graph_metric_{metric}_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=FEATURE_DISTANCE_DIRECTION_PATTERN,
        metric=PATIENT_REFERENCED_SIMILARITY_METRIC_PATTERN
    shell:
        r'''
        mkdir -p {params.out_base}
        Rscript {params.script} \
          --config {input.cfg} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --out_base {params.out_base} \
          --similarity_metric {params.similarity_metric} \
          --consensus_threshold {params.consensus_threshold} \
          --similarity_quantile {params.similarity_quantile} \
          > {log} 2>&1
        '''

rule resolve_dsmz_graph_neighbours_metric:
    """
    Re-runs the existing resolved-neighbour procedure on the metric-specific
    representation graphs so Pearson and Jaccard remain directly comparable.
    """
    input:
        cfg_file = CFGFILE_ABS,
        winners_tsv = os.path.join(P_CONS_ALL_DIR, "p_consensus_winners_by_frac_ge_thr.tsv"),
        dir_summary_tsv = os.path.join(P_CONS_ALL_DIR, "p_consensus_direction_summary.tsv"),
        graph_edges = lambda wc: expand(
            os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, wc.metric, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        resolved_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, "{metric}", "resolved_dsmz_neighbours.tsv"),
        validation_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, "{metric}", "resolved_graph_input_validation.tsv")
    params:
        script = os.path.join(BASE, "scripts", "resolve_dsmz_graph_neighbours.R"),
        config = CFGFILE_ABS,
        graph_root = lambda wc: os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT_ABS, wc.metric),
        # Resolution-eligible representations are enumerated by orchestration
        # from the effective configuration and passed explicitly, so the script
        # never reconstructs the representation universe.
        directions = ",".join(CONS_DIRECTIONS)
    log: os.path.join(LOGROOT, "resolve_dsmz_graph_neighbours_metric_{metric}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        metric=PATIENT_REFERENCED_SIMILARITY_METRIC_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.resolved_tsv})
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --winners_tsv {input.winners_tsv} \
          --direction_summary_tsv {input.dir_summary_tsv} \
          --graph_root {params.graph_root} \
          --directions {params.directions} \
          --output_tsv {output.resolved_tsv} \
          --validation_tsv {output.validation_tsv} \
          > {log} 2>&1
        test -s {output.resolved_tsv} || (echo "ERROR: missing {output.resolved_tsv}" >&2; exit 1)
        test -s {output.validation_tsv} || (echo "ERROR: missing {output.validation_tsv}" >&2; exit 1)
        '''

rule build_multi_representation_majority_threshold_consensus_network_metric:
    input:
        node_universe = P_CONS_SHORTNAMES_TSV,
        graph_edges = lambda wc: expand(
            os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, wc.metric, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        edges = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, "{metric}", "multi_representation_majority_threshold_edges.tsv"),
        edge_support = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, "{metric}", "multi_representation_majority_threshold_edge_support.tsv")
    params:
        script = os.path.join(BASE, "scripts", "build_consensus_from_direction_edgefiles.py"),
        tumour_nh_dir = lambda wc: os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT_ABS, wc.metric),
        directions = ",".join(CONS_DIRECTIONS),
        min_support = SIMILARITY_CONSENSUS_MIN_SUPPORT,
        cohort = profile_name
    log: os.path.join(LOGROOT, "build_multi_representation_majority_threshold_consensus_network_metric_{metric}.log")
    conda: CONDA_ENV_PY
    wildcard_constraints:
        metric=PATIENT_REFERENCED_SIMILARITY_METRIC_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.edges})
        python {params.script} \
          --tumour_nh_dir {params.tumour_nh_dir} \
          --out_edges {output.edges} \
          --out_support {output.edge_support} \
          --mode majority_threshold \
          --cohort {params.cohort} \
          --min_support {params.min_support} \
          --directions {params.directions} \
          --name_map {input.node_universe} \
          --node_universe {input.node_universe} \
          --require_short \
          > {log} 2>&1
        '''

rule build_multi_representation_union_supported_edges_network_metric:
    input:
        node_universe = P_CONS_SHORTNAMES_TSV,
        graph_edges = lambda wc: expand(
            os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, wc.metric, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        edges = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, "{metric}", "multi_representation_union_edges.tsv")
    params:
        script = os.path.join(BASE, "scripts", "build_consensus_from_direction_edgefiles.py"),
        tumour_nh_dir = lambda wc: os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT_ABS, wc.metric),
        directions = ",".join(CONS_DIRECTIONS),
        cohort = profile_name
    log: os.path.join(LOGROOT, "build_multi_representation_union_supported_edges_network_metric_{metric}.log")
    conda: CONDA_ENV_PY
    wildcard_constraints:
        metric=PATIENT_REFERENCED_SIMILARITY_METRIC_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.edges})
        python {params.script} \
          --tumour_nh_dir {params.tumour_nh_dir} \
          --out_edges {output.edges} \
          --mode union_supported_edges \
          --cohort {params.cohort} \
          --min_support 1 \
          --directions {params.directions} \
          --name_map {input.node_universe} \
          --node_universe {input.node_universe} \
          --require_short \
          > {log} 2>&1
        '''

rule audit_clustering_method_consensus_resolution:
    input:
        cfg = CFGFILE_ABS,
        consensus_tsvs = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "Final_consensus_tumour_neighbourhoods_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        audit_tsv = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "clustering_method_consensus_resolution.tsv")
    params:
        script = os.path.join(BASE, "scripts", "audit_clustering_method_consensus_resolution.py"),
        directions = ",".join(CONS_DIRECTIONS)
    log: os.path.join(LOGROOT, "audit_clustering_method_consensus_resolution.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR_ABS}
        python {params.script} \
          --config {input.cfg} \
          --profile "{profile_name}" \
          --directions {params.directions} \
          --out_tsv {output.audit_tsv} \
          > {log} 2>&1
        '''

rule compare_patient_referenced_similarity_metrics:
    input:
        pearson_pairs = expand(
            os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, PATIENT_REFERENCED_METRIC_A, "{direction}", "final_consensus",
                         "cell_line_similarity_pairs_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        ),
        jaccard_pairs = expand(
            os.path.join(PATIENT_REFERENCED_SIMILARITY_ROOT, PATIENT_REFERENCED_METRIC_B, "{direction}", "final_consensus",
                         "cell_line_similarity_pairs_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        ),
        pearson_resolved = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_A, "resolved_dsmz_neighbours.tsv"),
        jaccard_resolved = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_B, "resolved_dsmz_neighbours.tsv"),
        pearson_majority = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_A, "multi_representation_majority_threshold_edges.tsv"),
        jaccard_majority = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_B, "multi_representation_majority_threshold_edges.tsv"),
        pearson_union = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_A, "multi_representation_union_edges.tsv"),
        jaccard_union = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, PATIENT_REFERENCED_METRIC_B, "multi_representation_union_edges.tsv"),
        consensus_resolution = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "clustering_method_consensus_resolution.tsv")
    output:
        representation_graphs = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_representation_graphs.tsv"),
        pairwise_similarity = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_pairwise_similarity.tsv"),
        edge_agreement = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_edge_agreement.tsv"),
        resolved_graph = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_graph_comparison.tsv"),
        resolved_neighbours = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_neighbours.tsv"),
        provenance = os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "similarity_graph_provenance.tsv")
    params:
        script = os.path.join(BASE, "scripts", "compare_patient_referenced_similarity_metrics.py"),
        metric_root = PATIENT_REFERENCED_SIMILARITY_ROOT_ABS,
        directions = ",".join(CONS_DIRECTIONS),
        # Passed so the script fails if it is ever handed a configured metric
        # set it does not implement, rather than silently comparing the wrong
        # products.
        metrics = ",".join(PATIENT_REFERENCED_SIMILARITY_METRICS)
    log: os.path.join(LOGROOT, "compare_patient_referenced_similarity_metrics.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR_ABS}
        python {params.script} \
          --metric_root {params.metric_root} \
          --directions {params.directions} \
          --metrics {params.metrics} \
          --resolved_pearson {input.pearson_resolved} \
          --resolved_jaccard {input.jaccard_resolved} \
          --majority_pearson {input.pearson_majority} \
          --majority_jaccard {input.jaccard_majority} \
          --union_pearson {input.pearson_union} \
          --union_jaccard {input.jaccard_union} \
          --out_representation {output.representation_graphs} \
          --out_pairwise {output.pairwise_similarity} \
          --out_edge_agreement {output.edge_agreement} \
          --out_resolved_graph {output.resolved_graph} \
          --out_resolved_neighbours {output.resolved_neighbours} \
          --out_provenance {output.provenance} \
          > {log} 2>&1
        '''

rule patient_referenced_similarity_metric_comparison:
    input:
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_representation_graphs.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_pairwise_similarity.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_edge_agreement.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_graph_comparison.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_neighbours.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "similarity_graph_provenance.tsv"),
        os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "clustering_method_consensus_resolution.tsv")

# Rule: summarize_p_consensus_all
# Method role: aggregates p-consensus outputs across configured feature-distance representations.
# Flow: per-direction p-consensus RDS files -> ranking, winner, weight, and cross-direction metric tables.
# Analysis role: selects best-supported representations for resolved-neighbour graph construction and validation.
rule summarize_p_consensus_all:
    input:
        consensus_rds = [nh_final_consensus_rds(d) for d in P_CONS_SUMMARY_DIRECTIONS],
        cfg = CFGFILE_ABS
    output:
        dir_summary = os.path.join(P_CONS_ALL_DIR, "p_consensus_direction_summary.tsv"),
        cell_dir_summary = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.tsv"),
        winner_map = os.path.join(P_CONS_ALL_DIR, "p_consensus_winner_map_per_cell_line.tsv"),
        winners_max_p = os.path.join(P_CONS_ALL_DIR, "p_consensus_winners_by_max_p.tsv"),
        winners_frac = os.path.join(P_CONS_ALL_DIR, "p_consensus_winners_by_frac_ge_thr.tsv"),
        ranked_best = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_ranked.tsv"),
        long_tbl = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.long.tsv"),
        composite_weights = os.path.join(P_CONS_ALL_DIR, "p_consensus_composite_weights.tsv"),
        winning_direction = os.path.join(P_CONS_ALL_DIR, "winning_direction.txt"),
        top_fraction = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_top_fraction.tsv"),
        top_score = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_top_score.tsv"),
        fig_pdf = os.path.join(P_CONS_ALL_DIR, "Fig_p_consensus_direction_comparison.pdf"),
        fig_scree = os.path.join(P_CONS_ALL_DIR, "Fig_p_consensus_composite_PCA_scree.pdf"),
        fig_weights = os.path.join(P_CONS_ALL_DIR, "Fig_p_consensus_composite_PCA_weights.pdf")
    params:
        script = os.path.join(BASE, "scripts", "summarize_p_consensus_all.R"),
        config = CFGFILE_ABS,
        threshold = lambda wc: PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD,
        # The feature-distance representations summarised here are enumerated
        # by orchestration from configuration and passed explicitly, so the
        # script never has to infer the representation universe itself.
        directions = ",".join(P_CONS_SUMMARY_DIRECTIONS)
    log: os.path.join(LOGROOT, "summarize_p_consensus_all.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p {P_CONS_ALL_DIR}
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --threshold {params.threshold} \
          --directions {params.directions} \
          > {log} 2>&1
        '''


# Rule: plot_per_cellline_feature_distance_cleveland
# Method role: final plotting rule for per-cell-line representation comparison.
# Flow: p-consensus fraction tables -> Cleveland-style PDF/PNG plot.
# Purpose: reporting visualisation only; it does not alter analysis targets.
rule plot_per_cellline_feature_distance_cleveland:
    input:
        cell_dir_summary = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.tsv"),
        winners_frac = os.path.join(P_CONS_ALL_DIR, "p_consensus_winners_by_frac_ge_thr.tsv"),
        cfg = CFGFILE_ABS
    output:
        pdf = P_CONS_PER_CELLLINE_FEATURE_DISTANCE_PREFIX + ".pdf",
        png = P_CONS_PER_CELLLINE_FEATURE_DISTANCE_PREFIX + ".png"
    params:
        script = os.path.join(BASE, "scripts", "plot_per_cellline_feature_distance_cleveland.R"),
        disease_label = P_CONS_DISEASE_LABEL,
        threshold = lambda wc: PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD
    log: os.path.join(LOGROOT, "plot_per_cellline_feature_distance_cleveland.log")
    conda: CONDA_ENV_R_BASE
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR} $(dirname {log})
        Rscript {params.script} \
          --summary-tsv {input.cell_dir_summary} \
          --winners-tsv {input.winners_frac} \
          --out-pdf {output.pdf} \
          --out-png {output.png} \
          --disease-label "{params.disease_label}" \
          --threshold {params.threshold} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.png} || (echo "ERROR: missing {output.png}" >&2; exit 1)
        '''


# Rule: plot_p_consensus_direction_comparison_dumbbell
# Method role: supplementary plotting rule visualising p-consensus direction metrics across cohorts.
# Flow: cohort direction-metric tables -> dumbbell figure and caption text.
# Purpose: cross-cohort reporting visualisation, not a modelling step.
rule plot_p_consensus_direction_comparison_dumbbell:
    input:
        brca = p_consensus_direction_summary_for_profile("brca"),
        nbl = p_consensus_direction_summary_for_profile("nbl"),
        rbl = p_consensus_direction_summary_for_profile("rbl")
    output:
        pdf = P_CONS_DIRECTION_DUMBBELL_PREFIX + ".pdf",
        png = P_CONS_DIRECTION_DUMBBELL_PREFIX + ".png",
        svg = P_CONS_DIRECTION_DUMBBELL_PREFIX + ".svg",
        caption = P_CONS_DIRECTION_DUMBBELL_PREFIX + "_caption.txt"
    params:
        script = os.path.join(BASE, "scripts", "plot_p_consensus_direction_dumbbell.R"),
        out_prefix = P_CONS_DIRECTION_DUMBBELL_PREFIX,
        threshold = lambda wc: PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD
    log: os.path.join(LOGROOT, "plot_p_consensus_direction_comparison_dumbbell.log")
    conda: CONDA_ENV_R_BASE
    shell:
        r'''
        mkdir -p $(dirname {output.pdf}) $(dirname {log})
        Rscript {params.script} \
          --brca {input.brca} \
          --nbl {input.nbl} \
          --rbl {input.rbl} \
          --out_prefix {params.out_prefix} \
          --caption {output.caption} \
          --threshold {params.threshold} \
          > {log} 2>&1
        '''

# =============================================================================
# GRAPH-BASED NEIGHBOUR RESOLUTION
# =============================================================================
# Stage role: resolves stable cell-line neighbours from patient-referenced graph evidence and best-direction evidence.

# Rule: resolve_dsmz_graph_neighbours
# Method role: graph-resolution rule that resolves stable cell-line neighbours using patient-referenced graph evidence.
# Flow: p-consensus fraction tables and per-direction graph edges -> resolved-neighbour table and audit sidecars.
# Required by final graph plotting, validation, DESeq2 grouping, and support-network construction.
rule resolve_dsmz_graph_neighbours:
    """
    Resolves per-cell-line DSMZ neighbours by intersecting the global best-direction
    graph with each cell line's winner-direction graph. Consumes the cross-direction
    consensus summary outputs and all per-direction DSMZ similarity graph edge files
    produced by cell_line_similarity_graph. The intersection strategy retains only
    neighbours present in both the global best direction and the cell-line-specific
    winner direction, providing a robust estimate of stable similarity relationships.
    Uses CONS_DIRECTIONS (not FEATURE_DISTANCE_DIRECTIONS) so only configured consensus directions
    are included as explicit DAG dependencies.
    """
    input:
        cfg_file        = CFGFILE_ABS,
        winners_tsv     = os.path.join(P_CONS_ALL_DIR, "p_consensus_winners_by_frac_ge_thr.tsv"),
        dir_summary_tsv = os.path.join(P_CONS_ALL_DIR, "p_consensus_direction_summary.tsv"),
        graph_edges     = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv"),
        validation_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_graph_input_validation.tsv")
    params:
        script = os.path.join(BASE, "scripts", "resolve_dsmz_graph_neighbours.R"),
        config = CFGFILE_ABS,
        # Resolution-eligible representations are enumerated by orchestration
        # from the effective configuration and passed explicitly. The script
        # validates the observed representation graph set against this exact
        # expected set and never reconstructs it from config sections,
        # directory contents, or built-in lists.
        directions = ",".join(CONS_DIRECTIONS)
    log: os.path.join(LOGROOT, "resolve_dsmz_graph_neighbours.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p {P_CONS_ALL_DIR}
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --winners_tsv {input.winners_tsv} \
          --direction_summary_tsv {input.dir_summary_tsv} \
          --graph_root {TUMOUR_NH_ROOT_ABS} \
          --directions {params.directions} \
          --output_tsv {output.resolved_tsv} \
          --validation_tsv {output.validation_tsv} \
          > {log} 2>&1
        test -s {output.resolved_tsv} || (echo "ERROR: missing {output.resolved_tsv}" >&2; exit 1)
        test -s {output.validation_tsv} || (echo "ERROR: missing {output.validation_tsv}" >&2; exit 1)
        '''

# =============================================================================
# MULTI-REPRESENTATION SUPPORT NETWORKS
# =============================================================================
# Stage role: builds support networks by aggregating edges across feature-distance representations.

# Rule: build_multi_representation_majority_threshold_consensus_network
# Method role: analysis rule that retains edges meeting the configured majority threshold across representations.
# Flow: per-representation graph edges -> majority-threshold support network and edge-support table.
# Provides support-network plots and post-resolution support stratification.
rule build_multi_representation_majority_threshold_consensus_network:
    """
    Computes the shared feature--distance representation edge-support table and
    retains only edges satisfying the majority support threshold.
    """
    input:
        node_universe = P_CONS_SHORTNAMES_TSV,
        graph_edges = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        edges = P_CONS_MULTI_REP_MAJORITY_EDGES_TSV,
        edge_support = P_CONS_MULTI_REP_EDGE_SUPPORT_TSV
    params:
        script = os.path.join(BASE, "scripts", "build_consensus_from_direction_edgefiles.py"),
        tumour_nh_dir = TUMOUR_NH_ROOT_ABS,
        directions = ",".join(CONS_DIRECTIONS),
        min_support = SIMILARITY_CONSENSUS_MIN_SUPPORT,
        cohort = profile_name
    log: os.path.join(LOGROOT, "build_multi_representation_majority_threshold_consensus_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {params.script} \
          --tumour_nh_dir {params.tumour_nh_dir} \
          --out_edges {output.edges} \
          --out_support {output.edge_support} \
          --mode majority_threshold \
          --cohort {params.cohort} \
          --min_support {params.min_support} \
          --directions {params.directions} \
          --name_map {input.node_universe} \
          --node_universe {input.node_universe} \
          --require_short \
          > {log} 2>&1
        test -s {output.edges} || (echo "ERROR: missing {output.edges}" >&2; exit 1)
        test -s {output.edge_support} || (echo "ERROR: missing {output.edge_support}" >&2; exit 1)
        '''

# Rule: build_multi_representation_union_supported_edges_network
# Method role: analysis rule that retains every edge observed in at least one feature-distance representation.
# Flow: per-representation graph edges -> union support network and component sidecars.
# Provides a sensitivity view of support-network topology.
rule build_multi_representation_union_supported_edges_network:
    """
    Computes the same feature--distance representation support calculation and
    retains every edge observed in at least one configured representation.
    """
    input:
        node_universe = P_CONS_SHORTNAMES_TSV,
        graph_edges = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        edges = P_CONS_MULTI_REP_UNION_EDGES_TSV
    params:
        script = os.path.join(BASE, "scripts", "build_consensus_from_direction_edgefiles.py"),
        tumour_nh_dir = TUMOUR_NH_ROOT_ABS,
        directions = ",".join(CONS_DIRECTIONS),
        cohort = profile_name
    log: os.path.join(LOGROOT, "build_multi_representation_union_supported_edges_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {params.script} \
          --tumour_nh_dir {params.tumour_nh_dir} \
          --out_edges {output.edges} \
          --mode union_supported_edges \
          --cohort {params.cohort} \
          --min_support 1 \
          --directions {params.directions} \
          --name_map {input.node_universe} \
          --node_universe {input.node_universe} \
          --require_short \
          > {log} 2>&1
        test -s {output.edges} || (echo "ERROR: missing {output.edges}" >&2; exit 1)
        '''

# Rule: build_patient_referenced_support_threshold_consensus_cell_line_similarity_network
# Method role: analysis rule that builds the patient-referenced support-threshold network from representation edges.
# Flow: per-direction similarity edges -> consensus edge table, node statistics, and display-name sidecars.
# Provides the support-threshold final graph and multicohort inspection plots.
rule build_patient_referenced_support_threshold_consensus_cell_line_similarity_network:
    """
    Aggregates per-direction DSMZ similarity graph edge files into a single
    consensus edge table for the support-threshold consensus network figure.
    """
    input:
        shortnames = P_CONS_SHORTNAMES_TSV,
        graph_edges = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_edges_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        edges = P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV
    params:
        script = os.path.join(BASE, "scripts", "build_consensus_from_direction_edgefiles.py"),
        tumour_nh_dir = TUMOUR_NH_ROOT_ABS,
        directions = ",".join(CONS_DIRECTIONS),
        # Majority support threshold m = max(2, floor(|R| / 2) + 1),
        # derived from the configured active feature-distance representation list.
        min_support = SIMILARITY_CONSENSUS_MIN_SUPPORT
    log: os.path.join(LOGROOT, "build_patient_referenced_support_threshold_consensus_cell_line_similarity_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {params.script} \
          --tumour_nh_dir {params.tumour_nh_dir} \
          --out_edges {output.edges} \
          --min_support {params.min_support} \
          --directions {params.directions} \
          --name_map {input.shortnames} \
          --require_short \
          > {log} 2>&1
        test -s {output.edges} || (echo "ERROR: missing {output.edges}" >&2; exit 1)
        '''

# =============================================================================
# FINAL GRAPH PLOTTING
# =============================================================================
# Stage role: renders final and supplementary graph views without changing graph topology.

# Rule: derive_resolved_graph_node_statistics
# Scientific purpose: materialise graph-node analytical statistics before any
# figure rendering.
# Unit of analysis: resolved biological cell-line node in the patient-referenced
# graph. Transformation: resolved-neighbour table -> edge table, component and
# isolate labels, and anchor-centrality audit.
rule derive_resolved_graph_node_statistics:
    input:
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv"),
        script = os.path.join(BASE, "scripts", "derive_resolved_graph_node_statistics.py"),
        visual_helper = os.path.join(BASE, "scripts", "visualize_resolved_dsmz_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        shortnames = P_CONS_SHORTNAMES_TSV,
        edges = P_CONS_RESOLVED_EDGES_TSV,
        node_stats = P_CONS_RESOLVED_NODE_STATS_TSV,
        anchor_audit = P_CONS_ANCHOR_AUDIT_TSV
    params:
        cohort = profile_name.upper()
    log:
        os.path.join(LOGROOT, "derive_resolved_graph_node_statistics.log")
    conda:
        CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_GRAPH_NODE_STATS_DIR}
        python {input.script} \
          --resolved-neighbours {input.resolved_tsv} \
          --cohort {params.cohort} \
          --out-shortnames {output.shortnames} \
          --out-edges {output.edges} \
          --out-node-stats {output.node_stats} \
          --out-anchor-audit {output.anchor_audit} \
          > {log} 2>&1
        test -s {output.shortnames} || (echo "ERROR: missing {output.shortnames}" >&2; exit 1)
        test -s {output.edges} || (echo "ERROR: missing {output.edges}" >&2; exit 1)
        test -s {output.node_stats} || (echo "ERROR: missing {output.node_stats}" >&2; exit 1)
        test -s {output.anchor_audit} || (echo "ERROR: missing {output.anchor_audit}" >&2; exit 1)
        '''

# Rule: plot_patient_referenced_resolved_cell_line_neighbourhood_graph
# Method role: final plotting rule for the resolved-neighbour patient-referenced graph.
# Flow: resolved graph edges, node statistics, and display names -> PDF/PNG/SVG graph outputs.
# Purpose: reporting visualisation only; topology is fixed by the resolution rule.
rule plot_patient_referenced_resolved_cell_line_neighbourhood_graph:
    """
    Patient-referenced resolved-neighbourhood graph figure for the DSMZ
    similarity network after graph-based neighbour resolution. The helper also
    writes short-name, edge-list, and node-stat TSVs into the plot directory for
    figure provenance.
    """
    input:
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv"),
        script = os.path.join(BASE, "scripts", "visualize_resolved_dsmz_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py"),
        shortnames = P_CONS_SHORTNAMES_TSV,
        edges = P_CONS_RESOLVED_EDGES_TSV,
        node_stats = P_CONS_RESOLVED_NODE_STATS_TSV,
        anchor_audit = P_CONS_ANCHOR_AUDIT_TSV
    output:
        png = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".png",
        pdf = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".pdf",
        svg = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".svg",
        legend_png = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_legend.png",
        legend_pdf = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_legend.pdf",
        legend_svg = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_legend.svg"
    params:
        out_prefix = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX,
        legend_out_prefix = P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + "_legend",
        label = profile_name.upper(),
        fig_w = graph_layout_value("resolved_patient_graph", "fig_w"),
        fig_h = graph_layout_value("resolved_patient_graph", "fig_h"),
        dpi = graph_layout_value("resolved_patient_graph", "dpi"),
        node_size = graph_layout_value("resolved_patient_graph", "node_size"),
        font_size = graph_layout_value("resolved_patient_graph", "font_size"),
        edge_width_min = graph_layout_value("resolved_patient_graph", "edge_width_min"),
        edge_width_max = graph_layout_value("resolved_patient_graph", "edge_width_max"),
        legend_mode = graph_layout_value("resolved_patient_graph", "legend_mode"),
        use_halo_labels_arg = "--use-halo-labels" if graph_layout_value("resolved_patient_graph", "use_halo_labels") else "",
        halo_linewidth = graph_layout_value("resolved_patient_graph", "halo_linewidth"),
        halo_color = graph_layout_value("resolved_patient_graph", "halo_color"),
        dense_component_min_size = graph_layout_value("resolved_patient_graph", "dense_component_min_size"),
        dense_component_min_density = graph_layout_value("resolved_patient_graph", "dense_component_min_density"),
        dense_component_expand_factor = graph_layout_value("resolved_patient_graph", "dense_component_expand_factor"),
        adaptive_component_layout_arg = "--adaptive-component-layout" if graph_layout_value("resolved_patient_graph", "adaptive_component_layout") else "",
        adaptive_component_min_nodes = graph_layout_value("resolved_patient_graph", "adaptive_component_min_nodes"),
        adaptive_component_min_density = graph_layout_value("resolved_patient_graph", "adaptive_component_min_density"),
        adaptive_component_min_expand = graph_layout_value("resolved_patient_graph", "adaptive_component_min_expand"),
        adaptive_component_max_expand = graph_layout_value("resolved_patient_graph", "adaptive_component_max_expand"),
        adaptive_component_density_weight = graph_layout_value("resolved_patient_graph", "adaptive_component_density_weight"),
        adaptive_component_size_weight = graph_layout_value("resolved_patient_graph", "adaptive_component_size_weight"),
        adaptive_component_label_weight = graph_layout_value("resolved_patient_graph", "adaptive_component_label_weight"),
        position_fill_x = graph_layout_value("resolved_patient_graph", "position_fill_x"),
        position_fill_y = graph_layout_value("resolved_patient_graph", "position_fill_y"),
        axis_margin_frac = graph_layout_value("resolved_patient_graph", "axis_margin_frac"),
        label_overlap_avoidance_arg = "--label-overlap-avoidance" if graph_layout_value("resolved_patient_graph", "label_overlap_avoidance") else "",
        label_padding_factor = graph_layout_value("resolved_patient_graph", "label_padding_factor"),
        label_bbox_margin = graph_layout_value("resolved_patient_graph", "label_bbox_margin"),
        max_label_overlap_iterations = graph_layout_value("resolved_patient_graph", "max_label_overlap_iterations"),
        min_label_separation = graph_layout_value("resolved_patient_graph", "min_label_separation"),
        avoid_isolate_overlap_arg = "--avoid-isolate-overlap" if graph_layout_value("resolved_patient_graph", "avoid_isolate_overlap") else "",
        min_component_isolate_gap = graph_layout_value("resolved_patient_graph", "min_component_isolate_gap"),
        isolate_panel_reserved_height = graph_layout_value("resolved_patient_graph", "isolate_panel_reserved_height"),
        adaptive_isolate_y_shift_arg = "--adaptive-isolate-y-shift" if graph_layout_value("resolved_patient_graph", "adaptive_isolate_y_shift") else "",
        component_isolate_gap_factor = graph_layout_value("resolved_patient_graph", "component_isolate_gap_factor"),
        resolved_component_label_gap = graph_layout_value("resolved_patient_graph", "component_label_gap"),
        resolved_grid_cols = graph_layout_value("resolved_patient_graph", "component_grid_cols"),
        resolved_grid_rows = graph_layout_value("resolved_patient_graph", "component_grid_rows"),
        resolved_cell_width = graph_layout_value("resolved_patient_graph", "component_cell_width"),
        resolved_cell_height = graph_layout_value("resolved_patient_graph", "component_cell_height"),
        resolved_row_gap = graph_layout_value("resolved_patient_graph", "component_row_gap"),
        resolved_col_gap = graph_layout_value("resolved_patient_graph", "component_col_gap"),
        resolved_legend_width = graph_layout_value("resolved_patient_graph", "legend_width"),
        resolved_isolate_height = graph_layout_value("resolved_patient_graph", "isolate_region_height"),
        resolved_isolate_spacing = graph_layout_value("resolved_patient_graph", "isolate_spacing"),
        resolved_isolate_panel_padding_x = graph_layout_value("resolved_patient_graph", "isolate_panel_padding_x"),
        resolved_isolate_panel_padding_y = graph_layout_value("resolved_patient_graph", "isolate_panel_padding_y"),
        resolved_isolate_label_band_frac = graph_layout_value("resolved_patient_graph", "isolate_label_band_frac"),
        resolved_isolate_max_per_row = graph_layout_value("resolved_patient_graph", "isolate_max_per_row"),
        resolved_isolate_panel_width_frac = graph_layout_value("resolved_patient_graph", "isolate_panel_width_frac"),
        resolved_footnote_height = graph_layout_value("resolved_patient_graph", "footnote_region_height"),
        resolved_component_cell_fill_x = graph_layout_value("resolved_patient_graph", "component_cell_fill_x"),
        resolved_component_cell_fill_y = graph_layout_value("resolved_patient_graph", "component_cell_fill_y"),
        resolved_dense_component_cell_fill_x = graph_layout_value("resolved_patient_graph", "dense_component_cell_fill_x"),
        resolved_dense_component_cell_fill_y = graph_layout_value("resolved_patient_graph", "dense_component_cell_fill_y"),
        resolved_c2_component_cell_fill_x = graph_layout_value("resolved_patient_graph", "c2_component_cell_fill_x"),
        resolved_c2_component_cell_fill_y = graph_layout_value("resolved_patient_graph", "c2_component_cell_fill_y"),
        c1_x_expand_factor = graph_layout_value("resolved_patient_graph", "c1_x_expand_factor"),
        c1_y_expand_factor = graph_layout_value("resolved_patient_graph", "c1_y_expand_factor"),
        c2_x_expand_factor = graph_layout_value("resolved_patient_graph", "c2_x_expand_factor"),
        c2_y_expand_factor = graph_layout_value("resolved_patient_graph", "c2_y_expand_factor"),
        isolate_gap = graph_layout_value("resolved_patient_graph", "isolate_gap"),
        isolate_box_height_factor = graph_layout_value("resolved_patient_graph", "isolate_box_height_factor"),
        isolate_label_y_frac = graph_layout_value("resolved_patient_graph", "isolate_label_y_frac"),
        isolate_label_left_pad = graph_layout_value("resolved_patient_graph", "isolate_label_left_pad"),
        compact_isolate_panel_arg = "--compact-isolate-panel" if graph_layout_value("resolved_patient_graph", "compact_isolate_panel") else "",
        isolate_panel_width_mode = graph_layout_value("resolved_patient_graph", "isolate_panel_width_mode"),
        isolate_panel_spacing = graph_layout_value("resolved_patient_graph", "isolate_panel_spacing"),
        min_isolate_label_spacing = graph_layout_value("resolved_patient_graph", "min_isolate_label_spacing"),
        legend_alignment_mode = graph_layout_value("resolved_patient_graph", "legend_alignment_mode")
    log: os.path.join(LOGROOT, "plot_patient_referenced_resolved_cell_line_neighbourhood_graph.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        test -s {input.shortnames} || (echo "ERROR: missing analytical short-name table {input.shortnames}" >&2; exit 1)
        test -s {input.edges} || (echo "ERROR: missing analytical resolved-edge table {input.edges}" >&2; exit 1)
        test -s {input.node_stats} || (echo "ERROR: missing analytical node-statistics table {input.node_stats}" >&2; exit 1)
        test -s {input.anchor_audit} || (echo "ERROR: missing analytical anchor audit {input.anchor_audit}" >&2; exit 1)
        python {input.script} \
          {input.resolved_tsv} \
          {params.out_prefix} \
          {params.label} \
          --fig-w {params.fig_w} \
          --fig-h {params.fig_h} \
          --dpi {params.dpi} \
          --node-size {params.node_size} \
          --font-size {params.font_size} \
          --edge-width-min {params.edge_width_min} \
          --edge-width-max {params.edge_width_max} \
          --legend-mode {params.legend_mode} \
          --legend-out-prefix {params.legend_out_prefix} \
          {params.use_halo_labels_arg} \
          --halo-linewidth {params.halo_linewidth} \
          --halo-color {params.halo_color} \
          --dense-component-min-size {params.dense_component_min_size} \
          --dense-component-min-density {params.dense_component_min_density} \
          --dense-component-expand-factor {params.dense_component_expand_factor} \
          {params.adaptive_component_layout_arg} \
          --adaptive-component-min-nodes {params.adaptive_component_min_nodes} \
          --adaptive-component-min-density {params.adaptive_component_min_density} \
          --adaptive-component-min-expand {params.adaptive_component_min_expand} \
          --adaptive-component-max-expand {params.adaptive_component_max_expand} \
          --adaptive-component-density-weight {params.adaptive_component_density_weight} \
          --adaptive-component-size-weight {params.adaptive_component_size_weight} \
          --adaptive-component-label-weight {params.adaptive_component_label_weight} \
          --position-fill-x {params.position_fill_x} \
          --position-fill-y {params.position_fill_y} \
          --axis-margin-frac {params.axis_margin_frac} \
          {params.label_overlap_avoidance_arg} \
          --label-padding-factor {params.label_padding_factor} \
          --label-bbox-margin {params.label_bbox_margin} \
          --max-label-overlap-iterations {params.max_label_overlap_iterations} \
          --min-label-separation {params.min_label_separation} \
          {params.avoid_isolate_overlap_arg} \
          --min-component-isolate-gap {params.min_component_isolate_gap} \
          --isolate-panel-reserved-height {params.isolate_panel_reserved_height} \
          {params.adaptive_isolate_y_shift_arg} \
          --component-isolate-gap-factor {params.component_isolate_gap_factor} \
          --component-label-gap {params.resolved_component_label_gap} \
          --component-grid-ncols {params.resolved_grid_cols} \
          --component-grid-nrows {params.resolved_grid_rows} \
          --component-cell-width {params.resolved_cell_width} \
          --component-cell-height {params.resolved_cell_height} \
          --component-row-gap {params.resolved_row_gap} \
          --component-col-gap {params.resolved_col_gap} \
          --legend-width {params.resolved_legend_width} \
          --isolate-region-height {params.resolved_isolate_height} \
          --footnote-region-height {params.resolved_footnote_height} \
          --component-cell-fill-x {params.resolved_component_cell_fill_x} \
          --component-cell-fill-y {params.resolved_component_cell_fill_y} \
          --dense-component-cell-fill-x {params.resolved_dense_component_cell_fill_x} \
          --dense-component-cell-fill-y {params.resolved_dense_component_cell_fill_y} \
          --c2-component-cell-fill-x {params.resolved_c2_component_cell_fill_x} \
          --c2-component-cell-fill-y {params.resolved_c2_component_cell_fill_y} \
          --c1-x-expand-factor {params.c1_x_expand_factor} \
          --c1-y-expand-factor {params.c1_y_expand_factor} \
          --c2-x-expand-factor {params.c2_x_expand_factor} \
          --c2-y-expand-factor {params.c2_y_expand_factor} \
          --isolate-gap {params.isolate_gap} \
          --isolate-spacing {params.resolved_isolate_spacing} \
          --isolate-panel-padding-x {params.resolved_isolate_panel_padding_x} \
          --isolate-panel-padding-y {params.resolved_isolate_panel_padding_y} \
          --isolate-label-band-frac {params.resolved_isolate_label_band_frac} \
          --isolate-max-per-row {params.resolved_isolate_max_per_row} \
          --isolate-panel-width-frac {params.resolved_isolate_panel_width_frac} \
          --isolate-box-height-factor {params.isolate_box_height_factor} \
          --isolate-label-y-frac {params.isolate_label_y_frac} \
          --isolate-label-left-pad {params.isolate_label_left_pad} \
          {params.compact_isolate_panel_arg} \
          --isolate-panel-width-mode {params.isolate_panel_width_mode} \
          --isolate-panel-spacing {params.isolate_panel_spacing} \
          --min-isolate-label-spacing {params.min_isolate_label_spacing} \
          --legend-alignment-mode {params.legend_alignment_mode} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.legend_pdf} || (echo "ERROR: missing {output.legend_pdf}" >&2; exit 1)
        test -s {output.legend_png} || (echo "ERROR: missing {output.legend_png}" >&2; exit 1)
        test -s {output.legend_svg} || (echo "ERROR: missing {output.legend_svg}" >&2; exit 1)
        '''

# Rule: plot_multi_representation_majority_threshold_consensus_network
# Method role: final plotting rule for the majority-threshold support network.
# Flow: support-network edges and node sidecars -> PDF/PNG/SVG graph outputs.
# Purpose: visualises representation support without changing support thresholds.
rule plot_multi_representation_majority_threshold_consensus_network:
    """
    Plots the majority-threshold consensus network and writes node/component
    sidecars from the full cohort node universe.
    """
    input:
        edges = P_CONS_MULTI_REP_MAJORITY_EDGES_TSV,
        node_universe = P_CONS_SHORTNAMES_TSV,
        script = os.path.join(BASE, "scripts", "plot_consensus_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        pdf = P_CONS_MULTI_REP_MAJORITY_PREFIX + ".pdf",
        png = P_CONS_MULTI_REP_MAJORITY_PREFIX + ".png",
        svg = P_CONS_MULTI_REP_MAJORITY_PREFIX + ".svg",
        legend_pdf = P_CONS_MULTI_REP_MAJORITY_PREFIX + "_legend.pdf",
        legend_png = P_CONS_MULTI_REP_MAJORITY_PREFIX + "_legend.png",
        legend_svg = P_CONS_MULTI_REP_MAJORITY_PREFIX + "_legend.svg",
        node_stats = P_CONS_MULTI_REP_MAJORITY_NODE_STATS_TSV,
        components = P_CONS_MULTI_REP_MAJORITY_COMPONENTS_TSV
    params:
        out_prefix = P_CONS_MULTI_REP_MAJORITY_PREFIX,
        legend_out_prefix = P_CONS_MULTI_REP_MAJORITY_PREFIX + "_legend",
        fig_w = graph_layout_value("multi_representation_majority_graph", "fig_w"),
        fig_h = graph_layout_value("multi_representation_majority_graph", "fig_h"),
        dpi = graph_layout_value("multi_representation_majority_graph", "dpi"),
        node_size = graph_layout_value("multi_representation_majority_graph", "node_size"),
        font_size = graph_layout_value("multi_representation_majority_graph", "font_size"),
        edge_width_min = graph_layout_value("multi_representation_majority_graph", "edge_width_min"),
        edge_width_max = graph_layout_value("multi_representation_majority_graph", "edge_width_max"),
        legend_mode = graph_layout_value("multi_representation_majority_graph", "legend_mode"),
        use_halo_labels_arg = "--use-halo-labels" if graph_layout_value("multi_representation_majority_graph", "use_halo_labels") else "",
        halo_linewidth = graph_layout_value("multi_representation_majority_graph", "halo_linewidth"),
        halo_color = graph_layout_value("multi_representation_majority_graph", "halo_color"),
        match_resolved_style_arg = "--match-resolved-style" if graph_layout_value("multi_representation_majority_graph", "match_resolved_style") else "",
        position_fill_x = graph_layout_value("multi_representation_majority_graph", "position_fill_x"),
        position_fill_y = graph_layout_value("multi_representation_majority_graph", "position_fill_y"),
        axis_margin_frac = graph_layout_value("multi_representation_majority_graph", "axis_margin_frac"),
        compact_isolate_panel_arg = "--compact-isolate-panel" if graph_layout_value("multi_representation_majority_graph", "compact_isolate_panel") else "",
        isolate_panel_width_mode = graph_layout_value("multi_representation_majority_graph", "isolate_panel_width_mode"),
        isolate_panel_spacing = graph_layout_value("multi_representation_majority_graph", "isolate_panel_spacing"),
        isolate_panel_padding_x = graph_layout_value("multi_representation_majority_graph", "isolate_panel_padding_x"),
        isolate_panel_padding_y = graph_layout_value("multi_representation_majority_graph", "isolate_panel_padding_y"),
        min_isolate_label_spacing = graph_layout_value("multi_representation_majority_graph", "min_isolate_label_spacing"),
        adaptive_component_layout_arg = "--adaptive-component-layout" if graph_layout_value("multi_representation_majority_graph", "adaptive_component_layout") else "",
        adaptive_component_min_nodes = graph_layout_value("multi_representation_majority_graph", "adaptive_component_min_nodes"),
        adaptive_component_min_density = graph_layout_value("multi_representation_majority_graph", "adaptive_component_min_density"),
        adaptive_component_min_expand = graph_layout_value("multi_representation_majority_graph", "adaptive_component_min_expand"),
        adaptive_component_max_expand = graph_layout_value("multi_representation_majority_graph", "adaptive_component_max_expand"),
        adaptive_component_density_weight = graph_layout_value("multi_representation_majority_graph", "adaptive_component_density_weight"),
        adaptive_component_size_weight = graph_layout_value("multi_representation_majority_graph", "adaptive_component_size_weight"),
        adaptive_component_label_weight = graph_layout_value("multi_representation_majority_graph", "adaptive_component_label_weight"),
        label_overlap_avoidance_arg = "--label-overlap-avoidance" if graph_layout_value("multi_representation_majority_graph", "label_overlap_avoidance") else "",
        label_padding_factor = graph_layout_value("multi_representation_majority_graph", "label_padding_factor"),
        label_bbox_margin = graph_layout_value("multi_representation_majority_graph", "label_bbox_margin"),
        max_label_overlap_iterations = graph_layout_value("multi_representation_majority_graph", "max_label_overlap_iterations"),
        min_label_separation = graph_layout_value("multi_representation_majority_graph", "min_label_separation"),
        avoid_isolate_overlap_arg = "--avoid-isolate-overlap" if graph_layout_value("multi_representation_majority_graph", "avoid_isolate_overlap") else "",
        min_component_isolate_gap = graph_layout_value("multi_representation_majority_graph", "min_component_isolate_gap"),
        isolate_panel_reserved_height = graph_layout_value("multi_representation_majority_graph", "isolate_panel_reserved_height"),
        adaptive_isolate_y_shift_arg = "--adaptive-isolate-y-shift" if graph_layout_value("multi_representation_majority_graph", "adaptive_isolate_y_shift") else "",
        component_isolate_gap_factor = graph_layout_value("multi_representation_majority_graph", "component_isolate_gap_factor"),
        label = profile_name.upper()
    log: os.path.join(LOGROOT, "plot_multi_representation_majority_threshold_consensus_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {input.script} \
          --edges {input.edges} \
          --nodes {input.node_universe} \
          --nodes-col short_id \
          --out {params.out_prefix} \
          --label {params.label} \
          --graph-mode majority_threshold \
          --fig-w {params.fig_w} \
          --fig-h {params.fig_h} \
          --dpi {params.dpi} \
          {params.use_halo_labels_arg} \
          --halo-linewidth {params.halo_linewidth} \
          --halo-color {params.halo_color} \
          {params.match_resolved_style_arg} \
          --position-fill-x {params.position_fill_x} \
          --position-fill-y {params.position_fill_y} \
          --axis-margin-frac {params.axis_margin_frac} \
          {params.adaptive_component_layout_arg} \
          --adaptive-component-min-nodes {params.adaptive_component_min_nodes} \
          --adaptive-component-min-density {params.adaptive_component_min_density} \
          --adaptive-component-min-expand {params.adaptive_component_min_expand} \
          --adaptive-component-max-expand {params.adaptive_component_max_expand} \
          --adaptive-component-density-weight {params.adaptive_component_density_weight} \
          --adaptive-component-size-weight {params.adaptive_component_size_weight} \
          --adaptive-component-label-weight {params.adaptive_component_label_weight} \
          {params.label_overlap_avoidance_arg} \
          --label-padding-factor {params.label_padding_factor} \
          --label-bbox-margin {params.label_bbox_margin} \
          --max-label-overlap-iterations {params.max_label_overlap_iterations} \
          --min-label-separation {params.min_label_separation} \
          {params.avoid_isolate_overlap_arg} \
          --min-component-isolate-gap {params.min_component_isolate_gap} \
          --isolate-panel-reserved-height {params.isolate_panel_reserved_height} \
          {params.adaptive_isolate_y_shift_arg} \
          --component-isolate-gap-factor {params.component_isolate_gap_factor} \
          --node-size {params.node_size} \
          --font-size {params.font_size} \
          --edge-width-min {params.edge_width_min} \
          --edge-width-max {params.edge_width_max} \
          {params.compact_isolate_panel_arg} \
          --isolate-panel-width-mode {params.isolate_panel_width_mode} \
          --isolate-panel-spacing {params.isolate_panel_spacing} \
          --isolate-panel-padding-x {params.isolate_panel_padding_x} \
          --isolate-panel-padding-y {params.isolate_panel_padding_y} \
          --min-isolate-label-spacing {params.min_isolate_label_spacing} \
          --legend-mode {params.legend_mode} \
          --legend-out-prefix {params.legend_out_prefix} \
          --node-stats-out {output.node_stats} \
          --components-out {output.components} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.png} || (echo "ERROR: missing {output.png}" >&2; exit 1)
        test -s {output.svg} || (echo "ERROR: missing {output.svg}" >&2; exit 1)
        test -s {output.legend_pdf} || (echo "ERROR: missing {output.legend_pdf}" >&2; exit 1)
        test -s {output.legend_png} || (echo "ERROR: missing {output.legend_png}" >&2; exit 1)
        test -s {output.legend_svg} || (echo "ERROR: missing {output.legend_svg}" >&2; exit 1)
        test -s {output.node_stats} || (echo "ERROR: missing {output.node_stats}" >&2; exit 1)
        test -s {output.components} || (echo "ERROR: missing {output.components}" >&2; exit 1)
        '''

# Rule: plot_multi_representation_union_supported_edges_network
# Method role: supplementary plotting rule for the union support network sensitivity view.
# Flow: union network edges and node sidecars -> PDF/PNG/SVG graph outputs.
# Purpose: visualises union-supported support-network structure alongside majority-threshold outputs.
rule plot_multi_representation_union_supported_edges_network:
    """
    Plots all supported edges. Edges supported by one feature--distance
    representation are dashed; edges supported by multiple representations are solid.
    """
    input:
        edges = P_CONS_MULTI_REP_UNION_EDGES_TSV,
        node_universe = P_CONS_SHORTNAMES_TSV,
        script = os.path.join(BASE, "scripts", "plot_consensus_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        pdf = P_CONS_MULTI_REP_UNION_PREFIX + ".pdf",
        png = P_CONS_MULTI_REP_UNION_PREFIX + ".png",
        svg = P_CONS_MULTI_REP_UNION_PREFIX + ".svg",
        legend_pdf = P_CONS_MULTI_REP_UNION_PREFIX + "_legend.pdf",
        legend_png = P_CONS_MULTI_REP_UNION_PREFIX + "_legend.png",
        legend_svg = P_CONS_MULTI_REP_UNION_PREFIX + "_legend.svg",
        node_stats = P_CONS_MULTI_REP_UNION_NODE_STATS_TSV,
        components = P_CONS_MULTI_REP_UNION_COMPONENTS_TSV,
        full_node_labels = P_CONS_MULTI_REP_UNION_FULL_NODE_LABELS_TSV
    params:
        out_prefix = P_CONS_MULTI_REP_UNION_PREFIX,
        legend_out_prefix = P_CONS_MULTI_REP_UNION_PREFIX + "_legend",
        fig_w = graph_layout_value("multi_representation_union_graph", "fig_w"),
        fig_h = graph_layout_value("multi_representation_union_graph", "fig_h"),
        dpi = graph_layout_value("multi_representation_union_graph", "dpi"),
        node_size = graph_layout_value("multi_representation_union_graph", "node_size"),
        font_size = graph_layout_value("multi_representation_union_graph", "font_size"),
        node_label_mode = graph_layout_value("multi_representation_union_graph", "node_label_mode"),
        selected_label_max_component_labels = graph_layout_value("multi_representation_union_graph", "selected_label_max_component_labels"),
        selected_label_min_component_nodes = graph_layout_value("multi_representation_union_graph", "selected_label_min_component_nodes"),
        selected_label_min_component_density = graph_layout_value("multi_representation_union_graph", "selected_label_min_component_density"),
        selected_label_min_distance = graph_layout_value("multi_representation_union_graph", "selected_label_min_distance"),
        edge_width_min = graph_layout_value("multi_representation_union_graph", "edge_width_min"),
        edge_width_max = graph_layout_value("multi_representation_union_graph", "edge_width_max"),
        legend_mode = graph_layout_value("multi_representation_union_graph", "legend_mode"),
        use_halo_labels_arg = "--use-halo-labels" if graph_layout_value("multi_representation_union_graph", "use_halo_labels") else "",
        halo_linewidth = graph_layout_value("multi_representation_union_graph", "halo_linewidth"),
        halo_color = graph_layout_value("multi_representation_union_graph", "halo_color"),
        match_resolved_style_arg = "--match-resolved-style" if graph_layout_value("multi_representation_union_graph", "match_resolved_style") else "",
        position_fill_x = graph_layout_value("multi_representation_union_graph", "position_fill_x"),
        position_fill_y = graph_layout_value("multi_representation_union_graph", "position_fill_y"),
        axis_margin_frac = graph_layout_value("multi_representation_union_graph", "axis_margin_frac"),
        compact_isolate_panel_arg = "--compact-isolate-panel" if graph_layout_value("multi_representation_union_graph", "compact_isolate_panel") else "",
        isolate_panel_width_mode = graph_layout_value("multi_representation_union_graph", "isolate_panel_width_mode"),
        isolate_panel_spacing = graph_layout_value("multi_representation_union_graph", "isolate_panel_spacing"),
        isolate_panel_padding_x = graph_layout_value("multi_representation_union_graph", "isolate_panel_padding_x"),
        isolate_panel_padding_y = graph_layout_value("multi_representation_union_graph", "isolate_panel_padding_y"),
        min_isolate_label_spacing = graph_layout_value("multi_representation_union_graph", "min_isolate_label_spacing"),
        adaptive_component_layout_arg = "--adaptive-component-layout" if graph_layout_value("multi_representation_union_graph", "adaptive_component_layout") else "",
        adaptive_component_min_nodes = graph_layout_value("multi_representation_union_graph", "adaptive_component_min_nodes"),
        adaptive_component_min_density = graph_layout_value("multi_representation_union_graph", "adaptive_component_min_density"),
        adaptive_component_min_expand = graph_layout_value("multi_representation_union_graph", "adaptive_component_min_expand"),
        adaptive_component_max_expand = graph_layout_value("multi_representation_union_graph", "adaptive_component_max_expand"),
        adaptive_component_density_weight = graph_layout_value("multi_representation_union_graph", "adaptive_component_density_weight"),
        adaptive_component_size_weight = graph_layout_value("multi_representation_union_graph", "adaptive_component_size_weight"),
        adaptive_component_label_weight = graph_layout_value("multi_representation_union_graph", "adaptive_component_label_weight"),
        label_overlap_avoidance_arg = "--label-overlap-avoidance" if graph_layout_value("multi_representation_union_graph", "label_overlap_avoidance") else "",
        label_padding_factor = graph_layout_value("multi_representation_union_graph", "label_padding_factor"),
        label_bbox_margin = graph_layout_value("multi_representation_union_graph", "label_bbox_margin"),
        max_label_overlap_iterations = graph_layout_value("multi_representation_union_graph", "max_label_overlap_iterations"),
        min_label_separation = graph_layout_value("multi_representation_union_graph", "min_label_separation"),
        avoid_isolate_overlap_arg = "--avoid-isolate-overlap" if graph_layout_value("multi_representation_union_graph", "avoid_isolate_overlap") else "",
        min_component_isolate_gap = graph_layout_value("multi_representation_union_graph", "min_component_isolate_gap"),
        isolate_panel_reserved_height = graph_layout_value("multi_representation_union_graph", "isolate_panel_reserved_height"),
        adaptive_isolate_y_shift_arg = "--adaptive-isolate-y-shift" if graph_layout_value("multi_representation_union_graph", "adaptive_isolate_y_shift") else "",
        component_isolate_gap_factor = graph_layout_value("multi_representation_union_graph", "component_isolate_gap_factor"),
        label = profile_name.upper()
    log: os.path.join(LOGROOT, "plot_multi_representation_union_supported_edges_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {input.script} \
          --edges {input.edges} \
          --nodes {input.node_universe} \
          --nodes-col short_id \
          --out {params.out_prefix} \
          --label {params.label} \
          --graph-mode union_supported_edges \
          --node-label-mode {params.node_label_mode} \
          --selected-label-max-component-labels {params.selected_label_max_component_labels} \
          --selected-label-min-component-nodes {params.selected_label_min_component_nodes} \
          --selected-label-min-component-density {params.selected_label_min_component_density} \
          --selected-label-min-distance {params.selected_label_min_distance} \
          --fig-w {params.fig_w} \
          --fig-h {params.fig_h} \
          --dpi {params.dpi} \
          {params.use_halo_labels_arg} \
          --halo-linewidth {params.halo_linewidth} \
          --halo-color {params.halo_color} \
          {params.match_resolved_style_arg} \
          --position-fill-x {params.position_fill_x} \
          --position-fill-y {params.position_fill_y} \
          --axis-margin-frac {params.axis_margin_frac} \
          {params.adaptive_component_layout_arg} \
          --adaptive-component-min-nodes {params.adaptive_component_min_nodes} \
          --adaptive-component-min-density {params.adaptive_component_min_density} \
          --adaptive-component-min-expand {params.adaptive_component_min_expand} \
          --adaptive-component-max-expand {params.adaptive_component_max_expand} \
          --adaptive-component-density-weight {params.adaptive_component_density_weight} \
          --adaptive-component-size-weight {params.adaptive_component_size_weight} \
          --adaptive-component-label-weight {params.adaptive_component_label_weight} \
          {params.label_overlap_avoidance_arg} \
          --label-padding-factor {params.label_padding_factor} \
          --label-bbox-margin {params.label_bbox_margin} \
          --max-label-overlap-iterations {params.max_label_overlap_iterations} \
          --min-label-separation {params.min_label_separation} \
          {params.avoid_isolate_overlap_arg} \
          --min-component-isolate-gap {params.min_component_isolate_gap} \
          --isolate-panel-reserved-height {params.isolate_panel_reserved_height} \
          {params.adaptive_isolate_y_shift_arg} \
          --component-isolate-gap-factor {params.component_isolate_gap_factor} \
          --node-size {params.node_size} \
          --font-size {params.font_size} \
          --edge-width-min {params.edge_width_min} \
          --edge-width-max {params.edge_width_max} \
          {params.compact_isolate_panel_arg} \
          --isolate-panel-width-mode {params.isolate_panel_width_mode} \
          --isolate-panel-spacing {params.isolate_panel_spacing} \
          --isolate-panel-padding-x {params.isolate_panel_padding_x} \
          --isolate-panel-padding-y {params.isolate_panel_padding_y} \
          --min-isolate-label-spacing {params.min_isolate_label_spacing} \
          --legend-mode {params.legend_mode} \
          --legend-out-prefix {params.legend_out_prefix} \
          --node-labels-out {output.full_node_labels} \
          --node-stats-out {output.node_stats} \
          --components-out {output.components} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.png} || (echo "ERROR: missing {output.png}" >&2; exit 1)
        test -s {output.svg} || (echo "ERROR: missing {output.svg}" >&2; exit 1)
        test -s {output.legend_pdf} || (echo "ERROR: missing {output.legend_pdf}" >&2; exit 1)
        test -s {output.legend_png} || (echo "ERROR: missing {output.legend_png}" >&2; exit 1)
        test -s {output.legend_svg} || (echo "ERROR: missing {output.legend_svg}" >&2; exit 1)
        test -s {output.node_stats} || (echo "ERROR: missing {output.node_stats}" >&2; exit 1)
        test -s {output.components} || (echo "ERROR: missing {output.components}" >&2; exit 1)
        test -s {output.full_node_labels} || (echo "ERROR: missing {output.full_node_labels}" >&2; exit 1)
        '''

# Rule: plot_patient_referenced_support_threshold_consensus_cell_line_similarity_network
# Method role: final plotting rule for the patient-referenced support-threshold consensus network.
# Flow: consensus edge table, display names, and metadata -> graph figure and legend outputs.
# Purpose: reporting visualisation only; network membership is defined by earlier rules.
rule plot_patient_referenced_support_threshold_consensus_cell_line_similarity_network:
    """
    Support-threshold consensus network figure for the DSMZ similarity network
    aggregated across configured similarity-network directions.

    The plotting interface calls every threshold-filtered graph
    ``majority_threshold``; ``support_threshold_consensus_graph`` is the
    separate layout-config namespace, not a valid ``--graph-mode`` value.
    """
    input:
        edges = P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV,
        shortnames = P_CONS_SHORTNAMES_TSV,
        script = os.path.join(BASE, "scripts", "plot_consensus_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        png = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + ".png",
        pdf = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + ".pdf",
        svg = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + ".svg",
        legend_png = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_legend.png",
        legend_pdf = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_legend.pdf",
        legend_svg = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_legend.svg"
    params:
        out_prefix = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX,
        legend_out_prefix = P_CONS_SUPPORT_CONSENSUS_NETWORK_PREFIX + "_legend",
        label = profile_name.upper(),
        graph_title = graph_layout_value("support_threshold_consensus_graph", "graph_title"),
        fig_w = graph_layout_value("support_threshold_consensus_graph", "fig_w"),
        fig_h = graph_layout_value("support_threshold_consensus_graph", "fig_h"),
        dpi = graph_layout_value("support_threshold_consensus_graph", "dpi"),
        node_size = graph_layout_value("support_threshold_consensus_graph", "node_size"),
        font_size = graph_layout_value("support_threshold_consensus_graph", "font_size"),
        edge_width_min = graph_layout_value("support_threshold_consensus_graph", "edge_width_min"),
        edge_width_max = graph_layout_value("support_threshold_consensus_graph", "edge_width_max"),
        component_label_gap = graph_layout_value("support_threshold_consensus_graph", "component_label_gap"),
        nbl_component_gap_fraction = graph_layout_value("support_threshold_consensus_graph", "nbl_component_gap_fraction"),
        nbl_component_min_gap = graph_layout_value("support_threshold_consensus_graph", "nbl_component_min_gap"),
        legend_anchor_x = graph_layout_value("support_threshold_consensus_graph", "legend_anchor_x"),
        legend_anchor_y_arg = "" if graph_layout_value("support_threshold_consensus_graph", "legend_anchor_y") is None else f"--legend-anchor-y {float(graph_layout_value('support_threshold_consensus_graph', 'legend_anchor_y'))}",
        legend_align_mode = graph_layout_value("support_threshold_consensus_graph", "legend_align_mode"),
        legend_gap_frac = graph_layout_value("support_threshold_consensus_graph", "legend_gap_frac"),
        legend_mode = graph_layout_value("support_threshold_consensus_graph", "legend_mode"),
        use_halo_labels_arg = "--use-halo-labels" if graph_layout_value("support_threshold_consensus_graph", "use_halo_labels") else "",
        halo_linewidth = graph_layout_value("support_threshold_consensus_graph", "halo_linewidth"),
        halo_color = graph_layout_value("support_threshold_consensus_graph", "halo_color"),
        match_resolved_style_arg = "--match-resolved-style" if graph_layout_value("support_threshold_consensus_graph", "match_resolved_style") else "",
        position_fill_x = graph_layout_value("support_threshold_consensus_graph", "position_fill_x"),
        position_fill_y = graph_layout_value("support_threshold_consensus_graph", "position_fill_y"),
        axis_margin_frac = graph_layout_value("support_threshold_consensus_graph", "axis_margin_frac"),
        compact_isolate_panel_arg = "--compact-isolate-panel" if graph_layout_value("support_threshold_consensus_graph", "compact_isolate_panel") else "",
        isolate_panel_width_mode = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_width_mode"),
        isolate_panel_spacing = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_spacing"),
        min_isolate_label_spacing = graph_layout_value("support_threshold_consensus_graph", "min_isolate_label_spacing"),
        isolate_label_y_frac = graph_layout_value("support_threshold_consensus_graph", "isolate_label_y_frac"),
        isolate_box_height_factor = graph_layout_value("support_threshold_consensus_graph", "isolate_box_height_factor"),
        isolate_gap = graph_layout_value("support_threshold_consensus_graph", "isolate_gap"),
        isolate_label_left_pad = graph_layout_value("support_threshold_consensus_graph", "isolate_label_left_pad"),
        isolate_spacing = graph_layout_value("support_threshold_consensus_graph", "isolate_spacing"),
        isolate_panel_padding_x = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_padding_x"),
        isolate_panel_padding_y = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_padding_y"),
        isolate_label_band_frac = graph_layout_value("support_threshold_consensus_graph", "isolate_label_band_frac"),
        isolate_max_per_row = graph_layout_value("support_threshold_consensus_graph", "isolate_max_per_row"),
        isolate_panel_width_frac = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_width_frac"),
        dense_component_min_size = graph_layout_value("support_threshold_consensus_graph", "dense_component_min_size"),
        dense_component_min_density = graph_layout_value("support_threshold_consensus_graph", "dense_component_min_density"),
        dense_component_expand_factor = graph_layout_value("support_threshold_consensus_graph", "dense_component_expand_factor"),
        adaptive_component_layout_arg = "--adaptive-component-layout" if graph_layout_value("support_threshold_consensus_graph", "adaptive_component_layout") else "",
        adaptive_component_min_nodes = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_min_nodes"),
        adaptive_component_min_density = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_min_density"),
        adaptive_component_min_expand = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_min_expand"),
        adaptive_component_max_expand = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_max_expand"),
        adaptive_component_density_weight = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_density_weight"),
        adaptive_component_size_weight = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_size_weight"),
        adaptive_component_label_weight = graph_layout_value("support_threshold_consensus_graph", "adaptive_component_label_weight"),
        label_overlap_avoidance_arg = "--label-overlap-avoidance" if graph_layout_value("support_threshold_consensus_graph", "label_overlap_avoidance") else "",
        label_padding_factor = graph_layout_value("support_threshold_consensus_graph", "label_padding_factor"),
        label_bbox_margin = graph_layout_value("support_threshold_consensus_graph", "label_bbox_margin"),
        max_label_overlap_iterations = graph_layout_value("support_threshold_consensus_graph", "max_label_overlap_iterations"),
        min_label_separation = graph_layout_value("support_threshold_consensus_graph", "min_label_separation"),
        avoid_isolate_overlap_arg = "--avoid-isolate-overlap" if graph_layout_value("support_threshold_consensus_graph", "avoid_isolate_overlap") else "",
        min_component_isolate_gap = graph_layout_value("support_threshold_consensus_graph", "min_component_isolate_gap"),
        isolate_panel_reserved_height = graph_layout_value("support_threshold_consensus_graph", "isolate_panel_reserved_height"),
        adaptive_isolate_y_shift_arg = "--adaptive-isolate-y-shift" if graph_layout_value("support_threshold_consensus_graph", "adaptive_isolate_y_shift") else "",
        component_isolate_gap_factor = graph_layout_value("support_threshold_consensus_graph", "component_isolate_gap_factor"),
        support_grid_cols = graph_layout_value("support_threshold_consensus_graph", "component_grid_cols"),
        support_grid_rows = graph_layout_value("support_threshold_consensus_graph", "component_grid_rows"),
        support_cell_width = graph_layout_value("support_threshold_consensus_graph", "component_cell_width"),
        support_cell_height = graph_layout_value("support_threshold_consensus_graph", "component_cell_height"),
        support_row_gap = graph_layout_value("support_threshold_consensus_graph", "component_row_gap"),
        support_col_gap = graph_layout_value("support_threshold_consensus_graph", "component_col_gap"),
        support_legend_width = graph_layout_value("support_threshold_consensus_graph", "legend_width"),
        support_isolate_height = graph_layout_value("support_threshold_consensus_graph", "isolate_region_height"),
        support_footnote_height = graph_layout_value("support_threshold_consensus_graph", "footnote_region_height"),
        support_component_cell_fill_x = graph_layout_value("support_threshold_consensus_graph", "component_cell_fill_x"),
        support_component_cell_fill_y = graph_layout_value("support_threshold_consensus_graph", "component_cell_fill_y"),
        support_dense_component_cell_fill_x = graph_layout_value("support_threshold_consensus_graph", "dense_component_cell_fill_x"),
        support_dense_component_cell_fill_y = graph_layout_value("support_threshold_consensus_graph", "dense_component_cell_fill_y"),
        support_c2_component_cell_fill_x = graph_layout_value("support_threshold_consensus_graph", "c2_component_cell_fill_x"),
        support_c2_component_cell_fill_y = graph_layout_value("support_threshold_consensus_graph", "c2_component_cell_fill_y"),
        c1_x_expand_factor = graph_layout_value("support_threshold_consensus_graph", "c1_x_expand_factor"),
        c1_y_expand_factor = graph_layout_value("support_threshold_consensus_graph", "c1_y_expand_factor"),
        c2_x_expand_factor = graph_layout_value("support_threshold_consensus_graph", "c2_x_expand_factor"),
        c2_y_expand_factor = graph_layout_value("support_threshold_consensus_graph", "c2_y_expand_factor")
    log: os.path.join(LOGROOT, "plot_patient_referenced_support_threshold_consensus_cell_line_similarity_network.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {input.script} \
          --edges {input.edges} \
          --nodes {input.shortnames} \
          --nodes-col short_id \
          --out {params.out_prefix} \
          --label {params.label} \
          --graph-mode majority_threshold \
          --graph-title "{params.graph_title}" \
          --fig-w {params.fig_w} \
          --fig-h {params.fig_h} \
          --dpi {params.dpi} \
          {params.use_halo_labels_arg} \
          --halo-linewidth {params.halo_linewidth} \
          --halo-color {params.halo_color} \
          {params.match_resolved_style_arg} \
          --position-fill-x {params.position_fill_x} \
          --position-fill-y {params.position_fill_y} \
          --axis-margin-frac {params.axis_margin_frac} \
          --node-size {params.node_size} \
          --font-size {params.font_size} \
          --edge-width-min {params.edge_width_min} \
          --edge-width-max {params.edge_width_max} \
          --component-label-gap {params.component_label_gap} \
          --nbl-component-gap-fraction {params.nbl_component_gap_fraction} \
          --nbl-component-min-gap {params.nbl_component_min_gap} \
          --legend-anchor-x {params.legend_anchor_x} \
          {params.legend_anchor_y_arg} \
          --legend-align-mode {params.legend_align_mode} \
          --legend-gap-frac {params.legend_gap_frac} \
          --legend-mode {params.legend_mode} \
          --legend-out-prefix {params.legend_out_prefix} \
          {params.compact_isolate_panel_arg} \
          --isolate-panel-width-mode {params.isolate_panel_width_mode} \
          --isolate-panel-spacing {params.isolate_panel_spacing} \
          --min-isolate-label-spacing {params.min_isolate_label_spacing} \
          --isolate-label-y-frac {params.isolate_label_y_frac} \
          --isolate-box-height-factor {params.isolate_box_height_factor} \
          --isolate-gap {params.isolate_gap} \
          --isolate-label-left-pad {params.isolate_label_left_pad} \
          --isolate-spacing {params.isolate_spacing} \
          --isolate-panel-padding-x {params.isolate_panel_padding_x} \
          --isolate-panel-padding-y {params.isolate_panel_padding_y} \
          --isolate-label-band-frac {params.isolate_label_band_frac} \
          --isolate-max-per-row {params.isolate_max_per_row} \
          --isolate-panel-width-frac {params.isolate_panel_width_frac} \
          --dense-component-min-size {params.dense_component_min_size} \
          --dense-component-min-density {params.dense_component_min_density} \
          --dense-component-expand-factor {params.dense_component_expand_factor} \
          {params.adaptive_component_layout_arg} \
          --adaptive-component-min-nodes {params.adaptive_component_min_nodes} \
          --adaptive-component-min-density {params.adaptive_component_min_density} \
          --adaptive-component-min-expand {params.adaptive_component_min_expand} \
          --adaptive-component-max-expand {params.adaptive_component_max_expand} \
          --adaptive-component-density-weight {params.adaptive_component_density_weight} \
          --adaptive-component-size-weight {params.adaptive_component_size_weight} \
          --adaptive-component-label-weight {params.adaptive_component_label_weight} \
          {params.label_overlap_avoidance_arg} \
          --label-padding-factor {params.label_padding_factor} \
          --label-bbox-margin {params.label_bbox_margin} \
          --max-label-overlap-iterations {params.max_label_overlap_iterations} \
          --min-label-separation {params.min_label_separation} \
          {params.avoid_isolate_overlap_arg} \
          --min-component-isolate-gap {params.min_component_isolate_gap} \
          --isolate-panel-reserved-height {params.isolate_panel_reserved_height} \
          {params.adaptive_isolate_y_shift_arg} \
          --component-isolate-gap-factor {params.component_isolate_gap_factor} \
          --component-grid-ncols {params.support_grid_cols} \
          --component-grid-nrows {params.support_grid_rows} \
          --component-cell-width {params.support_cell_width} \
          --component-cell-height {params.support_cell_height} \
          --component-row-gap {params.support_row_gap} \
          --component-col-gap {params.support_col_gap} \
          --legend-width {params.support_legend_width} \
          --isolate-region-height {params.support_isolate_height} \
          --footnote-region-height {params.support_footnote_height} \
          --component-cell-fill-x {params.support_component_cell_fill_x} \
          --component-cell-fill-y {params.support_component_cell_fill_y} \
          --dense-component-cell-fill-x {params.support_dense_component_cell_fill_x} \
          --dense-component-cell-fill-y {params.support_dense_component_cell_fill_y} \
          --c2-component-cell-fill-x {params.support_c2_component_cell_fill_x} \
          --c2-component-cell-fill-y {params.support_c2_component_cell_fill_y} \
          --c1-x-expand-factor {params.c1_x_expand_factor} \
          --c1-y-expand-factor {params.c1_y_expand_factor} \
          --c2-x-expand-factor {params.c2_x_expand_factor} \
          --c2-y-expand-factor {params.c2_y_expand_factor} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.png} || (echo "ERROR: missing {output.png}" >&2; exit 1)
        test -s {output.svg} || (echo "ERROR: missing {output.svg}" >&2; exit 1)
        test -s {output.legend_pdf} || (echo "ERROR: missing {output.legend_pdf}" >&2; exit 1)
        test -s {output.legend_png} || (echo "ERROR: missing {output.legend_png}" >&2; exit 1)
        test -s {output.legend_svg} || (echo "ERROR: missing {output.legend_svg}" >&2; exit 1)
        '''

if IS_MULTICOHORT_PROFILE:

    # Rule: plot_pan_cancer_resolved_graph_inspection
    # Method role: supplementary plotting rule for dense multicohort resolved-graph inspection.
    # Flow: resolved graph sidecars -> component panels, interactive HTML, and node-label tables.
    # Used by: audit and visual inspection; it does not alter resolved neighbours.
    rule plot_pan_cancer_resolved_graph_inspection:
        """
        Component-panel and interactive inspection outputs for the dense
        multicohort resolved graph. These are generated from the audited graph
        sidecars and do not alter graph topology or anchor selection.
        """
        input:
            edges = P_CONS_RESOLVED_EDGES_TSV,
            node_stats = P_CONS_RESOLVED_NODE_STATS_TSV,
            display_names = P_CONS_SHORTNAMES_TSV,
            metadata = P_CONS_DSMZ_META_CSV,
            script = os.path.join(BASE, "scripts", "plot_pan_cancer_graph_inspection.py")
        output:
            panels_pdf = P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".pdf",
            panels_png = P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".png",
            panels_svg = P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX + ".svg",
            interactive_html = P_CONS_RESOLVED_INTERACTIVE_HTML,
            component_summary = P_CONS_RESOLVED_COMPONENT_SUMMARY_TSV,
            node_labels = P_CONS_RESOLVED_NODE_LABELS_TSV
        params:
            title_label = "Resolved multicohort cell-line neighbourhood graph"
        log: os.path.join(LOGROOT, "plot_pan_cancer_resolved_graph_inspection.log")
        conda: CONDA_ENV_PY
        shell:
            r'''
            mkdir -p {P_CONS_PLOTS_DIR}
            python {input.script} \
              --graph-type resolved \
              --edges {input.edges} \
              --node-stats {input.node_stats} \
              --display-names {input.display_names} \
              --metadata {input.metadata} \
              --panels-pdf {output.panels_pdf} \
              --panels-png {output.panels_png} \
              --panels-svg {output.panels_svg} \
              --interactive-html {output.interactive_html} \
              --component-summary-out {output.component_summary} \
              --node-labels-out {output.node_labels} \
              --label "{params.title_label}" \
              > {log} 2>&1
            test -s {output.panels_pdf} || (echo "ERROR: missing {output.panels_pdf}" >&2; exit 1)
            test -s {output.panels_png} || (echo "ERROR: missing {output.panels_png}" >&2; exit 1)
            test -s {output.panels_svg} || (echo "ERROR: missing {output.panels_svg}" >&2; exit 1)
            test -s {output.interactive_html} || (echo "ERROR: missing {output.interactive_html}" >&2; exit 1)
            test -s {output.component_summary} || (echo "ERROR: missing {output.component_summary}" >&2; exit 1)
            test -s {output.node_labels} || (echo "ERROR: missing {output.node_labels}" >&2; exit 1)
            '''


    # Rule: plot_pan_cancer_support_threshold_graph_inspection
    # Method role: supplementary plotting rule for dense multicohort support-threshold graph inspection.
    # Flow: support-threshold graph sidecars -> component panels, interactive HTML, and node-label tables.
    # Used by: audit and visual inspection; it does not alter support-network edges.
    rule plot_pan_cancer_support_threshold_graph_inspection:
        """
        Component-panel and interactive inspection outputs for the dense
        multicohort support-threshold consensus graph. These are generated from
        the audited consensus edge sidecar and short-name table.
        """
        input:
            edges = P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV,
            display_names = P_CONS_SHORTNAMES_TSV,
            metadata = P_CONS_DSMZ_META_CSV,
            script = os.path.join(BASE, "scripts", "plot_pan_cancer_graph_inspection.py")
        output:
            panels_pdf = P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".pdf",
            panels_png = P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".png",
            panels_svg = P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".svg",
            interactive_html = P_CONS_SUPPORT_INTERACTIVE_HTML,
            component_summary = P_CONS_SUPPORT_COMPONENT_SUMMARY_TSV,
            node_labels = P_CONS_SUPPORT_NODE_LABELS_TSV
        params:
            title_label = "Multicohort support-threshold consensus network"
        log: os.path.join(LOGROOT, "plot_pan_cancer_support_threshold_graph_inspection.log")
        conda: CONDA_ENV_PY
        shell:
            r'''
            mkdir -p {P_CONS_PLOTS_DIR}
            python {input.script} \
              --graph-type support_threshold \
              --edges {input.edges} \
              --display-names {input.display_names} \
              --metadata {input.metadata} \
              --panels-pdf {output.panels_pdf} \
              --panels-png {output.panels_png} \
              --panels-svg {output.panels_svg} \
              --interactive-html {output.interactive_html} \
              --component-summary-out {output.component_summary} \
              --node-labels-out {output.node_labels} \
              --label "{params.title_label}" \
              > {log} 2>&1
            test -s {output.panels_pdf} || (echo "ERROR: missing {output.panels_pdf}" >&2; exit 1)
            test -s {output.panels_png} || (echo "ERROR: missing {output.panels_png}" >&2; exit 1)
            test -s {output.panels_svg} || (echo "ERROR: missing {output.panels_svg}" >&2; exit 1)
            test -s {output.interactive_html} || (echo "ERROR: missing {output.interactive_html}" >&2; exit 1)
            test -s {output.component_summary} || (echo "ERROR: missing {output.component_summary}" >&2; exit 1)
            test -s {output.node_labels} || (echo "ERROR: missing {output.node_labels}" >&2; exit 1)
            '''


# =============================================================================
# STUDY-DESIGN MATERIALISATION
# =============================================================================
# Stage role: exports configured study-design metadata used by validation and reporting.

# Rule: materialize_study_design
# Method role: metadata rule that exports the configured study design for reporting and validation.
# Flow: project config and study-design YAML -> cohort, label, inference, and endpoint manifests.
# Documents the design assumptions used by dependent analyses.
rule materialize_study_design:
    input:
        cfg = CFGFILE_ABS,
        study_design = STUDY_DESIGN_FILE
    output:
        question_txt = os.path.join(STUDY_OUTPUT_DIR_REL, "study_question.txt"),
        cohort_manifest = os.path.join(STUDY_OUTPUT_DIR_REL, "cohort_manifest.tsv"),
        labels_manifest = os.path.join(STUDY_OUTPUT_DIR_REL, "cohort_labels.tsv"),
        inference_manifest = os.path.join(STUDY_OUTPUT_DIR_REL, "candidate_inference.tsv"),
        endpoint_manifest = os.path.join(STUDY_OUTPUT_DIR_REL, "endpoint_manifest.tsv")
    params:
        script = os.path.join(BASE, "scripts", "materialize_study_design.R"),
        config = CFGFILE_ABS,
        study_design = STUDY_DESIGN_FILE,
        # The strong-support threshold reported in the candidate-inference
        # manifest is the same configured p-consensus fraction threshold the
        # graph and reporting stages apply. It is transmitted explicitly so the
        # study-design manifest cannot become a second owner of the value.
        p_consensus_threshold = PATIENT_REFERENCED_P_CONSENSUS_THRESHOLD,
        # The representation grid documented in the manifest is the profile's
        # configured representation set, enumerated here rather than rebuilt
        # inside the script.
        directions = ",".join(FEATURE_DISTANCE_DIRECTIONS)
    log: os.path.join(LOGROOT, "materialize_study_design.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.question_txt})
        Rscript {params.script}           --config {params.config}           --study-design {params.study_design}           --profile "{profile_name}"           --p-consensus-threshold {params.p_consensus_threshold}           --directions {params.directions}           --out-question {output.question_txt}           --out-cohorts {output.cohort_manifest}           --out-labels {output.labels_manifest}           --out-inference {output.inference_manifest}           --out-endpoints {output.endpoint_manifest}           > {log} 2>&1
        '''

# =============================================================================
# VALIDATION ANALYSES
# =============================================================================
# Stage role: runs model-selection, permutation, baseline, and separation metrics without redefining targets.

# Rule: model_selection_summary
# Method role: validation rule evaluating selected representations and model-selection evidence.
# Flow: p-consensus fraction tables and graph annotations -> model-selection table, plot, and notes.
# Used by: reporting; it does not choose new workflow targets.
rule model_selection_summary:
    input:
        cfg = CFGFILE_ABS,
        ranked_best = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_ranked.tsv"),
        long_tbl = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.long.tsv"),
        winning_direction = os.path.join(P_CONS_ALL_DIR, "winning_direction.txt"),
        graph_nodes = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_node_annotations_{direction}.tsv"),
            direction=FEATURE_DISTANCE_DIRECTIONS
        )
    output:
        summary_tsv = os.path.join(VALIDATION_OUTPUT_DIR_REL, "model_selection_summary.tsv"),
        plot_pdf = os.path.join(VALIDATION_OUTPUT_DIR_REL, "model_selection_plot.pdf"),
        notes_txt = os.path.join(VALIDATION_OUTPUT_DIR_REL, "model_selection_notes.txt")
    params:
        script = os.path.join(BASE, "validation", "04_model_selection_summary.R"),
        config = CFGFILE_ABS
    log: os.path.join(LOGROOT, "model_selection_summary.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.summary_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --out-tsv {output.summary_tsv}           --out-plot {output.plot_pdf}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''


# Rule: neighbourhood_permutation_validation
# Method role: validation rule evaluating observed neighbourhood metrics against permutation expectations.
# Flow: p-consensus ranking tables -> permutation metrics, plot, and notes.
# Used by: sensitivity evidence for neighbourhood signal, not a replacement analysis.
rule neighbourhood_permutation_validation:
    input:
        cfg = CFGFILE_ABS,
        ranked_best = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_ranked.tsv"),
        long_tbl = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.long.tsv")
    output:
        summary_tsv = os.path.join(VALIDATION_OUTPUT_DIR_REL, "neighbourhood_permutation_summary.tsv"),
        plot_pdf = os.path.join(VALIDATION_OUTPUT_DIR_REL, "neighbourhood_permutation_plot.pdf"),
        notes_txt = os.path.join(VALIDATION_OUTPUT_DIR_REL, "neighbourhood_permutation_notes.txt")
    params:
        script = os.path.join(BASE, "validation", "01_permutation_test_neighbourhood.R"),
        config = CFGFILE_ABS,
        seed = PIPELINE_SEED
    log: os.path.join(LOGROOT, "neighbourhood_permutation_validation.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.summary_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --seed {params.seed}           --out-tsv {output.summary_tsv}           --out-plot {output.plot_pdf}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''


# Rule: random_baseline_comparison
# Method role: validation rule evaluating resolved rankings against a random-baseline distribution.
# Flow: ranked p-consensus cell-line table -> baseline metrics, plot, and notes.
# Purpose: contextualises ranking behaviour without changing graph outputs.
rule random_baseline_comparison:
    input:
        cfg = CFGFILE_ABS,
        ranked_best = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_ranked.tsv")
    output:
        summary_tsv = os.path.join(VALIDATION_OUTPUT_DIR_REL, "random_baseline_summary.tsv"),
        plot_pdf = os.path.join(VALIDATION_OUTPUT_DIR_REL, "random_baseline_plot.pdf"),
        notes_txt = os.path.join(VALIDATION_OUTPUT_DIR_REL, "random_baseline_notes.txt")
    params:
        script = os.path.join(BASE, "validation", "02_random_baseline_comparison.R"),
        config = CFGFILE_ABS,
        seed = PIPELINE_SEED
    log: os.path.join(LOGROOT, "random_baseline_comparison.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.summary_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --seed {params.seed}           --out-tsv {output.summary_tsv}           --out-plot {output.plot_pdf}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''


# Rule: silhouette_report
# Method role: validation rule reporting separation metrics for configured representations.
# Flow: profile config -> silhouette/separation report and notes.
# Purpose: reporting only.
rule silhouette_report:
    input:
        cfg = CFGFILE_ABS
    output:
        report_tsv = os.path.join(VALIDATION_OUTPUT_DIR_REL, "silhouette_report.tsv"),
        notes_txt = os.path.join(VALIDATION_OUTPUT_DIR_REL, "silhouette_notes.txt")
    params:
        script = os.path.join(BASE, "validation", "03_silhouette_report.R"),
        config = CFGFILE_ABS
    log: os.path.join(LOGROOT, "silhouette_report.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.report_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --out-tsv {output.report_tsv}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''



# Rule: community_stability_analysis
# Method role: sensitivity analysis evaluating community assignments across feature-distance representations.
# Flow: resolved graph, display names, and direction-specific annotations -> stability tables and figures.
# Used by: supplementary evidence on community robustness, not final community selection.
rule community_stability_analysis:
    """
    Cross-direction cell-line graph community stability assessment.

    Leiden/Louvain assignments are direction-specific assessments. This rule
    compares pairwise same-community co-assignment frequencies across directions
    against the final resolved graph components. It does not create or interpret
    a final consensus Leiden/Louvain community plot.
    """
    input:
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv"),
        shortnames = P_CONS_SHORTNAMES_TSV,
        node_annots = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_node_annotations_{direction}.tsv"),
            direction=CONS_DIRECTIONS
        )
    output:
        assignments = COMMUNITY_STABILITY_ASSIGNMENTS,
        matrix_leiden = os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_matrix_Leiden.tsv"),
        matrix_louvain = os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_matrix_Louvain.tsv"),
        nvalid_leiden = os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_n_valid_Leiden.tsv"),
        nvalid_louvain = os.path.join(COMMUNITY_STABILITY_DIR, "community_coassignment_n_valid_Louvain.tsv"),
        component_summary = COMMUNITY_STABILITY_COMPONENT_SUMMARY,
        global_summary = COMMUNITY_STABILITY_GLOBAL_SUMMARY,
        heatmap_leiden = os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_coassignment_heatmap_Leiden.pdf"),
        heatmap_louvain = os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_coassignment_heatmap_Louvain.pdf"),
        within_between_leiden = os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_stability_within_between_Leiden.pdf"),
        within_between_louvain = os.path.join(COMMUNITY_STABILITY_DIR, "Fig_community_stability_within_between_Louvain.pdf"),
        validation = COMMUNITY_STABILITY_VALIDATION
    params:
        script = os.path.join(BASE, "scripts", "analyze_community_stability.py"),
        tumour_nh_root = TUMOUR_NH_ROOT_ABS,
        outdir = COMMUNITY_STABILITY_DIR,
        directions = ",".join(CONS_DIRECTIONS),
        cohort = profile_name
    log: os.path.join(LOGROOT, "community_stability_analysis.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {params.outdir}
        python {params.script} \
          --cohort "{params.cohort}" \
          --tumour-nh-root {params.tumour_nh_root} \
          --directions "{params.directions}" \
          --resolved-tsv {input.resolved_tsv} \
          --shortnames-tsv {input.shortnames} \
          --outdir {params.outdir} \
          > {log} 2>&1
        test -s {output.assignments} || (echo "ERROR: missing {output.assignments}" >&2; exit 1)
        test -s {output.global_summary} || (echo "ERROR: missing {output.global_summary}" >&2; exit 1)
        '''


# Rule: post_resolution_edge_support_stratification
# Method role: supplementary analysis stratifying resolved graph edges by representation support.
# Flow: resolved/support network sidecars for one cohort -> edge-support stratification table and plots.
# Used by: cohort-level audit of post-resolution support.
rule post_resolution_edge_support_stratification:
    """
    Audit of edge support after resolved-neighbour graph construction.

    This rule does not modify the union-supported, support-threshold consensus,
    or resolved-neighbour graphs. It reads their existing edge/node sidecars and
    classifies each focal cohort cell line by the highest post-resolution edge
    evidence available.
    """
    input:
        union_edges = _post_resolution_union_edges,
        consensus_edges = _post_resolution_consensus_edges,
        resolved_edges = _post_resolution_resolved_edges,
        nodes = _post_resolution_nodes
    output:
        table = os.path.join(
            POST_RESOLUTION_TABLE_DIR,
            "{cohort}_post_resolution_edge_support_stratification.tsv",
        )
    wildcard_constraints:
        cohort = "NBL|BRCA|RBL"
    params:
        script = os.path.join(BASE, "scripts", "post_resolution_edge_support_stratification.py")
    log:
        os.path.join(LOGROOT, "post_resolution_edge_support_stratification", "{cohort}.log")
    conda:
        CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {POST_RESOLUTION_TABLE_DIR} $(dirname {log})
        python {params.script} classify \
          --cohort {wildcards.cohort} \
          --union-edges {input.union_edges} \
          --consensus-edges {input.consensus_edges} \
          --resolved-edges {input.resolved_edges} \
          --nodes {input.nodes} \
          --out {output.table} \
          > {log} 2>&1
        test -s {output.table} || (echo "ERROR: missing {output.table}" >&2; exit 1)
        '''


# Rule: combine_post_resolution_edge_support_stratification
# Method role: collector target combining post-resolution support stratifications across cohorts.
# Flow: cohort stratification tables -> combined support-stratification table.
# Used by: cross-cohort supplementary table for reporting.
rule combine_post_resolution_edge_support_stratification:
    input:
        tables = POST_RESOLUTION_STRATIFICATION_TABLES
    output:
        table = POST_RESOLUTION_COMBINED_STRATIFICATION_TSV
    params:
        script = os.path.join(BASE, "scripts", "post_resolution_edge_support_stratification.py")
    log:
        os.path.join(LOGROOT, "post_resolution_edge_support_stratification", "combined.log")
    conda:
        CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {POST_RESOLUTION_TABLE_DIR} $(dirname {log})
        python {params.script} combine \
          --out {output.table} \
          {input.tables} \
          > {log} 2>&1
        test -s {output.table} || (echo "ERROR: missing {output.table}" >&2; exit 1)
        '''


# =============================================================================
# AUDIT AND SUPPLEMENTARY GRAPH ANALYSES
# =============================================================================
# Stage role: computes graph audit tables, stratifications, and supplementary community reports.

# Rule: compute_multicohort_cancer_communities
# Method role: supplementary graph analysis assigning communities in the multicohort resolved graph.
# Flow: resolved graph sidecars -> community assignments, modularity, layout, and figure outputs.
# Used by: audit/supplementary community report for the multicohort graph.
rule compute_multicohort_cancer_communities:
    """
    Unweighted Leiden community detection on the final resolved MULTICOHORT_CANCER
    DSMZ cell-line graph (56 nodes, 155 edges, 9 components).
    Produces per-node community assignments, community summary, modularity table,
    layout coordinates, and a two-panel cancer-type/community PDF+PNG figure.
    Uses igraph::cluster_leiden() — does NOT fall back to Louvain.
    Runs only for the multicohort_cancer profile.
    """
    input:
        edges         = P_CONS_RESOLVED_EDGES_TSV,
        node_stats    = P_CONS_RESOLVED_NODE_STATS_TSV,
        anchor_audit  = P_CONS_ANCHOR_AUDIT_TSV,
        shortnames    = P_CONS_SHORTNAMES_TSV,
        metadata      = MC_META,
    output:
        communities       = P_CONS_LEIDEN_COMMUNITIES_TSV,
        community_summary = P_CONS_LEIDEN_SUMMARY_TSV,
        modularity        = P_CONS_LEIDEN_MODULARITY_TSV,
        layout            = P_CONS_LEIDEN_LAYOUT_TSV,
        pdf               = P_CONS_LEIDEN_FIG_PREFIX + ".pdf",
        png               = P_CONS_LEIDEN_FIG_PREFIX + ".png",
    params:
        script        = os.path.join(BASE, "scripts",
                            "compute_and_plot_multicohort_cancer_communities.R"),
        seed          = PIPELINE_SEED,
    log: os.path.join(LOGROOT, "compute_multicohort_cancer_communities.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.communities})
        mkdir -p $(dirname {log})
        Rscript {params.script} \
          --edges                   {input.edges} \
          --node-stats              {input.node_stats} \
          --anchor-audit            {input.anchor_audit} \
          --shortnames              {input.shortnames} \
          --metadata                {input.metadata} \
          --out-communities         {output.communities} \
          --out-community-summary   {output.community_summary} \
          --out-modularity          {output.modularity} \
          --out-layout              {output.layout} \
          --out-pdf                 {output.pdf} \
          --out-png                 {output.png} \
          --seed                    {params.seed} \
          > {log} 2>&1
        test -s {output.communities}       || (echo "ERROR: missing {output.communities}" >&2; exit 1)
        test -s {output.community_summary} || (echo "ERROR: missing {output.community_summary}" >&2; exit 1)
        test -s {output.modularity}        || (echo "ERROR: missing {output.modularity}" >&2; exit 1)
        test -s {output.layout}            || (echo "ERROR: missing {output.layout}" >&2; exit 1)
        test -s {output.pdf}               || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        test -s {output.png}               || (echo "ERROR: missing {output.png}" >&2; exit 1)
        '''

# =============================================================================
# DESEQ2 MARKER ANALYSIS
# =============================================================================
# Stage role: derives graph-informed marker contrasts from resolved isolates and component anchors.
# Gated to brca, nbl, and rbl only — the profiles whose disease-specific
# raw-count source tables are declared under defaults.deseq2_inputs and whose
# metadata schema is compatible with prepare_deseq2_inputs.R. heme and
# pan_cancer are deferred.

DESEQ2_PROFILE_CFG = cfg.get("deseq2", {})
DESEQ2_CFG = deep_merge(config.get("defaults", {}).get("deseq2", {}), DESEQ2_PROFILE_CFG)
DESEQ2_MARKER_PRIORITISATION_PROFILE_CFG = cfg.get("deseq2_marker_prioritisation", {})
DESEQ2_MARKER_PRIORITISATION_CFG = deep_merge(
    config.get("defaults", {}).get("deseq2_marker_prioritisation", {}),
    DESEQ2_MARKER_PRIORITISATION_PROFILE_CFG,
)
DESEQ2_ISOLATE_CONTRAST_CFG = DESEQ2_MARKER_PRIORITISATION_CFG.get("isolate_contrast", {})
DESEQ2_ANCHOR_CONTRAST_CFG = DESEQ2_MARKER_PRIORITISATION_CFG.get("anchor_contrast", {})
DESEQ2_COMPONENT_CONTRAST_CFG = DESEQ2_MARKER_PRIORITISATION_CFG.get("component_contrast", {})
DESEQ2_RETAINED_MARKER_RECURRENCE_THRESHOLD = 2
DESEQ2_INPUT_DIR     = os.path.join(UNSUP_REL, "deseq2_inputs")
DESEQ2_INPUT_DIR_ABS = os.path.join(UNSUP,     "deseq2_inputs")
DESEQ2_MARKER_OUTDIR_NAME = DESEQ2_CFG.get("marker_outdir_name", "deseq2_markers")
DESEQ2_ISOLATE_DIR     = os.path.join(UNSUP_REL, DESEQ2_MARKER_OUTDIR_NAME)
DESEQ2_ISOLATE_DIR_ABS = os.path.join(UNSUP,     DESEQ2_MARKER_OUTDIR_NAME)
DESEQ2_COMP_DIR     = os.path.join(UNSUP_REL, "deseq2", "component_vs_rest")
DESEQ2_COMP_DIR_ABS = os.path.join(UNSUP,     "deseq2", "component_vs_rest")
DESEQ2_INPUTS_PROFILE_CFG = cfg.get("deseq2_inputs", {})
DESEQ2_INPUTS_CFG = deep_merge(config.get("defaults", {}).get("deseq2_inputs", {}), DESEQ2_INPUTS_PROFILE_CFG)
# The profiles that run DESeq2 marker analysis are a configured scientific
# choice (defaults.deseq2.enabled_profiles); there is no built-in profile list
# to fall back to.
_enabled_profiles = DESEQ2_CFG.get("enabled_profiles")
if not _enabled_profiles:
    raise ValueError("deseq2.enabled_profiles is required and must be a non-empty list")
DESEQ2_ENABLED_PROFILES = {str(p) for p in _enabled_profiles}
_enabled_override = DESEQ2_PROFILE_CFG.get("enabled")
if _enabled_override is None:
    DESEQ2_ENABLED = profile_name in DESEQ2_ENABLED_PROFILES
else:
    DESEQ2_ENABLED = bool(_enabled_override)

def require_config_value(mapping, key, config_path):
    """Fail during workflow initialisation when an active canonical parameter is absent."""
    value = mapping.get(key)
    if value is None or value == "":
        raise ValueError(f"{config_path}.{key} is required")
    return value


if DESEQ2_ENABLED:
    for _key in ("minimum_total_gene_count", "minimum_base_mean", "dispersion_fit_type", "lfc_shrinkage_method"):
        require_config_value(DESEQ2_MARKER_PRIORITISATION_CFG, _key, "deseq2_marker_prioritisation")
    for _key in ("sample_id_column", "cell_line_column", "component_column"):
        require_config_value(DESEQ2_MARKER_PRIORITISATION_CFG, _key, "deseq2_marker_prioritisation")
    for _key in ("adjusted_p_value_threshold", "minimum_absolute_shrunken_log2fc", "maximum_markers_per_contrast"):
        require_config_value(DESEQ2_ISOLATE_CONTRAST_CFG, _key, "deseq2_marker_prioritisation.isolate_contrast")
        require_config_value(DESEQ2_ANCHOR_CONTRAST_CFG, _key, "deseq2_marker_prioritisation.anchor_contrast")
    for _key in ("adjusted_p_value_threshold", "minimum_absolute_shrunken_log2fc", "maximum_markers_per_direction"):
        require_config_value(DESEQ2_COMPONENT_CONTRAST_CFG, _key, "deseq2_marker_prioritisation.component_contrast")
    # Staged-input declarations are config-owned in the same way, so the rules
    # below can transmit them without restating a default.
    for _key in ("source_count_kind", "sample_population"):
        require_config_value(DESEQ2_INPUTS_CFG, _key, "deseq2_inputs")


def deseq2_profile_input_path(key):
    """Resolve a per-cohort source input path for canonical DESeq2 preparation."""
    by_profile = DESEQ2_INPUTS_CFG.get(key, {})
    if by_profile is None:
        by_profile = {}
    if not isinstance(by_profile, dict):
        raise ValueError(f"deseq2_inputs.{key} must map profile names to paths")
    rel = by_profile.get(profile_name)
    if DESEQ2_ENABLED and not rel:
        raise ValueError(
            f"deseq2_inputs.{key}.{profile_name} is required for canonical DESeq2 preparation"
        )
    return rel


def deseq2_expected_samples():
    """Return the configured prepared sample count for the active profile, if set."""
    by_profile = DESEQ2_INPUTS_CFG.get("expected_samples", {})
    if not isinstance(by_profile, dict):
        return ""
    value = by_profile.get(profile_name, "")
    return "" if value is None else str(value)


def looks_like_transformed_count_path(path):
    """Catch common transformed-expression path tokens before DESeq2 consumes them."""
    name = os.path.basename(str(path)).lower()
    tokens = ("vst", "tpm", "fpkm", "normalised", "normalized", "scaled", "batch", "corrected")
    return any(re.search(rf"(^|[^a-z0-9]){token}([^a-z0-9]|$)", name) for token in tokens)


DESEQ2_PREPARED_INPUT_DIR = os.path.join(DESEQ2_INPUT_DIR, "prepared")
DESEQ2_PREPARED_INPUT_DIR_ABS = os.path.join(DESEQ2_INPUT_DIR_ABS, "prepared")
DESEQ2_SOURCE_COUNTS = deseq2_profile_input_path("source_counts")
DESEQ2_SOURCE_METADATA = deseq2_profile_input_path("source_metadata")
DESEQ2_COUNTS_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "counts.tsv")
DESEQ2_METADATA_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "metadata.tsv")
DESEQ2_META_COMP_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "metadata_with_components.tsv")
DESEQ2_SAMPLE_MAPPING_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "sample_mapping.tsv")
DESEQ2_INPUT_PROVENANCE_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "input_provenance.tsv")
DESEQ2_ISOLATE_CSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "isolate_list.csv")
DESEQ2_ANCHOR_LIST_CSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "anchor_list.csv")
DESEQ2_ANCHOR_COMPONENTS_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "anchor_components.tsv")
DESEQ2_COMPONENTS_LIST_TXT = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "components_list.txt")
DESEQ2_INPUT_VALIDATION_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "deseq2_input_validation.tsv")
DESEQ2_INPUT_METRICS_TSV = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "deseq2_input_metrics.tsv")
DESEQ2_INPUT_VALIDATION_REQUIRED_TSV = DESEQ2_INPUT_VALIDATION_TSV
DESEQ2_EXPECTED_SAMPLE_COUNT = deseq2_expected_samples()

DESEQ2_NODE_STATS_TSV = P_CONS_RESOLVED_NODE_STATS_TSV

if DESEQ2_ENABLED:

    DESEQ2_SOURCE_COUNTS_ABS = abspath(DESEQ2_SOURCE_COUNTS)
    DESEQ2_SOURCE_METADATA_ABS = abspath(DESEQ2_SOURCE_METADATA)
    for _key in ("source_sample_id_column", "source_cell_line_column", "source_count_kind", "sample_population"):
        require_config_value(DESEQ2_INPUTS_CFG, _key, "deseq2_inputs")
    if looks_like_transformed_count_path(DESEQ2_SOURCE_COUNTS):
        raise ValueError(
            "Unsafe DESeq2 count source: deseq2_inputs.source_counts appears to reference "
            f"a transformed expression matrix ({DESEQ2_SOURCE_COUNTS}). Configure "
            "deseq2_inputs.source_counts for this profile to raw gene-level counts."
        )
    DESEQ2_CELL_LINE_VERIFICATION_CFG = DESEQ2_INPUTS_CFG.get("cell_line_verification", {})

    # Rule: validate_deseq2_inputs
    # Method role: validation rule checking canonical prepared DESeq2 input tables before marker analysis.
    # Flow: prepared count, metadata, graph contrast, and provenance files -> validation and metric tables.
    # Used by: guards DESeq2 marker rules against malformed prepared inputs.
    rule validate_deseq2_inputs:
        """
        Checks prepared count, metadata, isolate, anchor, and component
        sidecars before marker DESeq2 rules are allowed to run.
        """
        input:
            counts_tsv            = DESEQ2_COUNTS_TSV,
            meta_comp             = DESEQ2_META_COMP_TSV,
            isolate_csv           = DESEQ2_ISOLATE_CSV,
            anchor_list_csv       = DESEQ2_ANCHOR_LIST_CSV,
            anchor_components_tsv = DESEQ2_ANCHOR_COMPONENTS_TSV,
            components_list       = DESEQ2_COMPONENTS_LIST_TXT,
            node_stats            = DESEQ2_NODE_STATS_TSV,
            anchor_audit          = P_CONS_ANCHOR_AUDIT_TSV,
            provenance            = DESEQ2_INPUT_PROVENANCE_TSV
        output:
            validation_tsv = DESEQ2_INPUT_VALIDATION_TSV,
            metrics_tsv = DESEQ2_INPUT_METRICS_TSV
        params:
            script = os.path.join(BASE, "scripts", "validate_deseq2_staged_inputs.R"),
            cohort = profile_name,
            sample_id_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("sample_id_column"),
            cell_line_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("cell_line_column"),
            component_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("component_column"),
            expected_samples = DESEQ2_EXPECTED_SAMPLE_COUNT,
            verification_mode = DESEQ2_CELL_LINE_VERIFICATION_CFG.get("mode", "source_population_declaration"),
            verification_column = DESEQ2_CELL_LINE_VERIFICATION_CFG.get("column", ""),
            accepted_values = ",".join(DESEQ2_CELL_LINE_VERIFICATION_CFG.get("accepted_values", ["Cell Line"]))
        log: os.path.join(LOGROOT, "validate_deseq2_inputs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript {params.script} \
              --cohort "{params.cohort}" \
              --counts {input.counts_tsv} \
              --metadata {input.meta_comp} \
              --isolate_list {input.isolate_csv} \
              --anchor_list {input.anchor_list_csv} \
              --anchor_components {input.anchor_components_tsv} \
              --components_list {input.components_list} \
              --node_stats {input.node_stats} \
              --anchor_audit {input.anchor_audit} \
              --provenance {input.provenance} \
              --sample_id_col {params.sample_id_col} \
              --cell_line_col {params.cell_line_col} \
              --component_col {params.component_col} \
              --expected_samples "{params.expected_samples}" \
              --cell_line_verification_mode "{params.verification_mode}" \
              --cell_line_verification_column "{params.verification_column}" \
              --accepted_cell_line_values "{params.accepted_values}" \
              --output {output.validation_tsv} \
              --metrics_output {output.metrics_tsv} \
              > {log} 2>&1
            test -s {output.validation_tsv} || (echo "ERROR: missing {output.validation_tsv}" >&2; exit 1)
            test -s {output.metrics_tsv} || (echo "ERROR: missing {output.metrics_tsv}" >&2; exit 1)
            if awk -F'\t' 'NR>1 && $2=="FAIL" {{ bad=1 }} END {{ exit bad ? 0 : 1 }}' {output.validation_tsv}; then
                echo "ERROR: DESeq2 input validation contains FAIL rows: {output.validation_tsv}" >&2
                awk -F'\t' 'NR==1 || $2=="FAIL"' {output.validation_tsv} >&2
                exit 1
            fi
            '''

    # Rule: prepare_deseq2_inputs
    # Method role: data-preparation rule that builds DESeq2 count and metadata tables from resolved graph labels.
    # Flow: source counts/metadata and resolved-graph node statistics -> prepared DESeq2 input tables.
    # Provides input for isolate, component, and anchor marker contrasts.
    rule prepare_deseq2_inputs:
        """
        Extracts cancer-type-specific raw count profiles and metadata from the
        configured source tables, subsetting to biological cell lines present in
        the resolved graph node statistics.
        """
        input:
            source_counts = DESEQ2_SOURCE_COUNTS,
            source_metadata = DESEQ2_SOURCE_METADATA,
            node_stats  = DESEQ2_NODE_STATS_TSV
        output:
            counts_tsv   = DESEQ2_COUNTS_TSV,
            metadata_tsv = DESEQ2_METADATA_TSV,
            metadata_with_components_tsv = DESEQ2_META_COMP_TSV,
            sample_mapping_tsv = DESEQ2_SAMPLE_MAPPING_TSV,
            input_provenance_tsv = DESEQ2_INPUT_PROVENANCE_TSV
        params:
            script = os.path.join(BASE, "scripts", "prepare_deseq2_inputs.R"),
            outdir = DESEQ2_PREPARED_INPUT_DIR_ABS,
            source_count_kind = DESEQ2_INPUTS_CFG.get("source_count_kind"),
            sample_population = DESEQ2_INPUTS_CFG.get("sample_population"),
            verification_mode = DESEQ2_CELL_LINE_VERIFICATION_CFG.get("mode", "source_population_declaration"),
            verification_column = DESEQ2_CELL_LINE_VERIFICATION_CFG.get("column", ""),
            accepted_values = ",".join(DESEQ2_CELL_LINE_VERIFICATION_CFG.get("accepted_values", ["Cell Line"])),
            declared_population = DESEQ2_CELL_LINE_VERIFICATION_CFG.get("population", "cell_line_only"),
            source_sample_id_col = DESEQ2_INPUTS_CFG.get("source_sample_id_column"),
            source_cell_line_col = DESEQ2_INPUTS_CFG.get("source_cell_line_column")
        log: os.path.join(LOGROOT, "prepare_deseq2_inputs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p {params.outdir}
            Rscript {params.script} \
              --profile "{profile_name}" \
              --dsmz_counts {input.source_counts} \
              --dsmz_meta {input.source_metadata} \
              --node_stats {input.node_stats} \
              --source_count_kind "{params.source_count_kind}" \
              --source_sample_id_col "{params.source_sample_id_col}" \
              --source_cell_line_col "{params.source_cell_line_col}" \
              --sample_population "{params.sample_population}" \
              --cell_line_verification_mode "{params.verification_mode}" \
              --cell_line_verification_column "{params.verification_column}" \
              --accepted_cell_line_values "{params.accepted_values}" \
              --declared_population "{params.declared_population}" \
              --outdir {params.outdir} \
              > {log} 2>&1
            test -s {output.counts_tsv}   || (echo "ERROR: missing {output.counts_tsv}"   >&2; exit 1)
            test -s {output.metadata_tsv} || (echo "ERROR: missing {output.metadata_tsv}" >&2; exit 1)
            test -s {output.metadata_with_components_tsv} || (echo "ERROR: missing {output.metadata_with_components_tsv}" >&2; exit 1)
            test -s {output.sample_mapping_tsv} || (echo "ERROR: missing {output.sample_mapping_tsv}" >&2; exit 1)
            test -s {output.input_provenance_tsv} || (echo "ERROR: missing {output.input_provenance_tsv}" >&2; exit 1)
            '''

    # Rule: derive_isolate_list
    # Method role: DESeq2 staging rule that extracts resolved isolate cell-lines.
    # Flow: component-annotated metadata -> isolate list CSV.
    # Analysis role: defines isolate-vs-REST marker contrasts.
    rule derive_isolate_list:
        """
        Extracts the comma-separated list of isolate cell lines (is_isolate == TRUE)
        from metadata_with_components.tsv for use as --isolate_list in
        deseq2_isolate_degs.R.  Also writes a TSV manifest so the graph-derived
        cell-line selection is preserved as a tracked, inspectable pipeline output.
        Fails clearly if required columns are absent.
        """
        input:
            meta_comp = DESEQ2_META_COMP_TSV
        output:
            isolate_csv = DESEQ2_ISOLATE_CSV,
            isolate_tsv = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "isolate_cell_lines.tsv")
        params:
            cell_line_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("cell_line_column")
        log: os.path.join(LOGROOT, "derive_isolate_list.log")
        conda: CONDA_ENV_SHELL
        shell:
            r'''
            awk -F'\t' -v cell_line_col="{params.cell_line_col}" 'BEGIN {{ OFS="\t" }}
              NR==1 {{
                for(i=1;i<=NF;i++) {{
                  if($i==cell_line_col) cl=i
                  if($i=="is_isolate") iso=i
                  if($i=="component") comp=i
                  if($i=="degree") deg=i
                  if($i=="betweenness") btw=i
                }}
                if(cl==0) {{ print "ERROR: cell-line column " cell_line_col " not found in " FILENAME > "/dev/stderr"; exit 1 }}
                if(iso==0) {{ print "ERROR: column is_isolate not found in " FILENAME > "/dev/stderr"; exit 1 }}
                print "cell_line","is_isolate","component","degree","betweenness"
              }}
              NR>1 && toupper($iso)=="TRUE" {{
                cell_line=$cl
                if(!(cell_line in seen)) {{
                  seen[cell_line]=1
                  print cell_line,$iso,(comp ? $comp : ""),(deg ? $deg : ""),(btw ? $btw : "")
                }}
              }}' "{input.meta_comp}" > "{output.isolate_tsv}" 2> "{log}"
            awk -F'\t' 'NR>1 {{ print $1 }}' "{output.isolate_tsv}" | sort -u | paste -sd',' > "{output.isolate_csv}"
            test -s "{output.isolate_csv}" || (echo "ERROR: no isolate cell lines were selected" >&2; exit 1)
            '''

    # Rule: derive_components_list
    # Method role: DESeq2 staging rule that extracts resolved graph components.
    # Flow: component-annotated metadata -> component list file.
    # Analysis role: defines component-vs-rest marker contrasts.
    rule derive_components_list:
        """
        Extracts unique non-singleton component IDs from metadata_with_components.tsv.
        Singleton graph nodes are handled by isolate-vs-REST contrasts and are not
        duplicated as component-vs-rest contrasts. Writes a component-size manifest
        for reproducibility and reporting.
        """
        input:
            meta_comp = DESEQ2_META_COMP_TSV
        output:
            comp_list = DESEQ2_COMPONENTS_LIST_TXT,
            comp_tsv  = os.path.join(DESEQ2_PREPARED_INPUT_DIR, "component_cell_line_summary.tsv")
        params:
            cell_line_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("cell_line_column")
        log: os.path.join(LOGROOT, "derive_components_list.log")
        conda: CONDA_ENV_SHELL
        shell:
            r'''
            awk -F'\t' -v cell_line_col="{params.cell_line_col}" 'BEGIN {{ OFS="\t" }}
              NR==1 {{
                for(i=1;i<=NF;i++) {{
                  if($i=="component") c=i
                  if($i==cell_line_col) cl=i
                }}
                if(!c) {{
                  print "ERROR: component header missing in " FILENAME > "/dev/stderr"
                  exit 1
                }}
                if(!cl) {{
                  print "ERROR: cell-line column " cell_line_col " missing in " FILENAME > "/dev/stderr"
                  exit 1
                }}
                next
              }}
              {{
                val=$c
                gsub(/^[ \t\r\n]+/, "", val)
                gsub(/[ \t\r\n]+$/, "", val)
                if(val=="") next
                val_l=tolower(val)
                if(val_l=="na") next
                if(val !~ /^-?[0-9]+$/) {{
                  print "ERROR: non-integer component value on row " NR ": \"" val "\"" > "/dev/stderr"
                  exit 1
                }}
                component_id=val+0
                if(component_id < 0) next

                cl_val=$cl
                gsub(/^[ \t\r\n]+/, "", cl_val)
                gsub(/[ \t\r\n]+$/, "", cl_val)
                if(cl_val=="") {{
                  print "ERROR: empty cell-line identifier on row " NR > "/dev/stderr"
                  exit 1
                }}

                key=component_id SUBSEP cl_val
                if(!(key in seen)) {{
                  seen[key]=1
                  counts[component_id]++
                }}
              }}
              END {{
                print "component","n_cell_lines","eligible_for_deseq"
                for(k in counts) print k+0,counts[k],(counts[k]>=2 ? "TRUE" : "FALSE")
              }}' "{input.meta_comp}" > "{output.comp_tsv}.tmp" 2> "{log}"
            head -n 1 "{output.comp_tsv}.tmp" > "{output.comp_tsv}"
            tail -n +2 "{output.comp_tsv}.tmp" | sort -t $'\t' -k1,1n >> "{output.comp_tsv}"
            rm -f "{output.comp_tsv}.tmp"
            awk -F'\t' 'NR>1 && $3=="TRUE" {{ print $1 }}' "{output.comp_tsv}" | sort -n > "{output.comp_list}"
            '''

    # Rule: derive_anchor_list
    # Method role: DESeq2 staging rule that identifies anchor cell-lines and their components.
    # Flow: resolved-graph anchor audit -> anchor list and component mapping.
    # Supports anchor-aware directional marker tables.
    rule derive_anchor_list:
        """
        Extracts selected graph anchors from the resolved-graph anchor audit.
        Selection logic is schema-aware:
          - use anchor_selected when present;
          - otherwise use the union of most_connected_selected and
            canonical_bridge_selected for compatibility with the older audit schema.
        Duplicate node IDs are removed from the comma-separated anchor list.
        Two outputs are written for deseq2_isolate_degs.R:
          anchor_list.csv        -- comma-separated selected anchor node IDs
          anchor_components.tsv  -- two-column TSV (anchor<TAB>component)
        Selection-reason columns remain in the source anchor audit; they are not
        copied into anchor_components.tsv.
        """
        input:
            anchor_audit = P_CONS_ANCHOR_AUDIT_TSV
        output:
            anchor_list_csv       = DESEQ2_ANCHOR_LIST_CSV,
            anchor_components_tsv = DESEQ2_ANCHOR_COMPONENTS_TSV
        log: os.path.join(LOGROOT, "derive_anchor_list.log")
        conda: CONDA_ENV_SHELL
        shell:
            r'''
            awk -F'\t' 'BEGIN {{ OFS="\t" }}
              NR==1 {{
                for(i=1;i<=NF;i++) {{
                  if($i=="node_id")                  nid=i
                  if($i=="component_id")             cid=i
                  if($i=="anchor_selected")          asel=i
                  if($i=="most_connected_selected")  mcs=i
                  if($i=="canonical_bridge_selected") cbs=i
                }}
                if(!nid) {{ print "ERROR: column node_id not found" > "/dev/stderr"; exit 1 }}
                if(!cid) {{ print "ERROR: column component_id not found" > "/dev/stderr"; exit 1 }}
                if(!asel && !mcs && !cbs) {{
                  print "ERROR: no anchor selection column found (expected anchor_selected, most_connected_selected, or canonical_bridge_selected)" > "/dev/stderr"
                  exit 1
                }}
                print "anchor","component"
                next
              }}
              NR>1 {{
                if(asel) {{
                  selected = (toupper($asel)=="TRUE")
                }} else {{
                  selected = (toupper($mcs)=="TRUE" || toupper($cbs)=="TRUE")
                }}
                if(selected) print $nid,$cid
              }}
            ' "{input.anchor_audit}" > "{output.anchor_components_tsv}" 2> "{log}"
            awk -F'\t' 'NR>1 {{ print $1 }}' "{output.anchor_components_tsv}" \
              | sort -u | paste -sd',' > "{output.anchor_list_csv}"
            test -s "{output.anchor_list_csv}" || (echo "ERROR: no anchor cell lines were selected" >&2; exit 1)
            echo "[derive_anchor_list] Schema: $(head -1 "{input.anchor_audit}" | tr '\t' '\n' | grep -n 'anchor_selected\|most_connected\|canonical_bridge' | head -5)" >> "{log}"
            echo "[derive_anchor_list] Anchors: $(cat "{output.anchor_list_csv}")" >> "{log}"
            echo "[derive_anchor_list] anchor_components rows: $(wc -l < "{output.anchor_components_tsv}")" >> "{log}"
            '''

    # Rule: derive_graph_based_deseq2_contrast_definitions
    # Scientific purpose: materialise graph-topology-derived isolate and anchor
    # contrast definitions before any count modelling.
    # Unit of analysis: cell-line profile and resolved graph component.
    # Transformation: resolved graph labels and anchor audit -> explicit
    # isolate list, anchor list, and anchor-to-component mapping.
    rule derive_graph_based_deseq2_contrast_definitions:
        input:
            isolate_csv = DESEQ2_ISOLATE_CSV,
            anchor_list_csv = DESEQ2_ANCHOR_LIST_CSV,
            anchor_components_tsv = DESEQ2_ANCHOR_COMPONENTS_TSV

    # Rule: run_isolate_and_anchor_deseq2_contrasts
    # Scientific purpose: fit focal-versus-reference DESeq2 Wald contrasts and
    # prioritise retained markers within each contrast.
    # Unit of analysis: one isolate or anchor contrast within the active cancer type.
    # Transformation: prepared raw counts, metadata, and graph-derived contrast
    # definitions -> retained per-contrast marker tables plus manifest.
    rule run_isolate_and_anchor_deseq2_contrasts:
        """
        Runs within-cancer-type DESeq2 marker contrasts derived from the resolved
        cell-line graph:
          - focal isolate vs all other same-cancer cell-line profiles;
          - focal anchor vs same-cancer profiles outside its resolved component.
        Produces all-gene DESeq2 tables, retained marker gene lists, retained
        marker tables, size-factor QC, and a contrast-level marker manifest.
        """
        input:
            validation            = DESEQ2_INPUT_VALIDATION_REQUIRED_TSV,
            counts_tsv            = DESEQ2_COUNTS_TSV,
            meta_comp             = DESEQ2_META_COMP_TSV,
            isolate_csv           = DESEQ2_ISOLATE_CSV,
            anchor_list_csv       = DESEQ2_ANCHOR_LIST_CSV,
            anchor_components_tsv = DESEQ2_ANCHOR_COMPONENTS_TSV
        output:
            size_factors = os.path.join(DESEQ2_ISOLATE_DIR, "qc", "size_factors.tsv"),
            contrast_manifest = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "contrast_level_marker_manifest.tsv"),
            session_info = os.path.join(DESEQ2_ISOLATE_DIR, "sessionInfo.txt")
        params:
            script        = os.path.join(BASE, "scripts", "deseq2_isolate_degs.R"),
            outdir        = DESEQ2_ISOLATE_DIR_ABS,
            cancer_type   = cfg.get("analysis", {}).get("cancer_type", profile_name).lower(),
            sample_id_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("sample_id_column"),
            cell_line_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("cell_line_column"),
            component_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("component_column"),
            isolate_adjusted_p_value_threshold = DESEQ2_ISOLATE_CONTRAST_CFG.get("adjusted_p_value_threshold"),
            isolate_minimum_absolute_shrunken_log2fc = DESEQ2_ISOLATE_CONTRAST_CFG.get("minimum_absolute_shrunken_log2fc"),
            isolate_maximum_markers_per_contrast = DESEQ2_ISOLATE_CONTRAST_CFG.get("maximum_markers_per_contrast"),
            anchor_adjusted_p_value_threshold = DESEQ2_ANCHOR_CONTRAST_CFG.get("adjusted_p_value_threshold"),
            anchor_minimum_absolute_shrunken_log2fc = DESEQ2_ANCHOR_CONTRAST_CFG.get("minimum_absolute_shrunken_log2fc"),
            anchor_maximum_markers_per_contrast = DESEQ2_ANCHOR_CONTRAST_CFG.get("maximum_markers_per_contrast"),
            minimum_base_mean = DESEQ2_MARKER_PRIORITISATION_CFG.get("minimum_base_mean"),
            minimum_total_gene_count = DESEQ2_MARKER_PRIORITISATION_CFG.get("minimum_total_gene_count"),
            dispersion_fit_type = DESEQ2_MARKER_PRIORITISATION_CFG.get("dispersion_fit_type"),
            lfc_shrinkage_method = DESEQ2_MARKER_PRIORITISATION_CFG.get("lfc_shrinkage_method")
        log: os.path.join(LOGROOT, "deseq2_isolate_degs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            if [ ! -s {input.validation} ]; then
                echo "ERROR: DESeq2 input validation file is missing: {input.validation}" >&2
                exit 1
            fi
            if awk -F'\t' 'NR>1 && $2!="PASS" {{ bad=1 }} END {{ exit bad ? 0 : 1 }}' {input.validation}; then
                echo "ERROR: DESeq2 inputs failed validation; see {input.validation}" >&2
                awk -F'\t' 'NR==1 || $2!="PASS"' {input.validation} >&2
                exit 1
            fi
            ISOLATE_LIST=$(cat {input.isolate_csv})
            ANCHOR_LIST=$(cat {input.anchor_list_csv})
            if [ -z "$ISOLATE_LIST" ] && [ -z "$ANCHOR_LIST" ]; then
                echo "No isolate or anchor contrasts were derived for {params.cancer_type}." > {log}
                echo "This graph-derived marker analysis requires at least one contrast and the current graph state requires methodological review." >> {log}
                exit 1
            fi
            STAGING="{params.outdir}.staging_current"
            rm -rf "$STAGING"
            mkdir -p "$STAGING"
            Rscript {params.script} \
              --counts {input.counts_tsv} \
              --meta {input.meta_comp} \
              --sample_id_col {params.sample_id_col} \
              --cell_line_col {params.cell_line_col} \
              --component_col {params.component_col} \
              --cancer_type {params.cancer_type} \
              --isolate_list "$ISOLATE_LIST" \
              --anchor_list "$ANCHOR_LIST" \
              --anchor_components {input.anchor_components_tsv} \
              --isolate_adjusted_p_value_threshold {params.isolate_adjusted_p_value_threshold} \
              --isolate_minimum_absolute_shrunken_log2fc {params.isolate_minimum_absolute_shrunken_log2fc} \
              --isolate_maximum_markers_per_contrast {params.isolate_maximum_markers_per_contrast} \
              --anchor_adjusted_p_value_threshold {params.anchor_adjusted_p_value_threshold} \
              --anchor_minimum_absolute_shrunken_log2fc {params.anchor_minimum_absolute_shrunken_log2fc} \
              --anchor_maximum_markers_per_contrast {params.anchor_maximum_markers_per_contrast} \
              --minimum_base_mean {params.minimum_base_mean} \
              --minimum_total_gene_count {params.minimum_total_gene_count} \
              --dispersion_fit_type {params.dispersion_fit_type} \
              --lfc_shrinkage_method {params.lfc_shrinkage_method} \
              --outdir "$STAGING" \
              > {log} 2>&1
            test -s "$STAGING/markers/contrast_level_marker_manifest.tsv" || (echo "ERROR: missing current canonical marker manifest" >&2; exit 1)
            rsync -a --delete "$STAGING"/ "{params.outdir}"/
            rm -rf "$STAGING"
            '''

    # Rule: build_enrichment_recurrence_sets
    # Scientific purpose: provide a retained-marker recurrence sidecar for
    # enrichment-query construction without making recurrence a DESeq2 output.
    # Unit of analysis: retained per-contrast marker list within one cancer type.
    # Transformation: contrast-level marker manifest -> recurrence count table
    # and configured recurrence-threshold gene list for enrichment utilities.
    rule build_enrichment_recurrence_sets:
        input:
            manifest = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "contrast_level_marker_manifest.tsv")
        output:
            recurrence = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "gene_recurrence_across_contrasts.tsv"),
            unique_set = os.path.join(
                DESEQ2_ISOLATE_DIR,
                "markers",
                f"unique_feature_set_recurrence_ge_{DESEQ2_RETAINED_MARKER_RECURRENCE_THRESHOLD}.txt"
            )
        params:
            deseq2_dir = DESEQ2_ISOLATE_DIR_ABS,
            recurrence_threshold = DESEQ2_RETAINED_MARKER_RECURRENCE_THRESHOLD
        log:
            os.path.join(LOGROOT, "build_enrichment_recurrence_sets.log")
        conda:
            CONDA_ENV_PY
        shell:
            r'''
            python - <<'PY' > "{log}" 2>&1
            from collections import Counter
            from pathlib import Path
            import pandas as pd

            manifest = pd.read_csv("{input.manifest}", sep="\t")
            deseq2_dir = Path("{params.deseq2_dir}")
            recurrence_threshold = int("{params.recurrence_threshold}")
            counts = Counter()
            for row in manifest.itertuples(index=False):
                marker_list = getattr(row, "marker_gene_list_path", "")
                if not marker_list:
                    continue
                path = Path(marker_list)
                if not path.is_absolute():
                    path = deseq2_dir / path
                with path.open() as handle:
                    genes = {{line.strip().split(".", 1)[0] for line in handle if line.strip()}}
                counts.update(genes)
            rows = [
                {{"gene_id": gene_id, "retained_marker_list_count": count, "freq": count}}
                for gene_id, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            ]
            pd.DataFrame(rows, columns=["gene_id", "retained_marker_list_count", "freq"]).to_csv(
                "{output.recurrence}", sep="\t", index=False
            )
            selected = [row["gene_id"] for row in rows if row["retained_marker_list_count"] >= recurrence_threshold]
            Path("{output.unique_set}").write_text(("\n".join(selected) + "\n") if selected else "")
            print(f"[Marker recurrence] retained_marker_list_count >= {{recurrence_threshold}}: {{len(selected)}} genes")
            PY
            '''

    DESEQ2_DIRECTIONAL_MARKER_DIR = os.path.join(DESEQ2_ISOLATE_DIR, "directional_markers")
    DESEQ2_DIRECTIONAL_MARKER_DIR_ABS = os.path.join(DESEQ2_ISOLATE_DIR_ABS, "directional_markers")

    # Rule: write_deseq2_directional_marker_tables
    # Method role: marker-reporting rule converting isolate and anchor DESeq2 outputs into direction-specific tables.
    # Flow: completed marker manifest, retained marker lists, and DESeq2 result tables -> directional marker tables.
    # Provides marker post-processing and enrichment inputs.
    rule write_deseq2_directional_marker_tables:
        """
        Organise retained DESeq2 marker genes into directional marker tables
        based on the sign of log2FoldChange. This is a post-DESeq2 marker step:
        it consumes the completed marker manifest, marker gene lists, and full
        DESeq2 tables without rerunning DESeq2.
        """
        input:
            manifest = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "contrast_level_marker_manifest.tsv"),
            markers_dir = os.path.join(DESEQ2_ISOLATE_DIR, "markers"),
            tables_dir = os.path.join(DESEQ2_ISOLATE_DIR, "tables"),
            deseq2_script = os.path.join(BASE, "scripts", "deseq2_isolate_degs.R")
        output:
            orientation = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "directional_marker_orientation.txt"),
            all_up = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "all_upregulated_markers.tsv"),
            all_down = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "all_downregulated_markers.tsv"),
            all_mixed = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "all_mixed_direction_markers.tsv"),
            counts = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "directional_marker_counts.tsv"),
            top_up = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "top_upregulated_markers.tsv"),
            top_down = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "top_downregulated_markers.tsv"),
            top_mixed = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "top_mixed_direction_markers.tsv")
        params:
            script = os.path.join(BASE, "scripts", "write_deseq2_directional_marker_tables.R"),
            cohort = profile_name,
            deseq2_dir = DESEQ2_ISOLATE_DIR_ABS,
            outdir = DESEQ2_DIRECTIONAL_MARKER_DIR_ABS
        log:
            os.path.join(LOGROOT, "write_deseq2_directional_marker_tables.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript {params.script} \
              --cohort {params.cohort} \
              --deseq2_dir {params.deseq2_dir} \
              --deseq2_script {input.deseq2_script} \
              --outdir {params.outdir} \
              > {log} 2>&1
            '''

    # Rule: write_all_deseq2_directional_marker_tables
    # Method role: collector target requiring all direction-specific marker tables.
    # Flow: expected marker-table outputs -> Snakemake dependency aggregation only.
    # Stage output: provides a marker target without running a new analysis.
    rule write_all_deseq2_directional_marker_tables:
        """
        Aggregate target for the active profile's directional marker files.
        """
        input:
            rules.write_deseq2_directional_marker_tables.output

    # Rule: deseq2_component_vs_rest_all
    # Method role: marker-analysis collector that completes component-vs-rest DESeq2 contrasts.
    # Flow: component marker outputs -> completion sentinel.
    # Provides marker-analysis completion signal for PIPELINE_TARGET assembly.
    rule deseq2_component_vs_rest_all:
        """
        Runs deseq2_component_vs_rest.R once per component in components_list.txt.
        A shell loop iterates over component IDs, running the R script for each.
        A sentinel .done file is written only after all components complete.
        """
        input:
            validation = DESEQ2_INPUT_VALIDATION_REQUIRED_TSV,
            counts_tsv = DESEQ2_COUNTS_TSV,
            meta_comp  = DESEQ2_META_COMP_TSV,
            comp_list  = DESEQ2_COMPONENTS_LIST_TXT
        output:
            done = touch(os.path.join(DESEQ2_COMP_DIR, ".done"))
        params:
            script        = os.path.join(BASE, "scripts", "deseq2_component_vs_rest.R"),
            outdir        = DESEQ2_COMP_DIR_ABS,
            component_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("component_column"),
            sample_id_col = DESEQ2_MARKER_PRIORITISATION_CFG.get("sample_id_column"),
            adjusted_p_value_threshold = DESEQ2_COMPONENT_CONTRAST_CFG.get("adjusted_p_value_threshold"),
            minimum_absolute_shrunken_log2fc = DESEQ2_COMPONENT_CONTRAST_CFG.get("minimum_absolute_shrunken_log2fc"),
            maximum_markers_per_direction = DESEQ2_COMPONENT_CONTRAST_CFG.get("maximum_markers_per_direction"),
            minimum_base_mean = DESEQ2_MARKER_PRIORITISATION_CFG.get("minimum_base_mean"),
            minimum_total_gene_count = DESEQ2_MARKER_PRIORITISATION_CFG.get("minimum_total_gene_count"),
            dispersion_fit_type = DESEQ2_MARKER_PRIORITISATION_CFG.get("dispersion_fit_type"),
            lfc_shrinkage_method = DESEQ2_MARKER_PRIORITISATION_CFG.get("lfc_shrinkage_method")
        log: os.path.join(LOGROOT, "deseq2_component_vs_rest_all.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            if [ ! -s {input.validation} ]; then
                echo "ERROR: DESeq2 input validation file is missing: {input.validation}" >&2
                exit 1
            fi
            if awk -F'\t' 'NR>1 && $2!="PASS" {{ bad=1 }} END {{ exit bad ? 0 : 1 }}' {input.validation}; then
                echo "ERROR: DESeq2 inputs failed validation; see {input.validation}" >&2
                awk -F'\t' 'NR==1 || $2!="PASS"' {input.validation} >&2
                exit 1
            fi
            STAGING="{params.outdir}.staging_current"
            rm -rf "$STAGING"
            mkdir -p "$STAGING"
            echo "[START] Component-vs-rest loop" > {log}
            if [ ! -s {input.comp_list} ]; then
                echo "[WARN] No components to process; sentinel written" >> {log}
                rsync -a --delete "$STAGING"/ "{params.outdir}"/
                exit 0
            fi
            while IFS= read -r COMP; do
                [ -z "$COMP" ] && continue
                echo "[INFO] Running component $COMP" >> {log}
                Rscript {params.script} \
                  --counts {input.counts_tsv} \
                  --meta {input.meta_comp} \
                  --component "$COMP" \
                  --component_col {params.component_col} \
                  --sample_id_col {params.sample_id_col} \
                  --outdir "$STAGING" \
                  --adjusted_p_value_threshold {params.adjusted_p_value_threshold} \
                  --minimum_absolute_shrunken_log2fc {params.minimum_absolute_shrunken_log2fc} \
                  --maximum_markers_per_direction {params.maximum_markers_per_direction} \
                  --minimum_base_mean {params.minimum_base_mean} \
                  --minimum_total_gene_count {params.minimum_total_gene_count} \
                  --dispersion_fit_type {params.dispersion_fit_type} \
                  --lfc_shrinkage_method {params.lfc_shrinkage_method} \
                  >> {log} 2>&1
            done < {input.comp_list}
            rsync -a --delete "$STAGING"/ "{params.outdir}"/
            rm -rf "$STAGING"
            echo "[DONE] All components processed" >> {log}
            '''


# =============================================================================
# ENRICHMENT ANALYSIS
# =============================================================================
# Stage role: prepares, executes, aggregates, and reports canonical graph-derived
# functional enrichment for the pan-cancer selected feature panel.
#
# The historical per-profile enrichment routes are intentionally not declared.
# The canonical rules are included from rules/graph_derived_functional_enrichment.smk
# after pan-cancer feature-interface variables have been defined.

FUNCTIONAL_ENRICHMENT_CFG = cfg.get("functional_enrichment", {})
FUNCTIONAL_ENRICHMENT_ENABLED = bool(FUNCTIONAL_ENRICHMENT_CFG.get("enabled", False))

FUNCTIONAL_ENRICHMENT_DIR_REL = os.path.join(
    "results",
    "unsupervised",
    "pan_cancer",
    FUNCTIONAL_ENRICHMENT_CFG.get("output_dir_name", "enrichment"),
)
FUNCTIONAL_ENRICHMENT_DIR_ABS = abspath(FUNCTIONAL_ENRICHMENT_DIR_REL)
FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL = os.path.join(FUNCTIONAL_ENRICHMENT_DIR_REL, "query_sets")
FUNCTIONAL_ENRICHMENT_QUERY_DIR_ABS = abspath(FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL)
FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_REL = os.path.join(FUNCTIONAL_ENRICHMENT_DIR_REL, "gprofiler")
FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_ABS = abspath(FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_REL)
FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL = os.path.join(FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_REL, "aggregate")
FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_ABS = abspath(FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL)
FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL = os.path.join(FUNCTIONAL_ENRICHMENT_DIR_REL, "figures")
FUNCTIONAL_ENRICHMENT_FIGURE_DIR_ABS = abspath(FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL)

FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST = os.path.join(
    FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, "enrichment_query_manifest.tsv"
)
FUNCTIONAL_ENRICHMENT_QUERY_INTEGRITY = os.path.join(
    FUNCTIONAL_ENRICHMENT_QUERY_DIR_REL, "enrichment_query_background_integrity.tsv"
)
FUNCTIONAL_ENRICHMENT_GPROFILER_BATCH_SUCCESS = os.path.join(
    FUNCTIONAL_ENRICHMENT_GPROFILER_DIR_REL,
    "batches",
    "batch_0001_999999",
    ".success",
)
FUNCTIONAL_ENRICHMENT_AGGREGATE_SUCCESS = os.path.join(
    FUNCTIONAL_ENRICHMENT_AGGREGATE_DIR_REL, ".success"
)
FUNCTIONAL_ENRICHMENT_PRIMARY_HEATMAP_PDF = os.path.join(
    FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_top_terms_heatmap.pdf"
)
FUNCTIONAL_ENRICHMENT_CONTRAST_SUPPORT_HEATMAP_PDF = os.path.join(
    FUNCTIONAL_ENRICHMENT_FIGURE_DIR_REL, "Fig_enrichment_contrast_term_support_heatmap.pdf"
)


# =============================================================================
# MARKER POST-PROCESSING
# =============================================================================
# Stage role: converts cohort marker tables into marker-derived feature panels and pan-cancer outputs.

MARKER_POST_CFG = config.get("defaults", {}).get("marker_postprocessing", {})
# Boundary note: pan_cancer is the first-class profile for marker-derived
# pan-cancer outputs. These rules are declared only for the pan_cancer profile.
MARKER_POST_ENABLED = bool(MARKER_POST_CFG.get("enabled", False)) and DECLARE_PAN_CANCER_RULES
# The feature-panel rule alone is additionally declared for any profile whose
# configured representation universe contains a PanCancerFeatureSet
# representation, so the panel is a build dependency of that profile's
# representation, graph and resolution stages instead of a pre-existing file.
# Only the panel is declared here: the downstream pan-cancer expression,
# alignment, mapping and enrichment rules stay gated on MARKER_POST_ENABLED so
# no additional pan-cancer targets enter another profile's DAG.
MARKER_POST_FEATURE_PANEL_ENABLED = bool(MARKER_POST_CFG.get("enabled", False)) and (
    DECLARE_PAN_CANCER_RULES or REQUIRES_PAN_CANCER_FEATURE_PANEL
)
if REQUIRES_PAN_CANCER_FEATURE_PANEL and not MARKER_POST_FEATURE_PANEL_ENABLED:
    raise ValueError(
        f"Profile '{profile_name}' declares PanCancerFeatureSet representation(s) "
        f"{PAN_CANCER_FEATURE_SET_DIRECTIONS}, but "
        "defaults.marker_postprocessing.enabled is false, so the required "
        "pan-cancer feature panel has no producing rule."
    )

if MARKER_POST_FEATURE_PANEL_ENABLED:
    CONSENSUS_CFG = MARKER_POST_CFG.get("consensus", {})
    PAN_CANCER_MP_CFG = MARKER_POST_CFG.get("pan_cancer", {})
    PAN_EXPR_CFG = MARKER_POST_CFG.get("expression_matrix", {})
    MAPPING_CFG = MARKER_POST_CFG.get("mapping", {})

    PAN_PROFILES = [str(p) for p in PAN_CANCER_MP_CFG.get("disease_profiles", [])]
    CONSENSUS_OUTDIR_NAME = CONSENSUS_CFG.get("outdir_name", "consensus_markers")

    def profile_unsup_root_abs(profile):
        pcfg = get_profile_cfg(profile)
        unsup = pcfg.get("paths", {}).get("unsup_root", f"results/unsupervised/{profile}")
        if not os.path.isabs(unsup):
            unsup = os.path.join(PIPE_ROOT, unsup)
        return unsup

    def profile_consensus_outdir(profile):
        return os.path.join(profile_unsup_root_abs(profile), "deseq2_markers", CONSENSUS_OUTDIR_NAME)

    def profile_consensus_summary_dir(profile):
        return os.path.join(profile_consensus_outdir(profile), "summary")

    def profile_contrast_marker_manifest_rel(profile):
        return os.path.join(
            "results", "unsupervised", profile,
            "deseq2_markers", "markers", "contrast_level_marker_manifest.tsv"
        )

    def profile_marker_session_info_rel(profile):
        return os.path.join(
            "results", "unsupervised", profile,
            "deseq2_markers", "sessionInfo.txt"
        )

    def profile_summary_dir(profile):
        return profile_consensus_summary_dir(profile)

    PAN_FEATURES_OUTDIR = PAN_CANCER_MP_CFG.get("outdir", "results/unsupervised/pan_cancer/feature_space")
    PAN_FEATURES_OUTDIR_ABS = abspath(PAN_FEATURES_OUTDIR)
    PAN_FEATURES_TSV = PAN_CANCER_MP_CFG.get("final_features_tsv", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.tsv"))
    PAN_FEATURES_CLEAN = PAN_CANCER_MP_CFG.get(
        "final_features_clean",
        os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features_clean.txt")
    )
    PAN_FEATURES_UP = PAN_CANCER_MP_CFG.get("final_features_up", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.UP.txt"))
    PAN_FEATURES_DOWN = PAN_CANCER_MP_CFG.get("final_features_down", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.DOWN.txt"))
    PAN_FEATURES_MIXED = PAN_CANCER_MP_CFG.get("final_features_mixed", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.MIXED.txt"))
    PAN_FEATURES_SUMMARY = PAN_CANCER_MP_CFG.get("build_summary_tsv", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_feature_build_summary.tsv"))
    PAN_FEATURES_GENE_EVIDENCE = os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_feature_gene_evidence.tsv")
    PAN_FEATURES_METHOD = PAN_CANCER_MP_CFG.get("method", "graph_derived_pan_cancer_feature_selection_v1_revised")
    PAN_FEATURES_PREFIX = PAN_CANCER_MP_CFG.get("audit_output_prefix", "ranked_marker_source_panel")
    PAN_FEATURES_QUANTILE_CFG = PAN_CANCER_MP_CFG.get("empirical_quantile_thresholds", {})
    _PAN_FEATURES_REQUIRED_QUANTILES = {
        "adjusted_p_value_quantile",
        "absolute_shrunken_log2fc_quantile",
        "expression_quantile",
    }
    _PAN_FEATURES_QUANTILE_KEYS = set(PAN_FEATURES_QUANTILE_CFG)
    if _PAN_FEATURES_QUANTILE_KEYS != _PAN_FEATURES_REQUIRED_QUANTILES:
        raise ValueError(
            "marker_postprocessing.pan_cancer.empirical_quantile_thresholds must define exactly: "
            + ", ".join(sorted(_PAN_FEATURES_REQUIRED_QUANTILES))
        )
    PAN_FEATURES_RANKED_BY_COHORT = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_by_cohort.tsv")
    PAN_FEATURES_RANKED_BY_MARKER_SOURCE_CLASS = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_by_marker_source_class.tsv")
    PAN_FEATURES_RANKED_BY_FEATURE_CLASS = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_by_feature_class.tsv")
    PAN_FEATURES_CANDIDATE_POOL_EVIDENCE = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_candidate_pool_evidence.tsv")
    PAN_FEATURES_EMPIRICAL_THRESHOLDS = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_empirical_quantile_thresholds.tsv")
    PAN_FEATURES_CANDIDATE_ACCEPTANCE = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_candidate_acceptance.tsv")
    PAN_FEATURES_SELECTED_EVIDENCE_ROWS = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_selected_evidence_rows.tsv")
    PAN_FEATURES_VALIDATION = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_validation.tsv")
    PAN_FEATURES_RANKED_MANIFEST = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_run_manifest.tsv")
    PAN_FEATURES_ACTIVE_MANIFEST = os.path.join(PAN_FEATURES_OUTDIR, f"{PAN_FEATURES_PREFIX}_active_directory_manifest.tsv")

    def _profile_dir_args():
        args = []
        for prof in PAN_PROFILES:
            args.append(f"--profile-dir {prof}={shlex.quote(profile_summary_dir(prof))}")
        return " ".join(args)

    def profile_contrast_marker_manifest_abs(profile):
        return os.path.join(
            profile_unsup_root_abs(profile),
            "deseq2_markers", "markers", "contrast_level_marker_manifest.tsv"
        )

    def _profile_marker_manifest_args():
        args = []
        for prof in PAN_PROFILES:
            args.append(
                f"--profile-marker-manifest {prof}="
                f"{shlex.quote(profile_contrast_marker_manifest_abs(prof))}"
            )
        return " ".join(args)

    # Rule: construct_pan_cancer_feature_panel
    # Scientific purpose: aggregate retained contrast-level marker evidence,
    # classify recurrent/singleton/non-recurrent rows, apply all-three empirical
    # candidate acceptance, and export F = R union S union N.
    # Unit of analysis: cancer type x graph-derived marker-evidence stratum x gene.
    # Transformation: contrast-level marker manifests and retained marker tables
    # -> selected pan-cancer feature panel plus TSV audit tables.
    rule construct_pan_cancer_feature_panel:
        """Build the graph-derived pan-cancer feature panel from retained contrast-level marker evidence."""
        input:
            cfg = CFGFILE_ABS,
            marker_manifests=[
                profile_contrast_marker_manifest_rel(p)
                for p in PAN_PROFILES
            ],
            marker_session_info=[
                profile_marker_session_info_rel(p)
                for p in PAN_PROFILES
            ]
        output:
            features_tsv = PAN_FEATURES_TSV,
            clean_txt   = PAN_FEATURES_CLEAN,
            up_txt      = PAN_FEATURES_UP,
            down_txt    = PAN_FEATURES_DOWN,
            mixed_txt   = PAN_FEATURES_MIXED,
            summary_tsv = PAN_FEATURES_SUMMARY,
            gene_evidence_tsv = PAN_FEATURES_GENE_EVIDENCE,
            ranked_by_cohort_tsv = PAN_FEATURES_RANKED_BY_COHORT,
            ranked_by_marker_source_class_tsv = PAN_FEATURES_RANKED_BY_MARKER_SOURCE_CLASS,
            ranked_by_feature_class_tsv = PAN_FEATURES_RANKED_BY_FEATURE_CLASS,
            candidate_pool_evidence_tsv = PAN_FEATURES_CANDIDATE_POOL_EVIDENCE,
            empirical_quantile_thresholds_tsv = PAN_FEATURES_EMPIRICAL_THRESHOLDS,
            candidate_acceptance_tsv = PAN_FEATURES_CANDIDATE_ACCEPTANCE,
            selected_evidence_rows_tsv = PAN_FEATURES_SELECTED_EVIDENCE_ROWS,
            validation_tsv = PAN_FEATURES_VALIDATION,
            ranked_manifest_tsv = PAN_FEATURES_RANKED_MANIFEST,
            active_directory_manifest_tsv = PAN_FEATURES_ACTIVE_MANIFEST,
            done_file   = os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features_done.txt")
        params:
            script=os.path.join(SCRIPTS_DIR, "build_pan_cancer_features.py"),
            outdir=PAN_FEATURES_OUTDIR_ABS,
            method=PAN_FEATURES_METHOD,
            profile_marker_manifests=_profile_marker_manifest_args(),
            adjusted_p_value_quantile=PAN_FEATURES_QUANTILE_CFG.get("adjusted_p_value_quantile"),
            abs_log2fc_quantile=PAN_FEATURES_QUANTILE_CFG.get("absolute_shrunken_log2fc_quantile"),
            expression_quantile=PAN_FEATURES_QUANTILE_CFG.get("expression_quantile"),
            audit_output_prefix=PAN_FEATURES_PREFIX,
            remove_ribo_mt="--remove-ribo-mt" if PAN_CANCER_MP_CFG.get("remove_ribo_mt", False) else "",
            gene_annot=("--gene-annotation-tsv " + abspath(PAN_CANCER_MP_CFG["gene_annotation_tsv"]))
                        if PAN_CANCER_MP_CFG.get("gene_annotation_tsv") else ""
        log: os.path.join(LOGROOT, "build_pan_cancer_features.log")
        conda: CONDA_ENV_PY
        shell:
            r'''
            STAGING="{params.outdir}.staging_current"
            rm -rf "$STAGING"
            mkdir -p "$STAGING"
            python "{params.script}" \
              {params.profile_marker_manifests} \
              --output-dir "$STAGING" \
              --method "{params.method}" \
              --adjusted-p-value-quantile {params.adjusted_p_value_quantile} \
              --absolute-shrunken-log2fc-quantile {params.abs_log2fc_quantile} \
              --expression-quantile {params.expression_quantile} \
              --audit-output-prefix "{params.audit_output_prefix}" \
              {params.remove_ribo_mt} \
              {params.gene_annot} \
              > "{log}" 2>&1
            test -s "$STAGING/pan_cancer_features.tsv" || (echo "ERROR: missing staged pan_cancer_features.tsv" >&2; exit 1)
            test -s "$STAGING/pan_cancer_features_clean.txt" || (echo "ERROR: missing staged pan_cancer_features_clean.txt" >&2; exit 1)
            rsync -a --delete "$STAGING"/ "{params.outdir}"/
            rm -rf "$STAGING"
            test -s "{output.features_tsv}" || (echo "ERROR: missing {output.features_tsv}" >&2; exit 1)
            '''

# The remaining marker-postprocessing rules build pan-cancer expression,
# alignment, mapping and ranking products. They belong to the pan_cancer
# profile only; profiles that merely consume the feature panel stop above.
if MARKER_POST_ENABLED:
    PAN_EXPR_RDS = PAN_EXPR_CFG.get("output_rds", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds")
    PAN_EXPR_META = PAN_EXPR_CFG.get("output_meta_tsv", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_metadata.tsv")
    PAN_EXPR_CELL_LINES_RDS = PAN_EXPR_CFG.get(
        "output_cell_lines_rds",
        "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds"
    )
    PAN_EXPR_CELL_LINES_META = PAN_EXPR_CFG.get(
        "output_cell_lines_meta_tsv",
        "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only_metadata.tsv"
    )

    # Precompute optional CLI fragments for the expression matrix rule.
    # Snakemake shell blocks only support simple {name} placeholders, NOT
    # inline Python expressions like {expr if cond else ""}.
    # HEME is cell-line-only in the current inputs, so it is excluded from the
    # main tumour/cell-line alignment matrix and reserved for the cell-line-only
    # matrix below.
    _PAN_EXPR_HEME_PATH = ""
    _PAN_EXPR_CELL_LINES_HEME_PATH = profile_vst_joint_abs("heme") if PAN_EXPR_CFG.get("include_heme_cell_lines_only", True) else ""
    _PAN_EXPR_CELL_LINES_HEME_INPUT = [_PAN_EXPR_CELL_LINES_HEME_PATH] if _PAN_EXPR_CELL_LINES_HEME_PATH else []
    _PAN_EXPR_METADATA_PATH = abspath(PAN_EXPR_CFG.get("metadata_fallback_csv", "")) if PAN_EXPR_CFG.get("metadata_fallback_csv") else ""
    _PAN_EXPR_HEME_ARG = f'--heme-vst "{_PAN_EXPR_HEME_PATH}"' if _PAN_EXPR_HEME_PATH else ""
    _PAN_EXPR_CELL_LINES_HEME_ARG = f'--heme-vst "{_PAN_EXPR_CELL_LINES_HEME_PATH}"' if _PAN_EXPR_CELL_LINES_HEME_PATH else ""
    _PAN_EXPR_METADATA_ARG = f'--metadata "{_PAN_EXPR_METADATA_PATH}"' if _PAN_EXPR_METADATA_PATH else ""
    _PAN_EXPR_OUTPUT_META_ARG = f'--output-metadata "{PAN_EXPR_META}"' if PAN_EXPR_META else ""
    _PAN_EXPR_CELL_LINES_OUTPUT_META_ARG = f'--output-metadata "{PAN_EXPR_CELL_LINES_META}"' if PAN_EXPR_CELL_LINES_META else ""

    # Precompute concrete input paths (avoid lambdas in input:).
    _PAN_EXPR_BRCA_VST = profile_vst_joint_abs("brca")
    _PAN_EXPR_NBL_VST = profile_vst_joint_abs("nbl")
    _PAN_EXPR_RBL_VST = profile_vst_joint_abs("rbl")

    # Main expression inputs are limited to cohorts with tumour samples.
    _PAN_EXPR_INPUTS = [_PAN_EXPR_BRCA_VST, _PAN_EXPR_NBL_VST, _PAN_EXPR_RBL_VST]
    if _PAN_EXPR_HEME_PATH:
        _PAN_EXPR_INPUTS.append(_PAN_EXPR_HEME_PATH)

    # Rule: build_pan_cancer_feature_expression_matrix
    # Method role: expression-construction rule that filters tumour/cell-line matrices to the marker-derived feature panel.
    # Flow: cohort VST inputs plus feature table -> pan-cancer feature-space expression matrix.
    # Provides input for UMAP, mapping, ranking, and pan-cancer graph analyses.
    rule build_pan_cancer_feature_expression_matrix:
        """Combine per-profile VST inputs into the feature-limited pan-cancer matrix."""
        input:
            genes = PAN_FEATURES_CLEAN,
            brca = _PAN_EXPR_BRCA_VST,
            nbl  = _PAN_EXPR_NBL_VST,
            rbl  = _PAN_EXPR_RBL_VST
        output:
            expr_rds = PAN_EXPR_RDS,
            meta_tsv = PAN_EXPR_META
        params:
            script = os.path.join(SCRIPTS_DIR, "build_pan_cancer_expression_matrix.R"),
            heme_arg = _PAN_EXPR_HEME_ARG,
            metadata_arg = _PAN_EXPR_METADATA_ARG,
            output_meta_arg = _PAN_EXPR_OUTPUT_META_ARG,
            expr_dir = os.path.abspath(os.path.dirname(PAN_EXPR_CFG.get("output_rds", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds")))
        log: os.path.join(LOGROOT, "build_pan_cancer_expression_matrix.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.expr_dir}"
            Rscript "{params.script}" \
              --brca-vst "{input.brca}" \
              --nbl-vst "{input.nbl}" \
              --rbl-vst "{input.rbl}" \
              {params.heme_arg} \
              {params.metadata_arg} \
              --genes "{input.genes}" \
              --output "{output.expr_rds}" \
              {params.output_meta_arg} \
              > "{log}" 2>&1
            '''

    # Rule: build_pan_cancer_feature_expression_matrix_cell_lines_only
    # Method role: expression-construction rule for cell-line-only pan-cancer feature-space analyses.
    # Flow: configured cell-line VST inputs plus feature table -> cell-line-only feature matrix.
    # Provides input for cell-line-only similarity network and community sensitivity analyses.
    rule build_pan_cancer_feature_expression_matrix_cell_lines_only:
        """Build a feature-limited pan-cancer matrix containing cell lines only."""
        input:
            genes = PAN_FEATURES_CLEAN,
            brca = _PAN_EXPR_BRCA_VST,
            nbl  = _PAN_EXPR_NBL_VST,
            rbl  = _PAN_EXPR_RBL_VST,
            heme = _PAN_EXPR_CELL_LINES_HEME_INPUT
        output:
            expr_rds = PAN_EXPR_CELL_LINES_RDS,
            meta_tsv = PAN_EXPR_CELL_LINES_META
        params:
            script = os.path.join(SCRIPTS_DIR, "build_pan_cancer_expression_matrix.R"),
            heme_arg = _PAN_EXPR_CELL_LINES_HEME_ARG,
            metadata_arg = _PAN_EXPR_METADATA_ARG,
            output_meta_arg = _PAN_EXPR_CELL_LINES_OUTPUT_META_ARG,
            expr_dir = os.path.abspath(os.path.dirname(PAN_EXPR_CELL_LINES_RDS))
        log: os.path.join(LOGROOT, "build_pan_cancer_expression_matrix_cell_lines_only.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.expr_dir}"
            Rscript "{params.script}" \
              --brca-vst "{input.brca}" \
              --nbl-vst "{input.nbl}" \
              --rbl-vst "{input.rbl}" \
              {params.heme_arg} \
              {params.metadata_arg} \
              --genes "{input.genes}" \
              --sample-type-filter cell_line \
              --output "{output.expr_rds}" \
              {params.output_meta_arg} \
              > "{log}" 2>&1
            '''

    PAN_ALIGNMENT_UMAP_CFG = MARKER_POST_CFG.get("tumour_cell_line_alignment_umap", {})
    PAN_ALIGNMENT_UMAP_OUTDIR = PAN_ALIGNMENT_UMAP_CFG.get(
        "output_dir",
        "results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap"
    )
    PAN_ALIGNMENT_UMAP_META = PAN_ALIGNMENT_UMAP_CFG.get(
        "metadata_tsv",
        PAN_EXPR_META
    )
    _PAN_ALIGNMENT_UMAP_METRICS_RAW = PAN_ALIGNMENT_UMAP_CFG.get("metrics", ["cosine", "euclidean"])
    if isinstance(_PAN_ALIGNMENT_UMAP_METRICS_RAW, str):
        PAN_ALIGNMENT_UMAP_METRICS = [
            m.strip() for m in _PAN_ALIGNMENT_UMAP_METRICS_RAW.split(",") if m.strip()
        ]
    else:
        PAN_ALIGNMENT_UMAP_METRICS = [str(m).strip() for m in _PAN_ALIGNMENT_UMAP_METRICS_RAW if str(m).strip()]

    PAN_ALIGNMENT_UMAP_SOURCE_META = PAN_ALIGNMENT_UMAP_CFG.get(
        "source_metadata_tsv",
        "results/unsupervised/multicohort_cancer/inputs/joint_metadata.tsv"
    )
    PAN_ALIGNMENT_UMAP_FIGURE_WIDTH = float(PAN_ALIGNMENT_UMAP_CFG.get("figure_width", 9.0))
    PAN_ALIGNMENT_UMAP_FIGURE_HEIGHT = float(PAN_ALIGNMENT_UMAP_CFG.get("figure_height", 7.0))
    PAN_ALIGNMENT_UMAP_SLIDE_WIDTH = float(PAN_ALIGNMENT_UMAP_CFG.get("slide_width", 13.33))
    PAN_ALIGNMENT_UMAP_SLIDE_HEIGHT = float(PAN_ALIGNMENT_UMAP_CFG.get("slide_height", 7.5))
    PAN_ALIGNMENT_UMAP_LEGEND_POSITION = PAN_ALIGNMENT_UMAP_CFG.get("legend_position", "bottom")
    PAN_ALIGNMENT_UMAP_SOURCE_CANCER_SLIDE_THRESHOLD = int(PAN_ALIGNMENT_UMAP_CFG.get("source_cancer_slide_threshold", 6))

    def alignment_umap_metric_label(metric):
        return re.sub(r"[^A-Za-z0-9]+", "_", metric)

    def alignment_umap_stem(out_stem, feature_label, metric, view=None):
        parts = [out_stem, feature_label]
        if view:
            parts.append(view)
        parts.append(alignment_umap_metric_label(metric))
        return "_".join(parts)

    def alignment_umap_outputs(outdir, prefix, out_stem, feature_label, metrics, ext, view=None):
        return [
            os.path.join(
                outdir,
                f"{prefix}_{alignment_umap_stem(out_stem, feature_label, m, view)}.{ext}"
            )
            for m in metrics
        ]

    def pan_alignment_umap_stem(metric, view=None):
        return alignment_umap_stem(
            "pan_cancer_tumour_cell_line_alignment_umap",
            "DEG_SET",
            metric,
            view
        )

    PAN_ALIGNMENT_UMAP_PDFS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m)}.pdf")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SVGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m)}.svg")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_PNGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m)}.png")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_COORDS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"coords_{pan_alignment_umap_stem(m)}.tsv")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_PDFS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE')}.pdf")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_SVGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE')}.svg")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_PNGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE')}.png")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PDFS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE_CANCER')}.pdf")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_CANCER_SVGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE_CANCER')}.svg")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PNGS = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"Fig_{pan_alignment_umap_stem(m, 'SOURCE_CANCER')}.png")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SOURCE_COMPOSITION = [
        os.path.join(PAN_ALIGNMENT_UMAP_OUTDIR, f"source_composition_{pan_alignment_umap_stem(m)}.tsv")
        for m in PAN_ALIGNMENT_UMAP_METRICS
    ]
    PAN_ALIGNMENT_UMAP_SUMMARY = os.path.join(
        PAN_ALIGNMENT_UMAP_OUTDIR,
        "summary_pan_cancer_tumour_cell_line_alignment_umap.tsv"
    )

    # Rule: plot_pan_cancer_tumour_cell_line_alignment_umap
    # Method role: plotting rule for tumour/cell-line alignment in the pan-cancer feature space.
    # Flow: feature-space expression matrix and metadata -> UMAP coordinates, composition tables, and figures.
    # Purpose: visualises alignment; it does not define mapping or ranking outputs.
    rule plot_pan_cancer_tumour_cell_line_alignment_umap:
        """Generate pan-cancer tumour/cell-line alignment UMAPs."""
        input:
            expr_rds = PAN_EXPR_RDS,
            meta_tsv = PAN_ALIGNMENT_UMAP_META,
            source_meta_tsv = PAN_ALIGNMENT_UMAP_SOURCE_META,
            genes = PAN_FEATURES_CLEAN,
            script = os.path.join(SCRIPTS_DIR, "plot_pan_cancer_tumour_cell_line_alignment_umap.R")
        output:
            pdf = PAN_ALIGNMENT_UMAP_PDFS,
            svg = PAN_ALIGNMENT_UMAP_SVGS,
            png = PAN_ALIGNMENT_UMAP_PNGS,
            source_pdf = PAN_ALIGNMENT_UMAP_SOURCE_PDFS,
            source_svg = PAN_ALIGNMENT_UMAP_SOURCE_SVGS,
            source_png = PAN_ALIGNMENT_UMAP_SOURCE_PNGS,
            source_cancer_pdf = PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PDFS,
            source_cancer_svg = PAN_ALIGNMENT_UMAP_SOURCE_CANCER_SVGS,
            source_cancer_png = PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PNGS,
            coords = PAN_ALIGNMENT_UMAP_COORDS,
            source_composition = PAN_ALIGNMENT_UMAP_SOURCE_COMPOSITION,
            summary = PAN_ALIGNMENT_UMAP_SUMMARY
        params:
            outdir = PAN_ALIGNMENT_UMAP_OUTDIR,
            metrics = ",".join(PAN_ALIGNMENT_UMAP_METRICS),
            page = PAN_ALIGNMENT_UMAP_CFG.get("page", "deg_set_alignment"),
            figure_width = PAN_ALIGNMENT_UMAP_FIGURE_WIDTH,
            figure_height = PAN_ALIGNMENT_UMAP_FIGURE_HEIGHT,
            slide_width = PAN_ALIGNMENT_UMAP_SLIDE_WIDTH,
            slide_height = PAN_ALIGNMENT_UMAP_SLIDE_HEIGHT,
            legend_position = PAN_ALIGNMENT_UMAP_LEGEND_POSITION,
            source_cancer_slide_threshold = PAN_ALIGNMENT_UMAP_SOURCE_CANCER_SLIDE_THRESHOLD,
            seed = PIPELINE_SEED
        threads: 8
        log: os.path.join(LOGROOT, "plot_pan_cancer_tumour_cell_line_alignment_umap.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{input.script}" \
              --pipe_root "{BASE}" \
              --expr_rds "{input.expr_rds}" \
              --meta_tsv "{input.meta_tsv}" \
              --source_meta_tsv "{input.source_meta_tsv}" \
              --feature_list "{input.genes}" \
              --outdir "{params.outdir}" \
              --dist_metrics "{params.metrics}" \
              --page "{params.page}" \
              --width "{params.figure_width}" \
              --height "{params.figure_height}" \
              --slide_width "{params.slide_width}" \
              --slide_height "{params.slide_height}" \
              --legend_position "{params.legend_position}" \
              --source_cancer_slide_threshold "{params.source_cancer_slide_threshold}" \
              --seed {params.seed} \
              > "{log}" 2>&1
            for f in {output.pdf} {output.svg} {output.png} {output.source_pdf} {output.source_svg} {output.source_png} {output.source_cancer_pdf} {output.source_cancer_svg} {output.source_cancer_png} {output.coords} {output.source_composition} "{output.summary}"; do
              test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
            done
            '''

    # ==========================================================================
    # OPTIONAL: raw all-gene tumour/cell-line UMAP
    # ==========================================================================
    # This rule is activated only when a true raw expression matrix is configured.
    # It intentionally does not fall back to the VST all-gene matrix, so the raw,
    # VST-all-gene, and pan-cancer feature-set spaces remain distinct.
    RAW_ALLGENE_UMAP_EXPR_CFG = (
        PAN_ALIGNMENT_UMAP_CFG.get("raw_all_gene_expr_rds")
        or PAN_ALIGNMENT_UMAP_CFG.get("raw_expr_rds")
        or ""
    )
    RAW_ALLGENE_UMAP_TARGETS = []
    if RAW_ALLGENE_UMAP_EXPR_CFG:
        RAW_ALLGENE_UMAP_OUTDIR = PAN_ALIGNMENT_UMAP_CFG.get(
            "raw_output_dir",
            "results/unsupervised/pan_cancer/raw_all_gene_tumour_cell_line_umap"
        )
        RAW_ALLGENE_UMAP_EXPR_RDS = abspath(RAW_ALLGENE_UMAP_EXPR_CFG)
        RAW_ALLGENE_UMAP_META_TSV = abspath(PAN_ALIGNMENT_UMAP_CFG.get(
            "raw_metadata_tsv",
            PAN_ALIGNMENT_UMAP_META
        ))
        RAW_ALLGENE_UMAP_SOURCE_META = abspath(PAN_ALIGNMENT_UMAP_CFG.get(
            "raw_source_metadata_tsv",
            PAN_ALIGNMENT_UMAP_SOURCE_META
        ))
        RAW_ALLGENE_UMAP_METRICS = PAN_ALIGNMENT_UMAP_METRICS
        RAW_ALLGENE_UMAP_LABEL = "RAW_ALL_GENES"
        RAW_ALLGENE_UMAP_STEM = "tumour_cell_line_alignment_umap"
        RAW_ALLGENE_UMAP_PDFS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "pdf"
        )
        RAW_ALLGENE_UMAP_SVGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "svg"
        )
        RAW_ALLGENE_UMAP_PNGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "png"
        )
        RAW_ALLGENE_UMAP_COORDS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "coords", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "tsv"
        )
        RAW_ALLGENE_UMAP_SOURCE_PDFS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "pdf", "SOURCE"
        )
        RAW_ALLGENE_UMAP_SOURCE_SVGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "svg", "SOURCE"
        )
        RAW_ALLGENE_UMAP_SOURCE_PNGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "png", "SOURCE"
        )
        RAW_ALLGENE_UMAP_SOURCE_CANCER_PDFS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "pdf", "SOURCE_CANCER"
        )
        RAW_ALLGENE_UMAP_SOURCE_CANCER_SVGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "svg", "SOURCE_CANCER"
        )
        RAW_ALLGENE_UMAP_SOURCE_CANCER_PNGS = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "Fig", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "png", "SOURCE_CANCER"
        )
        RAW_ALLGENE_UMAP_SOURCE_COMPOSITION = alignment_umap_outputs(
            RAW_ALLGENE_UMAP_OUTDIR, "source_composition", RAW_ALLGENE_UMAP_STEM,
            RAW_ALLGENE_UMAP_LABEL, RAW_ALLGENE_UMAP_METRICS, "tsv"
        )
        RAW_ALLGENE_UMAP_SUMMARY = os.path.join(
            RAW_ALLGENE_UMAP_OUTDIR,
            f"summary_{RAW_ALLGENE_UMAP_STEM}_{RAW_ALLGENE_UMAP_LABEL}.tsv"
        )
        RAW_ALLGENE_UMAP_TARGETS = (
            RAW_ALLGENE_UMAP_PDFS + RAW_ALLGENE_UMAP_SVGS + RAW_ALLGENE_UMAP_PNGS +
            RAW_ALLGENE_UMAP_SOURCE_PDFS + RAW_ALLGENE_UMAP_SOURCE_SVGS + RAW_ALLGENE_UMAP_SOURCE_PNGS +
            RAW_ALLGENE_UMAP_SOURCE_CANCER_PDFS + RAW_ALLGENE_UMAP_SOURCE_CANCER_SVGS + RAW_ALLGENE_UMAP_SOURCE_CANCER_PNGS +
            RAW_ALLGENE_UMAP_COORDS + RAW_ALLGENE_UMAP_SOURCE_COMPOSITION + [RAW_ALLGENE_UMAP_SUMMARY]
        )

        # Rule: plot_raw_all_gene_tumour_cell_line_alignment_umap
        # Method role: optional plotting rule for raw all-gene tumour/cell-line alignment.
        # Flow: configured raw all-gene matrix -> UMAP coordinates and source-composition tables.
        # Used by: context visualisation only; it is gated by configured raw inputs.
        rule plot_raw_all_gene_tumour_cell_line_alignment_umap:
            """Optional raw-expression all-gene tumour/cell-line alignment UMAPs."""
            input:
                expr_rds = RAW_ALLGENE_UMAP_EXPR_RDS,
                meta_tsv = RAW_ALLGENE_UMAP_META_TSV,
                source_meta_tsv = RAW_ALLGENE_UMAP_SOURCE_META,
                script = os.path.join(SCRIPTS_DIR,
                    "plot_pan_cancer_tumour_cell_line_alignment_umap.R")
            output:
                pdf = RAW_ALLGENE_UMAP_PDFS,
                svg = RAW_ALLGENE_UMAP_SVGS,
                png = RAW_ALLGENE_UMAP_PNGS,
                source_pdf = RAW_ALLGENE_UMAP_SOURCE_PDFS,
                source_svg = RAW_ALLGENE_UMAP_SOURCE_SVGS,
                source_png = RAW_ALLGENE_UMAP_SOURCE_PNGS,
                source_cancer_pdf = RAW_ALLGENE_UMAP_SOURCE_CANCER_PDFS,
                source_cancer_svg = RAW_ALLGENE_UMAP_SOURCE_CANCER_SVGS,
                source_cancer_png = RAW_ALLGENE_UMAP_SOURCE_CANCER_PNGS,
                coords = RAW_ALLGENE_UMAP_COORDS,
                source_composition = RAW_ALLGENE_UMAP_SOURCE_COMPOSITION,
                summary = RAW_ALLGENE_UMAP_SUMMARY
            params:
                outdir = RAW_ALLGENE_UMAP_OUTDIR,
                metrics = ",".join(RAW_ALLGENE_UMAP_METRICS),
                feature_label = RAW_ALLGENE_UMAP_LABEL,
                out_stem = RAW_ALLGENE_UMAP_STEM,
                summary_basename = os.path.basename(RAW_ALLGENE_UMAP_SUMMARY),
                page = PAN_ALIGNMENT_UMAP_CFG.get("raw_page", "all_gene_alignment"),
                seed = PIPELINE_SEED
            threads: 8
            log: os.path.join(LOGROOT, "plot_raw_all_gene_tumour_cell_line_alignment_umap.log")
            conda: CONDA_ENV_R
            shell:
                r'''
                mkdir -p "{params.outdir}"
                export TMPDIR="{params.outdir}/.umap_tmp"
                mkdir -p "$TMPDIR"
                Rscript "{input.script}" \
                  --pipe_root "{BASE}" \
                  --expr_rds "{input.expr_rds}" \
                  --meta_tsv "{input.meta_tsv}" \
                  --source_meta_tsv "{input.source_meta_tsv}" \
                  --feature_mode all_genes \
                  --feature_label "{params.feature_label}" \
                  --out_stem "{params.out_stem}" \
                  --summary_basename "{params.summary_basename}" \
                  --outdir "{params.outdir}" \
                  --dist_metrics "{params.metrics}" \
                  --page "{params.page}" \
                  --seed {params.seed} \
                  > "{log}" 2>&1
                rm -rf "$TMPDIR"
                for f in {output.pdf} {output.svg} {output.png} {output.source_pdf} {output.source_svg} {output.source_png} {output.source_cancer_pdf} {output.source_cancer_svg} {output.source_cancer_png} {output.coords} {output.source_composition} "{output.summary}"; do
                  test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
                done
                '''

    # ==========================================================================
    # VST all-gene tumour/cell-line UMAP (BRCA + NBL + RBL)
    # ==========================================================================
    # Independent of the curated pan-cancer DEG / feature-set workflow above.
    # This unsupervised runs UMAP on the full VST gene space (after
    # NA + zero-variance filtering on the merged matrix), with no class label
    # used during fitting and no pan-cancer feature subsetting. Filenames use
    # the VST_ALL_GENE label to keep this run's outputs cleanly separable from
    # the DEG-set alignment outputs above.
    VST_ALLGENE_UMAP_OUTDIR = "results/unsupervised/pan_cancer/vst_all_gene_tumour_cell_line_umap"
    VST_ALLGENE_UMAP_EXPR_RDS = os.path.join(
        VST_ALLGENE_UMAP_OUTDIR,
        "vst_all_gene_brca_nbl_rbl_expr.rds"
    )
    VST_ALLGENE_UMAP_META_TSV = os.path.join(
        VST_ALLGENE_UMAP_OUTDIR,
        "vst_all_gene_metadata.tsv"
    )
    VST_ALLGENE_UMAP_MERGE_LOG = os.path.join(
        VST_ALLGENE_UMAP_OUTDIR,
        "vst_all_gene_metadata.tsv.merge_log.tsv"
    )
    VST_ALLGENE_UMAP_METRICS = ["cosine", "euclidean"]
    VST_ALLGENE_UMAP_LABEL   = "VST_ALL_GENE"
    VST_ALLGENE_UMAP_STEM    = "tumour_cell_line_alignment_umap"

    VST_ALLGENE_UMAP_PDFS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m)}.pdf")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SVGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m)}.svg")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_PNGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m)}.png")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_COORDS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"coords_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m)}.tsv")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_PDFS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE')}.pdf")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_SVGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE')}.svg")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_PNGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE')}.png")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_CANCER_PDFS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE_CANCER')}.pdf")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_CANCER_SVGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE_CANCER')}.svg")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_CANCER_PNGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"Fig_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m, 'SOURCE_CANCER')}.png")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SOURCE_COMPOSITION = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, f"source_composition_{alignment_umap_stem(VST_ALLGENE_UMAP_STEM, VST_ALLGENE_UMAP_LABEL, m)}.tsv")
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SUMMARY = os.path.join(
        VST_ALLGENE_UMAP_OUTDIR,
        f"summary_{VST_ALLGENE_UMAP_STEM}_{VST_ALLGENE_UMAP_LABEL}.tsv"
    )

    # Direct, profile-independent paths to the three per-cohort joint VST RDS
    # files (samples x genes after orientation handling in the merge script).
    VST_ALLGENE_BRCA_VST = profile_vst_joint_abs("brca")
    VST_ALLGENE_NBL_VST  = profile_vst_joint_abs("nbl")
    VST_ALLGENE_RBL_VST  = profile_vst_joint_abs("rbl")

    # Joint metadata produced by the multicohort_cancer profile (sample_id,
    # cancer_type, sample_type[, cohort]). Used here purely to annotate
    # samples AFTER UMAP fitting.
    VST_ALLGENE_JOINT_META = os.path.join(
        BASE, "results", "unsupervised", "multicohort_cancer", "inputs",
        "joint_metadata.tsv"
    )

    # Rule: merge_vst_all_genes_brca_nbl_rbl
    # Method role: data-preparation rule that merges cohort VST matrices without feature subsetting.
    # Flow: BRCA/NBL/RBL VST matrices and multicohort metadata -> all-gene VST expression matrix.
    # Provides the VST all-gene UMAP.
    rule merge_vst_all_genes_brca_nbl_rbl:
        """
        Merge BRCA + NBL + RBL joint VST tumour/cell-line matrices on the
        intersection of Ensembl gene IDs (versions stripped), keeping the
        first occurrence of duplicated IDs, then drop NA-bearing and
        zero-variance genes on the merged matrix. No class label is used
        during the merge or filtering.
        """
        input:
            brca = VST_ALLGENE_BRCA_VST,
            nbl  = VST_ALLGENE_NBL_VST,
            rbl  = VST_ALLGENE_RBL_VST,
            joint_meta = VST_ALLGENE_JOINT_META,
            script = os.path.join(SCRIPTS_DIR, "merge_vst_all_genes_brca_nbl_rbl.R")
        output:
            expr_rds  = VST_ALLGENE_UMAP_EXPR_RDS,
            meta_tsv  = VST_ALLGENE_UMAP_META_TSV,
            merge_log = VST_ALLGENE_UMAP_MERGE_LOG
        log: os.path.join(LOGROOT, "merge_vst_all_genes_brca_nbl_rbl.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "$(dirname "{output.expr_rds}")"
            Rscript "{input.script}" \
              --brca_rds  "{input.brca}" \
              --nbl_rds   "{input.nbl}" \
              --rbl_rds   "{input.rbl}" \
              --joint_meta "{input.joint_meta}" \
              --out_expr  "{output.expr_rds}" \
              --out_meta  "{output.meta_tsv}" \
              --out_log   "{output.merge_log}" \
              > "{log}" 2>&1
            for f in "{output.expr_rds}" "{output.meta_tsv}" "{output.merge_log}"; do
              test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
            done
            '''

    # Rule: plot_vst_all_gene_tumour_cell_line_alignment_umap
    # Method role: plotting rule for all-gene VST tumour/cell-line alignment.
    # Flow: merged all-gene VST matrix -> UMAP coordinates and source-composition tables.
    # Used by: context visualisation against the marker-derived feature space.
    rule plot_vst_all_gene_tumour_cell_line_alignment_umap:
        """
        Unsupervised UMAP on the BRCA+NBL+RBL VST all-gene merged
        matrix. Uses the same plotting script as the pan-cancer DEG-set
        alignment figure, invoked with --feature_mode=all_genes and the
        VST_ALL_GENE label so output filenames remain disjoint from the
        DEG-set alignment outputs.
        """
        input:
            expr_rds = VST_ALLGENE_UMAP_EXPR_RDS,
            meta_tsv = VST_ALLGENE_UMAP_META_TSV,
            source_meta_tsv = VST_ALLGENE_JOINT_META,
            script   = os.path.join(SCRIPTS_DIR,
                "plot_pan_cancer_tumour_cell_line_alignment_umap.R")
        output:
            pdf     = VST_ALLGENE_UMAP_PDFS,
            svg     = VST_ALLGENE_UMAP_SVGS,
            png     = VST_ALLGENE_UMAP_PNGS,
            source_pdf = VST_ALLGENE_UMAP_SOURCE_PDFS,
            source_svg = VST_ALLGENE_UMAP_SOURCE_SVGS,
            source_png = VST_ALLGENE_UMAP_SOURCE_PNGS,
            source_cancer_pdf = VST_ALLGENE_UMAP_SOURCE_CANCER_PDFS,
            source_cancer_svg = VST_ALLGENE_UMAP_SOURCE_CANCER_SVGS,
            source_cancer_png = VST_ALLGENE_UMAP_SOURCE_CANCER_PNGS,
            coords  = VST_ALLGENE_UMAP_COORDS,
            source_composition = VST_ALLGENE_UMAP_SOURCE_COMPOSITION,
            summary = VST_ALLGENE_UMAP_SUMMARY
        params:
            outdir           = VST_ALLGENE_UMAP_OUTDIR,
            metrics          = ",".join(VST_ALLGENE_UMAP_METRICS),
            feature_label    = VST_ALLGENE_UMAP_LABEL,
            out_stem         = VST_ALLGENE_UMAP_STEM,
            summary_basename = os.path.basename(VST_ALLGENE_UMAP_SUMMARY),
            page             = "all_gene_alignment",
            seed             = PIPELINE_SEED
        threads: 8
        log: os.path.join(LOGROOT, "plot_vst_all_gene_tumour_cell_line_alignment_umap.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            # uwot's Annoy backend writes its NN index to TMPDIR; the cluster
            # node's /tmp is small, so redirect to a per-rule scratch dir on the
            # work filesystem before invoking Rscript.
            export TMPDIR="{params.outdir}/.umap_tmp"
            mkdir -p "$TMPDIR"
            Rscript "{input.script}" \
              --pipe_root "{BASE}" \
              --expr_rds "{input.expr_rds}" \
              --meta_tsv "{input.meta_tsv}" \
              --source_meta_tsv "{input.source_meta_tsv}" \
              --feature_mode all_genes \
              --feature_label "{params.feature_label}" \
              --out_stem "{params.out_stem}" \
              --summary_basename "{params.summary_basename}" \
              --outdir "{params.outdir}" \
              --dist_metrics "{params.metrics}" \
              --page "{params.page}" \
              --seed {params.seed} \
              > "{log}" 2>&1
            rm -rf "$TMPDIR"
            for f in {output.pdf} {output.svg} {output.png} {output.source_pdf} {output.source_svg} {output.source_png} {output.source_cancer_pdf} {output.source_cancer_svg} {output.source_cancer_png} {output.coords} {output.source_composition} "{output.summary}"; do
              test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
            done
            '''

    MAPPING_OUTDIR = MAPPING_CFG.get("output_dir", "results/unsupervised/pan_cancer/tumour_mapping")

    def mapping_metrics_summary():
        return os.path.join(MAPPING_OUTDIR, "metrics_summary_group_level.tsv")

    # Precompute optional CLI fragment for gene-symbol-map.
    _MAPPING_GENE_MAP_RAW = MAPPING_CFG.get("gene_symbol_map", "")
    _MAPPING_GENE_MAP_ARG = f'--gene-symbol-map "{_MAPPING_GENE_MAP_RAW}"' if _MAPPING_GENE_MAP_RAW else ""

    # Rule: score_tumour_cellline_mapping
    # Method role: mapping/ranking rule that scores tumour-to-cell-line similarity in the pan-cancer feature space.
    # Flow: feature-space expression matrix and metadata -> mapping metrics and ranking tables.
    # Provides input for pan-cancer ranking evaluations and ECDF plots.
    rule score_tumour_cellline_mapping:
        """Score tumour-to-cell-line mappings using the combined expression object."""
        input:
            expr_rds = PAN_EXPR_RDS,
            genes = PAN_FEATURES_CLEAN
        output:
            metrics_summary = mapping_metrics_summary(),
            tumour_rankings = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_to_cellline_group_rankings.tsv"),
            tumour_summary = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_mapping_group_summary.tsv"),
            tumour_scores = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_group_scores.rds"),
            tumour_scores_long = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_group_scores_long.tsv.gz"),
            tumour_metrics = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "metrics_summary_group_level.tsv"),
            cellline_rankings = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_to_tumour_rankings.tsv"),
            cellline_summary = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_mapping_summary.tsv"),
            cellline_metrics = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_metrics_overall.tsv"),
            cellline_metrics_summary = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "metrics_summary.tsv"),
            cellline_scores = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_tumour_scores.rds"),
            cellline_scores_long = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_tumour_scores_long.tsv.gz")
        params:
            script=os.path.join(SCRIPTS_DIR, "score_tumour_cellline_mapping.R"),
            outdir=MAPPING_OUTDIR,
            top_k=MAPPING_CFG.get("top_k", 10),
            similarity=MAPPING_CFG.get("similarity", "spearman"),
            zscore_value="TRUE" if MAPPING_CFG.get("zscore_cohort", True) else "FALSE",
            drop_global_value="TRUE" if MAPPING_CFG.get("drop_global_genes", True) else "FALSE",
            global_regex=MAPPING_CFG.get("global_gene_regex", "^RPL|^RPS|^MT-|^HIST"),
            gene_map_arg=_MAPPING_GENE_MAP_ARG
        log: os.path.join(LOGROOT, "score_tumour_cellline_mapping.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --pan-cancer-expr "{input.expr_rds}" \
              --genes "{input.genes}" \
              --output-dir "{params.outdir}" \
              --top-k {params.top_k} \
              --similarity {params.similarity} \
              --global-gene-regex "{params.global_regex}" \
              --zscore-cohort {params.zscore_value} \
              --drop-global-genes {params.drop_global_value} \
              {params.gene_map_arg} \
              > "{log}" 2>&1
            test -s "{output.metrics_summary}" || (echo "ERROR: missing {output.metrics_summary}" >&2; exit 1)
            '''


    PAN_CANCER_DIR = "results/unsupervised/pan_cancer"
    PAN_GRAPH_DIR = os.path.join(PAN_CANCER_DIR, "graph")
    PAN_FIG_DIR = os.path.join(PAN_CANCER_DIR, "figures")
    PAN_COR_RDS = os.path.join(PAN_CANCER_DIR, "pan_cancer_cor.rds")

    # Rule: compute_pan_cancer_correlation
    # Method role: similarity analysis rule computing cell-line/tumour correlations in the pan-cancer feature space.
    # Flow: feature-space expression matrix -> correlation matrix.
    # Feeds graph construction and tumour-cell-line mapping tables.
    rule compute_pan_cancer_correlation:
        input: expr = PAN_EXPR_RDS
        output: cor = PAN_COR_RDS
        params: script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_correlation.R"), method = "spearman"
        log: os.path.join(LOGROOT, "compute_pan_cancer_correlation.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {PAN_CANCER_DIR} && Rscript {params.script} --input {input.expr} --output {output.cor} --method {params.method} > {log} 2>&1 && test -s {output.cor}"

    # Rule: build_pan_cancer_graph
    # Method role: graph-construction rule that constructs transcriptomic similarity in the pan-cancer feature space.
    # Flow: pan-cancer correlation matrix -> graph edge and node sidecars.
    # Provides input for community analysis, inspection, and graph plotting.
    rule build_pan_cancer_graph:
        input: cor = PAN_COR_RDS
        output:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "build_pan_cancer_graph.R"), outdir = PAN_GRAPH_DIR, k = 20, min_cor = 0.0
        log: os.path.join(LOGROOT, "build_pan_cancer_graph.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --cor-matrix {input.cor} --k {params.k} --min-cor {params.min_cor} --output-dir {params.outdir} > {log} 2>&1 && test -s {output.edges} && test -s {output.components}"

    # Rule: compute_pan_cancer_communities
    # Method role: graph analysis rule assigning communities on the pan-cancer similarity graph.
    # Flow: graph edges and metadata -> community assignments and metrics.
    # Provides annotations for graph plots and ranking evaluations.
    rule compute_pan_cancer_communities:
        input:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            meta = PAN_EXPR_RDS
        output:
            communities = os.path.join(PAN_GRAPH_DIR, "pan_cancer_communities.tsv"),
            summary = os.path.join(PAN_GRAPH_DIR, "community_validation", "community_summary.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_communities.R"), outdir = os.path.join(PAN_GRAPH_DIR, "community_validation"), seed = PIPELINE_SEED
        log: os.path.join(LOGROOT, "compute_pan_cancer_communities.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --edges {input.edges} --meta {input.meta} --out {output.communities} --outdir {params.outdir} --method both --seed {params.seed} > {log} 2>&1 && test -s {output.communities}"

    # Rule: inspect_pan_cancer_graph
    # Method role: audit rule for pan-cancer graph structure and metadata coverage.
    # Flow: graph sidecars and community assignments -> inspection tables.
    # Supports graph QA before final plotting.
    rule inspect_pan_cancer_graph:
        input:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            communities = os.path.join(PAN_GRAPH_DIR, "pan_cancer_communities.tsv"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv")
        output:
            summary = os.path.join(PAN_GRAPH_DIR, "inspection", "INSPECTION_SUMMARY.md"),
            bridges = os.path.join(PAN_GRAPH_DIR, "inspection", "isolated_or_bridge_nodes.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "inspect_pan_cancer_graph.R"), graph_dir = PAN_GRAPH_DIR, inspection_dir = os.path.join(PAN_GRAPH_DIR, "inspection")
        log: os.path.join(LOGROOT, "inspect_pan_cancer_graph.log")
        conda: CONDA_ENV_R
        shell: "Rscript {params.script} --graph-dir {params.graph_dir} --inspection-dir {params.inspection_dir} > {log} 2>&1 && test -s {output.bridges}"

    # Rule: plot_pan_cancer_graph
    # Method role: plotting rule for the pan-cancer feature-space similarity graph.
    # Flow: graph sidecars, communities, and layout inputs -> final pan-cancer graph figure.
    # Purpose: reporting visualisation only.
    rule plot_pan_cancer_graph:
        input:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv"),
            meta = PAN_EXPR_RDS
        output:
            pdf = os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_graph.pdf"),
            components_pdf = os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_graph_components.pdf"),
            size_pdf = os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_component_size.pdf"),
            provenance = os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_graph_provenance.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "plot_pan_cancer_graph.R"), outdir = PAN_FIG_DIR
        log: os.path.join(LOGROOT, "plot_pan_cancer_graph.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --edges {input.edges} --components {input.components} --meta {input.meta} --outdir {params.outdir} --expected-gene-count 0 > {log} 2>&1 && printf 'figure_name\tscript\tcommand\tgit_commit\ttimestamp\tinput_files\toutput_files\tupstream_tables\tkey_parameters\tsoftware_versions\tfigure_type\tsource_pipeline_root\tcopied_to_figure_export_path\tlegacy_source_path\tnotes\nFig_pan_cancer_graph.pdf\tscripts/plot_pan_cancer_graph.R\tRscript scripts/plot_pan_cancer_graph.R\tunavailable_not_git_worktree\tNA\t{input.edges};{input.components};{input.meta}\t{output.pdf};{output.components_pdf};{output.size_pdf}\t{input.edges};{input.components}\texpected_gene_count=infer_from_rds_genes\tR/igraph/ggplot2\tpan_cancer\t{PIPE_ROOT}\t\t\tTranscriptomic similarity network for prioritisation and neighbourhood assignment only\n' > {output.provenance} && test -s {output.pdf}"

    # Rule: cellline_precision_at_k
    # Method role: ranking evaluation computing cell-line-centred same-cancer-type ranking metrics.
    # Flow: mapping/ranking tables -> precision-at-k, rank-metric, and provenance outputs.
    # Supports cautious ranking-based agreement assessment.
    rule cellline_precision_at_k:
        """
        Cell-line-centred tumour ranking analysis with standalone retained plots.

        Replicate profiles are collapsed to biological cell-line groups by
        arithmetic mean Spearman correlation before ranking. Retained plots are
        standalone PDF/PNG outputs with enlarged fonts; confidence-margin and
        top-10 score-distributions are not retained plot targets.
        """
        input:
            c2t_long = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity",
                                    "cellline_tumour_scores_long.tsv.gz"),
            t2c_long = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity",
                                    "tumour_cellline_group_scores_long.tsv.gz"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv")
        output:
            top1_pdf = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top1_lineage_agreement.pdf"),
            top1_png = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top1_lineage_agreement.png"),
            precision_pdf = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_precision_at_k.pdf"),
            precision_png = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_precision_at_k.png"),
            rank_percentile_pdf = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_same_lineage_rank_percentile.pdf"),
            rank_percentile_png = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_same_lineage_rank_percentile.png"),
            top50_lineage_pdf = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top50_lineage_composition.pdf"),
            top50_lineage_png = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top50_lineage_composition.png"),
            top50_component_pdf = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top50_component_composition.pdf"),
            top50_component_png = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_to_tumour_top50_component_composition.png"),
            top1_accuracy = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_top1_accuracy.tsv"),
            topk_metrics  = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_topk_metrics.tsv"),
            topk_lineage_composition = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_topk_lineage_composition.tsv"),
            component_summary = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_component_mapping_summary.tsv"),
            rank_summary  = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_centred_rank_summary.tsv"),
            replicate_mapping = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "replicate_collapse_mapping.tsv"),
            reciprocal_summary = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "reciprocal_mapping_summary.tsv"),
            reciprocal_per_cellline = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "reciprocal_per_cell_line.tsv"),
            qc_log        = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_centred_QC_log.txt"),
            provenance    = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_similarity_precision_provenance.tsv")
        params:
            script      = os.path.join(SCRIPTS_DIR,
                              "cellline_centred_precision_at_k.R"),
            outdir      = os.path.join(MAPPING_OUTDIR,
                              "cellline_similarity_precision_bootstrap"),
            mapping_dir = MAPPING_OUTDIR
        log: os.path.join(LOGROOT, "cellline_precision_at_k.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --mapping-dir "{params.mapping_dir}" \
              --outdir "{params.outdir}" \
              --tumour-components "{input.components}" \
              > "{log}" 2>&1
            test -s "{output.top1_pdf}" || (echo "ERROR: missing {output.top1_pdf}" >&2; exit 1)
            test -s "{output.top1_png}" || (echo "ERROR: missing {output.top1_png}" >&2; exit 1)
            test -s "{output.precision_pdf}" || (echo "ERROR: missing {output.precision_pdf}" >&2; exit 1)
            test -s "{output.precision_png}" || (echo "ERROR: missing {output.precision_png}" >&2; exit 1)
            test -s "{output.rank_percentile_pdf}" || (echo "ERROR: missing {output.rank_percentile_pdf}" >&2; exit 1)
            test -s "{output.rank_percentile_png}" || (echo "ERROR: missing {output.rank_percentile_png}" >&2; exit 1)
            test -s "{output.top50_lineage_pdf}" || (echo "ERROR: missing {output.top50_lineage_pdf}" >&2; exit 1)
            test -s "{output.top50_lineage_png}" || (echo "ERROR: missing {output.top50_lineage_png}" >&2; exit 1)
            test -s "{output.top50_component_pdf}" || (echo "ERROR: missing {output.top50_component_pdf}" >&2; exit 1)
            test -s "{output.top50_component_png}" || (echo "ERROR: missing {output.top50_component_png}" >&2; exit 1)
            test -s "{output.topk_metrics}" || (echo "ERROR: missing {output.topk_metrics}" >&2; exit 1)
            test -s "{output.component_summary}" || (echo "ERROR: missing {output.component_summary}" >&2; exit 1)
            # Provenance record
            printf 'field\tvalue\n' > {output.provenance}
            printf 'script\t{params.script}\n' >> {output.provenance}
            printf 'mapping_dir\t{params.mapping_dir}\n' >> {output.provenance}
            printf 'outdir\t{params.outdir}\n' >> {output.provenance}
            printf 'component_input\t{input.components}\n' >> {output.provenance}
            printf 'timestamp\t%s\n' "$(date -Iseconds)" >> {output.provenance}
            printf 'source_pipeline\t{PIPE_ROOT}\n' >> {output.provenance}
            printf 'bootstrap_resamples\t2000\n' >> {output.provenance}
            printf 'bootstrap_statistic\tmean_precision_at_k\n' >> {output.provenance}
            printf 'bootstrap_bounds\t2.5th_97.5th_percentile\n' >> {output.provenance}
            printf 'resampling_unit\tbiological_cell_line_group\n' >> {output.provenance}
            printf 'replicate_collapse\tmean_score_across_replicate_profiles\n' >> {output.provenance}
            printf 'bootstrap_seed\t20260603\n' >> {output.provenance}
            '''


    # Rule: plot_tumour_to_cellline_rank_ecdf_top10_fraction
    # Method role: plotting rule for tumour-to-cell-line rank distributions and top-10 fraction.
    # Flow: tumour-to-cell-line ranking outputs -> ECDF and top-10 fraction figures.
    # Purpose: visualises ranking behaviour without asserting biological equivalence.
    rule plot_tumour_to_cellline_rank_ecdf_top10_fraction:
        input:
            metrics_summary = mapping_metrics_summary(),
            scores = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_group_scores.rds"),
            meta = PAN_EXPR_RDS
        output:
            rank_ecdf_pdf = os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.pdf"),
            rank_ecdf_png = os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.png"),
            top10_fraction_pdf = os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.pdf"),
            top10_fraction_png = os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.png"),
            summary = os.path.join(PAN_FIG_DIR, "ecdf_plots", "model_prioritisation_rank_summary.tsv"),
            label_map = os.path.join(PAN_FIG_DIR, "ecdf_plots", "model_prioritisation_display_label_mapping.tsv"),
            provenance = os.path.join(PAN_FIG_DIR, "ecdf_plots", "tumour_to_cellline_rank_ecdf_top10_fraction_provenance.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "plot_ecdf_rank_combined_pub.R"), outdir = os.path.join(PAN_FIG_DIR, "ecdf_plots"), pipeline_dir = PAN_CANCER_DIR, top_k = MAPPING_CFG.get("top_k", 10), n_show = 4
        log: os.path.join(LOGROOT, "plot_tumour_to_cellline_rank_ecdf_top10_fraction.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --pipeline-dir {params.pipeline_dir} --score-rds {input.scores} --meta-rds {input.meta} --outdir {params.outdir} --n-show {params.n_show} --top-k {params.top_k} > {log} 2>&1 && test -s {output.rank_ecdf_pdf} && test -s {output.rank_ecdf_png} && test -s {output.top10_fraction_pdf} && test -s {output.top10_fraction_png} && test -s {output.summary} && test -s {output.provenance}"


    PAN_RANKING_EVALUATION_DIR = os.path.join(PAN_CANCER_DIR, "ranking", "evaluation")
    PAN_T2C_MRR_AT10_BY_TUMOUR = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "tumour_to_cellline_mrr_at10_by_tumour.tsv",
    )
    PAN_T2C_MRR_AT10_BY_COHORT = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "tumour_to_cellline_mrr_at10_by_cohort.tsv",
    )
    PAN_T2C_MRR_AT10_PDF = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_tumour_to_cellline_mrr_at10_distribution.pdf",
    )
    PAN_T2C_MRR_AT10_PNG = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_tumour_to_cellline_mrr_at10_distribution.png",
    )
    PAN_C2T_COMPONENT_SAME_TABLE = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "cellline_to_tumour_same_cancer_top50_component_composition.tsv",
    )
    PAN_C2T_COMPONENT_ALL_TABLE = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "cellline_to_tumour_all_top50_component_composition.tsv",
    )
    PAN_C2T_COMPONENT_SUMMARY_TABLE = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "cellline_to_tumour_component_composition_summary.tsv",
    )
    PAN_C2T_COMPONENT_SAME_PDF = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_cellline_to_tumour_same_cancer_top50_component_composition.pdf",
    )
    PAN_C2T_COMPONENT_SAME_PNG = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_cellline_to_tumour_same_cancer_top50_component_composition.png",
    )
    PAN_C2T_COMPONENT_ALL_PDF = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_cellline_to_tumour_all_top50_component_composition.pdf",
    )
    PAN_C2T_COMPONENT_ALL_PNG = os.path.join(
        PAN_RANKING_EVALUATION_DIR,
        "Fig_cellline_to_tumour_all_top50_component_composition.png",
    )

    # Rule: plot_tumour_to_cellline_mrr_at10_distribution
    # Method role: plotting rule for tumour-to-cell-line MRR@10 metrics.
    # Flow: ranking metrics -> cohort-level MRR@10 tables and figures.
    # Used by: supplements ranking-based agreement assessment.
    rule plot_tumour_to_cellline_mrr_at10_distribution:
        """Build tumour-level RR@10 tables and a report-ready MRR@10 distribution plot."""
        input:
            script = os.path.join(SCRIPTS_DIR, "plot_tumour_to_cellline_mrr_at10_distribution.R"),
            tumour_rankings = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_to_cellline_group_rankings.tsv"),
            tumour_metrics = os.path.join(MAPPING_OUTDIR, "metrics_summary_group_level.tsv")
        output:
            by_tumour = PAN_T2C_MRR_AT10_BY_TUMOUR,
            by_cohort = PAN_T2C_MRR_AT10_BY_COHORT,
            pdf = PAN_T2C_MRR_AT10_PDF,
            png = PAN_T2C_MRR_AT10_PNG
        params:
            outdir = PAN_RANKING_EVALUATION_DIR,
            top_k = MAPPING_CFG.get("top_k", 10),
            expected_groups = 56
        log: os.path.join(LOGROOT, "plot_tumour_to_cellline_mrr_at10_distribution.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{input.script}" \
              --rankings "{input.tumour_rankings}" \
              --metrics "{input.tumour_metrics}" \
              --by-tumour "{output.by_tumour}" \
              --by-cohort "{output.by_cohort}" \
              --pdf "{output.pdf}" \
              --png "{output.png}" \
              --top-k {params.top_k} \
              --expected-groups {params.expected_groups} \
              > "{log}" 2>&1
            test -s "{output.by_tumour}" || (echo "ERROR: missing {output.by_tumour}" >&2; exit 1)
            test -s "{output.by_cohort}" || (echo "ERROR: missing {output.by_cohort}" >&2; exit 1)
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.png}" || (echo "ERROR: missing {output.png}" >&2; exit 1)
            '''

    # Rule: plot_cellline_to_tumour_same_cancer_component_composition
    # Method role: plotting rule for same-cancer component composition in cell-line-to-tumour rankings.
    # Flow: cell-line-centred ranking tables -> same-cancer composition table and figure.
    # Reports ranking composition in the configured feature space.
    rule plot_cellline_to_tumour_same_cancer_component_composition:
        """Build same-cancer-type top-50 component composition for cell-line-to-tumour ranks."""
        input:
            script = os.path.join(SCRIPTS_DIR, "plot_cellline_to_tumour_component_composition.R"),
            scores_long = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_tumour_scores_long.tsv.gz"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv")
        output:
            table = PAN_C2T_COMPONENT_SAME_TABLE,
            pdf = PAN_C2T_COMPONENT_SAME_PDF,
            png = PAN_C2T_COMPONENT_SAME_PNG
        params:
            outdir = PAN_RANKING_EVALUATION_DIR,
            top_n = 50
        log: os.path.join(LOGROOT, "plot_cellline_to_tumour_same_cancer_component_composition.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{input.script}" \
              --mode same \
              --scores-long "{input.scores_long}" \
              --components "{input.components}" \
              --top-n {params.top_n} \
              --same-table "{output.table}" \
              --same-pdf "{output.pdf}" \
              --same-png "{output.png}" \
              > "{log}" 2>&1
            test -s "{output.table}" || (echo "ERROR: missing {output.table}" >&2; exit 1)
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.png}" || (echo "ERROR: missing {output.png}" >&2; exit 1)
            '''

    # Rule: plot_cellline_to_tumour_all_top50_component_composition
    # Method role: plotting rule for all top-50 component composition in cell-line-to-tumour rankings.
    # Flow: cell-line-centred ranking tables -> all-cancer-type top-50 composition table and figure.
    # Analysis role: complements same-cancer composition results.
    rule plot_cellline_to_tumour_all_top50_component_composition:
        """Build all-top-50 component composition for cell-line-to-tumour ranks."""
        input:
            script = os.path.join(SCRIPTS_DIR, "plot_cellline_to_tumour_component_composition.R"),
            scores_long = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_tumour_scores_long.tsv.gz"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv"),
            same_table = PAN_C2T_COMPONENT_SAME_TABLE
        output:
            table = PAN_C2T_COMPONENT_ALL_TABLE,
            summary = PAN_C2T_COMPONENT_SUMMARY_TABLE,
            pdf = PAN_C2T_COMPONENT_ALL_PDF,
            png = PAN_C2T_COMPONENT_ALL_PNG
        params:
            outdir = PAN_RANKING_EVALUATION_DIR,
            top_n = 50
        log: os.path.join(LOGROOT, "plot_cellline_to_tumour_all_top50_component_composition.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{input.script}" \
              --mode all \
              --scores-long "{input.scores_long}" \
              --components "{input.components}" \
              --top-n {params.top_n} \
              --same-table "{input.same_table}" \
              --all-table "{output.table}" \
              --summary-table "{output.summary}" \
              --all-pdf "{output.pdf}" \
              --all-png "{output.png}" \
              > "{log}" 2>&1
            test -s "{output.table}" || (echo "ERROR: missing {output.table}" >&2; exit 1)
            test -s "{output.summary}" || (echo "ERROR: missing {output.summary}" >&2; exit 1)
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.png}" || (echo "ERROR: missing {output.png}" >&2; exit 1)
            '''

    # Rule: pan_cancer_ranking_plot_outputs
    # Method role: collector target for pan-cancer ranking evaluation plots and tables.
    # Flow: expected evaluation output paths -> Snakemake dependency aggregation only.
    # Stage output: provides a plotting target without running a new analysis.
    rule pan_cancer_ranking_plot_outputs:
        """Requested retained ranking plot outputs for correlation-based transcriptomic similarity analyses."""
        input:
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.pdf"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.png"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.pdf"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top1_lineage_agreement.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top1_lineage_agreement.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_precision_at_k.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_precision_at_k.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_same_lineage_rank_percentile.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_same_lineage_rank_percentile.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top50_lineage_composition.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top50_lineage_composition.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top50_component_composition.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap", "Fig_cellline_to_tumour_top50_component_composition.png"),
            PAN_T2C_MRR_AT10_PDF,
            PAN_T2C_MRR_AT10_PNG,
            PAN_C2T_COMPONENT_SAME_TABLE,
            PAN_C2T_COMPONENT_ALL_TABLE,
            PAN_C2T_COMPONENT_SUMMARY_TABLE,
            PAN_C2T_COMPONENT_SAME_PDF,
            PAN_C2T_COMPONENT_SAME_PNG,
            PAN_C2T_COMPONENT_ALL_PDF,
            PAN_C2T_COMPONENT_ALL_PNG

    # Rule: build_pan_cancer_bidirectional_ranking_evaluation
    # Method role: ranking evaluation rule checking bidirectional tumour-cell-line ranking tables.
    # Flow: tumour-to-cell-line and cell-line-to-tumour outputs -> crosscheck table and inspection tables.
    # Used by: QA for pan-cancer ranking reports.
    rule build_pan_cancer_bidirectional_ranking_evaluation:
        """Build ECDF and top-1 cancer-type agreement evaluation for bidirectional ranking."""
        input:
            script = os.path.join(SCRIPTS_DIR, "build_bidirectional_ranking_evaluation.R"),
            tumour_scores = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_group_scores.rds"),
            cellline_scores = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_tumour_scores.rds"),
            tumour_rankings = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_to_cellline_group_rankings.tsv"),
            cellline_rankings = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_to_tumour_rankings.tsv"),
            tumour_summary = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_mapping_group_summary.tsv"),
            cellline_summary = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_mapping_summary.tsv"),
            tumour_metrics = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "metrics_summary_group_level.tsv"),
            cellline_metrics = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "metrics_summary.tsv")
        output:
            tumour_ecdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_first_same_lineage_ecdf.tsv"),
            tumour_ecdf_pdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_first_same_lineage_ecdf.pdf"),
            tumour_ecdf_png = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_first_same_lineage_ecdf.png"),
            tumour_first_ranks = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_first_match_ranks.tsv"),
            tumour_confusion_counts = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_top1_confusion_matrix_counts.tsv"),
            tumour_confusion_fraction = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_top1_confusion_matrix_row_fraction.tsv"),
            tumour_confusion_pdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_top1_confusion_matrix.pdf"),
            tumour_confusion_png = os.path.join(PAN_RANKING_EVALUATION_DIR, "tumour_to_cell_line_top1_confusion_matrix.png"),
            cellline_ecdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_first_same_lineage_ecdf.tsv"),
            cellline_ecdf_pdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_first_same_lineage_ecdf.pdf"),
            cellline_ecdf_png = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_first_same_lineage_ecdf.png"),
            cellline_first_ranks = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_first_match_ranks.tsv"),
            cellline_confusion_counts = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_top1_confusion_matrix_counts.tsv"),
            cellline_confusion_fraction = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_top1_confusion_matrix_row_fraction.tsv"),
            cellline_confusion_pdf = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_top1_confusion_matrix.pdf"),
            cellline_confusion_png = os.path.join(PAN_RANKING_EVALUATION_DIR, "cell_line_to_tumour_top1_confusion_matrix.png"),
            crosscheck = os.path.join(PAN_RANKING_EVALUATION_DIR, "ranking_metric_crosscheck.tsv"),
            manifest = os.path.join(PAN_RANKING_EVALUATION_DIR, "ranking_evaluation_manifest.tsv")
        params:
            outdir = PAN_RANKING_EVALUATION_DIR
        log: os.path.join(LOGROOT, "build_pan_cancer_bidirectional_ranking_evaluation.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{input.script}" \
              --tumour-score-rds "{input.tumour_scores}" \
              --cellline-score-rds "{input.cellline_scores}" \
              --tumour-rankings "{input.tumour_rankings}" \
              --cellline-rankings "{input.cellline_rankings}" \
              --tumour-summary "{input.tumour_summary}" \
              --cellline-summary "{input.cellline_summary}" \
              --tumour-metrics "{input.tumour_metrics}" \
              --cellline-metrics "{input.cellline_metrics}" \
              --outdir "{params.outdir}" \
              > "{log}" 2>&1
            test -s "{output.crosscheck}" || (echo "ERROR: missing {output.crosscheck}" >&2; exit 1)
            test -s "{output.manifest}" || (echo "ERROR: missing {output.manifest}" >&2; exit 1)
            '''


# =============================================================================
# PAN-CANCER CELL-LINE-ONLY SIMILARITY NETWORK AND SENSITIVITY ANALYSES
# =============================================================================
# Stage role: analyses cell-line-only similarity and community structure in configured pan-cancer feature spaces.
# Pan-cancer cell line similarity network -- BUILD pipeline
# Enabled via config key: pan_cancer_cell_line_similarity
# ---------------------------------------------------------------------------
_PAN_CELL_LINE_SIM_CFG = config.get("pan_cancer_cell_line_similarity", {})

if _PAN_CELL_LINE_SIM_CFG:
    _CL_SIM_EXPR   = _PAN_CELL_LINE_SIM_CFG.get(
        "expr_rds",
        "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds"
    )
    _CL_SIM_DIR    = _PAN_CELL_LINE_SIM_CFG.get("output_dir", "results/unsupervised/pan_cancer/cell_line_similarity")
    _CL_SIM_COR    = _PAN_CELL_LINE_SIM_CFG.get("correlation", "spearman")
    _CL_SIM_K      = _PAN_CELL_LINE_SIM_CFG.get("k", 20)
    _CL_SIM_ALLOW_EMPTY_GRAPH = bool(_PAN_CELL_LINE_SIM_CFG.get("allow_empty_graph", False))
    _CL_SIM_ALLOW_EMPTY_GRAPH_ARG = "--allow-empty-graph" if _CL_SIM_ALLOW_EMPTY_GRAPH else ""
    _CL_SIM_FEATURES_CLEAN = (
        config.get("defaults", {})
        .get("marker_postprocessing", {})
        .get("pan_cancer", {})
        .get(
            "final_features_clean",
            "results/unsupervised/pan_cancer/feature_space/pan_cancer_features_clean.txt",
        )
    )

    def _infer_pan_cancer_feature_count(path, fallback_path=None):
        candidates = [abspath(path)]
        if fallback_path:
            candidates.append(abspath(fallback_path))
        feature_path = next((candidate for candidate in candidates if os.path.exists(candidate)), None)
        if feature_path is None:
            return 0
        with open(feature_path, "r", encoding="utf-8") as handle:
            return sum(1 for line in handle if line.strip())

    def _resolve_cell_line_expected_genes(raw_value, feature_path, fallback_path=None):
        if raw_value is None:
            return _infer_pan_cancer_feature_count(feature_path, fallback_path)
        if isinstance(raw_value, str) and raw_value.strip().lower() in (
            "", "auto", "infer", "dynamic"
        ):
            return _infer_pan_cancer_feature_count(feature_path, fallback_path)
        return int(raw_value)

    _CL_SIM_EXPECTED_GENES = _resolve_cell_line_expected_genes(
        _PAN_CELL_LINE_SIM_CFG.get("expected_genes", "auto"),
        _CL_SIM_FEATURES_CLEAN,
        _PAN_CELL_LINE_SIM_CFG.get("expected_genes_fallback_file"),
    )
    _CL_SIM_EXPECTED_NODES = _PAN_CELL_LINE_SIM_CFG.get("expected_nodes", 167)
    _CL_SIM_SEED   = PIPELINE_SEED
    _CL_SIM_LEIDEN_RESOLUTION = _PAN_CELL_LINE_SIM_CFG.get("leiden_resolution", 1.0)
    _CL_SIM_LEIDEN_RESOLUTION_SWEEP = _PAN_CELL_LINE_SIM_CFG.get(
        "leiden_resolution_sweep",
        [0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2,
         0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 1.0],
    )
    _CL_SIM_LOUVAIN_RESOLUTION_SWEEP = _PAN_CELL_LINE_SIM_CFG.get(
        "louvain_resolution_sweep",
        [0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2,
         0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 1.0],
    )
    _CL_SIM_FULL_CFG = _PAN_CELL_LINE_SIM_CFG.get("full_expression", {})

    def _format_cell_line_resolution_sweep(raw_value):
        if isinstance(raw_value, str):
            return raw_value
        if isinstance(raw_value, (int, float)):
            return str(raw_value)
        return ",".join(str(v) for v in raw_value)

    _CL_SIM_LEIDEN_RESOLUTION_SWEEP_ARG = _format_cell_line_resolution_sweep(
        _CL_SIM_LEIDEN_RESOLUTION_SWEEP
    )
    _CL_SIM_LOUVAIN_RESOLUTION_SWEEP_ARG = _format_cell_line_resolution_sweep(
        _CL_SIM_LOUVAIN_RESOLUTION_SWEEP
    )
    _CL_SIM_EDGES  = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_graph_edges.tsv")
    _CL_SIM_META   = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_node_metadata.tsv")
    _CL_SIM_COMMS  = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_communities.tsv")
    _CL_SIM_LEIDEN_COMMS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_leiden_communities.tsv")
    _CL_SIM_COMMUNITY_METRICS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_community_metrics.tsv")
    _CL_SIM_CANCER_TYPE_DISCORDANT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_cancer_type_discordant_profiles.tsv")
    _CL_SIM_LEIDEN_SWEEP_ASSIGNMENTS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_leiden_resolution_sweep_assignments.tsv")
    _CL_SIM_LEIDEN_SWEEP_SUMMARY = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_leiden_resolution_sweep_summary.tsv")
    _CL_SIM_LEIDEN_SWEEP_PLOT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_leiden_resolution_sweep.pdf")
    _CL_SIM_LOUVAIN_SWEEP_ASSIGNMENTS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_louvain_resolution_sweep_assignments.tsv")
    _CL_SIM_LOUVAIN_SWEEP_SUMMARY = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_louvain_resolution_sweep_summary.tsv")
    _CL_SIM_LOUVAIN_SWEEP_PLOT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_louvain_resolution_sweep.pdf")
    _CL_SIM_FULL_DIR = abspath(_CL_SIM_FULL_CFG.get(
        "output_dir",
        os.path.join(_CL_SIM_DIR, "dsmz_joint_expression")
    ))
    if not _CL_SIM_FULL_CFG.get("dsmz_vst_rds"):
        raise WorkflowError(
            "pan_cancer_cell_line_similarity.full_expression.dsmz_vst_rds "
            "must explicitly name a transformed DSMZ VST expression object"
        )
    _CL_SIM_FULL_DSMZ_VST = abspath(_CL_SIM_FULL_CFG["dsmz_vst_rds"])
    _CL_SIM_FULL_DSMZ_META = abspath(_CL_SIM_FULL_CFG.get(
        "dsmz_metadata_tsv",
        "results/unsupervised/multicohort_cancer/inputs/joint_metadata.tsv"
    ))
    _CL_SIM_FULL_HEME_VST = abspath(_CL_SIM_FULL_CFG.get(
        "heme_vst_rds",
        "data/heme/heme_joint_expr_vst_batch_corrected.rds"
    ))
    _CL_SIM_FULL_EXPECTED_NODES = _CL_SIM_FULL_CFG.get("expected_nodes", 167)
    _CL_SIM_FULL_EXPR = os.path.join(_CL_SIM_FULL_DIR, "dsmz_joint_expression_cell_line_expr.rds")
    _CL_SIM_FULL_META_OBJECT = os.path.join(_CL_SIM_FULL_DIR, "dsmz_joint_expression_cell_line_metadata.tsv")
    _CL_SIM_FULL_EDGES = os.path.join(_CL_SIM_FULL_DIR, "pan_cancer_cell_line_graph_edges.tsv")
    _CL_SIM_FULL_META = os.path.join(_CL_SIM_FULL_DIR, "pan_cancer_cell_line_node_metadata.tsv")
    _CL_SIM_FULL_LOUVAIN_SWEEP_ASSIGNMENTS = os.path.join(_CL_SIM_FULL_DIR, "dsmz_joint_expression_cell_line_louvain_resolution_sweep_assignments.tsv")
    _CL_SIM_FULL_LOUVAIN_SWEEP_SUMMARY = os.path.join(_CL_SIM_FULL_DIR, "dsmz_joint_expression_cell_line_louvain_resolution_sweep_summary.tsv")
    _CL_SIM_FULL_LOUVAIN_SWEEP_PLOT = os.path.join(_CL_SIM_FULL_DIR, "dsmz_joint_expression_cell_line_louvain_resolution_sweep.pdf")
    _CL_SIM_LAYOUT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_layout.tsv")
    _CL_SIM_VAL_DIR = os.path.join(_CL_SIM_DIR, "community_validation")

    # Rule: build_pan_cancer_cell_line_similarity_graph
    # Method role: graph-construction rule for cell-line-only similarity in the marker-derived feature space.
    # Flow: cell-line-only feature matrix -> cell-line graph edges and metadata.
    # Provides input for cell-line community analysis and resolution sweeps.
    rule build_pan_cancer_cell_line_similarity_graph:
        input:
            expr_rds = _CL_SIM_EXPR
        output:
            edges    = _CL_SIM_EDGES,
            metadata = _CL_SIM_META
        params:
            script  = os.path.join(SCRIPTS_DIR, "build_pan_cancer_cell_line_similarity_graph.R"),
            outdir  = _CL_SIM_DIR,
            k       = _CL_SIM_K,
            correlation = _CL_SIM_COR,
            expected_genes = _CL_SIM_EXPECTED_GENES,
            expected_nodes = _CL_SIM_EXPECTED_NODES,
            allow_empty_graph = _CL_SIM_ALLOW_EMPTY_GRAPH_ARG
        log: os.path.join(LOGROOT, "build_pan_cancer_cell_line_similarity_graph.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --expr-rds "{input.expr_rds}" \
              --output-dir "{params.outdir}" \
              --k {params.k} \
              --correlation "{params.correlation}" \
              --expected-genes {params.expected_genes} \
              --expected-nodes {params.expected_nodes} \
              {params.allow_empty_graph} \
              > "{log}" 2>&1
            test -s "{output.edges}"
            test -s "{output.metadata}"
            '''

    # Rule: compute_pan_cancer_cell_line_communities
    # Method role: graph analysis rule assigning cell-line communities in the marker-derived feature space.
    # Flow: cell-line graph edges and metadata -> community assignments, metrics, and discordance tables.
    # Provides annotations for cell-line-only network validation and plotting.
    rule compute_pan_cancer_cell_line_communities:
        input:
            edges    = _CL_SIM_EDGES,
            metadata = _CL_SIM_META
        output:
            communities = _CL_SIM_COMMS,
            leiden_communities = _CL_SIM_LEIDEN_COMMS,
            community_metrics = _CL_SIM_COMMUNITY_METRICS,
            cancer_type_discordant = _CL_SIM_CANCER_TYPE_DISCORDANT
        params:
            script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_cell_line_communities.R"),
            seed   = _CL_SIM_SEED,
            leiden_resolution = _CL_SIM_LEIDEN_RESOLUTION
        log: os.path.join(LOGROOT, "compute_pan_cancer_cell_line_communities.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript "{params.script}" \
              --edges "{input.edges}" \
              --meta "{input.metadata}" \
              --out "{output.communities}" \
              --leiden-out "{output.leiden_communities}" \
              --community-metrics-out "{output.community_metrics}" \
              --cancer-type-discordant-out "{output.cancer_type_discordant}" \
              --seed {params.seed} \
              --leiden-resolution {params.leiden_resolution} \
              > "{log}" 2>&1
            test -s "{output.communities}"
            test -s "{output.leiden_communities}"
            test -s "{output.community_metrics}"
            test -s "{output.cancer_type_discordant}"
            '''

    # Rule: compute_pan_cancer_cell_line_leiden_resolution_sweep
    # Method role: sensitivity analysis over Leiden resolutions for the cell-line-only graph.
    # Flow: cell-line graph/community outputs -> resolution-sweep assignments, metrics, and sweep plot.
    # Analysis role: assesses community sensitivity without replacing configured community outputs.
    rule compute_pan_cancer_cell_line_leiden_resolution_sweep:
        input:
            edges    = _CL_SIM_EDGES,
            metadata = _CL_SIM_META,
            communities = _CL_SIM_COMMS,
            leiden_communities = _CL_SIM_LEIDEN_COMMS,
            community_metrics = _CL_SIM_COMMUNITY_METRICS,
            cancer_type_discordant = _CL_SIM_CANCER_TYPE_DISCORDANT
        output:
            leiden_sweep_assignments = _CL_SIM_LEIDEN_SWEEP_ASSIGNMENTS,
            leiden_sweep_summary = _CL_SIM_LEIDEN_SWEEP_SUMMARY,
            leiden_sweep_plot = _CL_SIM_LEIDEN_SWEEP_PLOT
        params:
            script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_cell_line_communities.R"),
            seed   = _CL_SIM_SEED,
            leiden_resolution = _CL_SIM_LEIDEN_RESOLUTION,
            leiden_resolution_sweep = _CL_SIM_LEIDEN_RESOLUTION_SWEEP_ARG
        log: os.path.join(LOGROOT, "compute_pan_cancer_cell_line_leiden_resolution_sweep.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript "{params.script}" \
              --edges "{input.edges}" \
              --meta "{input.metadata}" \
              --out "{input.communities}" \
              --leiden-out "{input.leiden_communities}" \
              --community-metrics-out "{input.community_metrics}" \
              --cancer-type-discordant-out "{input.cancer_type_discordant}" \
              --seed {params.seed} \
              --leiden-resolution {params.leiden_resolution} \
              --leiden-resolution-sweep "{params.leiden_resolution_sweep}" \
              --leiden-sweep-assignments-out "{output.leiden_sweep_assignments}" \
              --leiden-sweep-summary-out "{output.leiden_sweep_summary}" \
              --leiden-sweep-plot-out "{output.leiden_sweep_plot}" \
              > "{log}" 2>&1
            test -s "{output.leiden_sweep_assignments}"
            test -s "{output.leiden_sweep_summary}"
            test -s "{output.leiden_sweep_plot}"
            '''

    # Rule: compute_pan_cancer_cell_line_louvain_resolution_sweep
    # Method role: sensitivity analysis over Louvain resolutions for the cell-line-only graph.
    # Flow: cell-line graph edges and metadata -> resolution-sweep assignments, metrics, and sweep plot.
    # Analysis role: compares community-resolution behaviour across algorithms.
    rule compute_pan_cancer_cell_line_louvain_resolution_sweep:
        input:
            edges    = _CL_SIM_EDGES,
            metadata = _CL_SIM_META
        output:
            louvain_sweep_assignments = _CL_SIM_LOUVAIN_SWEEP_ASSIGNMENTS,
            louvain_sweep_summary = _CL_SIM_LOUVAIN_SWEEP_SUMMARY,
            louvain_sweep_plot = _CL_SIM_LOUVAIN_SWEEP_PLOT
        params:
            script = os.path.join(SCRIPTS_DIR, "compute_cell_line_louvain_resolution_sweep.R"),
            seed   = _CL_SIM_SEED,
            louvain_resolution_sweep = _CL_SIM_LOUVAIN_RESOLUTION_SWEEP_ARG,
            feature_source = "marker_derived_pan_cancer_feature_space"
        log: os.path.join(LOGROOT, "compute_pan_cancer_cell_line_louvain_resolution_sweep.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript "{params.script}" \
              --edges "{input.edges}" \
              --meta "{input.metadata}" \
              --feature-source "{params.feature_source}" \
              --resolution-sweep "{params.louvain_resolution_sweep}" \
              --seed {params.seed} \
              --assignments-out "{output.louvain_sweep_assignments}" \
              --summary-out "{output.louvain_sweep_summary}" \
              --plot-out "{output.louvain_sweep_plot}" \
              > "{log}" 2>&1
            test -s "{output.louvain_sweep_assignments}"
            test -s "{output.louvain_sweep_summary}"
            test -s "{output.louvain_sweep_plot}"
            '''

    # Rule: build_dsmz_joint_expression_cell_line_object
    # Method role: data-preparation rule building a full DSMZ cell-line expression object for sensitivity analysis.
    # Flow: DSMZ and HEME VST inputs plus metadata -> joint cell-line expression object and metadata table.
    # Provides full-expression cell-line similarity sensitivity analysis.
    rule build_dsmz_joint_expression_cell_line_object:
        input:
            dsmz_vst = _CL_SIM_FULL_DSMZ_VST,
            dsmz_metadata = _CL_SIM_FULL_DSMZ_META,
            heme_vst = _CL_SIM_FULL_HEME_VST
        output:
            expr_rds = _CL_SIM_FULL_EXPR,
            metadata = _CL_SIM_FULL_META_OBJECT
        params:
            script = os.path.join(SCRIPTS_DIR, "build_dsmz_joint_expression_cell_line_object.R"),
            outdir = _CL_SIM_FULL_DIR
        log: os.path.join(LOGROOT, "build_dsmz_joint_expression_cell_line_object.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --dsmz-vst "{input.dsmz_vst}" \
              --dsmz-meta "{input.dsmz_metadata}" \
              --heme-vst "{input.heme_vst}" \
              --output "{output.expr_rds}" \
              --output-metadata "{output.metadata}" \
              > "{log}" 2>&1
            test -s "{output.expr_rds}"
            test -s "{output.metadata}"
            '''

    # Rule: build_dsmz_joint_expression_cell_line_similarity_graph
    # Method role: graph-construction rule for full-expression DSMZ cell-line similarity sensitivity analysis.
    # Flow: full DSMZ joint expression object -> similarity graph edges and metadata.
    # Provides full-expression Louvain resolution sweep.
    rule build_dsmz_joint_expression_cell_line_similarity_graph:
        input:
            expr_rds = _CL_SIM_FULL_EXPR
        output:
            edges    = _CL_SIM_FULL_EDGES,
            metadata = _CL_SIM_FULL_META
        params:
            script  = os.path.join(SCRIPTS_DIR, "build_pan_cancer_cell_line_similarity_graph.R"),
            outdir  = _CL_SIM_FULL_DIR,
            k       = _CL_SIM_K,
            correlation = _CL_SIM_COR,
            expected_nodes = _CL_SIM_FULL_EXPECTED_NODES,
            allow_empty_graph = _CL_SIM_ALLOW_EMPTY_GRAPH_ARG
        log: os.path.join(LOGROOT, "build_dsmz_joint_expression_cell_line_similarity_graph.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --expr-rds "{input.expr_rds}" \
              --output-dir "{params.outdir}" \
              --k {params.k} \
              --correlation "{params.correlation}" \
              --expected-genes 0 \
              --expected-nodes {params.expected_nodes} \
              {params.allow_empty_graph} \
              > "{log}" 2>&1
            test -s "{output.edges}"
            test -s "{output.metadata}"
            '''

    # Rule: compute_dsmz_joint_expression_cell_line_louvain_resolution_sweep
    # Method role: sensitivity analysis over Louvain resolutions for the full-expression DSMZ graph.
    # Flow: full-expression graph sidecars -> resolution-sweep assignments, metrics, and sweep plot.
    # Analysis role: compares feature-space and full-expression community behaviour.
    rule compute_dsmz_joint_expression_cell_line_louvain_resolution_sweep:
        input:
            edges    = _CL_SIM_FULL_EDGES,
            metadata = _CL_SIM_FULL_META
        output:
            louvain_sweep_assignments = _CL_SIM_FULL_LOUVAIN_SWEEP_ASSIGNMENTS,
            louvain_sweep_summary = _CL_SIM_FULL_LOUVAIN_SWEEP_SUMMARY,
            louvain_sweep_plot = _CL_SIM_FULL_LOUVAIN_SWEEP_PLOT
        params:
            script = os.path.join(SCRIPTS_DIR, "compute_cell_line_louvain_resolution_sweep.R"),
            seed   = _CL_SIM_SEED,
            louvain_resolution_sweep = _CL_SIM_LOUVAIN_RESOLUTION_SWEEP_ARG,
            feature_source = "full_dsmz_joint_expression"
        log: os.path.join(LOGROOT, "compute_dsmz_joint_expression_cell_line_louvain_resolution_sweep.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript "{params.script}" \
              --edges "{input.edges}" \
              --meta "{input.metadata}" \
              --feature-source "{params.feature_source}" \
              --resolution-sweep "{params.louvain_resolution_sweep}" \
              --seed {params.seed} \
              --assignments-out "{output.louvain_sweep_assignments}" \
              --summary-out "{output.louvain_sweep_summary}" \
              --plot-out "{output.louvain_sweep_plot}" \
              > "{log}" 2>&1
            test -s "{output.louvain_sweep_assignments}"
            test -s "{output.louvain_sweep_summary}"
            test -s "{output.louvain_sweep_plot}"
            '''

    # Rule: compute_pan_cancer_cell_line_validation
    # Method role: validation rule computing modularity and assortativity for cell-line communities.
    # Flow: cell-line graph, communities, and metadata -> validation tables.
    # Reports graph/community metrics without changing community assignments.
    rule compute_pan_cancer_cell_line_validation:
        input:
            edges       = _CL_SIM_EDGES,
            communities = _CL_SIM_COMMS,
            leiden_communities = _CL_SIM_LEIDEN_COMMS,
            metadata    = _CL_SIM_META
        output:
            modularity    = os.path.join(_CL_SIM_VAL_DIR, "validation_modularity.tsv"),
            assortativity = os.path.join(_CL_SIM_VAL_DIR, "validation_assortativity.tsv")
        params:
            script  = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_cell_line_validation.R"),
            out_dir = _CL_SIM_VAL_DIR
        log: os.path.join(LOGROOT, "compute_pan_cancer_cell_line_validation.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.out_dir}"
            Rscript "{params.script}" \
              --edges "{input.edges}" \
              --communities "{input.communities}" \
              --leiden-communities "{input.leiden_communities}" \
              --metadata "{input.metadata}" \
              --out-dir "{params.out_dir}" \
              > "{log}" 2>&1
            test -s "{output.modularity}"
            test -s "{output.assortativity}"
            '''

    # Rule: compute_pan_cancer_cell_line_layout
    # Method role: layout rule computing fixed coordinates for the cell-line-only graph.
    # Flow: cell-line graph edges -> layout table.
    # Provides input for manual two-panel cell-line network plotting.
    rule compute_pan_cancer_cell_line_layout:
        input:
            edges = _CL_SIM_EDGES
        output:
            layout = _CL_SIM_LAYOUT
        params:
            script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_cell_line_layout.R"),
            seed   = _CL_SIM_SEED
        log: os.path.join(LOGROOT, "compute_pan_cancer_cell_line_layout.log")
        conda: CONDA_ENV_R
        shell: "Rscript {params.script} --edges {input.edges} --out {output.layout} --seed {params.seed} > {log} 2>&1 && test -s {output.layout}"


# =============================================================================
# OPTIONAL / MANUAL PLOTTING RULES
# =============================================================================
# Stage role: contains configured reporting plots that can be added to the automated target list.
# Pan-cancer cell line similarity network -- two-panel figure
# Cell-line-to-cell-line transcriptomic similarity (not cell-line-to-tumour).
# Enabled and configured via config key: pan_cancer_cell_line_plot
#
# Required config keys:
#   pan_cancer_cell_line_plot:
#     edges_tsv:       path to edge list TSV  (columns: from, to, weight)
#     communities_tsv: path to community TSV  (columns: sample, component, cancer_type)
#     layout_tsv:      path to layout TSV     (columns: sample, x, y)
#     output_dir:      directory for output PDF and PNG
#
# Manual target command: snakemake plot_pan_cancer_cell_line_two_panel --cores 1
# Included in pan-cancer marker-postprocessing targets when configured.
# ---------------------------------------------------------------------------
_PAN_CELL_LINE_PLOT_CFG = config.get("pan_cancer_cell_line_plot", {})

if _PAN_CELL_LINE_PLOT_CFG:
    _PLOT_EDGES  = abspath(_PAN_CELL_LINE_PLOT_CFG["edges_tsv"])
    _PLOT_COMMS  = abspath(_PAN_CELL_LINE_PLOT_CFG["communities_tsv"])
    _PLOT_LAYOUT = abspath(_PAN_CELL_LINE_PLOT_CFG["layout_tsv"])
    _PLOT_OUTDIR = abspath(_PAN_CELL_LINE_PLOT_CFG["output_dir"])
    _PLOT_PDF = os.path.join(
        _PLOT_OUTDIR,
        "Fig_pan_cancer_cell_line_similarity_network_cancer_type_community.pdf",
    )
    _PLOT_PNG = os.path.join(
        _PLOT_OUTDIR,
        "Fig_pan_cancer_cell_line_similarity_network_cancer_type_community.png",
    )
    _PLOT_LOUVAIN_ALIAS_PDF = os.path.join(
        _PLOT_OUTDIR,
        _PAN_CELL_LINE_PLOT_CFG.get(
            "louvain_alias_pdf",
            "Fig_pan_cancer_cell_line_louvain.pdf",
        ),
    )
    _PLOT_FIGURE_WIDTH = float(_PAN_CELL_LINE_PLOT_CFG.get("figure_width", 11.2))
    _PLOT_FIGURE_HEIGHT = float(_PAN_CELL_LINE_PLOT_CFG.get("figure_height", 5.0))
    _PLOT_PANEL_WIDTH = float(_PAN_CELL_LINE_PLOT_CFG.get("panel_width", 4.75))
    _PLOT_LEGEND_WIDTH = float(_PAN_CELL_LINE_PLOT_CFG.get("legend_width", 1.70))
    _PLOT_LEGEND_MODE = _PAN_CELL_LINE_PLOT_CFG.get("legend_mode", "separate")

    # Rule: plot_pan_cancer_cell_line_two_panel
    # Method role: plotting rule for a two-panel cell-line-only similarity network figure.
    # Flow: edge list, communities, layout, and config-driven dimensions -> cancer-type/community PDF and PNG.
    # Purpose: reporting figure and stable Louvain/community PDF alias.
    rule plot_pan_cancer_cell_line_two_panel:
        """Two-panel figure: pan-cancer cell line transcriptomic similarity
        network (cell-line-to-cell-line, not cell-line-to-tumour).
        Panel A: cancer type coloured by Okabe-Ito palette.
        Panel B: Louvain communities (Dark2) with convex hulls and centroid labels.
        Accepts any edge list, community table, and layout produced by this
        or any compatible pipeline run."""
        input:
            edges  = _PLOT_EDGES,
            comms  = _PLOT_COMMS,
            layout = _PLOT_LAYOUT
        output:
            pdf = _PLOT_PDF,
            png = _PLOT_PNG,
            louvain_pdf = _PLOT_LOUVAIN_ALIAS_PDF
        params:
            script = os.path.join(SCRIPTS_DIR, "plot_pan_cancer_two_panel.R"),
            outdir = _PLOT_OUTDIR,
            figure_width = _PLOT_FIGURE_WIDTH,
            figure_height = _PLOT_FIGURE_HEIGHT,
            panel_width = _PLOT_PANEL_WIDTH,
            legend_width = _PLOT_LEGEND_WIDTH,
            legend_mode = _PLOT_LEGEND_MODE
        log: os.path.join(LOGROOT, "plot_pan_cancer_cell_line_two_panel.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p "{params.outdir}"
            Rscript "{params.script}" \
              --edges       "{input.edges}" \
              --communities "{input.comms}" \
              --layout      "{input.layout}" \
              --out-dir     "{params.outdir}" \
              --width       {params.figure_width} \
              --height      {params.figure_height} \
              --panel-width {params.panel_width} \
              --legend-width {params.legend_width} \
              --legend-mode "{params.legend_mode}" \
              > "{log}" 2>&1
            [ "{output.pdf}" = "{output.louvain_pdf}" ] || cp "{output.pdf}" "{output.louvain_pdf}"
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.png}" || (echo "ERROR: missing {output.png}" >&2; exit 1)
            test -s "{output.louvain_pdf}" || (echo "ERROR: missing {output.louvain_pdf}" >&2; exit 1)
            '''

# =============================================================================
# FINAL PIPELINE_TARGET ASSEMBLY
# =============================================================================
# Stage role: assembles the default target list after all optional section variables are defined.

# Method helper: assembles the default target list after optional analysis sections have declared their variables.
def build_pipeline_targets():
    targets = []

    if not IS_PAN_CANCER_PROFILE:
        targets.extend([
            os.path.join(UNSUP_REL, "tumour_neighbourhoods", "final_consensus_all", "resolved_dsmz_neighbours.tsv"),
            os.path.join(UNSUP_REL, "tumour_neighbourhoods", "final_consensus_all", "resolved_graph_input_validation.tsv"),
            P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".pdf",
            P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".png",
            P_CONS_RESOLVED_NEIGHBOURHOOD_GRAPH_PREFIX + ".svg",
            P_CONS_RESOLVED_EDGES_TSV,
            P_CONS_RESOLVED_NODE_STATS_TSV,
            P_CONS_ANCHOR_AUDIT_TSV,
            # P_CONS_RESOLVED_NODE_LABELS_TSV is produced by the
            # multicohort-only inspection rule and is already included for
            # that profile via PAN_CANCER_GRAPH_INSPECTION_TARGETS below;
            # listing it unconditionally breaks DAG construction for cohort
            # profiles, which have no producing rule.
            P_CONS_PER_CELLLINE_FEATURE_DISTANCE_PREFIX + ".pdf",
            P_CONS_PER_CELLLINE_FEATURE_DISTANCE_PREFIX + ".png",
            P_CONS_MULTI_REP_MAJORITY_PREFIX + ".pdf",
            P_CONS_MULTI_REP_MAJORITY_PREFIX + ".png",
            P_CONS_MULTI_REP_MAJORITY_PREFIX + ".svg",
            P_CONS_MULTI_REP_MAJORITY_EDGES_TSV,
            P_CONS_MULTI_REP_MAJORITY_NODE_STATS_TSV,
            P_CONS_MULTI_REP_MAJORITY_COMPONENTS_TSV,
            P_CONS_MULTI_REP_EDGE_SUPPORT_TSV,
            POST_RESOLUTION_COMBINED_STRATIFICATION_TSV,
            P_CONS_MULTI_REP_UNION_PREFIX + ".pdf",
            P_CONS_MULTI_REP_UNION_PREFIX + ".png",
            P_CONS_MULTI_REP_UNION_PREFIX + ".svg",
            P_CONS_MULTI_REP_UNION_EDGES_TSV,
            P_CONS_MULTI_REP_UNION_NODE_STATS_TSV,
            P_CONS_MULTI_REP_UNION_COMPONENTS_TSV,
            P_CONS_MULTI_REP_UNION_FULL_NODE_LABELS_TSV,
            *PAN_CANCER_GRAPH_INSPECTION_TARGETS,
            *COMMUNITY_STABILITY_TARGETS,
            os.path.join(STUDY_OUTPUT_DIR_REL, "study_question.txt"),
            os.path.join(STUDY_OUTPUT_DIR_REL, "cohort_manifest.tsv"),
            os.path.join(STUDY_OUTPUT_DIR_REL, "cohort_labels.tsv"),
            os.path.join(STUDY_OUTPUT_DIR_REL, "candidate_inference.tsv"),
            os.path.join(STUDY_OUTPUT_DIR_REL, "endpoint_manifest.tsv"),
            os.path.join(VALIDATION_OUTPUT_DIR_REL, "model_selection_summary.tsv"),
            os.path.join(VALIDATION_OUTPUT_DIR_REL, "neighbourhood_permutation_summary.tsv"),
            os.path.join(VALIDATION_OUTPUT_DIR_REL, "random_baseline_summary.tsv"),
            os.path.join(VALIDATION_OUTPUT_DIR_REL, "silhouette_report.tsv"),
        ])

        # Configured similarity metrics are part of default graph behaviour.
        # Every metric declared in patient_referenced_graph.similarity_metrics
        # produces its own per-representation graphs, resolved-neighbour table,
        # input-validation report and support networks, so the default DAG can
        # no longer ignore the configured metric set. The primary branch above
        # uses patient_referenced_graph.primary_similarity_metric, which is
        # validated to be a member of the same configured list.
        for _metric in PATIENT_REFERENCED_SIMILARITY_METRICS:
            _metric_dir = os.path.join(PATIENT_REFERENCED_SIMILARITY_ALL_DIR, _metric)
            targets.extend([
                os.path.join(_metric_dir, "resolved_dsmz_neighbours.tsv"),
                os.path.join(_metric_dir, "resolved_graph_input_validation.tsv"),
                os.path.join(_metric_dir, "multi_representation_majority_threshold_edges.tsv"),
                os.path.join(_metric_dir, "multi_representation_majority_threshold_edge_support.tsv"),
                os.path.join(_metric_dir, "multi_representation_union_edges.tsv"),
            ])

        if PATIENT_REFERENCED_METRIC_COMPARISON_ENABLED:
            targets.extend([
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_representation_graphs.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_pairwise_similarity.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_edge_agreement.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_graph_comparison.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "pearson_vs_jaccard_resolved_neighbours.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "similarity_graph_provenance.tsv"),
                os.path.join(PATIENT_REFERENCED_SIMILARITY_COMPARE_DIR, "clustering_method_consensus_resolution.tsv"),
            ])

        if DESEQ2_ENABLED:
            targets.append(os.path.join(DESEQ2_COMP_DIR, ".done"))

    if IS_PAN_CANCER_PROFILE and MARKER_POST_ENABLED:
        targets.append(PAN_EXPR_RDS)
        targets.append(mapping_metrics_summary())
        if FUNCTIONAL_ENRICHMENT_ENABLED:
            targets.extend([
                FUNCTIONAL_ENRICHMENT_QUERY_MANIFEST,
                FUNCTIONAL_ENRICHMENT_QUERY_INTEGRITY,
                FUNCTIONAL_ENRICHMENT_AGGREGATE_SUCCESS,
                FUNCTIONAL_ENRICHMENT_PRIMARY_HEATMAP_PDF,
                FUNCTIONAL_ENRICHMENT_CONTRAST_SUPPORT_HEATMAP_PDF,
            ])
        targets.extend(
            PAN_ALIGNMENT_UMAP_PDFS + PAN_ALIGNMENT_UMAP_SVGS + PAN_ALIGNMENT_UMAP_PNGS +
            PAN_ALIGNMENT_UMAP_SOURCE_PDFS + PAN_ALIGNMENT_UMAP_SOURCE_SVGS + PAN_ALIGNMENT_UMAP_SOURCE_PNGS +
            PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PDFS + PAN_ALIGNMENT_UMAP_SOURCE_CANCER_SVGS + PAN_ALIGNMENT_UMAP_SOURCE_CANCER_PNGS +
            PAN_ALIGNMENT_UMAP_COORDS + PAN_ALIGNMENT_UMAP_SOURCE_COMPOSITION + [PAN_ALIGNMENT_UMAP_SUMMARY] +
            VST_ALLGENE_UMAP_PDFS + VST_ALLGENE_UMAP_SVGS + VST_ALLGENE_UMAP_PNGS +
            VST_ALLGENE_UMAP_SOURCE_PDFS + VST_ALLGENE_UMAP_SOURCE_SVGS + VST_ALLGENE_UMAP_SOURCE_PNGS +
            VST_ALLGENE_UMAP_SOURCE_CANCER_PDFS + VST_ALLGENE_UMAP_SOURCE_CANCER_SVGS + VST_ALLGENE_UMAP_SOURCE_CANCER_PNGS +
            VST_ALLGENE_UMAP_COORDS + VST_ALLGENE_UMAP_SOURCE_COMPOSITION + [VST_ALLGENE_UMAP_SUMMARY] +
            RAW_ALLGENE_UMAP_TARGETS
        )
        if _PAN_CELL_LINE_PLOT_CFG:
            targets.extend([_PLOT_PDF, _PLOT_PNG, _PLOT_LOUVAIN_ALIAS_PDF])

        targets.extend([
            os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_graph.pdf"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.pdf"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_rank_ecdf.png"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.pdf"),
            os.path.join(PAN_FIG_DIR, "ecdf_plots", "Fig_tumour_to_cellline_top10_fraction.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top1_lineage_agreement.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top1_lineage_agreement.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_precision_at_k.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_precision_at_k.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_same_lineage_rank_percentile.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_same_lineage_rank_percentile.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top50_lineage_composition.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top50_lineage_composition.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top50_component_composition.pdf"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "Fig_cellline_to_tumour_top50_component_composition.png"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "cellline_topk_metrics.tsv"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "cellline_centred_rank_summary.tsv"),
            os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                         "cellline_similarity_precision_provenance.tsv"),
            os.path.join(PAN_RANKING_EVALUATION_DIR,
                         "ranking_metric_crosscheck.tsv"),
            PAN_T2C_MRR_AT10_BY_TUMOUR,
            PAN_T2C_MRR_AT10_BY_COHORT,
            PAN_T2C_MRR_AT10_PDF,
            PAN_T2C_MRR_AT10_PNG,
            PAN_C2T_COMPONENT_SAME_TABLE,
            PAN_C2T_COMPONENT_ALL_TABLE,
            PAN_C2T_COMPONENT_SUMMARY_TABLE,
            PAN_C2T_COMPONENT_SAME_PDF,
            PAN_C2T_COMPONENT_SAME_PNG,
            PAN_C2T_COMPONENT_ALL_PDF,
            PAN_C2T_COMPONENT_ALL_PNG,
        ])

    if IS_MULTICOHORT_PROFILE:
        targets.extend([
            P_CONS_LEIDEN_COMMUNITIES_TSV,
            P_CONS_LEIDEN_SUMMARY_TSV,
            P_CONS_LEIDEN_MODULARITY_TSV,
            P_CONS_LEIDEN_LAYOUT_TSV,
            P_CONS_LEIDEN_FIG_PREFIX + ".pdf",
            P_CONS_LEIDEN_FIG_PREFIX + ".png",
        ])
    return targets


PIPELINE_TARGET = build_pipeline_targets()
if not PIPELINE_TARGET:
    raise WorkflowError(
        f"No pipeline targets resolved for profile '{profile_name}'. "
        "Check marker_postprocessing.enabled and functional_enrichment.enabled."
    )

# =============================================================================
# INCLUDE FILES
# =============================================================================
# Stage role: loads external rule files after their required variables and target paths are defined.
# Canonical graph-derived functional-enrichment workflow.
include: "rules/graph_derived_functional_enrichment.smk"

# The synthetic test scaffold is opt-in and cannot collide with production rules by default.
if bool(config.get("minimal_example_test", {}).get("enabled", False)):
    include: "rules/minimal_example_test.smk"
