#!/usr/bin/env snakemake

# Unsupervised Clustering and Tumour Neighbourhood Analysis Pipeline

import json
import os
import re
import shlex
import sys
from pathlib import Path
from snakemake.shell import shell
from snakemake.io import expand, directory

# -----------------------------------------------------------------------------
# PATHS / PROFILE SETUP (MUST BE DEFINED BEFORE RULES)
# -----------------------------------------------------------------------------
# Anchors all path resolution to the main Snakefile location, preventing
# inadvertent resolution relative to Snakemake's site-packages directory when
# --directory or the working directory differs from the pipeline root.
_snakefile_path = getattr(workflow, "main_snakefile", None) or __file__
BASE = os.path.dirname(os.path.abspath(os.path.realpath(_snakefile_path)))
SCRIPTS_DIR = os.path.join(BASE, "scripts")
PIPE_ROOT = BASE

# Centralises conda environment paths to avoid hardcoded absolute strings in rules.
CONDA_ENV_R = os.path.join(PIPE_ROOT, "envs", "tcga-r-env.yaml")
CONDA_ENV_PY = os.path.join(PIPE_ROOT, "envs", "python-graph-env.yaml")
CONDA_ENV_QC = os.path.join(PIPE_ROOT, "envs", "tumour_nh_qc.yaml")

# Resolves the config file path once at pipeline load time.
CFGFILE_ABS = os.path.join(PIPE_ROOT, "config", "config.yaml")

# Conda environments — single source of truth; avoids hardcoded absolute paths in every rule.
CONDA_ENV_R  = os.path.join(PIPE_ROOT, "envs", "tcga-r-env.yaml")
CONDA_ENV_PY = os.path.join(PIPE_ROOT, "envs", "python-graph-env.yaml")
CONDA_ENV_QC = os.path.join(PIPE_ROOT, "envs", "tumour_nh_qc.yaml")

# Instructs Snakemake to load configuration from the resolved absolute path.
configfile: CFGFILE_ABS

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
IS_DRYRUN = any(arg in ("-n", "--dry-run", "--dryrun") for arg in sys.argv)

# Defines a safe default for PIPELINE_TARGET so DAG parsing never fails before
# profile-specific targets are resolved later in the Snakefile.
PIPELINE_TARGET = []

# -----------------------------------------------------------------------------
# DEFAULT TARGET (FIRST RULE)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# RULE: all
#
# This is the default entry point of the whole workflow.
# When you run Snakemake without naming a specific rule, Snakemake tries to build
# whatever appears here. In this pipeline, the final target is chosen later based
# on the active profile, so this rule acts like a profile-aware 'run everything'.
# -----------------------------------------------------------------------------
rule all:
    input:
        # Lazily evaluated so the target list is not resolved until after profile
        # setup and PIPELINE_TARGET is overwritten with profile-specific outputs.
        lambda wc: PIPELINE_TARGET,
        # Pan-cancer UMAP (optional, can be run independently)
        # Uncomment to include in default target:
        # os.path.join(BASE, "results", "unsupervised", "pan_cancer", "umap", "pan_cancer_mx_umap_coords.tsv")


# Utility function: converts a config-relative path to an absolute path.
def abspath(p):
    return p if os.path.isabs(p) else os.path.abspath(os.path.join(BASE, p))


STUDY_DESIGN_REL = config.get("study_design", {}).get("file", "config/study_design.yaml")
STUDY_DESIGN_FILE = abspath(STUDY_DESIGN_REL)

def deep_merge(base, override):
    """Recursively merges two configuration dictionaries without mutating either input."""
    merged = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(merged.get(k), dict):
            merged[k] = deep_merge(merged[k], v)
        else:
            merged[k] = v
    return merged


def get_profile_cfg(profile):
    """Returns the merged configuration for a named profile, used for cross-profile pan-cancer access."""
    default_cfg = config.get("defaults", {})
    profiles = config.get("profiles", {})
    if profile not in profiles:
        raise ValueError(f"pan_cancer profile '{profile}' not found in config.profiles")
    return deep_merge(default_cfg, profiles[profile])


def profile_vst_joint_abs(profile):
    """Resolves the absolute path to vst_joint_rds for a given profile."""
    cfgp = get_profile_cfg(profile)
    p = cfgp.get("paths", {}).get("vst_joint_rds")
    if not p:
        raise KeyError(f"paths.vst_joint_rds missing for profile '{profile}'")
    return abspath(p)


if "profiles" in config:
    profiles = config.get("profiles", {})
    # Reads the mandatory pipeline_profile flag; raises an informative error if absent.
    profile_name = config.get("pipeline_profile")

    # Enforces strict profile selection: pipeline_profile must be supplied via --config.
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


# Beginner note:
# From here down, the Snakefile converts config values into Python variables
# that make the rest of the workflow easier to read and reuse.

analysis_cfg = cfg.get("analysis", {})
# Controls whether PAM50 gene-set directions are included in the analysis.
USE_PAM50 = analysis_cfg.get("use_pam50", False)

# Multicohort benchmark flag that disables PCA pre-processing steps across all
# profiles to ensure comparability when merging cross-cancer datasets.
MULTICOHORT_CFG = config.get("multicohort_cancer", {})
DISABLE_PCA_EVERYWHERE = bool(MULTICOHORT_CFG.get("disable_pca_everywhere", False))
vprint(f"[Snakefile] DISABLE_PCA_EVERYWHERE = {DISABLE_PCA_EVERYWHERE}")

# Resolves the list of feature selection methods and distance metrics from config,
# with sensible defaults if neither feature_sets nor features keys are present.
FEATURE_METHODS = cfg.get("feature_sets", {}).get(
    "methods",
    cfg.get("features", {}).get(
        "methods", ["Variance", "MAD", "MeanAbsDev", "Entropy", "PCA", "Spearman", "MX", "kTotal", "HVG"]
    ),
)
DISTANCES = cfg.get("feature_sets", {}).get("distances", ["euc", "corr"])

# Per-method top-N gene counts.
# Defaults: Variance/MAD/MeanAbsDev/Entropy/PCA → 3000; Spearman/MX/kTotal → 500.
# Override per-method in config under feature_selection.method_topn.<method>.
_TOPN_DEFAULTS = {
    "Variance":   3000,
    "MAD":        3000,
    "MeanAbsDev": 3000,
    "Entropy":    3000,
    "PCA":        3000,
    "Spearman":   500,
    "MX":         500,
    "kTotal":     500,
    "HVG":        3000,
}
_topn_cfg = cfg.get("feature_selection", {}).get("method_topn", {})
METHOD_TOPN = {m: int(_topn_cfg.get(m, _TOPN_DEFAULTS.get(m, 500))) for m in FEATURE_METHODS}

def method_topn(m):
    """Return the top-N gene count for feature method m."""
    return METHOD_TOPN.get(m, 500)

TOPN = method_topn("MX")


def build_directions():
    """
    Constructs the full list of analysis directions as method × distance combinations.
    Config-explicit direction lists are authoritative for the active profile.
    If no explicit list is present, fall back to method × distance combinations
    and append PAM50 only when the legacy use_pam50 flag requests it.
    """
    feature_dirs = [f"{m}_{d}" for m in FEATURE_METHODS for d in DISTANCES]
    pam50_dirs = ["pam50_euc", "pam50_corr"] if USE_PAM50 else []
    explicit = cfg.get("tumour_neighbourhoods", {}).get("directions")
    if explicit is not None:
        return explicit
    return feature_dirs + pam50_dirs


AGN_DIRECTIONS = build_directions()
HAS_PAM50_DIRECTIONS = any(d.startswith("pam50") for d in AGN_DIRECTIONS)

if HAS_PAM50_DIRECTIONS and not cfg.get("pam50", {}).get("ensembl_gene_list") and not cfg.get("features", {}).get("pam50_ensembl_gene_list"):
    raise ValueError("PAM50 directions are configured, but pam50.ensembl_gene_list is not set in config.")


def cfgrel(*keys):
    """
    Traverses the merged config by key path and returns the value as-is (relative).
    Used for rule output declarations to keep paths portable across environments.
    """
    d = cfg
    for k in keys:
        d = d[k]
    return d


def cfgabs(*keys):
    """
    Traverses the merged config by key path and returns the value as an absolute path.
    Used for rule inputs and params where scripts require fully resolved paths.
    """
    return abspath(cfgrel(*keys))


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


# Sets strict shell execution flags and pins thread counts for linear algebra libraries
# to prevent unintended CPU over-subscription on shared HPC nodes.
shell.prefix("set -euo pipefail; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ")

# Resolves the unsupervised output root: relative for rule output declarations
# (Snakemake DAG portability), absolute for params/shell arguments (script compatibility).
UNSUP_BASE = os.path.normpath(cfgget_path_rel("results/unsupervised", "paths", "unsup_root"))
vprint(f"[Snakefile] UNSUP_BASE (from config) = {UNSUP_BASE}")

# Converts absolute config-sourced paths back to PIPE_ROOT-relative for output declarations.
if os.path.isabs(UNSUP_BASE):
    UNSUP_BASE = os.path.relpath(UNSUP_BASE, PIPE_ROOT)

# Ensures the profile name is embedded as the final directory component of the output root,
# so results from different profiles never collide under the same parent directory.
# Beginner note:
# The variables below define where outputs live and how directions are named.
# Think of them as the bookkeeping layer that keeps many clustering runs organised.

# Append the profile only when it is not already present in the configured
# path. This avoids doubled roots such as
# results/unsupervised/multicohort_cancer/unsupervised/multicohort_cancer.
if profile_name not in os.path.normpath(UNSUP_BASE).split(os.sep):
    UNSUP_REL = os.path.join(UNSUP_BASE, profile_name)
else:
    UNSUP_REL = UNSUP_BASE

# Absolute counterpart used only in params, inputs, and script arguments.
UNSUP = os.path.join(PIPE_ROOT, UNSUP_REL)
vprint(f"[Snakefile] UNSUP_REL (final, relative) = {UNSUP_REL}")
vprint(f"[Snakefile] UNSUP (final, absolute) = {UNSUP}")

# Resolves the joint VST-normalised expression matrix (absolute; external input).
VST_JOINT = cfgabs("paths", "vst_joint_rds")

# Multicohort input construction (only active when running the multicohort_cancer profile).
IS_MULTICOHORT_PROFILE = profile_name == "multicohort_cancer"
MC_PROFILES    = [str(p) for p in MULTICOHORT_CFG.get("profiles", [])]
MC_OUTDIR_REL  = MULTICOHORT_CFG.get("outdir", "results/unsupervised/multicohort_cancer")
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

# Feature list paths.
# _REL variant → used in rule output: declarations (relative, for portable DAG tracking).
# Non-suffixed absolute variant → used in rule input:, params:, and shell: references.
# Directory no longer encodes a single TOPN since each method may use a different value.
# File naming: genes_top{N}_{method}.txt  (N is method-specific).
FEATURESETS_DIR_REL = os.path.join(UNSUP_REL, "feature_selection_unsupervised", "feature_sets")
FEATURESETS_DIR     = os.path.join(UNSUP,     "feature_selection_unsupervised", "feature_sets")
FEATURELIST_FILES   = [
    os.path.join(FEATURESETS_DIR_REL, f"genes_top{method_topn(m)}_{m}.txt")
    for m in FEATURE_METHODS
    if m != "MX"
]

# Defines relative and absolute roots for agnostic clustering outputs and
# tumour neighbourhood input staging directories.
FEATURE_METHOD_OUTDIR_REL = os.path.join(UNSUP_REL, "feature_selection_unsupervised", "featuresets")
FEATURE_METHOD_OUTDIR     = os.path.join(UNSUP,     "feature_selection_unsupervised", "featuresets")
TUMOUR_NH_INPUT_ROOT_REL  = os.path.join(UNSUP_REL, "tumour_neighbourhoods_input")
TUMOUR_NH_INPUT_ROOT      = os.path.join(UNSUP,     "tumour_neighbourhoods_input")

# Initial pipeline target (Stage 6).  Overridden at end of file after
# Stage 7 variables are defined so that rule all: includes DESeq2 completion.
PIPELINE_TARGET = os.path.join(UNSUP_REL, "tumour_neighbourhoods", "final_consensus_all", "resolved_dsmz_neighbours.tsv")


# =============================================================================
# STAGE 0: FEATURE SELECTION
# =============================================================================

# -----------------------------------------------------------------------------
# RULE: feature_selection_unsupervised
#
# Purpose:
#   Start the analysis by creating gene lists that define different 'views' of the
#   same data. Each feature-selection method picks genes in a different way.
#
# Why this matters:
#   All downstream clustering rules depend on these gene sets. A direction such as
#   'Variance_euc' means: use the Variance-selected genes, then cluster with
#   Euclidean distance.
# -----------------------------------------------------------------------------
rule feature_selection_unsupervised:
    """
    Selects the top N genes from the joint VST expression matrix using multiple
    statistical methods (Variance, MAD, Entropy, PCA loadings, MX, HVG, etc.).
    Each method produces an independent ranked gene list that defines one
    feature direction propagated through all downstream clustering rules.
    The MX gene list is declared explicitly as mx_list so downstream rules
    can reference it as the canonical MX top-500 set without duplicating the path.
    """
    input:
        vst = VST_JOINT,
        cfg_file = CFGFILE_ABS
    output:
        # Relative paths for portable DAG tracking; R script writes via absolute --outdir.
        feature_lists = FEATURELIST_FILES,
        # Canonical MX gene list (MX method, top-500) — named output for downstream rules.
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
# STAGE 1: AGNOSTIC CLUSTERING — CONFIGURATION
# =============================================================================

# Resolves cell line and tumour VST expression matrices (absolute; external inputs).
# Falls back to vst_joint_rds for profiles (e.g. brca, nbl) that do not define
# separate cell_vst_rds / tumour_vst_rds keys.
CELL_VST   = cfgget_path_abs(cfgrel("paths", "vst_joint_rds"), "paths", "cell_vst_rds")
TUMOUR_VST = cfgget_path_abs(cfgrel("paths", "vst_joint_rds"), "paths", "tumour_vst_rds")

if profile_name in ("brca", "nbl", "rbl"):
    rule split_profile_joint_vst:
        input:
            joint = VST_JOINT
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

# Resolves agnostic clustering output root.
# Config override (agnostic_cluster_root) must be a relative path; if absent,
# the default is derived from the profile-scoped UNSUP_REL directory.
_agn_cfg_raw = cfg["paths"].get("agnostic_cluster_root")
if _agn_cfg_raw and not os.path.isabs(_agn_cfg_raw):
    AGN_ROOT_REL = _agn_cfg_raw
else:
    AGN_ROOT_REL = os.path.join(UNSUP_REL, "agnostic_clustering")
AGN_ROOT     = AGN_ROOT_REL                  # Backwards-compatibility alias.
AGN_ROOT_ABS = os.path.join(PIPE_ROOT, AGN_ROOT_REL) if not os.path.isabs(AGN_ROOT_REL) else AGN_ROOT_REL
AGN_BASE_FUN = abspath(cfg["paths"].get("agnostic_base_functions_dir", cfg["paths"].get("base_functions_dir", os.path.join(UNSUP, "base_functions"))))

# Reads agnostic clustering hyperparameters from config with sensible defaults.
AGN_CFG      = cfg.get("agnostic_clustering", {})
AGN_N_PCS    = AGN_CFG.get("n_pcs", 20)
AGN_MAX_K_HC = AGN_CFG.get("max_k_hc", 8)
AGN_KMIN     = AGN_CFG.get("k_min", 2)
AGN_KMAX     = AGN_CFG.get("k_max", 8)
AGN_SEED     = AGN_CFG.get("seed", 42)

# Reads and filters the consensus direction list independently of AGN_DIRECTIONS.
CONS_DIRECTIONS = cfg.get("tumour_neighbourhoods", {}).get("directions", AGN_DIRECTIONS)

# Builds wildcard constraint patterns from configured directions.
AGN_DIRECTION_PATTERN = "|".join(map(re.escape, AGN_DIRECTIONS))
AGN_EUC_DIRECTIONS    = [d for d in AGN_DIRECTIONS if d.endswith("_euc")]
# Falls back to an impossible pattern ('a^') when no Euclidean directions exist,
# preventing k-means rules from matching any wildcard.
AGN_EUC_PATTERN       = "|".join(map(re.escape, AGN_EUC_DIRECTIONS)) if AGN_EUC_DIRECTIONS else r"a^"

CONS_DIRECTION_PATTERN = "|".join(map(re.escape, CONS_DIRECTIONS))
CONS_EUC_DIRECTIONS    = [d for d in CONS_DIRECTIONS if d.endswith("_euc")]
CONS_EUC_PATTERN       = "|".join(map(re.escape, CONS_EUC_DIRECTIONS)) if CONS_EUC_DIRECTIONS else r"a^"
PAN_CANCER_FEATURE_SET_NAME = cfg.get("features", {}).get("pan_cancer_feature_set_name", "PanCancerFeatureSet")
PAN_CANCER_FEATURE_SET_GENE_LIST = cfg.get("features", {}).get("pan_cancer_feature_set_gene_list")
EXTERNAL_GENE_LIST_FEATURES = {PAN_CANCER_FEATURE_SET_NAME} if PAN_CANCER_FEATURE_SET_GENE_LIST else set()


def _resolve_majority_support_threshold():
    """
    Resolves the majority-style support threshold m = max(2, ceil(|R| / 2))
    for the support-threshold consensus cell-line similarity network, where
    |R| is the number of representation-specific graphs for this cohort.

    The value is recomputed from CONS_DIRECTIONS at Snakefile load time. If a
    cohort also declares
        profiles.<cohort>.tumour_neighbourhoods.similarity_consensus_min_support
    in config/config.yaml, the declared value must match the recomputed value.
    """
    nh_cfg = cfg.get("tumour_neighbourhoods", {})
    computed = max(2, (len(CONS_DIRECTIONS) + 1) // 2)
    if "similarity_consensus_min_support" not in nh_cfg:
        return computed
    m = nh_cfg["similarity_consensus_min_support"]
    if not isinstance(m, int) or m < 2:
        raise ValueError(
            "\n[Snakefile] Invalid value for "
            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support: "
            f"{m!r}. The majority-style support threshold must be an integer >= 2.\n"
        )
    if m != computed:
        raise ValueError(
            "\n[Snakefile] Configured support threshold does not match the "
            "configured representation count for profile "
            f"'{profile_name}'.\n"
            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support: "
            f"{m!r}; expected {computed} from |R|={len(CONS_DIRECTIONS)} "
            "(m = max(2, ceil(|R| / 2))).\n"
        )
    return computed


SIMILARITY_CONSENSUS_MIN_SUPPORT = _resolve_majority_support_threshold()

# Defines the six HC and six k-means clustering kind identifiers.
# PCA-prefixed kinds are optionally removed for pan-cancer runs where
# cross-cancer PCA spaces are not comparable.
HC_KINDS = ["pca_hc_cell", "pca_hc_tumour", "pca_hc_cell_tumour", "hc_cell", "hc_tumour", "hc_cell_tumour"]
KM_KINDS = ["pca_kmeans_cell", "pca_kmeans_tumour", "pca_kmeans_cell_tumour", "kmeans_cell", "kmeans_tumour", "kmeans_cell_tumour"]

if DISABLE_PCA_EVERYWHERE:
    HC_KINDS = [k for k in HC_KINDS if not k.startswith("pca_")]
    KM_KINDS = [k for k in KM_KINDS if not k.startswith("pca_")]
    vprint(f"[Snakefile] Filtered PCA kinds: HC_KINDS={HC_KINDS}, KM_KINDS={KM_KINDS}")


def dir_to_feature(direction):
    """Extracts the feature method name from a direction string (e.g., 'Variance_euc' → 'Variance')."""
    if direction.startswith("pam50"):
        return "pam50"
    else:
        if direction.endswith("_euc") or direction.endswith("_corr"):
            return direction.rsplit("_", 1)[0]
        else:
            return direction


def dir_to_dist(direction):
    """Maps a direction suffix to its R-compatible distance method string."""
    return "correlation" if direction.endswith("corr") else "euclidean"


def dir_to_gene_list(direction):
    """
    Returns the **relative** path to the gene list file appropriate for a given
    direction.  The path style MUST match the producer rule
    (feature_selection_unsupervised), which declares its outputs under the
    relative FEATURESETS_DIR_REL.  Using the absolute FEATURESETS_DIR here would
    cause a MissingInputException because Snakemake tracks files by the exact
    path string used in the output: declaration.

    PAM50 gene lists are external resource files (not rule-produced), so they
    use absolute paths via cfgabs().

    Raises ValueError if the direction contains an unknown feature method.
    """
    feature = dir_to_feature(direction)

    if feature == "pam50":
        # PAM50 gene lists are external resources, not produced by a rule.
        # Absolute paths are correct here.
        if "pam50" in cfg and "ensembl_gene_list" in cfg["pam50"]:
            return cfgabs("pam50", "ensembl_gene_list")
        return cfgabs("features", "pam50_ensembl_gene_list")
    elif feature in EXTERNAL_GENE_LIST_FEATURES:
        # External feature-set directions reuse a fixed gene list rather than
        # a feature_selection_unsupervised output.
        return cfgabs("features", "pan_cancer_feature_set_gene_list")
    elif feature in FEATURE_METHODS:
        # Feature-set directions: filename encodes the method-specific top-N.
        # Uses FEATURESETS_DIR_REL to match the producer rule output.
        return os.path.join(FEATURESETS_DIR_REL, f"genes_top{method_topn(feature)}_{feature}.txt")
    else:
        # Unknown feature: fail loudly to prevent silent reuse of a default gene list.
        raise ValueError(
            f"Unknown feature '{feature}' extracted from direction '{direction}'. "
            f"Known features: {sorted(set(FEATURE_METHODS) | {'pam50'} | EXTERNAL_GENE_LIST_FEATURES)}"
        )


def dir_to_geneset(direction):
    """Returns a human-readable gene set label for annotation and reporting purposes."""
    feature = dir_to_feature(direction)
    if feature == "pam50":
        return "PAM50"
    elif feature in EXTERNAL_GENE_LIST_FEATURES:
        return feature
    else:
        return f"{feature}_top{method_topn(feature)}"


def validate_direction_to_gene_list(direction):
    """
    Validates that a direction string correctly resolves to its gene-list file.
    Three checks:
      1. The feature method parsed from the direction is recognized.
      2. The resolved filename encodes the correct method name.
      3. Non-PAM50 paths are relative (matching the producer rule output style),
         preventing the absolute-vs-relative MissingInputException.

    Runs at Snakefile parse time to catch bugs before any job executes.
    """
    feature = dir_to_feature(direction)

    # 1. Check that the extracted feature is recognized.
    known_features = {"pam50"} | set(FEATURE_METHODS) | EXTERNAL_GENE_LIST_FEATURES
    assert feature in known_features, \
        f"Unknown feature '{feature}' extracted from direction '{direction}'. " \
        f"Known features: {sorted(known_features)}"

    # 2. Verify the gene-list file path is well-formed and method-specific.
    gene_list_path = dir_to_gene_list(direction)
    filename = os.path.basename(gene_list_path)

    # For feature-set directions (not pam50), verify filename format.
    if feature not in ("pam50",) and feature not in EXTERNAL_GENE_LIST_FEATURES:
        if not filename.startswith("genes_top"):
            raise AssertionError(
                f"Gene list file '{gene_list_path}' has unexpected name format. "
                f"Expected 'genes_top<N>_<method>.txt' pattern."
            )

    # Verify that the filename method matches the direction-derived feature.
    if feature not in ("pam50",) and feature not in EXTERNAL_GENE_LIST_FEATURES:
        if not filename.endswith(f"_{feature}.txt"):
            raise AssertionError(
                f"Gene list for direction '{direction}' (method '{feature}') "
                f"resolved to file '{filename}', which does not end with '_{feature}.txt'. "
                f"This indicates a misalignment in the direction-to-gene-list mapping."
            )

    # 3. Feature-selection paths must be relative to match the producer rule
    #    (feature_selection_unsupervised) output declarations.
    #    Absolute paths would cause Snakemake MissingInputException.
    if feature not in ("pam50",) and feature not in EXTERNAL_GENE_LIST_FEATURES and os.path.isabs(gene_list_path):
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
for _d in AGN_DIRECTIONS:
    validate_direction_to_gene_list(_d)
vprint(f"[Snakefile] Direction→gene-list validation passed for {len(AGN_DIRECTIONS)} directions")


def agn_outdir_dir(direction, kind):
    """Returns the absolute output directory path for a given direction–kind combination."""
    return os.path.join(AGN_ROOT_ABS, direction, kind)


# =============================================================================
# STAGE 1: AGNOSTIC CLUSTERING — RULES
# =============================================================================

# -----------------------------------------------------------------------------
# RULE: agnostic_cluster_build_inputs
#
# Purpose:
#   For one analysis direction, subset the cell-line and tumour expression matrices
#   to only the genes relevant for that direction.
#
# Output meaning:
#   These filtered matrices are the clean starting inputs for all clustering jobs
#   for that direction.
# -----------------------------------------------------------------------------
rule agnostic_cluster_build_inputs:
    """
    Subsets the cell line and tumour VST expression matrices to the gene list
    defined by the given direction, producing filtered RDS objects used as
    input to all agnostic clustering rules for that direction.
    """
    input:
        cell   = CELL_VST,
        tumour = TUMOUR_VST,
        genes  = lambda wc: dir_to_gene_list(wc.direction)
    output:
        cell_filt   = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour_filt = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    log: os.path.join(LOGROOT, "agnostic_cluster_build_inputs_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p "$(dirname "{output.cell_filt}")" "$(dirname "{output.tumour_filt}")"
        Rscript "{SCRIPTS_DIR}/build_agnostic_direction_mats.R" \
          --cell_rds   "{input.cell}" \
          --tumour_rds "{input.tumour}" \
          --genes      "{input.genes}" \
          --out_cell   "{output.cell_filt}" \
          --out_tumour "{output.tumour_filt}" \
          > "{log}" 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: agnostic_cluster_pca_hc_cell
#
# Purpose:
#   Cluster only the cell-line matrix using hierarchical clustering after PCA.
#
# Beginner note:
#   'pca_hc' means the data are first reduced to principal components, then
#   hierarchical clustering is performed.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_hc_cell:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the cell line expression matrix. Optimal k is selected by silhouette
    or gap statistic within [2, max_k]. Produces a cluster assignment RDS
    consumed by downstream tumour neighbourhood and consensus rules.
    """
    input:
        cell = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_hc_cell", "pca_hc_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_hc_cell"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_hc_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_pca_hc_tumour
#
# Same idea as the previous rule, but now only tumour samples are clustered.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_hc_tumour:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the tumour expression matrix. Mirrors agnostic_cluster_pca_hc_cell but
    operates on patient tumour samples independently to capture tumour-intrinsic
    structure.
    """
    input:
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_hc_tumour", "pca_hc_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_hc_tumour"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_hc_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_pca_hc_joint
#
# Purpose:
#   Cluster the combined cell-line + tumour matrix together after PCA.
#
# Why this is useful:
#   Joint clustering lets you see whether cell lines and tumours occupy similar
#   regions of transcriptomic space.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_hc_joint:
    """
    Applies PCA dimensionality reduction followed by hierarchical clustering
    to the concatenated cell line and tumour expression matrix. Joint embedding
    captures shared transcriptional structure across both sample types,
    enabling direct cross-space cluster comparison.
    """
    input:
        cell   = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_hc_cell_tumour", "pca_hc_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_hc_cell_tumour"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_hc_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_pca_kmeans_cell
#
# Purpose:
#   Run k-means on PCA-reduced cell-line data for one direction.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_kmeans_cell:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the cell line expression matrix. Constrained to Euclidean directions only,
    as k-means minimises Euclidean inertia and is undefined under correlation distance.
    """
    input:
        cell = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_kmeans_cell", "pca_kmeans_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_kmeans_cell"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_kmeans_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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
# RULE: agnostic_cluster_pca_kmeans_tumour
#
# Purpose:
#   Run k-means on PCA-reduced tumour data for one direction.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_kmeans_tumour:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the tumour expression matrix. Euclidean directions only; mirrors
    agnostic_cluster_pca_kmeans_cell but operates on patient samples independently.
    """
    input:
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_kmeans_tumour", "pca_kmeans_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_kmeans_tumour"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_kmeans_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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
# RULE: agnostic_cluster_pca_kmeans_joint
#
# Purpose:
#   Run k-means on the PCA-reduced joint cell-line + tumour matrix.
# -----------------------------------------------------------------------------
rule agnostic_cluster_pca_kmeans_joint:
    """
    Applies PCA dimensionality reduction followed by k-means clustering to
    the joint cell line and tumour expression matrix. Euclidean directions only.
    Joint embedding allows k-means to assign both sample types to shared clusters.
    """
    input:
        cell   = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "pca_kmeans_cell_tumour", "pca_kmeans_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "pca_kmeans_cell_tumour"),
        base_dir = AGN_BASE_FUN,
        n_pcs    = AGN_N_PCS,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_pca_kmeans_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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
# RULE: agnostic_cluster_hc_cell
#
# Purpose:
#   Run hierarchical clustering directly on the filtered cell-line matrix
#   without PCA. This keeps the analysis in the original selected-gene space.
# -----------------------------------------------------------------------------
rule agnostic_cluster_hc_cell:
    """
    Applies hierarchical clustering directly to the raw gene-subsetted cell line
    expression matrix without PCA pre-processing. Distance metric is determined
    by the direction suffix (Euclidean or correlation). Provides a non-PCA
    baseline for comparison with agnostic_cluster_pca_hc_cell under the same
    feature set.
    """
    input:
        cell = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "hc_cell", "hc_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "hc_cell"),
        base_dir = AGN_BASE_FUN,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_hc_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_hc_tumour
#
# Purpose:
#   Run hierarchical clustering directly on the filtered tumour matrix.
# -----------------------------------------------------------------------------
rule agnostic_cluster_hc_tumour:
    """
    Applies hierarchical clustering directly to the raw gene-subsetted tumour
    expression matrix without PCA pre-processing. Provides a non-PCA baseline
    for tumour-intrinsic cluster structure under each feature direction.
    """
    input:
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "hc_tumour", "hc_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "hc_tumour"),
        base_dir = AGN_BASE_FUN,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_hc_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_hc_joint
#
# Purpose:
#   Run hierarchical clustering directly on the filtered joint matrix.
# -----------------------------------------------------------------------------
rule agnostic_cluster_hc_joint:
    """
    Applies hierarchical clustering to the joint cell line and tumour expression
    matrix without PCA pre-processing. Supports both Euclidean and correlation
    distances, enabling direct assessment of co-clustering fidelity across
    sample types in the raw feature space.
    """
    input:
        cell   = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "hc_cell_tumour", "hc_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "hc_cell_tumour"),
        base_dir = AGN_BASE_FUN,
        max_k    = AGN_MAX_K_HC,
        dist     = lambda wc: dir_to_dist(wc.direction),
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_hc_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
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
# RULE: agnostic_cluster_kmeans_cell
#
# Purpose:
#   Run k-means directly on the filtered cell-line matrix without PCA.
# -----------------------------------------------------------------------------
rule agnostic_cluster_kmeans_cell:
    """
    Applies k-means clustering directly to the raw cell line expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    Serves as a non-PCA counterpart to agnostic_cluster_pca_kmeans_cell for
    sensitivity assessment.
    """
    input:
        cell = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "kmeans_cell", "kmeans_cell_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "kmeans_cell"),
        base_dir = AGN_BASE_FUN,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_kmeans_cell.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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
# RULE: agnostic_cluster_kmeans_tumour
#
# Purpose:
#   Run k-means directly on the filtered tumour matrix without PCA.
# -----------------------------------------------------------------------------
rule agnostic_cluster_kmeans_tumour:
    """
    Applies k-means clustering directly to the raw tumour expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    """
    input:
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "kmeans_tumour", "kmeans_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "kmeans_tumour"),
        base_dir = AGN_BASE_FUN,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_kmeans_tumour.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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
# RULE: agnostic_cluster_kmeans_joint
#
# Purpose:
#   Run k-means directly on the filtered joint cell-line + tumour matrix.
# -----------------------------------------------------------------------------
rule agnostic_cluster_kmeans_joint:
    """
    Applies k-means clustering to the joint cell line and tumour expression matrix
    without PCA pre-processing. Constrained to Euclidean directions only.
    Joint clustering provides a raw-space baseline for cross-type co-assignment.
    """
    input:
        cell   = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "cell_expr.rds"),
        tumour = os.path.join(AGN_ROOT_REL, "{direction}", "inputs", "tumour_expr.rds")
    output:
        clusters = os.path.join(AGN_ROOT_REL, "{direction}", "kmeans_cell_tumour", "kmeans_cell_tumour_clusters_optimal.rds")
    params:
        outdir   = lambda wc: agn_outdir_dir(wc.direction, "kmeans_cell_tumour"),
        base_dir = AGN_BASE_FUN,
        kmin     = AGN_KMIN,
        kmax     = AGN_KMAX,
        seed     = AGN_SEED
    log: os.path.join(LOGROOT, "agnostic_cluster_{direction}_kmeans_joint.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_EUC_PATTERN
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


# -----------------------------------------------------------------------------
# RULE: agnostic_cluster_all
#
# Purpose:
#   This is a collector rule. It does not perform a new analysis itself.
#   Instead, it asks Snakemake to build every clustering result that belongs to
#   the active direction set and clustering kinds.
#
# Why collector rules are useful:
#   They provide one clean target that stands for 'finish this whole stage'.
# -----------------------------------------------------------------------------
rule agnostic_cluster_all:
    """
    Convenience aggregation rule that collects all HC and k-means cluster outputs
    across every direction. Allows the full Stage 1 suite to be run independently
    via 'snakemake agnostic_cluster_all' without proceeding to consensus.
    """
    input:
        expand(os.path.join(AGN_ROOT_REL, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=AGN_DIRECTIONS, kind=HC_KINDS),
        expand(os.path.join(AGN_ROOT_REL, "{direction}", "{kind}", "{kind}_clusters_optimal.rds"), direction=AGN_EUC_DIRECTIONS, kind=KM_KINDS)


# =============================================================================
# STAGE 2: CONSENSUS CLUSTERING
# =============================================================================

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
# Purpose:
#   Take the clustering results from the previous stage and summarize them into
#   consensus-clustering outputs for one direction.
#
# Intuition:
#   Instead of trusting a single clustering view, consensus clustering measures
#   how stable the grouping is across related clustering runs.
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
        alg    = lambda wc: kind_to_alg(wc.kind)
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
          --outdir {params.outdir} \
          --cluster_rds {output.cluster_rds} > {log} 2>&1
        '''


# -----------------------------------------------------------------------------
# RULE: consensus_cluster_all
#
# Purpose:
#   Another collector rule. It makes Snakemake build all consensus outputs across
#   the configured directions for the active profile.
# -----------------------------------------------------------------------------
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
# STAGE 3: TUMOUR NEIGHBOURHOOD INPUT CONSTRUCTION
# =============================================================================

TUMOUR_NH_ROOT     = os.path.join(UNSUP_REL, "tumour_neighbourhoods")
# Absolute counterpart used in shell find commands to ensure correctness on HPC
# nodes where the job working directory may differ from PIPE_ROOT.
TUMOUR_NH_ROOT_ABS = os.path.join(UNSUP,     "tumour_neighbourhoods")

# Resolves PAM50 neighbourhood input paths only when a configured direction uses PAM50.
if HAS_PAM50_DIRECTIONS:
    TUMOUR_NH_EXPR_PAM50    = os.path.join(TUMOUR_NH_INPUT_ROOT, "expr_pam50.rds")
    TUMOUR_NH_MAP_RDS_PAM50 = os.path.join(TUMOUR_NH_INPUT_ROOT, "cell_line_to_original_sample_id_pam50.rds")
else:
    TUMOUR_NH_EXPR_PAM50    = None
    TUMOUR_NH_MAP_RDS_PAM50 = None

# This section shifts from clustering to tumour-neighbourhood analysis.
# Upstream clustering results are now treated as evidence that helps define
# neighbourhood relationships between cell lines and tumours.


def nh_expr_rds_pam50():
    """Returns the PAM50 expression matrix path. Only valid for PAM50 directions."""
    if TUMOUR_NH_EXPR_PAM50:
        return TUMOUR_NH_EXPR_PAM50
    raise ValueError("PAM50 expression path missing for configured PAM50 direction")


def nh_map_rds_pam50():
    """Returns the PAM50 cell-line mapping path. Only valid for PAM50 directions."""
    if TUMOUR_NH_MAP_RDS_PAM50:
        return TUMOUR_NH_MAP_RDS_PAM50
    raise ValueError("PAM50 map path missing for configured PAM50 direction")


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


def is_pam50_direction(direction):
    return dir_to_feature(direction) == "pam50"


def nh_staged_inputs_for_direction(direction):
    if is_pam50_direction(direction):
        return [nh_expr_rds_pam50(), nh_map_rds_pam50()]
    return []


def nh_expr_arg(direction):
    return f"--expr-rds {shlex.quote(nh_expr_rds_pam50())}" if is_pam50_direction(direction) else ""


def nh_map_arg(direction):
    return f"--mapping-rds {shlex.quote(nh_map_rds_pam50())}" if is_pam50_direction(direction) else ""


def nh_consensus_dependency(direction):
    sentinel = ".tumour_neighbourhoods_km_done" if direction.endswith("_euc") else ".tumour_neighbourhoods_done"
    return os.path.join(TUMOUR_NH_ROOT, direction, sentinel)


if HAS_PAM50_DIRECTIONS:
    # Resolves PAM50 expression matrix paths before instantiating the optional
    # build rule.
    # Prefers profile-specific keys (tcga_brca_pam50_expr, dsmz_bcc_pam50_expr);
    # falls back to generic keys. If source paths are absent, the staged PAM50
    # files are treated as external inputs.
    _pam50_paths = cfg.get("paths", {})
    if "tcga_brca_pam50_expr" in _pam50_paths:
        _TCGA_PAM50_PATH = cfgabs("paths", "tcga_brca_pam50_expr")
    elif "tcga_pam50_expr" in _pam50_paths:
        _TCGA_PAM50_PATH = cfgabs("paths", "tcga_pam50_expr")
    else:
        _TCGA_PAM50_PATH = None
    if "dsmz_bcc_pam50_expr" in _pam50_paths:
        _DSMZ_PAM50_PATH = cfgabs("paths", "dsmz_bcc_pam50_expr")
    elif "dsmz_pam50_expr" in _pam50_paths:
        _DSMZ_PAM50_PATH = cfgabs("paths", "dsmz_pam50_expr")
    else:
        _DSMZ_PAM50_PATH = None

if HAS_PAM50_DIRECTIONS and _TCGA_PAM50_PATH and _DSMZ_PAM50_PATH:
    rule tumour_nh_build_inputs_pam50:
        """
        Constructs the joint PAM50 expression matrix and cell-line-to-sample mapping
        file from pre-computed PAM50-scored TCGA and DSMZ expression tables.
        Only instantiated when the configured directions include PAM50.
        """
        input:
            tcga_pam50 = _TCGA_PAM50_PATH,
            dsmz_pam50 = _DSMZ_PAM50_PATH,
            dsmz_meta  = cfgget_path_abs(cfgget_path_rel("data/dsmz/DSMZ_metadata.csv", "paths", "meta_tsv"), "paths", "dsmz_meta_csv")
        output:
            expr_rds = TUMOUR_NH_EXPR_PAM50,
            map_rds  = TUMOUR_NH_MAP_RDS_PAM50
        log: os.path.join(LOGROOT, "tumour_nh_build_inputs_pam50.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p $(dirname {output.expr_rds}) $(dirname {output.map_rds})
            Rscript {SCRIPTS_DIR}/build_dsmz_tcga_pam50matrix.R --config {CFGFILE_ABS} --profile "{profile_name}" > {log} 2>&1
            test -s {output.expr_rds} || (echo "ERROR: missing {output.expr_rds}" >&2; exit 1)
            test -s {output.map_rds}  || (echo "ERROR: missing {output.map_rds}"  >&2; exit 1)
            '''



# =============================================================================
# STAGE 3: TUMOUR NEIGHBOURHOOD ANALYSIS
# =============================================================================

# -----------------------------------------------------------------------------
# RULE: tumour_nh_hc
#
# Purpose:
#   Run tumour-neighbourhood analysis for all directions using hierarchical
#   clustering outputs. PAM50 directions receive pre-built staged expression
#   and mapping files via params; feature-set directions resolve from config.
# -----------------------------------------------------------------------------
rule tumour_nh_hc:
    """
    Computes per-cell-line tumour neighbourhoods for all directions using
    hierarchical clustering outputs. Direction-driven: PAM50 directions pass
    staged files via params; non-PAM50 directions derive inputs from config.
    """
    input:
        cfg_file    = CFGFILE_ABS,
        cluster_rds = lambda wc: expand(
            os.path.join(CONS_ROOT, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=CONS_HC_KINDS
        ),
        staged = lambda wc: nh_staged_inputs_for_direction(wc.direction)
    output:
        done = touch(os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_done"))
    params:
        script   = os.path.join(BASE, "scripts", "compute_tumour_neighbourhoods.R"),
        expr_arg = lambda wc: nh_expr_arg(wc.direction),
        map_arg  = lambda wc: nh_map_arg(wc.direction),
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
          {params.expr_arg} \
          {params.map_arg} \
          --unsup-root {UNSUP} \
          --cluster-family hc \
          > {log} 2>&1
        test $(find {TUMOUR_NH_ROOT_ABS}/{wildcards.direction} -maxdepth 2 -type f -name "Top_m_long_*.csv" | wc -l) -gt 0 || (echo "ERROR: No Top_m_long_*.csv" >&2 && exit 1)
        touch {output.done}
        '''


rule tumour_nh_km:
    """
    k-means neighbourhood pass for all Euclidean directions. PAM50 directions
    pass staged files via params; feature-set directions derive inputs from config.
    """
    input:
        hc_done     = os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_done"),
        cfg_file    = CFGFILE_ABS,
        cluster_rds = lambda wc: expand(
            os.path.join(CONS_ROOT, wc.direction, "{kind}", "{kind}_clusters_optimal.rds"),
            kind=CONS_KM_KINDS
        ),
        staged = lambda wc: nh_staged_inputs_for_direction(wc.direction)
    output:
        done = touch(os.path.join(TUMOUR_NH_ROOT, "{direction}", ".tumour_neighbourhoods_km_done"))
    params:
        script   = os.path.join(BASE, "scripts", "compute_tumour_neighbourhoods.R"),
        expr_arg = lambda wc: nh_expr_arg(wc.direction),
        map_arg  = lambda wc: nh_map_arg(wc.direction),
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
          {params.expr_arg} \
          {params.map_arg} \
          --unsup-root {UNSUP} \
          --cluster-family km \
          > {log} 2>&1
        test $(find {TUMOUR_NH_ROOT_ABS}/{wildcards.direction} -maxdepth 2 -type f -name "Top_m_long_*_KM_*.csv" | wc -l) -gt 0 || (echo "ERROR: No KM-derived Top_m_long_*_KM_*.csv for {wildcards.direction}" >&2 && exit 1)
        touch {output.done}
        '''


# =============================================================================
# STAGE 4 PROBABILISTIC NEIGHBOURHOOD CONSENSUS
# =============================================================================

rule tumour_nh_consensus:
    """
    Probabilistic neighbourhood consensus for every configured direction.
    Euclidean directions depend on the KM sentinel so HC and KM neighbourhood
    files are present; correlation directions depend on HC only.
    """
    input:
        done = lambda wc: nh_consensus_dependency(wc.direction)
    output:
        consensus_rds = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "Final_consensus_tumour_neighbourhoods_{direction}.rds"),
        consensus_tsv = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "Final_consensus_tumour_neighbourhoods_{direction}.tsv")
    params:
        script = os.path.join(BASE, "scripts", "tumour_neighbourhood_p_consensus.R"),
        config = os.path.join(BASE, "config", "config.yaml")
    log: os.path.join(LOGROOT, "tumour_nh_consensus_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=CONS_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.consensus_rds})
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          --out_rds {output.consensus_rds} \
          --out_tsv {output.consensus_tsv} \
          > {log} 2>&1
        '''


# =============================================================================
# STAGE 5: Cell-line Similarity Graph
# =============================================================================

rule cell_line_similarity_graph:
    input:
        consensus_rds = lambda wc: nh_final_consensus_rds(wc.direction)
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
        graph_minimal = os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus", "plots", "Fig_cell_line_similarity_graph_minimal_{direction}.pdf")
    params:
        script = os.path.join(BASE, "scripts", "compute_cell_line_similarity.R"),
        config = os.path.join(BASE, "config", "config.yaml")
    log: os.path.join(LOGROOT, "cell_line_similarity_graph_{direction}.log")
    conda: CONDA_ENV_R
    wildcard_constraints:
        direction=AGN_DIRECTION_PATTERN
    shell:
        r'''
        mkdir -p $(dirname {output.sim_mat_rds})
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --direction {wildcards.direction} \
          > {log} 2>&1
        '''

P_CONS_ALL_DIR = os.path.join(TUMOUR_NH_ROOT, "final_consensus_all")
P_CONS_PLOTS_DIR = os.path.join(P_CONS_ALL_DIR, "plots")
P_CONS_THESIS_COHORTS = ["brca", "nbl", "rbl"]
P_CONS_THESIS_COHORT_LABELS = {
    "brca": "BRCA",
    "nbl": "NBL",
    "rbl": "RBL",
}
P_CONS_THESIS_FIG_DIR = abspath(
    os.path.join(
        config.get("multicohort_cancer", {}).get("outdir", "results/unsupervised/multicohort_cancer"),
        "thesis_figures",
    )
)
P_CONS_DIRECTION_DUMBBELL_PREFIX = os.path.join(
    P_CONS_THESIS_FIG_DIR,
    "Fig_p_consensus_direction_comparison_dumbbell",
)

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

_COHORT_UPPER = profile_name.upper()
# Naming-prefix selector: patient-referenced cohort graphs (BRCA/NBL/RBL etc.) are
# the cohort-specific outputs derived from tumour/patient-referenced similarity;
# the multicohort_cancer profile produces the pan-cancer integrated graph and
# carries the pan_cancer prefix instead. Keeping these two namings distinct is a
# hard requirement of the patient_referenced naming pass.
_PR_PREFIX = "pan_cancer_" if IS_MULTICOHORT_PROFILE else "patient_referenced_"
P_CONS_THESIS_RESOLUTION_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_{_PR_PREFIX}resolved_cell_line_neighbourhood_graph",
)
P_CONS_THESIS_CONSENSUS_PREFIX = os.path.join(
    P_CONS_PLOTS_DIR,
    f"Fig_{_COHORT_UPPER}_{_PR_PREFIX}support_threshold_consensus_cell_line_similarity_network",
)
P_CONS_SHORTNAMES_TSV = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}cell_line_display_names.tsv")
P_CONS_RESOLVED_EDGES_TSV = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_edges.tsv")
P_CONS_RESOLVED_NODE_STATS_TSV = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_graph_node_stats.tsv")
P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}support_threshold_consensus_cell_line_similarity_edges.tsv")
P_CONS_ANCHOR_AUDIT_TSV            = os.path.join(P_CONS_PLOTS_DIR, f"{_PR_PREFIX}resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv")
P_CONS_DSMZ_META_CSV = cfgget_path_abs(
    cfgget_path_rel("data/dsmz/DSMZ_metadata.csv", "paths", "meta_tsv"),
    "paths", "dsmz_meta_csv",
)
P_CONS_RESOLVED_COMPONENT_PANELS_PREFIX = P_CONS_THESIS_RESOLUTION_PREFIX + "_component_panels"
P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX = P_CONS_THESIS_CONSENSUS_PREFIX + "_component_panels"
P_CONS_RESOLVED_INTERACTIVE_HTML = P_CONS_THESIS_RESOLUTION_PREFIX + "_interactive.html"
P_CONS_SUPPORT_INTERACTIVE_HTML = P_CONS_THESIS_CONSENSUS_PREFIX + "_interactive.html"
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
        P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".pdf",
        P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".png",
        P_CONS_SUPPORT_COMPONENT_PANELS_PREFIX + ".svg",
        P_CONS_SUPPORT_INTERACTIVE_HTML,
        P_CONS_SUPPORT_COMPONENT_SUMMARY_TSV,
        P_CONS_SUPPORT_NODE_LABELS_TSV,
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
P_CONS_LEIDEN_FIG_PREFIX           = os.path.join(P_CONS_PLOTS_DIR, "Fig_MULTICOHORT_CANCER_lineage_community")
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

rule summarize_p_consensus_all:
    input:
        consensus_rds = [nh_final_consensus_rds(d) for d in AGN_DIRECTIONS],
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
        config = os.path.join(BASE, "config", "config.yaml"),
        threshold = lambda wc: cfg.get("tumour_neighbourhoods", {}).get("p_consensus_threshold", 0.7)
    log: os.path.join(LOGROOT, "summarize_p_consensus_all.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p {P_CONS_ALL_DIR}
        Rscript {params.script} \
          --config {params.config} \
          --profile "{profile_name}" \
          --threshold {params.threshold} \
          > {log} 2>&1
        '''


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
        threshold = lambda wc: cfg.get("tumour_neighbourhoods", {}).get("p_consensus_threshold", 0.7)
    log: os.path.join(LOGROOT, "plot_p_consensus_direction_comparison_dumbbell.log")
    conda: CONDA_ENV_R
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
        study_design = STUDY_DESIGN_FILE
    log: os.path.join(LOGROOT, "materialize_study_design.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.question_txt})
        Rscript {params.script}           --config {params.config}           --study-design {params.study_design}           --profile "{profile_name}"           --out-question {output.question_txt}           --out-cohorts {output.cohort_manifest}           --out-labels {output.labels_manifest}           --out-inference {output.inference_manifest}           --out-endpoints {output.endpoint_manifest}           > {log} 2>&1
        '''


rule model_selection_summary:
    input:
        cfg = CFGFILE_ABS,
        ranked_best = os.path.join(P_CONS_ALL_DIR, "p_consensus_best_cell_lines_ranked.tsv"),
        long_tbl = os.path.join(P_CONS_ALL_DIR, "p_consensus_cellline_direction_summary.long.tsv"),
        winning_direction = os.path.join(P_CONS_ALL_DIR, "winning_direction.txt"),
        graph_nodes = expand(
            os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                         "cell_line_similarity_graph_node_annotations_{direction}.tsv"),
            direction=AGN_DIRECTIONS
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
        config = CFGFILE_ABS
    log: os.path.join(LOGROOT, "neighbourhood_permutation_validation.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.summary_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --out-tsv {output.summary_tsv}           --out-plot {output.plot_pdf}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''


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
        config = CFGFILE_ABS
    log: os.path.join(LOGROOT, "random_baseline_comparison.log")
    conda: CONDA_ENV_R
    shell:
        r'''
        mkdir -p $(dirname {output.summary_tsv})
        Rscript {params.script}           --config {params.config}           --profile "{profile_name}"           --out-tsv {output.summary_tsv}           --out-plot {output.plot_pdf}           --out-notes {output.notes_txt}           > {log} 2>&1
        '''


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


# =============================================================================
# STAGE 6: GRAPH-BASED CONSENSUS NEIGHBOUR RESOLUTION
# =============================================================================

rule resolve_dsmz_graph_neighbours:
    """
    Resolves per-cell-line DSMZ neighbours by intersecting the global best-direction
    graph with each cell line's winner-direction graph. Consumes the cross-direction
    consensus summary outputs and all per-direction DSMZ similarity graph edge files
    produced by cell_line_similarity_graph. The intersection strategy retains only
    neighbours present in both the global best direction and the cell-line-specific
    winner direction, providing a robust estimate of stable similarity relationships.
    Uses CONS_DIRECTIONS (not AGN_DIRECTIONS) so only configured consensus directions
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
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv")
    params:
        script = os.path.join(BASE, "scripts", "resolve_dsmz_graph_neighbours.R"),
        config = os.path.join(BASE, "config", "config.yaml")
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
          --output_tsv {output.resolved_tsv} \
          > {log} 2>&1
        test -s {output.resolved_tsv} || (echo "ERROR: missing {output.resolved_tsv}" >&2; exit 1)
        '''


rule plot_patient_referenced_resolved_cell_line_neighbourhood_graph:
    """
    Thesis figure for the resolved DSMZ similarity network after graph-based
    neighbour resolution. The helper also writes short-name, edge-list, and
    node-stat TSVs into the plot directory for figure provenance.
    """
    input:
        resolved_tsv = os.path.join(P_CONS_ALL_DIR, "resolved_dsmz_neighbours.tsv"),
        script = os.path.join(BASE, "scripts", "visualize_resolved_dsmz_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        png = P_CONS_THESIS_RESOLUTION_PREFIX + ".png",
        pdf = P_CONS_THESIS_RESOLUTION_PREFIX + ".pdf",
        svg = P_CONS_THESIS_RESOLUTION_PREFIX + ".svg",
        shortnames = P_CONS_SHORTNAMES_TSV,
        edges = P_CONS_RESOLVED_EDGES_TSV,
        node_stats = P_CONS_RESOLVED_NODE_STATS_TSV,
        anchor_audit = P_CONS_ANCHOR_AUDIT_TSV
    params:
        out_prefix = P_CONS_THESIS_RESOLUTION_PREFIX,
        label = profile_name.upper(),
        dense_component_min_size = 8,
        dense_component_min_density = 0.30,
        dense_component_expand_factor = 2.00,
        resolved_component_label_gap = 2.0,
        resolved_grid_cols = 2,
        resolved_grid_rows = 2,
        resolved_cell_width = 96.0,
        resolved_cell_height = 60.0,
        resolved_row_gap = 16.0,
        resolved_col_gap = 18.0,
        resolved_legend_width = 68.0,
        resolved_isolate_height = 18.0,
        resolved_isolate_spacing = 18.0,
        resolved_isolate_panel_padding_x = 12.0,
        resolved_isolate_panel_padding_y = 1.4,
        resolved_isolate_label_band_frac = 0.34,
        resolved_isolate_max_per_row = 6,
        resolved_isolate_panel_width_frac = 0.62,
        resolved_footnote_height = 7.0,
        resolved_component_cell_fill_x = 0.78,
        resolved_component_cell_fill_y = 0.72,
        resolved_dense_component_cell_fill_x = 0.92,
        resolved_dense_component_cell_fill_y = 0.84,
        resolved_c2_component_cell_fill_x = 0.88,
        resolved_c2_component_cell_fill_y = 0.80,
        c1_x_expand_factor = 3.00,
        c1_y_expand_factor = 3.00,
        c2_x_expand_factor = 1.35,
        c2_y_expand_factor = 1.05,
        isolate_gap = 8.0,
        isolate_box_height_factor = 2.00,
        isolate_label_y_frac = 0.010,
        isolate_label_left_pad = 14.0,
        legend_alignment_mode = "component-band"
    log: os.path.join(LOGROOT, "plot_patient_referenced_resolved_cell_line_neighbourhood_graph.log")
    conda: CONDA_ENV_PY
    shell:
        r'''
        mkdir -p {P_CONS_PLOTS_DIR}
        python {input.script} \
          {input.resolved_tsv} \
          {params.out_prefix} \
          {params.label} \
          --dense-component-min-size {params.dense_component_min_size} \
          --dense-component-min-density {params.dense_component_min_density} \
          --dense-component-expand-factor {params.dense_component_expand_factor} \
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
          --legend-alignment-mode {params.legend_alignment_mode} \
          > {log} 2>&1
        test -s {output.pdf} || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
        '''


rule build_patient_referenced_support_threshold_consensus_cell_line_similarity_network:
    """
    Aggregates per-direction DSMZ similarity graph edge files into a single
    consensus edge table for the thesis consensus-network figure.
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
        # Majority-style support threshold m = max(2, ceil(|R| / 2)); resolved
        # at Snakefile load time from the configured representation list; any
        # explicit tumour_neighbourhoods.similarity_consensus_min_support value
        # must match the recomputed value (see
        # _resolve_majority_support_threshold above).
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


rule plot_patient_referenced_support_threshold_consensus_cell_line_similarity_network:
    """
    Thesis figure for the consensus DSMZ similarity network aggregated across
    configured similarity-network directions.
    """
    input:
        edges = P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV,
        shortnames = P_CONS_SHORTNAMES_TSV,
        script = os.path.join(BASE, "scripts", "plot_consensus_graph.py"),
        style = os.path.join(BASE, "scripts", "graph_plot_style.py")
    output:
        png = P_CONS_THESIS_CONSENSUS_PREFIX + ".png",
        pdf = P_CONS_THESIS_CONSENSUS_PREFIX + ".pdf",
        svg = P_CONS_THESIS_CONSENSUS_PREFIX + ".svg"
    params:
        out_prefix = P_CONS_THESIS_CONSENSUS_PREFIX,
        label = profile_name.upper(),
        component_label_gap = SUPPORT_THRESHOLD_COMPONENT_LABEL_GAP,
        nbl_component_gap_fraction = SUPPORT_THRESHOLD_NBL_COMPONENT_GAP_FRACTION,
        nbl_component_min_gap = SUPPORT_THRESHOLD_NBL_COMPONENT_MIN_GAP,
        legend_anchor_x = SUPPORT_THRESHOLD_LEGEND_ANCHOR_X,
        legend_anchor_y_arg = SUPPORT_THRESHOLD_LEGEND_ANCHOR_Y_ARG,
        legend_align_mode = SUPPORT_THRESHOLD_LEGEND_ALIGN_MODE,
        legend_gap_frac = SUPPORT_THRESHOLD_LEGEND_GAP_FRAC,
        isolate_label_y_frac = SUPPORT_THRESHOLD_ISOLATE_LABEL_Y_FRAC,
        isolate_box_height_factor = SUPPORT_THRESHOLD_ISOLATE_BOX_HEIGHT_FACTOR,
        isolate_gap = SUPPORT_THRESHOLD_ISOLATE_GAP,
        isolate_label_left_pad = 28.0 if _COHORT_UPPER == "BRCA" else 0.0,
        isolate_spacing = SUPPORT_THRESHOLD_ISOLATE_SPACING,
        isolate_panel_padding_x = SUPPORT_THRESHOLD_ISOLATE_PANEL_PADDING_X,
        isolate_panel_padding_y = SUPPORT_THRESHOLD_ISOLATE_PANEL_PADDING_Y,
        isolate_label_band_frac = SUPPORT_THRESHOLD_ISOLATE_LABEL_BAND_FRAC,
        isolate_max_per_row = SUPPORT_THRESHOLD_ISOLATE_MAX_PER_ROW,
        isolate_panel_width_frac = SUPPORT_THRESHOLD_ISOLATE_PANEL_WIDTH_FRAC,
        dense_component_min_size = SUPPORT_THRESHOLD_DENSE_COMPONENT_MIN_SIZE,
        dense_component_min_density = SUPPORT_THRESHOLD_DENSE_COMPONENT_MIN_DENSITY,
        dense_component_expand_factor = SUPPORT_THRESHOLD_DENSE_COMPONENT_EXPAND_FACTOR,
        support_grid_cols = SUPPORT_THRESHOLD_COMPONENT_GRID_NCOLS,
        support_grid_rows = SUPPORT_THRESHOLD_COMPONENT_GRID_NROWS,
        support_cell_width = SUPPORT_THRESHOLD_COMPONENT_CELL_WIDTH,
        support_cell_height = SUPPORT_THRESHOLD_COMPONENT_CELL_HEIGHT,
        support_row_gap = SUPPORT_THRESHOLD_COMPONENT_ROW_GAP,
        support_col_gap = SUPPORT_THRESHOLD_COMPONENT_COL_GAP,
        support_legend_width = SUPPORT_THRESHOLD_LEGEND_WIDTH,
        support_isolate_height = SUPPORT_THRESHOLD_ISOLATE_REGION_HEIGHT,
        support_footnote_height = SUPPORT_THRESHOLD_FOOTNOTE_REGION_HEIGHT,
        support_component_cell_fill_x = SUPPORT_THRESHOLD_COMPONENT_CELL_FILL_X,
        support_component_cell_fill_y = SUPPORT_THRESHOLD_COMPONENT_CELL_FILL_Y,
        support_dense_component_cell_fill_x = SUPPORT_THRESHOLD_DENSE_COMPONENT_CELL_FILL_X,
        support_dense_component_cell_fill_y = SUPPORT_THRESHOLD_DENSE_COMPONENT_CELL_FILL_Y,
        support_c2_component_cell_fill_x = SUPPORT_THRESHOLD_C2_COMPONENT_CELL_FILL_X,
        support_c2_component_cell_fill_y = SUPPORT_THRESHOLD_C2_COMPONENT_CELL_FILL_Y,
        c1_x_expand_factor = SUPPORT_THRESHOLD_C1_X_EXPAND_FACTOR,
        c1_y_expand_factor = SUPPORT_THRESHOLD_C1_Y_EXPAND_FACTOR,
        c2_x_expand_factor = SUPPORT_THRESHOLD_C2_X_EXPAND_FACTOR,
        c2_y_expand_factor = SUPPORT_THRESHOLD_C2_Y_EXPAND_FACTOR
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
          --component-label-gap {params.component_label_gap} \
          --nbl-component-gap-fraction {params.nbl_component_gap_fraction} \
          --nbl-component-min-gap {params.nbl_component_min_gap} \
          --legend-anchor-x {params.legend_anchor_x} \
          {params.legend_anchor_y_arg} \
          --legend-align-mode {params.legend_align_mode} \
          --legend-gap-frac {params.legend_gap_frac} \
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
        '''


if IS_MULTICOHORT_PROFILE:

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
            title_label = "Pan-cancer resolved cell-line neighbourhood graph"
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
            title_label = "Pan-cancer support-threshold consensus cell-line similarity network"
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


rule community_stability_analysis:
    """
    Cross-direction cell-line graph community stability diagnostic.

    Leiden/Louvain assignments are direction-specific diagnostics. This rule
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


# =============================================================================
# STAGE 6b: MULTICOHORT_CANCER LEIDEN COMMUNITY DETECTION
# =============================================================================

rule compute_multicohort_cancer_communities:
    """
    Unweighted Leiden community detection on the final resolved MULTICOHORT_CANCER
    DSMZ cell-line graph (56 nodes, 155 edges, 9 components).
    Produces per-node community assignments, community summary, modularity table,
    layout coordinates, and a two-panel lineage/community PDF+PNG figure.
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
        helper_annots = os.path.join(
            TUMOUR_NH_ROOT, "Entropy_corr", "final_consensus",
            "cell_line_similarity_graph_node_annotations_Entropy_corr.tsv"
        ),
        seed          = 42,
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
          --helper-node-annotations {params.helper_annots} \
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
# STAGE 7: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================
# Gated to brca, nbl, and rbl only — the three profiles that have
# dsmz_counts_rds + dsmz_meta_csv configured and whose metadata schema is
# compatible with prepare_deseq2_inputs.R.  heme and pan_cancer are deferred.

DESEQ2_PROFILE_CFG = cfg.get("deseq2", {})
DESEQ2_CFG = deep_merge(config.get("defaults", {}).get("deseq2", {}), DESEQ2_PROFILE_CFG)
DESEQ2_INPUT_DIR     = os.path.join(UNSUP_REL, "deseq2_inputs")
DESEQ2_INPUT_DIR_ABS = os.path.join(UNSUP,     "deseq2_inputs")
DESEQ2_MARKER_OUTDIR_NAME = DESEQ2_CFG.get("marker_outdir_name", "deseq2_markers")
DESEQ2_ISOLATE_DIR     = os.path.join(UNSUP_REL, DESEQ2_MARKER_OUTDIR_NAME)
DESEQ2_ISOLATE_DIR_ABS = os.path.join(UNSUP,     DESEQ2_MARKER_OUTDIR_NAME)
DESEQ2_COMP_DIR     = os.path.join(UNSUP_REL, "deseq2", "component_vs_rest")
DESEQ2_COMP_DIR_ABS = os.path.join(UNSUP,     "deseq2", "component_vs_rest")

_enabled_profiles = DESEQ2_CFG.get("enabled_profiles") or ["brca", "nbl", "rbl"]
DESEQ2_ENABLED_PROFILES = {str(p) for p in _enabled_profiles}
_enabled_override = DESEQ2_PROFILE_CFG.get("enabled")
if _enabled_override is None:
    DESEQ2_ENABLED = profile_name in DESEQ2_ENABLED_PROFILES
else:
    DESEQ2_ENABLED = bool(_enabled_override)

# Cross-direction aggregated graph node stats (produced by aggregate_graph_node_stats)
NODE_STATS_TSV = os.path.join(
    P_CONS_ALL_DIR,
    f"{_PR_PREFIX}aggregated_cell_line_similarity_graph_node_stats.tsv",
)
DESEQ2_NODE_STATS_TSV = P_CONS_RESOLVED_NODE_STATS_TSV

if DESEQ2_ENABLED:

    # -------------------------------------------------------------------------
    # RULE: aggregate_graph_node_stats
    #
    # Picks the best-direction graph node annotations and writes a single
    # patient_referenced_aggregated_cell_line_similarity_graph_node_stats.tsv
    # (or pan_cancer_aggregated_… for the multicohort profile) in P_CONS_ALL_DIR.
    # Best direction is selected at runtime from p_consensus_direction_summary.tsv.
    # -------------------------------------------------------------------------
    rule aggregate_graph_node_stats:
        input:
            dir_summary_tsv = os.path.join(P_CONS_ALL_DIR, "p_consensus_direction_summary.tsv"),
            node_annots     = expand(
                os.path.join(TUMOUR_NH_ROOT, "{direction}", "final_consensus",
                             "cell_line_similarity_graph_node_annotations_{direction}.tsv"),
                direction=CONS_DIRECTIONS
            )
        output:
            node_stats_tsv = NODE_STATS_TSV
        log: os.path.join(LOGROOT, "aggregate_graph_node_stats.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript -e '
            suppressPackageStartupMessages(library(readr))
            dir_sum <- read_tsv("{input.dir_summary_tsv}", show_col_types = FALSE)
            best_dir <- dir_sum[order(-dir_sum$frac_ge_thr), ][1, "direction", drop = TRUE]
            cat("Best overall direction:", best_dir, "\n")
            annot_path <- file.path("{TUMOUR_NH_ROOT_ABS}", best_dir, "final_consensus",
                                    paste0("cell_line_similarity_graph_node_annotations_", best_dir, ".tsv"))
            annot <- read_tsv(annot_path, show_col_types = FALSE)
            out <- data.frame(
                cell_line      = annot$cell_line,
                sample_id      = if ("sample_id" %in% colnames(annot)) annot$sample_id else annot$cell_line,
                cell_line_display = if ("cell_line_display" %in% colnames(annot)) annot$cell_line_display else annot$cell_line,
                component      = annot$component,
                is_isolate     = as.logical(annot$is_outlier),
                degree         = annot$degree,
                betweenness    = annot$betweenness,
                community_louv = annot$community_louv,
                community_leid = annot$community_leid,
                stringsAsFactors = FALSE
            )
            write_tsv(out, "{output.node_stats_tsv}")
            cat("Wrote", nrow(out), "rows to {output.node_stats_tsv}\n")
            ' > {log} 2>&1
            test -s {output.node_stats_tsv} || (echo "ERROR: missing {output.node_stats_tsv}" >&2; exit 1)
            '''


    # Resolved paths for DSMZ raw counts and metadata used by prepare_deseq2_inputs.R.
    DSMZ_COUNTS_RDS = cfgget_path_abs("data/dsmz/DSMZ_count_gene.rds", "paths", "dsmz_counts_rds")
    DSMZ_META_CSV   = cfgget_path_abs(
        cfgget_path_rel("data/dsmz/DSMZ_metadata.csv", "paths", "meta_tsv"),
        "paths", "dsmz_meta_csv"
    )

    # -----------------------------------------------------------------------------
    # RULE: prepare_deseq2_inputs
    # -----------------------------------------------------------------------------
    rule prepare_deseq2_inputs:
        """
        Extracts disease-specific raw counts and metadata from DSMZ archives,
        subsetting to the cell lines present in the aggregated graph node stats.
        Produces counts.tsv and metadata.tsv consumed by all downstream DESeq2 rules.
        """
        input:
            dsmz_counts = DSMZ_COUNTS_RDS,
            dsmz_meta   = DSMZ_META_CSV,
            node_stats  = DESEQ2_NODE_STATS_TSV
        output:
            counts_tsv   = os.path.join(DESEQ2_INPUT_DIR, "counts.tsv"),
            metadata_tsv = os.path.join(DESEQ2_INPUT_DIR, "metadata.tsv")
        params:
            script = os.path.join(BASE, "scripts", "prepare_deseq2_inputs.R"),
            outdir = DESEQ2_INPUT_DIR_ABS
        log: os.path.join(LOGROOT, "prepare_deseq2_inputs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p {params.outdir}
            Rscript {params.script} \
              --profile "{profile_name}" \
              --dsmz_counts {input.dsmz_counts} \
              --dsmz_meta {input.dsmz_meta} \
              --node_stats {input.node_stats} \
              --outdir {params.outdir} \
              > {log} 2>&1
            test -s {output.counts_tsv}   || (echo "ERROR: missing {output.counts_tsv}"   >&2; exit 1)
            test -s {output.metadata_tsv} || (echo "ERROR: missing {output.metadata_tsv}" >&2; exit 1)
            '''


    # -----------------------------------------------------------------------------
    # RULE: add_component_to_metadata
    # -----------------------------------------------------------------------------
    rule add_component_to_metadata:
        """
        Joins graph-derived component and is_isolate columns onto the staged
        metadata produced by prepare_deseq2_inputs, using the aggregated node
        stats as the source of component assignments.
        """
        input:
            meta       = os.path.join(DESEQ2_INPUT_DIR, "metadata.tsv"),
            node_stats = DESEQ2_NODE_STATS_TSV
        output:
            meta_comp = os.path.join(DESEQ2_INPUT_DIR, "metadata_with_components.tsv")
        params:
            script       = os.path.join(BASE, "scripts", "add_component_to_metadata.R"),
            cell_line_col = DESEQ2_CFG.get("cell_line_col", "cell_line")
        log: os.path.join(LOGROOT, "add_component_to_metadata.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            Rscript {params.script} \
              --node_stats {input.node_stats} \
              --meta {input.meta} \
              --cell_line_col {params.cell_line_col} \
              --output {output.meta_comp} \
              > {log} 2>&1
            test -s {output.meta_comp} || (echo "ERROR: missing {output.meta_comp}" >&2; exit 1)
            '''


    # -----------------------------------------------------------------------------
    # RULE: derive_isolate_list
    # -----------------------------------------------------------------------------
    rule derive_isolate_list:
        """
        Extracts the comma-separated list of isolate cell lines (is_isolate == TRUE)
        from metadata_with_components.tsv for use as --isolate_list in
        deseq2_isolate_degs.R.  Also writes a TSV manifest so the graph-derived
        cell-line selection is preserved as a tracked, inspectable pipeline output.
        Fails clearly if required columns are absent.
        """
        input:
            meta_comp = os.path.join(DESEQ2_INPUT_DIR, "metadata_with_components.tsv")
        output:
            isolate_csv = os.path.join(DESEQ2_INPUT_DIR, "isolate_list.csv"),
            isolate_tsv = os.path.join(DESEQ2_INPUT_DIR, "isolate_cell_lines.tsv")
        log: os.path.join(LOGROOT, "derive_isolate_list.log")
        shell:
            r'''
            awk -F'\t' 'BEGIN {{ OFS="\t" }}
              NR==1 {{
                for(i=1;i<=NF;i++) {{
                  if($i=="cell_line") cl=i
                  if($i=="is_isolate") iso=i
                  if($i=="component") comp=i
                  if($i=="degree") deg=i
                  if($i=="betweenness") btw=i
                }}
                if(cl==0) {{ print "ERROR: column cell_line not found in " FILENAME > "/dev/stderr"; exit 1 }}
                if(iso==0) {{ print "ERROR: column is_isolate not found in " FILENAME > "/dev/stderr"; exit 1 }}
                print "cell_line","is_isolate","component","degree","betweenness"
              }}
              NR>1 && $iso=="TRUE" {{
                print $cl,$iso,(comp ? $comp : ""),(deg ? $deg : ""),(btw ? $btw : "")
              }}' {input.meta_comp} > {output.isolate_tsv} 2> {log}
            awk -F'\t' 'NR>1 {{ print $1 }}' {output.isolate_tsv} | sort -u | paste -sd',' > {output.isolate_csv}
            '''


    # -----------------------------------------------------------------------------
    # RULE: derive_components_list
    # -----------------------------------------------------------------------------
    rule derive_components_list:
        """
        Extracts unique non-singleton component IDs from metadata_with_components.tsv.
        Components with fewer than 2 members are excluded since component-vs-rest
        DESeq2 requires at least one sample on each side of the contrast.  Writes
        a component-size manifest for reproducibility and reporting.
        """
        input:
            meta_comp = os.path.join(DESEQ2_INPUT_DIR, "metadata_with_components.tsv")
        output:
            comp_list = os.path.join(DESEQ2_INPUT_DIR, "components_list.txt"),
            comp_tsv  = os.path.join(DESEQ2_INPUT_DIR, "component_cell_line_summary.tsv")
        log: os.path.join(LOGROOT, "derive_components_list.log")
        shell:
            r'''
            awk -F'\t' 'BEGIN {{ OFS="\t" }}
              NR==1 {{
                for(i=1;i<=NF;i++) if($i=="component") {{ c=i; break }}
                if(!c) {{
                  print "ERROR: component header missing in " FILENAME > "/dev/stderr"
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
                if(val !~ /^[0-9]+$/) {{
                  print "ERROR: non-integer component value on row " NR ": \"" val "\"" > "/dev/stderr"
                  exit 1
                }}
                counts[val+0]++
              }}
              END {{
                print "component","n_cell_lines","eligible_for_deseq"
                for(k in counts) print k+0,counts[k],(counts[k]>=2 ? "TRUE" : "FALSE")
              }}' {input.meta_comp} > {output.comp_tsv}.tmp 2> {log}
            head -n 1 {output.comp_tsv}.tmp > {output.comp_tsv}
            tail -n +2 {output.comp_tsv}.tmp | sort -t '	' -k1,1n >> {output.comp_tsv}
            rm -f {output.comp_tsv}.tmp
            awk -F'\t' 'NR>1 && $3=="TRUE" {{ print $1 }}' {output.comp_tsv} | sort -n > {output.comp_list}
            '''


    # -----------------------------------------------------------------------------
    # RULE: derive_anchor_list
    # -----------------------------------------------------------------------------
    rule derive_anchor_list:
        """
        Extracts canonical anchor cell lines from
        patient_referenced_resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv
        (or the pan_cancer_… variant for the multicohort profile; produced by
        plot_patient_referenced_resolved_cell_line_neighbourhood_graph).  The anchor set is
        the union of:
          - degree_anchor_selected     : highest-degree node per component
            (alias for most_connected_selected)
          - bridge_betweenness_selected : maximum unweighted unnormalised
            betweenness node per component
            (alias for canonical_bridge_selected)
        Selection logic:
          Primary  -- use anchor_selected column if present (schema v2+)
          Fallback -- union of most_connected_selected OR canonical_bridge_selected
            (schema v1 backward compatibility)
        Duplicate nodes (selected by both criteria) are de-duplicated; both
        reasons are retained in anchor_components.tsv via the source audit file.
        Two output files are written for use in deseq2_isolate_degs.R:
          anchor_list.csv        -- comma-separated node IDs of canonical anchors
          anchor_components.tsv  -- two-column TSV (anchor<TAB>component)
        """
        input:
            anchor_audit = P_CONS_ANCHOR_AUDIT_TSV
        output:
            anchor_list_csv       = os.path.join(DESEQ2_INPUT_DIR, "anchor_list.csv"),
            anchor_components_tsv = os.path.join(DESEQ2_INPUT_DIR, "anchor_components.tsv")
        log: os.path.join(LOGROOT, "derive_anchor_list.log")
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
                  selected = ($asel=="True" || $asel=="TRUE")
                }} else {{
                  selected = ($mcs=="True" || $mcs=="TRUE" || $cbs=="True" || $cbs=="TRUE")
                }}
                if(selected) print $nid,$cid
              }}
            ' {input.anchor_audit} > {output.anchor_components_tsv} 2> {log}
            awk -F'\t' 'NR>1 {{ print $1 }}' {output.anchor_components_tsv} \
              | sort -u | paste -sd',' > {output.anchor_list_csv}
            echo "[derive_anchor_list] Schema: $(head -1 {input.anchor_audit} | tr '\t' '\n' | grep -n 'anchor_selected\|most_connected\|canonical_bridge' | head -5)" >> {log}
            echo "[derive_anchor_list] Anchors: $(cat {output.anchor_list_csv})" >> {log}
            echo "[derive_anchor_list] anchor_components rows: $(wc -l < {output.anchor_components_tsv})" >> {log}
            '''


    # -----------------------------------------------------------------------------
    # RULE: deseq2_isolate_degs
    # -----------------------------------------------------------------------------
    rule deseq2_isolate_degs:
        """
        Runs DESeq2 contrasts for each graph-derived group identified in the
        cell-line similarity graph:
          - isolate-vs-REST: each degree-zero cell line vs all others
          - anchor-vs-outside-component: each canonical anchor vs cell lines
            outside its component (union of most-connected and bridge anchors,
            de-duplicated by the derive_anchor_list rule)
        Produces per-contrast DEG tables, filtered marker gene lists, size-factor
        QC, and a recurrence-based unique feature set.
        """
        input:
            counts_tsv            = os.path.join(DESEQ2_INPUT_DIR, "counts.tsv"),
            meta_comp             = os.path.join(DESEQ2_INPUT_DIR, "metadata_with_components.tsv"),
            isolate_csv           = os.path.join(DESEQ2_INPUT_DIR, "isolate_list.csv"),
            isolate_tsv           = os.path.join(DESEQ2_INPUT_DIR, "isolate_cell_lines.tsv"),
            anchor_list_csv       = os.path.join(DESEQ2_INPUT_DIR, "anchor_list.csv"),
            anchor_components_tsv = os.path.join(DESEQ2_INPUT_DIR, "anchor_components.tsv")
        output:
            size_factors = os.path.join(DESEQ2_ISOLATE_DIR, "qc", "size_factors.tsv"),
            manifest     = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "marker_sets_manifest.tsv"),
            recurrence   = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "gene_recurrence_across_contrasts.tsv"),
            unique_set   = os.path.join(DESEQ2_ISOLATE_DIR, "markers",
                               f"unique_feature_set_recurrence_ge_{DESEQ2_CFG.get('recurrence_k', 2)}.txt"),
            session_info = os.path.join(DESEQ2_ISOLATE_DIR, "sessionInfo.txt")
        params:
            script        = os.path.join(BASE, "scripts", "deseq2_isolate_degs.R"),
            outdir        = DESEQ2_ISOLATE_DIR_ABS,
            sample_id_col = DESEQ2_CFG.get("staged_sample_id_col", "sample_id"),
            cell_line_col = DESEQ2_CFG.get("cell_line_col", "cell_line"),
            component_col = DESEQ2_CFG.get("component_col", "component"),
            fdr           = DESEQ2_CFG.get("fdr_isolate", 0.01),
            lfc           = DESEQ2_CFG.get("lfc_isolate", 1.5),
            topN          = DESEQ2_CFG.get("topN_isolate", 50),
            fdr_anchor    = DESEQ2_CFG.get("fdr_anchor",  0.05),
            lfc_anchor    = DESEQ2_CFG.get("lfc_anchor",  1.0),
            topN_anchor   = DESEQ2_CFG.get("topN_anchor", 200),
            min_baseMean  = DESEQ2_CFG.get("min_baseMean", 10),
            recurrence_k  = DESEQ2_CFG.get("recurrence_k", 2)
        log: os.path.join(LOGROOT, "deseq2_isolate_degs.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            ISOLATE_LIST=$(cat {input.isolate_csv})
            ANCHOR_LIST=$(cat {input.anchor_list_csv})
            if [ -z "$ISOLATE_LIST" ] && [ -z "$ANCHOR_LIST" ]; then
                echo "[WARN] No isolates and no anchors found; creating empty outputs" > {log}
                mkdir -p {params.outdir}/qc {params.outdir}/markers
                echo -e "sample_id\tsize_factor" > {output.size_factors}
                echo -e "contrast\tmarker_file\ttable_file\tn_markers" > {output.manifest}
                echo -e "gene_id\tn_contrasts" > {output.recurrence}
                touch {output.unique_set}
                touch {output.session_info}
                exit 0
            fi
            mkdir -p {params.outdir}
            ANCHOR_ARGS=""
            if [ -n "$ANCHOR_LIST" ]; then
                ANCHOR_ARGS="--anchor_list $ANCHOR_LIST --anchor_components {input.anchor_components_tsv} --fdr_anchor {params.fdr_anchor} --lfc_anchor {params.lfc_anchor} --topN_anchor {params.topN_anchor}"
            fi
            Rscript {params.script} \
              --counts {input.counts_tsv} \
              --meta {input.meta_comp} \
              --sample_id_col {params.sample_id_col} \
              --cell_line_col {params.cell_line_col} \
              --component_col {params.component_col} \
              --isolate_list "$ISOLATE_LIST" \
              --outdir {params.outdir} \
              --fdr_isolate {params.fdr} \
              --lfc_isolate {params.lfc} \
              --topN_isolate {params.topN} \
              --min_baseMean {params.min_baseMean} \
              --recurrence_k {params.recurrence_k} \
              $ANCHOR_ARGS \
              > {log} 2>&1
            '''

    DESEQ2_DIRECTIONAL_MARKER_DIR = os.path.join(DESEQ2_ISOLATE_DIR, "directional_markers")
    DESEQ2_DIRECTIONAL_MARKER_DIR_ABS = os.path.join(DESEQ2_ISOLATE_DIR_ABS, "directional_markers")

    # -----------------------------------------------------------------------------
    # RULE: write_deseq2_directional_marker_tables
    # -----------------------------------------------------------------------------
    rule write_deseq2_directional_marker_tables:
        """
        Organise retained DESeq2 marker genes into directional marker tables
        based on the sign of log2FoldChange. This is a downstream-only step:
        it consumes the completed marker manifest, marker gene lists, and full
        DESeq2 tables without rerunning DESeq2.
        """
        input:
            manifest = ancient(os.path.join(DESEQ2_ISOLATE_DIR, "markers", "marker_sets_manifest.tsv")),
            markers_dir = ancient(os.path.join(DESEQ2_ISOLATE_DIR, "markers")),
            tables_dir = ancient(os.path.join(DESEQ2_ISOLATE_DIR, "tables")),
            deseq2_script = os.path.join(BASE, "scripts", "deseq2_isolate_degs.R")
        output:
            orientation = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "directional_marker_orientation.txt"),
            all_up = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "all_upregulated_markers.tsv"),
            all_down = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "all_downregulated_markers.tsv"),
            counts = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "directional_marker_counts.tsv"),
            top_up = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "top_upregulated_markers.tsv"),
            top_down = os.path.join(DESEQ2_DIRECTIONAL_MARKER_DIR, "cohort_level", "top_downregulated_markers.tsv")
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

    rule write_all_deseq2_directional_marker_tables:
        """
        Aggregate target for the active profile's directional marker files.
        """
        input:
            rules.write_deseq2_directional_marker_tables.output

    # -----------------------------------------------------------------------------
    # RULE: deseq2_component_vs_rest_all
    # -----------------------------------------------------------------------------
    rule deseq2_component_vs_rest_all:
        """
        Runs deseq2_component_vs_rest.R once per component in components_list.txt.
        A shell loop iterates over component IDs, running the R script for each.
        A sentinel .done file is written only after all components complete.
        """
        input:
            counts_tsv = os.path.join(DESEQ2_INPUT_DIR, "counts.tsv"),
            meta_comp  = os.path.join(DESEQ2_INPUT_DIR, "metadata_with_components.tsv"),
            comp_list  = os.path.join(DESEQ2_INPUT_DIR, "components_list.txt"),
            comp_tsv   = os.path.join(DESEQ2_INPUT_DIR, "component_cell_line_summary.tsv"),
            isolate_manifest = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "marker_sets_manifest.tsv")
        output:
            done = touch(os.path.join(DESEQ2_COMP_DIR, ".done"))
        params:
            script        = os.path.join(BASE, "scripts", "deseq2_component_vs_rest.R"),
            outdir        = DESEQ2_COMP_DIR_ABS,
            component_col = DESEQ2_CFG.get("component_col", "component"),
            sample_id_col = DESEQ2_CFG.get("staged_sample_id_col", "sample_id"),
            fdr           = DESEQ2_CFG.get("fdr_component", 0.05),
            lfc           = DESEQ2_CFG.get("lfc_component", 1.0),
            topN          = DESEQ2_CFG.get("topN_component", 500)
        log: os.path.join(LOGROOT, "deseq2_component_vs_rest_all.log")
        conda: CONDA_ENV_R
        shell:
            r'''
            mkdir -p {params.outdir}
            echo "[START] Component-vs-rest loop" > {log}
            if [ ! -s {input.comp_list} ]; then
                echo "[WARN] No components to process; sentinel written" >> {log}
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
                  --outdir {params.outdir} \
                  --fdr {params.fdr} \
                  --lfc {params.lfc} \
                  --topN {params.topN} \
                  >> {log} 2>&1
            done < {input.comp_list}
            echo "[DONE] All components processed" >> {log}
            '''


# =============================================================================
# STAGE 7B: GENE SET ENRICHMENT ANALYSIS
# =============================================================================

ENRICH_PROFILE_CFG = cfg.get("enrichment", {})
ENRICH_CFG = deep_merge(config.get("defaults", {}).get("enrichment", {}), ENRICH_PROFILE_CFG)

ENRICH_ENABLED = bool(ENRICH_CFG.get("enabled", False)) and DESEQ2_ENABLED

ENRICH_DIR_REL = os.path.join(UNSUP_REL, ENRICH_CFG.get("output_dir_name", "enrichment"))
ENRICH_DIR_ABS = os.path.join(UNSUP, ENRICH_CFG.get("output_dir_name", "enrichment"))

_ENRICH_DISEASE_PROFILES = ENRICH_CFG.get(
    "disease_profiles",
    config.get("defaults", {}).get("marker_postprocessing", {}).get("pan_cancer", {}).get(
        "disease_profiles",
        sorted(DESEQ2_ENABLED_PROFILES)
    )
)
_ENRICH_DISEASE_PROFILES = [str(p) for p in _ENRICH_DISEASE_PROFILES]

if ENRICH_ENABLED:

    rule build_enrichment_query_sets:
        """
        Build enrichment query gene sets and matched custom backgrounds from
        DESeq2 isolate and component outputs. The script discovers contrasts,
        components, and disease-level recurrence sets from manifests and files.
        """
        input:
            isolate_manifest = os.path.join(DESEQ2_ISOLATE_DIR, "markers", "marker_sets_manifest.tsv"),
            isolate_recurrence = os.path.join(
                DESEQ2_ISOLATE_DIR,
                "markers",
                f"unique_feature_set_recurrence_ge_{DESEQ2_CFG.get('recurrence_k', 2)}.txt"
            ),
            component_done = os.path.join(DESEQ2_COMP_DIR, ".done")
        output:
            query_manifest = os.path.join(ENRICH_DIR_REL, "query_sets", "query_manifest.tsv"),
            skipped_queries = os.path.join(ENRICH_DIR_REL, "query_sets", "skipped_queries.tsv"),
            done = touch(os.path.join(ENRICH_DIR_REL, "query_sets", ".done"))
        params:
            script = os.path.join(SCRIPTS_DIR, "build_enrichment_query_sets.R"),
            isolate_dir = DESEQ2_ISOLATE_DIR_ABS,
            component_dir = DESEQ2_COMP_DIR_ABS,
            outdir = os.path.join(ENRICH_DIR_ABS, "query_sets"),
            profile = profile_name,
            min_query_genes = ENRICH_CFG.get("min_query_genes", 5),
            recurrence_k = DESEQ2_CFG.get("recurrence_k", 2),
            component_marker_pooling = ENRICH_CFG.get("component_marker_pooling", "union"),
            component_recurrence_min = ENRICH_CFG.get("component_recurrence_min", 1),
            component_rank_stat = ENRICH_CFG.get("component_rank_stat", "median_abs_log2fc"),
            operative_feature_set_gene_file = abspath(ENRICH_CFG.get("operative_feature_set_gene_file", "")) if ENRICH_CFG.get("operative_feature_set_gene_file", "") else "",
            operative_feature_set_rank_tsv = abspath(ENRICH_CFG.get("operative_feature_set_rank_tsv", "")) if ENRICH_CFG.get("operative_feature_set_rank_tsv", "") else "",
            strict_union_primary = "TRUE" if ENRICH_CFG.get("strict_union_primary", True) else "FALSE",
            operative_feature_set_primary = "TRUE" if ENRICH_CFG.get("operative_feature_set_primary", True) else "FALSE",
            compare_strict_vs_operative = "TRUE" if ENRICH_CFG.get("compare_strict_vs_operative", True) else "FALSE",
            ordered_strict_union = "TRUE" if ENRICH_CFG.get("ordered_strict_union", True) else "FALSE",
            ordered_operative_feature_set = "TRUE" if ENRICH_CFG.get("ordered_operative_feature_set", False) else "FALSE",
            disease_profiles = ",".join(_ENRICH_DISEASE_PROFILES)
        log:
            os.path.join(LOGROOT, "build_enrichment_query_sets.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            mkdir -p "{params.outdir}"

            Rscript "{params.script}" \
              --profile "{params.profile}" \
              --isolate-dir "{params.isolate_dir}" \
              --component-dir "{params.component_dir}" \
              --outdir "{params.outdir}" \
              --min-query-genes {params.min_query_genes} \
              --recurrence-k {params.recurrence_k} \
              --component-marker-pooling "{params.component_marker_pooling}" \
              --component-recurrence-min {params.component_recurrence_min} \
              --component-rank-stat "{params.component_rank_stat}" \
              --operative-feature-set-gene-file "{params.operative_feature_set_gene_file}" \
              --operative-feature-set-rank-tsv "{params.operative_feature_set_rank_tsv}" \
              --strict-union-primary {params.strict_union_primary} \
              --operative-feature-set-primary {params.operative_feature_set_primary} \
              --compare-strict-vs-operative {params.compare_strict_vs_operative} \
              --ordered-strict-union {params.ordered_strict_union} \
              --ordered-operative-feature-set {params.ordered_operative_feature_set} \
              --disease-profiles "{params.disease_profiles}" \
              > "{log}" 2>&1

            test -s "{output.query_manifest}" || (echo "ERROR: missing {output.query_manifest}" >&2; exit 1)
            '''

    rule run_gprofiler_enrichment:
        """
        Run g:Profiler enrichment for all non-skipped query sets listed in the
        query manifest. Produces primary IEA-included results and optional IEA
        excluded sensitivity summaries.
        """
        input:
            query_manifest = os.path.join(ENRICH_DIR_REL, "query_sets", "query_manifest.tsv"),
            query_sets_done = os.path.join(ENRICH_DIR_REL, "query_sets", ".done")
        output:
            corpus_manifest = os.path.join(ENRICH_DIR_REL, "gprofiler", "corpus_manifest.tsv"),
            top_terms = os.path.join(ENRICH_DIR_REL, "gprofiler", "top_terms.tsv"),
            iea_sensitivity_summary = os.path.join(ENRICH_DIR_REL, "gprofiler", "iea_sensitivity_summary.tsv"),
            version_info = os.path.join(ENRICH_DIR_REL, "gprofiler", "gprofiler_version.tsv"),
            done = touch(os.path.join(ENRICH_DIR_REL, "gprofiler", ".done"))
        params:
            script = os.path.join(SCRIPTS_DIR, "run_gprofiler_from_manifest.R"),
            outdir = os.path.join(ENRICH_DIR_ABS, "gprofiler"),
            organism = ENRICH_CFG.get("organism", "hsapiens"),
            sources = ",".join(ENRICH_CFG.get(
                "sources",
                ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", "TF", "HPA", "CORUM"]
            )),
            alpha = ENRICH_CFG.get("alpha", 0.05),
            correction_method = ENRICH_CFG.get("correction_method", "g_SCS"),
            archive_url = ENRICH_CFG.get("archive_url", ""),
            require_archive = "TRUE" if ENRICH_CFG.get("require_archive", True) else "FALSE",
            run_iea_sensitivity = "TRUE" if ENRICH_CFG.get("run_iea_sensitivity", True) else "FALSE",
            as_short_link = "TRUE" if ENRICH_CFG.get("as_short_link", False) else "FALSE",
            top_terms_per_source = ENRICH_CFG.get("top_terms_per_source", 5)
        log:
            os.path.join(LOGROOT, "run_gprofiler_enrichment.log")
        conda:
            CONDA_ENV_R
        shell:
            r'''
            set -euo pipefail
            mkdir -p "{params.outdir}"

            Rscript "{params.script}" \
              --query-manifest "{input.query_manifest}" \
              --outdir "{params.outdir}" \
              --organism "{params.organism}" \
              --sources "{params.sources}" \
              --alpha {params.alpha} \
              --correction-method "{params.correction_method}" \
              --archive-url "{params.archive_url}" \
              --require-archive {params.require_archive} \
              --run-iea-sensitivity {params.run_iea_sensitivity} \
              --as-short-link {params.as_short_link} \
              --top-terms-per-source {params.top_terms_per_source} \
              > "{log}" 2>&1

            test -s "{output.corpus_manifest}" || (echo "ERROR: missing {output.corpus_manifest}" >&2; exit 1)
            test -s "{output.top_terms}" || (echo "ERROR: missing {output.top_terms}" >&2; exit 1)
            '''


    rule build_enrichment_summary_top_terms:
        input:
            top_terms = os.path.join(ENRICH_DIR_REL, "gprofiler", "top_terms.tsv"),
            query_manifest = os.path.join(ENRICH_DIR_REL, "query_sets", "query_manifest.tsv")
        output:
            summary = os.path.join(ENRICH_DIR_REL, "enrichment_summary_top_terms.tsv"),
            provenance = os.path.join(ENRICH_DIR_REL, "enrichment_summary_top_terms_provenance.tsv")
        params:
            script = os.path.join(SCRIPTS_DIR, "build_enrichment_summary_top_terms.R")
        log: os.path.join(LOGROOT, "build_enrichment_summary_top_terms.log")
        conda: CONDA_ENV_R
        shell:
            "Rscript {params.script} --input {input.top_terms} --query-manifest {input.query_manifest} --output {output.summary} > {log} 2>&1 && test -s {output.summary} && test -s {output.provenance}"

    rule plot_enrichment_top_terms_heatmap:
        input:
            summary = os.path.join(ENRICH_DIR_REL, "enrichment_summary_top_terms.tsv")
        output:
            pdf = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap.pdf"),
            png = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap.png"),
            matrix = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap_matrix.tsv"),
            selected = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap_selected_terms.tsv"),
            excluded = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv"),
            provenance = os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap_provenance.tsv")
        params:
            script = os.path.join(SCRIPTS_DIR, "plot_enrichment_top_terms_heatmap.R"),
            outdir = os.path.join(ENRICH_DIR_ABS, "figures")
        log: os.path.join(LOGROOT, "plot_enrichment_top_terms_heatmap.log")
        conda: CONDA_ENV_R
        shell:
            "mkdir -p {params.outdir} && Rscript {params.script} --input {input.summary} --outdir {params.outdir} > {log} 2>&1 && test -s {output.pdf} && test -s {output.provenance}"


# =============================================================================
# STAGE 8: MARKER POST-PROCESSING (CONSENSUS → FEATURES → EXPRESSION → MAPPING)
# =============================================================================

MARKER_POST_CFG = config.get("defaults", {}).get("marker_postprocessing", {})
MARKER_POST_ENABLED = bool(MARKER_POST_CFG.get("enabled", False)) and profile_name == "multicohort_cancer"

if MARKER_POST_ENABLED:
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

    def profile_consensus_tables_dir(profile):
        return os.path.join("results", "unsupervised", profile, "deseq2_markers", "tables")

    def profile_consensus_outdir(profile):
        return os.path.join(profile_unsup_root_abs(profile), "deseq2_markers", CONSENSUS_OUTDIR_NAME)

    def profile_consensus_summary_dir(profile):
        return os.path.join(profile_consensus_outdir(profile), "summary")

    def profile_consensus_summary_rel(profile):
        return os.path.join(
            "results", "unsupervised", profile,
            "deseq2_markers", CONSENSUS_OUTDIR_NAME, "summary"
        )

    def profile_marker_manifest_rel(profile):
        return os.path.join(
            "results", "unsupervised", profile,
            "deseq2_markers", "markers", "marker_sets_manifest.tsv"
        )

    def profile_marker_tables_rel(profile):
        return os.path.join("results", "unsupervised", profile, "deseq2_markers", "tables")

    rule ensure_profile_deseq2_markers:
        """
        Build graph-derived isolate-vs-rest DESeq2 marker outputs for a disease
        profile by invoking the same Snakefile with pipeline_profile set to that
        disease. This lets the multicohort pan-cancer target depend on tracked
        disease-level DESeq2 outputs without hardcoded cell-line or gene lists.
        """
        output:
            tables_dir = directory(os.path.join(
                "results", "unsupervised", "{profile}", "deseq2_markers", "tables"
            )),
            marker_manifest = os.path.join(
                "results", "unsupervised", "{profile}", "deseq2_markers",
                "markers", "marker_sets_manifest.tsv"
            ),
            session_info = os.path.join(
                "results", "unsupervised", "{profile}", "deseq2_markers",
                "sessionInfo.txt"
            )
        params:
            snakemake_bin = os.path.join(BASE, "envs", ".conda", "smk", "bin", "snakemake"),
            snakefile = os.path.join(BASE, "Snakefile"),
            configfile = os.path.join(BASE, "config", "config.yaml")
        log:
            os.path.join(LOGROOT, "ensure_profile_deseq2_markers_{profile}.log")
        shell:
            r'''
            "{params.snakemake_bin}" \
              --snakefile "{params.snakefile}" \
              --configfile "{params.configfile}" \
              --config pipeline_profile="{wildcards.profile}" \
              --use-conda \
              --nolock \
              --cores 1 \
              deseq2_isolate_degs \
              > "{log}" 2>&1
            test -d "{output.tables_dir}" || (echo "ERROR: missing {output.tables_dir}" >&2; exit 1)
            test -s "{output.marker_manifest}" || (echo "ERROR: missing {output.marker_manifest}" >&2; exit 1)
            test -s "{output.session_info}" || (echo "ERROR: missing {output.session_info}" >&2; exit 1)
            '''

    rule build_component_consensus_markers:
        """Build per-profile consensus markers and summary exports."""
        input:
            tables_dir=lambda wc: profile_marker_tables_rel(wc.profile),
            marker_manifest=lambda wc: profile_marker_manifest_rel(wc.profile)
        output:
            summary_done=os.path.join(
                "results", "unsupervised", "{profile}",
                "deseq2_markers", CONSENSUS_OUTDIR_NAME,
                "summary", "consensus_export_done.txt"
            )
        params:
            script=os.path.join(SCRIPTS_DIR, "build_component_consensus_markers_v2.py"),
            outdir=lambda wc: profile_consensus_outdir(wc.profile),
            padj=CONSENSUS_CFG.get("padj", 0.05),
            basemean=CONSENSUS_CFG.get("basemean", 1.0),
            lfc=CONSENSUS_CFG.get("lfc", 0.5),
            core_support_min=CONSENSUS_CFG.get("core_support_min", 2),
            core_support_frac=CONSENSUS_CFG.get("core_support_frac", 0.30),
            min_anchor_genes=CONSENSUS_CFG.get("min_anchor_genes", 25),
            fallback_topn=CONSENSUS_CFG.get("fallback_topn", 200),
            flip_policy=CONSENSUS_CFG.get("flip_policy", "remove"),
            recurrence_k=CONSENSUS_CFG.get("recurrence_k", 2)
        log:
            os.path.join(LOGROOT, "build_component_consensus_markers_{profile}.log")
        conda: CONDA_ENV_PY
        shell:
            r'''
            mkdir -p "{params.outdir}"
            python "{params.script}" \
              --tables-dir "{input.tables_dir}" \
              --profile-name "{wildcards.profile}" \
              --outdir "{params.outdir}" \
              --padj {params.padj} \
              --basemean {params.basemean} \
              --lfc {params.lfc} \
              --core-support-min {params.core_support_min} \
              --core-support-frac {params.core_support_frac} \
              --min-anchor-genes {params.min_anchor_genes} \
              --fallback-topn {params.fallback_topn} \
              --flip-policy {params.flip_policy} \
              --recurrence-k {params.recurrence_k} \
              > "{log}" 2>&1
            test -s "{output.summary_done}" || (echo "ERROR: missing {output.summary_done}" >&2; exit 1)
            '''

    def consensus_summary_done(profile):
        return os.path.join(profile_consensus_summary_dir(profile), "consensus_export_done.txt")

    def profile_summary_dir(profile):
        return profile_consensus_summary_dir(profile)

    PAN_FEATURES_OUTDIR = PAN_CANCER_MP_CFG.get("outdir", "results/unsupervised/pan_cancer/feature_space")
    PAN_FEATURES_OUTDIR_ABS = abspath(PAN_FEATURES_OUTDIR)
    PAN_FEATURES_TSV = abspath(PAN_CANCER_MP_CFG.get("final_features_tsv", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.tsv")))
    PAN_FEATURES_CLEAN = abspath(PAN_CANCER_MP_CFG.get("final_features_clean", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features_clean.txt")))
    PAN_FEATURES_UP = abspath(PAN_CANCER_MP_CFG.get("final_features_up", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.UP.txt")))
    PAN_FEATURES_DOWN = abspath(PAN_CANCER_MP_CFG.get("final_features_down", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_features.DOWN.txt")))
    PAN_FEATURES_SUMMARY = abspath(PAN_CANCER_MP_CFG.get("build_summary_tsv", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_feature_build_summary.tsv")))
    PAN_FEATURES_REPORT = abspath(PAN_CANCER_MP_CFG.get("build_report_md", os.path.join(PAN_FEATURES_OUTDIR, "pan_cancer_feature_build_report.md")))

    def _profile_dir_args():
        args = []
        for prof in PAN_PROFILES:
            args.append(f"--profile-dir {prof}={shlex.quote(profile_summary_dir(prof))}")
        return " ".join(args)

    def _profile_marker_dir_args():
        args = []
        for prof in PAN_PROFILES:
            marker_dir = os.path.join(profile_unsup_root_abs(prof), "deseq2_markers", "markers")
            args.append(f"--profile-marker-dir {prof}={shlex.quote(marker_dir)}")
        return " ".join(args)

    def profile_marker_manifest_abs(profile):
        return os.path.join(
            profile_unsup_root_abs(profile),
            "deseq2_markers", "markers", "marker_sets_manifest.tsv"
        )

    def _profile_marker_manifest_args():
        args = []
        for prof in PAN_PROFILES:
            args.append(
                f"--profile-marker-manifest {prof}="
                f"{shlex.quote(profile_marker_manifest_abs(prof))}"
            )
        return " ".join(args)

    rule build_pan_cancer_features:
        """Merge per-profile summaries into a pan-cancer feature panel."""
        input:
            marker_manifests=[
                ancient(profile_marker_manifest_rel(p))
                for p in PAN_PROFILES
            ]
        output:
            features_tsv = PAN_FEATURES_TSV,
            clean_txt   = PAN_FEATURES_CLEAN,
            up_txt      = PAN_FEATURES_UP,
            down_txt    = PAN_FEATURES_DOWN,
            summary_tsv = PAN_FEATURES_SUMMARY,
            report_md   = PAN_FEATURES_REPORT,
            done_file   = os.path.join(PAN_FEATURES_OUTDIR_ABS, "pan_cancer_features_done.txt")
        params:
            script=os.path.join(SCRIPTS_DIR, "build_pan_cancer_features.py"),
            outdir=PAN_FEATURES_OUTDIR_ABS,
            profile_marker_dirs=_profile_marker_dir_args(),
            profile_marker_manifests=_profile_marker_manifest_args(),
            cap_isolate=PAN_CANCER_MP_CFG.get("cap_isolate", 0),
            remove_ribo_mt="--remove-ribo-mt" if PAN_CANCER_MP_CFG.get("remove_ribo_mt", False) else "",
            gene_annot=("--gene-annotation-tsv " + abspath(PAN_CANCER_MP_CFG["gene_annotation_tsv"]))
                        if PAN_CANCER_MP_CFG.get("gene_annotation_tsv") else ""
        log: os.path.join(LOGROOT, "build_pan_cancer_features.log")
        conda: CONDA_ENV_PY
        shell:
            r'''
            mkdir -p "{params.outdir}"
            python "{params.script}" \
              {params.profile_marker_dirs} \
              {params.profile_marker_manifests} \
              --output-dir "{params.outdir}" \
              --cap-isolate {params.cap_isolate} \
              {params.remove_ribo_mt} \
              {params.gene_annot} \
              > "{log}" 2>&1
            test -s "{output.features_tsv}" || (echo "ERROR: missing {output.features_tsv}" >&2; exit 1)
            '''

    PAN_EXPR_RDS = abspath(PAN_EXPR_CFG.get("output_rds", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds"))
    PAN_EXPR_META = abspath(PAN_EXPR_CFG.get("output_meta_tsv", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_metadata.tsv"))
    PAN_EXPR_CELL_LINES_RDS = abspath(PAN_EXPR_CFG.get(
        "output_cell_lines_rds",
        "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only.rds"
    ))
    PAN_EXPR_CELL_LINES_META = abspath(PAN_EXPR_CFG.get(
        "output_cell_lines_meta_tsv",
        "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr_cell_lines_only_metadata.tsv"
    ))

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

    # Main expression inputs are restricted to cohorts with tumour samples.
    _PAN_EXPR_INPUTS = [_PAN_EXPR_BRCA_VST, _PAN_EXPR_NBL_VST, _PAN_EXPR_RBL_VST]
    if _PAN_EXPR_HEME_PATH:
        _PAN_EXPR_INPUTS.append(_PAN_EXPR_HEME_PATH)

    rule build_pan_cancer_feature_expression_matrix:
        """Combine per-profile VST inputs into the feature-restricted pan-cancer matrix."""
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

    rule build_pan_cancer_feature_expression_matrix_cell_lines_only:
        """Build a feature-restricted pan-cancer matrix containing cell lines only."""
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
    PAN_ALIGNMENT_UMAP_OUTDIR = abspath(PAN_ALIGNMENT_UMAP_CFG.get(
        "output_dir",
        "results/unsupervised/pan_cancer/tumour_cell_line_alignment_umap"
    ))
    PAN_ALIGNMENT_UMAP_META = abspath(PAN_ALIGNMENT_UMAP_CFG.get(
        "metadata_tsv",
        PAN_EXPR_META
    ))
    _PAN_ALIGNMENT_UMAP_METRICS_RAW = PAN_ALIGNMENT_UMAP_CFG.get("metrics", ["cosine", "euclidean"])
    if isinstance(_PAN_ALIGNMENT_UMAP_METRICS_RAW, str):
        PAN_ALIGNMENT_UMAP_METRICS = [
            m.strip() for m in _PAN_ALIGNMENT_UMAP_METRICS_RAW.split(",") if m.strip()
        ]
    else:
        PAN_ALIGNMENT_UMAP_METRICS = [str(m).strip() for m in _PAN_ALIGNMENT_UMAP_METRICS_RAW if str(m).strip()]

    def pan_alignment_umap_tag(metric):
        return "DEG_SET_" + re.sub(r"[^A-Za-z0-9]+", "_", metric)

    def pan_alignment_umap_stem(metric):
        return "pan_cancer_tumour_cell_line_alignment_umap_" + pan_alignment_umap_tag(metric)

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
    PAN_ALIGNMENT_UMAP_SUMMARY = os.path.join(
        PAN_ALIGNMENT_UMAP_OUTDIR,
        "summary_pan_cancer_tumour_cell_line_alignment_umap.tsv"
    )

    rule plot_pan_cancer_tumour_cell_line_alignment_umap:
        """Generate thesis-facing pan-cancer tumour/cell-line alignment UMAPs."""
        input:
            expr_rds = PAN_EXPR_RDS,
            meta_tsv = PAN_ALIGNMENT_UMAP_META,
            genes = PAN_FEATURES_CLEAN,
            script = os.path.join(SCRIPTS_DIR, "plot_pan_cancer_tumour_cell_line_alignment_umap.R")
        output:
            pdf = PAN_ALIGNMENT_UMAP_PDFS,
            svg = PAN_ALIGNMENT_UMAP_SVGS,
            png = PAN_ALIGNMENT_UMAP_PNGS,
            coords = PAN_ALIGNMENT_UMAP_COORDS,
            summary = PAN_ALIGNMENT_UMAP_SUMMARY
        params:
            outdir = PAN_ALIGNMENT_UMAP_OUTDIR,
            metrics = ",".join(PAN_ALIGNMENT_UMAP_METRICS),
            page = PAN_ALIGNMENT_UMAP_CFG.get("page", "thesis")
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
              --deg_set "{input.genes}" \
              --outdir "{params.outdir}" \
              --dist_metrics "{params.metrics}" \
              --page "{params.page}" \
              > "{log}" 2>&1
            for f in {output.pdf} {output.svg} {output.png} {output.coords} "{output.summary}"; do
              test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
            done
            '''

    # ==========================================================================
    # DIAGNOSTIC: VST all-gene tumour/cell-line UMAP (BRCA + NBL + RBL)
    # ==========================================================================
    # Independent of the curated pan-cancer DEG / feature-set workflow above.
    # This unsupervised diagnostic runs UMAP on the full VST gene space (after
    # NA + zero-variance filtering on the merged matrix), with no class label
    # used during fitting and no pan-cancer feature subsetting. Filenames use
    # the VST_ALL_GENE label to keep this run's outputs cleanly separable from
    # the thesis-facing DEG-set outputs above.
    VST_ALLGENE_UMAP_OUTDIR = abspath(
        "results/unsupervised/pan_cancer/vst_all_gene_tumour_cell_line_umap"
    )
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

    def _vst_allgene_basename(prefix, metric, ext):
        return f"{prefix}_{VST_ALLGENE_UMAP_STEM}_{VST_ALLGENE_UMAP_LABEL}_{metric}.{ext}"

    VST_ALLGENE_UMAP_PDFS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, _vst_allgene_basename("Fig", m, "pdf"))
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SVGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, _vst_allgene_basename("Fig", m, "svg"))
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_PNGS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, _vst_allgene_basename("Fig", m, "png"))
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_COORDS = [
        os.path.join(VST_ALLGENE_UMAP_OUTDIR, _vst_allgene_basename("coords", m, "tsv"))
        for m in VST_ALLGENE_UMAP_METRICS
    ]
    VST_ALLGENE_UMAP_SUMMARY = os.path.join(
        VST_ALLGENE_UMAP_OUTDIR,
        f"summary_{VST_ALLGENE_UMAP_STEM}_{VST_ALLGENE_UMAP_LABEL}.tsv"
    )

    # Direct, profile-independent paths to the three per-cohort joint VST RDS
    # files (samples x genes after orientation handling in the merge script).
    VST_ALLGENE_BRCA_VST = os.path.join(BASE, "data", "brca", "brca_vst_joint.rds")
    VST_ALLGENE_NBL_VST  = os.path.join(BASE, "data", "nbl",  "nbl_vst_joint.rds")
    VST_ALLGENE_RBL_VST  = os.path.join(BASE, "data", "rbl",  "rbl_vst_joint.rds")

    # Joint metadata produced by the multicohort_cancer profile (sample_id,
    # cancer_type, sample_type[, cohort]). Used here purely to annotate
    # samples AFTER UMAP fitting.
    VST_ALLGENE_JOINT_META = os.path.join(
        BASE, "results", "unsupervised", "multicohort_cancer", "inputs",
        "joint_metadata.tsv"
    )

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

    rule plot_vst_all_gene_tumour_cell_line_alignment_umap:
        """
        Unsupervised diagnostic UMAP on the BRCA+NBL+RBL VST all-gene merged
        matrix. Uses the same plotting script as the pan-cancer DEG-set
        thesis figure, invoked with --feature_mode=all_genes and the
        VST_ALL_GENE label so output filenames remain disjoint from the
        thesis-facing DEG-set outputs.
        """
        input:
            expr_rds = VST_ALLGENE_UMAP_EXPR_RDS,
            meta_tsv = VST_ALLGENE_UMAP_META_TSV,
            script   = os.path.join(SCRIPTS_DIR,
                "plot_pan_cancer_tumour_cell_line_alignment_umap.R")
        output:
            pdf     = VST_ALLGENE_UMAP_PDFS,
            svg     = VST_ALLGENE_UMAP_SVGS,
            png     = VST_ALLGENE_UMAP_PNGS,
            coords  = VST_ALLGENE_UMAP_COORDS,
            summary = VST_ALLGENE_UMAP_SUMMARY
        params:
            outdir           = VST_ALLGENE_UMAP_OUTDIR,
            metrics          = ",".join(VST_ALLGENE_UMAP_METRICS),
            feature_label    = VST_ALLGENE_UMAP_LABEL,
            out_stem         = VST_ALLGENE_UMAP_STEM,
            summary_basename = os.path.basename(VST_ALLGENE_UMAP_SUMMARY),
            page             = "thesis"
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
              --feature_mode all_genes \
              --feature_label "{params.feature_label}" \
              --out_stem "{params.out_stem}" \
              --summary_basename "{params.summary_basename}" \
              --outdir "{params.outdir}" \
              --dist_metrics "{params.metrics}" \
              --page "{params.page}" \
              > "{log}" 2>&1
            rm -rf "$TMPDIR"
            for f in {output.pdf} {output.svg} {output.png} {output.coords} "{output.summary}"; do
              test -s "$f" || (echo "ERROR: missing or empty $f" >&2; exit 1)
            done
            '''

    MAPPING_OUTDIR = abspath(MAPPING_CFG.get("output_dir", "results/unsupervised/pan_cancer/tumour_mapping"))

    def mapping_metrics_summary():
        return os.path.join(MAPPING_OUTDIR, "metrics_summary.tsv")

    # Precompute optional CLI fragment for gene-symbol-map.
    _MAPPING_GENE_MAP_RAW = MAPPING_CFG.get("gene_symbol_map", "")
    _MAPPING_GENE_MAP_ARG = f'--gene-symbol-map "{_MAPPING_GENE_MAP_RAW}"' if _MAPPING_GENE_MAP_RAW else ""

    rule score_tumour_cellline_mapping:
        """Score tumour-to-cell-line mappings using the combined expression object."""
        input:
            expr_rds = PAN_EXPR_RDS,
            genes = PAN_FEATURES_CLEAN
        output:
            metrics_summary = mapping_metrics_summary(),
            tumour_rankings = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_to_cellline_rankings.tsv"),
            tumour_summary = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_mapping_summary.tsv"),
            tumour_scores = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_scores.rds"),
            cellline_rankings = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_to_tumour_rankings.tsv"),
            cellline_summary = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_mapping_summary.tsv"),
            cellline_metrics = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_metrics_overall.tsv")
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


    PAN_CANCER_DIR = abspath("results/unsupervised/pan_cancer")
    PAN_GRAPH_DIR = os.path.join(PAN_CANCER_DIR, "graph")
    PAN_FIG_DIR = os.path.join(PAN_CANCER_DIR, "figures")
    PAN_COR_RDS = os.path.join(PAN_CANCER_DIR, "pan_cancer_cor.rds")

    rule compute_pan_cancer_correlation:
        input: expr = PAN_EXPR_RDS
        output: cor = PAN_COR_RDS
        params: script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_correlation.R"), method = "spearman"
        log: os.path.join(LOGROOT, "compute_pan_cancer_correlation.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {PAN_CANCER_DIR} && Rscript {params.script} --input {input.expr} --output {output.cor} --method {params.method} > {log} 2>&1 && test -s {output.cor}"

    rule build_pan_cancer_graph:
        input: cor = PAN_COR_RDS
        output:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            components = os.path.join(PAN_GRAPH_DIR, "pan_cancer_components.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "build_pan_cancer_graph.R"), outdir = PAN_GRAPH_DIR, k = 20, min_cor = 0.0
        log: os.path.join(LOGROOT, "build_pan_cancer_graph.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --cor-matrix {input.cor} --k {params.k} --min-cor {params.min_cor} --output-dir {params.outdir} > {log} 2>&1 && test -s {output.edges} && test -s {output.components}"

    rule compute_pan_cancer_communities:
        input:
            edges = os.path.join(PAN_GRAPH_DIR, "pan_cancer_graph_edges.tsv"),
            meta = PAN_EXPR_RDS
        output:
            communities = os.path.join(PAN_GRAPH_DIR, "pan_cancer_communities.tsv"),
            summary = os.path.join(PAN_GRAPH_DIR, "community_validation", "community_summary.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "compute_pan_cancer_communities.R"), outdir = os.path.join(PAN_GRAPH_DIR, "community_validation"), seed = 1
        log: os.path.join(LOGROOT, "compute_pan_cancer_communities.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --edges {input.edges} --meta {input.meta} --out {output.communities} --outdir {params.outdir} --method both --seed {params.seed} > {log} 2>&1 && test -s {output.communities}"

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
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --edges {input.edges} --components {input.components} --meta {input.meta} --outdir {params.outdir} --expected-gene-count 0 > {log} 2>&1 && printf 'figure_name\tscript\tcommand\tgit_commit\ttimestamp\tinput_files\toutput_files\tupstream_tables\tkey_parameters\tsoftware_versions\tfigure_type\tsource_pipeline_root\tcopied_to_thesis_path\tlegacy_source_path\tnotes\nFig_pan_cancer_graph.pdf\tscripts/plot_pan_cancer_graph.R\tRscript scripts/plot_pan_cancer_graph.R\tunavailable_not_git_worktree\tNA\t{input.edges};{input.components};{input.meta}\t{output.pdf};{output.components_pdf};{output.size_pdf}\t{input.edges};{input.components}\texpected_gene_count=infer_from_rds_genes\tR/igraph/ggplot2\tpan_cancer\t/work/ugbogu/pipeline\t\t\tTranscriptomic similarity network for prioritisation and neighbourhood assignment only\n' > {output.provenance} && test -s {output.pdf}"

    rule plot_cellline_centred_figures:
        """
        Legacy bar-chart summaries for cell-line-centred similarity mapping.
        Uses simple top-k fraction statistics from score_tumour_cellline_mapping.R;
        does NOT compute Precision@k bootstrap confidence intervals.

        THIS IS NOT the four-panel Precision@k / bootstrap figure.
        The verified Precision@k bootstrap figure is produced by the
        cellline_precision_at_k rule (script: cellline_centred_precision_at_k.R).

        These outputs are legacy products retained for compatibility.
        """
        input:
            metrics_summary = mapping_metrics_summary(),
            tumour_rank = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_to_cellline_rankings.tsv"),
            cell_rank = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_to_tumour_rankings.tsv"),
            cell_metrics = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity", "cellline_metrics_overall.tsv")
        output:
            topk_bar_chart = os.path.join(MAPPING_OUTDIR, "cellline_centred", "Fig_cellline_centred_tumour_retrieval.pdf"),
            confidence = os.path.join(MAPPING_OUTDIR, "cellline_centred", "Fig_cellline_centred_confidence_scores.pdf"),
            component = os.path.join(MAPPING_OUTDIR, "cellline_centred", "Fig_cellline_centred_component_mapping.pdf"),
            replicate = os.path.join(MAPPING_OUTDIR, "cellline_centred", "Fig_cellline_centred_replicate_sensitivity.pdf"),
            provenance = os.path.join(MAPPING_OUTDIR, "cellline_centred", "cellline_centred_figure_provenance.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "plot_cellline_centred_figures.R"), outdir = os.path.join(MAPPING_OUTDIR, "cellline_centred"), top_k = MAPPING_CFG.get("top_k", 10), top_n = 50
        log: os.path.join(LOGROOT, "plot_cellline_centred_figures.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --mapping-dir {MAPPING_OUTDIR} --outdir {params.outdir} --top-k {params.top_k} --top-n {params.top_n} > {log} 2>&1 && test -s {output.topk_bar_chart} && test -s {output.provenance}"

    rule cellline_precision_at_k:
        """
        Cell-line-centred Precision@k bootstrap analysis.

        For each biological cell-line group (replicate profiles collapsed by
        arithmetic mean Spearman correlation before ranking), computes Precision@k
        across multiple top-k thresholds.  Bootstrap confidence intervals use
        B = 2000 resamples (with replacement) over biological cell-line groups,
        with 2.5th and 97.5th empirical percentiles as bounds.

        Lineage-stratified Precision@k summaries are computed separately within
        each lineage (BRCA, NBL, RBL) and as an Overall summary.

        Produces the four-panel figure (Fig_cellline_centred_tumour_class_similarity.pdf):
          A  Top-1 lineage confusion matrix
          B  Mean Precision@k with 95% bootstrap CI ribbons per cohort
          C  Same-lineage tumour rank percentile distributions per cell line
          D  Lineage composition of the top-k tumour neighbourhood

        This is the verified Drive pipeline computation.  Bootstrap seed = 20260603.
        tumour_components.tsv is optional (Panel B unaffected when absent).

        Expected QC values:
          raw profile-level observations : 58
          biological cell-line groups    : 56 (after replicate profile aggregation)
          BRCA groups : 29 | NBL groups : 18 | RBL groups : 9
          tumour samples               : 1,128
        """
        input:
            c2t_long = os.path.join(MAPPING_OUTDIR, "cellline_to_tumour_similarity",
                                    "cellline_tumour_scores_long.tsv.gz"),
            t2c_long = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity",
                                    "tumour_cellline_scores_long.tsv.gz")
        output:
            precision_fig = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_centred_tumour_class_similarity.pdf"),
            precision_at_k_fig = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "Fig_cellline_centred_precision_at_k.pdf"),
            topk_metrics  = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_topk_metrics.tsv"),
            rank_summary  = os.path.join(MAPPING_OUTDIR,
                                "cellline_similarity_precision_bootstrap",
                                "cellline_centred_rank_summary.tsv"),
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
              > "{log}" 2>&1
            test -s "{output.precision_fig}" || \
              (echo "ERROR: missing {output.precision_fig}" >&2; exit 1)
            test -s "{output.precision_at_k_fig}" || \
              (echo "ERROR: missing {output.precision_at_k_fig}" >&2; exit 1)
            test -s "{output.topk_metrics}" || \
              (echo "ERROR: missing {output.topk_metrics}" >&2; exit 1)
            # Provenance record
            printf 'field\tvalue\n' > {output.provenance}
            printf 'script\t{params.script}\n' >> {output.provenance}
            printf 'mapping_dir\t{params.mapping_dir}\n' >> {output.provenance}
            printf 'outdir\t{params.outdir}\n' >> {output.provenance}
            printf 'timestamp\t%s\n' "$(date -Iseconds)" >> {output.provenance}
            printf 'source_pipeline\tJarvis:/work/ugbogu/pipeline\n' >> {output.provenance}
            printf 'bootstrap_resamples\t2000\n' >> {output.provenance}
            printf 'bootstrap_statistic\tmean_precision_at_k\n' >> {output.provenance}
            printf 'bootstrap_bounds\t2.5th_97.5th_percentile\n' >> {output.provenance}
            printf 'resampling_unit\tbiological_cell_line_group\n' >> {output.provenance}
            printf 'replicate_collapse\tmean_score_across_replicate_profiles\n' >> {output.provenance}
            printf 'bootstrap_seed\t20260603\n' >> {output.provenance}
            '''


    rule plot_ecdf_model_prioritisation:
        input:
            metrics_summary = mapping_metrics_summary(),
            scores = os.path.join(MAPPING_OUTDIR, "tumour_to_cellline_similarity", "tumour_cellline_scores.rds"),
            meta = PAN_EXPR_RDS
        output:
            pdf = os.path.join(PAN_FIG_DIR, "ecdf_plots", "ecdf_model_prioritisation_combined.pdf"),
            png = os.path.join(PAN_FIG_DIR, "ecdf_plots", "ecdf_model_prioritisation_combined.png"),
            svg = os.path.join(PAN_FIG_DIR, "ecdf_plots", "ecdf_model_prioritisation_combined.svg"),
            summary = os.path.join(PAN_FIG_DIR, "ecdf_plots", "model_prioritisation_rank_summary.tsv"),
            provenance = os.path.join(PAN_FIG_DIR, "ecdf_plots", "ecdf_model_prioritisation_combined_provenance.tsv")
        params: script = os.path.join(SCRIPTS_DIR, "plot_ecdf_rank_combined_pub.R"), outdir = os.path.join(PAN_FIG_DIR, "ecdf_plots"), pipeline_dir = PAN_CANCER_DIR, top_k = MAPPING_CFG.get("top_k", 10), n_show = 4
        log: os.path.join(LOGROOT, "plot_ecdf_model_prioritisation.log")
        conda: CONDA_ENV_R
        shell: "mkdir -p {params.outdir} && Rscript {params.script} --pipeline-dir {params.pipeline_dir} --score-rds {input.scores} --meta-rds {input.meta} --outdir {params.outdir} --n-show {params.n_show} --top-k {params.top_k} > {log} 2>&1 && test -s {output.pdf} && test -s {output.provenance}"


# ---------------------------------------------------------------------------
# Pan-cancer cell line similarity network -- BUILD pipeline
# Enabled via config key: pan_cancer_cell_line_similarity
# ---------------------------------------------------------------------------
_PAN_CELL_LINE_SIM_CFG = config.get("pan_cancer_cell_line_similarity", {})

if _PAN_CELL_LINE_SIM_CFG:
    _CL_SIM_EXPR   = abspath(_PAN_CELL_LINE_SIM_CFG.get("expr_rds", "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds"))
    _CL_SIM_DIR    = abspath(_PAN_CELL_LINE_SIM_CFG.get("output_dir", "results/unsupervised/pan_cancer/cell_line_similarity"))
    _CL_SIM_COR    = _PAN_CELL_LINE_SIM_CFG.get("correlation", "spearman")
    _CL_SIM_K      = _PAN_CELL_LINE_SIM_CFG.get("k", 20)
    _CL_SIM_EXPECTED_GENES = _PAN_CELL_LINE_SIM_CFG.get("expected_genes", 171)
    _CL_SIM_EXPECTED_NODES = _PAN_CELL_LINE_SIM_CFG.get("expected_nodes", 167)
    _CL_SIM_SEED   = _PAN_CELL_LINE_SIM_CFG.get("seed", 1)
    _CL_SIM_LEIDEN_RESOLUTION = _PAN_CELL_LINE_SIM_CFG.get("leiden_resolution", 1.0)
    _CL_SIM_EDGES  = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_graph_edges.tsv")
    _CL_SIM_META   = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_node_metadata.tsv")
    _CL_SIM_COMMS  = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_communities.tsv")
    _CL_SIM_LEIDEN_COMMS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_leiden_communities.tsv")
    _CL_SIM_COMMUNITY_METRICS = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_community_metrics.tsv")
    _CL_SIM_LINEAGE_DISCORDANT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_lineage_discordant_profiles.tsv")
    _CL_SIM_LAYOUT = os.path.join(_CL_SIM_DIR, "pan_cancer_cell_line_layout.tsv")
    _CL_SIM_VAL_DIR = os.path.join(_CL_SIM_DIR, "community_validation")

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
            expected_nodes = _CL_SIM_EXPECTED_NODES
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
              > "{log}" 2>&1
            test -s "{output.edges}"
            test -s "{output.metadata}"
            '''

    rule compute_pan_cancer_cell_line_communities:
        input:
            edges    = _CL_SIM_EDGES,
            metadata = _CL_SIM_META
        output:
            communities = _CL_SIM_COMMS,
            leiden_communities = _CL_SIM_LEIDEN_COMMS,
            community_metrics = _CL_SIM_COMMUNITY_METRICS,
            lineage_discordant = _CL_SIM_LINEAGE_DISCORDANT
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
              --lineage-discordant-out "{output.lineage_discordant}" \
              --seed {params.seed} \
              --leiden-resolution {params.leiden_resolution} \
              > "{log}" 2>&1
            test -s "{output.communities}"
            test -s "{output.leiden_communities}"
            test -s "{output.community_metrics}"
            test -s "{output.lineage_discordant}"
            '''

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


# ---------------------------------------------------------------------------
# Pan-cancer cell line similarity network -- two-panel figure
# Cell-line-to-cell-line transcriptomic similarity (not cell-line-to-tumour).
# Enabled and configured via config key: pan_cancer_cell_line_plot
#
# Required config keys:
#   pan_cancer_cell_line_plot:
#     edges_tsv:       path to edge list TSV  (columns: from, to, weight)
#     communities_tsv: path to community TSV  (columns: sample, component, lineage)
#     layout_tsv:      path to layout TSV     (columns: sample, x, y)
#     output_dir:      directory for output PDF and PNG
#
# Run manually:  snakemake plot_pan_cancer_cell_line_two_panel --cores 1
# Not part of the automated pipeline; do NOT add to PIPELINE_TARGET.
# ---------------------------------------------------------------------------
_PAN_CELL_LINE_PLOT_CFG = config.get("pan_cancer_cell_line_plot", {})

if _PAN_CELL_LINE_PLOT_CFG:
    _PLOT_EDGES  = abspath(_PAN_CELL_LINE_PLOT_CFG["edges_tsv"])
    _PLOT_COMMS  = abspath(_PAN_CELL_LINE_PLOT_CFG["communities_tsv"])
    _PLOT_LAYOUT = abspath(_PAN_CELL_LINE_PLOT_CFG["layout_tsv"])
    _PLOT_OUTDIR = abspath(_PAN_CELL_LINE_PLOT_CFG["output_dir"])

    rule plot_pan_cancer_cell_line_two_panel:
        """Two-panel figure: pan-cancer cell line transcriptomic similarity
        network (cell-line-to-cell-line, not cell-line-to-tumour).
        Panel A: lineage coloured by Okabe-Ito palette.
        Panel B: Louvain communities (Dark2) with convex hulls and centroid labels.
        Accepts any edge list, community table, and layout produced by this
        or any compatible pipeline run."""
        input:
            edges  = _PLOT_EDGES,
            comms  = _PLOT_COMMS,
            layout = _PLOT_LAYOUT
        output:
            pdf = os.path.join(_PLOT_OUTDIR, "Fig_pan_cancer_cell_line_similarity_network_lineage_community.pdf"),
            png = os.path.join(_PLOT_OUTDIR, "Fig_pan_cancer_cell_line_similarity_network_lineage_community.png")
        params:
            script = os.path.join(SCRIPTS_DIR, "plot_pan_cancer_two_panel.R"),
            outdir = _PLOT_OUTDIR
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
              > "{log}" 2>&1
            test -s "{output.pdf}" || (echo "ERROR: missing {output.pdf}" >&2; exit 1)
            test -s "{output.png}" || (echo "ERROR: missing {output.png}" >&2; exit 1)
            '''


# Final pipeline terminal target list.
PIPELINE_TARGET = [
    os.path.join(UNSUP_REL, "tumour_neighbourhoods", "final_consensus_all", "resolved_dsmz_neighbours.tsv"),
    P_CONS_THESIS_RESOLUTION_PREFIX + ".pdf",
    P_CONS_THESIS_RESOLUTION_PREFIX + ".png",
    P_CONS_THESIS_CONSENSUS_PREFIX + ".pdf",
    P_CONS_THESIS_CONSENSUS_PREFIX + ".png",
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
]
if DESEQ2_ENABLED:
    PIPELINE_TARGET.append(NODE_STATS_TSV)
    PIPELINE_TARGET.append(os.path.join(DESEQ2_COMP_DIR, ".done"))

if ENRICH_ENABLED:
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "query_sets", "query_manifest.tsv"))
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "query_sets", "skipped_queries.tsv"))
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "gprofiler", "corpus_manifest.tsv"))
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "gprofiler", "top_terms.tsv"))
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "gprofiler", "iea_sensitivity_summary.tsv"))
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "gprofiler", "gprofiler_version.tsv"))

if MARKER_POST_ENABLED:
    PIPELINE_TARGET.append(PAN_EXPR_RDS)
    PIPELINE_TARGET.append(mapping_metrics_summary())

if ENRICH_ENABLED:
    PIPELINE_TARGET.append(os.path.join(ENRICH_DIR_REL, "figures", "Fig_enrichment_top_terms_heatmap.pdf"))
if MARKER_POST_ENABLED:
    PIPELINE_TARGET.extend([
        os.path.join(PAN_FIG_DIR, "Fig_pan_cancer_graph.pdf"),
        os.path.join(PAN_FIG_DIR, "ecdf_plots", "ecdf_model_prioritisation_combined.pdf"),
        # Verified four-panel Precision@k bootstrap figure (cellline_precision_at_k rule).
        # Legacy bar-chart outputs in cellline_centred/ are retained by plot_cellline_centred_figures
        # rule but removed from PIPELINE_TARGET; they are not the verified computation.
        os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                     "Fig_cellline_centred_tumour_class_similarity.pdf"),
        os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                     "cellline_topk_metrics.tsv"),
        os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                     "cellline_centred_rank_summary.tsv"),
        os.path.join(MAPPING_OUTDIR, "cellline_similarity_precision_bootstrap",
                     "cellline_similarity_precision_provenance.tsv"),
    ])

if IS_MULTICOHORT_PROFILE:
    PIPELINE_TARGET.extend([
        P_CONS_LEIDEN_COMMUNITIES_TSV,
        P_CONS_LEIDEN_SUMMARY_TSV,
        P_CONS_LEIDEN_MODULARITY_TSV,
        P_CONS_LEIDEN_LAYOUT_TSV,
        P_CONS_LEIDEN_FIG_PREFIX + ".pdf",
        P_CONS_LEIDEN_FIG_PREFIX + ".png",
    ])
