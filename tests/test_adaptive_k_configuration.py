"""Regression tests for adaptive tumour-neighbourhood sizing.

The adaptive neighbourhood rule is

    k_i = min( maximum, max( minimum, ceil(fraction * n_i) ), n_i )

where n_i is the number of candidate tumours for cell line i. The parameters
must come from configuration (tumour_neighbourhoods.adaptive_k), not from
values hard-coded at the call site.

These tests read a temporary, deliberately non-production configuration with
the same reader the pipeline uses, feed the resulting values into the real
neighbourhood function, and check that the computed neighbourhood sizes follow
the configured values. Production configuration is never modified.
"""

import os
import subprocess
import textwrap
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
NH_CORE = REPO_ROOT / "R" / "base_functions" / "tumour_neighbourhood.R"
LIB_CONFIG = REPO_ROOT / "scripts" / "lib_config.R"
COMPUTE_SCRIPT = REPO_ROOT / "scripts" / "compute_tumour_neighbourhoods.R"
PRODUCTION_CONFIG = REPO_ROOT / "config" / "config.yaml"

RSCRIPT = os.environ.get("PIPELINE_RSCRIPT", "Rscript")

# One cell line and 40 tumours, all in a single cluster, so the candidate set
# size n_i is exactly 40 and k depends only on the configured parameters.
N_TUMOURS = 40

PREAMBLE = f"""
suppressPackageStartupMessages({{
  library(dplyr)
  library(yaml)
}})
source("{LIB_CONFIG.as_posix()}")
source("{NH_CORE.as_posix()}")

make_inputs <- function(n_tumour = {N_TUMOURS}) {{
  set.seed(42)
  ids <- c("CELL_1", paste0("TUM_", seq_len(n_tumour)))
  m <- matrix(rnorm(length(ids) * 10), nrow = length(ids),
              dimnames = list(ids, paste0("g", 1:10)))
  cluster <- setNames(rep("C1", length(ids)), ids)
  dataset <- setNames(c("DSMZ", rep("TUMOUR", n_tumour)), ids)
  list(m = m, cluster = cluster, dataset = dataset)
}}

k_for <- function(cfg_path, n_tumour = {N_TUMOURS}) {{
  cfg <- read_profiled_config(cfg_path, profile_override = "brca")
  a <- cfg$tumour_neighbourhoods$adaptive_k
  # the same three reads the active call site performs
  top_frac  <- as.numeric(a$fraction)
  top_n_min <- as.integer(a$minimum)
  top_n_max <- as.integer(a$maximum)
  inp <- make_inputs(n_tumour)
  res <- compute_tumour_neighbourhoods(
    emb_mat = inp$m, cluster_m = inp$cluster, dataset = inp$dataset,
    top_frac = top_frac, top_n_min = top_n_min, top_n_max = top_n_max,
    method_id = "test", distance = "euclidean"
  )
  length(res$neighbourhoods[["CELL_1"]])
}}
"""


def write_temp_config(tmp_path: Path, fraction, minimum, maximum) -> Path:
    """A minimal profiled config carrying only what the reader needs."""
    cfg = textwrap.dedent(f"""\
        defaults:
          tumour_neighbourhoods:
            adaptive_k:
              fraction: {fraction}
              minimum: {minimum}
              maximum: {maximum}
        profiles:
          brca:
            analysis:
              cancer_type: BRCA
        """)
    p = tmp_path / "adaptive_k_test_config.yaml"
    p.write_text(cfg)
    return p


def run_r(code: str) -> list[str]:
    result = subprocess.run(
        [RSCRIPT, "-e", PREAMBLE + code],
        check=True, capture_output=True, text=True,
    )
    return [l.strip() for l in result.stdout.splitlines()
            if l.strip() and not l.startswith("[INFO]")]


def test_configured_values_propagate_into_computed_k(tmp_path=None):
    tmp_path = tmp_path or Path("/tmp")

    # Non-production: fraction 0.50, minimum 5, maximum 7.
    # ceil(0.50 * 40) = 20 -> max(5, 20) = 20 -> min(20, 7, 40) = 7
    non_prod = write_temp_config(tmp_path, 0.50, 5, 7)
    # Production: fraction 0.10, minimum 30, maximum 200.
    # ceil(0.10 * 40) = 4 -> max(30, 4) = 30 -> min(30, 200, 40) = 30
    prod_dir = tmp_path / "prod"
    prod_dir.mkdir(exist_ok=True)
    prod = write_temp_config(prod_dir, 0.10, 30, 200)

    lines = run_r(f'cat(k_for("{non_prod.as_posix()}"), k_for("{prod.as_posix()}"), sep = "\\n")')
    k_non_prod, k_prod = int(lines[0]), int(lines[1])

    assert k_non_prod == 7, (
        f"non-production config (0.50/5/7) must give k=7 for n=40, got {k_non_prod}"
    )
    assert k_prod == 30, (
        f"production config (0.10/30/200) must give k=30 for n=40, got {k_prod}"
    )
    assert k_non_prod != k_prod, "configuration change must alter the computed k"


def test_k_is_clamped_by_candidate_count(tmp_path=None):
    """With fewer candidates than the configured minimum, k is the candidate
    count: k = min(maximum, max(minimum, ceil(fraction*n)), n)."""
    tmp_path = tmp_path or Path("/tmp")
    d = tmp_path / "clamp"
    d.mkdir(exist_ok=True)
    prod = write_temp_config(d, 0.10, 30, 200)
    lines = run_r(f'cat(k_for("{prod.as_posix()}", n_tumour = 20), sep = "\\n")')
    assert int(lines[0]) == 20, f"k must clamp to the 20 available tumours, got {lines[0]}"


def _scaffold_repo(root: Path, fraction, minimum, maximum) -> Path:
    """Minimal repo layout so the script follows its production config path.

    compute_tumour_neighbourhoods.R resolves scripts/lib_config.R and the
    helper directory relative to the config file's parent, so those are
    symlinked to the real repository.
    """
    root.mkdir(parents=True, exist_ok=True)
    (root / "config").mkdir(exist_ok=True)
    for name in ("scripts", "R"):
        link = root / name
        if not link.exists():
            link.symlink_to(REPO_ROOT / name)

    cfg = textwrap.dedent(f"""\
        defaults:
          paths:
            tumour_nh_base_functions: R/base_functions
            unsup_root: results/unsupervised
          tumour_neighbourhoods:
            adaptive_k:
              fraction: {fraction}
              minimum: {minimum}
              maximum: {maximum}
        profiles:
          brca:
            analysis:
              cancer_type: BRCA
        """)
    cfg_path = root / "config" / "config.yaml"
    cfg_path.write_text(cfg)
    return cfg_path


def test_invalid_configuration_is_rejected_by_the_active_call_site(tmp_path=None):
    """An out-of-range fraction must abort the real script with an explicit
    message, proving the active call site reads and validates configuration
    rather than using hard-coded values."""
    tmp_path = tmp_path or Path("/tmp")
    bad = _scaffold_repo(tmp_path / "invalid_fraction", 1.5, 30, 200)

    result = subprocess.run(
        [RSCRIPT, COMPUTE_SCRIPT.as_posix(),
         "--config", bad.as_posix(), "--profile", "brca",
         "--direction", "Variance_euc", "--cluster-family", "hc"],
        capture_output=True, text=True,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0, "an out-of-range fraction must abort the run"
    assert "adaptive_k.fraction must satisfy" in combined, (
        f"expected the explicit adaptive_k validation error, got:\n{combined[-2000:]}"
    )


def test_missing_adaptive_k_section_is_rejected(tmp_path=None):
    """Removing the configuration section must abort rather than silently
    falling back to built-in defaults."""
    tmp_path = tmp_path or Path("/tmp")
    root = tmp_path / "missing_section"
    _scaffold_repo(root, 0.10, 30, 200)
    cfg_path = root / "config" / "config.yaml"
    cfg_path.write_text(textwrap.dedent("""\
        defaults:
          paths:
            tumour_nh_base_functions: R/base_functions
            unsup_root: results/unsupervised
        profiles:
          brca:
            analysis:
              cancer_type: BRCA
        """))

    result = subprocess.run(
        [RSCRIPT, COMPUTE_SCRIPT.as_posix(),
         "--config", cfg_path.as_posix(), "--profile", "brca",
         "--direction", "Variance_euc", "--cluster-family", "hc"],
        capture_output=True, text=True,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0, "a missing adaptive_k section must abort the run"
    assert "tumour_neighbourhoods.adaptive_k" in combined, (
        f"expected an explicit missing-section error, got:\n{combined[-2000:]}"
    )


def test_production_values_are_unchanged():
    """Guard: the production configuration must still carry 0.10 / 30 / 200."""
    import yaml
    cfg = yaml.safe_load(PRODUCTION_CONFIG.read_text())
    a = cfg["defaults"]["tumour_neighbourhoods"]["adaptive_k"]
    assert float(a["fraction"]) == 0.10, a
    assert int(a["minimum"]) == 30, a
    assert int(a["maximum"]) == 200, a


if __name__ == "__main__":
    base = Path("/tmp/adaptive_k_test")
    base.mkdir(parents=True, exist_ok=True)
    test_configured_values_propagate_into_computed_k(base)
    print("OK configured adaptive_k values propagate into computed k (7 vs 30 for n=40)")
    test_k_is_clamped_by_candidate_count(base)
    print("OK k clamps to the available candidate count")
    test_invalid_configuration_is_rejected_by_the_active_call_site(base)
    print("OK active call site rejects an out-of-range configured fraction")
    test_missing_adaptive_k_section_is_rejected(base)
    print("OK active call site rejects a missing adaptive_k section")
    test_production_values_are_unchanged()
    print("OK production values remain 0.10 / 30 / 200")
    print("PASS")
