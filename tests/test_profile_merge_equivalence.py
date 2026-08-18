"""Profile-merge equivalence between the Snakefile (Python) and lib_config.R.

Both read config/config.yaml and both feed the same analysis, so a profile must
resolve to one effective configuration. The two implementations are compared on
a fixture covering the cases that distinguish them:

  * scalar override
  * nested mapping override (sibling keys must survive)
  * sequence override (must replace, not be merged or ignored)
  * key present only in defaults (must survive)
  * key present only in the profile (must appear)
  * unknown profile (must fail in both)

The sequence case is the one that previously diverged: a YAML sequence parses to
an unnamed R list, so recursing on it iterated over NULL names and silently kept
the defaults value while Python replaced it. That would let orchestration and
implementation run with different scientific universes without raising.
"""

import json
import os
import subprocess
import textwrap
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
LIB_CONFIG = REPO_ROOT / "scripts" / "lib_config.R"

RSCRIPT = os.environ.get("PIPELINE_RSCRIPT", "Rscript")

FIXTURE = textwrap.dedent("""\
    defaults:
      scalar_key: 0.70
      defaults_only: kept
      nested:
        untouched: base_value
        overridden: base_value
      methods:
        - Variance
        - MAD
        - Entropy
      k_grid: [2, 3, 4, 5, 6, 7, 8]
    profiles:
      demo:
        scalar_key: 0.85
        profile_only: added
        nested:
          overridden: profile_value
        methods:
          - PCA
          - HVG
        k_grid: [2, 3]
      other:
        scalar_key: 0.10
    """)


# Mirrors the Snakefile's deep_merge() exactly.
def python_deep_merge(base, override):
    merged = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(merged.get(k), dict):
            merged[k] = python_deep_merge(merged[k], v)
        else:
            merged[k] = v
    return merged


def python_effective(config_path: Path, profile: str):
    cfg = yaml.safe_load(config_path.read_text())
    profiles = cfg.get("profiles", {})
    if profile not in profiles:
        raise KeyError(profile)
    merged = python_deep_merge(cfg.get("defaults", {}), profiles[profile])
    merged["profile"] = profile
    return merged


def r_effective(config_path: Path, profile: str):
    code = f'''
      suppressPackageStartupMessages({{ library(yaml); library(jsonlite) }})
      source("{LIB_CONFIG.as_posix()}")
      cfg <- read_profiled_config("{config_path.as_posix()}", profile_override = "{profile}")
      cat(jsonlite::toJSON(cfg, auto_unbox = TRUE, digits = NA))
    '''
    res = subprocess.run([RSCRIPT, "-e", code], capture_output=True, text=True, check=True)
    return json.loads(res.stdout)


def r_effective_expect_error(config_path: Path, profile: str) -> str:
    code = f'''
      suppressPackageStartupMessages({{ library(yaml) }})
      source("{LIB_CONFIG.as_posix()}")
      read_profiled_config("{config_path.as_posix()}", profile_override = "{profile}")
    '''
    res = subprocess.run([RSCRIPT, "-e", code], capture_output=True, text=True)
    assert res.returncode != 0, "unknown profile must fail"
    return res.stderr


def write_fixture(tmp_path: Path) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / "merge_fixture.yaml"
    p.write_text(FIXTURE)
    return p


def normalise(value):
    """JSON round-trip so R and Python values compare structurally."""
    return json.loads(json.dumps(value))


def test_python_and_r_agree_on_effective_config(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    py = python_effective(cfg_path, "demo")
    r = r_effective(cfg_path, "demo")

    for key in ("scalar_key", "defaults_only", "profile_only", "nested", "methods", "k_grid"):
        assert normalise(py[key]) == normalise(r[key]), (
            f"key '{key}' differs: python={py[key]!r} r={r[key]!r}"
        )


def test_scalar_override(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    assert python_effective(cfg_path, "demo")["scalar_key"] == 0.85
    assert r_effective(cfg_path, "demo")["scalar_key"] == 0.85


def test_nested_mapping_override_keeps_siblings(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    for effective in (python_effective(cfg_path, "demo"), r_effective(cfg_path, "demo")):
        assert effective["nested"]["overridden"] == "profile_value"
        assert effective["nested"]["untouched"] == "base_value", (
            "a nested mapping override must not drop sibling keys"
        )


def test_sequence_override_replaces(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    for label, effective in (("python", python_effective(cfg_path, "demo")),
                             ("R", r_effective(cfg_path, "demo"))):
        assert normalise(effective["methods"]) == ["PCA", "HVG"], (
            f"{label}: a profile sequence must replace the defaults sequence, "
            f"got {effective['methods']!r}"
        )
        assert normalise(effective["k_grid"]) == [2, 3], (
            f"{label}: k_grid must be replaced, got {effective['k_grid']!r}"
        )


def test_missing_profile_key_falls_back_to_defaults(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    for effective in (python_effective(cfg_path, "other"), r_effective(cfg_path, "other")):
        assert effective["scalar_key"] == 0.10
        assert normalise(effective["methods"]) == ["Variance", "MAD", "Entropy"]
        assert effective["defaults_only"] == "kept"


def test_unknown_profile_fails_in_both(tmp_path=None):
    cfg_path = write_fixture(tmp_path or Path("/tmp/merge_equivalence"))
    try:
        python_effective(cfg_path, "does_not_exist")
        raise AssertionError("python must reject an unknown profile")
    except KeyError:
        pass
    err = r_effective_expect_error(cfg_path, "does_not_exist")
    assert "Profile not found" in err, err


def test_real_config_agrees_for_every_profile():
    """The shipped configuration must resolve identically in both implementations."""
    real = REPO_ROOT / "config" / "config.yaml"
    cfg = yaml.safe_load(real.read_text())
    for profile in cfg.get("profiles", {}):
        py = python_effective(real, profile)
        r = r_effective(real, profile)
        for key in ("patient_referenced_graph", "tumour_neighbourhoods",
                    "hclust_kmeans", "clustering", "feature_sets",
                    "feature_selection"):
            if key in py or key in r:
                assert normalise(py.get(key)) == normalise(r.get(key)), (
                    f"profile '{profile}' key '{key}' differs between "
                    f"Snakemake and R merge"
                )


if __name__ == "__main__":
    base = Path("/tmp/merge_equivalence")
    test_python_and_r_agree_on_effective_config(base)
    print("OK python and R produce the same effective configuration")
    test_scalar_override(base)
    print("OK scalar override")
    test_nested_mapping_override_keeps_siblings(base)
    print("OK nested mapping override keeps sibling keys")
    test_sequence_override_replaces(base)
    print("OK sequence override replaces rather than being merged or ignored")
    test_missing_profile_key_falls_back_to_defaults(base)
    print("OK key absent from the profile falls back to defaults")
    test_unknown_profile_fails_in_both(base)
    print("OK unknown profile fails in both implementations")
    test_real_config_agrees_for_every_profile()
    print("OK shipped config resolves identically for every profile")
    print("PASS")
