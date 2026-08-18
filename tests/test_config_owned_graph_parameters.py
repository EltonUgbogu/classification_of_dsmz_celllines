"""Regression tests: config/config.yaml owns the graph-stage scientific choices.

Every parameter checked here was previously owned twice — once in
config/config.yaml and once as a literal or default inside orchestration or an
implementation script. A second owner is invisible: the configuration keeps
declaring a value while the DAG runs on the embedded one, so changing
configuration changes nothing.

The tests therefore assert *effectiveness*, not presence: an override applied on
top of the production configuration must change what the workflow plans to run.
Snakemake merges an extra ``--configfile`` over the workflow's own configfile,
so the production configuration is read but never modified.

The planned commands are read from a dry run (``-n -p``), which needs no conda
environments and executes nothing.
"""

import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = os.environ.get("PIPELINE_SNAKEMAKE", "snakemake")

pytestmark = pytest.mark.skipif(
    shutil.which(SNAKEMAKE) is None,
    reason=(
        "snakemake is not on PATH; set PIPELINE_SNAKEMAKE to the interpreter "
        "that runs this workflow"
    ),
)

# Dry runs traverse the whole DAG for a cohort profile; they are slow enough to
# deserve marking, but they are the only way to observe what the workflow would
# actually execute.
slow = pytest.mark.slow


def dry_run(tmp_path, profile, overrides=None, expect_success=True):
    """Plan `profile` with optional configuration overrides; return the output."""
    command = [
        SNAKEMAKE,
        "-n",
        "-p",
        "--directory",
        str(REPO_ROOT),
        "--snakefile",
        str(REPO_ROOT / "Snakefile"),
    ]
    if overrides is not None:
        override_file = tmp_path / "override.yaml"
        override_file.write_text(textwrap.dedent(overrides))
        command += ["--configfile", str(override_file)]
    command += ["--config", f"pipeline_profile={profile}"]

    completed = subprocess.run(
        command,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    if expect_success:
        assert completed.returncode == 0, output
    else:
        assert completed.returncode != 0, (
            "expected the workflow to refuse this configuration, but it planned "
            "successfully:\n" + output
        )
    return output


# ---------------------------------------------------------------------------
# similarity_quantile
# ---------------------------------------------------------------------------

@slow
def test_similarity_quantile_reaches_the_default_graph_branch(tmp_path):
    """The default graph rule must thread the configured quantile, not 0.9.

    The default branch previously called quantile(..., 0.9) inside
    compute_cell_line_similarity.R, so patient_referenced_graph.similarity_quantile
    had no effect on the graph the pipeline actually built.
    """
    output = dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            similarity_quantile: 0.55
        """,
    )
    assert "--similarity_quantile 0.55" in output
    assert "--similarity_quantile 0.9" not in output


@slow
def test_similarity_quantile_outside_the_unit_interval_is_refused(tmp_path):
    """An inadmissible quantile must fail at parse time, not silently."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            similarity_quantile: 1.5
        """,
        expect_success=False,
    )


def test_no_hard_coded_graph_quantile_remains():
    """compute_cell_line_similarity.R must not own a graph quantile literal."""
    source = (REPO_ROOT / "scripts" / "compute_cell_line_similarity.R").read_text()
    assert "quantile(sim_long$similarity, similarity_quantile" in source
    assert "quantile(sim_long$similarity, 0.9" not in source


# ---------------------------------------------------------------------------
# p_consensus_threshold
# ---------------------------------------------------------------------------

@slow
def test_p_consensus_threshold_reaches_every_active_reporting_path(tmp_path):
    """Reporting stages must receive the configured threshold, not a literal 0.7."""
    output = dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            p_consensus_threshold: 0.55
        """,
    )
    assert "--threshold 0.55" in output
    assert "--consensus_threshold 0.55" in output
    assert "--p-consensus-threshold 0.55" in output
    assert "--threshold 0.7" not in output
    assert "--consensus_threshold 0.7" not in output


def test_no_active_script_emits_a_threshold_encoding_column():
    """frac_ge_0_7 encodes a threshold configuration is allowed to change."""
    active_sources = [
        REPO_ROOT / "scripts" / "tumour_neighbourhood_p_consensus.R",
        REPO_ROOT / "scripts" / "summarize_p_consensus_all.R",
        REPO_ROOT / "scripts" / "compute_cell_line_similarity.R",
    ]
    offenders = [
        path.name for path in active_sources if "frac_ge_0_7" in path.read_text()
    ]
    assert offenders == [], f"threshold-encoding column name still present in {offenders}"


def test_study_design_manifest_does_not_restate_the_threshold():
    """study_design.yaml must not be a second owner of the consensus threshold."""
    study_design = (REPO_ROOT / "config" / "study_design.yaml").read_text()
    assert "strong_support_threshold:" not in study_design


# ---------------------------------------------------------------------------
# similarity_metrics
# ---------------------------------------------------------------------------

@slow
def test_default_targets_expand_over_configured_similarity_metrics(tmp_path):
    """Both configured metrics must appear in the default DAG.

    The default branch used to be Pearson-only while config declared
    pearson and jaccard, so the configured metric list was inert.
    """
    output = dry_run(tmp_path, "brca")
    assert "cell_line_similarity_graph_metric" in output
    assert "patient_referenced_similarity_metrics/pearson" in output
    assert "patient_referenced_similarity_metrics/jaccard" in output


@slow
def test_primary_similarity_metric_must_be_a_configured_metric(tmp_path):
    """The primary branch metric is selected by config and validated against it."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            primary_similarity_metric: spearman
        """,
        expect_success=False,
    )


@slow
def test_primary_metric_is_transmitted_to_the_default_graph_rule(tmp_path):
    """Pearson must be an explicit configured selection, not a script default."""
    output = dry_run(tmp_path, "brca")
    assert "--similarity_metric pearson" in output


# ---------------------------------------------------------------------------
# Feature-method and clustering-formulation universes
# ---------------------------------------------------------------------------

@slow
def test_enabled_feature_method_without_a_top_n_is_refused(tmp_path):
    """Every enabled method must declare its top-N; there is no built-in table."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        profiles:
          brca:
            feature_sets:
              methods:
                - Variance
                - MAD
                - MeanAbsDev
                - Entropy
                - PCA
                - Spearman
                - MX
                - kTotal
                - HVG
                - UndeclaredMethod
        """,
        expect_success=False,
    )


@slow
def test_removing_a_configured_representation_changes_expected_targets(tmp_path):
    """The representation universe is a configuration declaration."""
    baseline = dry_run(tmp_path, "brca")
    reduced = dry_run(
        tmp_path,
        "brca",
        overrides="""
        profiles:
          brca:
            tumour_neighbourhoods:
              directions:
                - Variance_euc
                - Variance_corr
                - MAD_euc
                - MAD_corr
              additional_summary_directions: []
        """,
    )
    assert "MX_corr" in baseline
    assert "MX_corr" not in reduced


@slow
def test_retired_clustering_formulation_identifiers_are_refused(tmp_path):
    """AGN_*/CCP_HC_*/CCP_KM_* must fail, not resolve to an unproduced directory."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            clustering_methods_by_distance:
              corr:
                - CCP_HC_expr_cell_tumour
        """,
        expect_success=False,
    )


@slow
def test_clustering_formulation_without_a_declared_product_is_refused(tmp_path):
    """A configured formulation must map onto a clustering product the DAG builds."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            clustering_methods_by_distance:
              corr:
                - HCLUST_expr_cell_only
        """,
        expect_success=False,
    )


@slow
def test_kmeans_formulation_under_correlation_distance_is_refused(tmp_path):
    """k-means is declared for Euclidean representations only."""
    dry_run(
        tmp_path,
        "brca",
        overrides="""
        defaults:
          patient_referenced_graph:
            clustering_methods_by_distance:
              corr:
                - HCLUST_expr_cell_tumour
                - KMEANS_expr_cell_tumour
        """,
        expect_success=False,
    )


# ---------------------------------------------------------------------------
# Multicohort representation completeness
# ---------------------------------------------------------------------------

@slow
def test_multicohort_plans_all_configured_representations(tmp_path):
    """multicohort_cancer declares 20 representations; all 20 must be planned."""
    output = dry_run(tmp_path, "multicohort_cancer")
    assert "PanCancerFeatureSet_euc" in output
    assert "PanCancerFeatureSet_corr" in output
    for line in output.splitlines():
        if line.strip().startswith("cell_line_similarity_graph "):
            assert line.split()[-1] == "20", line
            break
    else:  # pragma: no cover - the rule is always planned for this profile
        pytest.fail("cell_line_similarity_graph absent from the multicohort plan")


@slow
def test_multicohort_builds_the_pan_cancer_feature_panel_first(tmp_path):
    """PanCancerFeatureSet representations must depend on the marker-derived panel.

    Without the producing rule the panel was an ordinary file: present on disk
    the DAG succeeded, absent it failed late, and either way nothing tied the
    representation to the panel it is made of.
    """
    output = dry_run(tmp_path, "multicohort_cancer")
    assert "construct_pan_cancer_feature_panel" in output
    assert "pan_cancer_features_clean.txt" in output


@slow
def test_resolution_receives_the_configured_representation_list(tmp_path):
    """resolve_dsmz_graph_neighbours must be handed the expected representations."""
    output = dry_run(tmp_path, "multicohort_cancer")
    directions_arguments = [
        line for line in output.splitlines() if "--directions" in line
    ]
    assert directions_arguments, "no --directions argument reached the resolution stage"
    assert any("PanCancerFeatureSet_euc" in line for line in directions_arguments)
    assert any("PanCancerFeatureSet_corr" in line for line in directions_arguments)


def test_resolution_script_cannot_reconstruct_the_representation_universe():
    """The script must not hold a fallback representation grid."""
    source = (REPO_ROOT / "scripts" / "resolve_dsmz_graph_neighbours.R").read_text()
    assert "--directions is required" in source
    # The retired fallback built directions from a built-in nine-method grid.
    assert 'c("Variance", "MAD", "MeanAbsDev"' not in source
    assert '"Variance_euc", "Variance_corr"' not in source
