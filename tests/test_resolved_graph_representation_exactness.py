"""Regression tests: resolved-neighbour construction validates its inputs exactly.

resolve_dsmz_graph_neighbours.R used to rebuild the representation universe from
config sections that the multicohort profile does not declare, fall back to a
built-in nine-method by two-distance grid, and downgrade a missing graph product
to a warning. The combination silently resolved the multicohort graph over 18
representations while configuration declared 20, dropping
PanCancerFeatureSet_euc and PanCancerFeatureSet_corr without any failure.

These tests drive the real script over a synthetic graph root:

  * the complete expected set resolves and writes a validation report;
  * a missing expected representation is an error, not a warning;
  * a representation present on disk but not expected is never consumed;
  * an absent --directions argument is an error, so no built-in universe exists.

Nothing here reads or writes production results.
"""

import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RESOLVE_SCRIPT = REPO_ROOT / "scripts" / "resolve_dsmz_graph_neighbours.R"
RSCRIPT = os.environ.get("PIPELINE_RSCRIPT", "Rscript")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        shutil.which(RSCRIPT) is None,
        reason="Rscript is not on PATH; set PIPELINE_RSCRIPT to the workflow R",
    ),
]

REPRESENTATIONS = [
    "Variance_euc",
    "Variance_corr",
    "PanCancerFeatureSet_euc",
    "PanCancerFeatureSet_corr",
]

CONFIG_FIXTURE = textwrap.dedent("""\
    defaults:
      paths:
        unsup_root: results/unsupervised/fixture
    profiles:
      fixture:
        paths:
          unsup_root: results/unsupervised/fixture
    """)


def write_graph_root(root, representations):
    """Create one canonical edge file per representation."""
    for representation in representations:
        directory = root / representation / "final_consensus"
        directory.mkdir(parents=True, exist_ok=True)
        edges = directory / f"cell_line_similarity_graph_edges_{representation}.tsv"
        edges.write_text(
            "cell_line1\tcell_line2\tsimilarity\n"
            "CL_A\tCL_B\t0.9\n"
            "CL_B\tCL_C\t0.8\n"
        )


def write_summaries(directory, representations):
    """Write the winners and direction-summary tables the script consumes."""
    directory.mkdir(parents=True, exist_ok=True)
    winners = directory / "p_consensus_winners_by_frac_ge_thr.tsv"
    winners.write_text(
        "cell_line\tbest_dir_frac_ge_thr\n"
        f"CL_A\t{representations[0]}\n"
        f"CL_B\t{representations[0]}\n"
        f"CL_C\t{representations[0]}\n"
    )
    direction_summary = directory / "p_consensus_direction_summary.tsv"
    rows = ["direction\tfrac_ge_thr\tmedian_p_consensus\tmean_p_consensus"]
    for index, representation in enumerate(representations):
        score = 0.9 - 0.01 * index
        rows.append(f"{representation}\t{score}\t{score}\t{score}")
    direction_summary.write_text("\n".join(rows) + "\n")
    return winners, direction_summary


def run_resolution(tmp_path, expected_representations, on_disk_representations,
                   pass_directions=True):
    graph_root = tmp_path / "graphs"
    write_graph_root(graph_root, on_disk_representations)
    summary_dir = tmp_path / "summaries"
    winners, direction_summary = write_summaries(summary_dir, expected_representations)

    config_file = tmp_path / "config" / "config.yaml"
    config_file.parent.mkdir(parents=True, exist_ok=True)
    config_file.write_text(CONFIG_FIXTURE)

    command = [
        RSCRIPT,
        str(RESOLVE_SCRIPT),
        "--config", str(config_file),
        "--profile", "fixture",
        "--winners_tsv", str(winners),
        "--direction_summary_tsv", str(direction_summary),
        "--graph_root", str(graph_root),
        "--output_tsv", str(tmp_path / "resolved_dsmz_neighbours.tsv"),
        "--validation_tsv", str(tmp_path / "resolved_graph_input_validation.tsv"),
    ]
    if pass_directions:
        command += ["--directions", ",".join(expected_representations)]

    return subprocess.run(command, capture_output=True, text=True)


def test_complete_representation_set_resolves(tmp_path):
    completed = run_resolution(tmp_path, REPRESENTATIONS, REPRESENTATIONS)
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output

    resolved = tmp_path / "resolved_dsmz_neighbours.tsv"
    assert resolved.exists() and resolved.stat().st_size > 0

    report = (tmp_path / "resolved_graph_input_validation.tsv").read_text()
    for representation in REPRESENTATIONS:
        assert representation in report
    assert "expected_present" in report
    assert "expected_missing" not in report


def test_missing_pan_cancer_feature_set_representation_is_an_error(tmp_path):
    """The exact failure the 18-representation fallback used to hide."""
    on_disk = [r for r in REPRESENTATIONS if not r.startswith("PanCancerFeatureSet")]
    completed = run_resolution(tmp_path, REPRESENTATIONS, on_disk)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "PanCancerFeatureSet_euc" in output
    assert "PanCancerFeatureSet_corr" in output
    assert not (tmp_path / "resolved_dsmz_neighbours.tsv").exists()


def test_unexpected_representation_is_recorded_but_not_consumed(tmp_path):
    on_disk = REPRESENTATIONS + ["StaleMethod_euc"]
    completed = run_resolution(tmp_path, REPRESENTATIONS, on_disk)
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output

    report = (tmp_path / "resolved_graph_input_validation.tsv").read_text()
    assert "StaleMethod_euc" in report
    assert "unexpected_present_not_consumed" in report

    resolved = (tmp_path / "resolved_dsmz_neighbours.tsv").read_text()
    assert "StaleMethod_euc" not in resolved


def test_absent_directions_argument_is_an_error(tmp_path):
    """There is no built-in representation universe to fall back to."""
    completed = run_resolution(
        tmp_path, REPRESENTATIONS, REPRESENTATIONS, pass_directions=False
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "--directions is required" in output
