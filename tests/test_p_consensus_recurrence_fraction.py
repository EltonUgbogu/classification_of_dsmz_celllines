"""Regression tests for the p_consensus recurrence fraction.

p_consensus(i, t | r) must equal the raw recurrence count of tumour t in the
method-specific tumour neighbourhoods of cell line i, divided by the exact
configured number of eligible clustering formulations |A_r|
(patient_referenced_graph.clustering_methods_by_distance):

    Euclidean directions:   recurrence / 8 == p_consensus
    Correlation directions: recurrence / 4 == p_consensus

Eligible formulations are the JOINT cell-line + tumour outputs of the
HC/k-means clustering branch (AGN_*) and of ConsensusClusterPlus (CCP_*).

Config-level tests run anywhere. Output-level tests reconstruct p_consensus
from the per-method Top_m_long_<method_id>.csv tables and compare against the
reported Final_consensus_tumour_neighbourhoods_<direction>.tsv; they are
skipped when the profile's results tree is not present (e.g. code-only
checkouts).

Environment overrides:
    P_CONSENSUS_TEST_PROFILE     profile to check (default: brca)
    P_CONSENSUS_TEST_DIRECTIONS  comma-separated directions
                                 (default: Variance_euc,Variance_corr)
"""

import csv
import os
import sys
from collections import Counter
from pathlib import Path

import yaml

try:
    import pytest
except ImportError:  # standalone execution without pytest
    pytest = None

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"

PROFILE = os.environ.get("P_CONSENSUS_TEST_PROFILE", "brca")
DIRECTIONS = [
    d.strip()
    for d in os.environ.get(
        "P_CONSENSUS_TEST_DIRECTIONS", "Variance_euc,Variance_corr"
    ).split(",")
    if d.strip()
]

EXPECTED_METHODS = {
    "corr": [
        "AGN_HC_expr_cell_tumour",
        "AGN_HC_pca_cell_tumour",
        "CCP_HC_expr_cell_tumour",
        "CCP_HC_pca_cell_tumour",
    ],
    "euc": [
        "AGN_HC_expr_cell_tumour",
        "AGN_HC_pca_cell_tumour",
        "AGN_KM_expr_cell_tumour",
        "AGN_KM_pca_cell_tumour",
        "CCP_HC_expr_cell_tumour",
        "CCP_HC_pca_cell_tumour",
        "CCP_KM_expr_cell_tumour",
        "CCP_KM_pca_cell_tumour",
    ],
}


def load_config():
    return yaml.safe_load(CONFIG_PATH.read_text())


def declared_methods_by_distance(cfg):
    graph_cfg = cfg["defaults"]["patient_referenced_graph"]
    return graph_cfg["clustering_methods_by_distance"]


def profile_unsup_root(cfg, profile):
    prof = cfg["profiles"][profile]
    return REPO_ROOT / prof["paths"]["unsup_root"]


def skip_or_fail(message):
    if pytest is not None:
        pytest.skip(message)
    else:
        print(f"SKIP: {message}")
        return True


# ---------------------------------------------------------------------------
# Config-level tests: exact configured method sets
# ---------------------------------------------------------------------------

def test_configured_method_sets_are_exact():
    declared = declared_methods_by_distance(load_config())
    assert sorted(declared) == ["corr", "euc"]
    assert declared["corr"] == EXPECTED_METHODS["corr"], (
        "corr must declare exactly the 4 JOINT HC formulations "
        "(HC/k-means and CCP, expression and PCA space)"
    )
    assert declared["euc"] == EXPECTED_METHODS["euc"], (
        "euc must declare exactly the 8 JOINT formulations "
        "(HC and k-means from both branches, expression and PCA space)"
    )


def test_configured_methods_are_joint_only_and_unique():
    declared = declared_methods_by_distance(load_config())
    for distance, methods in declared.items():
        assert len(methods) == len(set(methods)), f"duplicate method ids in {distance}"
        for method_id in methods:
            assert method_id.endswith("_cell_tumour"), (
                f"{method_id}: only JOINT cell-line + tumour formulations are eligible"
            )
            assert method_id.startswith(("AGN_", "CCP_")), (
                f"{method_id}: formulations must come from the HC/k-means (AGN_*) "
                "or ConsensusClusterPlus (CCP_*) branch"
            )


def test_correlation_directions_have_no_kmeans():
    declared = declared_methods_by_distance(load_config())
    km_in_corr = [m for m in declared["corr"] if "_KM_" in m]
    assert km_in_corr == [], (
        "k-means minimises Euclidean inertia and is undefined under correlation "
        f"distance; found {km_in_corr} declared for corr"
    )


# ---------------------------------------------------------------------------
# Output-level tests: reconstruct p_consensus from method-specific tables
# ---------------------------------------------------------------------------

def _read_in_top_pairs(csv_path):
    """Return the set of (cell_tech_id, tumour_id) with in_top TRUE."""
    pairs = set()
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["in_top"].strip().upper() == "TRUE":
                pairs.add((row["cell_tech_id"], row["tumour_id"]))
    return pairs


def _reconstruct_and_compare(direction):
    cfg = load_config()
    declared = declared_methods_by_distance(cfg)
    distance = direction.rsplit("_", 1)[-1]
    declared_set = list(declared[distance])

    nh_dir = profile_unsup_root(cfg, PROFILE) / "tumour_neighbourhoods" / direction
    final_tsv = (
        nh_dir / "final_consensus"
        / f"Final_consensus_tumour_neighbourhoods_{direction}.tsv"
    )
    if not final_tsv.exists():
        return skip_or_fail(f"no generated outputs for {PROFILE}/{direction}: {final_tsv}")

    # Discovered method directories must equal the declared set exactly, with
    # exactly one correctly named neighbourhood table per formulation.
    discovered = sorted(
        p.name
        for p in nh_dir.iterdir()
        if p.is_dir() and p.name != "final_consensus"
        and list(p.glob("Top_m_long_*.csv"))
    )
    assert discovered == sorted(declared_set), (
        f"{direction}: discovered method dirs {discovered} != declared {sorted(declared_set)}"
    )
    for method_id in declared_set:
        csv_files = sorted((nh_dir / method_id).glob("Top_m_long_*.csv"))
        assert [p.name for p in csv_files] == [f"Top_m_long_{method_id}.csv"], (
            f"{direction}/{method_id}: expected exactly one Top_m_long_{method_id}.csv, "
            f"found {[p.name for p in csv_files]}"
        )

    n_methods_total = len(declared_set)
    expected_n = {"euc": 8, "corr": 4}[distance]
    assert n_methods_total == expected_n, (
        f"{direction}: configured |A_r| is {n_methods_total}, expected {expected_n}"
    )

    # Raw recurrence counts across the configured formulations.
    recurrence = Counter()
    for method_id in declared_set:
        method_pairs = _read_in_top_pairs(nh_dir / method_id / f"Top_m_long_{method_id}.csv")
        recurrence.update(method_pairs)

    # Reported p_consensus values.
    reported = {}
    with final_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["cell_tech_id"], row["tumour_id"])
            reported[key] = (int(row["n_methods"]), float(row["p_consensus"]))

    assert set(reported) == set(recurrence), (
        f"{direction}: pair universe mismatch between reconstruction "
        f"({len(recurrence)} pairs) and reported table ({len(reported)} pairs)"
    )

    for key, count in recurrence.items():
        rep_count, rep_p = reported[key]
        assert 1 <= count <= n_methods_total, (
            f"{direction} {key}: recurrence count {count} outside [1, {n_methods_total}] "
            "(a formulation was counted twice or the pair universe is corrupt)"
        )
        assert rep_count == count, (
            f"{direction} {key}: reported n_methods {rep_count} != reconstructed {count}"
        )
        assert abs(rep_p - count / n_methods_total) < 1e-9, (
            f"{direction} {key}: reported p_consensus {rep_p} != "
            f"{count}/{n_methods_total} = {count / n_methods_total}"
        )

    # Both branches must actually contribute somewhere in this direction.
    for branch in ("AGN_", "CCP_"):
        branch_methods = [m for m in declared_set if m.startswith(branch)]
        branch_pairs = set()
        for method_id in branch_methods:
            branch_pairs |= _read_in_top_pairs(
                nh_dir / method_id / f"Top_m_long_{method_id}.csv"
            )
        assert branch_pairs, f"{direction}: no neighbourhood members from branch {branch}*"

    print(
        f"OK {PROFILE}/{direction}: {len(recurrence)} pairs reconstructed, "
        f"|A_r| = {n_methods_total}, max recurrence = {max(recurrence.values())}"
    )
    return False


if pytest is not None:
    @pytest.mark.parametrize("direction", DIRECTIONS)
    def test_p_consensus_equals_recurrence_over_configured_methods(direction):
        _reconstruct_and_compare(direction)


def main():
    print(f"config: {CONFIG_PATH}")
    test_configured_method_sets_are_exact()
    test_configured_methods_are_joint_only_and_unique()
    test_correlation_directions_have_no_kmeans()
    print("OK config-level method-set checks")
    any_skipped = False
    for direction in DIRECTIONS:
        skipped = _reconstruct_and_compare(direction)
        any_skipped = any_skipped or bool(skipped)
    print("PASS (with skips)" if any_skipped else "PASS")


if __name__ == "__main__":
    sys.exit(main())
