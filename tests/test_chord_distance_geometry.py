"""Regression tests for the Ward.D2 correlation geometry.

Ward.D2 minimises variance in a Euclidean space, so correlation-distance
hierarchical clustering must use a Euclidean-compatible dissimilarity. The
chord distance

    d(x, y) = sqrt(2 * (1 - r(x, y)))

is the exact Euclidean distance between the centred, unit-norm forms of x and
y. These tests exercise the real code path used by the HC/k-means clustering
branch (`correlation_chord_dist()` in R/hclust_kmeans_utils.R), not a
reimplementation.

Reference values:
    r =  1  ->  d = 0
    r =  0  ->  d = sqrt(2)
    r = -1  ->  d = 2

A separate guard asserts that the tumour/cell-line ranking dissimilarity in
R/base_functions/tumour_neighbourhood.R remains 1 - r: that path ranks
neighbours and is not fed to Ward, so it must not be converted to chord
distance.
"""

import math
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
UTILS = REPO_ROOT / "R" / "hclust_kmeans_utils.R"
NH_CORE = REPO_ROOT / "R" / "base_functions" / "tumour_neighbourhood.R"

TOL = 1e-12

# Rscript from the pipeline R environment when the caller supplies one.
RSCRIPT = os.environ.get("PIPELINE_RSCRIPT", "Rscript")


def run_r(code: str) -> list[str]:
    result = subprocess.run(
        [RSCRIPT, "-e", f'suppressPackageStartupMessages(source("{UTILS.as_posix()}"));\n{code}'],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def chord_distances_for_reference_matrix() -> dict[str, float]:
    """Return chord distances for a 3-sample matrix with exact r = 1, 0, -1.

    Rows are samples; correlations are computed between samples.
        s1 = (1, 2, 3)
        s2 = (2, 4, 6)   -> r(s1, s2) =  1
        s3 = (1, -2, 1)  -> r(s1, s3) =  0   (centred vectors are orthogonal)
        s4 = (3, 2, 1)   -> r(s1, s4) = -1
    """
    lines = run_r(
        """
        X <- rbind(
          s1 = c(1, 2, 3),
          s2 = c(2, 4, 6),
          s3 = c(1, -2, 1),
          s4 = c(3, 2, 1)
        )
        d <- as.matrix(correlation_chord_dist(X))
        cat(sprintf("%.17g", c(d["s1", "s2"], d["s1", "s3"],
                               d["s1", "s4"], d["s1", "s1"])), sep = "\\n")
        """
    )
    values = [float(x) for x in lines]
    return {
        "r_plus_1": values[0],
        "r_zero": values[1],
        "r_minus_1": values[2],
        "self": values[3],
    }


def test_chord_distance_reference_values():
    d = chord_distances_for_reference_matrix()
    assert abs(d["r_plus_1"] - 0.0) < TOL, f"r=1 must give chord distance 0, got {d['r_plus_1']}"
    assert abs(d["r_zero"] - math.sqrt(2)) < TOL, f"r=0 must give sqrt(2), got {d['r_zero']}"
    assert abs(d["r_minus_1"] - 2.0) < TOL, f"r=-1 must give 2, got {d['r_minus_1']}"
    assert abs(d["self"] - 0.0) < TOL, "self-distance must be 0"


def test_chord_distance_is_euclidean_realisable():
    """The chord distance must equal the Euclidean distance between the
    centred, unit-norm profiles — the property Ward.D2 requires."""
    lines = run_r(
        """
        set.seed(1)
        X <- matrix(rnorm(5 * 40), nrow = 5,
                    dimnames = list(paste0("s", 1:5), NULL))
        chord <- as.matrix(correlation_chord_dist(X))
        U <- t(apply(X, 1, function(v) { v <- v - mean(v); v / sqrt(sum(v^2)) }))
        euc <- as.matrix(dist(U, method = "euclidean"))
        cat(sprintf("%.17g", max(abs(chord - euc))), sep = "\\n")
        """
    )
    max_abs_diff = float(lines[0])
    assert max_abs_diff < 1e-10, (
        "chord distance must equal Euclidean distance between centred unit-norm "
        f"profiles; max |difference| = {max_abs_diff}"
    )


def test_zero_variance_profiles_do_not_produce_nan():
    """A constant profile has undefined correlation; it is treated as r = 0
    rather than propagating NaN into hclust."""
    lines = run_r(
        """
        X <- rbind(s1 = c(1, 2, 3), s2 = c(5, 5, 5), s3 = c(3, 2, 1))
        d <- correlation_chord_dist(X)
        cat(any(!is.finite(d)), sprintf("%.17g", max(d)), sep = "\\n")
        """
    )
    assert lines[0] == "FALSE", "chord distances must all be finite"
    assert abs(float(lines[1]) - 2.0) < TOL


def test_ranking_dissimilarity_remains_one_minus_r():
    """The tumour-neighbourhood ranking dissimilarity is not a Ward input and
    must remain 1 - r; converting it would change neighbourhood membership."""
    source = NH_CORE.read_text()
    assert "as.numeric(1 - cors)" in source, (
        "tumour_neighbourhood.R correlation_dist must remain 1 - r; it ranks "
        "neighbours and is not fed to Ward.D2"
    )
    assert "correlation_chord_dist" not in source, (
        "the chord transform must not leak into the neighbourhood ranking path"
    )


def test_chord_transform_is_applied_exactly_once_in_the_hc_path():
    """Guard against a double transformation: hc_optimal must delegate to the
    helper, and must not also inline a chord or 1 - r construction."""
    source = UTILS.read_text()
    hc_body = source.split("hc_optimal <- function", 1)[1]
    assert hc_body.count("correlation_chord_dist(X)") == 1, (
        "hc_optimal must build the correlation dissimilarity exactly once"
    )
    assert "as.dist(1 - cm)" not in hc_body, (
        "the non-Euclidean 1 - r dissimilarity must not be used with Ward.D2"
    )


if __name__ == "__main__":
    test_chord_distance_reference_values()
    d = chord_distances_for_reference_matrix()
    print(f"OK reference values: r=1 -> {d['r_plus_1']}, "
          f"r=0 -> {d['r_zero']:.15g} (sqrt(2)={math.sqrt(2):.15g}), "
          f"r=-1 -> {d['r_minus_1']}")
    test_chord_distance_is_euclidean_realisable()
    print("OK chord distance equals Euclidean distance on centred unit-norm profiles")
    test_zero_variance_profiles_do_not_produce_nan()
    print("OK zero-variance profiles yield finite distances")
    test_ranking_dissimilarity_remains_one_minus_r()
    print("OK neighbourhood ranking dissimilarity still 1 - r")
    test_chord_transform_is_applied_exactly_once_in_the_hc_path()
    print("OK chord transform applied exactly once in the HC path")
    print("PASS")
