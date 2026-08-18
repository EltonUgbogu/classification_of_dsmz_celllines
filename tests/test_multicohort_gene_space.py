"""Regression tests for the multicohort common-gene-space construction.

The multicohort profile concatenates per-cohort matrices whose gene sets can
differ. These tests exercise R/multicohort_gene_space.R directly and assert
that it:

  - intersects gene sets on identifiers (never on row/column position);
  - returns an identical gene order for every cohort;
  - drops non-shared genes rather than zero-filling them;
  - rejects duplicate or blank gene identifiers;
  - refuses an implausibly small intersection;
  - detects cached feature lists that are no longer subsets of G_common.

The shuffled-order case is the important one: a builder that aligned genes by
position would silently mismatch expression values across cohorts while still
producing a matrix of the right shape.
"""

import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
HELPER = REPO_ROOT / "R" / "multicohort_gene_space.R"

RSCRIPT = os.environ.get("PIPELINE_RSCRIPT", "Rscript")


def run_r(code: str) -> list[str]:
    result = subprocess.run(
        [RSCRIPT, "-e", f'suppressPackageStartupMessages(source("{HELPER.as_posix()}"));\n{code}'],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def run_r_expect_error(code: str) -> str:
    result = subprocess.run(
        [RSCRIPT, "-e", f'suppressPackageStartupMessages(source("{HELPER.as_posix()}"));\n{code}'],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0, f"expected an error, got success:\n{result.stdout}"
    return result.stderr


# Three cohorts. Gene sets deliberately differ and are deliberately ordered
# differently, and the values encode (cohort, gene) so that any positional
# misalignment changes the recovered numbers.
FIXTURE = """
    mk <- function(genes, tag) {
      m <- matrix(
        as.numeric(outer(c(1, 2), seq_along(genes), function(s, g) tag * 1000 + g)),
        nrow = 2, byrow = FALSE,
        dimnames = list(paste0(tag, "_s", 1:2), genes)
      )
      # encode the gene identity itself so misalignment is detectable
      for (j in seq_along(genes)) m[, j] <- as.numeric(sub("^G", "", genes[j]))
      m
    }
    A <- mk(c("G1","G2","G3","G4","G5"), 1)          # superset
    B <- mk(c("G5","G3","G2","G1"),      2)          # shared four, reversed order
    C <- mk(c("G2","G1","G3","G5","G9"), 3)          # shared four plus a private gene
    X <- list(brca = A, nbl = B, rbl = C)
"""


def test_intersection_is_by_identifier_and_order_is_identical():
    lines = run_r(
        FIXTURE
        + """
        res <- harmonise_gene_space(X, min_common_genes = 1L)
        cat(paste(res$genes_common, collapse = ","), sep = "\\n")
        cat(paste(colnames(res$X_list$brca), collapse = ","), sep = "\\n")
        cat(paste(colnames(res$X_list$nbl),  collapse = ","), sep = "\\n")
        cat(paste(colnames(res$X_list$rbl),  collapse = ","), sep = "\\n")
        """
    )
    common, a, b, c = lines[0], lines[1], lines[2], lines[3]
    # G4 is brca-only, G9 is rbl-only; both must be dropped.
    assert common == "G1,G2,G3,G5", f"unexpected common gene space: {common}"
    assert a == b == c == common, (
        f"gene order must be identical across cohorts: brca={a}, nbl={b}, rbl={c}"
    )


def test_values_follow_gene_identity_not_position():
    """Each cell holds the gene's numeric id, so identifier-based alignment
    yields the same value in a column across cohorts. Positional alignment
    would not, because the cohorts were built with different gene orders."""
    lines = run_r(
        FIXTURE
        + """
        res <- harmonise_gene_space(X, min_common_genes = 1L)
        vals <- sapply(res$X_list, function(M) M[1, ])
        cat(paste(rownames(vals), collapse = ","), sep = "\\n")
        cat(paste(apply(vals, 1, function(r) length(unique(r))), collapse = ","), sep = "\\n")
        cat(paste(vals[, "brca"], collapse = ","), sep = "\\n")
        """
    )
    genes, uniques, brca_vals = lines[0], lines[1], lines[2]
    assert genes == "G1,G2,G3,G5"
    assert uniques == "1,1,1,1", (
        "each gene must carry one value across cohorts; differing values mean "
        "genes were aligned by position rather than identifier"
    )
    assert brca_vals == "1,2,3,5", f"values must match gene identity, got {brca_vals}"


def test_non_shared_genes_are_dropped_not_zero_filled():
    lines = run_r(
        FIXTURE
        + """
        res <- harmonise_gene_space(X, min_common_genes = 1L)
        cat(as.integer(any(c("G4","G9") %in% res$genes_common)), sep = "\\n")
        cat(as.integer(any(sapply(res$X_list, function(M) any(M == 0)))), sep = "\\n")
        cat(paste(res$audit$n_genes_before, collapse = ","), sep = "\\n")
        cat(paste(res$audit$n_genes_after,  collapse = ","), sep = "\\n")
        cat(paste(res$audit$n_genes_removed, collapse = ","), sep = "\\n")
        """
    )
    assert lines[0] == "0", "cohort-private genes must not appear in G_common"
    assert lines[1] == "0", "dropped genes must not be zero-filled"
    assert lines[2] == "5,4,5", f"n_genes_before wrong: {lines[2]}"
    assert lines[3] == "4,4,4", f"n_genes_after wrong: {lines[3]}"
    assert lines[4] == "1,0,1", f"n_genes_removed wrong: {lines[4]}"


def test_duplicate_gene_identifiers_are_rejected():
    err = run_r_expect_error(
        """
        A <- matrix(1:6, nrow = 2, dimnames = list(c("s1","s2"), c("G1","G2","G1")))
        B <- matrix(1:4, nrow = 2, dimnames = list(c("t1","t2"), c("G1","G2")))
        harmonise_gene_space(list(brca = A, nbl = B), min_common_genes = 1L)
        """
    )
    assert "duplicate gene identifier" in err, err


def test_blank_identifiers_are_rejected():
    err = run_r_expect_error(
        """
        A <- matrix(1:4, nrow = 2, dimnames = list(c("s1","s2"), c("G1","")))
        B <- matrix(1:4, nrow = 2, dimnames = list(c("t1","t2"), c("G1","G2")))
        harmonise_gene_space(list(brca = A, nbl = B), min_common_genes = 1L)
        """
    )
    assert "missing or blank gene identifiers" in err, err


def test_small_intersection_fails_explicitly():
    err = run_r_expect_error(
        FIXTURE + "harmonise_gene_space(X, min_common_genes = 50L)"
    )
    assert "are shared by all cohorts" in err, err


def test_feature_list_validation_detects_stale_entries():
    lines = run_r(
        """
        common <- c("G1","G2","G3")
        ok    <- validate_feature_list(c("G1","G3"), common, "ok_list")
        stale <- validate_feature_list(c("G1","G7","G8"), common, "stale_list")
        cat(as.integer(ok$ok), ok$n_missing, sep = "\\n")
        cat(as.integer(stale$ok), stale$n_missing,
            paste(stale$missing, collapse = ","), sep = "\\n")
        """
    )
    assert lines[0] == "1" and lines[1] == "0"
    assert lines[2] == "0" and lines[3] == "2" and lines[4] == "G7,G8"


if __name__ == "__main__":
    test_intersection_is_by_identifier_and_order_is_identical()
    print("OK intersection by identifier, identical order across cohorts")
    test_values_follow_gene_identity_not_position()
    print("OK values follow gene identity, not row position")
    test_non_shared_genes_are_dropped_not_zero_filled()
    print("OK non-shared genes dropped, never zero-filled; audit counts correct")
    test_duplicate_gene_identifiers_are_rejected()
    print("OK duplicate gene identifiers rejected")
    test_blank_identifiers_are_rejected()
    print("OK blank gene identifiers rejected")
    test_small_intersection_fails_explicitly()
    print("OK implausibly small intersection fails explicitly")
    test_feature_list_validation_detects_stale_entries()
    print("OK stale feature-list entries detected")
    print("PASS")
