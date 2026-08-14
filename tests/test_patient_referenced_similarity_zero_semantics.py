import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
HELPER_PATH = REPO_ROOT / "R" / "patient_referenced_similarity_utils.R"


def run_r(code: str) -> list[str]:
    result = subprocess.run(
        ["Rscript", "-e", f'source("{HELPER_PATH.as_posix()}");\n{code}'],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def test_threshold_includes_exact_boundary():
    lines = run_r(
        """
        mat <- matrix(c(0.70, 0.69, 1.00, 0.00), nrow = 2, byrow = TRUE)
        thr <- threshold_restrict_mean_profiles(mat, 0.70)
        cat(thr[1, 1], thr[1, 2], thr[2, 1], thr[2, 2], sep = "\\n")
        """
    )
    assert lines == ["0.7", "0", "1", "0"]


def test_joint_zero_invariance_for_pearson_and_jaccard():
    lines = run_r(
        """
        x1 <- c(1, 0)
        y1 <- c(0, 1)
        x2 <- c(1, 0, 0, 0, 0, 0)
        y2 <- c(0, 1, 0, 0, 0, 0)
        p1 <- compute_pairwise_active_union_similarity(x1, y1, metric = "pearson")
        p2 <- compute_pairwise_active_union_similarity(x2, y2, metric = "pearson")
        j1 <- compute_pairwise_active_union_similarity(x1, y1, metric = "jaccard")
        j2 <- compute_pairwise_active_union_similarity(x2, y2, metric = "jaccard")
        cat(p1$similarity, p2$similarity, j1$similarity, j2$similarity, sep = "\\n")
        """
    )
    assert lines == ["-1", "-1", "0", "0"]


def test_unilateral_zeros_remain_and_similarity_is_symmetric():
    lines = run_r(
        """
        x <- c(1.00, 0.75, 0, 0)
        y <- c(0, 0.75, 1.00, 0)
        p_xy <- compute_pairwise_active_union_similarity(x, y, metric = "pearson")
        p_yx <- compute_pairwise_active_union_similarity(y, x, metric = "pearson")
        j_xy <- compute_pairwise_active_union_similarity(x, y, metric = "jaccard")
        j_yx <- compute_pairwise_active_union_similarity(y, x, metric = "jaccard")
        cat(
          p_xy$n_pairwise_active_tumours,
          p_xy$n_shared_selected_tumours,
          p_xy$similarity == p_yx$similarity,
          j_xy$similarity == j_yx$similarity,
          sep = "\\n"
        )
        """
    )
    assert lines == ["3", "1", "TRUE", "TRUE"]


def test_empty_union_and_zero_variance_return_na_with_reason():
    lines = run_r(
        """
        j <- compute_pairwise_active_union_similarity(c(0, 0, 0), c(0, 0, 0), metric = "jaccard")
        p <- compute_pairwise_active_union_similarity(c(1, 1, 0), c(1, 1, 0), metric = "pearson")
        cat(is.na(j$similarity), j$undefined_similarity_reason, is.na(p$similarity), p$undefined_similarity_reason, sep = "\\n")
        """
    )
    assert lines == [
        "TRUE",
        "empty_pairwise_active_tumour_union",
        "TRUE",
        "zero_variance_active_profile",
    ]
