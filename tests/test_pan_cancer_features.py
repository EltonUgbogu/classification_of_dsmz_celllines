import importlib.util
from pathlib import Path

import math
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_pan_cancer_features.py"
spec = importlib.util.spec_from_file_location("build_pan_cancer_features", MODULE_PATH)
pcf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pcf)

# The empirical acceptance quantiles are configuration-owned
# (marker_postprocessing.pan_cancer.empirical_quantile_thresholds) and reach the
# implementation as arguments. The module no longer carries copies of them, so
# these tests supply the same values configuration currently declares as an
# explicit fixture rather than importing constants back out of the module.
ADJUSTED_P_VALUE_QUANTILE = 0.25
ABSOLUTE_SHRUNKEN_LOG2FC_QUANTILE = 0.75
EXPRESSION_QUANTILE = 0.50


def evidence_row(
    gene_id,
    cancer_type="brca",
    stratum="anchor_associated",
    count=1,
    total=1,
    padj=0.01,
    lfc=2.0,
    base_mean=100.0,
    rank=1,
    direction="consistently_upregulated",
):
    marker_source_class = pcf.marker_source_class_from_evidence_stratum(stratum)
    return {
        "method": pcf.METHOD_NAME,
        "profile": cancer_type,
        "cancer_type": cancer_type,
        "cohort": cancer_type,
        "marker_evidence_stratum": stratum,
        "marker_source_class": marker_source_class,
        "gene_id": gene_id,
        "gene_id_raw": gene_id,
        "retained_marker_list_count": count,
        "retained_marker_list_total": total,
        "source_count": count,
        "marker_source_class_source_count": total,
        "source_contrast": f"{stratum}_{rank}",
        "source_file": f"{gene_id}.tsv",
        "source_marker_ranks": f"{stratum}_{rank}:{rank}",
        "minimum_adjusted_p_value": padj,
        "median_adjusted_p_value": padj,
        "minimum_p_value": padj / 2,
        "median_absolute_shrunken_log2fc": lfc,
        "maximum_absolute_shrunken_log2fc": lfc,
        "median_shrunken_log2fc": lfc,
        "mean_baseMean": base_mean,
        "median_base_mean": base_mean,
        "max_baseMean": base_mean,
        "minimum_contrast_marker_rank": rank,
        "n_positive_effects": 1 if direction == "consistently_upregulated" else 0,
        "n_negative_effects": 1 if direction == "consistently_downregulated" else 0,
        "n_zero_effects": 0,
        "direction_consistency_class": direction,
    }


def classify_one(count, total):
    df = pd.DataFrame([evidence_row("ENSG1", count=count, total=total)])
    return pcf.recurrence_classification(df).loc[0, "candidate_pool_type"]


def test_recurrence_classification_invariants():
    assert classify_one(count=1, total=1) == "singleton_candidate"
    assert classify_one(count=1, total=3) == "non_recurrent_candidate"
    assert classify_one(count=2, total=3) == "recurrent"
    assert classify_one(count=4, total=5) == "recurrent"


def test_empirical_threshold_boundaries_and_candidate_acceptance():
    pool = pcf.candidate_pool_construction(
        pcf.recurrence_classification(
            pd.DataFrame([evidence_row("ENSG1", padj=0.03, lfc=2.0, base_mean=20.0)])
        )
    )
    pool, thresholds = pcf.empirical_quantile_threshold_calculation(
        pool,
        ADJUSTED_P_VALUE_QUANTILE,
        ABSOLUTE_SHRUNKEN_LOG2FC_QUANTILE,
        EXPRESSION_QUANTILE,
    )
    pool = pcf.evaluate_candidate_acceptance(pool)
    threshold = thresholds.iloc[0]
    boundary_row = pool.loc[pool["gene_id"] == "ENSG1"].iloc[0]
    assert math.isclose(boundary_row["minimum_adjusted_p_value"], threshold["adjusted_p_value_quantile_threshold"])
    assert math.isclose(boundary_row["median_absolute_shrunken_log2fc"], threshold["absolute_shrunken_log2fc_quantile_threshold"])
    assert math.isclose(boundary_row["median_base_mean"], threshold["expression_quantile_threshold"])
    assert bool(boundary_row["passes_candidate_acceptance"])

    rejected_pool = pool.copy()
    rejected_pool.loc[:, "gene_id"] = "ENSG2"
    rejected_pool.loc[:, "median_absolute_shrunken_log2fc"] = 1.0
    rejected_pool.loc[:, "absolute_shrunken_log2fc_quantile_threshold"] = 2.0
    rejected_pool = pcf.evaluate_candidate_acceptance(rejected_pool)
    rejected = rejected_pool.iloc[0]
    assert bool(rejected["passes_statistical_evidence"])
    assert not bool(rejected["passes_effect_magnitude"])
    assert bool(rejected["passes_expression_evidence"])
    assert not bool(rejected["passes_candidate_acceptance"])


def test_recurrent_row_is_retained_without_candidate_flags():
    recurrent = pd.DataFrame([evidence_row("ENSG1", count=2, total=3, padj=1.0, lfc=0.1, base_mean=1.0)])
    pool = pcf.candidate_pool_construction(pcf.recurrence_classification(recurrent))
    pool, thresholds = pcf.empirical_quantile_threshold_calculation(
        pool,
        ADJUSTED_P_VALUE_QUANTILE,
        ABSOLUTE_SHRUNKEN_LOG2FC_QUANTILE,
        EXPRESSION_QUANTILE,
    )
    assert thresholds.empty
    pool = pcf.effect_direction_consistency_classification(pool)
    pool = pcf.evaluate_candidate_acceptance(pool)
    selected = pcf.select_pan_cancer_evidence_rows(pool)
    assert selected["gene_id"].tolist() == ["ENSG1"]
    assert selected.loc[0, "selection_basis"] == "recurrent"
    assert selected[[
        "passes_statistical_evidence",
        "passes_effect_magnitude",
        "passes_expression_evidence",
        "passes_candidate_acceptance",
    ]].isna().all(axis=None)


def test_cross_context_dedup_preserves_provenance_and_priority():
    rows = [
        evidence_row("ENSG1", cancer_type="brca", stratum="anchor_associated", count=2, total=3, rank=1),
        evidence_row("ENSG1", cancer_type="nbl", stratum="isolate_associated", count=1, total=1, rank=2),
        evidence_row("ENSG2", cancer_type="rbl", stratum="anchor_associated", count=1, total=1, rank=3),
    ]
    pool = pcf.candidate_pool_construction(pcf.recurrence_classification(pd.DataFrame(rows)))
    pool, _ = pcf.empirical_quantile_threshold_calculation(
        pool,
        ADJUSTED_P_VALUE_QUANTILE,
        ABSOLUTE_SHRUNKEN_LOG2FC_QUANTILE,
        EXPRESSION_QUANTILE,
    )
    pool = pcf.effect_direction_consistency_classification(pool)
    pool = pcf.evaluate_candidate_acceptance(pool)
    selected = pcf.select_pan_cancer_evidence_rows(pool)
    gene_evidence = pcf.build_gene_evidence_table(selected)
    features, duplicate_count, removed = pcf.collapse_gene_evidence_across_cancer_types_and_strata(
        ["brca", "nbl", "rbl"],
        gene_evidence,
        annotation={},
        remove_ribo_mt=False,
    )
    assert removed == set()
    assert duplicate_count == 1
    assert features["gene_id"].tolist().count("ENSG1") == 1
    row = features.loc[features["gene_id"] == "ENSG1"].iloc[0]
    assert row["feature_class"] == "recurrent"
    assert set(row["selection_basis_classes"].split(";")) == {"recurrent", "accepted_singleton"}
    assert row["has_anchor_associated_marker_evidence"]
    assert row["has_isolate_associated_marker_evidence"]


def test_non_empty_candidate_context_has_single_finite_threshold_row():
    pool = pcf.candidate_pool_construction(
        pcf.recurrence_classification(pd.DataFrame([evidence_row("ENSG1", count=1, total=2)]))
    )
    _, thresholds = pcf.empirical_quantile_threshold_calculation(
        pool,
        ADJUSTED_P_VALUE_QUANTILE,
        ABSOLUTE_SHRUNKEN_LOG2FC_QUANTILE,
        EXPRESSION_QUANTILE,
    )
    assert len(thresholds) == 1
    assert thresholds[[
        "adjusted_p_value_quantile_threshold",
        "absolute_shrunken_log2fc_quantile_threshold",
        "expression_quantile_threshold",
    ]].map(math.isfinite).all(axis=None)
