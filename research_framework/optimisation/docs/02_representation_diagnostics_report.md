# Representation diagnostics report

## Redundancy diagnostics

The redundancy step compared representation-specific p-consensus value matrices within each cohort. Pairwise matrix agreement was estimated by Spearman correlation after missing pair values were set to zero. Nearest-neighbour overlap was estimated by Jaccard overlap of cell-line--tumour pairs with p-consensus values at or above the configured threshold.

## Assumptions

- Strong neighbourhood threshold: 0.7.
- Missing cell-line--tumour pairs were treated as observed p-consensus value 0 for matrix comparisons.
- Feature-set overlap was only recorded where direct feature membership tables were available; otherwise it remains unestimated.

## Redundancy classes

- complementary: 73 representation(s).
- partly_redundant: 3 representation(s).

## Outputs

- `research_framework/optimisation/results/representation_diagnostics/representation_redundancy.tsv`
- `research_framework/optimisation/results/representation_diagnostics/representation_distance_matrix_correlations.tsv`
- `research_framework/optimisation/results/representation_diagnostics/nearest_neighbour_overlap.tsv`

## Manual inspection required

- Inspect highly redundant representations before treating them as independent evidence sources.
- Inspect complementary representations with low matrix correlation and low nearest-neighbour overlap because they may encode distinct tumour-neighbourhood structure or noisy evidence.

## Stability diagnostics

The stability step estimated representation-level stability from observed p-consensus values, edge similarity when graph edges were available, and sensitivity to the configured strong-neighbourhood threshold.

## Assumptions

- Strong neighbourhood threshold: 0.7.
- Threshold sensitivity window: 0.05.
- Lightweight bootstrap enabled: FALSE.
- If no bootstrap artefacts were available and lightweight bootstrap was disabled, stability was estimated from existing thesis outputs only.

## Stability classes

- moderately_stable: 7 representation(s).
- stable: 69 representation(s).

## Output

- `research_framework/optimisation/results/representation_diagnostics/representation_stability.tsv`

## Manual inspection required

- Inspect representations with low ranking stability scores because they have many values close to the threshold.
- Inspect representations with missing edge stability scores because graph-edge evidence was unavailable or unreadable. 

## Representation reliability weights

Representation reliability weights were estimated as the product of stability, non-redundancy, batch-independence, and cancer-type preservation scores. Weights were normalised within cohort when configured.

## Highest-weighted representations

- `rbl::MeanAbsDev_corr`: weight 0.06939.
- `rbl::Spearman_corr`: weight 0.06925.
- `nbl::MX_corr`: weight 0.06685.
- `rbl::kTotal_euc`: weight 0.065.
- `rbl::kTotal_corr`: weight 0.06467.
- `nbl::MX_euc`: weight 0.06367.
- `rbl::Spearman_euc`: weight 0.06211.
- `nbl::Entropy_corr`: weight 0.06198.
- `rbl::MX_corr`: weight 0.06173.
- `brca::MAD_corr`: weight 0.06159.

## Lowest-weighted representations

- `brca::pam50_corr`: weight 0.02127.
- `brca::pam50_euc`: weight 0.02181.
- `brca::kTotal_corr`: weight 0.03577.
- `multicohort_cancer::Variance_euc`: weight 0.04142.
- `multicohort_cancer::MeanAbsDev_euc`: weight 0.04234.
- `nbl::kTotal_euc`: weight 0.04263.
- `nbl::kTotal_corr`: weight 0.04402.
- `multicohort_cancer::PCA_corr`: weight 0.04416.
- `multicohort_cancer::PCA_euc`: weight 0.04421.
- `rbl::MeanAbsDev_euc`: weight 0.04457.

## Output

- `research_framework/optimisation/results/representation_weights/representation_weights.tsv`

## Manual inspection required

- Inspect low-weight representations before excluding them from scientific interpretation.
- Inspect representations with neutral batch-independence scores because metadata associations were not estimated from available files. 
