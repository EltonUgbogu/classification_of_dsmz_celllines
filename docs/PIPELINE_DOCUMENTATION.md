# Unsupervised Clustering and Patient-Referenced Tumour Neighbourhood Pipeline

## Overview

This document describes the computational pipeline used to quantify transcriptomic
relationships between patient tumour samples and cancer cell lines. The workflow
combines multiple feature-selection strategies, Euclidean and Pearson-correlation
dissimilarities, clustering formulations, adaptive tumour neighbourhoods,
tumour-wise `p_consensus` recurrence, patient-referenced cell-line graphs,
multi-representation graph resolution, graph-derived marker contrasts, pan-cancer
marker aggregation, and functional enrichment analysis.

The patient-referenced graph stage evaluates two similarity definitions in
parallel:

1. Pearson correlation of threshold-restricted continuous `p_consensus` fractions
   over the pairwise active-tumour union.
2. Jaccard similarity derived from the same threshold-restricted continuous
   `p_consensus` fractions after conversion to binary tumour membership.

Both definitions use the production `p_consensus` threshold of 0.70 and a
representation-specific similarity quantile of 0.90. Tumours unselected by both
cell lines do not contribute to either pairwise similarity calculation.

> **Implementation status (2026-08-14):** the tumour-neighbourhood and
> `p_consensus` path consumes BOTH eligible clustering formulations — the
> HC/k-means clustering branch (`AGN_*`) and ConsensusClusterPlus (`CCP_*`) —
> using JOINT cell-line–tumour outputs only. The contributing formulation set
> is declared in `patient_referenced_graph.clustering_methods_by_distance`
> and validated exactly at run time: Euclidean representations use
> `n_methods = 8`, correlation representations `n_methods = 4`. The
> `p_consensus` denominator is the configured formulation count, never
> inferred from directories present on disk; missing or unexpected
> formulations are explicit errors.

## Table of Contents

1. Scientific Background and Objectives
2. Pipeline Architecture
3. Data Sources and Inputs
4. Feature Selection
5. Dissimilarity Measures and Clustering
6. Adaptive Tumour Neighbourhoods and `p_consensus`
7. Patient-Referenced Cell-Line Graphs
8. Multi-Representation Graph Resolution
9. Marker Analysis
10. Multicohort and Pan-Cancer Analysis
11. Key Output Files
12. Configuration
13. Reproducibility and Validation
14. Glossary

---

## 1. Scientific Background and Objectives

### 1.1 Patient-Referenced Cell-Line Analysis

Cancer cell lines are widely used experimental models. Prolonged culture and
adaptation can produce transcriptomic states that differ from those observed in
patient tumours. The pipeline therefore evaluates cell lines relative to patient
tumour expression profiles rather than assuming that a labelled cancer type alone
establishes transcriptomic proximity.

### 1.2 Main Objectives

The pipeline asks:

1. Which patient tumours are transcriptomically close to each cell line under
   different feature-distance representations?
2. Which cell-line–tumour neighbourhood relationships recur across clustering
   formulations within a fixed representation?
3. Which cell-line relationships recur across multiple feature-distance
   representations?
4. Which cell lines become isolates or central nodes after graph resolution?
5. How sensitive is the patient-referenced graph to the definition of
   cell-line similarity?
6. Which genes distinguish graph-defined cell-line roles, and which biological
   programmes are represented by the resulting pan-cancer marker panel?

### 1.3 Tumour-Neighbourhood Concept

For a given feature-distance representation and clustering formulation, each cell
line is associated with an adaptive neighbourhood of patient tumours. Candidate
tumours are identified from the relevant joint clustering solution and ranked by
their dissimilarity to the cell line.

The neighbourhood is therefore a method-defined set of nearby patient tumours. Its
size is governed by the adaptive neighbourhood rule and should not be interpreted
as a stand-alone measure of how representative a cell line is.

Recurrence of a specific cell-line–tumour relationship across clustering
formulations is quantified by the `p_consensus` fraction.

---

## 2. Pipeline Architecture

The analysis is organised into the following layers.

### Layer 1: Data Preprocessing and Feature Selection

- Batch correction and variance-stabilising transformation (VST)
- Purity filtering of patient tumour profiles
- Six feature-ranking methods, each selecting a top-3,000 gene set
- Three connectivity-based methods, each selecting a top-500 gene set

### Layer 2: Feature-Distance Representations

- Nine feature-selection methods
- Euclidean distance
- Pearson correlation distance
- Eighteen cancer-specific feature-distance representations

### Layer 3: Clustering

- HC/k-means clustering in expression and PCA-reduced spaces
- ConsensusClusterPlus clustering in expression and PCA-reduced spaces
- Cell-only, tumour-only, and joint cell-line–tumour sample scopes where implemented
- Only joint cell-line–tumour clustering solutions are eligible to define
  cell-line tumour neighbourhoods

### Layer 4: Adaptive Tumour-Neighbourhood Computation

- Candidate tumours are obtained from the relevant joint clustering solution
- Candidate tumours are ranked by tumour–cell-line dissimilarity
- Adaptive top-k tumour neighbourhoods are constructed
- `p_consensus` is calculated tumour-wise within each feature-distance representation

The two clustering formulations of Layer 3 are parallel branches: neither
consumes the other's outputs, and their method-specific cell-line tumour
neighbourhoods are brought together at the `p_consensus` stage.

```text
                Feature-distance representation
                           ↓
                  clustering formulations
                    ↙             ↘
             HC / k-means      Consensus
              clustering       clustering
                    ↘             ↙
             method-specific cell-line
                tumour neighbourhoods
                           ↓
                 p_consensus fraction
```

### Layer 5: Representation-Specific Patient-Referenced Graphs

- Profile-level `p_consensus` fractions are mean-pooled for biological replicates
- A threshold of `p_consensus >= 0.70` defines the threshold-restricted continuous
  tumour-association profile
- Pearson and Jaccard cell-line similarity are calculated independently
- A representation- and metric-specific 90th-percentile threshold selects graph edges

### Layer 6: Multi-Representation Graph Resolution

- Edge recurrence is evaluated across active feature-distance representations
- A majority recurrence threshold is derived from the number of active
  representations
- Graph-resolved neighbour sets are obtained through the global–local neighbour
  intersection

### Layer 7: Marker Prioritisation and Functional Analysis

- Graph-defined isolate and anchor contrasts
- Differential expression analysis with DESeq2
- Marker selection under declared statistical and expression criteria
- Pan-cancer marker aggregation
- Functional enrichment analysis

### Layer 8: Multicohort and Pan-Cancer Analysis

- Cancer-type annotation of multicohort UMAP and graph visualisations
- Pan-cancer marker-derived analyses
- Qualitative examination of multicohort transcriptomic organisation

### 2.1 Feature-Distance Representation Strategy

A central feature of this pipeline is the systematic evaluation of multiple
feature-distance representations, defined by combinations of feature-selection
methods and dissimilarity measures. This allows transcriptomic similarity to be
evaluated across several representations of tumour–cell-line relationships.

The standard cancer-specific grid contains:

- 9 feature-selection methods
- 2 dissimilarity measures

This yields 18 feature-distance representations per cancer type, subject to
profile-specific configuration.

PAM50 is breast-cancer-specific and is treated separately from the generic
cancer-specific feature-distance grid.

---

## 3. Data Sources and Inputs

### 3.1 Patient Tumour Cohorts

Current patient tumour profiles include:

- **Breast cancer:** TCGA-BRCA
- **Neuroblastoma:** TARGET and additional public neuroblastoma cohorts
- **Retinoblastoma:** public retinoblastoma cohorts

The current purity-filtered joint VST matrices contain:

- **Breast cancer:** 773 patient tumours + 29 cell-line profiles
- **Neuroblastoma:** 235 patient tumours + 18 cell-line profiles
- **Retinoblastoma:** 68 patient tumours + 9 cell-line profiles

For retinoblastoma, profile-level replicates are combined to 9 biological cell
lines before patient-referenced graph construction.

### 3.2 Cell-Line Inputs

The analysed cell lines are obtained from DSMZ. The pipeline uses cell-line
expression profiles and metadata alongside patient tumour data.

Expression preprocessing is designed to place tumour and cell-line profiles in a
common VST expression coordinate system for the patient-referenced analyses.

### 3.3 Required Inputs

Profile-specific inputs include:

1. Joint tumour–cell-line VST expression matrix
2. Tumour-only VST expression matrix where required
3. Cell-line-only VST expression matrix where required
4. Cell-line metadata
5. Study-design and source metadata
6. Direction-specific gene lists
7. Cell-line-only raw integer count matrices for DESeq2 marker analysis
8. PAM50 gene list where breast-cancer-specific annotation is required

Clustering and graph construction use VST expression data. DESeq2 marker analysis
uses raw integer counts rather than VST values.

---

## 4. Feature Selection

### 4.1 Six Top-3,000 Feature-Ranking Methods

Six methods each select a top-3,000 gene set:

- **Variance:** ranks genes by variance across samples
- **HVG residual variance:** ranks genes by residual variation after accounting
  for the mean–variance relationship
- **Median absolute deviation (MAD):** ranks genes by median-based dispersion
- **Mean absolute deviation (MeanAbsDev):** ranks genes by mean absolute dispersion
- **Shannon entropy:** ranks genes by diversity of expression values across samples
- **PCA loadings:** ranks genes by their contribution to leading principal components

These methods do not all estimate the same statistical property. They provide
alternative feature rankings based on dispersion, expression-value diversity, or
contribution to multivariate variation.

For entropy, a high score indicates a broader distribution of observations across
the discretised expression states. Entropy does not specifically identify
bimodality or multimodality.

### 4.2 Three Top-500 Connectivity-Based Methods

The union of the six top-3,000 gene lists defines the candidate feature set for
three connectivity-based rankings:

- **Spearman connectivity:** mean absolute Spearman correlation
- **MX score:** Spearman connectivity weighted by scaled variance
- **WGCNA kTotal:** total unsigned, soft-thresholded connectivity

Each method selects its top 500 genes.

Conceptually:

```text
Spearman connectivity:
K_S(g) = mean_h |rho_S(g,h)|

MX score:
MX(g) = K_S(g) * scaled_variance(g)

WGCNA kTotal:
a_gh = |cor(g,h)|^beta
kTotal(g) = sum_h a_gh
```

### 4.3 Feature-Distance Representations

Each selected gene set defines a sample representation. Crossing the nine
feature-selection methods with the two dissimilarity measures gives:

```text
9 feature-selection methods × 2 dissimilarity measures = 18 representations
```

The two dissimilarity measures are:

- Euclidean distance
- Pearson correlation distance

---

## 5. Dissimilarity Measures and Clustering

### 5.1 Euclidean Distance

For two samples \(x\) and \(y\) in a selected \(p\)-gene expression space:

\[
d_E(x,y)
=
\sqrt{\sum_{g=1}^{p}(x_g-y_g)^2}
\]

Euclidean distance measures absolute separation in the selected VST expression
coordinates.

### 5.2 Pearson Correlation Distance

\[
d_{\mathrm{corr}}(x,y)
=
1-
ho(x,y)
\]

where \(
ho(x,y)\) is the Pearson correlation coefficient across the selected
genes.

Pearson correlation distance measures disagreement in relative expression
patterns. It is invariant to additive shifts and positive rescaling of an entire
sample profile.

`1 - Pearson correlation` is a dissimilarity measure and is not generally a
mathematical metric.

### 5.3 HC/K-Means Clustering

The workflow contains a single-run HC/k-means branch applied separately to
configured feature-distance representations.

It includes:

- hierarchical clustering in expression space
- hierarchical clustering in PCA-reduced space
- k-means in expression space for Euclidean representations
- k-means in PCA-reduced space for Euclidean representations

The branch is evaluated over:

- cell-line-only profiles
- tumour-only profiles
- joint cell-line–tumour profiles

Only joint cell-line–tumour clustering solutions can directly define
cell-line tumour neighbourhoods.

The hierarchical-clustering \(k\) grid is evaluated from \(k=2\) to \(k=8\)
where applicable. K-means is restricted to Euclidean representations.

Hierarchical clustering uses Ward.D2. Ward.D2 is Euclidean in its standard
formulation, so correlation-distance HC uses an explicitly Euclidean-compatible
construction: the chord distance `sqrt(2 * (1 - r))`, which is the exact
Euclidean distance between centred, unit-norm sample vectors. This matches the
correlation-geometry transform (centre + unit norm + Euclidean distance) applied
by the ConsensusClusterPlus branch, so both branches operate on the same
correlation geometry and both are eligible tumour-neighbourhood formulations.

### 5.4 Expression-Space and PCA-Space Clustering

Two input spaces are evaluated:

- **Expression space:** clustering is applied to the selected-gene expression
  coordinates
- **PCA-reduced space:** clustering is applied after projection onto the
  configured leading principal components

PCA provides a lower-dimensional coordinate system derived from the selected
expression matrix. It should not be interpreted as automatically removing
biological or technical noise.

### 5.5 ConsensusClusterPlus Clustering

ConsensusClusterPlus provides a separate resampling-based clustering stage. Its
joint cell-line–tumour formulations contribute to tumour-neighbourhood
construction alongside the HC/k-means formulations:

For correlation representations:

- `CCP_HC_expr_cell_tumour`
- `CCP_HC_pca_cell_tumour`

For Euclidean representations:

- `CCP_HC_expr_cell_tumour`
- `CCP_HC_pca_cell_tumour`
- `CCP_KM_expr_cell_tumour`
- `CCP_KM_pca_cell_tumour`

ConsensusClusterPlus does not consume the HC/k-means branch outputs. The two
branches operate as separate clustering formulations over the relevant
representation.

## 6. Adaptive Tumour Neighbourhoods and `p_consensus`

### 6.1 Candidate Tumour Set

For cell line \(i\), feature-distance representation \(r\), and clustering
formulation \(a\), candidate tumours are obtained from the relevant joint
cell-line–tumour clustering solution.

The candidate tumours are then ranked by tumour–cell-line dissimilarity within
representation \(r\).

### 6.2 Adaptive Neighbourhood Size

If \(n_i\) is the number of candidate tumours for cell line \(i\), the adaptive
neighbourhood size is:

\[
k_i
=
\min\left(
200,\;
\max\left(30,\left\lceil0.10n_i
ight
ceil
ight),\;
n_i

ight)
\]

Thus the neighbourhood:

- scales with 10% of the candidate set
- has a target minimum of 30 tumours when at least 30 candidates are available
- is capped at 200 tumours
- never exceeds the available candidate count

The method-specific tumour neighbourhood is denoted
\(\mathcal{N}_{i,r,a}\).

### 6.3 Tumour-Wise `p_consensus` Fraction

For cell line \(i\), tumour \(t\), representation \(r\), and configured
clustering-method set \(\mathcal{A}_r\):

\[
p_{\mathrm{consensus}}(i,t\mid r)
=

rac{1}{|\mathcal{A}_r|}
\sum_{a\in\mathcal{A}_r}
\mathbf{1}
\left[
t\in\mathcal{N}_{i,r,a}

ight]
\]

The numerator is the number of eligible clustering formulations in which tumour
\(t\) belongs to the adaptive tumour neighbourhood of cell line \(i\).

`p_consensus` is a recurrence fraction across the configured clustering
formulations.

### 6.5 Biological Replicate Pooling and Thresholding

For biological cell line \(s\), let \(U_s\) be its set of profile-level
measurements. The pooled tumour-wise fraction is:

\[
ar p_{s,t,r}
=

rac{1}{|U_s|}
\sum_{u\in U_s}
p_{\mathrm{consensus}}(u,t\mid r)
\]

The patient-referenced graph uses the production threshold:

\[
	heta=0.70
\]

and defines the threshold-restricted continuous profile:

\[
z_{s,t,r}
=
egin{cases}
ar p_{s,t,r}, & ar p_{s,t,r}\ge	heta\
0, & 	ext{otherwise}
\end{cases}
\]

Both Pearson and Jaccard graph definitions begin from this same
threshold-restricted continuous profile.

---

## 7. Patient-Referenced Cell-Line Graphs

### 7.1 Pairwise Active-Tumour Union

For biological cell lines \(i\) and \(j\) in representation \(r\), define:

\[
A_{ij,r}
=
\left\{
t:
z_{i,t,r}>0
\;\lor\;
z_{j,t,r}>0

ight\}
\]

This is the pairwise active-tumour union.

For each tumour:

- positive for both cell lines: included
- positive for only one cell line: included as disagreement
- zero for both cell lines: excluded

Joint zeros therefore do not contribute evidence of cell-line similarity.

### 7.2 Pearson Similarity

Pearson uses the threshold-restricted continuous `p_consensus` fractions over
the pairwise active-tumour union:

\[
S^{(P)}_r(i,j)
=
\operatorname{cor}_{t\in A_{ij,r}}
\left(
z_{i,t,r},
z_{j,t,r}

ight)
\]

It quantifies linear concordance between the two threshold-restricted tumour
association profiles among tumours selected by at least one of the two cell
lines.

Pearson similarity is undefined when the restricted vectors do not permit a
valid correlation, including zero-variance cases. Undefined similarities are
recorded and excluded from graph-threshold estimation.

### 7.3 Jaccard Similarity

Jaccard starts from the same threshold-restricted continuous profile and converts
it to binary tumour membership:

\[
b_{i,t,r}
=
\mathbf{1}[z_{i,t,r}>0]
\]

Then:

\[
J_r(i,j)
=

rac{
\sum_t b_{i,t,r}b_{j,t,r}
}{
\sum_t
\mathbf{1}
[b_{i,t,r}=1\lor b_{j,t,r}=1]
}
\]

Equivalently:

\[
J_r(i,j)
=

rac{|N_{i,r}\cap N_{j,r}|}
{|N_{i,r}\cup N_{j,r}|}
\]

where \(N_{i,r}\) is the set of tumours selected for cell line \(i\) after the
0.70 threshold.

Jaccard therefore measures agreement in binary tumour membership. Qualifying
continuous values such as 0.75 and 1.00 both become binary membership value 1.

If neither cell line selects any tumour, the union is empty and Jaccard is
undefined (`NA`).

### 7.4 Representation-Specific Graph Threshold

For metric \(M\in\{\mathrm{Pearson},\mathrm{Jaccard}\}\), collect all defined
non-self unordered cell-line-pair similarities in representation \(r\):

\[
\mathcal{S}_{r,M}
=
\left\{
S_{r,M}(i,j):
i<j,\;
S_{r,M}(i,j)	ext{ defined}

ight\}
\]

The production edge threshold is:

\[
	au_{r,M}
=
Q_{0.90}(\mathcal{S}_{r,M})
\]

An undirected edge is selected when:

\[
\{i,j\}\in E_{r,M}
\iff
S_{r,M}(i,j)\ge	au_{r,M}
\]

The 0.90 quantile is computed separately for every metric and representation.
The Pearson threshold is not reused for Jaccard.

The `>=` edge rule includes all ties at the percentile threshold. The selected
edge fraction can therefore exceed 10%, especially for discrete Jaccard
similarities.

### 7.5 Representation-Specific Graph

For each representation and similarity metric:

\[
G_{r,M}
=
(S,E_{r,M})
\]

where \(S\) is the set of biological cell lines and \(E_{r,M}\) is the
threshold-selected edge set.

The graph is a patient-referenced cell-line similarity graph. Its edges encode
similarity derived from tumour-neighbourhood recurrence profiles and should not
be interpreted as direct biological interactions.

### 7.6 Community Structure

Community detection can be applied to graph outputs to describe groups of nodes
with relatively dense internal connectivity.

Community assignments describe graph structure under the selected graph
definition. They are not assumed to correspond to predefined biological classes,
common driver mutations, or drug-response classes unless independent evidence is
provided.

---

## 8. Multi-Representation Graph Resolution

### 8.1 Edge Recurrence Across Representations

For a fixed similarity metric \(M\) and unordered cell-line edge \(e\):

\[
a_M(e)
=
\sum_{r\in\mathcal{R}}
\mathbf{1}[e\in E_{r,M}]
\]

The edge-recurrence fraction is:

\[
f_M(e)
=

rac{a_M(e)}{|\mathcal{R}|}
\]

Pearson-derived and Jaccard-derived representation-specific graphs are processed
independently. Their edge-recurrence counts are not pooled.

### 8.2 Majority Recurrence Threshold

For \(|\mathcal{R}|\) active representations, the majority threshold is:

\[
m
=
\left\lfloor

rac{|\mathcal{R}|}{2}

ight
floor
+1
\]

This threshold is derived from the number of active representations rather than
being independently fixed for each cancer type.

### 8.3 Global–Local Neighbour Intersection

The graph-resolution procedure operates independently for each similarity
metric.

For cell line \(i\):

\[
N_{\mathrm{global}}(i)
=
N_{r^*}(i)
\]

\[
N_{\mathrm{local}}(i)
=
igcup_{r\in T_i}
N_r(i)
\]

\[
N_{\mathrm{final}}(i)
=
N_{\mathrm{global}}(i)
\cap
N_{\mathrm{local}}(i)
\]

where \(r^*\) is the selected global representation and \(T_i\) is the relevant
local representation set for cell line \(i\).

The resolved graph is defined by the resulting graph-resolved neighbour sets. It
should not be described as simple pruning of the edge-recurrence network.

### 8.4 Isolates and Central Nodes

A graph isolate is a biological cell line with degree zero in the final resolved
graph.

Isolation means that no neighbour is selected by the graph-resolution rule. It
does not mean that the cell line has no transcriptomic proximity to any patient
tumour.

Central graph roles are described using graph quantities such as:

- degree
- betweenness centrality
- connected-component membership

Anchor candidates are drawn from highest-degree and highest-betweenness nodes
within connected components.

---

## 9. Marker Analysis

### 9.1 Graph-Derived Cell-Line Contrasts

Graph-defined cell-line roles provide the focal and reference groups used for
marker analysis.

The primary contrast classes include:

- isolate-based contrasts
- anchor-based contrasts

These are cell-line contrasts. They are not tumour-cluster versus cell-line
community contrasts.

### 9.2 Differential Expression Analysis

DESeq2 is applied to cell-line-only raw integer counts.

For gene \(j\) in sample \(i\):

\[
K_{ij}
\sim
\mathrm{NB}(\mu_{ij},lpha_j)
\]

with:

\[
\log(\mu_{ij})
=
\log(s_i)
+
eta_{j0}
+
eta_{j1}x_i
\]

where:

- \(s_i\) is the sample-specific size factor
- \(lpha_j\) is the gene-specific dispersion parameter
- \(x_i\) encodes the two-group graph-derived contrast

Reported log-fold changes are expressed on the \(\log_2\) scale.

### 9.3 Marker Selection

Genes passing the declared differential-expression and expression filters are
selected for marker analysis.

Selection uses quantities including:

- adjusted \(p\)-value
- absolute log2 fold change
- base expression requirements
- test-sample expression requirements
- recurrence across graph-derived contrasts

The exact criteria differ between isolate- and anchor-based contrasts and are
defined by the marker-analysis configuration and scripts.

---

## 10. Multicohort and Pan-Cancer Analysis

### 10.1 Pan-Cancer Marker Aggregation

Cancer-specific marker evidence is aggregated across breast cancer,
neuroblastoma, and retinoblastoma to construct the pan-cancer marker panel.

For cancer type \(c\), marker-source class \(s\), contrast \(j\), and gene \(g\),
let \(M_{c,s,j}\) denote the marker list passing the declared filters.

Within-cancer recurrence is:

\[
r_{c,s}(g)
=
\sum_j
\mathbf{1}
[g\in M_{c,s,j}]
\]

Genes recurring across multiple contrasts are distinguished from non-recurrent
candidates. Candidate genes are evaluated using the declared statistical,
fold-change, and expression criteria before construction of the final marker
panel.

The pan-cancer marker panel is therefore constructed from cancer-specific
DESeq2 marker evidence. It is not produced by merging all cohort matrices and
performing a new joint feature-selection step.

### 10.2 Functional Enrichment Analysis

Functional over-representation analysis is applied to the pan-cancer marker
panel using g:Profiler and an eligible background gene set.

The analysis evaluates biological processes and pathways represented among the
selected marker genes. Multiple-testing correction is applied using the
configured g:Profiler procedure.

### 10.3 Multicohort Visualisation

Multicohort UMAP and graph visualisations can be annotated by cancer type to
examine the organisation of cell-line and tumour profiles across cohorts.

UMAP is used for qualitative visual interpretation. Distances and apparent
clusters in the two-dimensional embedding are not treated as quantitative
evidence of discrete biological classes.

---

## 11. Key Output Files

### 11.1 Tumour-Neighbourhood Outputs

`Final_consensus_tumour_neighbourhoods_{direction}.tsv`

- contains cell-line–tumour neighbourhood records
- contains the corresponding `p_consensus` fraction
- records tumour-wise recurrence across the contributing clustering formulations

Interpretation:

A higher `p_consensus` means that the same cell-line–tumour neighbourhood
relationship occurs in a larger fraction of the configured clustering
formulations for that feature-distance representation.

### 11.2 Representation-Specific Patient-Referenced Graph Outputs

`cell_line_similarity_pairs_{direction}.tsv`

- pairwise cell-line similarities derived from threshold-restricted
  `p_consensus` profiles

`cell_line_similarity_graph_edges_{direction}.tsv`

- threshold-selected edges for the representation-specific patient-referenced
  graph

`cell_line_similarity_graph_node_summary_{direction}.tsv`

- node-level graph quantities for the representation-specific graph

`cell_line_similarity_undefined_pairs_{direction}.tsv`

- pairs with undefined similarity values
- reason for the undefined similarity where recorded

### 11.3 Graph-Resolution and Metric-Evaluation Outputs

Current metric-evaluation outputs include:

- `pearson_vs_jaccard_representation_graphs.tsv`
- `pearson_vs_jaccard_pairwise_similarity.tsv`
- `pearson_vs_jaccard_edge_agreement.tsv`
- `pearson_vs_jaccard_resolved_graph_comparison.tsv`
- `pearson_vs_jaccard_resolved_neighbours.tsv`
- `similarity_graph_provenance.tsv`
- `clustering_method_consensus_resolution.tsv`

These files report representation-specific graph structure, edge agreement,
resolved-neighbour agreement, graph provenance, and clustering-method
composition.

### 11.4 Marker and Functional-Analysis Outputs

The marker-analysis stage produces:

- DESeq2 results for graph-derived cell-line contrasts
- marker lists passing the declared filters
- cancer-specific marker recurrence results
- pan-cancer marker-panel outputs
- functional enrichment results for the marker panel

Exact filenames are defined by the active marker-analysis rules and should not be
replaced by legacy tumour-cluster or representativeness-scoring filenames.

---

## 12. Configuration

### 12.1 Profile-Based Configuration

The workflow uses YAML configuration with shared defaults and cancer-specific
profile sections.

Profile-specific values override inherited defaults through deep merging.

### 12.2 Feature-Distance Representations

Active cancer-specific directions are defined through configuration rather than
by scanning whatever output directories happen to exist.

The standard feature-selection grid contains:

- Variance
- HVG residual variance
- MAD
- MeanAbsDev
- Entropy
- PCA loadings
- Spearman connectivity
- MX score
- WGCNA kTotal

crossed with:

- Euclidean distance
- Pearson correlation distance

### 12.3 Patient-Referenced Graph Parameters

Production graph parameters are:

```yaml
patient_referenced_graph:
  p_consensus_threshold: 0.70
  similarity_quantile: 0.90
  similarity_metrics:
    - pearson
    - jaccard
```

Values such as 0.75 for the `p_consensus` threshold and 0.85 for the similarity
quantile were used only in configuration-perturbation tests and are not
production choices.

### 12.4 Clustering-Method Configuration

`patient_referenced_graph.clustering_methods_by_distance` declares the exact
clustering formulations contributing to the tumour-neighbourhood recurrence
calculation. `AGN_*` identifies the HC/k-means clustering branch (retained
internal compatibility identifier); `CCP_*` identifies ConsensusClusterPlus.
Identifiers match the `tumour_neighbourhoods/<direction>/<method_id>/` output
directories exactly.

Correlation (`n_methods = 4`):

```text
AGN_HC_expr_cell_tumour
AGN_HC_pca_cell_tumour
CCP_HC_expr_cell_tumour
CCP_HC_pca_cell_tumour
```

Euclidean (`n_methods = 8`):

```text
AGN_HC_expr_cell_tumour
AGN_HC_pca_cell_tumour
AGN_KM_expr_cell_tumour
AGN_KM_pca_cell_tumour
CCP_HC_expr_cell_tumour
CCP_HC_pca_cell_tumour
CCP_KM_expr_cell_tumour
CCP_KM_pca_cell_tumour
```

The `p_consensus` denominator is determined from these exact configured method
identities, never inferred from directories present on disk. Declared methods
whose outputs are missing, and discovered outputs that are not declared, both
raise explicit errors (`scripts/compute_tumour_neighbourhoods.R`,
`scripts/tumour_neighbourhood_p_consensus.R`). k-means formulations are
declared for Euclidean directions only, since k-means minimises Euclidean
inertia and is undefined under correlation distance.

### 12.5 Adaptive-Neighbourhood Configuration

The methodological values are:

```text
fraction = 0.10
minimum = 30
maximum = 200
```

A previous implementation audit found that these values were still hard-coded at
an active call site despite corresponding configuration entries. This wiring
should be revalidated after the current refactor. Documentation should describe
the configuration as the execution source only after that check passes.

---

## 13. Reproducibility and Validation

### 13.1 Random Seeds

Stochastic procedures use fixed random seeds where applicable. Seed values are
defined through configuration or analysis scripts.

### 13.2 Environment Management

Conda environments define software dependencies for the pipeline stages,
including:

- R-based transcriptomic and clustering analysis
- tumour-neighbourhood analysis
- Python-based graph analysis
- visualisation procedures

### 13.3 Logging

Pipeline rules write execution logs to rule- and profile-specific locations.
Logs record standard output, error output, and workflow failures needed for
reproducibility and debugging.

### 13.4 Graph Provenance

Metric-specific graph provenance should record at minimum:

- cancer profile
- feature-distance representation
- similarity metric
- `p_consensus` threshold
- similarity quantile
- computed similarity threshold
- number and identities of contributing clustering formulations
- number of candidate cell-line pairs
- number of defined and undefined similarities
- number of pairs above the threshold
- number of pairs tied at the threshold
- selected edge count
- selected edge fraction
- graph density

Pairwise similarity outputs should also record:

- `n_pairwise_active_tumours`
- shared selected-tumour count
- selected-tumour count for each cell line
- `undefined_similarity_reason`

### 13.5 Undefined Similarities

Undefined Pearson and Jaccard values are excluded from the similarity-quantile
population and cannot become graph edges.

For Pearson, a small pairwise active-tumour count and zero-variance restricted
vectors require explicit audit.

For Jaccard, an empty active-tumour union is recorded as undefined (`NA`).

### 13.6 Percentile Ties

The graph rule uses:

```text
similarity >= computed_similarity_threshold
```

All pairs tied at the configured percentile threshold are therefore selected.
The graph is not an exact top-10%-of-pairs graph.

### 13.7 Current Implementation Checks

The following points were validated on the regenerated outputs (2026-08-14):

1. Joint HC/k-means clustering outputs are connected to tumour-neighbourhood
   generation and `p_consensus`: the `tumour_nh_hc` / `tumour_nh_km` rules
   declare the JOINT `AGN_*` cluster files as inputs alongside the `CCP_*`
   files, and both branches produce method-specific neighbourhood tables.
2. `n_methods = 4` for correlation and `n_methods = 8` for Euclidean
   representations, confirmed for every configured representation by
   `scripts/audit_clustering_method_consensus_resolution.py`
   (`method_file_validation = exact_match` for all directions).
3. Cell-only and tumour-only clustering outputs never contribute: only
   `*_cell_tumour` identifiers are declarable, and undeclared method
   directories raise explicit errors.
4. The exact configured clustering-method identities are validated against
   realised inputs at both the neighbourhood-generation and `p_consensus`
   stages.
5. Missing declared formulations and unexpected discovered formulations are
   explicit errors (verified by fault-injection runs).
6. Adaptive-neighbourhood configuration wiring is unchanged
   (`adaptive_k`: fraction 0.10, minimum 30, maximum 200) and revalidated by
   the regenerated runs.
7. Ward.D2 compatibility for correlation-distance HC is established by the
   chord-distance construction `sqrt(2 * (1 - r))` (Section 5.3), matching
   the ConsensusClusterPlus correlation-geometry transform.
8. `tests/test_p_consensus_recurrence_fraction.py` reconstructs every
   `p_consensus` value of a Euclidean and a correlation representation
   directly from the method-specific tumour memberships
   (`recurrence / n_methods == p_consensus`, all pairs).

---

## 14. Glossary

**Active-tumour union:**
The pair-specific set of patient tumours selected by at least one of two cell
lines after the `p_consensus` threshold is applied.

**Biological cell line:**
A unique cell-line identity after profile-level replicates are combined where
applicable.

**Consensus clustering:**
Resampling-based clustering implemented with ConsensusClusterPlus.

**Edge-recurrence count:**
The number of active feature-distance representations in which a specific
cell-line edge is threshold-selected.

**Edge-recurrence fraction:**
The edge-recurrence count divided by the number of active feature-distance
representations.

**Feature-distance representation:**
A transcriptomic representation defined by a feature-selection method and a
dissimilarity measure.

**HC/k-means clustering:**
The workflow branch containing hierarchical clustering and k-means applied in
expression or PCA-reduced spaces. The term describes the algorithms in that
branch and does not imply a distinct biological class of clustering.

**HVG residual variance:**
A feature-ranking statistic based on residual variation after accounting for the
mean–variance relationship.

**Jaccard similarity:**
Intersection divided by union after threshold-restricted continuous
`p_consensus` fractions are converted to binary patient-tumour membership.

**PAM50:**
A breast-cancer-specific 50-gene expression signature used for molecular subtype
annotation.

**Patient-referenced cell-line graph:**
A graph whose nodes are biological cell lines and whose edges are selected from
cell-line similarity values derived from patient tumour-neighbourhood profiles.

**Pearson correlation distance:**
The dissimilarity \(1-
ho\), where \(
ho\) is the Pearson correlation
coefficient across selected genes.

**`p_consensus`:**
A tumour-wise recurrence fraction describing how often a particular
cell-line–tumour neighbourhood relationship occurs across the configured
clustering formulations within a fixed feature-distance representation.

**Threshold-restricted continuous `p_consensus` profile:**
The biological-cell-line mean-pooled tumour-wise `p_consensus` profile after
values below the production threshold are set to zero while qualifying values
remain continuous.

**Tumour neighbourhood:**
The adaptive set of patient tumours associated with a cell line under a
particular feature-distance representation and clustering formulation.

**VST (Variance-Stabilising Transformation):**
A transformation of count data that reduces the dependence of variance on mean
expression level.
