# Concepts and interpretation

This document defines terminology used in the transcriptomic similarity workflow. The terms connect biological interpretation, graph-based computation, and high-dimensional expression geometry.

## Biological concepts

The workflow evaluates transcriptomic similarity between DSMZ cancer cell lines and patient tumour RNA-seq cohorts. The central biological question is whether each cell line occupies a patient-referenced expression neighbourhood concordant with its annotated cancer type.

Key biological concepts include:

* transcriptomic heterogeneity
* transcriptomic substructure
* tumour-type identity
* disease-specific expression states
* cancer-type-specific expression programmes
* breast cancer transcriptomic heterogeneity
* neuroblastoma transcriptional states
* retinoblastoma transcriptional states
* tumour purity
* immune and stromal admixture
* tumour microenvironmental signal
* marker-gene expression
* differential gene expression

## Computational concepts

The workflow uses computational methods to structure comparisons between DSMZ cancer cell lines and patient tumour RNA-seq cohorts. Expression profiles are represented through multiple feature-distance settings and carried through dimensionality reduction, high-dimensional neighbour search, tumour-neighbourhood construction, consensus clustering, graph construction, graph resolution, community detection, and reciprocal ranking analyses.

Key computational concepts include:

* high-dimensional transcriptomic data
* feature representation
* feature-distance representation
* feature-representation and distance-metric
* dimensionality reduction
* UMAP-based qualitative embedding
* high-dimensional nearest-neighbour search
* adaptive tumour-neighbourhood construction
* consensus clustering
* multi-representation consensus
* patient-referenced graph construction
* graph-based neighbourhood resolution
* global-local neighbour intersection
* connected-component analysis
* isolate detection
* bridge-like anchor detection
* centrality-based graph characterisation
* community detection
* Louvain community detection
* Leiden community detection
* reciprocal tumour-cell-line ranking
* rank-based transcriptomic proximity
* correlation-based similarity scoring
* distance-metric sensitivity
* workflow reproducibility
* workflow provenance

The patient-referenced graph should be treated as an inference layer over tumour-neighbourhood evidence, not as a direct biological interaction network.

Isolates indicate lack of stable resolved neighbours under the declared graph-resolution criteria. They do not indicate absence of tumour similarity.

## Statistical and mathematical concepts

After feature sets are defined, the analysis depends on distance metrics, graph thresholds, and high-dimensional nearest-neighbour behaviour. These factors should be considered when interpreting graph-resolved neighbours, isolates, ranking metrics, and enrichment results.

Key statistical and mathematical concepts include:

* threshold selection
* consensus threshold
* edge-recurrence threshold
* adaptive neighbourhood size
* distance concentration
* high-dimensional expression geometry
* normalisation
* variance stabilisation
* batch correction
* edge recurrence
* graph sparsity
* majority-style thresholding
* global-local evidence integration
* balanced accuracy
* Wilson confidence interval
* multiple-testing correction
* g:SCS correction
* Benjamini-Hochberg adjustment

## Methodological considerations and future directions

Important methodological considerations include:

* threshold selection
* feature-set selection
* distance-metric selection
* high-dimensional nearest-neighbour uncertainty
* threshold-constrained graph-resolution rules
* limited pathway-level enrichment signal in some marker sets
* need for independent validation
* need for subtype harmonisation
* need for experimental validation

## Interpretive framework

The main interpretation of the workflow can be organised around the following points:

1. Patient-referenced transcriptomic similarity can recover broad cancer-type identity from RNA-seq expression profiles.
2. Cell-line model prioritisation can vary within cancer types, reflecting tumour heterogeneity and transcriptomic substructure.
3. Graph-resolved neighbourhoods provide an interpretable edge-recurrence-thresholded graph representation of patient-neighbourhood evidence.
4. Feature-distance choices can alter which patient-tumour expression states are emphasised in graph-resolved neighbourhoods, influencing how cell-line cancer-type concordance and model prioritisation are interpreted.
5. Marker-gene and functional analyses provide candidate biological follow-up.
6. Threshold selection, feature-set selection, distance-metric selection, and high-dimensional nearest-neighbour geometry remain important methodological considerations that can be further optimised.
7. Isolates and bridge-like anchors should first be interpreted as graph-resolution outcomes rather than direct biological absence-or-presence claims.
