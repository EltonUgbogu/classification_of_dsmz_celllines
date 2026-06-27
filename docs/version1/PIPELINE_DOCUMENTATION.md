# =============================================================================
# UNSUPERVISED CLUSTERING AND TUMOUR NEIGHBOURHOOD ANALYSIS PIPELINE
# =============================================================================
#
# SUMMARY
#
# This document provides comprehensive documentation of a computational
# bioinformatics pipeline designed to investigate transcriptomic relationships between
# primary tumour samples and cancer cell line models. The pipeline addresses the
# fundamental question of whether established cancer cell lines faithfully recapitulate
# the molecular characteristics of their tissues of origin.
#
# -----------------------------------------------------------------------------
# TABLE OF CONTENTS
#
# 1. Scientific Background and Rationale
# 2. Pipeline Architecture Overview
# 3. Data Sources and Inputs
# 4. Feature Selection Methodology
# 5. Clustering Approaches
# 6. Tumour Neighbourhood Analysis
# 7. Cell Line Similarity Networks
# 8. Biological Characterisation
# 9. Pan-Cancer Integration
# 10. Key Output Files and Interpretation
# 11. Configuration System
# 12. Reproducibility Considerations
#
# -----------------------------------------------------------------------------
# 1. SCIENTIFIC BACKGROUND AND RATIONALE
#
# 1.1 The Cell Line Model Problem
#
# Cancer cell lines have served as the primary in vitro models for cancer research
# and drug development for over half a century. However, prolonged culture may
# lead to genetic drift, causing cell lines to accumulate alterations that
# distinguish them from their tissues of origin. This divergence raises concerns
# about the translational relevance of findings derived from cell line studies.
#
# 1.2 Research Objectives
#
# The pipeline addresses several key questions:
#
# 1. Which cell lines most faithfully represent primary tumour transcriptomic
#    profiles within each cancer type?
#
# 2. Are there cell lines that have diverged significantly from any known
#    tumour molecular state?
#
# 3. Can cell lines be grouped into functionally related communities based
#    on their shared tumour neighbourhood characteristics?
#
# 4. How robust are cell line-tumour relationships across different
#    analytical methodologies?
#
# 1.3 The Tumour Neighbourhood Concept
#
# For each cell line, the pipeline computes a "tumour neighbourhood" consisting
# of primary tumour samples that exhibit similar gene expression profiles. The
# size and composition of this neighbourhood provides a quantitative measure
# of cell line representativeness. Cell lines with large neighbourhoods spanning
# many patients are considered highly representative, while those with small
# or empty neighbourhoods may have limited utility as disease models.
#
# -----------------------------------------------------------------------------
# 2. PIPELINE ARCHITECTURE OVERVIEW
#
# 2.1 Layered Design
#
# The pipeline is organised into sequential analytical layers:
#
# Layer 0: Data Preprocessing and Feature Selection
# - Variance-stabilised transformation (VST) of RNA-seq count data
# - Multi-method highly variable gene (HVG) identification
# - Optional PAM50 molecular subtype gene signature extraction
#
# Layer 1: Agnostic Clustering
# - Hierarchical clustering with Ward linkage
# - K-means clustering across k=2 to k=8
# - Both PCA-reduced and full expression space approaches
# - Three sample views: cell-only, tumour-only, joint
#
# Layer 2: Consensus Clustering
# - ConsensusClusterPlus-based stability assessment
# - Bootstrap resampling for robustness evaluation
# - Optimal cluster number determination
#
# Layer 3: Tumour Neighbourhood Computation
# - Per-cell-line similarity scoring against all tumours
# - p_consensus calculation across parameter combinations
# - Cross-direction winner determination
#
# Layer 4: Similarity Network Analysis
# - DSMZ cell line similarity graph construction
# - Louvain and Leiden community detection
# - Component characterisation
#
# Layer 5: Biological Interpretation
# - Differential expression analysis (DESeq2)
# - Representativeness scoring
#
# Layer 6: Pan-Cancer Integration
# - Cross-disease comparative analysis
# - Disease-aware visualisation
#
# 2.2 Multi-Direction Ensemble Strategy
#
# A central innovation of this pipeline is the systematic evaluation of multiple
# "directions" - combinations of feature selection methods and distance metrics.
# This ensemble approach distinguishes robust biological findings from
# method-specific artefacts.
#
# The pipeline evaluates 8 feature selection methods:
# - Variance: Standard variance across samples
# - MAD: Median Absolute Deviation (robust to outliers)
# - MeanAbsDev: Mean Absolute Deviation
# - Entropy: Information-theoretic measure of distribution complexity
# - PCA: Principal component loadings
# - Spearman: Rank correlation-based selection
# - MX: Maximum expression-based selection
# - kTotal: Network connectivity-based selection
#
# Combined with 2 distance metrics:
# - Euclidean: Appropriate for VST-transformed data
# - Correlation: Captures relative expression patterns
#
# This yields 16 direction combinations per cancer type (plus optional PAM50
# directions for breast cancer).
#
# -----------------------------------------------------------------------------
# 3. DATA SOURCES AND INPUTS
#
# 3.1 Tumour Data: 
#
# The pipeline utilises RNA-seq expression profiles from TCGA,and GEO (Gene Expression Omnibus)
# the largest publicly available cancer genomics resource. The following cohorts are
# supported:
#
# - BRCA: Breast invasive carcinoma (n approximately 1100 samples)
# - NBL: Neuroblastoma (available through TARGET)
# - RBL: Retinoblastoma
#
# 3.2 Cell Line Data: DSMZ Collection
#
# The Leibniz Institute DSMZ maintains one of the world's largest collections
# of authenticated human and animal cell lines. Expression profiles are
# processed identically to Patient Tumour data to enable direct comparison.
#
# 3.3 Input File Requirements
#
# The pipeline requires:
#
# 1. vst_joint_rds: Combined VST-normalised expression matrix containing
#    both tumour and cell line samples (genes x samples)
#
# 2. cell_vst_rds: Cell line-only expression matrix
#
# 3. tumour_vst_rds: Tumour-only expression matrix
#
# 4. dsmz_meta_csv: Cell line metadata including disease type, tissue
#    of origin, and authentication status
#
# 5. PAM50 gene list (optional, breast cancer only): 50-gene signature
#    for molecular subtyping
#
# -----------------------------------------------------------------------------
# 4. FEATURE SELECTION METHODOLOGY
#
# 4.1 Rationale for Feature Selection
#
# With approximately 20,000 protein-coding genes in the human genome,
# dimensionality reduction through feature selection is essential for:
#
# 1. Removing uninformative genes (housekeeping, low expression)
# 2. Reducing computational burden
# 3. Improving signal-to-noise ratio
# 4. Focusing on biologically variable genes
#
# 4.2 Variance-Based Methods
#
# Variance: The simplest approach computes the variance of each gene's
# expression across all samples. Genes with high variance are likely to
# capture meaningful biological differences, though this metric can be
# dominated by highly expressed genes.
#
# MAD (Median Absolute Deviation): A robust alternative to variance that
# is less sensitive to outliers:
#
#     MAD(x) = median(|x - median(x)|)
#
# This metric is preferred when samples may contain technical artefacts or
# extreme biological outliers.
#
# 4.3 Information-Theoretic Methods
#
# Entropy: Quantifies the complexity of a gene's expression distribution
# using Shannon entropy. Genes with bimodal or multimodal distributions
# (characteristic of regulatory switches) receive high entropy scores.
#
# 4.4 Dimensionality Reduction Methods
#
# PCA Loadings: Principal Component Analysis identifies linear combinations
# of genes that capture maximal variance. Genes with high loadings on the
# dominant principal components are selected.
#
# 4.5 The Top-500 Gene Strategy
#
# After applying each method, the pipeline retains the top 500 genes by rank.
# This number represents an empirical balance between capturing biological
# signal and avoiding noise from lowly-expressed or invariant genes.
#
# -----------------------------------------------------------------------------
# 5. CLUSTERING APPROACHES
#
# 5.1 Hierarchical Clustering
#
# Hierarchical clustering constructs a dendrogram representing sample
# relationships through iterative merging of clusters. The pipeline employs
# Ward's minimum variance method, which minimises within-cluster variance
# at each merge step.
#
# Key Parameters:
# - Distance metric: Euclidean or correlation-based
# - Linkage method: Ward
# - Number of clusters (k): Evaluated from k=2 to k=8
#
# Advantages:
# - Produces interpretable dendrograms
# - No need to pre-specify k
# - Reveals hierarchical structure
#
# 5.2 K-Means Clustering
#
# K-means partitioning assigns samples to k clusters by minimising the
# sum of squared distances to cluster centroids. The algorithm iterates
# between:
#
# 1. Assignment step: Assign each sample to nearest centroid
# 2. Update step: Recompute centroids as cluster means
#
# Note: K-means requires Euclidean distance and is therefore only
# applied to "_euc" directions.
#
# 5.3 PCA-Reduced vs. Direct Clustering
#
# The pipeline evaluates both approaches:
#
# PCA-Reduced: Expression data is first projected onto the top 20
# principal components before clustering. This:
# - Reduces noise by discarding low-variance components
# - Improves computational efficiency
# - May lose gene-level interpretability
#
# Direct: Clustering is performed on the full (filtered) expression
# matrix, preserving gene-level information but potentially including noise.
#
# 5.4 Sample Views
#
# Three clustering configurations are evaluated:
#
# 1. Cell-only: Clusters cell lines independently to identify intrinsic
#    groupings
#
# 2. Tumour-only: Clusters tumours independently to identify molecular
#    subtypes
#
# 3. Joint (cell_tumour): Clusters cell lines and tumours together to
#    assess integration
#
# -----------------------------------------------------------------------------
# 6. TUMOUR NEIGHBOURHOOD ANALYSIS
#
# 6.1 Core Algorithm
#
# For each cell line, the tumour neighbourhood is computed as follows:
#
# 1. Identify which tumour samples cluster together with the cell line
#    across multiple clustering parameter combinations
#
# 2. Count the frequency with which each tumour appears in the same
#    cluster as the cell line
#
# 3. Retain tumours exceeding a frequency threshold as neighbourhood
#    members
#
# 6.2 The p_consensus Metric
#
# The p_consensus value quantifies the proportion of clustering analyses
# in which a tumour-cell line relationship is observed:
#
#     p_consensus = (number of analyses where tumour T clusters with cell line C) /
#                   (total number of analyses)
#
# High p_consensus (e.g., >0.7) indicates robust relationships that persist
# across methodological choices.
#
# 6.3 Winner Determination
#
# After computing p_consensus across all directions, the pipeline identifies
# the "winning" direction - the combination of feature selection and distance
# metric that produces the most robust tumour neighbourhoods. Selection
# criteria include:
#
# 1. Fraction of cell lines with p_consensus > threshold
# 2. Concordance across directions
# 3. PCA-derived composite scoring
#
# -----------------------------------------------------------------------------
# 7. CELL LINE SIMILARITY NETWORKS
#
# 7.1 Graph Construction
#
# Cell lines sharing similar tumour neighbourhoods may represent related
# biology. The pipeline constructs similarity networks where:
#
# - Nodes: Individual cell lines
# - Edges: Weighted by tumour neighbourhood overlap (Jaccard similarity)
#
#     Jaccard(A,B) = |TN_A intersection TN_B| / |TN_A union TN_B|
#
# where TN denotes tumour neighbourhood.
#
# 7.2 Community Detection
#
# Louvain Algorithm: A greedy modularity optimisation algorithm that
# iteratively merges nodes into communities to maximise:
#
#     Q = (1/2m) sum_ij [A_ij - (k_i * k_j)/(2m)] delta(c_i, c_j)
#
# where Q is modularity, A is the adjacency matrix, and c denotes community
# assignment.
#
# Leiden Algorithm: An improved variant that guarantees well-connected
# communities and addresses the resolution limit of Louvain.
#
# 7.3 Biological Interpretation
#
# Cell line communities identified through network analysis may correspond to:
# - Shared tissue of origin
# - Common driver mutations
# - Similar drug sensitivity profiles
#
# -----------------------------------------------------------------------------
# 8. BIOLOGICAL CHARACTERISATION
#
# 8.1 Differential Expression Analysis
#
# Using DESeq2, the pipeline identifies genes differentially expressed
# between:
#
# 1. Tumour clusters identified by consensus clustering
# 2. Cell line communities identified by network analysis
#
# DESeq2 employs a negative binomial generalised linear model:
#
#     K_ij ~ NB(mu_ij, alpha_j)
#     log2(mu_ij) = beta_j0 + beta_j1 * x_i + ...
#
# 8.2 Functional Enrichment
#
# Marker gene lists are submitted to enrichment analysis using:
#
# - Gene Ontology (GO): Biological Process, Molecular Function,
#   Cellular Component
#
# - Pathway databases: KEGG, Reactome
#
# Enrichment significance is assessed using hypergeometric tests with
# Benjamini-Hochberg correction for multiple testing.
#
# 8.3 Representativeness Scoring
#
# A composite representativeness score is computed for each cell line:
#
#     Score = p_consensus * specificity_weight
#
# where specificity weights reward cell lines with focused (cluster-specific)
# neighbourhoods versus diffuse (pan-cluster) neighbourhoods.
#
# -----------------------------------------------------------------------------
# 9. PAN-CANCER INTEGRATION
#
# 9.1 Rationale
#
# Pan-cancer analysis extends the single-disease framework to identify:
#
# 1. Cell lines that model biology common to multiple cancer types
# 2. Disease-specific versus shared transcriptomic programmes
# 3. Cross-disease cell line relationships
#
# 9.2 Joint Feature Space
#
# Expression matrices from multiple disease cohorts are merged, and
# joint feature selection is performed. This ensures that comparisons
# occur in a unified analytical space.
#
# 9.3 Disease-Aware Visualisation
#
# UMAP projections and network visualisations are annotated with
# disease type, enabling assessment of disease-specific clustering
# versus cross-disease integration.
#
# -----------------------------------------------------------------------------
# 10. KEY OUTPUT FILES AND INTERPRETATION
#
# 10.1 Tumour Neighbourhood Results
#
# Final_consensus_tumour_neighbourhoods_{direction}.tsv
# - Columns: cell_line, tumour_id, p_consensus
# - Interpretation: Higher p_consensus indicates stronger cell line-tumour
#   relationship
#
# p_consensus_best_cell_lines_ranked.tsv
# - Cell lines ranked by overall representativeness
# - Top-ranked lines are the best tumour models
#
# 10.2 Similarity Network Results
#
# DSMZ_DSMZ_graph_edges_{direction}.tsv
# - Network edge list with similarity weights
# - For import into network visualisation software
#
# DSMZ_DSMZ_Louvain_vs_Leiden_community_table_{direction}.tsv
# - Community assignments for each cell line
# - Comparison of Louvain and Leiden results
#
# 10.3 Characterisation Results
#
# cell_line_characterisation.tsv
# - Per-cell-line summary including cluster assignment and p_consensus
# - Primary reference for cell line selection
#
# cluster_DE_markers.rds
# - Differential expression results for tumour clusters
# - R data object for downstream analysis
#
# -----------------------------------------------------------------------------
# 11. CONFIGURATION SYSTEM
#
# 11.1 Profile-Based Configuration
#
# The pipeline uses YAML configuration files organised by disease profile:
#
# profiles:
#   brca:
#     analysis:
#       cancer_type: "BRCA"
#       use_pam50: true
#     paths:
#       vst_joint_rds: "data/brca/vst_joint.rds"
#   nbl:
#     analysis:
#       cancer_type: "NBL"
#       use_pam50: false
#     paths:
#       vst_joint_rds: "data/nbl/vst_joint.rds"
#
# 11.2 Default Inheritance
#
# Profile-specific settings override defaults through deep merging,
# allowing compact profile definitions that inherit common parameters.
#
# -----------------------------------------------------------------------------
# 12. REPRODUCIBILITY CONSIDERATIONS
#
# 12.1 Random Seeds
#
# All stochastic algorithms (k-means, UMAP) use fixed random seeds
# (default: 42) specified in the configuration file.
#
# 12.2 Environment Management
#
# Conda environment YAML files specify exact package versions:
#
# - tcga-r-env.yaml: R packages for clustering and differential expression
# - tumour_nh_qc.yaml: UMAP and visualisation tools
# - python-graph-env.yaml: Network analysis packages
#
# 12.3 Logging
#
# Comprehensive logs are written to profile-specific directories:
#
#     logs/{profile_name}/{rule_name}.log
#
# Logs capture all stdout/stderr from analysis scripts for debugging
# and audit purposes.
#
# -----------------------------------------------------------------------------
# APPENDIX: GLOSSARY OF KEY TERMS
#
# Agnostic Clustering: Unsupervised clustering without reference to
# known biological labels.
#
# Consensus Clustering: Stability assessment through repeated clustering
# of resampled data.
#
# Direction: A specific combination of feature selection method and
# distance metric.
#
# HVG (Highly Variable Genes): Genes with high expression variability
# across samples.
#
# p_consensus: Proportion of analyses supporting a specific cell line-
# tumour relationship.
#
# PAM50: 50-gene signature for breast cancer molecular subtyping.
#
# Tumour Neighbourhood: Set of primary tumours transcriptomically similar
# to a given cell line.
#
# VST (Variance Stabilising Transformation): Normalisation method that
# stabilises variance across the expression range.
#
# -----------------------------------------------------------------------------
# CITATION AND REFERENCES
#
# This pipeline integrates methodology from:
#
# 1. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold
#    change and dispersion for RNA-seq data with DESeq2. Genome Biology.
#
# 2. Wilkerson MD, Hayes DN. (2010). ConsensusClusterPlus: a class discovery
#    tool with confidence assessments and item tracking. Bioinformatics.
#
# 3. Parker JS, et al. (2009). Supervised risk predictor of breast cancer
#    based on intrinsic subtypes. Journal of Clinical Oncology.
#
