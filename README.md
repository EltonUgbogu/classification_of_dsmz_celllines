# Classification of DSMZ Cell Lines to Cancer Types Based on Transcriptomic Analysis

This repository contains the Snakemake workflow used to prioritise DSMZ cancer cell-line models by their transcriptomic proximity to patient tumour RNA-seq cohorts. The workflow maps DSMZ cell lines to patient tumours by transcriptomic proximity; these mappings provide computational evidence of relative transcriptomic similarity and should be interpreted alongside biological context and experimental validation.

The workflow constructs feature-distance representations, infers patient tumour neighbourhoods, builds multi-representation consensus networks, resolves cell-line neighbours, derives graph-informed marker evidence, constructs a pan-cancer feature panel, and performs reciprocal tumour--cell-line ranking. These outputs provide relative transcriptomic model prioritisation within the configured data and feature space.

The repository is maintained as a code-focused workflow repository. Raw sequencing data, large expression matrices, reference indexes, complete result trees, runtime logs, and most generated outputs are excluded from the main Git history. A tiny synthetic example is included so that repository structure and input schemas can be tested without access to the real datasets.

**Navigation:** [Key terms](#key-terms) | [Scientific aim](#scientific-aim) | [Workflow overview](#workflow-overview) | [Current analytical design](#current-analytical-design) | [Repository structure](#repository-structure) | [Quick start](#quick-start) | [Installation](#installation) | [Synthetic minimal example](#synthetic-minimal-example) | [Data and outputs](#data-and-outputs) | [Configuration](#configuration) | [Running the full workflow](#running-the-full-workflow) | [Expected output classes](#expected-output-classes) | [Reproducibility checks](#reproducibility-checks) | [Interpretation guide](#interpretation-guide) | [Scope assumptions and validation requirements](#scope-assumptions-and-validation-requirements) | [Citation and contact](#citation-and-contact)

## Key terms

- **Feature-distance representation:** a selected gene set, expression space, and distance or similarity metric considered together as one representation.
- **Tumour neighbourhood:** the adaptive set of retained patient tumour samples assigned near a focal cell-line profile within a feature-distance representation.
- **p-consensus fraction:** the tumour-wise recurrence fraction, within a fixed feature-distance representation, of the configured eligible clustering formulations in which a tumour sample occurs in the adaptive neighbourhood of a focal cell-line profile. Eligible formulations are the JOINT cell-line + tumour outputs of HC/k-means clustering and of ConsensusClusterPlus: 8 formulations for Euclidean directions and 4 for correlation directions, declared in `patient_referenced_graph.clustering_methods_by_distance`.
- **HC/k-means clustering:** hierarchical clustering and k-means applied to each feature-distance representation in either expression space or PCA-reduced space (internal identifiers `AGN_*`, rules `agnostic_cluster_*`, outputs under `agnostic_clustering/`). A clustering formulation parallel to ConsensusClusterPlus consensus clustering; the two meet at the p-consensus stage via their method-specific tumour neighbourhoods.
- **Patient-referenced cell-line similarity graph:** a weighted, undirected graph constructed separately for each feature-distance representation from Pearson correlation between tumour-wise p-consensus profiles. It reduces tumour-mediated cell-line similarity to a graph representation for graph-based comparison and resolved-neighbour analysis.
- **Multi-representation consensus network:** the primary cross-representation graph layer, which selects an edge when it recurs in at least `m = max(2, floor(|R| / 2) + 1)` configured feature-distance representations. A separate union network includes every edge observed in at least one representation as a sensitivity view.
- **Resolved graph:** a patient-referenced cell-line graph with resolved edge set `E_resolved`. For each focal cell-line profile `v`, retained incident edges are defined by the intersection of the global-best and local-best incident edge sets:

  ```text
  E_resolved(v) = E_G(v) ∩ E_L(v)
  ```
- **Resolved neighbour:** a cell-line profile connected to a focal cell-line profile in the resolved graph.
- **Global-best representation:** the cohort-level feature–distance representation with the largest overall fraction of `p_consensus` values meeting the declared threshold across the configured cell-line profiles.
- **Local-best representation set:** the one or more feature–distance representations with the largest fraction of `p_consensus` values meeting the declared threshold for an individual cell-line profile.
- **Component anchor:** a profile selected from a non-isolate resolved component through highest-degree or highest unnormalised-betweenness criteria.
- **Marker-evidence stratum:** the graph-derived evidence category used during marker aggregation: anchor-associated or isolate-associated marker evidence.
- **Pan-cancer feature panel:** the graph-informed marker-derived feature set used for pan-cancer cell-line similarity analysis and reciprocal tumour--cell-line ranking.
- **Model prioritisation:** relative ordering of candidate cell-line models according to patient-referenced transcriptomic proximity in the configured feature space.

## Scientific aim

Cancer cell lines are widely used as laboratory models; however, selecting a cell line by cancer-type label alone does not establish how well its transcriptome represents patient tumours. Patient tumours and cell-line profiles can vary substantially within the same cancer type, and this variation can affect which models appear most representative.

This workflow evaluates DSMZ cell-line profiles in a patient-referenced expression space. It constructs unsupervised feature–distance representations and quantifies patient-referenced tumour-neighbourhood proximity between DSMZ cell-line profiles and retained tumour samples across those representations. Resolved-graph construction then selects cell-line adjacencies that satisfy both cohort-level and cell-line-specific representation criteria, while graph-derived marker contrasts connect resolved graph structure to pan-cancer feature construction.

Reciprocal ranking provides label-free DSMZ model ranking within a cohort-informed, graph-derived marker feature space. Cancer-type labels are not used to compute tumour–cell-line similarity scores or rank orderings; BRCA, NBL, and RBL annotations are applied afterwards to evaluate post hoc cancer-type annotation concordance.

## Workflow overview

The workflow performs the following stages:

1. Validate raw patient-tumour count matrices and gene identifiers.
2. Estimate tumour purity from a tumour-only expression representation derived inside the purity stage, then apply the configured threshold to the original raw tumour counts.
3. Join only purity-retained raw tumour counts with raw DSMZ counts, harmonise genes, apply ComBat-seq, and generate joint DESeq2 variance-stabilised matrices.
4. Derive profile-specific unsupervised feature-selection gene lists from joint VST expression using variance/HVG, MAD/MeanAbsDev, entropy, PCA-loading, and Spearman-, MX-, and WGCNA-derived connectivity criteria.
5. Combine feature sets, expression spaces, and distance or similarity metrics into feature-distance representations.
6. Apply the eligible clustering formulations to each representation: HC/k-means clustering (hierarchical clustering and k-means in expression space and PCA-reduced space) and ConsensusClusterPlus consensus clustering, each producing a JOINT cell-line + tumour cluster assignment.
7. Infer adaptive patient tumour neighbourhoods for focal cell-line profiles, one neighbourhood per eligible clustering formulation.
8. Aggregate the tumour-wise p-consensus recurrence fraction across the configured formulation set within each representation (Euclidean directions: 8 formulations; correlation directions: 4).
9. Construct representation-specific patient-referenced cell-line graphs.
10. Construct a majority-threshold multi-representation consensus network and a union sensitivity network.
11. Resolve cell-line neighbours by intersecting the global-best neighbour set with the pooled local-best neighbour set for the cell-line-specific tied representations, then render retained adjacencies as undirected edges.
12. Identify graph-derived isolate and anchor profiles for differential-expression analysis.
13. Run DESeq2 isolate and anchor contrasts using cancer-type-specific cell-line count matrices.
14. Aggregate contrast-level marker evidence into pan-cancer feature panel.
15. Perform reciprocal tumour--cell-line ranking in the configured pan-cancer feature space.
16. Construct the pan-cancer cell-line-only \(k\)-nearest-neighbour similarity network.
17. Compute Louvain and Leiden community assignments and evaluate both algorithms across configured resolution sweeps.
18. Construct and evaluate functional-enrichment query and background sets for marker interpretation.

## Current analytical design

### Patient-referenced graph resolution

Tumour-neighbourhood assignments are calculated separately within each feature-distance representation, once per configured clustering formulation (HC/k-means clustering and ConsensusClusterPlus, JOINT outputs only), and the p-consensus fraction records each pair's recurrence across that formulation set. The representation-specific patient-referenced cell-line similarity graph is then constructed from Pearson correlation between tumour-wise p-consensus profiles. The primary multi-representation consensus network selects edges recurring in at least the configured majority threshold across representations, using \(m = \max(2, \lfloor |R|/2 \rfloor + 1)\). The union network instead includes every edge observed in at least one representation.

Let \(E_r\) denote the undirected edge set of representation \(r\). The global-best representation \(g\) maximises cohort-level `frac_ge_thr`; the resolver orders ties by median `p_consensus` and, when supplied as a distinct field, mean `p_consensus`. Its focal edge set is \(E_G(v) = \{\{v,u\} \in E_g\}\). For each focal cell-line profile \(v\), the local-best set \(L(v)\) contains every representation tied for its highest per-cell-line `frac_ge_thr`, and its pooled local edge set is \(E_L(v) = \bigcup_{r \in L(v)} \{\{v,u\} \in E_r\}\). The resolver retains \(E_G(v) \cap E_L(v)\). The final graph takes the union of these retained focal sets over all \(v\) and renders each retained pair as one undirected edge.

### Graph-derived DESeq2 contrasts

Differential-expression analysis uses cell-line-only raw integer count matrices within each cancer type. Patient tumour samples are not included in the DESeq2 contrasts. Before model fitting, genes with a total raw count below 10 across the relevant cell-line count matrix are removed. Marker retention additionally requires at least 10 DESeq2-normalised counts in the focal isolate or anchor profile.

| Contrast type | Focal profile | Reference profiles | Marker filter |
| --- | --- | --- | --- |
| isolate | one degree-zero resolved-graph profile | all other same-cancer-type cell-line profiles | adjusted \(p \leq 0.01\), \(|\log_2\mathrm{FC}| \geq 1.5\), `baseMean >= 10`, maximum 50 markers |
| anchor | one graph-selected anchor profile | same-cancer-type profiles outside the focal resolved component | adjusted \(p \leq 0.05\), \(|\log_2\mathrm{FC}| \geq 1.0\), `baseMean >= 10`, maximum 200 markers |

For anchor contrasts, non-anchor profiles from the focal component are excluded from the reference set. Marker-list caps are maxima, not required list sizes.

### Pan-cancer feature selection

Marker recurrence is counted separately within each cancer type and marker-evidence stratum.

- Genes observed in at least two filtered contrast-level marker lists are included directly as recurrent evidence.
- Genes observed once are separated into singleton candidates or non-recurrent candidates according to the number of available contrast lists in the relevant cancer type and marker-evidence stratum.
- Singleton and non-recurrent candidate genes are accepted only when they pass all three empirical evidence criteria:
  - minimum adjusted \(p\)-value at or below the candidate-pool first quartile;
  - median absolute \(\log_2\) fold change at or above the candidate-pool third quartile;
  - median `baseMean` at or above the candidate-pool median.
- Empirical thresholds are estimated within cancer type, marker-evidence stratum, and candidate-pool type.
- The final panel is the union of recurrent genes, accepted singleton candidates, and accepted non-recurrent candidates.

## Repository structure

The public repository contains workflow source code, configuration, small resources, environment definitions, documentation, tests, and the synthetic minimal example.

```text
.
├── README.md                          # project overview and usage notes
├── Snakefile                          # main Snakemake workflow
├── config/                            # workflow configuration and study-design files
├── rules/                             # included Snakemake rule modules
├── scripts/                           # R, Python, and shell scripts called by rules
├── R/                                 # shared R helper functions
├── envs/                              # Conda environment specifications
├── resources/                         # small curated resources and identifier maps
├── preprocessing_and_quality_control/ # preprocessing workflows for tumour-purity analysis, batch correction, and normalisation
├── profiles/                          # local and SLURM Snakemake profiles
├── docs/                              # project documentation
├── examples/minimal/                  # tiny synthetic schema-compatible inputs
└── tests/                             # repository checks and test-related files
```

## Quick start

Clone the repository and create the Snakemake driver environment:

```bash
git clone --depth 1 https://github.com/EltonUgbogu/classification_of_dsmz_celllines.git
cd classification_of_dsmz_celllines
conda env create -f envs/smk.yaml
conda activate smk
```

List the available rules:

```bash
snakemake --list --config pipeline_profile=brca
```

Run the bundled synthetic validation target:

```bash
snakemake -n --use-conda test_minimal_inputs --config pipeline_profile=brca
snakemake --use-conda --cores 1 test_minimal_inputs --config pipeline_profile=brca
```

The synthetic target checks the expected matrix and metadata schemas and the inclusion of the rule in the main workflow. It does not execute the full biological analysis or produce scientifically interpretable results.

## Installation

The driver Snakemake environment is defined by [`envs/smk.yaml`](envs/smk.yaml) and uses the Conda environment name `smk`:

```bash
conda env create -f envs/smk.yaml
conda activate smk
```

The workflow uses per-rule Conda environments through `--use-conda`. Tracked environment definitions include:

- [`envs/tcga-r-env.yaml`](envs/tcga-r-env.yaml): the principal R analysis environment, including DESeq2, clustering, ranking, enrichment, and plotting dependencies;
- [`envs/python-graph-env.yaml`](envs/python-graph-env.yaml): Python graph construction and plotting utilities;
- [`envs/tumour_nh_qc.yaml`](envs/tumour_nh_qc.yaml): tumour-neighbourhood quality control and UMAP dependencies;
- [`envs/r-base.yaml`](envs/r-base.yaml): lightweight R rules where the full analysis environment is unnecessary.

Additional preprocessing environments are stored under [`preprocessing_and_quality_control/envs/`](preprocessing_and_quality_control/envs/).

Install Conda, Mamba, or an equivalent compatible environment manager separately, then construct the tracked environments from the YAML definitions.

## Synthetic minimal example

The repository includes tiny synthetic BRCA, NBL, and RBL inputs under [`examples/minimal/`](examples/minimal/). They contain:

- synthetic Ensembl-like gene identifiers;
- synthetic cell-line and tumour sample identifiers;
- small numeric VST-like matrices stored as RDS files;
- metadata tables describing the synthetic samples.

The matrices reproduce the object class expected by the validation target: a numeric R `matrix`/`array` with genes as rows and samples as columns. They are not biologically meaningful.

Regenerate the synthetic inputs:

```bash
Rscript examples/minimal/make_synthetic_minimal_inputs.R
```

Validate them directly:

```bash
Rscript scripts/validate_minimal_example_inputs.R examples/minimal
```

Validate them through Snakemake:

```bash
snakemake -n --use-conda test_minimal_inputs --config pipeline_profile=brca
snakemake --use-conda --cores 1 test_minimal_inputs --config pipeline_profile=brca
```

The Snakemake rule creates `results/tests/minimal_inputs.ok` as a generated touch file. It should not be committed.

## Data and outputs

Raw FASTQ, BAM/CRAM, SRA files, reference genomes, STAR indexes, large count matrices, complete VST matrices, full result trees, runtime logs, and most generated outputs are not redistributed in the main repository.

Abbreviations used by the project include TCGA, GDC, GEO, SRA, ENA, and TARGET.

## Configuration

The main workflow configuration is [`config/config.yaml`](config/config.yaml). Cohort design, candidate-inference settings, and study-level metadata are recorded in [`config/study_design.yaml`](config/study_design.yaml). Paths are interpreted relative to the repository root unless explicitly provided as absolute paths.

The configuration separates biological study design, data interfaces, analysis parameters, and runtime paths. Before a full run:

- configure profile-specific inputs under `profiles.<profile>.paths` in [`config/config.yaml`](config/config.yaml);
- make the required count matrices, expression matrices, metadata tables, gene-identifier references, tumour-purity outputs, and DSMZ cell-line annotations available;
- confirm that preprocessing outputs for BRCA, NBL, and RBL point to the intended retained tumour and DSMZ cell-line expression objects;
- ensure that cancer-type-specific DESeq2 inputs use raw, non-negative integer count matrices rather than transformed expression values;
- verify that feature-selection, tumour-neighbourhood, graph-resolution, ranking, network, and enrichment parameters match the intended analysis profile;
- use repository-relative paths for reproducible workflow inputs wherever possible;
- use absolute paths or an untracked local configuration override only for machine-specific storage locations;


`pipeline_profile` selects the active analysis profile:

| profile | scope | notes |
| --- | --- | --- |
| `brca` | breast cancer | patient-referenced BRCA workflow using retained TCGA-BRCA tumour samples and DSMZ BRCA cell-line profiles |
| `nbl` | neuroblastoma | patient-referenced NBL workflow using retained NBL tumour samples and DSMZ NBL cell-line profiles |
| `rbl` | retinoblastoma | patient-referenced RBL workflow using retained RBL tumour samples and DSMZ RBL cell-line profiles |
| `heme` | standalone haematological malignancy analysis | uses configured LAML and CLL tumour-expression inputs; only its cell-line profiles provide auxiliary context to the BRCA/NBL/RBL pan-cancer cell-line network |
| `multicohort_cancer` | combined BRCA, NBL, and RBL | multicohort patient-referenced analysis and compatibility profile for explicit pan-cancer targets |
| `pan_cancer` | marker-derived pan-cancer analysis | selectable profile for shared feature-space, ranking, enrichment, and cell-line network products under `results/unsupervised/pan_cancer` |

## Running the full workflow

Run commands from the repository root. Full execution requires prepared inputs, configured paths, adequate storage, and sufficient CPU and memory.

After the required inputs are available, inspect the directed acyclic graph for each profile with a dry-run:

```bash
snakemake -n --use-conda --config pipeline_profile=brca
snakemake -n --use-conda --config pipeline_profile=nbl
snakemake -n --use-conda --config pipeline_profile=rbl
snakemake -n --use-conda --config pipeline_profile=heme
snakemake -n --use-conda --config pipeline_profile=multicohort_cancer
snakemake -n --use-conda --config pipeline_profile=pan_cancer
```

Execute the required profile locally with an explicit core count:

```bash
snakemake --use-conda --cores 8 --config pipeline_profile=brca
snakemake --use-conda --cores 8 --config pipeline_profile=nbl
snakemake --use-conda --cores 8 --config pipeline_profile=rbl
snakemake --use-conda --cores 8 --config pipeline_profile=heme
snakemake --use-conda --cores 8 --config pipeline_profile=multicohort_cancer
snakemake --use-conda --cores 8 --config pipeline_profile=pan_cancer
```

On a SLURM cluster, inspect and adapt [`profiles/slurm/config.yaml`](profiles/slurm/config.yaml) before execution:

```bash
snakemake --profile profiles/slurm --config pipeline_profile=multicohort_cancer
```

Specific rules or output targets may be requested directly. For example:

```bash
snakemake --use-conda --cores 1 test_minimal_inputs --config pipeline_profile=brca
```

## Expected output classes

When the required real inputs are available, the workflow can generate:

- source-aware count and expression objects;
- tumour-purity-filtered patient sample sets;
- unsupervised feature-selection outputs;
- feature-distance representation outputs;
- patient tumour-neighbourhood and p-consensus fraction tables;
- representation-specific patient-referenced graphs;
- majority-threshold consensus-network and union-network outputs;
- resolved cell-line neighbours and node summaries;
- graph-derived DESeq2 marker tables;
- pan-cancer feature-panel tables;
- reciprocal tumour--cell-line ranking metrics;
- pan-cancer cell-line similarity-network outputs;
- Louvain and Leiden community and sensitivity reports;
- functional-enrichment query, background, result, and plotting tables;
- selected figures, manifests, and run reports.

## Reproducibility checks

Repository-level checks that do not require the real datasets include:

```bash
git status --short --branch
snakemake --list --config pipeline_profile=brca
Rscript scripts/validate_minimal_example_inputs.R examples/minimal
snakemake -n --use-conda test_minimal_inputs --config pipeline_profile=brca
snakemake --use-conda --cores 1 test_minimal_inputs --config pipeline_profile=brca
```

Optional rule-graph visualisation can be generated when Graphviz is installed:

```bash
snakemake --rulegraph --config pipeline_profile=brca | dot -Tpdf > rulegraph_brca.pdf
```

Full-analysis releases should provide checksums and manifests alongside the corresponding release assets or archival deposit.

## Interpretation guide

### Patient tumour neighbourhoods and graph outputs

The p-consensus fraction describes how consistently retained tumour samples are assigned near a focal cell-line profile within the configured representations. Patient-referenced edges and resolved neighbours reflect threshold-selected tumour-mediated cell-line similarity under the declared representation rules.

An isolate is therefore a profile without a retained resolved cell-line neighbour. It may still show strong similarity to patient tumours in ranking analyses.

### Reciprocal tumour--cell-line ranking

Reciprocal ranking provides **relative transcriptomic model prioritisation** among the configured candidate cell-line profiles and retained patient tumour samples. A favourable rank is consistent with prioritisation within that analysis space; it is not a universal declaration that the model is suitable for every biological question.

Post hoc cancer-type annotation-concordance metrics describe whether highly ranked tumour--cell-line relationships share the annotated cancer type after the label-free and graph-informed stages. They are evaluation measures, not training objectives.

### Pan-cancer cell-line network

The cell-line-only network describes transcriptomic proximity among DSMZ profiles in the selected pan-cancer feature space. Louvain and Leiden assignments, together with their configured resolution sweeps, assess community structure and its sensitivity to the algorithm and resolution setting. Co-assignment indicates network proximity under those settings; it does not by itself establish a shared biological programme, molecular subtype, or functional equivalence.

### Functional enrichment

Functional-enrichment results provide biological context for graph-derived marker sets. Their interpretation is conditional on the query definition, eligible background, annotation source, database version, and multiple-testing procedure.

## Scope, assumptions, and validation requirements

The workflow provides computational transcriptomic model prioritisation. Its results should be interpreted under the following assumptions:

- The unsupervised patient-referenced workflow is defined by the tumour samples retained after ESTIMATE-based purity filtering. Retained cohort composition, cancer-subtype representation, and tumour microenvironmental admixture may influence tumour-neighbourhood structure, graph resolution, and ranking results.

- Expression-space proximity is affected by RNA-seq source effects, library-preparation differences, count processing, normalisation, batch correction, and systematic differences between patient tumour and DSMZ cell-line samples.

- The inferred patient-referenced graph depends on the configured feature-selection components, expression transformation, distance or similarity metric, the tumour-wise `p_consensus` fractions, and the cross-representation edge-recurrence rule. These choices may affect graph topology, resolved-neighbour status, isolate or anchor assignment, and final rank ordering.

- Tumour–cell-line ranking was performed in a graph-derived pan-cancer marker feature space constructed from recurrent and empirically retained isolate- and anchor-associated DESeq2 markers, rather than across the full transcriptome. This reduces, but does not remove, high-dimensional nearest-neighbour instability and distance-concentration effects; rankings should therefore be interpreted cautiously when the highest-ranked DSMZ models have similar Spearman similarity scores.

- Reciprocal ranking prioritises DSMZ cell-line models relative to the configured candidate set and retained patient-reference cohorts. It does not establish universal model validity, causal equivalence, or suitability for every experimental endpoint.

- Functional-enrichment findings depend on the retained marker lists, enrichment-query construction, eligible gene background, query size, annotation coverage, database version, and multiple-testing correction.

- The computational classification and model prioritisation performed by the workflow require independent biological, perturbational, or pharmacological validation before experimental or translational conclusions are drawn.

## Citation and contact

Repository: [classification_of_dsmz_celllines](https://github.com/EltonUgbogu/classification_of_dsmz_celllines)
