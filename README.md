# Classification of DSMZ Cell Lines to Cancer Types Based on Transcriptomic Analysis

This repository contains the Snakemake workflow used to prioritise DSMZ cancer cell-line models by their transcriptomic proximity to patient tumour RNA-seq cohorts. The central question is not whether a cell line carries the same cancer-type label as a patient cohort, but how closely its expression profile relates to retained patient tumour samples in a configured, graph-informed feature space.

The patient-referenced analyses use breast cancer (BRCA), neuroblastoma (NBL), and retinoblastoma (RBL) cohorts. An auxiliary haematological malignancy (HEME) profile is included for cell-line-only network context, but it is not a source cohort for the active BRCA/NBL/RBL pan-cancer marker panel or the corresponding patient-referenced ranking analyses.

The workflow constructs feature-distance representations, infers patient tumour neighbourhoods, builds multi-representation consensus networks, resolves cell-line neighbours, derives graph-informed marker evidence, constructs a revised Version 1 pan-cancer feature panel, and performs reciprocal tumour--cell-line ranking. These outputs provide relative transcriptomic model prioritisation within the configured data and feature space.

The repository is maintained as a **code-focused workflow repository**. Raw sequencing data, large expression matrices, reference indexes, complete result trees, runtime logs, and most generated outputs are intentionally excluded from the main Git history. A tiny synthetic example is included so that repository structure and input schemas can be tested without access to the real datasets.

**Navigation:** [Key terms](#key-terms) | [Scientific aim](#scientific-aim) | [Workflow overview](#workflow-overview) | [Current analytical design](#current-analytical-design) | [Repository structure](#repository-structure) | [Quick start](#quick-start) | [Installation](#installation) | [Synthetic minimal example](#synthetic-minimal-example) | [Data and outputs](#data-and-outputs) | [Configuration](#configuration) | [Running the full workflow](#running-the-full-workflow) | [Expected output classes](#expected-output-classes) | [Reproducibility checks](#reproducibility-checks) | [Interpretation guide](#interpretation-guide) | [Scope assumptions and validation requirements](#scope-assumptions-and-validation-requirements) | [Citation and contact](#citation-and-contact)

## Key terms

- **Feature-distance representation:** a selected gene set, expression space, and distance or similarity metric considered together as one representation.
- **Tumour neighbourhood:** the adaptive set of retained patient tumour samples assigned near a focal cell-line profile within a feature-distance representation.
- **p-consensus tumour-neighbourhood support:** the fraction of available tumour-neighbourhood clustering outputs within a fixed feature-distance representation that support assignment of a tumour sample to a focal cell-line neighbourhood.
- **Patient-referenced graph:** a cell-line graph derived from shared patient tumour-neighbourhood evidence. It is not a direct molecular or biological interaction network.
- **Multi-representation consensus network:** a graph layer that retains cell-line relationships recurring across the configured feature-distance representations under the selected cross-representation support rule.
- **Resolved graph:** the graph obtained after combining cohort-level representation selection with cell-line-specific tied-best representation support.
- **Resolved neighbour:** a cell-line neighbour retained in the resolved graph.
- **Global-best representation:** the cohort-level feature-distance representation selected from the configured tumour-neighbourhood performance criteria.
- **Local-best representation set:** the tied-best or co-optimal representation set selected for an individual cell-line profile.
- **Isolate:** a cell-line profile with degree zero in the resolved graph. This means no stable cell-line neighbour was retained under the configured graph-resolution rules; it does not mean that the cell line lacks similarity to patient tumours.
- **Component anchor:** a profile selected from a non-isolate resolved component through highest-degree or highest unnormalised-betweenness criteria, with deterministic tie handling.
- **Marker-evidence stratum:** the graph-derived evidence category used during marker aggregation: anchor-associated or isolate-associated marker evidence.
- **Pan-cancer feature panel:** the graph-informed marker-derived feature set used for pan-cancer cell-line similarity analysis and reciprocal tumour--cell-line ranking.
- **Model prioritisation:** relative ordering of candidate cell-line models according to patient-referenced transcriptomic proximity in the configured feature space.

## Scientific aim

Cancer cell lines are widely used as laboratory models, but matching a cell line to a cancer-type label does not establish how well its transcriptome represents patient tumours. Patient tumours and cell-line profiles can vary substantially within the same cancer type, and that heterogeneity can change which models appear most representative.

This workflow therefore evaluates DSMZ cell-line profiles in patient-referenced expression space. It constructs unsupervised feature-distance representations and asks whether cell-line profiles and retained patient tumour samples remain close across those representations. Graph resolution summarises stable shared tumour-neighbourhood evidence, while graph-derived marker contrasts connect resolved graph structure to pan-cancer feature construction.

Reciprocal ranking then evaluates the relative position of tumour and cell-line profiles in the marker-derived feature space. Cancer-type agreement is assessed after the unsupervised and graph-informed stages; it is not used as a supervised training objective.

## Workflow overview

At a high level, the workflow performs the following stages:

1. Prepare source-aware RNA-seq count and expression objects.
2. Apply tumour-purity filtering to patient tumour cohorts where configured.
3. Generate variance-stabilised expression matrices for patient tumour and DSMZ cell-line profiles.
4. Derive unsupervised feature sets using variance-, distribution-, principal-component-, and network-based criteria.
5. Combine feature sets, expression spaces, and distance or similarity metrics into feature-distance representations.
6. Infer adaptive patient tumour neighbourhoods for focal cell-line profiles.
7. Aggregate p-consensus tumour-neighbourhood support within each representation.
8. Construct representation-specific patient-referenced cell-line graphs.
9. Retain recurrent cell-line relationships in a multi-representation consensus network.
10. Resolve cell-line neighbours using cohort-level representation selection and cell-line-specific tied-best representation support.
11. Identify graph-derived isolate and anchor profiles for differential-expression analysis.
12. Run DESeq2 isolate and anchor contrasts using cancer-type-specific cell-line count matrices.
13. Aggregate contrast-level marker evidence into the revised Version 1 pan-cancer feature panel.
14. Perform reciprocal tumour--cell-line ranking in the configured pan-cancer feature space.
15. Construct the pan-cancer cell-line-only \(k\)-nearest-neighbour similarity network.
16. Evaluate Louvain community structure and Leiden sensitivity results.
17. Construct and evaluate functional-enrichment query and background sets for marker interpretation.

## Current analytical design

### Patient-referenced graph resolution

Tumour-neighbourhood evidence is calculated separately within each feature-distance representation. Strong within-representation p-consensus evidence is used to construct representation-specific cell-line graphs. Cell-line relationships recurring across representations form the multi-representation consensus network.

The resolved graph uses a cohort-level selected representation together with each cell line's tied-best representation set. A resolved neighbour therefore represents stable support under the configured graph-resolution rules, rather than a direct claim of biological interaction or functional equivalence.

### Graph-derived DESeq2 contrasts

Differential-expression analysis uses cell-line-only raw integer count matrices within each cancer type. Patient tumour samples are not included in the DESeq2 contrasts.

| Contrast type | Focal profile | Reference profiles | Marker filter |
| --- | --- | --- | --- |
| isolate | one degree-zero resolved-graph profile | all other same-cancer-type cell-line profiles | adjusted \(p \leq 0.01\), \(|\log_2\mathrm{FC}| \geq 1.5\), `baseMean >= 10`, maximum 50 markers |
| anchor | one graph-selected anchor profile | same-cancer-type profiles outside the focal resolved component | adjusted \(p \leq 0.05\), \(|\log_2\mathrm{FC}| \geq 1.0\), `baseMean >= 10`, maximum 200 markers |

For anchor contrasts, non-anchor profiles from the focal component are excluded from the reference set. Marker-list caps are maxima, not required list sizes.

### Revised Version 1 pan-cancer feature selection

Marker recurrence is counted separately within each cancer type and marker-evidence stratum.

- Genes observed in at least two filtered contrast-level marker lists are included directly as recurrent evidence.
- Genes observed once are separated into singleton candidates or non-recurrent candidates according to the number of available contrast lists in the relevant cancer type and marker-evidence stratum.
- Candidate genes are accepted only when they pass all three empirical evidence criteria:
  - minimum adjusted \(p\)-value at or below the candidate-pool first quartile;
  - median absolute \(\log_2\) fold change at or above the candidate-pool third quartile;
  - median `baseMean` at or above the candidate-pool median.
- Empirical thresholds are estimated within cancer type, marker-evidence stratum, and candidate-pool type.
- The final panel is the union of recurrent genes, accepted singleton candidates, and accepted non-recurrent candidates.

This revised Version 1 procedure replaces the obsolete minimum-panel-size and sequential top-up rules used in earlier development states.

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
├── preprocessing_and_quality_control/ # upstream preprocessing workflows and configuration
├── profiles/                          # local and SLURM Snakemake profiles
├── docs/                              # project documentation
├── examples/minimal/                  # tiny synthetic schema-compatible inputs
└── tests/                             # repository checks and test-related files
```

The following paths contain local data, references, generated outputs, provenance, or runtime state and are excluded from the public code-only main branch unless deliberately packaged as release assets:

```text
data/
reference_data/
results/
reports/
figures/
supplementary_data/
logs/
benchmarks/
archive/
derived_results/
sync_logs/
.snakemake/
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

The synthetic target validates the expected matrix and metadata schemas and confirms that the rule is included in the main workflow. It does not execute the full biological analysis or produce scientifically interpretable results.

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

No machine-specific Conda bootstrap script is required. Install Conda, Mamba, or an equivalent compatible environment manager separately, then construct the tracked environments from the YAML definitions.

## Synthetic minimal example

The repository includes tiny synthetic BRCA, NBL, and RBL inputs under [`examples/minimal/`](examples/minimal/). They contain:

- fake Ensembl-like gene identifiers;
- fake cell-line and tumour sample identifiers;
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

Large or controlled-access data must be obtained from their original repositories or supplied as equivalent prepared inputs. Source abbreviations used by the project include TCGA, GDC, GEO, SRA, ENA, and TARGET.

| Repository item | Public-main policy | Reason |
| --- | --- | --- |
| `Snakefile`, `config/`, `rules/`, `scripts/`, `R/` | tracked | workflow source and reproducible method definition |
| `envs/*.yaml` | tracked | portable Conda environment specifications |
| `examples/minimal/` | tracked | tiny synthetic schema-compatible validation inputs |
| `resources/` | tracked when small and curated | lightweight resources required by code |
| raw sequencing files | excluded | bulky, controlled-access, or reconstructable from original repositories |
| reference genomes and indexes | excluded | bulky external resources |
| `data/`, `reference_data/` | excluded | execution-environment inputs and references |
| `results/`, `logs/`, `.snakemake/`, `benchmarks/` | excluded | generated runtime state |
| `figures/`, `supplementary_data/`, `reports/` | excluded from code-only main | generated outputs, release assets, or provenance material |

Curated figures, result tables, checksums, and file manifests should be distributed through versioned release assets or an external archival repository rather than restored as large generated trees in the main branch.

Historical outputs from obsolete feature-panel definitions are provenance material only and must not be presented as the active analysis state.

## Configuration

The main workflow configuration is [`config/config.yaml`](config/config.yaml). The study design is recorded in [`config/study_design.yaml`](config/study_design.yaml). Paths are interpreted relative to the repository root unless they are explicitly absolute.

Before a full run:

- configure profile-specific inputs under `profiles.<profile>.paths` in `config/config.yaml`;
- make the required count, expression, metadata, and reference inputs available;
- configure NBL and RBL preprocessing paths in their respective preprocessing configuration files;
- use absolute paths or an untracked local configuration override for machine-specific storage locations;
- do not commit workstation- or cluster-specific paths.

`pipeline_profile` selects the analysis profile:

| profile | scope | notes |
| --- | --- | --- |
| `brca` | breast cancer | patient-referenced BRCA workflow |
| `nbl` | neuroblastoma | patient-referenced NBL workflow |
| `rbl` | retinoblastoma | patient-referenced RBL workflow |
| `heme` | haematological malignancy cell-line profiles | auxiliary cell-line context; not a BRCA/NBL/RBL patient-reference cohort |
| `multicohort_cancer` | combined BRCA, NBL, and RBL | shared multicohort and pan-cancer stages |

`pan_cancer` is an output namespace, not a selectable Snakemake profile.

## Running the full workflow

Run commands from the repository root. Full execution requires prepared inputs, configured paths, adequate storage, and sufficient CPU and memory.

After the required inputs are available, inspect the directed acyclic graph with a dry-run:

```bash
snakemake -n --use-conda --config pipeline_profile=brca
snakemake -n --use-conda --config pipeline_profile=nbl
snakemake -n --use-conda --config pipeline_profile=rbl
snakemake -n --use-conda --config pipeline_profile=multicohort_cancer
snakemake -n --use-conda --config pipeline_profile=heme
```

Execute locally with an explicit core count:

```bash
snakemake --use-conda --cores 8 --config pipeline_profile=brca
snakemake --use-conda --cores 8 --config pipeline_profile=nbl
snakemake --use-conda --cores 8 --config pipeline_profile=rbl
snakemake --use-conda --cores 8 --config pipeline_profile=multicohort_cancer
snakemake --use-conda --cores 8 --config pipeline_profile=heme
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
- patient tumour-neighbourhood and p-consensus support tables;
- representation-specific patient-referenced graphs;
- multi-representation consensus-network outputs;
- resolved cell-line neighbours and node summaries;
- graph-derived DESeq2 marker tables;
- revised Version 1 pan-cancer feature-panel tables;
- reciprocal tumour--cell-line ranking metrics;
- pan-cancer cell-line similarity-network outputs;
- Louvain and Leiden community and sensitivity summaries;
- functional-enrichment query, background, result, and plotting tables;
- selected figures, manifests, and run reports.

Generated outputs are excluded from GitHub by default. A clean Git working tree means the tracked code tree is clean; it does not mean all workflow outputs exist locally.

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

Patient tumour-neighbourhood support describes how consistently retained tumour samples are assigned near a focal cell-line profile within the configured representations. Patient-referenced edges and resolved neighbours summarise shared tumour-neighbourhood evidence under the selected thresholds and representation rules.

An isolate is therefore a profile without a retained resolved cell-line neighbour. It may still show strong similarity to patient tumours in ranking analyses.

### Reciprocal tumour--cell-line ranking

Reciprocal ranking provides **relative transcriptomic model prioritisation** among the configured candidate cell-line profiles and retained patient tumour samples. A favourable rank supports prioritisation within that analysis space; it is not a universal declaration that the model is suitable for every biological question.

Cancer-type agreement metrics describe whether highly ranked tumour--cell-line relationships agree with the annotated cancer type after the unsupervised and graph-informed stages. They are evaluation measures, not training objectives.

### Pan-cancer cell-line network

The cell-line-only network describes transcriptomic proximity among DSMZ profiles in the selected pan-cancer feature space. Louvain or Leiden co-assignment indicates network proximity under the chosen graph construction and resolution settings; it does not by itself establish a shared biological programme, molecular subtype, or functional equivalence.

### Functional enrichment

Functional-enrichment results provide biological context for graph-derived marker sets. Their interpretation is conditional on the query definition, eligible background, annotation source, database version, and multiple-testing procedure.

## Scope, assumptions, and validation requirements

The workflow provides computational, transcriptomic model prioritisation. Its results should be interpreted under the following assumptions and boundaries:

- Patient-referenced proximity is conditional on the tumour samples retained after purity filtering. Cohort composition, cancer-subtype representation, source distribution, and sample availability can influence tumour neighbourhoods, graph resolution, and ranking results.
- RNA-seq source effects, library preparation, preprocessing choices, and differences between patient tumours and cell-line samples can influence expression-space proximity.
- Feature-set choice, expression transformation, distance metric, p-consensus threshold, and cross-representation support rule can alter graph topology, resolved-neighbour status, and rank ordering.
- Restricting the analysis to a controlled marker-derived feature space and using rank-based comparisons mitigates, but does not eliminate, high-dimensional nearest-neighbour and distance-concentration effects. Interpretation should be cautious when rank margins are small.
- Reciprocal ranking prioritises models relative to the configured candidate set and retained patient-reference cohorts. It does not establish universal model validity, causal equivalence, or suitability for every experimental endpoint.
- Study-specific claims about mechanism, perturbation response, drug sensitivity, or therapeutic relevance require independent biological, perturbational, or pharmacological validation. This requirement qualifies downstream claims; it does not negate the computational prioritisation performed by the workflow.
- Functional-enrichment findings depend on query size, candidate selection, eligible background, annotation coverage, database version, and multiple-testing correction.
- Full reproduction requires access to the original data sources or equivalent prepared inputs because controlled-access and bulky datasets are not redistributed.

## Citation and contact

Repository: [classification_of_dsmz_celllines](https://github.com/EltonUgbogu/classification_of_dsmz_celllines)

Until a versioned release DOI is available, cite the repository URL together with the Git commit hash or release tag used for the analysis. Add formal citation metadata and archival DOI information when the workflow is deposited as a versioned software release.
