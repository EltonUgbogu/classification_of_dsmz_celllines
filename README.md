# Classification of DSMZ Cell Lines to Cancer Types Based on Transcriptomic Analysis

Cancer cell lines are standard laboratory models, however matching a cell line to a disease label does not guarantee that its transcriptome resembles patient tumours from that disease. This project asks a more specific question: in gene-expression space, how close are DSMZ (German Collection of Microorganisms and Cell Cultures) cancer cell lines to patient tumour RNA-seq cohorts, and does that proximity agree with the expected cancer context? The analysis constructs feature-distance representations from expression data and evaluates cancer-type agreement afterwards. The resulting scores describe feature-space transcriptomic proximity and require biological or experimental follow-up.

This repository provides the Snakemake workflow used for the analysis. It combines tumour-purity filtering, source-aware RNA-seq preprocessing, unsupervised feature selection, feature-distance representations, adaptive tumour-neighbourhood inference, patient-referenced graph resolution, graph-derived marker selection, pan-cancer feature-panel construction, and reciprocal tumour-cell-line ranking across breast cancer (BRCA), neuroblastoma (NBL), retinoblastoma (RBL), and a configured haematological malignancy (HEME) profile.

Curated outputs are organised under [`supplementary_data/`](supplementary_data/); the file-level output catalogue belongs in [`supplementary_data/README.md`](supplementary_data/README.md).

**Navigation:** [Key terms](#key-terms) | [Scientific rationale](#scientific-rationale) | [Method overview](#method-overview) | [Repository structure](#repository-structure) | [Installation](#installation) | [Data availability](#data-availability) | [Configuration](#configuration) | [Running the workflow](#running-the-workflow) | [Outputs](#outputs) | [Interpretation notes](#interpretation-notes) | [Limitations](#limitations) | [Citation and contact](#citation-and-contact)

## Key terms

- **Feature-distance representation:** a selected gene set paired with an expression space and distance or similarity metric.
- **Tumour neighbourhood:** the adaptive set of patient tumour samples assigned near a focal cell line within a clustering representation.
- **p-consensus tumour-neighbourhood support:** the fraction of available clustering outputs within a feature-distance representation that assign a tumour sample to a cell-line neighbourhood. Values range from 0 to 1; configured strong-support thresholds identify stable neighbourhood assignments.
- **Patient-referenced graph:** a support-filtered cell-line graph derived from shared patient tumour-neighbourhood evidence, not a direct biological interaction network.
- **Resolved neighbour:** a cell-line neighbour retained after intersecting cohort-level representation support with cell-line-specific representation support.
- **Global-best representation:** the cohort-level representation selected from p-consensus tumour-neighbourhood metrics across cell lines.
- **Local-best representation:** the cell-line-specific representation set selected from co-optimal or best-supported evidence for an individual cell line.
- **Isolate:** a cell line with degree zero in the final resolved graph. This means no stable graph-resolved cell-line neighbour was retained under the threshold-constrained graph-resolution rules, not absence of tumour similarity.
- **Component anchor:** a representative cell line selected within a non-isolate resolved component, including most-connected and bridge-like graph anchors used for marker-contrast construction.
- **Marker-source class:** the graph-derived source category recording whether retained marker evidence came from isolate-derived or anchor/component-derived contrasts. Membership can be non-exclusive.

## Scientific rationale

Cancer cell lines are widely used to study tumour biology and test experimental hypotheses. Label-level cancer matching is insufficient because patient tumours and cell lines can show intra-cancer-type transcriptomic heterogeneity.

Patient-referenced tumour-neighbourhood evidence gives the workflow a structured comparison space. Instead of relying on label-only matching, the analysis evaluates whether cell lines and patient tumours appear near one another across unsupervised feature-distance representations.

Graph resolution aggregates stable shared tumour-neighbourhood support into a threshold-constrained graph-resolution layer. Graph-derived marker contrasts connect graph structure to pan-cancer feature-panel construction. Reciprocal ranking evaluates post hoc cancer-type agreement in the configured marker-derived feature space.

## Method overview

The workflow proceeds through the following stages:

1. Builds source-aware RNA-seq count and metadata objects for patient tumours and DSMZ cell lines.
2. Applies tumour-purity filtering to patient tumour samples.
3. Constructs filtered, source-aware corrected expression objects and variance-stabilised transformation (VST) matrices.
4. Derives unsupervised feature-selection outputs from variance, distributional, principal-component, and network-based criteria.
5. Constructs feature-distance representations from selected genes, expression spaces, and distance or similarity metrics.
6. Computes adaptive tumour-neighbourhood evidence for cell lines and patient tumours.
7. Aggregates p-consensus tumour-neighbourhood values within each feature-distance representation.
8. Performs patient-referenced graph construction from shared tumour-neighbourhood evidence.
9. Resolves cell-line neighbours through global-local representation intersection.
10. Computes DESeq2 graph-derived marker contrasts for isolates, component anchors, and graph components.
11. Derives the ranked marker-source pan-cancer feature panel.
12. Performs reciprocal tumour-cell-line similarity ranking in the configured feature space.
13. Constructs the pan-cancer cell-line-only k-nearest-neighbour similarity network and evaluates Louvain and Leiden community structure.
14. Evaluates enrichment outputs and quality-control visualisations.

## Repository structure

```text
.
├── README.md                         # project overview
├── Snakefile                         # main Snakemake workflow
├── config/                           # workflow configuration and study-design manifest
├── rules/                            # included Snakemake rule modules
├── scripts/                          # R, Python, and shell scripts called by rules
├── R/                                # shared R helper functions
├── envs/                             # Conda environment specifications
├── resources/                        # small curated resources, such as gene lists and symbol maps
├── preprocessing_and_quality_control/ # upstream preprocessing workflows and environment specifications
├── profiles/                         # Snakemake profile configuration
├── docs/                             # project documentation
├── reports/                          # generated or local reporting artifacts
├── validation/                       # validation scripts and run records
├── supplementary_data/               # curated supplementary data
├── figures/                          # selected curated figures
├── data/                             # external input data; excluded from GitHub
├── results/                          # generated workflow outputs; excluded by default except curated exports
├── logs/                             # runtime logs; excluded from GitHub
├── reference_data/                   # external reference and index resources; excluded from GitHub
└── .snakemake/                       # Snakemake runtime state; excluded from GitHub
```

The repository excludes raw data, large generated outputs, runtime state, and local environment artifacts through ignore policies. Selected curated outputs may be versioned when they support reproducibility or release documentation.

## Quick start

No bundled toy dataset is currently provided. Full execution requires prepared expression/count inputs, reference resources, and profile-specific configuration. Repository-level checks, configuration inspection, and Snakemake dry-runs can still be used to validate the workflow structure after cloning.

Typical first checks are:

```bash
git status
snakemake --list --config pipeline_profile=brca
snakemake -n --use-conda --config pipeline_profile=brca
```

Runtime and memory requirements depend on the selected profile, available input matrices, and enabled optional stages. Lightweight checks such as `snakemake -n` and `snakemake --list` can be run locally. Full RNA-seq preprocessing, clustering grids, graph construction, and enrichment stages are intended for a configured execution environment with sufficient CPU, memory, and storage.

Verified execution-profile settings include:

| profile file | configured setting |
| --- | --- |
| `profiles/local/config.yaml` | `cores: 8`; `use-conda: true` |
| `profiles/slurm/config.yaml` | `jobs: 200`; default `mem_mb=4000`, `runtime=60`, `threads=1`; rule-specific overrides for selected pan-cancer rules |

## Installation

Clone the repository with HTTPS:

```bash
git clone https://github.com/EltonUgbogu/classification_of_dsmz_celllines.git
cd classification_of_dsmz_celllines
git status
```

An SSH clone is also valid when GitHub SSH keys are configured:

```bash
git clone git@github.com:EltonUgbogu/classification_of_dsmz_celllines.git
cd classification_of_dsmz_celllines
git status
```

The driver Snakemake environment is defined by [`envs/smk.yaml`](envs/smk.yaml) and uses the Conda environment name `smk`:

```bash
conda env create -f envs/smk.yaml
conda activate smk
```

The workflow uses per-rule Conda environments through `--use-conda`. Main rule environments include:

- [`envs/tcga-r-env.yaml`](envs/tcga-r-env.yaml): R analysis rules, including DESeq2, clustering, ranking, enrichment, and plotting.
- [`envs/python-graph-env.yaml`](envs/python-graph-env.yaml): Python graph construction and plotting utilities.
- [`envs/tumour_nh_qc.yaml`](envs/tumour_nh_qc.yaml): tumour-neighbourhood quality-control and UMAP (Uniform Manifold Approximation and Projection) rules.

Additional upstream preprocessing environments are stored under [`preprocessing_and_quality_control/envs/`](preprocessing_and_quality_control/envs/).

## Data availability

Raw FASTQ, BAM/CRAM, reference genome, index, and large intermediate files are not stored in GitHub. Controlled-access or bulky data must be obtained from the original repositories or made available in the execution environment before running the workflow.

Source abbreviations used by the project include TCGA (The Cancer Genome Atlas), GDC (Genomic Data Commons), GEO (Gene Expression Omnibus), SRA (Sequence Read Archive), ENA (European Nucleotide Archive), and TARGET (Therapeutically Applicable Research to Generate Effective Treatments). Curated data-description files record transformed inputs, retained tumour identifiers, cell-line metadata, unavailable sources, and checksum information without redistributing restricted raw sequencing data.

| Repository item | GitHub policy | Reason |
| --- | --- | --- |
| `Snakefile`, `config/`, `rules/`, `scripts/`, `R/` | tracked | workflow source and reproducible method definition |
| `envs/*.yaml` | tracked | portable Conda environment specifications |
| Raw FASTQ/BAM/CRAM/SRA files | excluded | controlled-access, bulky, or reconstructable from original repositories |
| Reference genomes, indexes, and large annotations | excluded | bulky external resources |
| `data/` and full `results/` trees | excluded by default | external inputs, intermediates, and generated outputs |
| Selected final tables and supplementary data | tracked when curated | reproducibility and release documentation |
| Selected final figures | tracked when curated | reporting outputs |
| `.snakemake/`, `logs/`, `benchmarks/`, Conda payloads | excluded | runtime state and local execution artifacts |

The archived legacy feature-panel analysis is retained only as provenance material under [`supplementary_data/archive/stale_171_gene_panel/`](supplementary_data/archive/stale_171_gene_panel/) and is not part of the current active feature space.

## Configuration

The main workflow configuration is [`config/config.yaml`](config/config.yaml). The study-design manifest is [`config/study_design.yaml`](config/study_design.yaml). Paths in the configuration are relative to the repository root unless explicitly absolute.

Before the first run, configure the local data locations:

- Main workflow inputs are under `profiles.<profile>.paths` in `config/config.yaml`. Relative values resolve from the repository root; absolute values may point to external storage.
- NBL and RBL preprocessing roots and reference files are configured in `preprocessing_and_quality_control/nbl/config/config.yaml` and `preprocessing_and_quality_control/rbl/config/config.yaml`. Set `data_root` and the `genome.star_index`, `genome.gtf`, and `genome.fasta` values for the execution environment.
- Standalone download and preprocessing helpers accept `BRCA_DATA_ROOT`, `NBL_DATA_ROOT`, `RBL_DATA_ROOT`, `TARGET_NBL_DIR`, `REFERENCE_ROOT`, and `STAR_INDEX_DIR` environment overrides. `CONDA_SH_PATH` may identify a non-standard Conda installation.
- `external_inputs.retained_source_contrasts` is optional. Set it to the local TSV path only when rebuilding the corresponding g:Profiler query sets; `GPROFILER_RETAINED_SOURCE_CONTRASTS` is the equivalent environment override.

Repository-relative defaults preserve the existing layout when data are stored under `data/`. Large data and references may instead remain outside the checkout by using absolute values in a local configuration override. Machine-specific paths should not be committed.

`pipeline_profile` selects the configured analysis profile:

| profile | biological scope | notes |
| --- | --- | --- |
| `brca` | breast cancer | full unsupervised grid plus PAM50 directions |
| `nbl` | neuroblastoma | full unsupervised grid without PAM50 |
| `rbl` | retinoblastoma | full unsupervised grid without PAM50 |
| `heme` | haematological malignancy profiles | configured for cell-line network analyses; not part of the BRCA/NBL/RBL active pan-cancer feature panel |
| `multicohort_cancer` | combined BRCA, NBL, and RBL | patient-referenced multicohort context for tumour-neighbourhood and graph outputs |
| `pan_cancer` | marker-derived pan-cancer analysis layer | builds pan-cancer feature-space, expression-matrix, tumour/cell-line mapping, ranking, enrichment-query, graph, and community outputs under `results/unsupervised/pan_cancer/` |

`pan_cancer` is both a selectable Snakemake profile and the output namespace for the marker-derived pan-cancer analysis layer. `multicohort_cancer` and `pan_cancer` are not interchangeable; explicit legacy pan-cancer targets may still resolve through `multicohort_cancer` for compatibility, but new pan-cancer commands should use `pipeline_profile=pan_cancer`.

## Running the workflow

Commands are run from the repository root. Dry-runs build the directed acyclic graph without executing jobs:

```bash
snakemake -n --use-conda --config pipeline_profile=brca
snakemake -n --use-conda --config pipeline_profile=nbl
snakemake -n --use-conda --config pipeline_profile=rbl
snakemake -n --use-conda --config pipeline_profile=multicohort_cancer
snakemake -n --use-conda --config pipeline_profile=pan_cancer
snakemake -n --use-conda --config pipeline_profile=heme
```

Execution uses the same profile selection with an explicit core count:

```bash
snakemake --use-conda --cores 8 --config pipeline_profile=brca
snakemake --use-conda --cores 8 --config pipeline_profile=nbl
snakemake --use-conda --cores 8 --config pipeline_profile=rbl
snakemake --use-conda --cores 8 --config pipeline_profile=multicohort_cancer
snakemake --use-conda --cores 8 --config pipeline_profile=pan_cancer
snakemake --use-conda --cores 8 --config pipeline_profile=heme
```

Specific rules or exact file targets can be dry-run or executed without requesting the default target. A verified rule-level example for the pan-cancer two-panel figure is:

```bash
snakemake -n --use-conda plot_pan_cancer_cell_line_two_panel --config pipeline_profile=pan_cancer
snakemake --use-conda --cores 4 plot_pan_cancer_cell_line_two_panel --config pipeline_profile=pan_cancer
```

A verified file target for the same rule is:

```bash
TARGET_FILE=figures/Fig_pan_cancer_cell_line_similarity_network_lineage_community.pdf
snakemake -n --use-conda "$TARGET_FILE" --config pipeline_profile=pan_cancer
snakemake --use-conda --cores 4 "$TARGET_FILE" --config pipeline_profile=pan_cancer
```

## Outputs

### Current pan-cancer feature panel

The current pan-cancer feature panel is graph-informed and DESeq2 marker-derived:

- method: `graph_derived_pan_cancer_feature_selection_v1_revised`
- selection route: recurrent direct retention plus all-three acceptance for singleton and non-recurrent candidates
- feature table: [`supplementary_data/feature_space/pan_cancer_features.tsv`](supplementary_data/feature_space/pan_cancer_features.tsv)
- clean Ensembl list: [`supplementary_data/feature_space/pan_cancer_features_clean.txt`](supplementary_data/feature_space/pan_cancer_features_clean.txt)

Marker-source-class membership records whether a retained marker had anchor-derived or isolate-derived graph-derived marker evidence within a cancer type. Membership is non-exclusive.

### Curated supplementary data

See [`supplementary_data/README.md`](supplementary_data/README.md) for file-level descriptions. The top-level categories are:

| directory | contents |
| --- | --- |
| `supplementary_data/feature_space/` | current pan-cancer feature panel and marker-source annotations |
| `supplementary_data/ranking/` | tumour-to-cell-line and cell-line-to-tumour ranking metrics |
| `supplementary_data/tumour_neighbourhood/` | patient-neighbourhood and resolved-neighbour outputs |
| `supplementary_data/pan_cancer_similarity/` | pan-cancer cell-line network and community outputs |
| `supplementary_data/model_prioritisation/` | model-prioritisation inputs and scores |
| `supplementary_data/enrichment/` | enrichment query, background, and result tables |
| `supplementary_data/data_description/` | input descriptions, retained samples, exclusions, and checksums |
| `supplementary_data/manifests/` | file manifests, checksums, and figure provenance |
| `supplementary_data/archive/` | provenance-only stale or historical outputs |

Selected curated figures are stored under [`figures/`](figures/), and figure provenance is recorded in [`supplementary_data/manifests/figure_manifest.tsv`](supplementary_data/manifests/figure_manifest.tsv).

## Reproducibility checks

Common repository checks include:

```bash
git status
snakemake --list --config pipeline_profile=brca
snakemake -n --use-conda --config pipeline_profile=brca
sha256sum -c supplementary_data/manifests/checksums_sha256.tsv
```

Optional rule-graph visualisation can be generated when Graphviz is installed:

```bash
snakemake --rulegraph --config pipeline_profile=brca | dot -Tpdf > rulegraph_brca.pdf
```

Generated workflow outputs are intentionally excluded from GitHub by default. A clean Git state does not imply that all generated outputs are present.

## Interpretation notes

High tumour-cell-line similarity indicates rank-based transcriptomic proximity in the configured feature space. It does not establish model relevance, causal similarity, or functional substitutability.

Resolved neighbours indicate stable patient-referenced support after applying the configured global-best and local-best representation rules. Isolates indicate lack of retained graph-resolved neighbours under the threshold-constrained graph-resolution rules; they do not indicate absence of tumour similarity.

The patient-referenced graph is a support-filtered graph representation of shared tumour-neighbourhood evidence, not a direct biological interaction network.

Marker-derived features reflect graph-informed contrast selection and filtering. Enrichment terms provide interpretability and must be interpreted with respect to the configured query sets, background definitions, source databases, and multiple-testing correction.

Ranking metrics assess agreement with annotated cancer type after unsupervised analysis. They are post hoc agreement metrics rather than training objectives.

## Limitations

- RNA-seq source effects and sample-type confounding can affect transcriptomic proximity.
- Feature-set dependence and distance-metric dependence influence graph topology and rank ordering.
- Threshold dependence affects p-consensus support calls and the graph-resolved neighbourhood.
- High-dimensional nearest-neighbour instability can arise from preprocessing, scaling, and sparse rank margins.
- Small rank-margin sensitivity can alter model-prioritisation order when candidates have similar scores.
- Tumour-purity uncertainty and cohort-composition dependence affect the patient-reference space.
- Functional-enrichment signal can be limited by query size, background definition, annotation coverage, and multiple-testing correction.
- Cell-line model prioritisation requires independent experimental, perturbational, or pharmacological validation beyond transcriptomic similarity.
- Controlled-access and bulky data are not redistributed, so full reproduction requires access to the original data sources or equivalent prepared inputs.

## Citation and contact

Repository: [classification_of_dsmz_celllines](https://github.com/EltonUgbogu/classification_of_dsmz_celllines)

Citation metadata and release DOI information will be added when a versioned software release is archived.
