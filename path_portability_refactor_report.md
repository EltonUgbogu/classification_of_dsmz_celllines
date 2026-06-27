# Path portability refactor report

## Outcome

Project-owned paths now resolve from the checkout rather than requiring the Jarvis or Mac filesystem layout. External data, references, indexes, and downloaded payload roots use repository-relative defaults and support explicit local overrides. Scientific logic, thresholds, feature-selection settings, graph-resolution criteria, marker filters, rule names, and output filenames were not changed. No full workflow was run, no generated results were edited, no files were deleted, and no commit was created.

## Pre-refactor audit

`path_portability_audit.tsv` records 711 path tokens across 172 active source/documentation files. The narrower machine-specific scan found 158 matching lines across 52 files.

| Classification | Entries |
| --- | ---: |
| `project_internal_absolute` | 45 |
| `external_data_absolute` | 96 |
| `generated_output_reference` | 6 |
| `documentation_only` | 18 |
| `uncertain` | 8 |
| `environment_variable_path` | 272 |
| `system_absolute_allowed` | 251 |
| `url_or_accession` | 15 |
| **Total** | **711** |

The audit manifest SHA-256 is `561696bb778b1709c2c1d1481c5b8bc5bfce43f85ea9d4eef0b2335f6cbb7032`.

## Files edited

Fifty-three source, configuration, launcher, and documentation files were edited:

- Root/workflow: `Snakefile`, `config/config.yaml`, `README.md`, `run_all_neighbourhoods.sh`, and `rules/pan_cancer_enrichment_marker_framework.smk`.
- Operational documentation: `docs/AGENT_TASK_PIPELINE_ANALYSIS.md`, `docs/DSMZ_GRAPH_RULE_DEBUG.md`, `docs/PIPELINE_PROFILES_AND_OUTPUT_NAMESPACES.md`, and `envs/CONDA_SETUP.md`.
- Shared preprocessing: `preprocessing_and_quality_control/build_star_index.sh` and `preprocessing_and_quality_control/brca/build_tcga_brca_star_counts_manifest.py`.
- NBL preprocessing: `preprocessing_and_quality_control/nbl/Snakefile`; `config/config.yaml`; `preprocessing_snakefile.sh`; `preprocessing_snakefile_count_only.sh`; `batch_corr_and_normalisation/Snakefile`; `batch_corr_and_normalisation/config/config.yaml`; `batch_corr_and_normalisation/config/config.test.yaml`; `batch_corr_and_normalisation/preprocess_batch_correction_normalisation.sh`; `tumour_purity_analysis/nbl_tumour_purity_analysis.R`; and 12 active files under `preprocessing_and_quality_control/nbl/scripts/` (`compare_target_nbl_id_sets.py`, `download_nbl_primary_pe.sh`, `download_nbl_tumour_sample_srr_ids.py`, `download_target_nbl_star_counts.sh`, `generate_GSE189367_comparison_nbl.py`, `generate_nbl_target_id_comparison_report.py`, `generate_target_nbl_star_manifest.py`, `list_nbl_sample_ids_by_cohort.R`, `merge_star_counts.R`, `merge_target_nbl_counts.R`, `run_download_gse100148_missing.sh`, and `run_download_nbl_srr_ids.sh`).
- RBL preprocessing: `preprocessing_and_quality_control/rbl/Snakefile`; `config/config.yaml`; `preprocessing_count_only_rbl.sh`; `preprocessing_snakefile.sh`; `preprocessing_snakefile_rbl_full_preprocess.sh`; `batch_corr_and_normalisation/preprocess_batch_correction_normalisation.sh`; five active files under `preprocessing_and_quality_control/rbl/scripts/` (`download_rbl_primary_pe.sh`, `download_rbl_tumour_sample_srr_ids.py`, `list_rbl_tumour_sample_ids_by_cohort.R`, `merge_star_counts.R`, and `run_download_rbl_srr_ids.sh`); and `preprocessing_and_quality_control/tumour_purity_analysis/rbl/rbl_tumour_purity_analysis.R`.
- Main analysis scripts: `scripts/build_enrichment_summary_top_terms.R`, `scripts/build_ranked_marker_source_panel_enrichment_queries.py`, `scripts/finalize_ranked_marker_source_panel_reports.py`, `scripts/make_cell_cycle_ensg_from_gmt.R`, `scripts/mycn_gene_presence.R`, `scripts/package_ranked_marker_source_panel_outputs.py`, `scripts/pan_agnostic_clustering.R`, `scripts/plot_ecdf_rank_combined_pub.R`, and `scripts/plot_enrichment_top_terms_heatmap.R`.

## Project-internal paths converted

All 45 audited project-internal absolute entries and all 6 embedded checkout/provenance references were converted. The main `Snakefile` now resolves `REPO_ROOT` from `workflow.basedir`, exposes `repo_path()` and `cfg_path()`, and exports the resolved root to rule processes. Nested preprocessing workflows resolve relative configuration values against the same repository root. Python scripts derive repository defaults from `Path(__file__).resolve()`, R scripts derive them from the `Rscript --file` path, and shell launchers derive them from `BASH_SOURCE`.

Provenance files now record the runtime checkout rather than a fixed Jarvis label. Operational commands and generated-report command examples use repository-relative targets. Intentional `$CONDA_PREFIX/bin/python` uses remain unchanged.

## External paths made configurable

All 96 audited external-data entries were converted to relative configuration values or explicit overrides. Supported controls include:

- YAML: profile input paths in `config/config.yaml`; NBL/RBL `data_root`; `genome.star_index`, `genome.gtf`, and `genome.fasta`; purity and batch-correction output/input paths; and optional `external_inputs.retained_source_contrasts`.
- Environment: `BRCA_DATA_ROOT`, `NBL_DATA_ROOT`, `RBL_DATA_ROOT`, `TARGET_NBL_DIR`, `REFERENCE_ROOT`, `STAR_INDEX_DIR`, `GPROFILER_RETAINED_SOURCE_CONTRASTS`, and `CONDA_SH_PATH`.
- CLI arguments retained by standalone comparison/report scripts.

Relative defaults reproduce the current Jarvis layout when the checkout is `/work/ugbogu/pipeline`, while absolute user overrides can point to external storage. The optional retained-source-contrast input is represented as an empty Snakemake input list when unset.

The generated 379-gene feature list previously caused parse-time failure when `results/` was absent. The workflow now falls back only when necessary to `supplementary_data/feature_space/pan_cancer_features_clean.txt`; the generated and curated files were verified byte-identical (`a72e94f5083248fbe7a7af035a39abd28a5c5c803ceec8edb738cd51f632e84d`). This preserves the configured `expected_genes: auto` criterion.

## Allowed paths retained

The final allowed-path scan found 182 lines containing valid interpreter, temporary, device, environment, or Conda paths. Retained categories include `/usr/bin/env`, `/bin/bash`, `/bin/sh`, `/dev/null`, `/dev/stderr`, `/tmp`, system certificate paths, URLs, `${HOME}`-derived Conda paths, and `$CONDA_PREFIX/bin/python` or `$CONDA_PREFIX/bin/Rscript`.

Three historical traceback lines in `docs/SUPPORT_THRESHOLD_PASS_REPORT.md` retain `/work/ugbogu/pipeline/Snakefile` as an execution record. Existing `.bak` files were found only inside the excluded `.cleanup_quarantine_deep_clean/` quarantine; this refactor created none.

## Uncertain paths not changed

Two debug-only files remain unchanged:

- `preprocessing_and_quality_control/nbl/batch_corr_and_normalisation/check.R`
- `preprocessing_and_quality_control/rbl/batch_corr_and_normalisation/check.R`

Each contains the same Mac Google Drive scratch path. They are not referenced by a Snakefile or launcher and were classified as uncertain rather than silently rewritten. The other six initially uncertain entries were literal user Conda paths and were converted to runtime discovery.

## Validation results

| Check | Result |
| --- | --- |
| Final forbidden-path scan | Five reviewed matches only: two debug `check.R` paths and three historical traceback lines; no active workflow/config/script match |
| Allowed system/environment scan | 182 lines retained |
| Bash syntax | Passed for all edited launchers/scripts |
| Python `py_compile` | Passed for active `scripts/**/*.py` and edited preprocessing Python files; bytecode redirected to `/tmp` |
| R parse | Passed for all edited R files using the existing Snakemake R environment |
| Main profile parse (`--list-rules`) | Passed for BRCA, NBL, RBL, multicohort cancer, and HEME |
| Nested workflow parse | Passed for NBL preprocessing, RBL preprocessing, and NBL batch correction |
| BRCA dry-run | Parsed; stopped at the pre-existing missing `patient_referenced_resolved_cell_line_neighbourhood_graph_full_node_labels.tsv` default target |
| NBL dry-run | Parsed; stopped at the equivalent pre-existing missing default target |
| RBL dry-run | Parsed; stopped at the equivalent pre-existing missing default target |
| HEME dry-run | Parsed; stopped at the equivalent pre-existing missing default target |
| Multicohort dry-run | Passed (status 0) |
| Temporary-root parse | Passed from `/tmp/tmp.jJBzT83229/portable_repo` with `.git`, `.snakemake`, `data`, and `results` excluded |
| Temporary-root BRCA dry-run | Reached the expected missing external `data/brca/brca_vst_joint.rds`; zero references to the original checkout and zero external-root symlinks |
| `git diff --check` | Unavailable: `/work/ugbogu/pipeline` is not a Git worktree (`git` status 129) |

The four default-target dry-run failures are not path-resolution failures: each workflow parsed and listed its rules successfully, then Snakemake reported an existing target with no available producer/input. The temporary-copy failure is expected because the test deliberately excluded required external data.

## Remaining portability risks

- The two debug `check.R` files remain machine-specific until their intended status is confirmed.
- Four profiles have an existing incomplete default-target dependency for the full-node-label TSV; this should be repaired separately from portability work.
- External datasets and references are not distributed with the repository. A clean clone parses, but execution requires their configured local paths.
- Relative `#SBATCH` log destinations are interpreted from the submission directory before shell code starts. Submit Slurm wrappers from the repository root or override `sbatch --output/--error` for site-specific logging.
- The remote directory lacks Git metadata, so this run could not prove change scope with Git or execute `git diff --check`.

## Configuration for a new machine

1. Clone the repository and enter its root.
2. In `config/config.yaml`, set required `profiles.<profile>.paths` inputs. Relative paths resolve from the clone; absolute paths may point to mounted data storage.
3. For upstream NBL/RBL processing, set `data_root` and the three `genome` paths in the corresponding preprocessing config.
4. For standalone helpers, export roots as needed, for example:

   ```bash
   export NBL_DATA_ROOT=/path/to/nbl_data
   export RBL_DATA_ROOT=/path/to/rbl_data
   export REFERENCE_ROOT=/path/to/reference_data
   ```

5. Set `external_inputs.retained_source_contrasts` only if rebuilding the marker-framework enrichment query sets with that optional audit artefact.
6. Activate the driver environment and parse/dry-run the required profile:

   ```bash
   conda activate smk
   snakemake --list-rules --config pipeline_profile=brca
   snakemake -n --config pipeline_profile=brca --cores 1
   ```

The current Jarvis layout continues to resolve to the same `data/`, `results/`, `resources/`, script, environment, and reference locations when launched from the repository checkout.
