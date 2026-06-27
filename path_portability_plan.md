# Path portability policy

## Repository-root resolution

The main `Snakefile` resolves the checkout root once from `workflow.basedir`:

```python
REPO_ROOT = Path(workflow.basedir).resolve()

def repo_path(*parts):
    return str(REPO_ROOT.joinpath(*parts))

def cfg_path(value):
    path = Path(str(value)).expanduser()
    return str(path if path.is_absolute() else REPO_ROOT / path)
```

`BASE` and `PIPE_ROOT` may remain as string aliases of `REPO_ROOT` where existing rules require strings. Nested cohort Snakefiles derive the same repository root from their own `workflow.basedir` and resolve relative configuration values against it. No workflow path depends on the caller's current working directory.

## Project-internal paths

Repository-owned scripts, rules, configuration, environments, resources, results, figures, logs, and supplementary outputs are represented as repository-relative paths in configuration and rule outputs. Absolute paths needed by scripts or shell commands are constructed with `repo_path(...)` or `cfg_path(...)`. Standalone Python and R scripts take paths from arguments, Snakemake, or environment variables; fallback repository paths are resolved from the script location. Shell launchers derive `SCRIPT_DIR` and `REPO_ROOT` before constructing paths, with all path variables quoted.

## External data and reference roots

Raw data, downloaded payloads, reference FASTA/GTF files, aligner indexes, and external audit artefacts are configuration inputs. Relative defaults such as `data/nbl` and `data/reference` resolve against the repository root; absolute local/HPC overrides remain supported through YAML keys, command-line arguments, or documented environment variables. User-specific paths are never the sole default.

## Allowed absolute paths

System interpreter paths (`/usr/bin/env`, `/bin/bash`, `/bin/sh`), device paths (`/dev/null`, `/dev/stderr`), temporary paths under `/tmp`, dynamically discovered system certificate paths, URLs, and intentional `$CONDA_PREFIX/bin/python` uses are allowed. Environment-derived paths such as `${HOME}` are allowed when they are portable runtime discovery rather than a literal user directory.

## Forbidden absolute paths

Active source and operational documentation must not require `/work/ugbogu/pipeline`, `/Users/...`, `/home/<user>/...`, `/mnt/...`, `/Volumes/...`, or any other checkout-specific path for repository-owned files. External resources may use absolute paths only as user-supplied overrides, never as committed machine-specific defaults.

## Files to edit

Edits are limited to manifest entries requiring action: the main `Snakefile`; active scripts under `scripts/`; root launchers; active cohort Snakefiles, YAML files, scripts, and launchers under `preprocessing_and_quality_control/`; `README.md`; operational documentation; and provenance writers that embed the checkout path. Configuration changes alter path representation only, not scientific values.

## Files intentionally left unchanged

Generated results, `.snakemake/`, raw data, archives (including the legacy `archieve` directory), deprecated scripts, quarantine/backup directories, and historical report tracebacks are excluded. Allowed system paths, URLs, accession identifiers, `$CONDA_PREFIX/bin/python`, and portable environment-derived paths are retained. Debug-only `batch_corr_and_normalisation/check.R` files remain unchanged unless they are confirmed as active; their Mac path is recorded as uncertain in the audit. The external `/work/ugbogu/local_audits/...` marker-audit input will be made optional/configurable without changing its current Jarvis fallback behaviour.

## Validation boundary

Validation consists of forbidden-path scans, allowed-system-path checks, syntax checks, Snakemake parsing/dry-runs for configured profiles, a data/results-excluded temporary-copy dry-run, and whitespace checks. The full workflow and environment creation are out of scope.
