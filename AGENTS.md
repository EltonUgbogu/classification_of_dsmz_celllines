# AGENTS.md

## Cursor Cloud specific instructions

### Overview

This is a **bioinformatics research pipeline** for classification of DSMZ cancer cell lines against TCGA primary tumour samples. It consists of:

1. **Snakemake pipeline** (`pipeline/`) — unsupervised clustering and tumour neighbourhood analysis
2. **Expression Transformer** (`dsmz/transformer/`) — PyTorch-based self-supervised gene expression model (GPU-only training)
3. **R analysis scripts** — correlation, preprocessing, QC

### Environment Setup

- **Conda**: Miniforge3 installed at `$HOME/miniforge3`. Always source conda before use:
  ```
  export PATH="$HOME/miniforge3/bin:$PATH"
  source $HOME/miniforge3/etc/profile.d/conda.sh
  ```
- **Snakemake**: Available in the `smk` conda env (`conda activate smk`). Contains R 4.5+ and Snakemake 9.x.
- **Python**: System Python 3.12 with PyTorch, NumPy, tqdm installed globally via pip.
- **No databases or web services** — all data is file-based (RDS, CSV, TSV, NPZ).

### Running the Pipeline

See `pipeline/README.md` and `pipeline/AGENT.md` for detailed pipeline documentation.

- **Entry point**: `pipeline/run_pipeline.sh --profile <name> --executor local --cores N`
- **Profiles**: `brca`, `nbl`, `rbl`, `heme`, `pan_cancer` (defined in `pipeline/config/config.yaml`)
- **Data requirement**: Input RDS/CSV files under `pipeline/data/` are **not** included in the repo (large bioinformatics data). The pipeline will fail at runtime with `MissingInputException` if data files are absent. This is expected.
- **Conda envs for rules**: Snakemake auto-creates conda environments from `pipeline/envs/*.yaml` when run with `--use-conda`.

### Known Gotchas

- The `brca` profile is missing `cell_vst_rds` and `tumour_vst_rds` config keys, causing a `KeyError` during Snakefile parsing. Use `rbl` or `heme` profiles for dry-run validation instead.
- The `smoke_test` rule has unescaped shell braces (`{...}`) that cause a Snakemake `NameError`. This is a known issue in the current Snakefile.
- Expression Transformer training (`train_mgm.py`) is **GPU-only** — it exits with error code 3 if CUDA is not available. Use `export_embeddings.py --help` or model architecture validation on CPU to verify the module.
- Some Snakefile rules reference absolute paths to conda envs (e.g., `/work/ugbogu/pipeline/envs/tcga-r-env.yaml`). These are flagged by `snakemake --lint` but don't affect local env setup since Snakemake resolves them at runtime.

### Lint and Testing

- **Python linting**: `ruff check dsmz/transformer/` (no formal config exists; use `--select E,F,W --ignore E501`)
- **Snakemake linting**: `snakemake --lint --snakefile pipeline/Snakefile --configfile pipeline/config/config.yaml --config pipeline_profile=rbl`
- **Python compilation check**: `python3 -m py_compile <file.py>`
- **Pipeline dry-run**: `cd pipeline && ./run_pipeline.sh --profile rbl --executor local --cores 2 -n`
- **No automated test suite** — this is a research pipeline; validation is done via dry-runs and output contract checks (see `pipeline/AGENT.md` §9).
