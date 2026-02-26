# Unsupervised Clustering and Tumour Neighbourhood Analysis Pipeline

A portable, HPC-independent pipeline for pan-cancer clustering and tumour neighbourhood analysis.

## Quick Start

### Prerequisites

- Python 3.8+ with Snakemake installed
- Conda/Mamba (for environment management)
- R 4.0+ with required packages (managed via conda)

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd pipeline

# Install Snakemake (if not already installed)
conda create -n smk -c conda-forge -c bioconda snakemake
conda activate smk
```

### Running the Pipeline

**Basic usage:**
```bash
# Local execution (laptop/single machine)
./run_pipeline.sh --profile brca --executor local --cores 8

# Slurm execution (HPC cluster)
./run_pipeline.sh --profile brca --executor slurm --jobs 200

# Run specific targets
./run_pipeline.sh --profile brca --executor local --cores 4 smoke_test
```

**Available profiles:**
Valid profiles are the keys under `profiles:` in `config/config.yaml`:
- `brca` - Breast cancer analysis
- `nbl` - Neuroblastoma analysis  
- `rbl` - Retinoblastoma analysis

*(The runner will list available profiles if you omit `--profile`)*

**Available executors:**
- `local` - Run on local machine (uses `--cores`)
- `slurm` - Submit to Slurm cluster (uses `--jobs`)

### Pan-Cancer Analysis

Pan-cancer analysis combines data from multiple disease profiles. **You don't use a separate `--profile pan_cancer`** - instead, you specify a pan-cancer target rule:

```bash
# Run full pan-cancer pipeline (clustering + consensus + graphs)
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_like_per_disease_targets

# Run only pan-cancer graphs (if clustering already done)
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_graph_targets

# Run pan-cancer joint benchmark
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_joint_benchmark
```

**Important:** 
- `--profile brca` (or any profile) is still required - it sets log paths and config resolution
- Pan-cancer **profiles included** come from `config.pan_cancer.profiles` (default: `["brca", "nbl", "rbl"]`)
- Pan-cancer outputs go to `results/pan_cancer/joint_benchmark/` (not profile-specific)

**To override which profiles are included** (without editing config):
```bash
./run_pipeline.sh --profile brca --executor local --cores 8 \
  --config "pan_cancer_profiles=['brca','rbl']" \
  pan_cancer_like_per_disease_targets
```

### Smoke Test

Validate your setup with a quick smoke test:
```bash
./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test
```

This runs a minimal end-to-end test that should complete in minutes.

## Project Structure

```
pipeline/
├── config/
│   └── config.yaml          # Main configuration (profiles, paths, methods)
├── profiles/
│   ├── local/               # Local execution profile
│   └── slurm/               # Slurm execution profile
├── scripts/                 # R scripts for analysis steps
├── workflow/
│   └── lib/                 # Python utilities
├── run/                     # Snakemake working directories (per profile)
├── results/                 # Pipeline outputs
├── logs/                    # Execution logs
└── run_pipeline.sh          # Main entry point
```

## Configuration

### Profile Selection

Profiles are selected via `--config pipeline_profile=<name>`:
- **Required**: Always specify `--profile <name>` when running
- **No defaults**: Pipeline will fail with clear error if missing
- **Available**: `brca`, `nbl`, `rbl`

### Paths

All paths in `config/config.yaml` are **relative** to the pipeline root:
- Inputs: `data/...`
- Outputs: `results/...`
- Logs: `logs/...`

The pipeline automatically converts these to absolute paths internally.

### Pan-Cancer Analysis

Pan-cancer analysis combines data from multiple profiles:
- Configured in `config/config.yaml` under `pan_cancer:`
- Uses `cosine`/`euclidean` distance metrics (not `euc`/`corr`)
- Excludes PCA by default (`disable_pca_everywhere: true`)

## Execution Profiles

### Local Profile (`profiles/local/`)

For running on a single machine:
- Uses `--cores` for parallelization
- Default: 8 cores (auto-detected from system)
- Resources: 4GB RAM, 60min runtime per job

### Slurm Profile (`profiles/slurm/`)

For HPC clusters:
- Uses `--jobs` for concurrent job submission
- Default: 200 jobs
- Resources: Configurable per rule
- Partition: `cpu` (customize in `profiles/slurm/config.yaml`)

## Troubleshooting

### "Profile not found" error

Ensure you're using `--profile <name>` and the profile exists in `config/config.yaml`:
```bash
grep -E "^  [a-z]+:" config/config.yaml
```

### "Missing pipeline_profile" error

Always specify `--profile`:
```bash
./run_pipeline.sh --profile brca --executor local --cores 8
```

### Conda environment issues

The pipeline uses conda environments automatically. If you see conda errors:
1. Ensure conda is in your PATH
2. Activate the `smk` environment: `conda activate smk`
3. Check that `--use-conda` is enabled (default)

### Path issues

All paths are relative to the pipeline root. If you see "file not found" errors:
1. Ensure you're running from the pipeline root directory
2. Check that `config/config.yaml` paths are correct
3. Verify input files exist in `data/...`

## Development

### Adding a New Profile

1. Add profile section to `config/config.yaml`:
```yaml
profiles:
  new_profile:
    analysis:
      cancer_type: NEW
    paths:
      unsup_root: "results/unsupervised/new_profile"
      # ... other paths
```

2. Test with smoke test:
```bash
./run_pipeline.sh --profile new_profile --executor local --cores 2 smoke_test
```

### Adding a New Rule

Follow the standard pattern:
```python
rule my_new_rule:
    input: ...
    output: ...
    log: os.path.join(LOGROOT, "my_new_rule.log")
    conda: os.path.join(BASE, "envs", "tcga-r-env.yaml")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"
        
        Rscript "{SCRIPTS_DIR}/my_script.R" \
          --config "{input.cfg}" \
          --profile "{profile_name}" \
          --workdir "{BASE}" \
          --input "{input}" \
          --output "{output}" \
          2>&1 | tee "{log}"
        '''
```

## License

[Add your license here]

## Citation

[Add citation information here]
