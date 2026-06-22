# Conda Setup Guide for Snakemake Pipeline

## Problem

Your current conda installation (from Spack) is broken:
- `CONDA_EXE` points to `/srv/software/spack/.../bin/conda`
- Running `python3 -c "import conda"` fails with `ModuleNotFoundError`
- Snakemake cannot create conda environments

## Solution: Install Your Own Miniforge

### Option 1: Automated Setup (Recommended)

Run the setup script:

```bash
cd classification_of_dsmz_celllines
./setup_conda.sh
```

This will:
1. Download Miniforge3 to `$HOME/miniforge3`
2. Install it
3. Initialize conda for your shell

After running, **restart your shell** or run:
```bash
source ~/.zshrc  # or ~/.bashrc if using bash
```

### Option 2: Manual Installation

```bash
# Download Miniforge
cd $HOME
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Install
bash Miniforge3.sh -b -p $HOME/miniforge3

# Initialize
$HOME/miniforge3/bin/conda init zsh  # or 'bash' if using bash

# Restart shell
exec zsh
```

## Verify Installation

```bash
# Check conda is in PATH
which conda
# Should show: $HOME/miniforge3/bin/conda

# Check version
conda --version

# Verify Python module works
python3 -c "import conda; print('OK')"
# Should print: OK
```

## Configure Snakemake to Use Your Conda

### Method 1: Use the `smk` Environment (Recommended)

The `unsupervised_pipeline.sh` script already handles this automatically:
- It creates/activates the environment from `envs/smk.yaml`
- This environment contains Snakemake

### Method 2: Manual Snakemake Setup

```bash
# Create Snakemake environment
conda env create -f envs/smk.yaml

# Activate it
conda activate smk

# Verify Snakemake works
snakemake --version
```

## Environment Files

Your pipeline has three environment files:

1. **`envs/smk.yaml`**: Snakemake + basic R packages
   - Used by `unsupervised_pipeline.sh` automatically
   - Contains: `snakemake-minimal`, `r-base`, `r-dplyr`, `r-matrixstats`, `r-ggplot2`

2. **`envs/r-base.yaml`**: Minimal R environment (NEW)
   - For scripts that only need `optparse`, `yaml`, `dplyr`
   - Use this for `rebuild_joint_vst.R` and similar simple scripts
   - Much faster to create than `tcga-r-env.yaml`

3. **`envs/tcga-r-env.yaml`**: Full R/Bioconductor environment
   - Used by most R rules in the Snakefile
   - Contains all Bioconductor packages, DESeq2, etc.
   - Large but comprehensive

## Using Different Environments Per Script

### Current Setup

All R rules in `Snakefile` currently use `tcga-r-env.yaml`:
```python
conda: os.path.join(BASE, "envs", "tcga-r-env.yaml")
```

### Recommended: Use Minimal Env for Simple Scripts

For scripts like `rebuild_joint_vst.R` that only need `optparse` and `yaml`, you can use the lighter `r-base.yaml`:

```python
rule rebuild_joint_vst:
    conda: os.path.join(BASE, "envs", "r-base.yaml")  # Instead of tcga-r-env.yaml
    shell:
        r"""
        Rscript --vanilla {SCRIPTS_DIR}/rebuild_joint_vst.R ...
        """
```

This will:
- Create environments faster
- Use less disk space
- Still work correctly

## Troubleshooting

### "conda: command not found"

Make sure conda is in your PATH:
```bash
export PATH="$HOME/miniforge3/bin:$PATH"
```

Or restart your shell after running `conda init`.

### "ModuleNotFoundError: No module named 'conda'"

You're still using the broken Spack conda. Fix it:

```bash
# Unset Spack conda
unset CONDA_EXE
hash -r

# Use your Miniforge conda
export PATH="$HOME/miniforge3/bin:$PATH"

# Verify
python3 -c "import conda; print('OK')"
```

### Snakemake Can't Create Environments

1. Verify conda works: `python3 -c "import conda"`
2. Make sure `--use-conda` flag is set
3. Check Snakemake version: `snakemake --version` (should be 5.10+)

### Manual Script Execution

When running R scripts manually (outside Snakemake), use the conda environment:

```bash
# Option 1: Use Snakemake-created env
"$CONDA_PREFIX/bin/Rscript" rebuild_joint_vst.R ...

# Option 2: Use your conda env
conda activate tcga-r-env  # or r-base
Rscript rebuild_joint_vst.R ...
```

## Summary

1. **Install Miniforge**: Run `./setup_conda.sh` or follow manual steps
2. **Restart shell**: `source ~/.zshrc` or `exec zsh`
3. **Verify**: `python3 -c "import conda"`
4. **Use Snakemake**: The `unsupervised_pipeline.sh` script handles everything automatically

Your pipeline will now work with proper conda environments!
