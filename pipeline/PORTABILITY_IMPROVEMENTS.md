# Portability Improvements Summary

This document summarizes the changes made to improve pipeline portability and HPC-independence.

## Changes Made

### 1. Portable Pipeline Runner (`run_pipeline.sh`)

**Created:** `run_pipeline.sh` - A single entrypoint that works everywhere

**Features:**
- Works on laptop, HPC, or cloud
- Supports `--executor local|slurm`
- Auto-detects cores/jobs
- No hardcoded paths
- Clear error messages

**Usage:**
```bash
./run_pipeline.sh --profile brca --executor local --cores 8
./run_pipeline.sh --profile brca --executor slurm --jobs 200
```

**Replaces:** Old `run_smk.sh` with HPC-specific paths

---

### 2. Snakemake Execution Profiles

**Created:**
- `profiles/local/config.yaml` - For local execution
- `profiles/slurm/config.yaml` - For Slurm HPC clusters

**Features:**
- Resource defaults (mem, runtime)
- Rule-specific resource overrides
- Clean separation of execution backend from pipeline logic

**Benefits:**
- No Slurm-specific code in Snakefile
- Easy to add new executors (SGE, PBS, etc.)
- Consistent resource management

---

### 3. Configuration Improvements

**Updated:** `config/config.yaml`

**Changes:**
- Removed top-level `profile: brca` field
- Made `pipeline_profile` mandatory (no silent fallbacks)
- Clear error messages when profile is missing

**Updated:** `Snakefile`

**Changes:**
- Strict validation: `pipeline_profile` is REQUIRED
- No fallback to config.yaml `profile:` field
- Clear error messages with available profiles listed

---

### 4. Configuration Utilities Module

**Created:** `workflow/lib/config_utils.py`

**Functions:**
- `abspath()` - Convert relative to absolute paths
- `deep_merge()` - Recursively merge configs
- `get_profile_cfg()` - Load profile configs
- `load_profile()` - Convenience wrapper

**Benefits:**
- Reusable config logic
- Testable independently
- Cleaner Snakefile

---

### 5. Improved Shell Wrappers

**Updated:** Example rule (`tumour_neighbourhood_p_consensus_pan_featureset`)

**Changes:**
- Added `tee` for dual output (terminal + log)
- Always create log directory
- Consistent error handling with `set -euo pipefail`

**Pattern:**
```python
shell:
    r'''
    set -euo pipefail
    mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"
    
    Rscript "{SCRIPTS_DIR}/script.R" \
      --args ... \
      2>&1 | tee "{log}"
    '''
```

---

### 6. Smoke Test Rule

**Added:** `rule smoke_test` to Snakefile

**Purpose:**
- Quick validation (<5 minutes)
- Tests config loading, path resolution
- Validates profile exists
- Can be run in CI

**Usage:**
```bash
./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test
```

---

### 7. Documentation

**Created:** `README.md`

**Contents:**
- Quick start guide
- Installation instructions
- Usage examples
- Troubleshooting section
- Development guidelines

---

## Migration Guide

### For Existing Users

1. **Update your run command:**
   ```bash
   # Old
   ./run_smk.sh --config profile=brca --cores 8
   
   # New
   ./run_pipeline.sh --profile brca --executor local --cores 8
   ```

2. **Always specify `--profile`:**
   - No more silent fallbacks
   - Pipeline will fail with clear error if missing

3. **Use execution profiles:**
   - Local: `--executor local`
   - Slurm: `--executor slurm`
   - Customize resources in `profiles/*/config.yaml`

### For New Users

1. Clone repository
2. Install Snakemake: `conda create -n smk -c conda-forge snakemake`
3. Run smoke test: `./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test`
4. Run full pipeline: `./run_pipeline.sh --profile brca --executor local --cores 8`

---

## Path Portability

### All Paths Are Relative

- **Inputs:** `data/...` (relative to pipeline root)
- **Outputs:** `results/...` (relative to pipeline root)
- **Logs:** `logs/...` (relative to pipeline root)

### No Hardcoded Paths

- No `/work/ugbogu/...` anywhere
- All paths derived from `workflow.basedir`
- Works from any directory

---

## Environment Portability

### Conda Strategy

- Environments defined in `envs/*.yaml`
- Snakemake manages environments automatically
- No manual conda activation required

### Future: Container Support

To add Singularity/Apptainer support:
1. Create `containers/` directory
2. Build images: `singularity build containers/r-base.sif r-base.def`
3. Add to rules: `singularity: "containers/r-base.sif"`

---

## Testing

### Smoke Test
```bash
./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test
```

### Dry Run
```bash
./run_pipeline.sh --profile brca --executor local --cores 8 --dry-run
```

### Lint
```bash
snakemake --snakefile Snakefile --lint
```

---

## Next Steps

1. **Add CI:** GitHub Actions to run smoke test on every push
2. **Add containers:** Singularity/Apptainer images for complete portability
3. **Add more executors:** SGE, PBS profiles if needed
4. **Documentation:** Add more examples and troubleshooting

---

## Benefits Summary

✅ **Portable:** Works on laptop, HPC, cloud  
✅ **Clear errors:** No silent failures  
✅ **Maintainable:** Clean separation of concerns  
✅ **Testable:** Smoke test validates setup  
✅ **Documented:** README and inline docs  
✅ **Professional:** Follows Snakemake best practices  

---

## Questions?

See `README.md` for usage examples and troubleshooting.

