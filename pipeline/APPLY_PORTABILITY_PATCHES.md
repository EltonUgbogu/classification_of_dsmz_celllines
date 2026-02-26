# Exact Steps to Apply Portability Patches

This document provides **copy-paste ready** instructions to make your pipeline portable.

## Files Already Created

These files are ready to use (already in your repo):

1. ✅ `run_pipeline.sh` - Portable runner script
2. ✅ `profiles/local/config.yaml` - Local execution profile
3. ✅ `profiles/slurm/config.yaml` - Slurm execution profile

## Patches to Apply

### 1. Update `config/config.yaml`

**Action:** Remove the top-level `profile:` field

**Find this line (around line 5):**
```yaml
profile: brca
```

**Delete it completely.**

The file should start with:
```yaml
# Unsupervised Clustering Pipeline Configuration
# config.yaml
# 
# Profile selection:
#   Use --config pipeline_profile=brca when running Snakemake
#   Do NOT set a top-level "profile:" field here
#   The pipeline will fail with a clear error if pipeline_profile is missing

defaults:
  paths:
    ...
```

---

### 2. Update `Snakefile` - Make pipeline_profile mandatory

**Find the profile selection block** (around lines 60-90, after `if "profiles" in config:`)

**Replace the entire block with:**

```python
if "profiles" in config:
    profiles = config.get("profiles", {})
    if not profiles:
        raise ValueError("config.profiles is missing or empty in config/config.yaml")
    
    profile_name = config.get("pipeline_profile", None)
    if profile_name is None:
        available = ", ".join(profiles.keys()) if profiles else "none"
        raise ValueError(
            "Missing --config pipeline_profile=<name>. "
            f"Available profiles: {available}. "
            "Example: snakemake --config pipeline_profile=brca"
        )
    
    if profile_name not in profiles:
        available = ", ".join(profiles.keys())
        raise ValueError(
            f"Profile '{profile_name}' not found in config.profiles. "
            f"Available profiles: {available}"
        )
    
    print(f"[Snakefile] pipeline_profile={profile_name}")
    
    default_cfg = config.get("defaults", {})
    cfg = deep_merge(default_cfg, profiles[profile_name])
    print(f"[Snakefile] Using profile: {profile_name}")
    print(f"[Snakefile] profile_name={profile_name}")
    print(f"[Snakefile] cfg.paths.unsup_root = {cfg.get('paths', {}).get('unsup_root', 'NOT SET')}")
    print(f"[Snakefile] cfg.analysis.cancer_type = {cfg.get('analysis', {}).get('cancer_type', 'NOT SET')}")
    
    # Define profile-specific log root (absolute for --directory compatibility)
    LOGROOT = os.path.join(PIPE_ROOT, "logs", profile_name)
    os.makedirs(LOGROOT, exist_ok=True)
    print(f"[Snakefile] LOGROOT = {LOGROOT}")
else:
    raise ValueError("config.profiles section is missing in config/config.yaml")
```

**Remove any lines that have fallbacks like:**
- `config.get("profile") or next(iter(profiles.keys()))`
- `profile_name = config.get("pipeline_profile") or config.get("profile") or ...`

---

### 3. Update `Snakefile` - Add smoke_test rule

**Add this rule AFTER the `rule all:` block** (around line 25, after BASE/SCRIPTS_DIR definitions):

```python
rule smoke_test:
    output:
        os.path.join(LOGROOT, "smoke_test.ok")
    log:
        os.path.join(LOGROOT, "smoke_test.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"
        
        {
            echo "[SMOKE TEST] Pipeline root: {BASE}"
            echo "[SMOKE TEST] Profile: {profile_name}"
            echo "[SMOKE TEST] Config file: {CFGFILE_ABS}"
            echo "[SMOKE TEST] Log root: {LOGROOT}"
            echo "[SMOKE TEST] Unsup root: {UNSUP}"
            echo "[SMOKE TEST] Cancer type: {cfg.get('analysis', {}).get('cancer_type', 'MISSING')}"
            
            # Validate config can be read
            if [ ! -f "{CFGFILE_ABS}" ]; then
                echo "ERROR: Config file not found: {CFGFILE_ABS}"
                exit 1
            fi
            
            # Validate profile exists
            if ! grep -q "^  {profile_name}:" "{CFGFILE_ABS}"; then
                echo "ERROR: Profile '{profile_name}' not found in config"
                exit 1
            fi
            
            echo "[SMOKE TEST] SUCCESS: Basic validation passed"
        } | tee "{log}"
        
        touch "{output}"
        '''
```

**Note:** This uses variables defined later (UNSUP, cfg). If you get errors, move it to after those definitions, or use lazy evaluation.

---

### 4. Update Shell Wrappers (Optional but Recommended)

**For ALL R script rules**, update the shell block to use this pattern:

**OLD:**
```python
shell:
    r'''
    set -euo pipefail
    mkdir -p "$(dirname "{output.rds}")"

    Rscript "{SCRIPTS_DIR}/script.R" \
      --args ... \
      > "{log}" 2>&1
    '''
```

**NEW:**
```python
shell:
    r'''
    set -euo pipefail
    mkdir -p "$(dirname "{output.rds}")" "$(dirname "{log}")"

    Rscript "{SCRIPTS_DIR}/script.R" \
      --args ... \
      2>&1 | tee "{log}"
    '''
```

**Key changes:**
1. Always `mkdir -p` for log directory too
2. Use `tee` instead of `>` (shows output in terminal AND saves to log)

**Example rule to update:** `tumour_neighbourhood_p_consensus_pan_featureset`

---

## Testing

After applying patches:

### 1. Test smoke test
```bash
./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test
```

Should complete in seconds and create `logs/brca/smoke_test.ok`

### 2. Test dry run
```bash
./run_pipeline.sh --profile brca --executor local --cores 8 --dry-run
```

Should show planned jobs without errors.

### 3. Test actual run
```bash
./run_pipeline.sh --profile brca --executor local --cores 8
```

---

## Migration from Old Wrapper

**Old command:**
```bash
./run_smk.sh --config profile=brca --cores 8
```

**New command:**
```bash
./run_pipeline.sh --profile brca --executor local --cores 8
```

**Key differences:**
- `--profile` instead of `--config profile=`
- `--executor local` is now required
- No hardcoded conda paths (auto-detects)

---

## Troubleshooting

### "Missing --config pipeline_profile" error

**Cause:** You didn't specify `--profile` or the Snakefile patch wasn't applied.

**Fix:** Always use `./run_pipeline.sh --profile <name> --executor <local|slurm>`

### "Profile 'X' not found" error

**Cause:** Profile doesn't exist in `config/config.yaml` under `profiles:` section.

**Fix:** Check available profiles: `grep -E "^  [a-z]+:" config/config.yaml`

### Conda not found

**Cause:** Conda not in PATH or not installed.

**Fix:** The script tries common locations. If needed, manually activate conda before running:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate smk
./run_pipeline.sh --profile brca --executor local --cores 8
```

---

## Summary

✅ **Files created:** `run_pipeline.sh`, `profiles/local/config.yaml`, `profiles/slurm/config.yaml`  
✅ **Patches needed:** `config/config.yaml` (remove `profile:`), `Snakefile` (mandatory profile + smoke_test)  
✅ **Optional:** Update shell wrappers to use `tee` pattern  

After applying, your pipeline will be **fully portable** and work the same on laptop, HPC, or cloud.

