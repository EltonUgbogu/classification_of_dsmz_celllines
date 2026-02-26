# Understanding Profile Selection vs Pan-Cancer Analysis

## Quick Answer

- **`--profile`**: Selects a **single disease profile** (brca, nbl, rbl) for per-disease analysis
- **Pan-cancer**: Uses **target rules** (not a separate profile) and combines data from multiple profiles

## Per-Disease Profiles

### What `--profile` does

Selects one disease profile from `config/config.yaml` → `profiles:`:

```yaml
profiles:
  brca: ...
  nbl: ...
  rbl: ...
```

### Usage

```bash
# BRCA analysis
./run_pipeline.sh --profile brca --executor local --cores 8

# NBL analysis  
./run_pipeline.sh --profile nbl --executor local --cores 8

# RBL analysis
./run_pipeline.sh --profile rbl --executor local --cores 8
```

### What it controls

- **Output directory**: `results/unsupervised/{profile}/`
- **Log directory**: `logs/{profile}/`
- **Config resolution**: Merges `defaults:` + `profiles.{profile}:`
- **Input paths**: Profile-specific data paths

### Available profiles

The runner will list them if you omit `--profile`:
```bash
$ ./run_pipeline.sh
[ERROR] --profile and --executor are required.

Available profiles (from config/config.yaml):
  - brca
  - nbl
  - rbl
```

---

## Pan-Cancer Analysis

### How to run pan-cancer

**You don't use `--profile pan_cancer`**. Instead, you specify a **pan-cancer target rule**:

```bash
# Full pan-cancer pipeline (clustering + consensus + graphs)
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_like_per_disease_targets

# Only pan-cancer graphs
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_graph_targets

# Pan-cancer joint benchmark
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_joint_benchmark
```

### Why `--profile` is still required

Even for pan-cancer, you must specify `--profile` (typically `brca`):

- Sets **log directory**: `logs/{profile}/pan_cancer/...`
- Used for **config path resolution**
- Can be any valid profile (doesn't affect which diseases are included)

### What controls pan-cancer disease inclusion

Pan-cancer **profiles included** come from `config.pan_cancer.profiles`:

```yaml
pan_cancer:
  profiles: ["brca", "nbl", "rbl"]  # These are included
  outdir: "results/pan_cancer/joint_benchmark"
  ...
```

**Default**: All profiles listed in `pan_cancer.profiles` are included.

### Pan-cancer outputs

- **Location**: `results/pan_cancer/joint_benchmark/` (shared, not profile-specific)
- **Logs**: `logs/{profile}/pan_cancer/...` (profile-specific for organization)

---

## Overriding Pan-Cancer Profiles

To run pan-cancer with only a subset of profiles (without editing config):

```bash
./run_pipeline.sh --profile brca --executor local --cores 8 \
  --config "pan_cancer_profiles=['brca','rbl']" \
  pan_cancer_like_per_disease_targets
```

**Note**: This requires the Snakefile to support `pan_cancer_profiles` override (currently uses `config.pan_cancer.profiles`).

---

## Examples Summary

### Per-disease analysis
```bash
# BRCA only
./run_pipeline.sh --profile brca --executor local --cores 8

# Outputs: results/unsupervised/brca/
# Logs:    logs/brca/
```

### Pan-cancer analysis
```bash
# Pan-cancer (includes brca, nbl, rbl by default)
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_like_per_disease_targets

# Outputs: results/pan_cancer/joint_benchmark/
# Logs:    logs/brca/pan_cancer/  (profile used for log organization)
```

### Both in one session
```bash
# Run BRCA analysis
./run_pipeline.sh --profile brca --executor local --cores 8

# Then run pan-cancer (uses outputs from all profiles)
./run_pipeline.sh --profile brca --executor local --cores 8 pan_cancer_like_per_disease_targets
```

---

## Key Takeaways

1. ✅ **`--profile`** = single disease profile (brca/nbl/rbl)
2. ✅ **Pan-cancer** = target rule (not a profile)
3. ✅ **`--profile` required for pan-cancer** = sets log paths (can be any profile)
4. ✅ **Pan-cancer diseases** = from `config.pan_cancer.profiles` (not from `--profile`)
5. ✅ **Pan-cancer outputs** = `results/pan_cancer/` (shared, not profile-specific)

---

## Common Confusion

**Q: Why do I need `--profile brca` for pan-cancer if it includes all diseases?**

**A:** The `--profile` argument controls:
- Where logs go (`logs/brca/pan_cancer/`)
- Config path resolution
- It does **NOT** control which diseases are included in pan-cancer

The diseases included come from `config.pan_cancer.profiles`, not from `--profile`.

**Q: Can I use `--profile nbl` for pan-cancer?**

**A:** Yes! It will work the same, just logs will go to `logs/nbl/pan_cancer/` instead. The pan-cancer outputs and included diseases are unchanged.

