# AGENT.md — Unsupervised Clustering + Tumour Neighbourhood Pipeline

This document is the **execution + output contract** for running the pipeline end-to-end and understanding the **expected outputs**.  
Goal: run one profile (BRCA/NBL/RBL), produce the full result tree, then proceed directly to downstream analysis.

---

## 1) What this pipeline does (high-level)

For a given profile (e.g. `brca`, `nbl`, `rbl`), the pipeline:

1. **Feature selection** on a joint VST matrix  
   - Produces HVG list and feature-set gene lists for multiple methods.

2. **Agnostic clustering** across multiple "directions"  
   - Directions = feature method × distance (+ PAM50 directions if enabled)
   - Algorithms = HC and PCA+HC for all; KMeans variants for Euclidean directions only.

3. **Consensus clustering** per direction and clustering-kind  
   - Produces "optimal cluster" RDS outputs for multiple consensus kinds.

4. **Tumour neighbourhood analysis**
   - Builds tumour-neighbourhood input matrices for:
     - HVG (always)
     - PAM50 (if enabled)
     - Feature-set methods (Variance, MAD, …) as separate submatrices
   - Computes tumour neighbourhoods (HC; plus KMeans for euclidean directions)
   - Computes p-consensus outputs per direction

5. **DSMZ–DSMZ similarity graph** per direction  
   - Builds similarity matrix, graph edges, community summaries, and figures

6. **Cross-direction summary** across all directions  
   - Aggregates p-consensus results and ranks best cell lines

7. **QC + integrated UMAP** per direction  
   - Produces QC summary tables and UMAP plots/tables

Final signal that a profile is complete:
- `results/unsupervised/<profile>/pipeline_targets.done`

---

## 2) Repository expectations (layout)

Expected repo structure:

```
pipeline/
├── Snakefile
├── AGENT.md
├── README.md
├── run_smk.sh
├── config/
│   └── config.yaml
├── envs/
│   ├── smk.yaml
│   ├── tcga-r-env.yaml
│   └── tumour_nh_qc.yaml
├── scripts/
│   ├── lib_config.R
│   ├── feature_selection_unsupervised.R
│   ├── build_agnostic_direction_mats.R
│   ├── pca_hc_cell.R
│   ├── ...
│   ├── build_tumour_neighbourhood_input.R
│   ├── build_dsmz_tcga_pam50matrix.R
│   ├── comp_tumour_neighbourhoods.R
│   ├── tumour_neighbourhood_p_consensus.R
│   ├── compute_dsmz_dsmz_similarity.R
│   ├── summarize_p_consensus_all.R
│   └── tumour_neighbourhood_qc_umap.R
├── R/
│   ├── base_functions/
│   └── ...
├── data/
│   ├── brca/
│   ├── nbl/
│   ├── rbl/
│   ├── dsmz/
│   └── resources/
└── results/
    └── unsupervised/
        ├── brca/
        ├── nbl/
        └── rbl/
```

Important conventions:
- **Inputs** live under `data/`
- **All outputs** live under `results/unsupervised/<profile>/...`
- **Logs** for Snakemake jobs are written under `./logs/`

---

## 3) Profiles and config model

The pipeline uses a **profile** to select a dataset and toggle features like PAM50.

**Important**: The actual source of truth is `config/config.yaml`. The profile determines which section (`profiles.<profile>`) is merged with `defaults` to produce the final config.

Typical configuration:
- `profile: brca` → typically has `use_pam50: true` in config (PAM50 inputs/outputs expected)
- `profile: nbl` → typically has `use_pam50: false` in config (PAM50 steps skipped)
- `profile: rbl` → typically has `use_pam50: false` in config (PAM50 steps skipped)

**Note**: PAM50-related outputs (e.g. `expr_pam50.rds`, `pam50_euc`/`pam50_corr` directions) are **only required when `use_pam50: true`** in the merged config. For profiles with `use_pam50: false`, completion criteria do not include PAM50 outputs.

Config file:
- `config/config.yaml`

Invocation overrides:
- `--config profile=<profile>`

To verify the active profile and merged config:

```bash
# Show available profiles in YAML
python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config/config.yaml"))
print("Available profiles:", ", ".join((cfg.get("profiles") or {}).keys()))
PY

# Verify the profile is recognized by Snakemake (source of truth)
./run_smk.sh -n --config profile=brca 2>&1 | grep -E "Using profile|profile="
```

---

## 4) One-time setup

### 4.1 Create logs directory
From repo root:

```bash
mkdir -p logs
```

### 4.2 Conda reliability on HPC

**Critical**: Never call `snakemake` directly on HPC. Always use `./run_smk.sh`.

If you run `snakemake` directly, you may see:
- `Conda environments: ignored` (conda envs won't be created/activated)
- `ModuleNotFoundError: No module named 'conda'` (broken conda on PATH)

#### run_smk.sh (wrapper)

`run_smk.sh` ensures we use the **working** Miniforge and then runs Snakemake with `--use-conda`.

Expected behavior:

* It puts `/work/ugbogu/software/miniforge3/bin` first on PATH
* Activates `smk` environment (where snakemake is installed)
* Executes `snakemake --use-conda --conda-frontend conda ...`

If `run_smk.sh` is not executable:

```bash
chmod +x run_smk.sh
```

### 4.3 Strict channel priority (recommended)

```bash
/work/ugbogu/software/miniforge3/bin/conda config --set channel_priority strict
```

---

## 4.5) Preflight checks (run before execution)

These checks prevent ~80% of common failures. **Run these before your first execution**:

### Verify wrapper and conda setup

```bash
# Show which conda + snakemake you will use
./run_smk.sh --version || true
./run_smk.sh -h | head
which conda
conda info | head -n 5

# Verify conda works (should not error)
conda info --json >/dev/null && echo "✓ conda OK" || echo "✗ conda broken"
```

### Verify profile and DAG structure

```bash
# Confirm the profile expands to the expected output root
./run_smk.sh -n --config profile=brca 2>&1 | head -n 30

# Check that profile is recognized
./run_smk.sh -n --config profile=brca 2>&1 | grep -i "Using profile"
```

### Verify required input files exist (adjust per profile)

**Important**: Each profile requires its own VST matrix path under `data/<profile>/...` as configured in `config/config.yaml`.

For BRCA:

```bash
# Joint VST matrix (required for all profiles)
ls -lh data/brca/bcc_analysis_results/vst_dsmz_bcc_tcga_brca_batch_corrected.rds

# PAM50 inputs (BRCA only, if use_pam50: true)
ls -lh data/brca/bcc_preprocessing_results/TCGA-BRCA_PAM50_expr_matrix.rds
ls -lh data/brca/bcc_preprocessing_results/DSMZ-BCC_PAM50_expr_matrix.rds

# DSMZ metadata (required for all profiles)
ls -lh data/dsmz/DSMZ_metadata.csv
```

For NBL/RBL, check corresponding `data/nbl/...` or `data/rbl/...` paths. The required joint VST matrix path differs per profile and must exist under `data/<profile>/...` as configured in the profile-specific config section.

### Verify config structure

```bash
# Show available profiles in YAML
python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config/config.yaml"))
print("Available profiles:", ", ".join((cfg.get("profiles") or {}).keys()))
PY

# Verify the profile is recognized by Snakemake (source of truth)
./run_smk.sh -n --config profile=brca 2>&1 | grep -E "Using profile|profile=" || echo "Profile check: no explicit 'Using profile' message (this is OK if DAG builds correctly)"
```

**If any preflight check fails, fix it before running the pipeline.**

---

## 5) Execution commands (end-to-end)

### 5.1 Dry-run

Always sanity-check DAG:

```bash
./run_smk.sh -n --config profile=brca
```

Optional: show summary:

```bash
./run_smk.sh -n --summary --config profile=brca | head -n 60
```

### 5.2 Run the full pipeline for a profile

The canonical "complete everything" target is:

* `results/unsupervised/<profile>/pipeline_targets.done`

**Recommended execution flags** (for robustness on HPC):

* `-j <cores>`: number of parallel jobs (adjust to your allocation)
* `--rerun-incomplete`: re-run jobs if outputs are incomplete
* `--keep-going`: continue even if some jobs fail (useful for large batches)
* `--latency-wait 60`: wait up to 60s for filesystem sync (important on shared filesystems)

Run BRCA:

```bash
./run_smk.sh --config profile=brca -j 8 \
  --rerun-incomplete --keep-going --latency-wait 60 \
  results/unsupervised/brca/pipeline_targets.done
```

Run NBL:

```bash
./run_smk.sh --config profile=nbl -j 8 \
  --rerun-incomplete --keep-going --latency-wait 60 \
  results/unsupervised/nbl/pipeline_targets.done
```

Run RBL:

```bash
./run_smk.sh --config profile=rbl -j 8 \
  --rerun-incomplete --keep-going --latency-wait 60 \
  results/unsupervised/rbl/pipeline_targets.done
```

**Resource recommendations**:

* For local execution: `-j 4` to `-j 8` (depends on CPU cores)
* For SLURM cluster: `-j 20` to `-j 32` (adjust based on partition limits)
* If conda env builds are throttled: use `--conda-prefix /path/to/shared/envs` to reuse environments across runs

### 5.3 Re-run only a subset (examples)

Feature selection only:

```bash
./run_smk.sh --config profile=brca -j 1 \
  results/unsupervised/brca/feature_selection_unsupervised/feature_sets/genes_top500_MX.txt
```

Tumour neighbourhood input PAM50 (BRCA only):

```bash
./run_smk.sh --config profile=brca -j 1 \
  results/unsupervised/brca/tumour_neighbourhoods_input/expr_pam50.rds
```

---

## 6) Expected outputs (by stage)

All outputs are under:

* `results/unsupervised/<profile>/...`

Below is the **expected result tree** and "what to look for".

### 6.1 Feature selection outputs

Directory:

* `results/unsupervised/<profile>/feature_selection_unsupervised/`

✅ **Declared outputs (must exist)**:

* Joint ranks table (exact filename depends on profile/dataset naming):

  * `*_joint_ranks_*.tsv`
* MX final list (top-500 network-ranked genes):

  * `feature_sets/genes_top500_MX.txt`
* HVG trend-corrected list (top-3000 HVGs):

  * `feature_sets/genes_top3000_HVG.txt`
* Feature-set gene lists (one per method):

  * `feature_sets/genes_top{N}_{METHOD}.txt` where `{N}` matches the method-specific 
    top-N configured in `feature_selection.method_topn` (e.g. 3000 for HVG, 500 for MX).
    * Methods: `Variance`, `MAD`, `MeanAbsDev`, `Entropy`, `PCA`, `Spearman`, `MX`, `kTotal`, `HVG`

Example (BRCA):

```
results/unsupervised/brca/feature_selection_unsupervised/
├── BRCA_TCGA-DSMZ_joint_ranks_kTotal_vs_MX.tsv
└── feature_sets/
    ├── genes_top3000_HVG.txt
    ├── genes_top500_MX.txt
    ├── genes_top500_Variance.txt
    ├── genes_top500_MAD.txt
    ├── genes_top500_MeanAbsDev.txt
    ├── genes_top500_Entropy.txt
    ├── genes_top500_PCA.txt
    ├── genes_top500_Spearman.txt
    └── genes_top500_kTotal.txt
```

**Completion criteria**:

* All declared gene lists exist and are non-empty (check with `test -s <file>`).

---

### 6.2 Agnostic clustering outputs

Directory:

* `results/unsupervised/<profile>/agnostic_clustering/<direction>/...`

✅ **Declared outputs (must exist)**:

Per direction:

* `inputs/cell_expr.rds`
* `inputs/tumour_expr.rds`
* Cluster RDS outputs for:

  * `pca_hc_cell/pca_hc_cell_clusters_optimal.rds`
  * `pca_hc_tumour/pca_hc_tumour_clusters_optimal.rds`
  * `pca_hc_cell_tumour/pca_hc_cell_tumour_clusters_optimal.rds`
  * `hc_cell/hc_cell_clusters_optimal.rds`
  * `hc_tumour/hc_tumour_clusters_optimal.rds`
  * `hc_cell_tumour/hc_cell_tumour_clusters_optimal.rds`
* Euclidean-only directions additionally include:

  * `pca_kmeans_cell/pca_kmeans_cell_clusters_optimal.rds`
  * `pca_kmeans_tumour/pca_kmeans_tumour_clusters_optimal.rds`
  * `pca_kmeans_cell_tumour/pca_kmeans_cell_tumour_clusters_optimal.rds`
  * `kmeans_cell/kmeans_cell_clusters_optimal.rds`
  * `kmeans_tumour/kmeans_tumour_clusters_optimal.rds`
  * `kmeans_cell_tumour/kmeans_cell_tumour_clusters_optimal.rds`

Example fragment:

```
results/unsupervised/brca/agnostic_clustering/Variance_euc/
├── inputs/
│   ├── cell_expr.rds
│   └── tumour_expr.rds
├── pca_hc_cell/pca_hc_cell_clusters_optimal.rds
├── pca_hc_tumour/pca_hc_tumour_clusters_optimal.rds
├── pca_hc_cell_tumour/pca_hc_cell_tumour_clusters_optimal.rds
├── hc_cell/hc_cell_clusters_optimal.rds
├── hc_tumour/hc_tumour_clusters_optimal.rds
└── hc_cell_tumour/hc_cell_tumour_clusters_optimal.rds
```

**Completion criteria**:

* Every direction has all required cluster RDS files (as declared in Snakefile outputs).

---

### 6.3 Consensus clustering outputs

Directory:

* `results/unsupervised/<profile>/consensus/<direction>/<kind>/...`

✅ **Declared outputs (must exist)**:

* `<kind>_clusters_optimal.rds` for all consensus kinds.

Consensus kinds include:
* `ccp_hc_expr_cell_tumour`, `ccp_hc_pca_cell_tumour`
* `ccp_hc_expr_cell`, `ccp_hc_pca_cell`
* `ccp_hc_expr_tumour`, `ccp_hc_pca_tumour`
* `ccp_kmeans_expr_cell_tumour`, `ccp_kmeans_pca_cell_tumour`
* `ccp_kmeans_expr_cell`, `ccp_kmeans_pca_cell`
* `ccp_kmeans_expr_tumour`, `ccp_kmeans_pca_tumour`

Example:

```
results/unsupervised/brca/consensus/MAD_euc/ccp_hc_expr_cell_tumour/
└── ccp_hc_expr_cell_tumour_clusters_optimal.rds
```

**Completion criteria**:

* All consensus kinds exist for all directions configured (as declared in Snakefile outputs).

---

### 6.4 Tumour neighbourhood inputs (direction-aware)

Directory:

* `results/unsupervised/<profile>/tumour_neighbourhoods_input/`

✅ **Declared outputs (must exist)**:

* PAM50 (if `use_pam50: true` in config):

  * `expr_pam50.rds`
  * `cell_line_to_original_sample_id_pam50.rds`
* For each feature-set method (including HVG and MX):

  * `tumour_neighbourhoods_input/cell_line_to_original_sample_id_<METHOD>.rds`
  * `feature_selection_unsupervised/featuresets/<METHOD>/expr_submatrix.rds`

🟡 **Derived artifacts (may exist depending on script logic)**:

* Additional per-method diagnostics may appear under 
  `feature_selection_unsupervised/feature_sets/` (e.g. unique gene tables).

**Completion criteria**:

* All declared `.rds` files exist and are non-empty (check with `test -s <file>`).
* **Note**: PAM50 files are only required when `use_pam50: true` in the merged config.

---

### 6.5 Tumour neighbourhood results + p-consensus

Directory:

* `results/unsupervised/<profile>/tumour_neighbourhoods/<direction>/...`

✅ **Declared outputs (must exist)**:

* "done marker":

  * `.tumour_neighbourhoods_done` (for HC-based tumour neighbourhoods)
  * `.tumour_neighbourhoods_km_done` (for KMeans-based tumour neighbourhoods, euclidean directions only)
* p-consensus:

  * `final_consensus/Final_consensus_tumour_neighbourhoods_<direction>.rds`
  * `final_consensus/Final_consensus_tumour_neighbourhoods_<direction>.tsv`

**Completion criteria**:

* For every direction, the `final_consensus/*Final_consensus*` files exist.
* **Note**: PAM50 directions (`pam50_euc`, `pam50_corr`) are only required when `use_pam50: true` in config.

---

### 6.6 DSMZ–DSMZ similarity graph outputs

Directory:

* `results/unsupervised/<profile>/tumour_neighbourhoods/<direction>/final_consensus/`

✅ **Declared outputs (must exist)**:

* Similarity matrix + derived tables:

  * `DSMZ_DSMZ_similarity_matrix_<direction>.rds`
  * `DSMZ_DSMZ_similarity_pairs_<direction>.tsv`
  * `DSMZ_DSMZ_graph_edges_<direction>.tsv`
  * `DSMZ_DSMZ_graph_node_summary_<direction>.tsv`
  * `DSMZ_DSMZ_graph_node_annotations_<direction>.tsv`
  * `DSMZ_DSMZ_graph_community_summary_<direction>.tsv`
  * `DSMZ_DSMZ_Louvain_vs_Leiden_community_table_<direction>.tsv`
* Figures:

  * `Fig_DSMZ_DSMZ_similarity_histogram_<direction>.pdf`
  * `Fig_DSMZ_p_consensus_cell_scatter_<direction>.pdf`
  * `Fig_DSMZ_DSMZ_Louvain_vs_Leiden_heatmap_<direction>.pdf`
  * `Fig_DSMZ_DSMZ_graph_Leiden_<direction>.pdf`
  * `Fig_DSMZ_DSMZ_graph_Louvain_<direction>.pdf`
  * `Fig_DSMZ_DSMZ_graph_minimal_<direction>.pdf`

**Completion criteria**:

* For every direction, all declared similarity tables + PDFs exist.
* **Note**: PAM50 directions are only required when `use_pam50: true` in config.

---

### 6.7 Cross-direction p-consensus summary outputs

Directory:

* `results/unsupervised/<profile>/tumour_neighbourhoods/final_consensus_all/`

✅ **Declared outputs (must exist)**:

* `p_consensus_direction_summary.tsv`
* `p_consensus_cellline_direction_summary.tsv`
* `p_consensus_winner_map_per_cell_line.tsv`
* `p_consensus_winners_by_max_p.tsv`
* `p_consensus_winners_by_frac_ge_thr.tsv`
* `p_consensus_best_cell_lines_ranked.tsv`
* `p_consensus_cellline_direction_summary.long.tsv`
* `p_consensus_composite_weights.tsv`
* `p_consensus_best_cell_lines_top_fraction.tsv`
* `p_consensus_best_cell_lines_top_score.tsv`
* Figures:

  * `Fig_p_consensus_direction_comparison.pdf`
  * `Fig_p_consensus_composite_PCA_scree.pdf`
  * `Fig_p_consensus_composite_PCA_weights.pdf`

**Completion criteria**:

* All declared outputs exist, especially the ranked list (`p_consensus_best_cell_lines_ranked.tsv`).

---

### 6.8 QC + UMAP outputs

Directory:

* `results/unsupervised/<profile>/tumour_neighbourhoods/<direction>/qc/`

✅ **Declared outputs (must exist)**:

* `nh_qc_summary.tsv`
* `nh_umap.tsv`
* `nh_umap.pdf`

**Completion criteria**:

* For all directions: QC TSV + UMAP TSV + UMAP PDF exist.
* **Note**: PAM50 directions are only required when `use_pam50: true` in config.

---

### 6.9 Final pipeline target

File:

* `results/unsupervised/<profile>/pipeline_targets.done`

**Completion criteria**:

* This file exists (touch output), which implies all aggregated outputs exist.

---

## 7) Logs and debugging

All Snakemake job logs are written to:

* `./logs/*.log`

Snakemake's internal run log:

* `.snakemake/log/<timestamp>.snakemake.log`

Fast debug steps:

1. Identify failing rule + job:

   * read Snakemake stderr output
2. Open its log:

   * `cat logs/<rule_or_step>.log`
3. Re-run only that output target with verbose output:

   * `./run_smk.sh --config profile=brca -j 1 --printshellcmds --reason <target_file>`
   * This shows the exact shell command and why the rule is running.

**Quick contract validation** (check if all declared outputs are "ok"):

```bash
./run_smk.sh --config profile=brca -n results/unsupervised/brca/pipeline_targets.done --summary 2>&1 | tail -n +2 | awk -F'\t' '$7!="ok"{print}'
```

This shows any outputs that are not "ok" (missing, outdated, etc.).

---

## 8) Known failure classes and fixes

### 8.1 Conda is broken / wrong conda on PATH

Symptom:

* `ModuleNotFoundError: No module named 'conda'`
  Fix:
* Always run via `./run_smk.sh ...`

### 8.2 R packages missing (optparse/yaml/etc.)

Symptom (from logs):

- `Error in library(optparse) : there is no package called 'optparse'`
- `Error in loadNamespace(x) : there is no package called 'yaml'`

**Important**: before editing any env YAML, verify whether the package is already declared:

```bash
grep -nE "r-(yaml|optparse)\b" envs/tcga-r-env.yaml envs/tumour_nh_qc.yaml || true
```

If `r-yaml` / `r-optparse` are already present in the env YAMLs and you *still* see missing-package errors, treat it as one of:

1. **Wrong Rscript is being used** (system `Rscript`, not the conda `Rscript`), or
2. **Stale/corrupted Snakemake conda env cache** (`.snakemake/conda/*`)

#### Level-1 safe remediation

**A) Confirm the `Rscript` being used inside a Snakemake job**

Run with printshellcmds and inspect the command:

```bash
./run_smk.sh --config profile=brca -j 1 --printshellcmds --reason <target_file>
```

Also inspect the job log and look for an absolute `Rscript` path:

```bash
grep -R "Rscript" -n logs/*.log | tail -n 20 || true
```

**B) Purge Snakemake conda env cache and rebuild**

This is safe: it deletes only cached envs, not results.

```bash
rm -rf .snakemake/conda
./run_smk.sh --config profile=brca -j 1 <target_file>
```

If the package is *not* present in the env YAMLs, then add the missing dependency (`r-yaml`, `r-optparse`, etc.) to the relevant env file, and rerun (Snakemake will rebuild a new hashed env).

#### Common gotcha: concatenated env YAMLs

Each environment must be a separate YAML file under `envs/` (one `name:` per file).  
If multiple env definitions are accidentally concatenated into one YAML, conda behavior becomes unpredictable.

Check:

```bash
for f in envs/*.yaml; do echo "== $f"; grep -n '^name:' "$f"; done
```

### 8.3 Profile config not applied inside R scripts

Symptom:

* `Error in readRDS(cfg$paths$tcga_brca_pam50_expr %||% cfg$paths$tcga_pam50_expr) : bad 'file' argument`
* R scripts read `cfg$paths$...` but config actually stores keys under `profiles:<name>:paths:...`
* The expression evaluates to `NULL` instead of a file path string

**Root cause**: R scripts using `yaml::read_yaml()` directly don't merge `defaults + profiles[profile]`, so profile-specific keys are missing.

**Fix**:

* Standardize config reading in R:

  * Use `scripts/lib_config.R` which provides `read_profiled_config()` that merges `defaults + profiles[profile]` (same logic as Snakefile),
  * OR pass resolved file paths from Snakemake to R scripts as explicit CLI args.

**Debug checklist**:

1. Confirm config keys exist under the selected profile:

   ```bash
   # Inspect config structure
   cat config/config.yaml | grep -A 20 "profiles:" | grep -A 10 "brca:"
   ```

2. Confirm R config loader merges `defaults + profiles[profile]`:

   * Check that the R script sources `scripts/lib_config.R` and calls `read_profiled_config(opt$config)` instead of `yaml::read_yaml(opt$config)`

3. Confirm those fields are *strings* and point to existing `.rds`:

   ```bash
   # Debug-only: dump merged config paths (simplified merge check)
   python3 - <<'PY'
   import yaml
   cfg = yaml.safe_load(open("config/config.yaml"))
   defaults = cfg.get("defaults", {}).get("paths", {})
   profile_cfg = cfg.get("profiles", {}).get("brca", {}).get("paths", {})
   merged = {**defaults, **profile_cfg}
   print("tcga_brca_pam50_expr:", merged.get("tcga_brca_pam50_expr", "NOT FOUND"))
   PY
   ```

4. Verify the file exists:

   ```bash
   # Check if the resolved path exists (adjust profile name as needed)
   python3 -c "import yaml; cfg=yaml.safe_load(open('config/config.yaml')); print(cfg['profiles']['brca']['paths']['tcga_brca_pam50_expr'])" | xargs ls -lh
```

### 8.4 Automated orchestrator (optional)

An automated orchestrator script (`agent/agent_orchestrator.py`) is available for Level-1 safe remediation:

- **Error classification** aligned to AGENT.md §8
- **Auto-fixes**: conda cache purge (when R packages are missing despite being in env YAMLs), latency-wait increases, targeted reruns
- **Controlled stops** for config-merge bugs, OOM/time limits, or broken conda
- **Final report**: `logs/agent_report_<date>.md` with per-profile artifact checklist

**Usage**:

```bash
python3 agent/agent_orchestrator.py --profiles brca nbl rbl --cores 8 --max-retries 3
```

**Safety**: The orchestrator only performs Level-1 safe operations (cache purge, env YAML append, reruns). It **must stop** on config-merge bugs, resource limits, or broken conda (these require manual fixes).

---

## 9) "Ready for analysis" checklist

A profile is ready for downstream analysis if:

* `results/unsupervised/<profile>/pipeline_targets.done` exists
* `.../tumour_neighbourhoods/final_consensus_all/p_consensus_best_cell_lines_ranked.tsv` exists
* Per direction:

  * `.../tumour_neighbourhoods/<direction>/final_consensus/Final_consensus_tumour_neighbourhoods_<direction>.tsv`
  * `.../tumour_neighbourhoods/<direction>/qc/nh_umap.pdf`

Recommended "analysis entry points":

1. `final_consensus_all/p_consensus_best_cell_lines_ranked.tsv`
2. `final_consensus_all/p_consensus_winner_map_per_cell_line.tsv`
3. per-direction similarity graph outputs (`DSMZ_DSMZ_graph_*`)
4. per-direction QC (`nh_qc_summary.tsv` + `nh_umap.tsv` + `nh_umap.pdf`)

---

## 10) Canonical one-liners (copy/paste)

Run BRCA end-to-end:

```bash
./run_smk.sh --config profile=brca -j 8 results/unsupervised/brca/pipeline_targets.done
```

Run NBL end-to-end:

```bash
./run_smk.sh --config profile=nbl -j 8 results/unsupervised/nbl/pipeline_targets.done
```

Run RBL end-to-end:

```bash
./run_smk.sh --config profile=rbl -j 8 results/unsupervised/rbl/pipeline_targets.done
```

---

## 11) Output contract (what scripts must guarantee)

For every rule output declared in the Snakefile:

* The file must be created deterministically
* Output paths must be **profile-scoped**
* Output files must be non-empty unless explicitly marked as `touch()`

Where appropriate, scripts should:

* `stop()` with a clear error if input files are missing
* write summary TSVs with stable column names
* write figures as PDF into the declared paths
* avoid writing outside `results/unsupervised/<profile>/...`

---

## 12) Next step: analysis layer

Once `pipeline_targets.done` exists, proceed to analysis notebooks/scripts that:

* load `final_consensus_all/*ranked*.tsv`
* load direction-level `Final_consensus_tumour_neighbourhoods_*.tsv`
* load DSMZ–DSMZ graph outputs to explore community structure
* use QC UMAP embeddings for visualization and sanity checks

End.

