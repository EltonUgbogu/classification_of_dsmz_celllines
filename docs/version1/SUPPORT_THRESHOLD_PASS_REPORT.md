# Support-threshold reproducibility pass — report

Date: 2026-05-31 (UTC)
Convention: timestamped `.bak_support_threshold_20260531` copies (no git in tree).
Scope: make the majority-style support threshold `m = max(2, ceil(|R| / 2))` explicit per cohort; remove the silent fallback; make missing thresholds fail-fast.

---

## 1. Files edited

| Path | Diff size |
|---|---:|
| `config/config.yaml` | 52 lines |
| `Snakefile` | 63 lines |

Diffs at repo root: `SUPPORT_THRESHOLD_PASS_diff_config.patch`, `SUPPORT_THRESHOLD_PASS_diff_Snakefile.patch`.

---

## 2. Backups

| Path | md5 |
|---|---|
| `config/config.bak_support_threshold_20260531.yaml` | `f068ab7a03425df84fb72f064fd123eb` |
| `Snakefile.bak_support_threshold_20260531` | `1d24a0ae0ec68ecbd13dcd559f18e1a6` |

---

## 3. Exact config keys and values added

Added under each cohort's `profiles.<cohort>.tumour_neighbourhoods` block in `config/config.yaml`, with a comment line stating the formula:

```yaml
profiles:
  multicohort_cancer:
    tumour_neighbourhoods:
      directions: [...]                                 # |R| = 18
      similarity_consensus_min_support: 9

  brca:
    tumour_neighbourhoods:
      directions: [..., pam50_euc, pam50_corr]          # |R| = 20 (incl. pam50)
      similarity_consensus_min_support: 10

  nbl:
    tumour_neighbourhoods:
      directions: [...]                                 # |R| = 18
      similarity_consensus_min_support: 9

  rbl:
    tumour_neighbourhoods:
      directions: [...]                                 # |R| = 18
      similarity_consensus_min_support: 9

  heme:
    tumour_neighbourhoods:
      directions: [...]                                 # |R| = 18
      similarity_consensus_min_support: 9
```

The three brief-mandated values (`brca=10, nbl=9, rbl=9`) plus matching values for `multicohort_cancer` and `heme` (so the fail-fast logic does not surprise those profiles on next use). Each insertion is preceded by a two-line comment naming the formula and reminding the maintainer that there is no silent fallback.

---

## 4. Exact Snakefile changes

### (a) New module-level helper + constant (after `CONS_EUC_PATTERN`, around L497)

```python
def _resolve_majority_support_threshold():
    """
    Resolves the majority-style support threshold m = max(2, ceil(|R| / 2))
    for the support-threshold consensus cell-line similarity network, where
    |R| is the number of representation-specific graphs for this cohort.

    The value MUST be set explicitly per cohort under
        profiles.<cohort>.tumour_neighbourhoods.similarity_consensus_min_support
    in config/config.yaml. There is no silent fallback — running with the
    threshold missing fails at Snakefile load time with a clear message.
    """
    nh_cfg = cfg.get("tumour_neighbourhoods", {})
    if "similarity_consensus_min_support" not in nh_cfg:
        recommended = max(2, (len(CONS_DIRECTIONS) + 1) // 2)
        raise ValueError(
            "\n[Snakefile] Missing required config key for profile "
            f"'{profile_name}': "
            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support.\n"
            "This is the majority-style support threshold "
            "(m = max(2, ceil(|R| / 2))). "
            f"For this cohort |R| = {len(CONS_DIRECTIONS)}, so the recommended "
            f"value is {recommended}. Add the key explicitly under the cohort's "
            "tumour_neighbourhoods block in config/config.yaml; no silent "
            "fallback is permitted.\n"
        )
    m = nh_cfg["similarity_consensus_min_support"]
    if not isinstance(m, int) or m < 2:
        raise ValueError(
            "\n[Snakefile] Invalid value for "
            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support: "
            f"{m!r}. The majority-style support threshold must be an integer >= 2.\n"
        )
    return m


SIMILARITY_CONSENSUS_MIN_SUPPORT = _resolve_majority_support_threshold()
```

### (b) Rule `build_patient_referenced_support_threshold_consensus_cell_line_similarity_network`'s `params:` block (line ~1946)

Before:
```python
        min_support = cfg.get(
            "tumour_neighbourhoods", {}
        ).get("similarity_consensus_min_support", max(2, (len(CONS_DIRECTIONS) + 1) // 2))
```

After:
```python
        # Majority-style support threshold m = max(2, ceil(|R| / 2)); resolved
        # at Snakefile load time from the explicit per-cohort config key
        # tumour_neighbourhoods.similarity_consensus_min_support (no fallback —
        # missing values fail with a clear error; see
        # _resolve_majority_support_threshold above).
        min_support = SIMILARITY_CONSENSUS_MIN_SUPPORT
```

No other Snakefile lines touched. Rule names, inputs, outputs, log paths, conda envs, shell bodies, and the `build_consensus_from_direction_edgefiles.py` script itself are unchanged.

---

## 5. Hidden fallback removal — confirmation

Confirmed by `grep`:

```
$ grep -n "max(2, (len(CONS_DIRECTIONS)" Snakefile
$ grep -n "similarity_consensus_min_support" Snakefile
514:    if "similarity_consensus_min_support" not in nh_cfg:
517:            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support.\n"
525:    m = nh_cfg["similarity_consensus_min_support"]
529:            f"profiles.{profile_name}.tumour_neighbourhoods.similarity_consensus_min_support: "
```

The only remaining references are inside the fail-fast helper. The previous fallback expression `max(2, (len(CONS_DIRECTIONS) + 1) // 2)` is gone from the build rule's `params:` block; it survives only as the *recommended value* reported in the error message, never as a silently-applied default.

---

## 6. Missing-threshold failure — confirmation

Test procedure:

1. Backed up `config/config.yaml` to `config/config.yaml.SAFE_LIVE_$$` (PID-suffixed).
2. Replaced live config with a copy where the **rbl** block's `similarity_consensus_min_support: 9` (and its two leading comment lines) had been removed.
3. Ran `snakemake -n --config pipeline_profile=rbl <support-threshold consensus target>`.
4. EXIT-trap restored the original config from the safe backup, removed the temp config.

Result:

```
Building DAG of jobs...
ValueError in file "<repo>/Snakefile", line 514:

[Snakefile] Missing required config key for profile 'rbl':
profiles.rbl.tumour_neighbourhoods.similarity_consensus_min_support.
This is the majority-style support threshold (m = max(2, ceil(|R| / 2))).
For this cohort |R| = 18, so the recommended value is 9. Add the key
explicitly under the cohort's tumour_neighbourhoods block in config/config.yaml;
no silent fallback is permitted.

  File "<repo>/Snakefile", line 535, in <module>
  File "<repo>/Snakefile", line 514, in _resolve_majority_support_threshold

snakemake exit code: 1
```

Live config restored verbatim by the trap; post-test md5 of `config/config.yaml` equals the post-edit md5 (`caf4cd51bf0ffa8344029e5a87ff8b5b`). No temp artefact left in the project tree.

---

## 7. Baseline checksums for support-threshold edge TSVs (pre-edit)

```
37679268f4c2905b1d275a8d4b57c4e2  results/unsupervised/brca/.../patient_referenced_support_threshold_consensus_cell_line_similarity_edges.tsv  (size 7514)
8d502e630e80f0c413671450cdc5f72a  results/unsupervised/nbl/.../patient_referenced_support_threshold_consensus_cell_line_similarity_edges.tsv  (size 3014)
2547508135267c623509ce931dfdf66d  results/unsupervised/rbl/.../patient_referenced_support_threshold_consensus_cell_line_similarity_edges.tsv  (size 888)
```

---

## 8. Post-edit checksums for support-threshold edge TSVs

```
37679268f4c2905b1d275a8d4b57c4e2  brca
8d502e630e80f0c413671450cdc5f72a  nbl
2547508135267c623509ce931dfdf66d  rbl
```

Plus independent content verification — invoked `scripts/build_consensus_from_direction_edgefiles.py` directly with the new explicit `--min_support 10/9/9` arguments, writing to `scratch/st_verify/<cohort>_rebuild.tsv`. md5s of those independent rebuilds:

```
37679268f4c2905b1d275a8d4b57c4e2  brca rebuild
8d502e630e80f0c413671450cdc5f72a  nbl  rebuild
2547508135267c623509ce931dfdf66d  rbl  rebuild
```

The build script reported, for each cohort:
```
BRCA: [OK] wrote scratch/st_verify/brca_rebuild.tsv with 31 edges from 20 directions
NBL:  [OK] wrote scratch/st_verify/nbl_rebuild.tsv  with 14 edges from 18 directions
RBL:  [OK] wrote scratch/st_verify/rbl_rebuild.tsv  with  4 edges from 18 directions
```

These edge counts match every prior pass's counts. Temp rebuilds removed after verification.

---

## 9. Edge TSVs identical before and after?

**Yes — byte-identical for all three cohorts** (the post-edit md5s, the on-disk on-jarvis md5s, and the independently-rebuilt md5s all agree). The explicit values exactly reproduce the previously-implicit fallback.

---

## 10. Plots regenerated?

**No.** No plotting rule was re-run. `snakemake --touch` was used after the edit to reconcile mtimes between the rule outputs and the edited `config/config.yaml` (which is an `input:` to `feature_selection_unsupervised`), but `--touch` updates mtimes without re-executing rule scripts. See §11 for the corresponding content check.

---

## 11. Plot graph content changed?

**No.** md5 comparison of the 9 consensus figure files (3 cohorts × {pdf, svg, png}) and a spot check of the 3 resolved PDFs against the pre-edit baseline:

| cohort | consensus.pdf | consensus.svg | consensus.png | resolved.pdf |
|---|---|---|---|---|
| BRCA | MATCH | MATCH | MATCH | MATCH |
| NBL  | MATCH | MATCH | MATCH | MATCH |
| RBL  | MATCH | MATCH | MATCH | MATCH |

All 15 files byte-identical to baseline. No graph content can have changed, because no rule re-executed.

---

## 12. Resolved-graph plots affected?

**No, and they don't depend on support annotations either.**

- Dependency check: `rule plot_patient_referenced_resolved_cell_line_neighbourhood_graph` reads only `resolved_dsmz_neighbours.tsv` (and indirectly the per-direction `cell_line_similarity_graph_edges_<dir>.tsv` files via the script's directory scan); it does NOT read `P_CONS_SIMILARITY_CONSENSUS_EDGES_TSV` or any support-annotation TSV.
- Initial dry-run (before `--touch`) showed `plot_patient_referenced_resolved_cell_line_neighbourhood_graph` in the cascade for all three cohorts, but only because its input chain reached back to `feature_selection_unsupervised` whose `cfg_file = CFGFILE_ABS` input got a new mtime when I edited `config/config.yaml`. The cascade was an mtime artefact — no content dependency on the threshold change.
- After `--touch`, the resolved plots are mtime-current and unchanged in content (md5 matches baseline — §11).

---

## 13. Dry-run command and result

Command (per cohort; ran for brca/nbl/rbl):
```bash
snakemake --snakefile Snakefile \
  --config pipeline_profile=<cohort> \
  --rerun-triggers mtime \
  --cores 1 -n \
  results/unsupervised/<cohort>/.../plots/Fig_<COHORT>_patient_referenced_support_threshold_consensus_cell_line_similarity_network.pdf \
  results/unsupervised/<cohort>/.../plots/patient_referenced_support_threshold_consensus_cell_line_similarity_edges.tsv
```

Result before `--touch`: cascade of 257 / 232 / 232 jobs (BRCA / NBL / RBL), root cause `feature_selection_unsupervised` triggered by the edited `config/config.yaml` mtime.

Result after `--touch`: **"Nothing to be done"** for all three cohorts. See §15.

---

## 14. Real Snakemake command and result

Used `snakemake --touch` (not a build):

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml \
  --config pipeline_profile=<cohort> --cores 1 --touch \
  <support-threshold + resolved targets + cohort-level aggregated node-stats>
```

- BRCA: 258 / 258 jobs marked current.
- NBL:  232 / 232 jobs marked current.
- RBL:  233 / 233 jobs marked current.

No file contents rewritten (verified by md5 comparison against baseline — see §11). No new figures, no new TSVs. Side effect: mtimes of touched outputs are now newer than the edited `config/config.yaml`, so future dry-runs do not see a stale chain.

---

## 15. Final dry-run result

```
$ snakemake -n --config pipeline_profile=brca --rerun-triggers mtime \
    <support-threshold and resolved targets>
Building DAG of jobs...
Nothing to be done (all requested files are present and up to date).

$ snakemake -n --config pipeline_profile=nbl ...
Nothing to be done (all requested files are present and up to date).

$ snakemake -n --config pipeline_profile=rbl ...
Nothing to be done (all requested files are present and up to date).
```

(NBL's run additionally notes `1 jobs have missing provenance/metadata so that it in part cannot be used to trigger re-runs. Rules with missing metadata: split_profile_joint_vst` — pre-existing, unrelated to this pass, present in every previous pass's final dry-run.)

---

## 16. Scientific logic unchanged

Confirmed:

- **Representation-specific graph construction** (`cell_line_similarity_graph` rule, `scripts/compute_cell_line_similarity.R`, `scripts/run_tumour_neighbourhood_analysis.R`, `R/base_functions/tumour_neighbourhood.R`): untouched.
- **Resolved graph logic** (`resolve_dsmz_graph_neighbours` rule + `scripts/resolve_dsmz_graph_neighbours.R`; `N_final(i) = N_global(i) ∩ N_local(i)`): untouched.
- **DESeq2 focal-group selection** (`derive_anchor_list`, `derive_isolate_list`, `prepare_deseq2_inputs`, `add_component_to_metadata`): inputs follow the renamed constants from the prior patient_referenced pass; column names unchanged; this pass touches none of the DESeq2 rules or their scripts.
- **Node centrality calculations** (in `scripts/visualize_resolved_dsmz_graph.py`'s `pack_non_isolate_components`): untouched.
- **Plotting aesthetics** (`scripts/graph_plot_style.py`, `scripts/plot_consensus_graph.py`, `scripts/visualize_resolved_dsmz_graph.py`): untouched.
- **Thesis files** (`methods.tex`, `methods.md`, `algorithm.md`, `README.md`, `AGENT.md`, `task.md`): untouched (`stat` mtime unchanged).
- **Patient-referenced naming** and **pan-cancer naming**: untouched (Snakefile constants unchanged from prior pass).
- **Generated figure/TSV output names**: identical to prior pass.
- **Graph topology** (nodes, edges, support counts, isolates, components, anchors): identical (md5 confirmation, §11).
- **`scripts/build_consensus_from_direction_edgefiles.py`**: untouched.

---

## 17. Unresolved issues requiring approval

None. The pass is complete, content-preserving, and the fail-safe is in place.

Two adjacent observations (not blockers; flagged in the spirit of the brief's "any unresolved issue" clause):

1. **`config/config.yaml` mtime cascade.** Editing `config.yaml` triggers six rules that declare it as `input:` (including `feature_selection_unsupervised`), which cascades to ~230 jobs per cohort. This pass handled it with `snakemake --touch`, which is the standard Snakemake idiom for "I know the content is the same — please update mtimes". Each future config edit that does not change rule semantics will hit the same cascade and the same `--touch` workaround. Not a problem to solve in this pass; just a flag.

2. **Older `.bak` files in `config/`.** I noticed (during the SAFE_LIVE check) that `config/` contains six older `.bak*` files from prior agent passes (`config.yaml.bak`, `config.yaml.bak_20260503_062313`, `…_064528`, `…_codex_20260510_…`, `…_codex_20260513_…`, `…_no_heme_multicohort_20260510`). They are not actively referenced. A future housekeeping pass could move them to `archive/` analogous to last session's root-level cleanup. Not touched this pass.

---

## Verification commands (post-edit, reproducible)

```bash
# 1. Confirm fallback gone:
grep -n "max(2, (len(CONS_DIRECTIONS)" Snakefile          # → empty
grep -nE "similarity_consensus_min_support:\s*[0-9]+" config/config.yaml
#   → 5 lines (one per cohort) with values 9/10/9/9/9

# 2. Sanity dry-run (no --set-threads etc.):
$SMK --snakefile Snakefile --config pipeline_profile=rbl \
  --rerun-triggers mtime --cores 1 -n \
  results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/plots/patient_referenced_support_threshold_consensus_cell_line_similarity_edges.tsv
#   → "Nothing to be done"

# 3. Reproduce the failure mode (DOES NOT modify the live config):
# (See §6 — the backup-replace-restore script with a trap is in
#  the turn's command history but does not need to be re-run.)
```
