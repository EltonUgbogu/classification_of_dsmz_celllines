# Why "MissingRuleException" for DSMZ_DSMZ_graph_edges_*.tsv

## Diagnostic results

### 1) Rule is **not** in the workflow Snakemake loads for RBL

```text
snakemake --configfile config/config.yaml --config pipeline_profile=rbl --list-rules | grep "dsmz_dsmz_similarity_graph"
→ RULE NOT LOADED (no match)
```

So for `pipeline_profile=rbl`, Snakemake does **not** list `dsmz_dsmz_similarity_graph`. That’s why requesting any `DSMZ_DSMZ_graph_edges_*.tsv` file gives **MissingRuleException**: no rule in the loaded workflow has that output pattern.

### 2) The rule **is** in the Snakefile

- **Rule:** `dsmz_dsmz_similarity_graph` (Snakefile, around line 1106).
- **Output pattern:**  
  `{TUMOUR_NH_ROOT}/{direction}/final_consensus/DSMZ_DSMZ_graph_edges_{direction}.tsv`
- **Wildcard constraint:** `direction=r"[a-zA-Z0-9_]+"` (already relaxed).

So the rule and path template are correct in the file; the issue is that this rule is **not active** when you run with `--config pipeline_profile=rbl`.

### 3) How the rule is used

- `rule pipeline_targets` (line ~1850) has in its `input:`  
  `expand(rules.dsmz_dsmz_similarity_graph.output.edges_tsv, direction=AGN_DIRECTIONS)` (and other outputs of the same rule).
- So if the Snakefile parses fully, `dsmz_dsmz_similarity_graph` is part of the same workflow as `pipeline_targets` and `all`.

So either:

- Parsing stops or fails **before** the rule is registered (e.g. exception when evaluating config/AGN_DIRECTIONS/inputs for RBL), or  
- Something (e.g. `--list-rules` behaviour or a different Snakefile) means the rule never appears in the loaded workflow for RBL.

---

## What to do next

### A. Confirm whether the rule is loaded (your run)

From pipeline root:

```bash
cd classification_of_dsmz_celllines

# Full list of rules (no grep) – check if dsmz_dsmz_similarity_graph appears at all
mkdir -p scratch
snakemake --configfile config/config.yaml --config pipeline_profile=rbl --list-rules 2>&1 | tee scratch/rbl_rules.txt
grep -i "dsmz\|similarity_graph" scratch/rbl_rules.txt || echo "RULE NOT IN LIST"
```

If you still see no `dsmz_dsmz_similarity_graph`, the rule is not in the loaded rule set for RBL.

### B. See if the DAG can be built (and if the rule appears there)

```bash
cd classification_of_dsmz_celllines

# Dry-run default target; look for dsmz_dsmz_similarity_graph in the job list
mkdir -p scratch
snakemake --configfile config/config.yaml --config pipeline_profile=rbl -n -p all 2>&1 | tee scratch/rbl_dag.txt
grep -i "dsmz_dsmz_similarity_graph\|DSMZ_DSMZ_graph_edges" scratch/rbl_dag.txt || echo "NOT IN DAG"
```

- If the rule **appears** here, the rule is loaded and the problem is only with how you’re requesting the edge files (e.g. path mismatch).
- If the DAG build **fails** (exception or error), the message will point to what breaks for RBL (e.g. config, `AGN_DIRECTIONS`, or an input function).
- If the rule **does not appear** and there’s no error, then the rule is not part of the workflow Snakemake builds for RBL (same as “rule not loaded”).

### C. Force Snakemake to consider only this rule (sanity check)

```bash
cd classification_of_dsmz_celllines

snakemake --configfile config/config.yaml --config pipeline_profile=rbl -n -p \
  --allowed-rules dsmz_dsmz_similarity_graph \
  results/unsupervised/rbl/tumour_neighbourhoods/Spearman_euc/final_consensus/DSMZ_DSMZ_graph_edges_Spearman_euc.tsv
```

- If this **succeeds** (no MissingRuleException), the rule and path are correct and the issue is that the rule is otherwise excluded from the normal DAG for RBL.
- If this still gives **MissingRuleException**, then even with `--allowed-rules`, no rule matches that path in the loaded Snakefile (e.g. wrong Snakefile or path still doesn’t match).

### D. Check for parse/execution errors with RBL config

```bash
cd classification_of_dsmz_celllines
snakemake --configfile config/config.yaml --config pipeline_profile=rbl --debug 2>&1 | head -150
```

Look for Python tracebacks or config/validation errors that could prevent the Snakefile from loading the rule.

---

## Likely cause and fix

**Likely cause:** For `pipeline_profile=rbl`, something during **parsing or DAG building** fails or is skipped (e.g. exception when evaluating `AGN_DIRECTIONS`, `pipeline_targets` inputs, or an input function), so Snakemake never registers or uses `dsmz_dsmz_similarity_graph`. That would explain both “RULE NOT LOADED” and MissingRuleException.

**Concrete checks:**

1. Run **B** and **D** and fix any **config/validation/input** error that appears for RBL (e.g. missing key, wrong type, or file path).
2. Ensure **RBL’s** `tumour_neighbourhoods.directions` in `config/config.yaml` is valid and that `AGN_DIRECTIONS` is non-empty and contains the directions you expect (e.g. `Spearman_euc`).
3. After the workflow loads cleanly, run **C**; if that works, then request all 17 edge TSVs again (or use `rebuild_rbl_direction_edges.sh`).

Once the rule is loaded and the DAG builds for RBL, the existing output pattern and rebuild script are correct; no path or rule-name guess is needed.
