# DESeq2 Re-run Recommendation: BRCA and RBL

**Date:** 2025-01-11  
**Status:** Recommendation to re-run for consistency

---

## Current Status

### NBL ✅
- **Status:** Re-run with improved script (Jan 11 05:03)
- **Improvements applied:** ✅ All
- **Contrasts:** 7 (3 isolates + 4 anchors)
- **Unique features:** 20 genes (recurrence ≥ 2)

### BRCA ⚠️
- **Status:** Old run (Jan 11 04:36) - before improvements
- **Improvements applied:** ❌ None
- **Contrasts:** 8 (3 isolates + 5 anchors)
- **Unique features:** 261 genes (recurrence ≥ 2)

### RBL ⚠️
- **Status:** Old run (Jan 11 04:39) - before improvements
- **Improvements applied:** ❌ None
- **Contrasts:** 6 (4 isolates + 2 anchors)
- **Unique features:** 22 genes (recurrence ≥ 2)

---

## Why Re-run BRCA and RBL?

### 1. **Consistency Across Diseases**
- NBL now uses improved script with:
  - Expression filtering (≥10 normalized counts)
  - Better DESeq2 parameters
  - Enhanced marker ranking
- BRCA and RBL should use the same methodology for fair comparison

### 2. **Improved Statistical Rigor**
- **Expression filtering:** Prevents selection of genes that are "low everywhere but stochastic"
- **Better DESeq2 parameters:** `fitType="local"`, `minReplicatesForReplace=Inf`
- **LFC shrinkage:** Reduces noise in effect size estimates (if apeglm available)

### 3. **Enhanced Reporting**
- Reports recurrence ≥ 1 (for context)
- Reports recurrence ≥ k (strict threshold)
- Better logging and diagnostics

### 4. **Expected Changes**
- **Marker counts may decrease** (due to expression filter)
- **Marker quality should improve** (better ranking, less noise)
- **Unique feature sets may change** (more stringent filtering)

---

## Impact Assessment

### BRCA (Current: 261 unique genes)
- **Expected:** May decrease due to expression filtering
- **Reason:** More stringent criteria will filter out low-expression markers
- **Benefit:** Higher quality, more reliable markers

### RBL (Current: 22 unique genes)
- **Expected:** May decrease slightly or stay similar
- **Reason:** Already has strict filtering; expression filter may remove a few more
- **Benefit:** Ensures consistency with NBL methodology

---

## Re-run Commands

### BRCA
```bash
ISO_BRC=$(paste -sd, results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/isolate_list.txt)
ANCHOR_BRC=$(cut -f1 results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv | tail -n +2 | paste -sd,)

nohup /work/ugbogu/.conda/envs/tcga-r-env/bin/Rscript scripts/deseq2_isolate_degs.R \
  --counts results/unsupervised/brca/deseq2_inputs/counts.tsv \
  --meta results/unsupervised/brca/deseq2_inputs/metadata_with_components.tsv \
  --sample_id_col sample_id \
  --cell_line_col cell_line \
  --component_col component \
  --isolate_list "$ISO_BRC" \
  --anchor_list "$ANCHOR_BRC" \
  --anchor_components results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv \
  --outdir results/unsupervised/brca/deseq2_markers \
  --fdr_isolate 0.01 --lfc_isolate 1.5 --topN_isolate 50 \
  --fdr_anchor 0.05 --lfc_anchor 1.0 --topN_anchor 200 \
  --recurrence_k 2 > results/unsupervised/brca/deseq2_markers/deseq2_run.log 2>&1 &
```

### RBL
```bash
ISO_RBL=$(paste -sd, results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/isolate_list.txt)
ANCHOR_RBL=$(cut -f1 results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv | tail -n +2 | paste -sd,)

nohup /work/ugbogu/.conda/envs/tcga-r-env/bin/Rscript scripts/deseq2_isolate_degs.R \
  --counts results/unsupervised/rbl/deseq2_inputs/counts.tsv \
  --meta results/unsupervised/rbl/deseq2_inputs/metadata_with_components.tsv \
  --sample_id_col sample_id \
  --cell_line_col cell_line \
  --component_col component \
  --isolate_list "$ISO_RBL" \
  --anchor_list "$ANCHOR_RBL" \
  --anchor_components results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv \
  --outdir results/unsupervised/rbl/deseq2_markers \
  --fdr_isolate 0.01 --lfc_isolate 1.5 --topN_isolate 50 \
  --fdr_anchor 0.05 --lfc_anchor 1.0 --topN_anchor 200 \
  --recurrence_k 2 > results/unsupervised/rbl/deseq2_markers/deseq2_run.log 2>&1 &
```

---

## Recommendation

**✅ Re-run both BRCA and RBL** to ensure:
1. Consistent methodology across all three diseases
2. Application of improved statistical methods
3. Fair comparison in pan-cancer analysis
4. Higher quality marker gene sets

---

## Expected Timeline

- **BRCA:** ~5-10 minutes (8 contrasts)
- **RBL:** ~3-5 minutes (6 contrasts)
- **Total:** ~10-15 minutes for both

---

## Monitoring

After starting the runs, monitor with:
```bash
# Check BRCA progress
tail -f results/unsupervised/brca/deseq2_markers/deseq2_run.log

# Check RBL progress
tail -f results/unsupervised/rbl/deseq2_markers/deseq2_run.log
```
