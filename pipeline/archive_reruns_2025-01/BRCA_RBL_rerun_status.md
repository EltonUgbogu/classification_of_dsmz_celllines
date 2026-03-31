# BRCA and RBL DESeq2 Re-run Status

**Date:** 2025-01-11  
**Status:** 🔄 Running in background

---

## Processes Started

Both BRCA and RBL DESeq2 analyses have been started in the background with the improved script.

---

## Monitor Progress

### BRCA
```bash
# View log
tail -f results/unsupervised/brca/deseq2_markers/deseq2_run.log

# Check if still running
ps aux | grep "deseq2.*brca" | grep -v grep
```

### RBL
```bash
# View log
tail -f results/unsupervised/rbl/deseq2_markers/deseq2_run.log

# Check if still running
ps aux | grep "deseq2.*rbl" | grep -v grep
```

---

## Expected Improvements Applied

1. ✅ **Expression filtering:** ≥10 normalized counts in test sample
2. ✅ **Better DESeq2 parameters:** `fitType="local"`, `minReplicatesForReplace=Inf`
3. ✅ **LFC shrinkage attempt:** apeglm if available
4. ✅ **Enhanced reporting:** Recurrence ≥ 1 and ≥ k

---

## Expected Results

### BRCA
- **Contrasts:** 8 (3 isolates + 5 anchors)
- **Previous unique features:** 261 genes
- **Expected:** May decrease due to expression filtering (higher quality)

### RBL
- **Contrasts:** 6 (4 isolates + 2 anchors)
- **Previous unique features:** 22 genes
- **Expected:** May decrease slightly or stay similar

---

## Check Completion

```bash
# Check if processes completed
ps aux | grep -E "deseq2.*brca|deseq2.*rbl" | grep -v grep

# If no processes, check final results
cat results/unsupervised/brca/deseq2_markers/markers/marker_sets_manifest.tsv
cat results/unsupervised/rbl/deseq2_markers/markers/marker_sets_manifest.tsv
```

---

## Timeline

- **BRCA:** ~5-10 minutes (8 contrasts)
- **RBL:** ~3-5 minutes (6 contrasts)
- **Total:** ~10-15 minutes for both
