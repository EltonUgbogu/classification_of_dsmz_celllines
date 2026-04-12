# Cell Line Canonicalization Fix Summary

**Date:** 2025-01-11  
**Issue:** LAN_5 and other LAN lines not matching in NBL metadata  
**Status:** ✅ Fixed

---

## Problem

NBL node stats used cell line names like `LAN_5`, but DSMZ metadata had normalized names like `LAN5` (no underscore). The prep script couldn't match them, resulting in:

- **Missing from metadata:** LAN_1, LAN_2, LAN_5, LAN_6
- **Anchor LAN_5 skipped** in DESeq2 analysis

---

## Root Cause

Three different naming conventions for the same biological cell line:

1. **Node stats:** `LAN_5` (underscore)
2. **DSMZ_Cell_line_norm:** `LAN5` (no separator)
3. **Sample ID:** `NG-A0611_LAN5_libLAB6087` (embedded in sample name)

The original matching logic only handled exact matches or simple regex patterns, missing the canonicalization step.

---

## Solution: Canonicalization Function

Added a `canon_cl()` function that:

- Converts to uppercase
- Removes all non-alphanumeric characters (underscores, hyphens, spaces, etc.)

**Examples:**
- `LAN_5` → `LAN5`
- `LAN-5` → `LAN5`
- `LAN5` → `LAN5`
- `CHP_126` → `CHP126`
- `CHP-126` → `CHP126`

---

## Implementation

### 1. Canonicalization Function

```r
canon_cl <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]", "", x)   # remove underscores, hyphens, spaces, etc.
  x
}
```

### 2. Metadata Canonicalization

```r
meta_df <- meta_df %>%
  mutate(
    cell_line_clean = ifelse(!is.na(DSMZ_Cell_line_norm) & DSMZ_Cell_line_norm != "", 
                             DSMZ_Cell_line_norm, 
                             Cell_Line),
    cell_line_canon = canon_cl(cell_line_clean),
    Cell_Line_canon = canon_cl(Cell_Line)
  )
```

### 3. Target Cell Line Canonicalization

```r
target_cell_lines <- unique(node_stats$cell_line)
target_canon <- canon_cl(target_cell_lines)
```

### 4. Matching by Canonical Key

```r
for (i in seq_along(target_cell_lines)) {
  cl <- target_cell_lines[i]
  cl_can <- target_canon[i]
  
  # Match by canonicalized name
  meta_matches <- meta_df %>%
    filter(cell_line_canon == cl_can | Cell_Line_canon == cl_can)
  # ... rest of matching logic
}
```

### 5. Unmatched Cell Line Reporting

```r
unmatched <- target_cell_lines[!(target_cell_lines %in% matched_cell_lines)]
if (length(unmatched) > 0) {
  message("[WARN] Unmatched cell lines from node stats:")
  message(paste("  ", unmatched, collapse = "\n"))
}
```

---

## Results

### Before Fix
- **Matched:** 14/18 cell lines
- **Missing:** LAN_1, LAN_2, LAN_5, LAN_6
- **Anchor LAN_5:** Skipped in DESeq2

### After Fix
- **Matched:** 18/18 cell lines ✅
- **All LAN lines present:** LAN_1, LAN_2, LAN_5, LAN_6 ✅
- **Anchor LAN_5:** Now available for DESeq2 analysis ✅

---

## Verification

```bash
# All node stats cell lines now in metadata
$ cut -f2 nbl/deseq2_inputs/metadata.tsv | tail -n +2 | sort -u
CCLF_PEDS_0051_T
CHP_126
CHP_134
GI_ME_N
IMR_32
KELLY
LAN_1      ← Now present
LAN_2      ← Now present
LAN_5      ← Now present
LAN_6      ← Now present
LS
MHH_NB_11
NBL_S
NGP
NMB
SH_SY5Y
SIMA
SK_N_BE_2
```

---

## Benefits

1. **Robust matching:** Handles all naming variations (underscore, hyphen, no separator)
2. **Preserves graph naming:** Uses node-stats labels (`LAN_5`) as `cell_line` in metadata
3. **Future-proof:** Fixes similar issues for other cell lines (e.g., CHP-126 vs CHP_126)
4. **Transparent:** Reports unmatched cell lines for debugging

---

## Next Steps

1. ✅ Re-run NBL prep script → All 18 cell lines matched
2. ✅ Update metadata with components → LAN_5 has component assignment
3. ⏭️ Re-run NBL DESeq2 → LAN_5 anchor will now be processed

---

## Files Modified

- `/work/ugbogu/pipeline/scripts/prepare_deseq2_inputs.R`
  - Added `canon_cl()` function
  - Updated matching logic to use canonicalization
  - Added unmatched cell line reporting

---

## Design Rationale

**Why preserve node-stats naming (`LAN_5`) as `cell_line`?**

- Keeps consistency with graph analysis outputs
- Isolate/anchor lists use node-stats names
- Downstream scripts expect these names
- Only the matching step uses canonicalization (internal)

**Why canonicalize for matching only?**

- Preserves original naming conventions in outputs
- Makes matching robust without changing data representation
- Allows different naming conventions to coexist
