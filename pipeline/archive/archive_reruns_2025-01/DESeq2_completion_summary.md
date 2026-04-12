# DESeq2 Analysis Completion Summary

**Date:** 2025-01-11  
**Status:** ✅ Complete for BRCA, NBL, and RBL

---

## Overview

Successfully completed DESeq2 differential expression analysis for anchor and isolate cell lines across three cancer types. The analysis identified marker genes using:

- **Isolate DEGs**: Each isolate vs all other cell lines (strict thresholds: FDR ≤ 0.01, |log2FC| ≥ 1.5, top 50 genes)
- **Anchor DEGs**: Each anchor vs cell lines outside its component (moderate thresholds: FDR ≤ 0.05, |log2FC| ≥ 1.0, top 200 genes)
- **Unique feature sets**: Recurrence-filtered genes appearing in ≥2 contrasts

---

## Results Summary

### 🟣 BRCA (Breast Cancer)

**Inputs:**
- 29 cell lines
- 5 anchors: BT_474, EFM_192C (Component 0), HDQ_P1 (Component 1), COLO_824 (Component 3), ETCC_006 (Component 5)
- 3 isolates: CAL_51, DU_4475, HCC_1937

**Outputs:**
- **8 contrasts** (3 isolate + 5 anchor)
- **Unique feature set: 261 genes** (recurrence ≥ 2)
- Location: `results/unsupervised/brca/deseq2_markers/`

**Marker counts per contrast:**
- isolate_CAL_51_vs_rest: 50 markers
- isolate_DU_4475_vs_rest: 50 markers
- isolate_HCC_1937_vs_rest: 50 markers
- anchor_BT_474_vs_outside_component_0: 200 markers
- anchor_EFM_192C_vs_outside_component_0: 200 markers
- anchor_HDQ_P1_vs_outside_component_1: 200 markers
- anchor_COLO_824_vs_outside_component_3: 200 markers
- anchor_ETCC_006_vs_outside_component_5: 200 markers

---

### 🔵 NBL (Neuroblastoma)

**Inputs:**
- 14 cell lines (matched from 18 target)
- 4 anchors: IMR_32, CHP_126 (Component 1), CCLF_PEDS_0051_T (Component 0), LAN_5 (Component 2)
- 3 isolates: LS, SH_SY5Y, SK_N_BE_2

**Outputs:**
- **6 contrasts** (3 isolate + 3 anchor; LAN_5 not found in data)
- **Unique feature set: 153 genes** (recurrence ≥ 2)
- Location: `results/unsupervised/nbl/deseq2_markers/`

**Marker counts per contrast:**
- isolate_LS_vs_rest: 50 markers
- isolate_SH_SY5Y_vs_rest: 39 markers
- isolate_SK_N_BE_2_vs_rest: 50 markers
- anchor_IMR_32_vs_outside_component_1: 200 markers
- anchor_CHP_126_vs_outside_component_1: 200 markers
- anchor_CCLF_PEDS_0051_T_vs_outside_component_0: 200 markers

**Note:** LAN_5 anchor was not found in the matched samples (likely a naming mismatch).

---

### 🟠 RBL (Retinoblastoma)

**Inputs:**
- 11 cell lines
- 2 anchors: NG-30919_RBL_18_lib626623_10098_1 (Component 0), NG-30919_RBL_15_lib626622_10098_2 (Component 2)
- 4 isolates: NG-30919_RBL_14_lib628470_10098_1, NG-30919_RBL_7_lib628468_10098_1, NG-30919_WERI_RB1_lib628472_10098_1, NG-30919_Y_79_lib626626_10098_1

**Outputs:**
- **6 contrasts** (4 isolate + 2 anchor)
- **Unique feature set: 22 genes** (recurrence ≥ 2)
- Location: `results/unsupervised/rbl/deseq2_markers/`

**Marker counts per contrast:**
- isolate_NG-30919_RBL_14_lib628470_10098_1_vs_rest: 50 markers
- isolate_NG-30919_RBL_7_lib628468_10098_1_vs_rest: 50 markers
- isolate_NG-30919_WERI_RB1_lib628472_10098_1_vs_rest: 50 markers
- isolate_NG-30919_Y_79_lib626626_10098_1_vs_rest: 50 markers
- anchor_NG-30919_RBL_18_lib626623_10098_1_vs_outside_component_0: 88 markers
- anchor_NG-30919_RBL_15_lib626622_10098_2_vs_outside_component_2: 200 markers

---

## Technical Notes

### Single-Sample-Per-Cell-Line Handling

The analysis successfully handled the constraint of **one sample per cell line** by:

1. Using **binary grouping designs** (isolate vs REST, anchor vs OUTSIDE_COMP)
2. Applying **local dispersion estimation** for single-sample-per-group contrasts
3. Maintaining proper size factor estimation across all samples

### File Structure

For each disease (`{profile}` = brca, nbl, rbl):

```
results/unsupervised/{profile}/deseq2_markers/
├── tables/
│   ├── isolate_{cell_line}_vs_rest.tsv
│   └── anchor_{cell_line}_vs_outside_component_{comp}.tsv
├── markers/
│   ├── isolate_{cell_line}_vs_rest_markers_top50.txt
│   ├── anchor_{cell_line}_vs_outside_component_{comp}_markers_top200.txt
│   ├── gene_recurrence_across_contrasts.tsv
│   ├── unique_feature_set_recurrence_ge_2.txt
│   └── marker_sets_manifest.tsv
└── qc/
    └── size_factors.tsv
```

---

## Next Steps

1. **Review marker gene lists** for biological relevance
2. **Tune recurrence threshold** (`recurrence_k`) if needed (see `gene_recurrence_across_contrasts.tsv`)
3. **Combine pan-cancer feature sets** across BRCA, NBL, and RBL
4. **Apply to tumour data** for signature scoring and classification

---

## Scripts Used

1. `scripts/prepare_deseq2_inputs.R` - Extract disease-specific counts and metadata
2. `scripts/add_component_to_metadata.R` - Add component assignments
3. `scripts/deseq2_isolate_degs.R` - Run DESeq2 analysis

All scripts are located in `/work/ugbogu/pipeline/scripts/`.
