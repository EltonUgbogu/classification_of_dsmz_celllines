# DESeq2 Analysis Plan: Anchor + Isolate DEGs

## Summary

This document defines the **final anchor and isolate selections** for BRCA, NBL, and RBL based on graph topology analysis (degree and betweenness centrality).

---

## 🟣 BRCA (Breast Cancer)

### Anchors (5 total)

| Anchor | Component | Rationale |
|--------|-----------|-----------|
| **BT_474** | 0 | Highest degree (6) - hub node |
| **EFM_192C** | 0 | Highest betweenness (0.127) - connector |
| **HDQ_P1** | 1 | Highest degree (5) + betweenness (0.011) |
| **COLO_824** | 3 | Only meaningful anchor (degree = 2) |
| **ETCC_006** | 5 | Representative of triangle micro-module |

**File:** `results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv`

### Isolates (3 total)

- **CAL_51** (component 2)
- **DU_4475** (component 4)
- **HCC_1937** (component 6)

**File:** `results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/isolate_list.txt`

### DESeq2 Contrasts

**Isolate DEGs (3 contrasts):**
- `CAL_51 vs (all other cell lines)`
- `DU_4475 vs (all other cell lines)`
- `HCC_1937 vs (all other cell lines)`

**Anchor DEGs (5 contrasts):**
- `BT_474 vs (cell lines NOT in component 0)`
- `EFM_192C vs (cell lines NOT in component 0)`
- `HDQ_P1 vs (cell lines NOT in component 1)`
- `COLO_824 vs (cell lines NOT in component 3)`
- `ETCC_006 vs (cell lines NOT in component 5)`

**Total:** 8 contrasts

---

## 🔵 NBL (Neuroblastoma)

### Anchors (4 total)

| Anchor | Component | Rationale |
|--------|-----------|-----------|
| **IMR_32** | 1 | Highest degree (4) |
| **CHP_126** | 1 | Highest betweenness (0.213) - bridge |
| **CCLF_PEDS_0051_T** | 0 | Representative of 2-node pair |
| **LAN_5** | 2 | Representative of 2-node pair |

**File:** `results/unsupervised/nbl/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv`

### Isolates (3 total)

- **LS** (component 3)
- **SH_SY5Y** (component 4)
- **SK_N_BE_2** (component 5)

**File:** `results/unsupervised/nbl/tumour_neighbourhoods/final_consensus_all/isolate_list.txt`

### DESeq2 Contrasts

**Isolate DEGs (3 contrasts):**
- `LS vs (all other cell lines)`
- `SH_SY5Y vs (all other cell lines)`
- `SK_N_BE_2 vs (all other cell lines)`

**Anchor DEGs (4 contrasts):**
- `IMR_32 vs (cell lines NOT in component 1)`
- `CHP_126 vs (cell lines NOT in component 1)`
- `CCLF_PEDS_0051_T vs (cell lines NOT in component 0)`
- `LAN_5 vs (cell lines NOT in component 2)`

**Total:** 7 contrasts

---

## 🟠 RBL (Retinoblastoma)

### Anchors (2 total)

| Anchor | Component | Rationale |
|--------|-----------|-----------|
| **NG-30919_RBL_18_lib626623_10098_1** | 0 | Highest degree (3) + betweenness (0.089) |
| **NG-30919_RBL_15_lib626622_10098_2** | 2 | Representative of 2-node replicate pair |

**File:** `results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv`

### Isolates (4 total)

- **NG-30919_RBL_14_lib628470_10098_1** (component 1)
- **NG-30919_RBL_7_lib628468_10098_1** (component 3)
- **NG-30919_WERI_RB1_lib628472_10098_1** (component 4)
- **NG-30919_Y_79_lib626626_10098_1** (component 5)

**File:** `results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/isolate_list.txt`

### DESeq2 Contrasts

**Isolate DEGs (4 contrasts):**
- `NG-30919_RBL_14_lib628470_10098_1 vs (all other cell lines)`
- `NG-30919_RBL_7_lib628468_10098_1 vs (all other cell lines)`
- `NG-30919_WERI_RB1_lib628472_10098_1 vs (all other cell lines)`
- `NG-30919_Y_79_lib626626_10098_1 vs (all other cell lines)`

**Anchor DEGs (2 contrasts):**
- `NG-30919_RBL_18_lib626623_10098_1 vs (cell lines NOT in component 0)`
- `NG-30919_RBL_15_lib626622_10098_2 vs (cell lines NOT in component 2)`

**Total:** 6 contrasts

---

## DESeq2 Script Usage

### Prerequisites

1. **Raw count matrix** (integer counts, not VST/TPM/logcounts)
2. **Metadata TSV** with columns:
   - `sample_id`: matches count matrix column names
   - `cell_line`: cell line identifier
   - `component`: component assignment (for anchor DEGs)

### Example Command (BRCA)

```bash
Rscript scripts/deseq2_isolate_degs.R \
  --counts data/brca/raw_counts.tsv \
  --meta data/brca/metadata_with_components.tsv \
  --sample_id_col sample_id \
  --cell_line_col cell_line \
  --component_col component \
  --isolate_list CAL_51,DU_4475,HCC_1937 \
  --anchor_list BT_474,EFM_192C,HDQ_P1,COLO_824,ETCC_006 \
  --anchor_components results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/anchor_components.tsv \
  --outdir results/degs/brca \
  --fdr_isolate 0.01 \
  --lfc_isolate 1.5 \
  --topN_isolate 50 \
  --fdr_anchor 0.05 \
  --lfc_anchor 1.0 \
  --topN_anchor 200 \
  --recurrence_k 2
```

### Output Structure

```
results/degs/{profile}/
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

## Biological Interpretation

### Component-level programmes (anchors)
- Extract **shared, stable transcriptional programmes** per component
- Robust to within-component variation
- Directly transferable to tumours

### Singleton states (isolates)
- Capture **rare / extreme / plastic states**
- No topology forcing
- Maintains biological diversity

### Pan-cancer feature set
- Recurrence-filtered unique genes (appearing in ≥k contrasts)
- Structured pool of DEG sets, not a flat union
- Enables multi-class pan-cancer alignment preserving topology

---

## Validation Checklist

- ✅ Graph topology respected
- ✅ Rare states included
- ✅ No forced reassignment
- ✅ No DEG explosion (controlled via thresholds)
- ✅ Pan-cancer ready
- ✅ Methods-section defensible
