#!/usr/bin/env python
"""
Merge STAR ReadsPerGene.out.tab count files across multiple NBL cohorts,
map Ensembl IDs to HGNC symbols, filter for protein-coding genes (optional),
apply CPM filtering, and output a unified count matrix in Parquet format.

Output:
    - final_nbl_count_matrix_hgnc.parquet : (samples x genes) int32 count matrix
    - final_nbl_sample_metadata.csv       : Sample metadata
    - final_nbl_gene_list.txt             : List of final HGNC gene symbols
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Set, Dict, List

# ==============================
# Configuration & Paths
# ==============================
ROOT_DIR = Path("/home/chu25/data/nbl")
COHORTS = ["GSE100148", "GSE189367", "GSE49711", "SRP409177", "TARGET-NBL"]

MAPPING_FILE = Path("/home/chu25/data/ensembl_to_hgnc_map.tsv")
OUTPUT_FILE = ROOT_DIR / "final_nbl_count_matrix_hgnc.parquet"
META_FILE = ROOT_DIR / "final_nbl_sample_metadata.csv"
GENE_LIST_FILE = OUTPUT_FILE.parent / "final_nbl_gene_list.txt"

# Optional: Biotype file for protein-coding gene filtering
BIOTYPE_MAP_PATH = Path("/home/chu25/data/ensembl_biotype.tsv")


# ==============================
# Load Gene Mapping (Ensembl → HGNC)
# ==============================
def load_gene_mapping() -> pd.Series:
    """
    Load Ensembl ID to HGNC symbol mapping.
    Returns a Series: Ensembl_ID (index) → HGNC_Symbol (value), with NaN dropped.
    """
    print(f"Loading gene mapping from: {MAPPING_FILE}")
    gene_map = pd.read_csv(MAPPING_FILE, sep="\t")
    gene_map = gene_map.set_index("Ensembl_ID")["HGNC_Symbol"].dropna()
    print(f"   → Loaded {len(gene_map):,} Ensembl → HGNC mappings")
    return gene_map


# ==============================
# Optional: Load Protein-Coding Filter
# ==============================
def load_protein_coding_ids() -> Optional[Set[str]]:
    """
    Load Ensembl biotype file and return set of protein-coding Ensembl IDs.
    Returns None if file doesn't exist.
    """
    if not BIOTYPE_MAP_PATH.exists():
        print("Biotype file not found. Skipping protein-coding filter.")
        return None

    print(f"Loading biotype mapping from: {BIOTYPE_MAP_PATH}")
    biotype_map = pd.read_csv(BIOTYPE_MAP_PATH, sep="\t").set_index("Ensembl_ID")["gene_biotype"]
    protein_coding_ids = set(biotype_map[biotype_map == "protein_coding"].index)
    print(f"   → Found {len(protein_coding_ids):,} protein-coding genes")
    return protein_coding_ids


# ==============================
# Auto-detect Strandedness in STAR Output
# ==============================
def detect_stranded_col(df_raw: pd.DataFrame) -> int:
    """
    Auto-detect strandedness in STAR ReadsPerGene.out.tab file.
    Analyzes columns 1, 2, 3 (0-indexed) which correspond to:
        col1: unstranded
        col2: forward stranded (+)
        col3: reverse stranded (-)

    Returns:
        2 (unstranded), 3 (forward), or 4 (reverse)
    """
    # Remove summary rows (e.g., N_unmapped, N_multimapping)
    df = df_raw[~df_raw.iloc[:, 0].astype(str).str.startswith("N_", na=False)].copy()

    if df.empty or df.shape[1] < 4:
        print("   → Not enough data/columns. Defaulting to unstranded.")
        return 2  # fallback

    # Sum counts in each stranded column
    s2 = pd.to_numeric(df.iloc[:, 1], errors="coerce").fillna(0).sum()  # unstranded
    s3 = pd.to_numeric(df.iloc[:, 2], errors="coerce").fillna(0).sum()  # forward
    s4 = pd.to_numeric(df.iloc[:, 3], errors="coerce").fillna(0).sum()  # reverse

    # Decision logic: one stranded count must dominate
    if s3 >= 1.5 * s4 and s3 >= 0.6 * s2:
        return 3  # forward stranded
    if s4 >= 1.5 * s3 and s4 >= 0.6 * s2:
        return 4  # reverse stranded
    return 2  # default: unstranded


# ==============================
# Load STAR Counts (ReadsPerGene.out.tab)
# ==============================
def load_star_counts_tsv(
    path: Path,
    sample_name: str,
    stranded_col: Optional[int] = None,
    protein_coding_ids: Optional[Set[str]] = None,
) -> pd.Series:
    """
    Load a single STAR ReadsPerGene.out.tab file.

    Args:
        path: Path to .tab file
        sample_name: Name for output column
        stranded_col: 2, 3, or 4. If None, auto-detect.
        protein_coding_ids: Optional set to filter protein-coding genes

    Returns:
        pd.Series: Ensembl_ID (index) → count (int32)
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        engine="python",
        dtype={0: "string"},
        na_filter=False,
    )

    if df.shape[1] < 4:
        raise ValueError(f"{path} has {df.shape[1]} columns; expected ≥4 for STAR output")

    # Auto-detect strandedness if not provided
    col = stranded_col if stranded_col is not None else detect_stranded_col(df)
    col_name_map = {2: "unstranded", 3: "forward/+", 4: "reverse/-"}
    print(f"   [{sample_name}] Using column {col} ({col_name_map.get(col, 'unknown')})")

    # Extract gene ID and selected count column
    sub = df.iloc[:, [0, col - 1]].copy()
    sub.columns = ["gene_id_ver", sample_name]
    sub["gene_id_ver"] = sub["gene_id_ver"].astype("string").str.strip()

    # Remove STAR summary rows
    sub = sub[~sub["gene_id_ver"].str.startswith("N_", na=False)]

    # Convert counts to int32
    sub[sample_name] = pd.to_numeric(sub[sample_name], errors="coerce").fillna(0).astype("int32")

    # Strip version from Ensembl ID (e.g., ENSG000001 → ENSG000001)
    sub["Ensembl_ID"] = sub["gene_id_ver"].str.replace(r"\.\d+$", "", regex=True)

    # Optional: filter to protein-coding genes
    if protein_coding_ids is not None:
        before = len(sub)
        sub = sub[sub["Ensembl_ID"].isin(protein_coding_ids)]
        print(f"   → Protein-coding filter: {before} → {len(sub)} genes")

    # Return as Series: Ensembl_ID → count
    return sub[["Ensembl_ID", sample_name]].set_index("Ensembl_ID")[sample_name]


# ==============================
# CPM-based Gene Filtering
# ==============================
def cpm_filter(
    counts_df: pd.DataFrame,
    min_cpm: float = 1.0,
    min_samples: int = 3,
) -> pd.DataFrame:
    """
    Filter genes by CPM (counts per million).

    Keeps genes with CPM >= min_cpm in at least min_samples samples.

    Args:
        counts_df: DataFrame (genes x samples), raw counts
        min_cpm: Minimum CPM threshold
        min_samples: Minimum number of samples passing threshold

    Returns:
        Filtered DataFrame
    """
    library_sizes = counts_df.sum(axis=0).replace(0, 1)  # avoid div-by-zero
    cpm = counts_df.div(library_sizes, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    filtered = counts_df.loc[keep]
    print(f"   → CPM ≥ {min_cpm} in ≥ {min_samples} samples: {len(counts_df)} → {len(filtered)} genes")
    return filtered


# ==============================
# Main Processing Loop
# ==============================
def main():
    # Load mappings
    gene_map = load_gene_mapping()
    protein_coding_ids = load_protein_coding_ids()

    all_cohort_matrices: List[pd.DataFrame] = []
    meta_rows: List[Dict[str, str]] = []

    # Process each cohort
    for cohort in COHORTS:
        cohort_path = ROOT_DIR / cohort / "results" / "valid_count"
        files = sorted(cohort_path.glob("*.tab"))

        if not files:
            print(f"\nWARNING: No .tab files found in {cohort_path}")
            continue

        print(f"\nProcessing cohort: {cohort} ({len(files)} samples)")
        sample_series_list: List[pd.Series] = []

        # Load each sample
        for fp in files:
            sample_name = f"{cohort}_{fp.stem}"
            try:
                counts_series = load_star_counts_tsv(
                    path=fp,
                    sample_name=sample_name,
                    stranded_col=None,
                    protein_coding_ids=protein_coding_ids,
                )
                sample_series_list.append(counts_series)

                # Record metadata
                meta_rows.append({
                    "sample": sample_name,
                    "cohort": cohort,
                    "file": str(fp),
                })
            except Exception as e:
                print(f"   ERROR loading {fp.name}: {e}")
                continue

        if not sample_series_list:
            print(f"   No valid samples loaded for {cohort}. Skipping.")
            continue

        # Merge samples in cohort
        cohort_counts = pd.concat(sample_series_list, axis=1).fillna(0).astype("int32")
        print(f"   Initial matrix: {cohort_counts.shape[1]} samples × {cohort_counts.shape[0]} genes")

        # CPM filtering (adaptive min_samples: at least 5% of samples or 3)
        lib_sizes = cohort_counts.sum(axis=0).replace(0, 1)
        min_samples = max(3, int(0.05 * cohort_counts.shape[1]))
        cohort_counts = cpm_filter(cohort_counts, min_cpm=1.0, min_samples=min_samples)

        # Map Ensembl → HGNC and sum duplicates
        hgnc_symbols = cohort_counts.index.to_series().map(gene_map)
        mapped = cohort_counts.loc[hgnc_symbols.notna()].copy()
        mapped.index = hgnc_symbols[hgnc_symbols.notna()].values
        mapped = mapped.groupby(level=0).sum().astype("int32")

        print(f"   After HGNC mapping & dedup: {mapped.shape}")
        all_cohort_matrices.append(mapped)

    # ==============================
    # Final Merge Across Cohorts
    # ==============================
    if not all_cohort_matrices:
        raise RuntimeError("No data processed from any cohort.")

    print("\nMerging all cohorts...")
    final = pd.concat(all_cohort_matrices, axis=1, join="outer").fillna(0).astype("int32")
    final = final.loc[~final.index.duplicated(keep='first')].sort_index()
    final_t = final.T  # Samples x Genes
    final_t = final_t.sort_index(axis=1)  # Sort genes alphabetically

    print(f"\nFinal matrix (Samples × Genes): {final_t.shape}")
    print(f"   Sample preview: {final_t.index[:5].tolist()}")
    print(f"   Gene preview : {final_t.columns[:5].tolist()} ...")

    # ==============================
    # Save Outputs
    # ==============================
    print(f"\nSaving count matrix → {OUTPUT_FILE}")
    final_t.to_parquet(OUTPUT_FILE, index=True, engine="pyarrow", compression="zstd")

    print(f"Saving metadata → {META_FILE}")
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(META_FILE, index=False)

    print(f"Saving gene list → {GENE_LIST_FILE}")
    final_t.columns.to_series().to_csv(GENE_LIST_FILE, index=False, header=False)
    print(f"   → {len(final_t.columns):,} genes saved")

    print("\nProcessing complete!")


# ==============================
# Entry Point
# ==============================
if __name__ == "__main__":
    main()
