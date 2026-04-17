#!/usr/bin/env python3
# Build expression matrix for TARGET NBL RNA-seq samples

from functools import reduce
from pathlib import Path
import pandas as pd

BASE = Path("/work/ugbogu/pipeline/data/nbl/target_nbl")
SAMPLE_LIST = BASE / "target_nbl_rnaseq_sample_ids.tsv"
COUNT_DIR = BASE / "count_data"
OUTFILE = BASE / "TARGET_expression_matrix.tsv"

mapping = pd.read_csv(SAMPLE_LIST, sep="\t")

if "sample_id" not in mapping.columns:
    raise ValueError("Expected column 'sample_id' in target_nbl_rnaseq_sample_ids.tsv")

dfs = []

for _, row in mapping.iterrows():
    sample_id = row["sample_id"]
    sample_dir = COUNT_DIR / sample_id

    if not sample_dir.exists():
        print(f"Warning: directory not found, skipping: {sample_dir}")
        continue

    tsv_files = list(sample_dir.glob("*.tsv"))
    if not tsv_files:
        print(f"Warning: no TSV file found in {sample_dir}, skipping")
        continue
    if len(tsv_files) > 1:
        print(f"Warning: multiple TSV files found in {sample_dir}, using first one: {tsv_files[0]}")

    file_path = tsv_files[0]

    try:
        df = pd.read_csv(file_path, sep="\t", comment="#")
    except Exception as e:
        print(f"Warning: failed to read {file_path}: {e}")
        continue

    df.columns = df.columns.str.strip()

    required_cols = ["gene_id", "unstranded"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Warning: missing columns {missing} in {file_path}, skipping")
        continue

    # Remove STAR summary rows like N_unmapped, N_multimapping, etc.
    df = df[~df["gene_id"].astype(str).str.startswith("N_")].copy()

    # Keep only gene_id and raw counts
    df = df[["gene_id", "unstranded"]].copy()
    df.columns = ["gene_id", sample_id]

    dfs.append(df)

if not dfs:
    raise ValueError("No valid count files were loaded.")

expr = reduce(lambda left, right: pd.merge(left, right, on="gene_id", how="inner"), dfs)

expr.to_csv(OUTFILE, sep="\t", index=False)

print(f"Saved expression matrix: {OUTFILE}")
print(f"Genes: {expr.shape[0]}")
print(f"Samples: {expr.shape[1] - 1}")