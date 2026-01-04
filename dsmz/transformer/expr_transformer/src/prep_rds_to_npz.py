#!/usr/bin/env python3
"""
# prep_rds_to_npz.py
==================

This script converts an R Data Serialization (RDS) file containing a gene
expression matrix into a Python-friendly NumPy archive (NPZ) format.

**Why this conversion?**

- RDS files are R's native binary format, efficient for R workflows
- NPZ files are NumPy's compressed archive format, ideal for Python/PyTorch


**What it does:**

1. Reads an RDS file containing a samples × genes expression matrix
2. Applies z-score normalization (standardization) per gene
3. Saves the normalized matrix as a compressed NPZ file
4. Saves sample/gene names as a JSON metadata file

**Z-score normalization:**

For each gene: z = (x - μ) / σ

This centers each gene at mean=0 with standard deviation=1, which:
- Puts all genes on the same scale (important for neural networks)
- Prevents genes with high expression from dominating
- Allows the model to focus on relative patterns, not absolute values

**Input:**
    - BRCA_TCGA-DSMZ_HVG500_samples_x_genes.rds (from R preprocessing pipeline)

**Outputs:**
    - BRCA_HVG500_joint.npz: Compressed array with X (data), mu, sd (for reverse transform)
    - BRCA_HVG500_joint.meta.json: Sample IDs and gene names
"""

import argparse
import hashlib
import json
import os
import sys


META_GENE_KEY = "genes"


def sha256_lines(values):
    """Compute SHA256 hash of a list of strings (one per line)."""
    joined = "\n".join(values).encode("utf-8")
    return hashlib.sha256(joined).hexdigest()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert expression matrix RDS to NPZ and metadata JSON",
    )
    parser.add_argument(
        "--rds",
        required=True,
        help="Input RDS file containing samples × genes matrix (pandas DataFrame after read)",
    )
    parser.add_argument(
        "--npz",
        required=True,
        help="Path to write compressed NPZ archive",
    )
    parser.add_argument(
        "--meta",
        required=True,
        help="Path to write metadata JSON",
    )
    parser.add_argument(
        "--cancer",
        required=True,
        help="Cancer label for this dataset (e.g., brca, nbl, rbl, pan)",
    )
    parser.add_argument(
        "--no-zscore",
        action="store_true",
        help="If set, do NOT z-score normalize; store raw matrix as X.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Heavy imports AFTER args are parsed (so --help works without deps)
    try:
        import numpy as np
        import pandas as pd
        import pyreadr  # Library for reading R data files in Python
    except ModuleNotFoundError as e:
        print(f"[ERROR] Missing dependency: {e.name}. Activate env or install it.", file=sys.stderr)
        return 2

    # =============================================================================
    # LOAD RDS FILE
    # =============================================================================

    # pyreadr.read_r() returns an OrderedDict where:
    #   - Keys are object names from R (or None for unnamed objects)
    #   - Values are the actual data (converted to pandas DataFrames or numpy arrays)
    #
    # An RDS file typically contains a single R object (unlike RData which can have many)
    print(f"Reading RDS file: {args.rds}")
    obj = pyreadr.read_r(args.rds)

    # Extract the first (and usually only) object from the RDS file
    # next(iter(...)) is a Pythonic way to get the first value from a dict
    mat = next(iter(obj.values()))

    # Validate that we got a DataFrame (expected for expression matrices)
    # RDS files can contain various R objects; we specifically need a matrix/data.frame
    if not isinstance(mat, pd.DataFrame):
        raise TypeError(f"Expected DataFrame from RDS, got: {type(mat)}")

    print(f"Loaded matrix: {mat.shape[0]} samples × {mat.shape[1]} genes")


    # =============================================================================
    # EXTRACT DATA AND METADATA
    # =============================================================================

    # Convert to NumPy array with 32-bit floats
    # - float32 uses half the memory of float64 (default)
    # - Sufficient precision for neural networks
    # - Faster computation, especially on GPUs
    X = mat.to_numpy(dtype=np.float32)

    # Extract sample IDs from DataFrame index (row names in R)
    # Examples: "TCGA-A1-A0SK-01A", "ACC-201", etc.
    # map(str, ...) ensures all IDs are strings (handles numeric indices)
    samples = list(map(str, mat.index))

    # Extract gene names from DataFrame columns
    # Examples: "ENSG00000141510", "TP53", etc.
    genes = list(map(str, mat.columns))

    print(f"Samples: {samples[:3]} ... ({len(samples)} total)")
    print(f"Genes: {genes[:3]} ... ({len(genes)} total)")

    # =============================================================================
    # Z-SCORE NORMALIZATION (STANDARDIZATION)
    # =============================================================================

    # Z-score formula: z = (x - μ) / σ
    #
    # **Why z-score normalize?**
    #
    # 1. Scale invariance: Genes with high baseline expression (e.g., housekeeping
    #    genes) won't dominate the model's attention over lowly-expressed genes.
    #
    # 2. Neural network training: Most neural networks work best when inputs are
    #    centered around 0 with unit variance. This helps with:
    #    - Gradient flow (avoids vanishing/exploding gradients)
    #    - Weight initialization assumptions
    #    - Faster convergence
    #
    # 3. Biological interpretation: We care about relative changes in expression,
    #    not absolute counts. Z-scores capture "how many standard deviations
    #    from normal" each measurement is.
    #
    # **Axis=0 means "per gene" (across all samples):**
    #
    #    Sample1:  gene1  gene2  gene3
    #    Sample2:  gene1  gene2  gene3
    #    Sample3:  gene1  gene2  gene3
    #              ↓      ↓      ↓
    #           normalize each column independently

    # Compute mean and std for each gene (axis=0 = along samples)
    # keepdims=True maintains shape (1, n_genes) for broadcasting
    mu = X.mean(axis=0, keepdims=True).astype(np.float32)
    sd = X.std(axis=0, keepdims=True).astype(np.float32)

    # Handle edge case: genes with zero variance (constant expression)
    # - Division by zero would produce NaN/Inf
    # - Setting sd=1 means z-score becomes just (x - mu), which equals 0 for constant values
    # - These genes carry no information anyway and could be filtered upstream
    sd[sd == 0] = 1.0

    if args.no_zscore:
        Xz = X
        print("\n[INFO] --no-zscore set: saving raw X (VST-scale) to NPZ")
    else:
        # Apply z-score transformation
        # Broadcasting: (n_samples, n_genes) - (1, n_genes) / (1, n_genes)
        Xz = (X - mu) / sd
        print(f"\nZ-score normalization applied:")
        print(f"  Before: mean={X.mean():.3f}, std={X.std():.3f}")
        print(f"  After:  mean={Xz.mean():.6f}, std={Xz.std():.3f}")

    # =============================================================================
    # SAVE OUTPUTS
    # =============================================================================

    # Create output directory if it doesn't exist
    # exist_ok=True prevents error if directory already exists
    os.makedirs(os.path.dirname(args.npz), exist_ok=True)

    # -----------------------------------------------------------------------------
    # Save compressed NPZ archive
    # -----------------------------------------------------------------------------
    # NPZ format stores multiple named arrays in a single file:
    #   - X: The z-scored expression matrix (n_samples, n_genes)
    #   - mu: Per-gene means (1, n_genes) - for reverse transformation
    #   - sd: Per-gene stds (1, n_genes) - for reverse transformation
    #
    # To reverse z-score later: X_original = X * sd + mu
    #
    # savez_compressed uses ZIP compression (typically 2-5x smaller than raw binary)
    # Always save mu and sd even if --no-zscore (useful for reference/debugging)
    np.savez_compressed(
        args.npz,
        X=Xz,
        mu=mu,
        sd=sd
    )

    # -----------------------------------------------------------------------------
    # Save metadata as JSON
    # -----------------------------------------------------------------------------
    # JSON format includes:
    # - samples: List of sample IDs
    # - genes: List of gene names (in order matching NPZ columns)
    # - gene_key: Which key contains the gene list (for consistency)
    # - gene_hash: SHA256 hash of gene list (for checkpoint validation)
    # - cancer: Cancer type label
    # - n_samples, n_genes: Counts for quick reference
    # - source_rds: Original RDS file path
    # - normalized: Whether z-score normalization was applied
    #
    # Why separate from NPZ?
    # - JSON is human-readable (can inspect with any text editor)
    # - Easy to load in any language (R, JavaScript, etc.)
    # - Sample/gene names can have complex characters that NPZ handles poorly
    gene_hash = sha256_lines(genes)

    meta_obj = {
        "samples": samples,
        META_GENE_KEY: genes,
        "gene_key": META_GENE_KEY,
        "gene_hash": gene_hash,
        "cancer": str(args.cancer),
        "n_samples": int(len(samples)),
        "n_genes": int(len(genes)),
        "source_rds": os.path.abspath(args.rds),
        "normalized": (not args.no_zscore),
    }

    os.makedirs(os.path.dirname(args.meta), exist_ok=True)
    with open(args.meta, "w") as f:
        json.dump(meta_obj, f, indent=2)

    # =============================================================================
    # SUMMARY
    # =============================================================================

    print(f"\n{'='*60}")
    print("CONVERSION COMPLETE")
    print(f"{'='*60}")
    print(f"NPZ file:  {args.npz}")
    print(f"Meta file: {args.meta}")
    print(f"Shape:     {Xz.shape[0]} samples × {Xz.shape[1]} genes")
    print(f"\nTo load in Python:")
    print(f"    data = np.load('{os.path.basename(args.npz)}')")
    print(f"    X = data['X']  # z-scored expression matrix")
    print(f"    # To reverse: X_original = X * data['sd'] + data['mu']")
    print(f"{'='*60}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
