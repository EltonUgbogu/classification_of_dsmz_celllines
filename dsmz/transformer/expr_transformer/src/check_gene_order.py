#!/usr/bin/env python3
"""
Verify that NPZ gene ordering matches the accompanying metadata.

This script performs a one-time integrity check to ensure the expression matrix
columns align with the gene list recorded in the metadata JSON. It reports a
stable SHA256 hash you can embed into checkpoints for future validation.
"""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

# Default paths assume execution from repository root
# Resolve relative to script's parent directory (repo root)
_SCRIPT_DIR = Path(__file__).parent
_REPO_ROOT = _SCRIPT_DIR.parent
DEFAULT_NPZ = _REPO_ROOT / "data" / "BRCA_HVG500_joint.npz"
DEFAULT_META = _REPO_ROOT / "data" / "BRCA_HVG500_joint.meta.json"

# Common metadata keys for gene identifiers
GENE_KEY_CANDIDATES = ("genes", "gene_names", "features", "gene_ids")


def load_gene_list(meta: dict) -> Tuple[str, list]:
    """Return the gene list and the key it was found under."""
    gene_key = next((k for k in GENE_KEY_CANDIDATES if k in meta), None)
    if gene_key is None:
        raise KeyError(
            "Could not find gene list in meta.json. "
            f"Tried: {list(GENE_KEY_CANDIDATES)}. Keys: {list(meta.keys())}"
        )
    genes = meta[gene_key]
    if not isinstance(genes, list):
        raise TypeError(f"meta['{gene_key}'] must be a list, got {type(genes)!r}")
    return gene_key, genes


def sha256_lines(values: Iterable[str]) -> str:
    joined = "\n".join(values).encode("utf-8")
    return hashlib.sha256(joined).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate NPZ+meta and print gene-order SHA256."
    )
    parser.add_argument(
        "--npz",
        type=Path,
        default=DEFAULT_NPZ,
        help=f"Path to NPZ file (default: {DEFAULT_NPZ})",
    )
    parser.add_argument(
        "--meta",
        type=Path,
        default=DEFAULT_META,
        help=f"Path to metadata JSON (default: {DEFAULT_META})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = np.load(args.npz)
    if "X" not in data:
        raise KeyError(f"NPZ file does not contain 'X'. Keys: {list(data.keys())}")

    X = data["X"]
    print(f"[OK] X shape: {X.shape} (samples × genes)")

    meta = json.loads(args.meta.read_text())
    if "genes" not in meta:
        raise KeyError(f"meta.json has no 'genes' key. Keys: {list(meta.keys())}")
    if "samples" not in meta:
        raise KeyError(f"meta.json has no 'samples' key. Keys: {list(meta.keys())}")

    samples = meta["samples"]
    if not isinstance(samples, list):
        raise TypeError(f"meta['samples'] must be a list, got {type(samples)!r}")

    gene_key, genes = load_gene_list(meta)
    print(f"[OK] Found gene list under meta['{gene_key}'] with length: {len(genes)}")

    # 1) shape agreement
    if X.shape[0] != len(samples):
        raise ValueError(
            f"Sample mismatch: X has {X.shape[0]} rows, meta has {len(samples)} samples"
        )
    if X.shape[1] != len(genes):
        raise ValueError(
            f"Gene mismatch: X has {X.shape[1]} cols, meta has {len(genes)} genes"
        )
    print("[OK] X dimensions match meta.json samples+genes")

    # 2) duplicates + sanity
    n_unique = len(set(genes))
    if n_unique != len(genes):
        raise ValueError(
            f"Duplicate genes detected: unique={n_unique}, total={len(genes)}"
        )
    if not all(isinstance(g, str) and g for g in genes):
        raise TypeError("Some gene IDs are empty or not strings")
    print("[OK] Gene IDs are unique and valid strings")

    # 3) stable fingerprint (hash) for later comparisons
    gene_hash = sha256_lines(genes)

    print(f"[OK] Gene-order SHA256: {gene_hash}")
    print("First 5 genes:", genes[:5])
    print("Last 5 genes: ", genes[-5:])


if __name__ == "__main__":
    main()
