#!/usr/bin/env python3
"""Create tiny placeholder VST RDS files for Snakemake dry-runs.

Reads config/config.yaml, inspects the selected profiles, and creates
minimal valid placeholder inputs at the configured external source paths
if they do not already exist.  Does not overwrite real files.

Usage:
    python scripts/create_dryrun_placeholders.py \
        --config config/config.yaml \
        --profiles brca,nbl,rbl,heme
"""
import argparse
import os
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML is required to run create_dryrun_placeholders.py")

try:
    import pandas as pd
    import pyreadr
    _HAS_PYREADR = True
except ImportError:
    _HAS_PYREADR = False

# Rscript fallback for environments without pyreadr
import subprocess
import tempfile


# ---------------------------------------------------------------------------
# Placeholder writers
# ---------------------------------------------------------------------------

def write_placeholder_rds(path: Path, profile: str) -> None:
    """Write a tiny valid RDS file containing a small data frame.
    Uses pyreadr if available, falls back to Rscript.
    """
    if _HAS_PYREADR:
        genes = ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"]
        samples = [f"{profile}_tumour_1", f"{profile}_tumour_2",
                   f"{profile}_cellline_1", f"{profile}_cellline_2"]
        # Build a numeric DataFrame that mimics a VST expression matrix
        import random
        data = {s: [random.uniform(0.0, 10.0) for _ in genes] for s in samples}
        df = pd.DataFrame(data, index=genes)
        df.index.name = "gene_id"
        pyreadr.write_rds(str(path), df)
    else:
        # Fallback: use Rscript
        genes = ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"]
        samples = [f"{profile}_tumour_1", f"{profile}_tumour_2",
                   f"{profile}_cellline_1", f"{profile}_cellline_2"]
        genes_vec = ", ".join(f'"{g}"' for g in genes)
        samples_vec = ", ".join(f'"{s}"' for s in samples)
        safe_path = str(path).replace("\\", "\\\\").replace('"', '\\"')
        r_code = (
            "genes <- c(" + genes_vec + ")\n"
            "samples <- c(" + samples_vec + ")\n"
            "mat <- matrix(runif(length(genes) * length(samples)),\n"
            "              nrow=length(genes), ncol=length(samples))\n"
            "rownames(mat) <- genes\n"
            "colnames(mat) <- samples\n"
            'saveRDS(mat, file="' + safe_path + '")\n'
        )
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as tmp:
            tmp.write(r_code)
            tmp_path = tmp.name
        try:
            subprocess.run(
                ["Rscript", tmp_path], check=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
        finally:
            os.unlink(tmp_path)


def write_placeholder_tsv(path: Path, profile: str) -> None:
    """Write a tiny valid metadata TSV/CSV."""
    lineage = profile.upper()
    path.write_text(
        "sample_id,lineage,type,cell_line,component\n"
        f"{profile}_tumour_1,{lineage},tumour,,\n"
        f"{profile}_tumour_2,{lineage},tumour,,\n"
        f"{profile}_cellline_1,{lineage},cell_line,CL1,comp_1\n"
        f"{profile}_cellline_2,{lineage},cell_line,CL2,comp_2\n"
    )


def write_placeholder_deseq_table(path: Path, profile: str) -> None:
    """Write a tiny valid DESeq2 results table."""
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n"
        "ENSG00000000001\t10\t1.0\t0.001\t0.01\n"
        "ENSG00000000002\t15\t-1.2\t0.002\t0.02\n"
    )


# ---------------------------------------------------------------------------
# Helper: ensure a placeholder at a given path
# ---------------------------------------------------------------------------

def ensure_rds_placeholder(path_str, profile, manifest_fh) -> str:
    """Create a placeholder RDS if the file does not exist. Returns action taken."""
    if not path_str:
        return "empty_path"
    path = Path(path_str)
    if path.exists():
        manifest_fh.write(f"{path}\t{profile}\texists\n")
        return "exists"
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        write_placeholder_rds(path, profile)
        manifest_fh.write(f"{path}\t{profile}\tcreated\n")
        print(f"[PLACEHOLDER] Created {path}")
        return "created"
    except subprocess.CalledProcessError as err:
        manifest_fh.write(f"{path}\t{profile}\tfailure:{err.returncode}\n")
        sys.stderr.write(err.stderr.decode() if err.stderr else "")
        return f"failure:{err.returncode}"


def ensure_meta_placeholder(path_str, profile, manifest_fh) -> str:
    """Create a placeholder metadata file if it does not exist."""
    if not path_str:
        return "empty_path"
    path = Path(path_str)
    if path.exists():
        manifest_fh.write(f"{path}\t{profile}\texists\n")
        return "exists"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_placeholder_tsv(path, profile)
    manifest_fh.write(f"{path}\t{profile}\tcreated\n")
    print(f"[PLACEHOLDER] Created {path}")
    return "created"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config.yaml")
    parser.add_argument("--profiles", default="brca,nbl,rbl,heme",
                        help="Comma-separated profile list")
    parser.add_argument("--manifest",
                        default="results/dryrun_placeholders/manifest.tsv",
                        help="Path for the manifest TSV")
    args = parser.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    profiles_cfg = cfg.get("profiles", {})
    selected = [p.strip() for p in args.profiles.split(",") if p.strip()]
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    # Write fresh manifest header (overwrite previous runs).
    manifest.write_text("path\tprofile\taction\treason\n")

    had_failure = False

    with manifest.open("a") as mf:
        for profile in selected:
            profile_cfg = profiles_cfg.get(profile)
            if not profile_cfg:
                mf.write(f"\t{profile}\tskipped\tunknown_profile\n")
                print(f"[SKIP] Unknown profile: {profile}")
                continue

            paths_cfg = profile_cfg.get("paths", {})

            # RDS input paths: vst_joint_rds, cell_vst_rds, tumour_vst_rds, dsmz_counts_rds
            rds_keys = ["vst_joint_rds", "cell_vst_rds", "tumour_vst_rds", "dsmz_counts_rds"]
            seen_rds = set()
            for key in rds_keys:
                rds_path = paths_cfg.get(key, "")
                if rds_path and rds_path not in seen_rds:
                    seen_rds.add(rds_path)
                    result = ensure_rds_placeholder(rds_path, profile, mf)
                    if result.startswith("failure"):
                        had_failure = True

            # Metadata paths: meta_tsv, dsmz_meta_csv
            meta_keys = ["meta_tsv", "dsmz_meta_csv"]
            seen_meta = set()
            for key in meta_keys:
                meta_path = paths_cfg.get(key, "")
                if meta_path and meta_path not in seen_meta:
                    seen_meta.add(meta_path)
                    ensure_meta_placeholder(meta_path, profile, mf)

            # DESeq2 placeholder tables directory
            unsup_root = paths_cfg.get("unsup_root") or f"results/unsupervised/{profile}"
            tables_dir = Path(unsup_root) / "deseq2_markers" / "tables"
            if not tables_dir.exists():
                tables_dir.mkdir(parents=True, exist_ok=True)
            placeholder_table = tables_dir / f"anchor_{profile}_placeholder_vs_outside_component_0001.tsv"
            if not placeholder_table.exists():
                write_placeholder_deseq_table(placeholder_table, profile)
                mf.write(f"{placeholder_table}\t{profile}\tcreated\tdeseq2_table\n")
                print(f"[PLACEHOLDER] Created {placeholder_table}")
            else:
                mf.write(f"{placeholder_table}\t{profile}\texists\tdeseq2_table\n")

    print(f"\n[MANIFEST] Written to {manifest}")

    if had_failure:
        print("[WARNING] Some placeholders failed to create. Check the manifest.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
