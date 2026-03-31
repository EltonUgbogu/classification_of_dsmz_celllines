#!/usr/bin/env python3
"""Create tiny placeholder VST RDS files for Snakemake dry-runs."""
import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML is required to run create_dryrun_placeholders.py")


def write_placeholder_rds(path: Path, profile: str) -> None:
    genes = ["ENSG00000000001", "ENSG00000000002"]
    samples = [f"{profile}_tumour", f"{profile}_cellline"]
    genes_vec = ", ".join(f'"{g}"' for g in genes)
    samples_vec = ", ".join(f'"{s}"' for s in samples)
    safe_path = str(path).replace("\\", "\\\\").replace('"', '\\"')
    r_code = (
        "genes <- c(" + genes_vec + ")\n"
        "samples <- c(" + samples_vec + ")\n"
        "mat <- matrix(seq_along(genes) * 0.01, nrow=length(genes), ncol=length(samples))\n"
        "rownames(mat) <- genes\n"
        "colnames(mat) <- samples\n"
        "saveRDS(mat, file=\"" + safe_path + "\")\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as tmp:
        tmp.write(r_code)
        tmp_path = tmp.name
    try:
        subprocess.run(["Rscript", tmp_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    finally:
        os.unlink(tmp_path)


def write_placeholder_tsv(path: Path, profile: str) -> None:
    lineage = profile.upper()
    path.write_text(
        "sample_id,lineage,type\n"
        f"{profile}_tumour,{lineage},tumour\n"
        f"{profile}_cellline,{lineage},cell_line\n"
    )


def write_placeholder_deseq_table(path: Path, profile: str) -> None:
    path.write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n"
        f"ENSG00000000001\t10\t1.0\t0.001\t0.01\n"
        f"ENSG00000000002\t15\t-1.2\t0.002\t0.02\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--profiles", default="brca,nbl,rbl,heme",
                        help="Comma-separated profile list")
    parser.add_argument("--manifest", default="results/dryrun_placeholders/manifest.tsv")
    args = parser.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    profiles_cfg = cfg.get("profiles", {})
    selected = [p.strip() for p in args.profiles.split(",") if p.strip()]
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    if not manifest.exists():
        manifest.write_text("path\tprofile\taction\n")

    with manifest.open("a") as mf:
        def ensure_placeholder(path_str: str, profile: str) -> None:
            if not path_str:
                return
            path = Path(path_str)
            if path.exists():
                mf.write(f"{path}\t{profile}\texists\n")
                return
            path.parent.mkdir(parents=True, exist_ok=True)
            try:
                write_placeholder_rds(path, profile)
                mf.write(f"{path}\t{profile}\tcreated\n")
                print(f"[PLACEHOLDER] Created {path}")
            except subprocess.CalledProcessError as err:
                mf.write(f"{path}\t{profile}\tfailure:{err.returncode}\n")
                sys.stderr.write(err.stderr.decode())
                sys.exit(err.returncode)

        for profile in selected:
            profile_cfg = profiles_cfg.get(profile)
            if not profile_cfg:
                mf.write(f"\t{profile}\tskipped_unknown_profile\n")
                continue
            paths_cfg = profile_cfg.get("paths", {})
            ensure_placeholder(paths_cfg.get("vst_joint_rds"), profile)
            ensure_placeholder(paths_cfg.get("dsmz_counts_rds"), profile)
            meta_path = paths_cfg.get("meta_tsv")
            if meta_path:
                meta = Path(meta_path)
                if meta.exists():
                    mf.write(f"{meta}\t{profile}\texists\n")
                else:
                    meta.parent.mkdir(parents=True, exist_ok=True)
                    write_placeholder_tsv(meta, profile)
                    mf.write(f"{meta}\t{profile}\tcreated\n")
                    print(f"[PLACEHOLDER] Created {meta}")

            unsup_root = paths_cfg.get("unsup_root") or f"results/unsupervised/{profile}"
            tables_dir = Path(unsup_root) / "deseq2_markers" / "tables"
            if not tables_dir.exists():
                tables_dir.mkdir(parents=True, exist_ok=True)
            placeholder_table = tables_dir / f"anchor_{profile}_placeholder_vs_outside_component_0001.tsv"
            if not placeholder_table.exists():
                write_placeholder_deseq_table(placeholder_table, profile)
                mf.write(f"{placeholder_table}\t{profile}\tcreated\n")
                print(f"[PLACEHOLDER] Created {placeholder_table}")


if __name__ == "__main__":
    main()
