#!/usr/bin/env python3
"""
Compare GSE189367 sample-export SRR IDs against SRRs present in FASTQ files.

Defaults:
- sample list: /work/ugbogu/pipeline/data/nbl/preprocessing_results/sample_id_exports/GSE189367_sample_ids.tsv
- fastq dir:   /work/ugbogu/pipeline/data/nbl/GSE189367/fastq
- output dir:  /work/ugbogu/pipeline/data/nbl

Outputs:
- gse189367_ids_in_sample_and_fastq.txt
- gse189367_ids_only_in_sample_export.txt
- gse189367_ids_only_in_fastq.txt
- gse189367_sample_vs_fastq_report.txt
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_SAMPLE_LIST = Path(
    "/work/ugbogu/pipeline/data/nbl/preprocessing_results/sample_id_exports/GSE189367_sample_ids.tsv"
)
DEFAULT_FASTQ_DIR = Path(
    "/work/ugbogu/pipeline/data/nbl/GSE189367/fastq"
)
DEFAULT_OUTPUT_DIR = Path("/work/ugbogu/pipeline/data/nbl")


def extract_srr_from_sample_export(path: Path) -> set[str]:
    """
    Read sample export lines such as:
        sample_id
        GSE189367_SRR17011027
    or:
        SRP409177_SRR22373268
    or:
        SRR17011027

    Return:
        {"SRR17011027", ...}
    """
    ids: set[str] = set()

    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        value = raw.strip()
        if not value or value == "sample_id":
            continue

        if "_SRR" in value:
            srr = value.split("_")[-1]
            if srr.startswith("SRR"):
                ids.add(srr)
        elif value.startswith("SRR"):
            ids.add(value)

    return ids


def load_srr_ids_from_fastq_dir(path: Path) -> set[str]:
    """
    Extract SRR IDs from FASTQ filenames in a directory.

    Only count an SRR as present if both paired files exist:
      SRRxxxxxxx_1.fastq.gz
      SRRxxxxxxx_2.fastq.gz
    """
    read1: set[str] = set()
    read2: set[str] = set()

    for fq in path.glob("*.fastq.gz"):
        name = fq.name
        if name.endswith("_1.fastq.gz"):
            read1.add(name[: -len("_1.fastq.gz")])
        elif name.endswith("_2.fastq.gz"):
            read2.add(name[: -len("_2.fastq.gz")])

    return read1 & read2


def write_lines(path: Path, values: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for item in values:
            handle.write(item + "\n")


def build_report(
    sample_path: Path,
    fastq_dir: Path,
    in_both: list[str],
    only_sample: list[str],
    only_fastq: list[str],
) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines: list[str] = []
    lines.append("GSE189367_sample_ids.tsv-Centered SRR Comparison Report")
    lines.append("========================================================")
    lines.append(f"Generated: {stamp}")
    lines.append("")
    lines.append("Input files")
    lines.append("-----------")
    lines.append(f"- Sample list : {sample_path}")
    lines.append(f"- FASTQ dir   : {fastq_dir}")
    lines.append("")
    lines.append("Comparison basis")
    lines.append("----------------")
    lines.append(
        "- Sample export IDs were normalised to raw SRR accessions before comparison."
    )
    lines.append(
        "- FASTQ-derived IDs were extracted from filenames and counted only when both _1.fastq.gz and _2.fastq.gz were present."
    )
    lines.append("")
    lines.append("Set comparison summary")
    lines.append("----------------------")
    lines.append(f"- Sample-export SRR IDs : {len(in_both) + len(only_sample)}")
    lines.append(f"- FASTQ-derived SRR IDs : {len(in_both) + len(only_fastq)}")
    lines.append(f"- IDs in both           : {len(in_both)}")
    lines.append(f"- Only in sample export : {len(only_sample)}")
    lines.append(f"- Only in FASTQ         : {len(only_fastq)}")
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append(
        "- IDs only in sample export are expected by the sample list but missing from the current paired FASTQ set."
    )
    lines.append(
        "- IDs only in FASTQ are present on disk but absent from the sample export."
    )
    lines.append("")
    lines.append("IDs only in sample export (missing in FASTQ)")
    lines.append("---------------------------------------------")
    if only_sample:
        lines.extend(f"- {x}" for x in only_sample)
    else:
        lines.append("- none")
    lines.append("")
    lines.append("IDs only in FASTQ")
    lines.append("-----------------")
    if only_fastq:
        lines.extend(f"- {x}" for x in only_fastq)
    else:
        lines.append("- none")
    lines.append("")
    lines.append("Conclusion")
    lines.append("----------")
    lines.append(
        "- This report identifies which SRR IDs from the GSE189367 sample export are represented in the current paired FASTQ directory and which are missing."
    )

    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-list", type=Path, default=DEFAULT_SAMPLE_LIST)
    parser.add_argument("--fastq-dir", type=Path, default=DEFAULT_FASTQ_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    if not args.sample_list.exists():
        raise FileNotFoundError(f"Sample list not found: {args.sample_list}")
    if not args.fastq_dir.exists():
        raise FileNotFoundError(f"FASTQ directory not found: {args.fastq_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = extract_srr_from_sample_export(args.sample_list)
    fastq_ids = load_srr_ids_from_fastq_dir(args.fastq_dir)

    in_both = sorted(sample_ids & fastq_ids)
    only_sample = sorted(sample_ids - fastq_ids)
    only_fastq = sorted(fastq_ids - sample_ids)

    in_both_path = args.output_dir / "gse189367_ids_in_sample_and_fastq.txt"
    only_sample_path = args.output_dir / "gse189367_ids_only_in_sample_export.txt"
    only_fastq_path = args.output_dir / "gse189367_ids_only_in_fastq.txt"
    report_path = args.output_dir / "gse189367_sample_vs_fastq_report.txt"

    write_lines(in_both_path, in_both)
    write_lines(only_sample_path, only_sample)
    write_lines(only_fastq_path, only_fastq)
    report_path.write_text(
        build_report(
            sample_path=args.sample_list,
            fastq_dir=args.fastq_dir,
            in_both=in_both,
            only_sample=only_sample,
            only_fastq=only_fastq,
        ),
        encoding="utf-8",
    )

    print(f"Wrote: {in_both_path}")
    print(f"Wrote: {only_sample_path}")
    print(f"Wrote: {only_fastq_path}")
    print(f"Wrote: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())