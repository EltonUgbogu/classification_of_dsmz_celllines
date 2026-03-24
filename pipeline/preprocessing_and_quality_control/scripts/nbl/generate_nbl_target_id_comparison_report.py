#!/usr/bin/env python3
"""
Compare TARGET-NBL sample/aliquot ID sets and write a report.

Defaults:
- sample list:  /work/ugbogu/pipeline/data/nbl/preprocessing_results/sample_id_exports/TARGET_sample_ids.tsv
- aliquot list: /work/ugbogu/pipeline/data/nbl/target_nbl/metadata/target_nbl_aliquot_ids_unique.txt
- output dir:   /work/ugbogu/pipeline/data/nbl

Outputs:
- target_nbl_ids_in_both.txt
- target_nbl_only_in_sample_list.txt
- target_nbl_only_in_aliquot_list.txt
- target_nbl_id_set_comparison_report.txt
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_SAMPLE_LIST = Path(
    "/work/ugbogu/pipeline/data/nbl/preprocessing_results/sample_id_exports/TARGET_sample_ids.tsv"
)
DEFAULT_ALIQUOT_LIST = Path(
    "/work/ugbogu/pipeline/data/nbl/target_nbl/metadata/target_nbl_aliquot_ids_unique.txt"
)
DEFAULT_OUTPUT_DIR = Path("/work/ugbogu/pipeline/data/nbl")


def load_target_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        value = raw.strip()
        if value.startswith("TARGET-"):
            ids.add(value)
    return ids


def write_lines(path: Path, values: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for item in values:
            handle.write(item + "\n")


def build_report(
    sample_path: Path,
    aliquot_path: Path,
    in_both: list[str],
    only_sample: list[str],
    only_aliquot: list[str],
) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    lines: list[str] = []
    lines.append("TARGET-NBL ID Set Comparison Report")
    lines.append("===================================")
    lines.append(f"Generated: {stamp}")
    lines.append("")
    lines.append("Input files")
    lines.append("-----------")
    lines.append(f"- Sample list : {sample_path}")
    lines.append(f"- Aliquot list: {aliquot_path}")
    lines.append("")
    lines.append("Set comparison summary")
    lines.append("----------------------")
    lines.append(f"- Sample IDs (TARGET_sample_ids.tsv style): {len(in_both) + len(only_sample)}")
    lines.append(f"- Aliquot IDs (raw metadata list): {len(in_both) + len(only_aliquot)}")
    lines.append(f"- IDs in both: {len(in_both)}")
    lines.append(f"- Only in sample list: {len(only_sample)}")
    lines.append(f"- Only in aliquot list: {len(only_aliquot)}")
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append(
        "- The sample-ID export represents the downstream analysis subset."
    )
    lines.append(
        "- The aliquot list represents the broader identifier set derived from TARGET-NBL metadata/count downloads."
    )
    lines.append(
        "- IDs present only in the aliquot list were not retained in the downstream sample-ID export."
    )
    lines.append(
        "- Their exclusion may reflect tumour-purity filtering and/or other preprocessing or analysis-subset selection steps."
    )
    lines.append("")
    lines.append("IDs only in sample list")
    lines.append("-----------------------")
    if only_sample:
        lines.extend(f"- {x}" for x in only_sample)
    else:
        lines.append("- none")
    lines.append("")
    lines.append("IDs only in aliquot list")
    lines.append("------------------------")
    if only_aliquot:
        lines.extend(f"- {x}" for x in only_aliquot)
    else:
        lines.append("- none")
    lines.append("")
    lines.append("Conclusion")
    lines.append("----------")
    lines.append(
        "- Raw metadata IDs and downstream analysis IDs were compared as related but not necessarily identical sets."
    )
    lines.append(
        "- IDs absent from TARGET_sample_ids.tsv should be interpreted as excluded from the downstream analysis subset unless a separate purity table confirms the exact exclusion reason."
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-list", type=Path, default=DEFAULT_SAMPLE_LIST)
    parser.add_argument("--aliquot-list", type=Path, default=DEFAULT_ALIQUOT_LIST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    if not args.sample_list.exists():
        raise FileNotFoundError(f"Sample list not found: {args.sample_list}")
    if not args.aliquot_list.exists():
        raise FileNotFoundError(f"Aliquot list not found: {args.aliquot_list}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = load_target_ids(args.sample_list)
    aliquot_ids = load_target_ids(args.aliquot_list)

    in_both = sorted(sample_ids & aliquot_ids)
    only_sample = sorted(sample_ids - aliquot_ids)
    only_aliquot = sorted(aliquot_ids - sample_ids)

    in_both_path = args.output_dir / "target_nbl_ids_in_both.txt"
    only_sample_path = args.output_dir / "target_nbl_only_in_sample_list.txt"
    only_aliquot_path = args.output_dir / "target_nbl_only_in_aliquot_list.txt"
    report_path = args.output_dir / "target_nbl_id_set_comparison_report.txt"

    write_lines(in_both_path, in_both)
    write_lines(only_sample_path, only_sample)
    write_lines(only_aliquot_path, only_aliquot)
    report_path.write_text(
        build_report(
            sample_path=args.sample_list,
            aliquot_path=args.aliquot_list,
            in_both=in_both,
            only_sample=only_sample,
            only_aliquot=only_aliquot,
        ),
        encoding="utf-8",
    )

    print(f"Wrote: {in_both_path}")
    print(f"Wrote: {only_sample_path}")
    print(f"Wrote: {only_aliquot_path}")
    print(f"Wrote: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
