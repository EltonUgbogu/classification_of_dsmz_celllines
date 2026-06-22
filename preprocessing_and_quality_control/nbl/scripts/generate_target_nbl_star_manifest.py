#!/usr/bin/env python3
"""
Generate a GDC manifest for TARGET-NBL primary-tumour RNA-seq STAR-count files.

Filters:
- project.project_id = TARGET-NBL
- data_category = Transcriptome Profiling
- data_type = Gene Expression Quantification
- analysis.workflow_type = STAR - Counts
- cases.samples.sample_type = Primary Tumor

Output:
- manifest.tsv under `NBL_DATA_ROOT/target_nbl` by default.

Usage:
    python3 generate_target_nbl_star_manifest.py
    python3 generate_target_nbl_star_manifest.py /path/to/output_manifest.tsv
"""

from __future__ import annotations

import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"

REPO_ROOT = Path(__file__).resolve().parents[3]
NBL_DATA_ROOT = Path(os.environ.get("NBL_DATA_ROOT", REPO_ROOT / "data" / "nbl")).expanduser()
DEFAULT_MANIFEST_PATH = NBL_DATA_ROOT / "target_nbl" / "manifest.tsv"


def build_filters() -> dict:
    return {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["TARGET-NBL"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": ["Transcriptome Profiling"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": ["Gene Expression Quantification"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": ["STAR - Counts"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.sample_type",
                    "value": ["Primary Tumor"],
                },
            },
        ],
    }


def main() -> int:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MANIFEST_PATH
    out_path.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "filters": json.dumps(build_filters()),
        "return_type": "manifest",
    }

    url = f"{GDC_FILES_ENDPOINT}?{urllib.parse.urlencode(params)}"

    print(f"[INFO] Requesting manifest from: {GDC_FILES_ENDPOINT}", file=sys.stderr)
    print(f"[INFO] Writing to: {out_path}", file=sys.stderr)

    req = urllib.request.Request(url, method="GET")
    try:
        with urllib.request.urlopen(req, timeout=300) as response:
            manifest_bytes = response.read()
    except Exception as exc:
        print(f"ERROR: Failed to retrieve manifest: {exc}", file=sys.stderr)
        return 1

    if not manifest_bytes.strip():
        print("ERROR: Empty manifest returned.", file=sys.stderr)
        return 1

    out_path.write_bytes(manifest_bytes)
    print(f"[INFO] Manifest written: {out_path}", file=sys.stderr)

    # Quick sanity check
    try:
        with out_path.open("r", encoding="utf-8", errors="replace") as handle:
            line_count = sum(1 for _ in handle)
        print(f"[INFO] Manifest lines: {line_count}", file=sys.stderr)
    except Exception as exc:
        print(f"[WARN] Could not count manifest lines: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
