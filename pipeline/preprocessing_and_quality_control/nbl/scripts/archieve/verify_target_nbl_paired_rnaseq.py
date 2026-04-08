#!/usr/bin/env python3
"""
verify_target_nbl_paired_rnaseq.py

Verifies whether TARGET-NBL primary tumour RNA-seq samples were sequenced as
paired-end reads by querying the NCI Genomic Data Commons (GDC) file metadata API.

The script proceeds in four sequential steps:
  1. Filter construction  — builds a query restricting results to TARGET-NBL
                            primary tumour RNA-seq sequencing-read files.
  2. Data retrieval       — POSTs the filter to the GDC /files endpoint and
                            collects raw file metadata records.
  3. Record flattening    — denormalises the nested file→case→sample structure
                            into flat rows, one per (file, sample) association.
  4. Pairing inference    — groups rows by sample and checks whether both mate 1
                            and mate 2 can be identified via explicit GDC
                            annotation or filename pattern matching. A sample is
                            classified as proven paired only when bilateral
                            evidence exists; otherwise it is marked ambiguous.

Output validity: the binary classification is conservative — a sample is
promoted to "proven paired" only when both mate designations are positively
found. Missing evidence is never equated with single-end sequencing.
"""

from __future__ import annotations

import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# GDC RESTful endpoint for file catalogue queries.
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"

# Metadata fields requested from the API. Dot-notation paths instruct the GDC
# to include nested case and sample annotations alongside each file record;
# these are required for per-sample grouping in Step 4.
FIELDS = [
    "file_id",
    "file_name",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "read_pair_number",           # Explicit GDC paired-end annotation (Step 4, strategy a)
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.project.project_id",
]

# Regex patterns for inferring mate-1 identity from a filename when
# read_pair_number is absent (Step 4, strategy b). Three families cover the
# most common FASTQ naming conventions: _R1_, _1_, and _1.fastq.gz suffixes.
R1_PATTERNS = [
    re.compile(r"(^|[._-])R1([._-]|$)", re.IGNORECASE),
    re.compile(r"(^|[._-])1([._-]|$)", re.IGNORECASE),
    re.compile(r"_1\.f(ast)?q(\.gz)?$", re.IGNORECASE),
]

# Corresponding patterns for mate-2 identity (Step 4, strategy b).
R2_PATTERNS = [
    re.compile(r"(^|[._-])R2([._-]|$)", re.IGNORECASE),
    re.compile(r"(^|[._-])2([._-]|$)", re.IGNORECASE),
    re.compile(r"_2\.f(ast)?q(\.gz)?$", re.IGNORECASE),
]


def post_json(url: str, payload: dict[str, Any], retries: int = 3) -> dict[str, Any]:
    """
    Serialises payload as JSON and POSTs it to the given URL, returning the
    decoded response dictionary.

    A retry loop with linear back-off (2 s × attempt index) tolerates
    transient HTTP or network failures; a RuntimeError is raised only once
    all attempts are exhausted.
    """
    data = json.dumps(payload).encode("utf-8")
    req = Request(
        url,
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    for attempt in range(1, retries + 1):
        try:
            with urlopen(req, timeout=180) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (HTTPError, URLError) as exc:
            if attempt == retries:
                raise RuntimeError(f"GDC API request failed: {exc}") from exc
            time.sleep(2 * attempt)
    raise RuntimeError("Unexpected API failure")


def build_filters() -> dict[str, Any]:
    """
    Constructs the Step 1 GDC filter as a four-way AND conjunction:
      - project = TARGET-NBL          (neuroblastoma cohort only)
      - sample type = Primary Tumor   (excludes normals, metastases, cell lines)
      - data category = Sequencing Reads  (raw reads, not derived alignments)
      - experimental strategy = RNA-Seq   (transcriptomic assays only)

    All four conditions must be satisfied simultaneously, ensuring that only
    the relevant file subset is retrieved before any local processing begins.
    """
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
                    "field": "cases.samples.sample_type",
                    "value": ["Primary Tumor"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": ["Sequencing Reads"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"],
                },
            },
        ],
    }


def fetch_all_hits() -> list[dict[str, Any]]:
    """
    Executes Step 2: assembles the complete API payload and retrieves file
    records from the GDC.

    A page size of 5 000 exceeds the expected number of TARGET-NBL primary
    tumour RNA-seq files, so the full result set is obtained in a single
    request. Returns the raw list of hit dictionaries from the GDC response 
    envelope.
    """
    payload = {
        "filters": build_filters(),
        "fields": ",".join(FIELDS),
        "format": "JSON",
        "size": "5000",
    }
    result = post_json(GDC_FILES_ENDPOINT, payload)
    return result.get("data", {}).get("hits", [])


def flatten_hits(hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Performs Step 3: denormalises the nested GDC structure
    (file → cases → samples) into a flat list of row dictionaries.

    Because one file may be linked to multiple cases and each case to multiple
    samples, nested loops emit one row per (file, case, sample) combination.
    This expansion is necessary so that Step 4 can group rows by sample
    identifier. Missing fields are stored as None rather than raising errors.
    """
    rows: list[dict[str, Any]] = []
    for hit in hits:
        cases = hit.get("cases", []) or []
        for case in cases:
            samples = case.get("samples", []) or []
            for sample in samples:
                rows.append(
                    {
                        "file_id": hit.get("file_id"),
                        "file_name": hit.get("file_name"),
                        "data_category": hit.get("data_category"),
                        "data_type": hit.get("data_type"),
                        "data_format": hit.get("data_format"),
                        "experimental_strategy": hit.get("experimental_strategy"),
                        "read_pair_number": hit.get("read_pair_number"),
                        "case_id": case.get("case_id"),
                        "case_submitter_id": case.get("submitter_id"),
                        "sample_id": sample.get("sample_id"),
                        "sample_submitter_id": sample.get("submitter_id"),
                        "sample_type": sample.get("sample_type"),
                        "project_id": (case.get("project") or {}).get("project_id"),
                    }
                )
    return rows


def normalise_pair_number(value: Any) -> int | None:
    """
    Converts the GDC read_pair_number field to an integer in {1, 2}, or
    returns None if the value is absent or does not encode a valid mate
    designation. The field may arrive as an integer, a string, or None.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s in {"1", "2"}:
        return int(s)
    return None


def infer_pair_from_filename(name: str | None) -> int | None:
    """
    Applies the Step 4 strategy-b fallback: scans the filename against the
    R1_PATTERNS list first; if a match is found, mate 1 is inferred. The
    R2_PATTERNS list is tested only when no R1 pattern matches. Returns None
    when neither set of patterns produces a match, leaving the file
    unclassified rather than guessing.
    """
    if not name:
        return None
    for pat in R1_PATTERNS:
        if pat.search(name):
            return 1
    for pat in R2_PATTERNS:
        if pat.search(name):
            return 2
    return None


def main() -> int:
    # --- Step 2 & 3: retrieve and flatten GDC records ---
    hits = fetch_all_hits()
    rows = flatten_hits(hits)

    if not rows:
        print("No TARGET-NBL primary tumour RNA-seq sequencing-read files were returned.")
        return 1

    # --- Step 4a: group all files by (case_id, sample_id) ---
    # Each key represents one unique biological sample; its value list holds
    # all associated file records. Records lacking either identifier are
    # discarded to prevent spurious grouping under a None key.
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        case_id = row.get("case_id")
        sample_id = row.get("sample_id")
        if case_id and sample_id:
            grouped[(case_id, sample_id)].append(row)

    proven_paired = []
    ambiguous = []

    # --- Step 4b: classify each sample as proven paired or ambiguous ---
    # For every file in a sample group, the mate designation is resolved by
    # first checking the explicit GDC annotation (strategy a) and falling back
    # to filename pattern matching (strategy b) only when the annotation is
    # absent. The resolved designations are collected in a set; a sample
    # reaches "proven paired" only when the set contains both 1 and 2.
    for key, recs in grouped.items():
        pair_nums = set()
        filenames = []

        for rec in recs:
            p = normalise_pair_number(rec.get("read_pair_number"))  # Strategy a
            if p is None:
                p = infer_pair_from_filename(rec.get("file_name"))  # Strategy b fallback
            if p in {1, 2}:
                pair_nums.add(p)
            filenames.append(rec.get("file_name"))

        if pair_nums == {1, 2}:
            # Both mate designations positively identified — sample is paired-end.
            proven_paired.append((key, filenames))
        else:
            # One or both designations missing — sample is not provably paired.
            ambiguous.append((key, filenames, sorted(pair_nums)))

    # --- Output: summary counts and diagnostic examples for ambiguous cases ---
    print("=== TARGET-NBL RNA-seq pairing verification ===")
    print(f"Unique primary tumour RNA-seq samples examined: {len(grouped)}")
    print(f"Samples proven paired-end: {len(proven_paired)}")
    print(f"Samples ambiguous / not provable as paired from returned metadata: {len(ambiguous)}")

    if ambiguous:
        print("\nExamples of ambiguous samples (first 10):")
        for (case_id, sample_id), filenames, pair_nums in ambiguous[:10]:
            print(f"- case={case_id} sample={sample_id} pair_markers={pair_nums}")
            for name in filenames[:4]:
                print(f"    {name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())