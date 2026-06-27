#!/usr/bin/env python3

import csv
import json
import os
import urllib.request
from pathlib import Path

PROJECT = "TCGA-BRCA"
REPO_ROOT = Path(__file__).resolve().parents[2]
BRCA_DATA_ROOT = Path(os.environ.get("BRCA_DATA_ROOT", REPO_ROOT / "data" / "brca" / "new_data")).expanduser()
OUTDIR = BRCA_DATA_ROOT / "TCGA_BRCA"
MANIFEST = OUTDIR / "manifests" / "gdc_manifest_tcga_brca_star_counts_primary_tumour.txt"
META = OUTDIR / "metadata" / "tcga_brca_star_counts_primary_tumour_files.tsv"

API = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT]}},
        {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
        {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
    ],
}

fields = ",".join([
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "state",
    "data_category",
    "data_type",
    "experimental_strategy",
    "analysis.workflow_type",
    "cases.submitter_id",
    "cases.case_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.sample_type",
    "cases.samples.portions.analytes.aliquots.submitter_id",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
])

def post_json(payload):
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        API,
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))

# First request gets total
payload = {
    "filters": filters,
    "fields": fields,
    "format": "JSON",
    "size": 1,
    "from": 0,
}
first = post_json(payload)
total = int(first["data"]["pagination"]["total"])
print(f"Total matching files: {total}")

rows = []
page_size = 500
for start in range(0, total, page_size):
    payload["size"] = page_size
    payload["from"] = start
    out = post_json(payload)
    rows.extend(out["data"]["hits"])

OUTDIR.joinpath("manifests").mkdir(parents=True, exist_ok=True)
OUTDIR.joinpath("metadata").mkdir(parents=True, exist_ok=True)

# GDC manifest format
with MANIFEST.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["id", "filename", "md5", "size", "state"])
    for r in rows:
        w.writerow([
            r.get("file_id", ""),
            r.get("file_name", ""),
            r.get("md5sum", ""),
            r.get("file_size", ""),
            r.get("state", ""),
        ])

def flatten_case_sample(row):
    case_submitter = ""
    sample_submitter = ""
    sample_type = ""
    aliquot_submitter = ""

    cases = row.get("cases") or []
    if cases:
        c = cases[0]
        case_submitter = c.get("submitter_id", "")
        samples = c.get("samples") or []
        if samples:
            s = samples[0]
            sample_submitter = s.get("submitter_id", "")
            sample_type = s.get("sample_type", "")
            portions = s.get("portions") or []
            if portions:
                analytes = portions[0].get("analytes") or []
                if analytes:
                    aliquots = analytes[0].get("aliquots") or []
                    if aliquots:
                        aliquot_submitter = aliquots[0].get("submitter_id", "")

    return case_submitter, sample_submitter, sample_type, aliquot_submitter

with META.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "file_id", "file_name", "md5sum", "file_size", "state",
        "data_category", "data_type", "experimental_strategy", "workflow_type",
        "case_submitter_id", "sample_submitter_id", "sample_type", "aliquot_submitter_id",
    ])
    for r in rows:
        case_submitter, sample_submitter, sample_type, aliquot_submitter = flatten_case_sample(r)
        workflow = ""
        if isinstance(r.get("analysis"), dict):
            workflow = r["analysis"].get("workflow_type", "")
        w.writerow([
            r.get("file_id", ""),
            r.get("file_name", ""),
            r.get("md5sum", ""),
            r.get("file_size", ""),
            r.get("state", ""),
            r.get("data_category", ""),
            r.get("data_type", ""),
            r.get("experimental_strategy", ""),
            workflow,
            case_submitter,
            sample_submitter,
            sample_type,
            aliquot_submitter,
        ])

print(f"Wrote manifest: {MANIFEST}")
print(f"Wrote metadata: {META}")
