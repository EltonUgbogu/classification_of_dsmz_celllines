#!/usr/bin/env python3

import json
import pandas as pd
from pathlib import Path

meta_file = Path("/work/ugbogu/pipeline/data/nbl/target_nbl/metadata.cart.2026-03-24.json")
out_file = Path("/work/ugbogu/pipeline/data/nbl/target_nbl/target_nbl_uuid_to_submitter.tsv")

with open(meta_file) as f:
    meta = json.load(f)

rows = []

for rec in meta:
    file_id = rec.get("id")
    file_name = rec.get("file_name")

    submitter_id = None

    # Best candidate: associated_entities
    for ent in rec.get("associated_entities", []):
        sid = ent.get("entity_submitter_id")
        if sid and sid.startswith("TARGET-"):
            submitter_id = sid
            break

    # Fallback: cases -> samples
    if submitter_id is None:
        for case in rec.get("cases", []):
            for sample in case.get("samples", []):
                sid = sample.get("submitter_id")
                if sid and sid.startswith("TARGET-"):
                    submitter_id = sid
                    break
            if submitter_id:
                break

    rows.append({
        "file_id": file_id,
        "file_name": file_name,
        "submitter_id": submitter_id
    })

df = pd.DataFrame(rows)
df.to_csv(out_file, sep="\t", index=False)

print(df.head(20).to_string(index=False))
print()
print(f"Saved: {out_file}")
print(f"Total rows: {len(df)}")
print(f"Mapped TARGET submitter IDs: {df['submitter_id'].notna().sum()}")