#!/usr/bin/env python3
import argparse
import gzip
import re
import pandas as pd

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def parse_attr(attr: str) -> dict:
    # GTF attributes look like: key "value"; key2 "value2";
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        out[m.group(1)] = m.group(2)
    return out

def strip_version(x: str) -> str:
    return x.split(".", 1)[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="Path to hg38 GTF (can be .gz)")
    ap.add_argument("--out", required=True, help="Output TSV mapping")
    args = ap.parse_args()

    rows = []
    with open_text(args.gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "gene":
                continue
            attrs = parse_attr(parts[8])
            gid = attrs.get("gene_id")
            gname = attrs.get("gene_name") or attrs.get("Name")
            gtype = attrs.get("gene_type") or attrs.get("gene_biotype")
            if not gid:
                continue
            rows.append((strip_version(gid), gname, gtype))

    df = pd.DataFrame(rows, columns=["ensembl_id", "symbol", "gene_type"]).drop_duplicates()
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote {args.out} ({len(df)} rows)")

if __name__ == "__main__":
    main()
