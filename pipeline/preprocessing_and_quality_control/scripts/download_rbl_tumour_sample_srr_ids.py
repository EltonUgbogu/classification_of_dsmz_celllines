#!/usr/bin/env python3
import csv
import gzip
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DatasetConfig:
    bioproject: str
    outdir: Path
    title_pattern: str
    runinfo_gsm_fields: tuple[str, ...]


DATASETS: dict[str, DatasetConfig] = {
    "GSE111168": DatasetConfig(
        bioproject="PRJNA436090",
        outdir=Path("/work/ugbogu/pipeline/data/rbl/GSE111168_primary_tumours"),
        title_pattern=r"tumor[0-9]+",
        runinfo_gsm_fields=("SampleName", "LibraryName"),
    ),
    "GSE196420": DatasetConfig(
        bioproject="PRJNA804875",
        outdir=Path("/work/ugbogu/pipeline/data/rbl/GSE196420_primary_tumours"),
        title_pattern=r".*Primary Retinoblastoma - RNA-Seq",
        runinfo_gsm_fields=("LibraryName", "SampleName"),
    ),
    "GSE268136": DatasetConfig(
        bioproject="PRJNA1114744",
        outdir=Path("/work/ugbogu/pipeline/data/rbl/GSE268136_primary_tumours"),
        title_pattern=r"RB_[0-9]+",
        runinfo_gsm_fields=("LibraryName", "SampleName"),
    ),
}


def run_shell(cmd: str) -> str:
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        text=True,
        capture_output=True,
        executable="/bin/bash",
    )
    return result.stdout


def geo_family_soft_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz"


def is_valid_runinfo_csv(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with path.open(newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
        return (
            "Run" in header
            and "LibraryLayout" in header
            and ("LibraryStrategy" in header or "Assay Type" in header)
        )
    except Exception:
        return False


def ensure_family_soft(gse: str, outdir: Path) -> Path:
    family_soft = outdir / f"{gse}_family.soft.gz"
    if family_soft.exists() and family_soft.stat().st_size > 0:
        return family_soft

    subprocess.run(
        ["wget", "-O", str(family_soft), geo_family_soft_url(gse)],
        check=True,
    )
    return family_soft


def ensure_runinfo(cfg: DatasetConfig, runinfo_path: Path) -> None:
    if is_valid_runinfo_csv(runinfo_path):
        print(f"[INFO] using existing valid RunInfo: {runinfo_path}")
        return

    print(f"[INFO] downloading RunInfo for {cfg.bioproject}")
    stdout = run_shell(f'esearch -db sra -query "{cfg.bioproject}" | efetch -format runinfo')
    runinfo_path.write_text(stdout)

    if not is_valid_runinfo_csv(runinfo_path):
        raise RuntimeError(f"Downloaded invalid RunInfo CSV: {runinfo_path}")


def extract_gsm_titles(family_soft: Path, gsm_titles_tsv: Path) -> dict[str, str]:
    gsm_titles: dict[str, str] = {}
    current_gsm: str | None = None

    with gzip.open(family_soft, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current_gsm = line.split(" = ", 1)[1]
            elif current_gsm and line.startswith("!Sample_title = "):
                gsm_titles[current_gsm] = line.split(" = ", 1)[1]

    with gsm_titles_tsv.open("w") as out:
        for gsm, title in gsm_titles.items():
            out.write(f"{gsm}\t{title}\n")

    return gsm_titles


def select_primary_gsms(
    gsm_titles: dict[str, str],
    title_pattern: str,
    primary_gsms_tsv: Path,
) -> set[str]:
    pattern = re.compile(title_pattern)
    selected = {
        gsm: title
        for gsm, title in gsm_titles.items()
        if pattern.fullmatch(title) or pattern.search(title)
    }

    with primary_gsms_tsv.open("w") as out:
        for gsm, title in selected.items():
            out.write(f"{gsm}\t{title}\n")

    return set(selected.keys())


def resolve_gsm(row: dict[str, str], candidate_fields: tuple[str, ...]) -> str:
    for field in candidate_fields:
        value = row.get(field, "").strip()
        if value:
            return value
    return ""


def get_assay(row: dict[str, str]) -> str:
    return row.get("LibraryStrategy", "").strip() or row.get("Assay Type", "").strip()


def filter_primary_rnaseq(
    runinfo_csv: Path,
    primary_gsms: set[str],
    runinfo_gsm_fields: tuple[str, ...],
    out_all_csv: Path,
    out_all_txt: Path,
    out_pe_csv: Path,
    out_pe_txt: Path,
) -> tuple[int, int]:
    all_hits: list[str] = []
    pe_hits: list[str] = []

    with runinfo_csv.open(newline="") as f, \
         out_all_csv.open("w", newline="") as all_csv_handle, \
         out_pe_csv.open("w", newline="") as pe_csv_handle:

        reader = csv.DictReader(f)
        writer_all = csv.DictWriter(all_csv_handle, fieldnames=reader.fieldnames)
        writer_pe = csv.DictWriter(pe_csv_handle, fieldnames=reader.fieldnames)
        writer_all.writeheader()
        writer_pe.writeheader()

        for row in reader:
            assay = get_assay(row)
            layout = row.get("LibraryLayout", "").strip()
            run = row.get("Run", "").strip()
            gsm = resolve_gsm(row, runinfo_gsm_fields)

            if assay != "RNA-Seq" or not run or gsm not in primary_gsms:
                continue

            all_hits.append(run)
            writer_all.writerow(row)

            if layout == "PAIRED":
                pe_hits.append(run)
                writer_pe.writerow(row)

    with out_all_txt.open("w") as out:
        for run in all_hits:
            out.write(run + "\n")

    with out_pe_txt.open("w") as out:
        for run in pe_hits:
            out.write(run + "\n")

    return len(all_hits), len(pe_hits)


def process_dataset(gse: str) -> int:
    if gse not in DATASETS:
        print(f"ERROR: unsupported dataset '{gse}'", file=sys.stderr)
        print(f"Supported: {', '.join(DATASETS)}", file=sys.stderr)
        return 1

    cfg = DATASETS[gse]
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    family_soft = ensure_family_soft(gse, cfg.outdir)
    runinfo_csv = cfg.outdir / "SraRunInfo.csv"

    gsm_titles_tsv = cfg.outdir / "gsm_titles.tsv"
    primary_gsms_tsv = cfg.outdir / "primary_gsms.tsv"

    out_all_txt = cfg.outdir / "all_primary_tumour_srr.txt"
    out_all_csv = cfg.outdir / "SraRunInfo_primary.csv"
    out_pe_txt = cfg.outdir / "all_primary_tumour_pe.txt"
    out_pe_csv = cfg.outdir / "SraRunInfo_primary_PE.csv"

    gsm_titles = extract_gsm_titles(family_soft, gsm_titles_tsv)
    primary_gsms = select_primary_gsms(gsm_titles, cfg.title_pattern, primary_gsms_tsv)

    if not primary_gsms:
        raise RuntimeError(f"No primary GSMs matched for {gse}")

    ensure_runinfo(cfg, runinfo_csv)

    total_rna, total_pe = filter_primary_rnaseq(
        runinfo_csv=runinfo_csv,
        primary_gsms=primary_gsms,
        runinfo_gsm_fields=cfg.runinfo_gsm_fields,
        out_all_csv=out_all_csv,
        out_all_txt=out_all_txt,
        out_pe_csv=out_pe_csv,
        out_pe_txt=out_pe_txt,
    )

    print(f"[INFO] {gse}")
    print(f"[INFO] primary GSMs matched: {len(primary_gsms)}")
    print(f"[INFO] tumour RNA-seq runs: {total_rna}")
    print(f"[INFO] paired-end tumour RNA-seq runs: {total_pe}")
    print(f"[INFO] wrote: {out_all_txt}")
    print(f"[INFO] wrote: {out_pe_txt}")

    return 0


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: python {Path(sys.argv[0]).name} GSE111168|GSE196420|GSE268136|all")
        return 1

    arg = sys.argv[1]

    if arg == "all":
        rc = 0
        for gse in DATASETS:
            print(f"\n=== {gse} ===")
            try:
                ret = process_dataset(gse)
                if ret != 0:
                    rc = ret
            except Exception as e:
                print(f"[ERROR] {gse}: {e}", file=sys.stderr)
                rc = 1
        return rc

    try:
        return process_dataset(arg)
    except Exception as e:
        print(f"[ERROR] {arg}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())