#!/usr/bin/env python3
"""
Check coverage for all BAM files in a directory using samtools idxstats.
Reports total mapped reads per BAM file to a TSV summary.
"""

import os
import subprocess
import glob

# ==== CONFIGURATION (Edit these paths) ====
bam_dir = "/home/chu25/dsmz/results/star"  # Directory containing BAM files
output_path = "/home/chu25/dsmz/results/strandness/bam_coverage.tsv"  # Output summary file
min_reads = 1_000_000  # Minimum acceptable mapped reads per sample

# ==== SCAN FOR BAM FILES ====
bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))

if not bam_files:
    print(f"[ERROR] No BAM files found in {bam_dir}")
    exit(1)

# ==== RUN COVERAGE CHECK ====
with open(output_path, "w") as out:
    out.write("sample\tmapped_reads\tstatus\n")  # Header row

    for bam in bam_files:
        sample = os.path.basename(bam).split(".")[0]

        try:
            # samtools idxstats returns one line per reference:
            # ref_name, ref_length, mapped_reads, unmapped_reads
            result = subprocess.run(["samtools", "idxstats", bam],
                                    capture_output=True, text=True, check=True)

            # Parse and sum mapped reads (3rd column)
            mapped = sum(int(line.split("\t")[2]) for line in result.stdout.strip().split("\n") if line)

            # Determine status based on threshold
            status = "OK" if mapped >= min_reads else "LOW_COVERAGE"

            # Write result to TSV
            out.write(f"{sample}\t{mapped}\t{status}\n")
            print(f"[INFO] {sample}: {mapped} mapped reads ({status})")

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] samtools idxstats failed for {sample}: {e}")
            out.write(f"{sample}\tNA\tFAILED\n")
        except Exception as ex:
            print(f"[ERROR] Unexpected error for {sample}: {ex}")
            out.write(f"{sample}\tNA\tERROR\n")

print(f"\n[INFO] Coverage summary written to: {output_path}")
