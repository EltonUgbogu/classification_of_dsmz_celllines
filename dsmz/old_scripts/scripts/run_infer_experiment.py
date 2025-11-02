#!/usr/bin/env python3
# infer_experiment.py pipeline: validates BAM inputs, runs strandness inference via RSeQC, and summarizes results.

import os
import subprocess
import glob
import shutil
import sys

# ==== CONFIGURATION SECTION ====

bam_dir = "/home/chu25/dsmz/results/star"  # Directory containing BAM files
bed = "/home/chu25/data/hg38_bed12/hg38.refseq.bed12"  # BED12 annotation file
output_dir = "/home/chu25/dsmz/results/strandness"  # Output directory for strandness inference
strandness_check = "/home/chu25/dsmz/logs/strandness"  # Directory for logging
summary_file = os.path.join(output_dir, "summary", "strandness_summary.tsv")  # Final TSV summary
log_file = os.path.join(strandness_check, "strandness_check.log")  # Log file

# ==== SETUP ====

# Ensure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.dirname(summary_file), exist_ok=True)
os.makedirs(strandness_check, exist_ok=True)

valid_bam_files = []

with open(log_file, "w") as log:

    # ==== BED FILE CHECKS ====
    if not os.path.exists(bed):
        log.write(f"[ERROR] BED file not found: {bed}\n")
        sys.exit(1)

    if os.stat(bed).st_size == 0:
        log.write(f"[ERROR] BED file is empty: {bed}\n")
        sys.exit(1)

    if "hg38" not in os.path.basename(bed).lower():
        log.write(f"[ERROR] BED file name does not suggest GRCh38: {bed}\n")
        sys.exit(1)

    # ==== RSeQC AVAILABILITY CHECK ====
    infer_tool = shutil.which("infer_experiment.py")
    if infer_tool is None:
        log.write("[ERROR] RSeQC's infer_experiment.py not found in PATH. Activate the correct conda environment.\n")
        sys.exit(1)
    log.write(f"[INFO] infer_experiment.py found at: {infer_tool}\n")

    # ==== BAM FILE DISCOVERY ====
    bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    if not bam_files:
        log.write(f"[ERROR] No BAM files found in {bam_dir}\n")
        sys.exit(1)

    log.write(f"[INFO] {len(bam_files)} BAM files found.\n")

    # ==== BAM FILE VALIDATION LOOP ====
    for bam in bam_files:
        sample = os.path.basename(bam).split(".")[0]
        bai = bam + ".bai"

        # Check for index file, create if missing
        if not os.path.exists(bai):
            log.write(f"[WARNING] Missing BAM index for {sample}. Attempting to create it with samtools...\n")
            result = subprocess.run(["samtools", "index", bam], stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                log.write(f"[ERROR] Failed to create BAM index for {sample}: {result.stderr}\n")
                continue
            log.write(f"[INFO] Index created for {sample}\n")

        # Check for coordinate sorting in BAM header
        result = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True)
        if "@HD" not in result.stdout or "SO:coordinate" not in result.stdout:
            log.write(f"[WARNING] BAM not coordinate-sorted: {sample}. Skipping.\n")
            continue

        # Passed all checks
        valid_bam_files.append(bam)
        log.write(f"[INFO] Valid BAM: {sample} ({bam})\n")

    # ==== LOG TOTAL VALID ====
    log.write(f"\n[INFO] Total valid BAMs for inference: {len(valid_bam_files)}\n")
    for valid_bam in valid_bam_files:
        log.write(f"  - {valid_bam}\n")

# ==== STRANDNESS INFERENCE ====

for bam in valid_bam_files:
    sample = os.path.basename(bam).split(".")[0]
    out_file = os.path.join(output_dir, f"{sample}.infer.txt")

    cmd = [infer_tool, "-i", bam, "-r", bed]
    print(f"[INFO] Running: {' '.join(cmd)}")

    with open(out_file, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"[ERROR] Inference failed for {sample}")
        with open(log_file, "a") as log:
            log.write(f"[ERROR] infer_experiment.py failed on {sample}: {result.stderr}\n")
    else:
        print(f"[INFO] Inference complete for {sample}")

# ==== SUMMARIZE RESULTS ====

summary_lines = []

# Collect strandness summaries
for f in sorted(glob.glob(os.path.join(output_dir, "*.infer.txt"))):
    sample = os.path.basename(f).replace(".infer.txt", "")
    with open(f) as fh:
        lines = fh.readlines()
        strand = lines[-1].strip() if lines else "Unknown"
    summary_lines.append(f"{sample}\t{strand}")

# Write summary to TSV
with open(summary_file, "w") as f:
    f.write("Sample\tStrandness\n")
    f.write("\n".join(summary_lines))

print(f"\n[INFO] Strandness summary saved to: {summary_file}")
