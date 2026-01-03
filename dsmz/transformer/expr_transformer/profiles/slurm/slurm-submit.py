#!/usr/bin/env python3
# slurm-submit.py - SLURM job submission script for the expr_transformer project
# Compatible with Snakemake 5.10 (classic cluster interface)
import os
import re
import shlex
import subprocess
import sys


def extract_int(key: str, text: str, default: int) -> int:
    """
    Extract integer values from Snakemake jobscript text.

    Supports:
      - key=123
      - key: 123
      - "key": 123
      - 'key': 123
    """
    patterns = [
        rf"{re.escape(key)}\s*=\s*([0-9]+)",
        rf"{re.escape(key)}\s*:\s*([0-9]+)",
        rf'"{re.escape(key)}"\s*:\s*([0-9]+)',
        rf"'{re.escape(key)}'\s*:\s*([0-9]+)",
    ]
    for pat in patterns:
        m = re.search(pat, text)
        if m:
            return int(m.group(1))
    return default


def main():
    if len(sys.argv) != 2:
        print("Usage: slurm-submit.py <jobscript>", file=sys.stderr)
        return 2

    jobscript = sys.argv[1]

    try:
        with open(jobscript, "r", encoding="utf-8") as fh:
            text = fh.read()
    except FileNotFoundError:
        print(f"ERROR: jobscript not found: {jobscript}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: failed to read jobscript: {e}", file=sys.stderr)
        return 1

    # Extract resources from jobscript
    gpu = extract_int("gpu", text, 0)
    threads = extract_int("threads", text, 1)  # Snakemake 5 uses 'threads' not 'cpus'
    mem_mb = extract_int("mem_mb", text, 8000)
    runtime = extract_int("runtime", text, 60)

    # Safety override: force GPU for known GPU rules (belt-and-suspenders)
    force_gpu = ("train_mgm" in text) or ("export_embeddings" in text)
    if force_gpu:
        gpu = max(gpu, 1)

    # Determine partition and QOS based on resources
    partition = "gpu" if gpu > 0 else "cpu"
    qos = "short" if runtime <= 720 else "normal"

    # Convert runtime from minutes to HH:MM:SS format for Slurm
    hours = runtime // 60
    mins = runtime % 60
    time_str = f"{hours:02d}:{mins:02d}:00"

    # Build sbatch command
    cmd = [
        "sbatch",
        "--parsable",  # Print only jobid (cleaner for Snakemake)
        f"--partition={partition}",
        f"--qos={qos}",
        f"--cpus-per-task={threads}",
        f"--mem={mem_mb}M",
        f"--time={time_str}",
    ]

    if gpu > 0:
        cmd.append("--gres=gpu:1")

    # Use explicit user instead of %u (more reliable)
    user = os.environ.get("USER", "unknown")
    cmd += [
        f"--output=/work/{user}/expr_transformer/logs/%x_%j.out",
        f"--error=/work/{user}/expr_transformer/logs/%x_%j.err",
        jobscript,
    ]

    # Submit job
    print("Submitting:", " ".join(shlex.quote(arg) for arg in cmd), file=sys.stderr)
    
    try:
        res = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            sys.stderr.write(res.stderr)
            return res.returncode
        
        # Snakemake expects the jobid on stdout
        jobid = res.stdout.strip()
        if jobid:
            sys.stdout.write(jobid + "\n")
            return 0
        else:
            print("ERROR: sbatch returned empty jobid", file=sys.stderr)
            return 1
    except Exception as e:
        print(f"ERROR: submission failed: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
