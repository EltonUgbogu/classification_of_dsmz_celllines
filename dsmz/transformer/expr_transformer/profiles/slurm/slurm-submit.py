#!/usr/bin/env python3
# slurm-submit.py - SLURM job submission script for the expr_transformer project
# Compatible with Snakemake 5.10 (classic cluster interface)
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from typing import Any, Dict

DEFAULTS: Dict[str, int] = {
    "gpu": 0,
    "cpus": 1,
    "mem_mb": 8000,
    "runtime": 60,  # minutes
}

PROPS_RE = re.compile(r"^#\s*properties\s*=\s*(\{.*\})\s*$")


def minutes_to_hms(minutes: int) -> str:
    """Convert minutes to HH:MM:SS format for Slurm."""
    hours, mins = divmod(max(int(minutes), 1), 60)
    return f"{hours:02d}:{mins:02d}:00"


def read_props_from_jobscript(jobscript: str) -> Dict[str, Any]:
    """Read the # properties = {...} line from the jobscript and parse as JSON."""
    try:
        with open(jobscript, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = PROPS_RE.match(line.strip())
                if m:
                    return json.loads(m.group(1))
    except Exception as e:
        print(f"ERROR: failed to read/parse jobscript properties: {e}", file=sys.stderr)
    return {}


def get_int(d: Dict[str, Any], key: str, default: int) -> int:
    """Safely extract integer from dict with fallback to default."""
    try:
        return int(d.get(key, default))
    except Exception:
        return int(default)


def main() -> None:
    print("### USING THIS slurm-submit.py:", __file__, file=sys.stderr)
    
    if len(sys.argv) != 2:
        print("Usage: slurm-submit.py <jobscript>", file=sys.stderr)
        sys.exit(2)

    jobscript = sys.argv[1]
    job_props = read_props_from_jobscript(jobscript)

    resources = job_props.get("resources", {}) or {}
    threads = get_int(job_props, "threads", 1)

    gpu = get_int(resources, "gpu", DEFAULTS["gpu"])
    mem_mb = get_int(resources, "mem_mb", DEFAULTS["mem_mb"])
    runtime = get_int(resources, "runtime", DEFAULTS["runtime"])

    # Prefer threads for cpus-per-task; fall back to resources.cpus if you use it
    cpus = max(threads, get_int(resources, "cpus", DEFAULTS["cpus"]))

    # Debug: show extracted resources
    resources_dict = {
        "gpu": gpu,
        "threads": threads,
        "cpus": cpus,
        "mem_mb": mem_mb,
        "runtime": runtime
    }
    print("### JOB RESOURCES:", resources_dict, file=sys.stderr)

    # Determine partition and QOS based on resources
    # QOS threshold: short <= 60 min, normal > 60 min (matches typical HPC cluster policies)
    partition = "gpu" if gpu > 0 else "cpu"
    qos = "short" if runtime <= 60 else "normal"

    # Build sbatch command
    cmd = [
        "sbatch",
        "--parsable",  # Print only jobid (cleaner for Snakemake)
        f"--partition={partition}",
        f"--qos={qos}",
        f"--cpus-per-task={cpus}",
        f"--mem={mem_mb}M",
        f"--time={minutes_to_hms(runtime)}",
    ]

    if gpu > 0:
        cmd.append(f"--gres=gpu:{gpu}")

    # Use explicit user path (more reliable than %u)
    user = "ugbogu"  # Hardcoded for this cluster setup
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
            sys.exit(res.returncode)
        
        # Snakemake expects the jobid on stdout
        jobid = res.stdout.strip()
        if jobid:
            # Handle parsable output (may have multiple fields, take last)
            print(jobid.split()[-1])
            sys.exit(0)
        else:
            print("ERROR: sbatch returned empty jobid", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"ERROR: submission failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
