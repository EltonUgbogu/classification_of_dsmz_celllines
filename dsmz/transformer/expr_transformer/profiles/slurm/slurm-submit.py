#!/usr/bin/env python3
import re
import shlex
import subprocess
import sys

jobscript = sys.argv[1]

with open(jobscript, "r", encoding="utf-8") as fh:
    text = fh.read()


def extract(key: str, default: int) -> int:
    match = re.search(rf"{key}\s*=\s*([0-9]+)", text)
    return int(match.group(1)) if match else default


gpu = extract("gpu", 0)
cpus = extract("cpus", 1)
mem_mb = extract("mem_mb", 8000)
runtime = extract("runtime", 60)

partition = "gpu" if gpu > 0 else "cpu"
qos = "short" if runtime <= 720 else "normal"

cmd = [
    "sbatch",
    f"--partition={partition}",
    f"--qos={qos}",
    f"--cpus-per-task={cpus}",
    f"--mem={mem_mb}M",
    f"--time={runtime}",
]

if gpu > 0:
    cmd.append("--gres=gpu:1")

cmd += [
    "--output=/work/%u/expr_transformer/logs/%x_%j.out",
    "--error=/work/%u/expr_transformer/logs/%x_%j.err",
    jobscript,
]

print("Submitting:", " ".join(shlex.quote(arg) for arg in cmd), file=sys.stderr)
response = subprocess.check_output(cmd).decode().strip()
print(response.split()[-1])

