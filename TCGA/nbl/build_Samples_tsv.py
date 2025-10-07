#!/usr/bin/env python3
# Build per-project samples.tsv files (PAIRED-END ONLY, strict)
# - Scans <data_dir>/<PROJECT>/fastq/
# - Requires matching R1/R2 for every sample (multi-lane supported)
# - If any sample is SINGLE or missing mates -> hard error and exit(1)

# Import standard library modules for file system operations, regular expressions, system interface, and data manipulation.
import os
import re
import sys
import pandas as pd
# Import defaultdict from collections for grouping files by sample base name.
from collections import defaultdict

# === CONFIG ===
# Define the root directory containing all project data.
data_dir = "/home/chu25/data/nbl"
# List of project IDs to process. This should be updated to include all relevant projects.
PROJECTS = ["GSE100148", "GSE189367", "GSE49711", "SRP409177"]

# --------- helpers ---------
# Define a natural sorting key function to handle alphanumeric filenames correctly (e.g., for lane numbers).
def natural_key(s):
    # Split the string by digits and non-digits, converting digits to integers for proper sorting.
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

# Define regex patterns for identifying paired-end FASTQ filenames from different sequencing platforms/naming conventions.
# Each tuple contains a compiled regex and a descriptive name for the pattern.
READ_PATTERNS = [
    # Illumina standard pattern: supports sample_S<set>_L<lane>_R<read>.[fastq|fq][.gz]
    (re.compile(r'(.+?)(_S\d+)?(_L\d{3})?_R([12])(_\d+)?\.(?:fastq|fq)(?:\.gz)?$', re.I), 'illumina_std'),
    # Pattern for filenames ending in .R1. or .R2.[fastq|fq][.gz]
    (re.compile(r'(.+?)\.R([12])\.(?:fastq|fq)(?:\.gz)?$', re.I), 'dotR'),
    # Pattern for filenames with _R1. or _R2.[fastq|fq][.gz]
    (re.compile(r'(.+?)_R([12])\.(?:fastq|fq)(?:\.gz)?$', re.I), 'underscoreR'),
    # Pattern for filenames with _1. or _2.[fastq|fq][.gz] (simple underscore numbering).
    (re.compile(r'(.+?)_([12])\.(?:fastq|fq)(?:\.gz)?$', re.I), 'underscore_num'),
]

# Define regex pattern for identifying single-end FASTQ filenames.
SE_PATTERNS = [
    # Basic pattern for any filename ending in .[fastq|fq][.gz] without read indicators.
    re.compile(r'^(.+?)\.(?:fastq|fq)(?:\.gz)?$', re.I),
]

# Parse a FASTQ filename to extract the sample base name and read identifier (R1, R2, or None for SE).
def parse_fastq_name(fname):
    # Iterate over paired-end patterns to match the filename.
    for pat, _ in READ_PATTERNS:
        # Attempt to match the pattern against the filename.
        m = pat.match(fname)
        if m:
            # Extract the base sample name (group 1).
            base = m.group(1)
            # Extract the read number (group 4 for most patterns, group 2 for dotR).
            read = m.group(4) if pat.groups >= 4 else m.group(2)
            # Return base name and read identifier as "R{read}".
            return base, f"R{read}"
    # If no paired-end match, check single-end patterns.
    for pat in SE_PATTERNS:
        # Attempt to match the single-end pattern.
        m = pat.match(fname)
        if m:
            # Return base name and None for single-end.
            return m.group(1), None
    # If no match, return None for both (unparseable file).
    return None, None

# Load an optional allowlist of sample prefixes from text files in the project root.
def load_allowlist(project_root):
    # Initialize an empty set for allowed sample prefixes.
    allow = set()
    # Iterate over all files in the project root.
    for fn in os.listdir(project_root):
        # Check if the file is a text file (case-insensitive).
        if fn.lower().endswith(".txt"):
            # Open and read the file.
            with open(os.path.join(project_root, fn)) as fh:
                # Process each line.
                for line in fh:
                    # Strip whitespace and skip empty lines or comments.
                    s = line.strip()
                    if s and not s.startswith("#"):
                        # Add the stripped string to the allowlist set.
                        allow.add(s)
    # Return the set of allowed prefixes.
    return allow

# Build the samples.tsv file for a single project, enforcing strict paired-end requirements.
def build_samples_tsv_for_project(project_id):
    # Construct the project root directory path.
    proj_root = os.path.join(data_dir, project_id)
    # Construct the FASTQ input directory path.
    fastq_dir = os.path.join(proj_root, "fastq")
    # Construct the config output directory path.
    cfg_dir   = os.path.join(proj_root, "config")
    # Create the config directory if it doesn't exist.
    os.makedirs(cfg_dir, exist_ok=True)

    # Check if the FASTQ directory exists; if not, log an error and return False.
    if not os.path.isdir(fastq_dir):
        print(f"[ERROR] Project {project_id}: missing fastq/ dir: {fastq_dir}", file=sys.stderr)
        return False

    # Load the allowlist for this project if any text files are present.
    allow = load_allowlist(proj_root)
    # Initialize a defaultdict to bucket FASTQ files by sample base name, with lists for R1, R2, and SE.
    buckets = defaultdict(lambda: {"R1": [], "R2": [], "SE": []})

    # Iterate over files in the FASTQ directory, sorted naturally for consistent ordering.
    for fn in sorted(os.listdir(fastq_dir), key=natural_key):
        # Skip files that do not end in .fastq, .fq, .fastq.gz, or .fq.gz (case-insensitive).
        if not re.search(r'\.(fastq|fq)(\.gz)?$', fn, flags=re.I):
            continue
        # Parse the filename to extract base and read type.
        base, read = parse_fastq_name(fn)
        # If parsing fails, log a warning and skip the file.
        if base is None:
            print(f"[WARN] {project_id}: Could not parse read from {fn}, skipping", file=sys.stderr)
            continue

        # Apply allowlist filtering: skip if base doesn't match or start with an allowed prefix.
        if allow and not (base in allow or any(base.startswith(srr) for srr in allow)):
            continue

        # Construct the full path to the FASTQ file.
        fpath = os.path.join(fastq_dir, fn)
        # Categorize and append the file path to the appropriate bucket list.
        if read == "R1":
            buckets[base]["R1"].append(fpath)
        elif read == "R2":
            buckets[base]["R2"].append(fpath)
        else:
            buckets[base]["SE"].append(fpath)

    # Initialize lists to collect errors and valid rows for the DataFrame.
    errors = []
    rows = []

    # Process each bucketed sample.
    for base, parts in buckets.items():
        # Sort file lists naturally for consistent multi-lane ordering.
        r1_list = sorted(parts["R1"], key=natural_key)
        r2_list = sorted(parts["R2"], key=natural_key)
        se_list = sorted(parts["SE"], key=natural_key)

        # Strict rule: reject samples that are purely single-end (no R1 or R2).
        if se_list and not (r1_list or r2_list):
            errors.append(f"{project_id}: {base} appears SINGLE-END ({len(se_list)} files).")
            continue

        # Strict rule: require both R1 and R2 lists to be non-empty.
        if not r1_list or not r2_list:
            errors.append(f"{project_id}: {base} missing mate files (R1:{len(r1_list)} / R2:{len(r2_list)}).")
            continue

        # Strict rule: ensure equal number of lanes/files for R1 and R2.
        if len(r1_list) != len(r2_list):
            errors.append(f"{project_id}: {base} lane count mismatch (R1:{len(r1_list)} vs R2:{len(r2_list)}).")
            continue

        # If all checks pass, create a row dictionary for the sample.
        rows.append({
            "sample": base,
            "r1": ",".join(r1_list),  # Comma-separated list of R1 file paths.
            "r2": ",".join(r2_list),  # Comma-separated list of R2 file paths.
            "layout": "PAIRED"         # Fixed layout indicator for paired-end.
        })

    # If any errors were found, print them and return False to indicate failure.
    if errors:
        print("\n[PAIRED-END ENFORCEMENT] The following issues were found:", file=sys.stderr)
        for e in errors:
            print(" - " + e, file=sys.stderr)
        print("\nAborting. Fix FASTQs or adjust allowlist.", file=sys.stderr)
        return False

    # If no valid rows, log an error and return False.
    if not rows:
        print(f"[ERROR] {project_id}: No valid paired-end samples found.", file=sys.stderr)
        return False

    # Create a DataFrame from valid rows and sort by sample name.
    df = pd.DataFrame(rows).sort_values("sample")
    # Define the output TSV path in the config directory.
    out_tsv = os.path.join(cfg_dir, "samples.tsv")
    # Write the DataFrame to TSV without index.
    df.to_csv(out_tsv, sep="\t", index=False)
    # Log success message with sample count.
    print(f"[OK] {project_id}: wrote {out_tsv} with {len(df)} paired samples")
    # Return True to indicate success.
    return True

# Define the main entry point function.
def main():
    # Initialize success flag.
    ok = True
    # Process each project in the list.
    for project in PROJECTS:
        # Update success flag based on individual project result.
        ok = build_samples_tsv_for_project(project) and ok
    # Exit with 0 on success, 1 on any failure.
    sys.exit(0 if ok else 1)

# Standard Python idiom to run main() if the script is executed directly.
if __name__ == "__main__":
    main()
