# Import necessary libraries
import pandas as pd                      # for handling tabular data and writing TSV
from collections import defaultdict      # for grouping R1 and R2 by sample prefix
import os                                # for file and path manipulation

# Define the root directory where all FASTQ files are stored
data_dir = "/home/cpo14/Dokumente/ngs/data/gatc"

# Path to a text file that contains relative paths to the FASTQ files
txt_path = "/home/chu25/dsmz/scripts/sorted_fastq_files.txt"

# Step 1: Read all file paths from the text file and make them absolute by prepending data_dir
with open(txt_path) as f:
    
    paths = [
        os.path.normpath(os.path.join(data_dir, line.strip().lstrip("/")))  # Remove leading slash and join with data_dir
        for line in f
        if line.strip().endswith(".fastq.gz")  # Only include FASTQ files
    ]

# Step 2: Prepare a dictionary to group R1 and R2 FASTQ files by sample prefix
# Each sample will map to a dictionary with 'r1' and 'r2' keys, each holding a list of paths
paired = defaultdict(lambda: {"r1": [], "r2": []})

# Step 3: Process each FASTQ path
for path in paths:
    fname = os.path.basename(path)           # Extract just the filename from the full path
    prefix = fname.rsplit("_", 1)[0] + "_"   # Get everything before the last underscore (e.g., "_1" or "_2")

    # Detect and group R1 files
    if fname.endswith("_1.fastq.gz") or fname.endswith("_R1.fastq.gz"):
        paired[prefix]["r1"].append(path)

    # Detect and group R2 files
    elif fname.endswith("_2.fastq.gz") or fname.endswith("_R2.fastq.gz"):
        paired[prefix]["r2"].append(path)

# Step 4: Build a table of valid paired-end samples
records = []
for prefix, files in paired.items():
    # Only keep samples that have both R1 and R2 files
    if files["r1"] and files["r2"]:
        records.append({
            "sample": prefix.rstrip("_"),                # Remove trailing underscore from sample ID
            "r1": ",".join(sorted(files["r1"])),         # Join multiple R1 paths with commas
            "r2": ",".join(sorted(files["r2"]))          # Join multiple R2 paths with commas
        })

# Step 5: Save the collected records as a TSV file
df = pd.DataFrame(records)                              # Convert list of dictionaries to a DataFrame
os.makedirs("config", exist_ok=True)                    # Ensure the output directory exists
df.to_csv("config/samples.tsv", sep="\t", index=False)  # Save as tab-separated values without row index

# Final message
print(f"Saved {len(df)} samples to config/samples.tsv")  # Print number of valid paired samples found
