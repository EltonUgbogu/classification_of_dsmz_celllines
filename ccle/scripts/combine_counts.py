import sys
import os
import pandas as pd

print("Starting combine_counts.py...")

# Extract input files and output file
files = sys.argv[1:-1]
outfile = sys.argv[-1]

print("Files:",files)
print(f"Number of input count files: {len(files)}")
print(f"Output will be written to: {outfile}")

dfs = []

for f in files:
    print(f"Reading file: {f}")
    try:
        df = pd.read_csv(f, comment='#', sep='\t')
        print(f"Successfully read {f} with shape {df.shape}")

        # Extract the last column (read count) and Geneid
        count_col = df.columns[-1]
        print(f"Using columns: ['Geneid', '{count_col}']")

        df = df[["Geneid", count_col]]

        # Rename count column to sample name
        sample_name = os.path.basename(f).replace(".counts.txt", "")
        df.columns = ["Geneid", sample_name]
        print(f"Renamed count column to: {sample_name}")

        dfs.append(df.set_index("Geneid"))
        print(f"Appended data from {f}")

    except Exception as e:
        print(f"Error processing {f}: {e}", file=sys.stderr)

if not dfs:
    raise ValueError("No valid count files were processed. Exiting.")

print("Concatenating all count tables...")
combined = pd.concat(dfs, axis=1)
print(f"Combined matrix shape: {combined.shape}")

print(f"Writing merged counts to {outfile}")
combined.to_csv(outfile, sep='\t')

print("Done.")
