import pandas as pd
import sys

# Inputs
count_files = snakemake.input.counts
lengths_file = snakemake.input.lengths
output_file = snakemake.output[0]

# Load gene lengths
gene_lengths = pd.read_csv(lengths_file, sep="\t", index_col=0, names=["gene_id", "length"])

# TPM calculation function
def compute_tpm(counts, gene_lengths):
    rpk = counts.div(gene_lengths["length"], axis=0)
    scaling_factors = rpk.sum(axis=0) / 1e6
    return rpk.div(scaling_factors, axis=1)

# Load and merge count tables
all_counts = []
sample_names = []

for f in count_files:
    df = pd.read_csv(f, sep="\t", comment="#", index_col=0)
    df = df[[df.columns[5]]]  # assumes 6th column is raw count
    sample_name = f.split("/")[-1].replace(".tsv", "")
    df.columns = [sample_name]
    all_counts.append(df)
    sample_names.append(sample_name)

# Merge all samples into one dataframe
counts_df = pd.concat(all_counts, axis=1)

# Align with gene lengths
counts_df = counts_df.loc[counts_df.index.isin(gene_lengths.index)]
gene_lengths = gene_lengths.loc[counts_df.index]

# Calculate TPM
tpm_df = compute_tpm(counts_df, gene_lengths)

# Save
tpm_df.to_csv(output_file, sep="\t")
