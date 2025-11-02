#!/usr/bin/env python3
import pandas as pd

# -----------------------------------------------
# Step 1: Read Input File Paths from Snakemake
# -----------------------------------------------
# These variables are passed in from the Snakemake workflow
counts_file = snakemake.input.counts        # Input: STAR gene count file (TSV with 4 columns: gene_id, unstranded, forward, reverse)
lengths_file = snakemake.input.lengths      # Input: gene lengths (TSV with columns: gene, length)
out_file = snakemake.output.tpm             # Output: calculated TPM file (TSV with columns: gene, TPM)

# -----------------------------------------------
# Step 2: Load Raw Gene Counts
# -----------------------------------------------
# STAR's ReadsPerGene.out.tab contains:
# Column 1: Gene ID
# Column 2: Unstranded read count
# Column 3: Forward-stranded read count
# Column 4: Reverse-stranded read count
# We only use the unstranded counts here (column 2), for general compatibility
# Lines starting with "_" or "__" are STAR summary lines — skipped using `comment='_'`
counts = pd.read_csv(
    counts_file,
    sep='\t',
    header=None,        # No header row in the file
    comment='_',        # Skip STAR metadata rows
    usecols=[0, 1],     # Only read gene_id and unstranded counts
    names=["gene", "count"]  # Rename columns for clarity
)

# -----------------------------------------------
# Step 3: Load Gene Lengths
# -----------------------------------------------
# Format: TSV with columns "gene" and "length" (length in base pairs)
lengths = pd.read_csv(lengths_file, sep='\t')

# -----------------------------------------------
# Step 4: Merge Counts and Lengths
# -----------------------------------------------
# We join on the "gene" column (inner join by default).
# This ensures that only genes present in both files are retained.
df = counts.merge(lengths, on='gene')

# -----------------------------------------------
# Step 5: Filter Out Invalid Entries
# -----------------------------------------------
# TPM calculation requires non-zero lengths (to avoid division by zero)
# Drop genes with length = 0 or missing
df = df[df['length'] > 0]

# -----------------------------------------------
# Step 6: Compute Normalized Read Count (CPK)
# -----------------------------------------------
# Mathematical explanation:
#   CPK = Count per kilobase
#   = raw_count / (length_in_bp / 1000)
#   = (C_g × 10^3) / L_g
#   where:
#     - C_g = raw count for gene g
#     - L_g = gene length in base pairs
#
# Why normalize by length?
#   Longer genes tend to have more reads just by chance.
#   Normalizing adjusts for gene size.
df['norm'] = df['count'] / (df['length'] / 1000)

# -----------------------------------------------
# Step 7: Calculate TPM (Transcripts Per Million)
# -----------------------------------------------
# TPM_g = (CPK_g / ΣCPK_all_genes) × 1e6
#
# This scales the normalized expression so the total expression per sample = 1,000,000
# Benefits of TPM:
#   - Comparable across samples
#   - Reflects relative abundance of transcripts
df['TPM'] = (df['norm'] / df['norm'].sum()) * 1e6

# -----------------------------------------------
# Step 8: Output Final TPM Table
# -----------------------------------------------
# Only keep gene ID and TPM columns
# Write as tab-separated file
df[['gene', 'TPM']].to_csv(out_file, sep='\t', index=False)
