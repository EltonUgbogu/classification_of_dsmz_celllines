import pandas as pd
from collections import defaultdict
import gzip

# Function to parse a GTF file and compute the total exon length per gene
def parse_gtf(gtf_file):
    # Dictionary to store exon intervals (start, end) for each gene_id
    gene_lengths = defaultdict(set)

    # Automatically handle both compressed (.gz) and uncompressed files
    open_func = gzip.open if gtf_file.endswith(".gz") else open

    with open_func(gtf_file, 'rt') as f:
        for line in f:
            # Skip comment lines (starting with '#')
            if line.startswith("#"):
                continue

            # Split the line into its 9 standard GTF columns
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue  # Skip malformed lines

            feature_type = fields[2]
            # Only process lines where the feature type is 'exon'
            if feature_type != "exon":
                continue

            # Get the start and end positions of the exon
            start = int(fields[3])
            end = int(fields[4])
            length = end - start + 1  # Not used here but could be useful

            # Parse the 9th column (attributes), which is a key-value string
            attrs = {}
            for attr in fields[8].split(";"):
                if attr.strip() == "":
                    continue
                key, value = attr.strip().split(" ", 1)
                attrs[key] = value.strip('"')  # Remove quotes around value

            # Extract the gene_id from attributes
            gene_id = attrs.get("gene_id")
            if gene_id:
                # Store the exon interval (start, end) for this gene
                # Using a set ensures each exon interval is unique per gene
                gene_lengths[gene_id].add((start, end))

    # Now compute the total non-overlapping exon length for each gene
    gene_length_sum = {}
    for gene, exons in gene_lengths.items():
        # Merge overlapping exon intervals to avoid double-counting
        merged_exons = merge_intervals(list(exons))
        # Sum the length of merged intervals
        gene_length_sum[gene] = sum(e - s + 1 for s, e in merged_exons)

    # Return the result as a Pandas DataFrame
    return pd.DataFrame.from_dict(gene_length_sum, orient="index", columns=["length"])

# Helper function to merge overlapping intervals
def merge_intervals(intervals):
    # Sort intervals by their start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []

    for interval in sorted_intervals:
        if not merged or interval[0] > merged[-1][1]:
            # No overlap with the previous interval
            merged.append(interval)
        else:
            # Overlap: merge with the last interval
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    
    return merged

# Main script logic when run from the command line
if __name__ == "__main__":
    import sys

    # Use input and output from Snakemake
    gtf_file = snakemake.input.gtf
    output_file = snakemake.output[0]

    # Parse the GTF file to get gene lengths
    gene_lengths_df = parse_gtf(gtf_file)

    # Set 'gene_id' as the index name (for clarity in output)
    gene_lengths_df.index.name = "gene_id"

    # Save the result to a tab-separated file
    gene_lengths_df.to_csv(output_file, sep="\t")
