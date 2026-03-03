import pandas as pd
from pysradb.sraweb import SRAweb  # Correct import for web API (no local DB needed)
import os

# Path to your master list (format: "GSE_DIR SRR" per line)
master_list_path = '/home/chu25/data/rbl/all_rbl_srrs.txt'

# Initialize SRA web API (no DB file required)
sra = SRAweb()

# Read your master list
df_list = pd.read_csv(master_list_path, sep='\s+', header=None, names=['GSE_DIR', 'SRR'])
srrs = df_list['SRR'].unique().tolist()  # Unique SRRs only

# Fetch metadata for each SRR using srr_to_srp (detailed=True for full metadata) and concatenate
dfs = []
for srr in srrs:
    try:
        df = sra.srr_to_srp(srr, detailed=True)
        dfs.append(df)
    except Exception as e:
        print(f"Error fetching {srr}: {e}")
dfs_concat = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# Function to classify based on keywords in title/description
def classify_sc(row):
    text = str(row.get('sample_title', '') + ' ' + str(row.get('description', ''))).lower()
    sc_keywords = ['scrna', 'single-cell', '10x', 'chromium', 'drop-seq']
    if any(kw in text for kw in sc_keywords):
        return 'scRNA-Seq'
    return 'Bulk RNA-Seq'

# Add classification
dfs_concat['Type'] = dfs_concat.apply(classify_sc, axis=1)

# Merge with your list (add GSE_DIR back; map by SRR to 'run_accession')
df_result = df_list.merge(dfs_concat[['run_accession', 'Type', 'sample_title']], left_on='SRR', right_on='run_accession', how='left')

# Output to CSV
output_csv = master_list_path.replace('.txt', '_classified.csv')
df_result.to_csv(output_csv, index=False)
print(f"Classification complete! Results saved to {output_csv}")
print(df_result[['SRR', 'GSE_DIR', 'Type', 'sample_title']].head(10))  # Preview first 10