import pandas as pd
import matplotlib.pyplot as plt
import time
from datetime import datetime

start_time = time.time()

with open(snakemake.log[0], "w") as log:
    try:
        log.write(f"[{datetime.now()}] START: plot_read_counts.py\n")
        log.write(f"Input file: {snakemake.input[0]}\n")
        log.write(f"Output file: {snakemake.output[0]}\n\n")

        log.write(f"[{datetime.now()}] Reading count matrix...\n")
        df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

        log.write(f"[{datetime.now()}] Summing read counts per sample...\n")
        sample_totals = df.sum(axis=0)

        log.write(f"[{datetime.now()}] Generating bar plot...\n")
        plt.figure(figsize=(10, 6))
        sample_totals.plot(kind='bar', color='skyblue')
        plt.title('Total Read Counts per Sample')
        plt.ylabel('Total Counts')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(snakemake.output[0], format="pdf")
        log.write(f"[{datetime.now()}] Plot saved to: {snakemake.output[0]}\n")

        duration = round(time.time() - start_time, 2)
        log.write(f"\n[{datetime.now()}] FINISHED in {duration} seconds\n")

    except Exception as e:
        log.write(f"\n[{datetime.now()}] ERROR: {str(e)}\n")
        raise
