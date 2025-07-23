import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import numpy as np
import time
from datetime import datetime

start_time = time.time()

with open(snakemake.log[0], "w") as log:
    try:
        log.write(f"[{datetime.now()}]START: plot_pca_heatmaps.py\n")
        log.write(f"Input file: {snakemake.input[0]}\n")
        log.write(f"Output files: {snakemake.output[0]}, {snakemake.output[1]}\n\n")

        log.write(f"[{datetime.now()}] Loading count matrix...\n")
        df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

        log.write(f"[{datetime.now()}] Log-transforming counts...\n")
        df_log = np.log1p(df)

        log.write(f"[{datetime.now()}] Selecting top 100 variable genes...\n")
        top_genes = df_log.var(axis=1).sort_values(ascending=False).head(100).index
        df_top = df_log.loc[top_genes]

        log.write(f"[{datetime.now()}] Plotting heatmap...\n")
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_top, cmap="viridis", xticklabels=True, yticklabels=False)
        plt.title("Top 100 Variable Genes (log-transformed)")
        plt.tight_layout()
        plt.savefig(snakemake.output[0], format="pdf")
        log.write(f"[{datetime.now()}] Heatmap saved to {snakemake.output[0]}\n")

        log.write(f"[{datetime.now()}] Running PCA on samples...\n")
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df_log.T)
        pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=df.columns)

        log.write(f"[{datetime.now()}] Plotting PCA scatterplot...\n")
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x="PC1", y="PC2", data=pca_df, s=100)
        for label in pca_df.index:
            plt.text(pca_df.loc[label, "PC1"] + 0.2, pca_df.loc[label, "PC2"], label)
        plt.title("PCA of Samples")
        plt.tight_layout()
        plt.savefig(snakemake.output[1], format="pdf")
        log.write(f"[{datetime.now()}] PCA plot saved to {snakemake.output[1]}\n")

        duration = round(time.time() - start_time, 2)
        log.write(f"\n[{datetime.now()}] FINISHED in {duration} seconds\n")

    except Exception as e:
        log.write(f"\n[{datetime.now()}] ERROR: {str(e)}\n")
        raise
