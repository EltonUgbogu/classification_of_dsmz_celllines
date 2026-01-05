
+This folder bundles the command-line utilities used to prepare expression matrices, train the masked gene modeling (MGM) Transformer, validate data integrity, and export embeddings. Run them from the repository root (`classification_of_dsmz_celllines`) so relative paths resolve correctly.
+
+## File map
+- **`prep_rds_to_npz.py`** – Reads an RDS matrix (samples × genes), z-scores per gene, and writes a compressed NPZ plus JSON metadata (sample IDs and gene list). Retains per-gene mean/SD for reversing the transform and tolerates zero-variance genes.
+- **`check_gene_order.py`** – Confirms that an NPZ + metadata pair has consistent dimensions, enforces unique gene IDs, and prints a SHA256 fingerprint of the gene order for checkpoint provenance. Accepts `--npz` and `--meta` overrides.
+- **`train_mgm.py`** – Trains the Expression Transformer for masked gene modeling: dataset masking, model construction (`ExprTransformer`), training/validation loop, checkpointing with gene-order hashes, and GPU-only enforcement. Defaults point to `data/BRCA_HVG500_joint.npz` and the matching metadata JSON, but you can override `--npz`, `--meta`, `--ckpt-dir`, `--ckpt-prefix`, `--epochs`, `--batch-size`, and `--log-every`.
+- **`train_exprtf.sh`** – Slurm wrapper that activates the Conda env, runs training, then exports embeddings. Adjust partition/`--gres` to your cluster.
+- **`export_embeddings.py`** – Loads a trained checkpoint, validates gene order, and writes CLS + mean-pooled embeddings to TSV. Keeps the full matrix on CPU and moves batches to GPU; execution is GPU-only.
+
+## Typical workflow
+1. Convert RDS expression data: `python dsmz/transformer/expr_transformer/src/prep_rds_to_npz.py --rds <input.rds> --npz data/BRCA_HVG500_joint.npz --meta data/BRCA_HVG500_joint.meta.json`.
+2. Verify integrity and capture the gene hash: `python dsmz/transformer/expr_transformer/src/check_gene_order.py --npz data/BRCA_HVG500_joint.npz --meta data/BRCA_HVG500_joint.meta.json`.
+3. Train the MGM Transformer: `python dsmz/transformer/expr_transformer/src/train_mgm.py --npz data/BRCA_HVG500_joint.npz --meta data/BRCA_HVG500_joint.meta.json --ckpt-prefix exprtf_brca_hvg500 --epochs 30`. On Slurm, prefer `bash dsmz/transformer/expr_transformer/src/train_exprtf.sh` after tweaking resource flags.
+4. Export embeddings: `python dsmz/transformer/expr_transformer/src/export_embeddings.py --npz data/BRCA_HVG500_joint.npz --meta data/BRCA_HVG500_joint.meta.json --ckpt checkpoints/latest.pt --out embeddings.tsv`.
+
+## GPU expectations
+Training and embedding export are GPU-only and will fail fast if CUDA is unavailable. Batches are moved to GPU on-demand while the full dataset stays on CPU to avoid exhausting device memory.
+
+## Snakemake automation
+The repository now includes a GPU-only Snakemake workflow to standardize validation, training, and embedding export across datasets:
+
+- **Config:** `config.yaml` defines datasets, checkpoint prefixes, default hyperparameters, and output locations.
+- **Workflow:** `Snakefile` orchestrates `check_gene_order.py`, `train_mgm.py`, and `export_embeddings.py` with logs/benchmarks under `logs/` and `benchmarks/`.
+- **Conda env:** `workflow/envs/exprtf.yml` (PyTorch + data deps) is used for all rules.
+- **SLURM profile:** `profiles/slurm/` requests GPUs where required and ships a submit helper.
+
+Run locally on a GPU node:
+
+```bash
+snakemake -j 4 --use-conda
+```
+
+Run on SLURM with GPU requests handled automatically:
+
+```bash
+snakemake --profile profiles/slurm
+```
+
+If you prefer to point Snakemake directly at the workflow that lives next to
+these scripts, use the bundled path explicitly (this works from any directory):
+
+```bash
+snakemake -s dsmz/transformer/expr_transformer/src/Snakefile --use-conda -j 4
+```
+
+Outputs land in `ckpt/` and `embeddings/` using the dataset-specific prefixes from `config.yaml` (e.g., `ckpt/exprtf_brca_hvg500_ep30.pt`, `embeddings/brca_hvg500_embeddings.tsv`).
+
+## Path hints
+- Defaults assume you run from the repo root. Override paths with CLI flags as needed.
+- `prep_rds_to_npz.py` auto-fills `--rds` when `$USER` is set (e.g., `/work/<user>/data/BRCA/...`).
