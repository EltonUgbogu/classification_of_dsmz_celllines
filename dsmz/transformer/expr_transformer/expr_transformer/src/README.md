# Expression Transformer (Masked Gene Modeling)

This directory contains the workflow components for training a **Transformer-based representation model on gene expression data** and exporting learned sample embeddings.

The implementation targets **transcriptomic datasets** (e.g. bulk RNA-seq from cell lines or tumour samples) and supports **self-supervised learning** via masked gene modeling. The resulting embeddings provide a compact representation of each sample and are intended for downstream analyses such as clustering, similarity assessment, and dataset integration.

---

This module applies **Transformer architectures** to gene expression data in order to learn **context-aware representations of samples** and explicitly capture higher-order gene–gene relationships, **without relying on phenotype labels, subtype annotations, or prior biological stratification**.

---

## Masked Gene Modeling (MGM)

Masked Gene Modeling (MGM) is a **self-supervised learning strategy** inspired by masked language modeling, adapted here to gene expression data.

From a biological perspective, MGM leverages the fact that gene expression profiles are not independent measurements, but are structured by **co-regulated pathways, regulatory programs, and cellular state**.

During training:

- Each **sample** is represented as a vector of gene expression values
- A random subset of genes is **masked**
- The model is trained to **reconstruct the masked expression values**
- Learning is driven by the model's ability to infer missing genes from the remaining expression context

In doing so, the model learns gene–gene dependencies and global expression structure that reflect underlying biological organisation rather than predefined labels.

No phenotype information, subtype labels, or external supervision are required.

---

## Scope of this module

This directory provides tools to:

1. Convert expression matrices into a model-ready format
2. Validate dataset integrity and gene order consistency
3. Train a Transformer using masked gene modeling
4. Export learned **sample-level embeddings** for downstream analyses

The focus of this module is **representation learning**, not end-to-end prediction or classification.

---

## Typical workflow

### 1. Expression preprocessing

Convert an RDS expression matrix (samples × genes) into NumPy format with accompanying metadata:

```bash
python dsmz/transformer/expr_transformer/src/prep_rds_to_npz.py \
  --rds <input_matrix.rds> \
  --npz dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.npz \
  --meta dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.meta.json
```

Processing includes per-gene standardisation and storage of reversible preprocessing statistics.

---

### 2. Dataset integrity validation

Before training, gene order and matrix dimensions are validated:

```bash
python dsmz/transformer/expr_transformer/src/check_gene_order.py \
  --npz dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.npz \
  --meta dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.meta.json
```

This step emits a **gene-order hash**, which is used to guarantee compatibility between datasets, trained checkpoints, and exported embeddings.

---

### 3. Model training (GPU recommended)

```bash
python dsmz/transformer/expr_transformer/src/train_mgm.py \
  --npz dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.npz \
  --meta dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.meta.json \
  --ckpt-dir dsmz/transformer/expr_transformer/ckpt \
  --ckpt-prefix exprtf_brca_hvg500 \
  --epochs 30 \
  --batch-size 64
```

Training uses masked gene modeling and stores checkpoints annotated with gene-order metadata to ensure reproducibility.

---

### 4. Embedding export

```bash
python dsmz/transformer/expr_transformer/src/export_embeddings.py \
  --npz dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.npz \
  --meta dsmz/transformer/expr_transformer/data/brca/BRCA_HVG500_joint.meta.json \
  --ckpt dsmz/transformer/expr_transformer/ckpt/<checkpoint>.pt \
  --out dsmz/transformer/expr_transformer/embeddings/brca_hvg500_embeddings.tsv
```

The output is a tab-separated file containing **one embedding vector per sample**, suitable for use in R or Python-based downstream analyses.

---

## Script summary

| Script                 | Purpose                                             |
| ---------------------- | --------------------------------------------------- |
| `prep_rds_to_npz.py`   | Converts RDS matrices to model-ready NumPy format   |
| `check_gene_order.py`  | Validates dimensional consistency and gene identity |
| `train_mgm.py`         | Trains a Transformer via masked gene modeling       |
| `export_embeddings.py` | Exports learned sample embeddings                   |
| `train_exprtf.sh`      | SLURM helper script for cluster execution           |

---

## Computational considerations

* Training is **GPU-recommended**
* Expression matrices remain resident in CPU memory
* Mini-batches are transferred to GPU during training and inference
* Large intermediate artifacts are intentionally excluded from version control

---

## Automation

A Snakemake workflow is provided for reproducible execution:

```bash
snakemake -s dsmz/transformer/expr_transformer/src/Snakefile --use-conda
```

Configuration parameters (datasets, checkpoints, and hyperparameters) are defined in the accompanying `config.yaml`.

---

## Relation to the broader repository

The embeddings generated are used in the
**DSMZ cell line–tumour alignment framework**, enabling unsupervised
comparison, clustering, and representation-level matching across datasets.
