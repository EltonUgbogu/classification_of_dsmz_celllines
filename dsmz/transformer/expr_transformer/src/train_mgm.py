#!/usr/bin/env python3
"""
train_mgm.py
============

This script trains a Masked Gene Model (MGM) on gene expression data using a
Transformer architecture - similar to how BERT learns language by predicting
masked words, but for genes instead.

**The Core Idea: Masked Gene Modeling**

Just as BERT learns to understand language by predicting [MASK] tokens:
    "The cat sat on the [MASK]" → predict "mat"

We learn gene expression patterns by predicting masked gene values:
    [TP53=2.1, BRCA1=[MASK], MYC=1.5, ...] → predict BRCA1's expression

The model learns gene-gene relationships because to predict a masked gene's
expression, it must understand how that gene correlates with visible genes.

**Why This Works for Biology:**

1. Genes don't act in isolation - they form regulatory networks
2. If TP53 is highly expressed, certain downstream genes will be too
3. The model learns these co-expression patterns implicitly
4. The learned representations capture biological signal (pathways, subtypes)

**Training Pipeline:**

1. Load z-scored expression matrix (samples × genes)
2. For each sample, randomly mask 15% of genes (set to 0)
3. Model predicts original values for masked genes only
4. Backpropagate loss, update weights
5. Save checkpoint after each epoch

**Output:**

Checkpoint files containing trained model weights that can be used for:
- Generating sample embeddings (for clustering, UMAP, etc.)
- Transfer learning to other expression datasets
- Predicting missing gene values

"""

import json
import hashlib
import os
import math
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm  # Progress bar library


# =============================================================================
# FILE PATHS
# =============================================================================

# Input: Preprocessed and z-scored expression matrix
NPZ_PATH = os.path.abspath("data/BRCA_HVG500_joint.npz")

# Input: Metadata (sample IDs, gene names) - not used during training but
# included for completeness
META_PATH = os.path.abspath("data/BRCA_HVG500_joint.meta.json")

# Metadata keys commonly used for gene lists
GENE_KEY_CANDIDATES = ("genes", "gene_names", "features", "gene_ids")

# Output: Directory for model checkpoints
CKPT_DIR = os.path.abspath("ckpt")
os.makedirs(CKPT_DIR, exist_ok=True)


def load_gene_list(meta):
    """Return the gene list and the metadata key used to store it."""

    gene_key = next((k for k in GENE_KEY_CANDIDATES if k in meta), None)
    if gene_key is None:
        raise KeyError(
            "Could not find gene list in meta.json. "
            f"Tried: {list(GENE_KEY_CANDIDATES)}. Keys: {list(meta.keys())}"
        )

    genes = meta[gene_key]
    if not isinstance(genes, list):
        raise TypeError(f"meta['{gene_key}'] must be a list, got {type(genes)!r}")

    return gene_key, genes


def sha256_lines(values):
    joined = "\n".join(values).encode("utf-8")
    return hashlib.sha256(joined).hexdigest()


# =============================================================================
# DATASET CLASS
# =============================================================================

class ExprDataset(Dataset):
    """
    PyTorch Dataset for masked gene modeling.
    
    This class handles the data loading and masking strategy. Each time we
    access a sample, it randomly masks a portion of genes and returns:
    - The masked input (what the model sees)
    - The original values (ground truth for loss calculation)
    - The mask itself (which positions were masked)
    
    **Why random masking each access?**
    
    Unlike a fixed train/test split, we randomly mask different genes each
    time we see a sample. This provides data augmentation - the same sample
    contributes to learning different gene relationships across epochs.
    
    Parameters
    ----------
    X : np.ndarray
        Expression matrix of shape (n_samples, n_genes), z-scored.
    mask_prob : float, default=0.15
        Fraction of genes to mask. 15% is a common choice from BERT.
        - Too low (5%): Model doesn't learn enough from each sample
        - Too high (50%): Not enough context to make good predictions
    """
    
    def __init__(self, X, mask_prob=0.15):
        # Convert to PyTorch tensor once, keep in memory
        # Shape: (N, G) where N=samples, G=genes
        self.X = torch.from_numpy(X)
        self.mask_prob = mask_prob
    
    def __len__(self):
        """Return number of samples in dataset."""
        return self.X.shape[0]
    
    def __getitem__(self, i):
        """
        Get a single sample with random masking applied.
        
        Returns
        -------
        x : torch.Tensor
            Masked input, shape (G,). Masked positions are set to 0.
        target : torch.Tensor
            Original values, shape (G,). Used for computing loss.
        mask : torch.Tensor
            Boolean mask, shape (G,). True where genes were masked.
        mask_f : torch.Tensor
            Float mask indicator (1.0 for masked positions) used for adding a
            learned mask embedding inside the model.
        """
        # Clone to avoid modifying the original data
        x = self.X[i].clone()  # Shape: (G,)
        target = x.clone()     # Keep original for loss computation
        
        # Generate random mask: each gene has mask_prob chance of being masked
        # torch.rand_like creates uniform random values in [0, 1)
        # Comparing to mask_prob gives us ~15% True values
        mask = torch.rand_like(x) < self.mask_prob
        
        # ---------------------------------------------------------------------
        # MASKING STRATEGY
        # ---------------------------------------------------------------------
        # Set masked positions to 0.0 (the mean in z-scored space)
        #
        # Alternative strategies (from BERT paper):
        # - 80% of time: replace with 0 (our [MASK] token equivalent)
        # - 10% of time: replace with random gene's value
        # - 10% of time: keep original (forces model to be robust)
        #
        # For simplicity, we always use 0. This works well because:
        # - In z-scored space, 0 = mean expression (uninformative)
        # - The model must infer the true value from context
        x[mask] = 0.0

        # Float mask (1.0 where masked) used to add a learned indicator vector
        mask_f = mask.float()

        return x, target, mask, mask_f


# =============================================================================
# MODEL DEFINITION
# =============================================================================

class ExprTransformer(nn.Module):
    """
    Expression Transformer for Masked Gene Modeling.
    
    This architecture treats each gene as a "token" and uses self-attention
    to learn relationships between genes. The key insight is that predicting
    masked genes requires understanding the gene regulatory network.
    
    **Architecture Diagram:**
    
    Input: [gene1_val, gene2_val, ..., geneG_val]  (G scalars)
           ↓
    Value Projection: each scalar → d_model vector
           ↓
    + Gene Embeddings: add learned gene identity vectors
           ↓
    Transformer Encoder (6 layers of self-attention)
           ↓
    Output Projection: each d_model vector → 1 scalar
           ↓
    Output: [pred1, pred2, ..., predG]  (G predicted values)
    
    **What each component learns:**
    
    - gene_emb: "Gene TP53 has these properties" (gene identity)
    - val_proj: "Expression value 2.5 means this" (value semantics)
    - encoder: "When TP53 is high, BRCA1 tends to be..." (gene relationships)
    - out: "Based on context, this gene should express at level X"
    
    Parameters
    ----------
    n_genes : int
        Number of genes (vocabulary size)
    d_model : int, default=256
        Hidden dimension throughout the model
    n_layers : int, default=6
        Number of Transformer encoder layers
    n_heads : int, default=8
        Number of attention heads per layer
    dropout : float, default=0.1
        Dropout rate for regularization
    """
    
    def __init__(self, n_genes, d_model=256, n_layers=6, n_heads=8, dropout=0.1):
        super().__init__()
        self.n_genes = n_genes
        self.d_model = d_model
        
        # ---------------------------------------------------------------------
        # GENE IDENTITY EMBEDDINGS
        # ---------------------------------------------------------------------
        # Each gene gets a unique learnable vector that represents its identity.
        # This is analogous to word embeddings in NLP.
        #
        # The model learns things like:
        # - "TP53 is a tumor suppressor involved in apoptosis"
        # - "MYC is an oncogene involved in proliferation"
        # (Encoded implicitly in the embedding vectors)
        self.gene_emb = nn.Embedding(n_genes, d_model)
        
        # ---------------------------------------------------------------------
        # EXPRESSION VALUE PROJECTION
        # ---------------------------------------------------------------------
        # Project scalar expression values to d_model dimensions.
        # This allows the model to understand "what does expression=2.5 mean?"
        #
        # Input: (B, G, 1) - one scalar per gene
        # Output: (B, G, d_model) - rich representation of expression level
        self.val_proj = nn.Linear(1, d_model)

        # ---------------------------------------------------------------------
        # MASK INDICATOR EMBEDDING
        # ---------------------------------------------------------------------
        # Learned vector added at positions that were masked during training.
        # Initialized to zeros so early training behaves like the previous
        # version before the model learns to leverage the indicator.
        self.mask_emb = nn.Parameter(torch.zeros(1, 1, d_model))

        # ---------------------------------------------------------------------
        # GLOBAL [CLS] TOKEN
        # ---------------------------------------------------------------------
        # Learnable token prepended to every sequence so the encoder is trained
        # with the same structure used at embedding export time.
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model))
        
        # ---------------------------------------------------------------------
        # TRANSFORMER ENCODER
        # ---------------------------------------------------------------------
        # Stack of self-attention layers that learn gene-gene relationships.
        #
        # Key hyperparameters:
        # - d_model=256: Hidden size (affects capacity and memory)
        # - nhead=8: Parallel attention heads (different relationship types)
        # - dim_feedforward=1024: MLP hidden size (4×d_model is standard)
        # - activation="gelu": Smooth activation (better than ReLU for Transformers)
        # - batch_first=True: Input shape (B, G, d) not (G, B, d)
        # - norm_first=True: Pre-LayerNorm (more stable training)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=4 * d_model,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True
        )
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=n_layers)
        
        # ---------------------------------------------------------------------
        # OUTPUT PROJECTION
        # ---------------------------------------------------------------------
        # Project from d_model back to scalar predictions.
        # One prediction per gene position.
        self.out = nn.Linear(d_model, 1)
        
        # ---------------------------------------------------------------------
        # GENE INDEX BUFFER
        # ---------------------------------------------------------------------
        # Pre-computed tensor [0, 1, 2, ..., n_genes-1] for embedding lookup.
        # register_buffer: moves with model to GPU, but not a learnable parameter
        # persistent=False: don't save to checkpoint (trivial to recreate)
        self.register_buffer(
            "gene_ids",
            torch.arange(n_genes, dtype=torch.long),
            persistent=False
        )
    
    def forward(self, x, mask_f=None):
        """
        Forward pass for training (predicts all gene values).

        Parameters
        ----------
        x : torch.Tensor
            Masked expression values, shape (B, G).
            Masked positions contain 0.0.
        mask_f : torch.Tensor, optional
            Float mask indicator of shape (B, G) where 1.0 marks masked
            positions. When provided, a learned mask embedding is added at
            those locations.
        
        Returns
        -------
        y : torch.Tensor
            Predicted expression values, shape (B, G).
            Predictions for ALL genes (loss computed only on masked ones).
        """
        B, G = x.shape
        
        # Expand gene indices for batch: (G,) → (B, G)
        gene_ids = self.gene_ids.unsqueeze(0).expand(B, G)
        
        # Project expression values: (B, G) → (B, G, 1) → (B, G, d_model)
        v = self.val_proj(x.unsqueeze(-1))
        
        # Look up gene embeddings: (B, G) → (B, G, d_model)
        g = self.gene_emb(gene_ids)
        
        # Combine value information with gene identity
        # This is the input representation for the Transformer
        h = v + g  # Shape: (B, G, d_model)

        if mask_f is not None:
            # mask_f: (B, G) -> (B, G, 1)
            h = h + self.mask_emb * mask_f.unsqueeze(-1)

        # Prepend [CLS] token so the encoder trains with the same context used
        # during embedding export.
        cls = self.cls_token.expand(B, -1, -1)  # (B, 1, d_model)
        h = torch.cat([cls, h], dim=1)  # (B, G + 1, d_model)

        # Pass through Transformer encoder
        # Self-attention allows each gene to aggregate information from all others
        h = self.encoder(h)  # Shape: (B, G + 1, d_model)

        # Drop the CLS token for per-gene predictions
        h_genes = h[:, 1:, :]  # (B, G, d_model)

        # Project to predictions: (B, G, d_model) → (B, G, 1) → (B, G)
        y = self.out(h_genes).squeeze(-1)
        
        return y
    
    @torch.no_grad()
    def embed(self, x):
        """
        Extract sample embeddings (for downstream tasks, not training).
        
        This method is used AFTER training to generate fixed-size
        representations of samples for clustering, visualization, etc.
        
        Parameters
        ----------
        x : torch.Tensor
            Full expression values (no masking), shape (B, G).
        
        Returns
        -------
        embeddings : torch.Tensor
            Sample embeddings, shape (B, 2 * d_model) when default d_model=256.
            Concatenates [CLS] token with mean-pooled gene tokens for a hybrid
            global representation.
        """
        B, G = x.shape
        gene_ids = self.gene_ids.unsqueeze(0).expand(B, G)
        
        v = self.val_proj(x.unsqueeze(-1))
        g = self.gene_emb(gene_ids)
        h = v + g  # (B, G, d_model)

        # Prepend [CLS] for consistency with training/inference forward pass
        cls = self.cls_token.expand(B, -1, -1)  # (B, 1, d_model)
        h = torch.cat([cls, h], dim=1)  # (B, G + 1, d_model)

        # Encode sequence
        h = self.encoder(h)  # (B, G + 1, d_model)

        # Hybrid embedding: CLS + mean over gene tokens
        cls_emb = h[:, 0, :]  # (B, d_model)
        mean_emb = h[:, 1:, :].mean(dim=1)  # (B, d_model)
        return torch.cat([cls_emb, mean_emb], dim=1)  # (B, 2 * d_model)


# =============================================================================
# TRAINING LOOP
# =============================================================================

def main():
    """
    Main training function.
    
    This implements a standard PyTorch training loop with:
    - Mixed precision training (faster on modern GPUs)
    - Gradient scaling (prevents underflow with FP16)
    - Progress bars and logging
    - Checkpoint saving after each epoch
    """
    
    # -------------------------------------------------------------------------
    # SETUP
    # -------------------------------------------------------------------------
    
    if not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA required. Submit this job to a GPU node (e.g. --partition=gpu --gres=gpu:1)."
        )
    device = "cuda"
    print("Training on: cuda:", torch.cuda.get_device_name(0))
    
    # Load preprocessed data
    data = np.load(NPZ_PATH)
    if "X" not in data:
        raise KeyError(f"NPZ file does not contain 'X'. Keys: {list(data.keys())}")

    X = data["X"]  # Shape: (N, G), dtype: float32, z-scored

    N, G = X.shape
    print(f"Data: {N} samples × {G} genes")

    # Validate metadata alignment with NPZ
    with open(META_PATH) as f:
        meta = json.load(f)

    if "samples" not in meta:
        raise KeyError(f"meta.json has no 'samples' key. Keys: {list(meta.keys())}")
    samples = meta["samples"]
    if not isinstance(samples, list):
        raise TypeError(f"meta['samples'] must be a list, got {type(samples)!r}")

    gene_key, genes = load_gene_list(meta)
    if len(genes) != G:
        raise ValueError(
            f"Gene count mismatch: NPZ has {G} genes but meta['{gene_key}'] has {len(genes)}"
        )

    if len(samples) != N:
        raise ValueError(
            f"Sample count mismatch: NPZ has {N} rows but meta['samples'] has {len(samples)}"
        )

    n_unique = len(set(genes))
    if n_unique != len(genes):
        raise ValueError(
            f"Duplicate genes detected: unique={n_unique}, total={len(genes)}"
        )

    if not all(isinstance(g, str) and g for g in genes):
        raise TypeError("Some gene IDs are empty or not strings")

    gene_hash = sha256_lines(genes)
    print(
        "Gene list validated: "
        f"first={genes[:5]} last={genes[-5:] if genes else []} sha256={gene_hash}"
    )
    
    # Create dataset with 15% masking probability
    ds = ExprDataset(X, mask_prob=0.15)
    
    # -------------------------------------------------------------------------
    # DATALOADER
    # -------------------------------------------------------------------------
    # batch_size: Number of samples per gradient update
    #   - Larger = more stable gradients, faster training, but more memory
    #   - GPU can handle larger batches than CPU
    #
    # shuffle=True: Randomize sample order each epoch (important for training)
    #
    # num_workers=2: Background processes for data loading (overlap with GPU)
    #
    # pin_memory=True: Use pinned (page-locked) memory for faster CPU→GPU transfer
    batch_size = 256
    dl = DataLoader(
        ds,
        batch_size=batch_size,
        shuffle=True,
        num_workers=2,
        pin_memory=True
    )
    
    # -------------------------------------------------------------------------
    # MODEL
    # -------------------------------------------------------------------------
    model = ExprTransformer(
        n_genes=G,
        d_model=256,   # Hidden dimension
        n_layers=6,    # Transformer depth
        n_heads=8,     # Attention heads
        dropout=0.1    # Regularization
    ).to(device)
    
    # Count parameters
    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {n_params:,}")
    
    # -------------------------------------------------------------------------
    # OPTIMIZER
    # -------------------------------------------------------------------------
    # AdamW: Adam with proper weight decay (L2 regularization)
    #   - lr=3e-4: Learning rate (how big each update step is)
    #   - weight_decay=0.01: Regularization to prevent overfitting
    opt = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=0.01)
    
    # -------------------------------------------------------------------------
    # MIXED PRECISION TRAINING
    # -------------------------------------------------------------------------
    # GradScaler helps with mixed precision training (FP16 + FP32):
    #   - Forward pass in FP16 (faster, less memory)
    #   - Gradients scaled up to prevent underflow
    #   - Optimizer step in FP32 (full precision for stability)
    #   - Only enabled on CUDA (CPU doesn't benefit)
    scaler = torch.amp.GradScaler("cuda", enabled=(device == "cuda"))
    
    # -------------------------------------------------------------------------
    # LOSS FUNCTION
    # -------------------------------------------------------------------------
    # MSE (Mean Squared Error) for regression
    # reduction="sum": Return sum of squared errors (we'll normalize manually)
    # This allows us to correctly weight by number of masked genes
    mse = nn.MSELoss(reduction="sum")
    
    # -------------------------------------------------------------------------
    # TRAINING HYPERPARAMETERS
    # -------------------------------------------------------------------------
    epochs = 30      # Number of passes through entire dataset
    log_every = 100  # Log progress every N batches
    global_step = 0  # Total batches processed (across all epochs)
    
    print(f"\nStarting training for {epochs} epochs...")
    print(f"Batch size: {batch_size}, Batches per epoch: {len(dl)}")
    
    # -------------------------------------------------------------------------
    # TRAINING LOOP
    # -------------------------------------------------------------------------
    for ep in range(1, epochs + 1):
        model.train()  # Enable dropout, etc.
        
        # Accumulators for computing average loss
        running = 0.0  # Sum of losses
        denom = 0.0    # Count of masked genes (for averaging)
        
        # Progress bar for this epoch
        pbar = tqdm(dl, desc=f"epoch {ep}/{epochs}", leave=True)
        
        for x, target, mask, mask_f in pbar:
            # -----------------------------------------------------------------
            # DATA TRANSFER
            # -----------------------------------------------------------------
            # Move data to GPU (non_blocking=True allows overlap with computation)
            x = x.to(device, non_blocking=True)
            target = target.to(device, non_blocking=True)
            mask = mask.to(device, non_blocking=True)
            mask_f = mask_f.to(device, non_blocking=True)
            
            # -----------------------------------------------------------------
            # FORWARD PASS
            # -----------------------------------------------------------------
            # Zero gradients (set_to_none=True is slightly faster than zero_grad())
            opt.zero_grad(set_to_none=True)
            
            # Automatic mixed precision context
            # autocast: automatically use FP16 where safe, FP32 where needed
            with torch.amp.autocast("cuda", enabled=(device == "cuda")):
                # Get predictions for all genes
                pred = model(x, mask_f=mask_f)  # Shape: (B, G)
                
                # Compute MSE only on masked positions
                # diff[mask] selects only the masked predictions
                diff = pred - target
                
                # MSE loss: sum of squared differences for masked genes
                # We compare to zeros because we want (pred - target)² = diff²
                loss = mse(diff[mask], torch.zeros_like(diff[mask]))
                
                # Normalize by number of masked entries
                # This keeps the loss scale consistent regardless of mask count
                # clamp(min=1) prevents division by zero edge case
                mcount = mask.sum().clamp(min=1)
                loss = loss / mcount
            
            # -----------------------------------------------------------------
            # BACKWARD PASS
            # -----------------------------------------------------------------
            # scale: Multiply loss by scale factor before backward
            # This prevents gradient underflow in FP16
            scaler.scale(loss).backward()
            
            # step: Unscale gradients, then call optimizer.step()
            scaler.step(opt)
            
            # update: Adjust scale factor for next iteration
            # (Increases if no overflow, decreases if overflow detected)
            scaler.update()
            
            # -----------------------------------------------------------------
            # LOGGING
            # -----------------------------------------------------------------
            # Accumulate loss (weighted by number of masked genes)
            running += loss.item() * float(mcount)
            denom += float(mcount)
            
            global_step += 1
            if global_step % log_every == 0:
                avg_loss = running / max(denom, 1.0)
                pbar.set_postfix({"masked_mse": f"{avg_loss:.4f}"})
        
        # ---------------------------------------------------------------------
        # EPOCH COMPLETE: SAVE CHECKPOINT
        # ---------------------------------------------------------------------
        avg_loss = running / max(denom, 1.0)
        
        # Checkpoint contains everything needed to resume or use the model:
        # - epoch: Training progress
        # - model_state: Learned weights
        # - G: Number of genes (needed to instantiate model)
        # - config: Architecture hyperparameters
        ckpt_path = os.path.join(CKPT_DIR, f"exprtf_brca_hvg500_ep{ep:02d}.pt")
        torch.save({
            "epoch": ep,
            "model_state": model.state_dict(),
            "G": G,
            "gene_hash": gene_hash,
            "gene_key": gene_key,
            "gene_examples": {
                "first": genes[:5],
                "last": genes[-5:]
            },
            "config": {
                "d_model": 256,
                "n_layers": 6,
                "n_heads": 8,
                "dropout": 0.1
            }
        }, ckpt_path)
        
        print(f"[OK] epoch={ep} masked_mse={avg_loss:.6f} saved={ckpt_path}")
    
    print("\nTraining complete!")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    main()

