#!/usr/bin/env python3
"""
export_embeddings.py
====================

This script loads a pretrained Expression Transformer model and uses it to
generate embeddings (learned representations) from gene expression data.

**What are embeddings?**
Embeddings are dense, lower-dimensional vector representations that capture
meaningful patterns in high-dimensional data. Here, we transform ~500 gene
expression values per sample into a 256-dimensional vector that captures
the "essence" of that sample's expression profile.

**Workflow:**
1. Load preprocessed gene expression data (NPZ format)
2. Load a pretrained Transformer model checkpoint
3. Pass all samples through the model to generate embeddings
4. Save embeddings as a TSV file for downstream analysis (clustering, UMAP, etc.)

"""

import os
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn


# =============================================================================
# FILE PATHS - Configure these for your environment
# =============================================================================

# Input: Preprocessed expression matrix (samples × genes) saved as NumPy archive
NPZ_PATH = os.path.abspath("data/BRCA_HVG500_joint.npz")

# Input: Metadata JSON containing sample IDs and gene names
META_PATH = os.path.abspath("data/BRCA_HVG500_joint.meta.json")

# Input: Trained model checkpoint containing weights and configuration
CKPT_PATH = os.path.abspath("ckpt/exprtf_brca_hvg500_ep30.pt")

# Output: TSV file where each row is a sample, columns are embedding dimensions
OUT_TSV = os.path.abspath("embeddings/BRCA_HVG500_joint_exprtf.tsv")


# =============================================================================
# MODEL DEFINITION
# =============================================================================

class ExprTransformer(nn.Module):
    """
    Expression Transformer: A Transformer-based model for gene expression data.
    
    **Architecture Overview:**
    
    The model treats each gene as a "token" (similar to words in NLP) and learns
    to understand relationships between genes using self-attention.
    
    For each sample:
        Input:  [gene1_expr, gene2_expr, ..., geneG_expr]  (G values)
        Output: [z1, z2, ..., z256]                        (d_model dimensions)
    
    **Key Components:**
    
    1. Gene Embeddings (gene_emb):
       - Each gene gets a learnable vector representing its "identity"
       - Similar to word embeddings in NLP
       - Shape: (n_genes, d_model)
    
    2. Value Projection (val_proj):
       - Projects scalar expression values into d_model dimensions
       - Allows the model to process expression magnitudes
    
    3. Transformer Encoder (encoder):
       - Stack of self-attention layers
       - Learns which genes "pay attention to" other genes
       - Captures complex gene-gene relationships
    
    4. Output Layer (out):
       - Used during training for reconstruction loss
       - Not used for embedding extraction
    
    Parameters
    ----------
    n_genes : int
        Number of genes in the expression matrix (vocabulary size)
    d_model : int, default=256
        Dimensionality of the embedding space
    n_layers : int, default=6
        Number of Transformer encoder layers (depth)
    n_heads : int, default=8
        Number of attention heads (parallelism in attention)
    dropout : float, default=0.1
        Dropout rate for regularization during training
    """
    
    def __init__(self, n_genes, d_model=256, n_layers=6, n_heads=8, dropout=0.1):
        super().__init__()
        
        # ---------------------------------------------------------------------
        # LEARNABLE GENE IDENTITY EMBEDDINGS
        # ---------------------------------------------------------------------
        # Each gene gets a unique d_model-dimensional vector that the model
        # learns during training. This allows the model to distinguish between
        # genes and learn gene-specific patterns.
        #
        # Think of it like: "Gene TP53 has this personality, BRCA1 has that one"
        self.gene_emb = nn.Embedding(n_genes, d_model)
        
        # ---------------------------------------------------------------------
        # EXPRESSION VALUE PROJECTION
        # ---------------------------------------------------------------------
        # Raw expression values are scalars (e.g., 5.2, -1.3).
        # We project them to d_model dimensions so they can be combined with
        # gene embeddings.
        #
        # This is like converting "expression level" into a rich representation
        # that can interact with gene identity information.
        self.val_proj = nn.Linear(1, d_model)
        
        # ---------------------------------------------------------------------
        # TRANSFORMER ENCODER
        # ---------------------------------------------------------------------
        # The core of the model: self-attention layers that learn gene-gene
        # relationships.
        #
        # Key hyperparameters:
        # - d_model: Hidden dimension throughout the model
        # - nhead: Number of parallel attention heads
        # - dim_feedforward: Size of the MLP inside each layer (typically 4×d_model)
        # - activation="gelu": Smooth activation function (better than ReLU)
        # - batch_first=True: Input shape is (batch, sequence, features)
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
        # OUTPUT PROJECTION (for training)
        # ---------------------------------------------------------------------
        # Projects back to scalar values for reconstruction loss during training.
        # Not used when extracting embeddings.
        self.out = nn.Linear(d_model, 1)
        
        # ---------------------------------------------------------------------
        # GENE INDEX BUFFER
        # ---------------------------------------------------------------------
        # Pre-computed tensor [0, 1, 2, ..., n_genes-1] for looking up gene
        # embeddings. Using register_buffer ensures this moves to GPU with the
        # model but isn't treated as a learnable parameter.
        #
        # persistent=False: Don't save this in checkpoints (it's trivial to recreate)
        self.register_buffer(
            "gene_ids",
            torch.arange(n_genes, dtype=torch.long),
            persistent=False
        )
    
    @torch.no_grad()  # Disable gradient computation for inference (saves memory)
    def embed(self, x):
        """
        Extract embeddings from expression data.
        
        This is the key method for generating representations we can use
        for downstream analysis (clustering, visualization, etc.).
        
        Parameters
        ----------
        x : torch.Tensor
            Expression matrix of shape (batch_size, n_genes).
            Each row is a sample, each column is a gene's expression value.
        
        Returns
        -------
        torch.Tensor
            Embeddings of shape (batch_size, d_model).
            Each row is the learned representation for one sample.
        
        How it works
        ------------
        1. For each gene position, look up its identity embedding
        2. Project the expression value to d_model dimensions
        3. ADD them together: final_repr = value_info + gene_identity
        4. Pass through Transformer to learn cross-gene relationships
        5. Mean-pool across all genes to get a single vector per sample
        """
        B, G = x.shape  # B = batch size, G = number of genes
        
        # Expand gene IDs for the entire batch
        # Shape: (1, G) -> (B, G)
        gene_ids = self.gene_ids.unsqueeze(0).expand(B, G)
        
        # Combine expression values with gene identities
        # x.unsqueeze(-1): (B, G) -> (B, G, 1) to allow projection
        # val_proj output: (B, G, d_model)
        # gene_emb output: (B, G, d_model)
        # h (combined): (B, G, d_model)
        h = self.val_proj(x.unsqueeze(-1)) + self.gene_emb(gene_ids)
        
        # Pass through Transformer encoder
        # Self-attention allows each gene to "see" all other genes
        # Output shape: (B, G, d_model)
        h = self.encoder(h)
        
        # Mean pooling: aggregate information across all genes
        # (B, G, d_model) -> (B, d_model)
        # This produces one fixed-size vector per sample regardless of G
        return h.mean(dim=1)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # DEVICE SELECTION
    # -------------------------------------------------------------------------
    # Use GPU if available (much faster), otherwise fall back to CPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # -------------------------------------------------------------------------
    # LOAD INPUT DATA
    # -------------------------------------------------------------------------
    # NPZ is a NumPy archive format that can store multiple arrays
    # X: Expression matrix, shape (n_samples, n_genes)
    data = np.load(NPZ_PATH)
    X = torch.from_numpy(data["X"]).to(device)
    print(f"Loaded expression data: {X.shape[0]} samples × {X.shape[1]} genes")
    
    # Load metadata (sample IDs for labeling output rows)
    with open(META_PATH) as f:
        meta = json.load(f)
    samples = meta["samples"]
    
    # -------------------------------------------------------------------------
    # LOAD PRETRAINED MODEL
    # -------------------------------------------------------------------------
    # The checkpoint contains:
    # - model_state: Learned weights from training
    # - config: Hyperparameters (d_model, n_layers, etc.)
    # - G: Number of genes the model was trained on
    ckpt = torch.load(CKPT_PATH, map_location=device, weights_only=False)
    
    cfg = ckpt["config"]  # Model architecture hyperparameters
    G = ckpt["G"]         # Number of genes (must match input data)
    
    # Instantiate model with saved configuration
    model = ExprTransformer(G, **cfg).to(device)
    
    # Load the trained weights
    # strict=True: Fail if there's any mismatch (good for catching errors)
    model.load_state_dict(ckpt["model_state"], strict=True)
    
    # Set to evaluation mode
    # - Disables dropout (use all neurons)
    # - Sets BatchNorm to use running statistics (not relevant here)
    model.eval()
    
    print(f"Loaded model: {G} genes, d_model={cfg.get('d_model', 256)}")
    
    # -------------------------------------------------------------------------
    # BATCH INFERENCE
    # -------------------------------------------------------------------------
    # Process samples in batches to avoid running out of GPU memory.
    # GPU can handle larger batches; CPU needs smaller ones.
    batch_size = 512 if device == "cuda" else 64
    
    embeddings_list = []
    n_samples = X.shape[0]
    
    print(f"Generating embeddings in batches of {batch_size}...")
    for i in range(0, n_samples, batch_size):
        # Extract batch
        batch = X[i:i + batch_size]
        
        # Generate embeddings for this batch
        # .detach(): Ensure no gradient tracking (redundant with @torch.no_grad but explicit)
        # .cpu(): Move back to CPU for numpy conversion
        # .numpy(): Convert to numpy array
        batch_embeddings = model.embed(batch).detach().cpu().numpy()
        embeddings_list.append(batch_embeddings)
        
        # Progress indicator
        if (i // batch_size) % 10 == 0:
            print(f"  Processed {min(i + batch_size, n_samples)}/{n_samples} samples")
    
    # Stack all batches into a single array
    # Shape: (n_samples, d_model)
    E = np.vstack(embeddings_list)
    
    # -------------------------------------------------------------------------
    # SAVE OUTPUT
    # -------------------------------------------------------------------------
    # Create a DataFrame with:
    # - Index: Sample IDs (e.g., "TCGA-A1-A0SK-01A", "ACC-201")
    # - Columns: Embedding dimensions (z1, z2, ..., z256)
    df = pd.DataFrame(
        E,
        index=samples,
        columns=[f"z{i+1}" for i in range(E.shape[1])]
    )
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    
    # Save as tab-separated values (TSV)
    # TSV is easy to load in R, Python, or Excel
    df.to_csv(OUT_TSV, sep="\t")
    
    print(f"\n{'='*60}")
    print(f"SUCCESS: Saved embeddings to {OUT_TSV}")
    print(f"Shape: {df.shape[0]} samples × {df.shape[1]} dimensions")
    print(f"{'='*60}")