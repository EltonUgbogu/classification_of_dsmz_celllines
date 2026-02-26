#!/usr/bin/env python3
"""
Build pan-cancer feature set from disease-specific marker sets.

Merge strategy:
1. PAN_CORE: genes present in ≥2 diseases (using within-disease recurrence ≥2)
2. Disease-specific: capped per-disease genes (balanced quotas)
3. Remove ribosomal/mitochondrial families (RPL/RPS/MT genes)
"""

import os
import sys
from pathlib import Path
from collections import defaultdict
import argparse

# Default parameters
DEFAULT_CAP_SPECIFIC = 300  # per disease
DEFAULT_CAP_ISOLATE = 50    # per disease
DEFAULT_MIN_CROSS_DISEASE = 2  # genes must appear in ≥2 diseases for PAN_CORE


def load_gene_list(filepath):
    """Load gene IDs from a file (one per line)."""
    genes = set()
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found, skipping", file=sys.stderr)
        return genes
    
    with open(filepath, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    
    return genes


def load_recurrence_table(filepath):
    """Load gene recurrence table (TSV with gene_id and freq columns)."""
    gene_freq = {}
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found, skipping", file=sys.stderr)
        return gene_freq
    
    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        if 'gene_id' not in header or 'freq' not in header:
            print(f"Warning: {filepath} doesn't have expected columns", file=sys.stderr)
            return gene_freq
        
        gene_idx = header.index('gene_id')
        freq_idx = header.index('freq')
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(gene_idx, freq_idx):
                gene_id = fields[gene_idx]
                try:
                    freq = int(fields[freq_idx])
                    gene_freq[gene_id] = freq
                except ValueError:
                    continue
    
    return gene_freq


def is_ribosomal_or_mitochondrial(gene_id):
    """Check if gene is ribosomal or mitochondrial."""
    # Check gene symbol patterns (we'll need to map or use annotation)
    # For now, check common patterns in ENSEMBL IDs or use gene symbol lookup
    # This is a simplified check - you may want to use a gene annotation file
    gene_lower = gene_id.lower()
    
    # Common ribosomal gene patterns in ENSEMBL annotations
    # Note: This is a heuristic. For production, use proper gene annotation.
    ribosomal_patterns = ['rpl', 'rps', 'rpm', 'rpn', 'mrpl', 'mrps']
    mitochondrial_patterns = ['mt-', 'mt_']
    
    # Since we have ENSEMBL IDs, we can't directly check gene symbols
    # We'll need to either:
    # 1. Load a gene annotation file (GTF/GFF)
    # 2. Use a gene ID to symbol mapping
    # 3. Skip this check and do it later with annotation
    
    # For now, return False - we'll filter later if needed
    return False


def filter_ribosomal_mitochondrial(genes, annotation_file=None):
    """
    Filter out ribosomal and mitochondrial genes.
    
    If annotation_file is provided, use it to map gene IDs to symbols.
    Otherwise, return genes as-is (filtering can be done later).
    """
    if annotation_file and os.path.exists(annotation_file):
        # Load annotation and filter
        # This would require parsing GTF/GFF or a gene ID mapping file
        pass
    
    # For now, we'll do a simple check on gene IDs
    # Most ribosomal genes in human have specific ENSEMBL patterns
    # But without annotation, we can't reliably filter
    # Return all genes for now - user can filter later with proper annotation
    return genes


def build_pan_cancer_features(
    brca_markers_dir,
    nbl_markers_dir,
    rbl_markers_dir,
    output_file,
    min_cross_disease=DEFAULT_MIN_CROSS_DISEASE,
    cap_specific=DEFAULT_CAP_SPECIFIC,
    cap_isolate=DEFAULT_CAP_ISOLATE,
    remove_ribo_mt=True
):
    """
    Build pan-cancer feature set.
    
    Parameters:
    -----------
    brca_markers_dir : str
        Path to BRCA DESeq2 markers directory
    nbl_markers_dir : str
        Path to NBL DESeq2 markers directory
    rbl_markers_dir : str
        Path to RBL DESeq2 markers directory
    output_file : str
        Output file path for pan-cancer features
    min_cross_disease : int
        Minimum number of diseases a gene must appear in for PAN_CORE
    cap_specific : int
        Maximum disease-specific genes per disease
    cap_isolate : int
        Maximum isolate-specific genes per disease
    remove_ribo_mt : bool
        Whether to remove ribosomal/mitochondrial genes
    """
    
    # Disease directories
    diseases = {
        'BRCA': brca_markers_dir,
        'NBL': nbl_markers_dir,
        'RBL': rbl_markers_dir
    }
    
    # Load unique feature sets (recurrence ≥2) for each disease
    disease_genes = {}
    disease_recurrence = {}
    
    for disease, markers_dir in diseases.items():
        unique_file = os.path.join(markers_dir, 'unique_feature_set_recurrence_ge_2.txt')
        recurrence_file = os.path.join(markers_dir, 'gene_recurrence_across_contrasts.tsv')
        
        print(f"Loading {disease} markers from {markers_dir}...", file=sys.stderr)
        
        # Load unique genes (recurrence ≥2)
        disease_genes[disease] = load_gene_list(unique_file)
        print(f"  Loaded {len(disease_genes[disease])} unique genes (recurrence ≥2)", file=sys.stderr)
        
        # Load recurrence frequencies for ranking
        disease_recurrence[disease] = load_recurrence_table(recurrence_file)
        print(f"  Loaded recurrence for {len(disease_recurrence[disease])} genes", file=sys.stderr)
    
    # Step 1: Build PAN_CORE (genes in ≥min_cross_disease diseases)
    print(f"\nBuilding PAN_CORE (genes in ≥{min_cross_disease} diseases)...", file=sys.stderr)
    
    gene_disease_count = defaultdict(int)
    for disease, genes in disease_genes.items():
        for gene in genes:
            gene_disease_count[gene] += 1
    
    pan_core = {
        gene for gene, count in gene_disease_count.items()
        if count >= min_cross_disease
    }
    
    print(f"  PAN_CORE: {len(pan_core)} genes", file=sys.stderr)
    
    # Step 2: Build disease-specific sets (capped)
    print(f"\nBuilding disease-specific sets (cap={cap_specific} per disease)...", file=sys.stderr)
    
    disease_specific = {}
    
    for disease, genes in disease_genes.items():
        # Get genes that are NOT in PAN_CORE
        specific_genes = genes - pan_core
        
        # Rank by recurrence frequency (higher = better)
        # If gene not in recurrence table, use frequency 1
        ranked = sorted(
            specific_genes,
            key=lambda g: disease_recurrence[disease].get(g, 1),
            reverse=True
        )
        
        # Cap at cap_specific
        disease_specific[disease] = set(ranked[:cap_specific])
        print(f"  {disease}-specific: {len(disease_specific[disease])} genes (from {len(specific_genes)} candidates)", file=sys.stderr)
    
    # Step 3: Load isolate markers (if available, capped separately)
    print(f"\nLoading isolate markers (cap={cap_isolate} per disease)...", file=sys.stderr)
    
    isolate_genes = {}
    selected_genes = set(pan_core)
    for genes in disease_specific.values():
        selected_genes.update(genes)
    for disease, markers_dir in diseases.items():
        isolate_files = [
            f for f in os.listdir(markers_dir)
            if f.startswith('isolate_') and f.endswith('_markers_top50.txt')
        ]
        
        disease_isolates = set()
        for isolate_file in isolate_files:
            filepath = os.path.join(markers_dir, isolate_file)
            genes = load_gene_list(filepath)
            disease_isolates.update(genes)
        
        # Remove genes already selected (PAN_CORE + all disease-specific + prior isolates)
        disease_isolates = disease_isolates - selected_genes
        
        # Cap at cap_isolate, ensuring global uniqueness
        ranked_isolates = sorted(
            disease_isolates,
            key=lambda g: disease_recurrence[disease].get(g, 1),
            reverse=True
        )
        chosen = []
        for gene in ranked_isolates:
            if gene in selected_genes:
                continue
            chosen.append(gene)
            if len(chosen) >= cap_isolate:
                break
        isolate_genes[disease] = set(chosen)
        selected_genes.update(isolate_genes[disease])
        print(f"  {disease} isolates: {len(isolate_genes[disease])} genes (from {len(disease_isolates)} candidates)", file=sys.stderr)
    
    # Step 4: Combine all features
    all_features = pan_core.copy()
    for genes in disease_specific.values():
        all_features.update(genes)
    for genes in isolate_genes.values():
        all_features.update(genes)
    
    print(f"\nTotal features before filtering: {len(all_features)}", file=sys.stderr)
    
    # Step 5: Filter ribosomal/mitochondrial (if requested)
    if remove_ribo_mt:
        print("\nFiltering ribosomal/mitochondrial genes...", file=sys.stderr)
        print("  Note: This requires gene annotation. Skipping for now.", file=sys.stderr)
        print("  You can filter later using a gene annotation file.", file=sys.stderr)
        # For now, we keep all genes
        # TODO: Add proper gene annotation filtering
    
    # Step 6: Write output
    print(f"\nWriting {len(all_features)} features to {output_file}...", file=sys.stderr)
    
    with open(output_file, 'w') as f:
        # Write PAN_CORE first
        f.write("# PAN_CORE (genes in ≥{} diseases): {} genes\n".format(
            min_cross_disease, len(pan_core)
        ))
        for gene in sorted(pan_core):
            f.write(f"{gene}\n")
        
        # Write disease-specific
        for disease in sorted(disease_specific.keys()):
            genes = disease_specific[disease]
            if genes:
                f.write(f"\n# {disease}-specific (capped at {cap_specific}): {len(genes)} genes\n")
                for gene in sorted(genes):
                    f.write(f"{gene}\n")
        
        # Write isolate genes
        for disease in sorted(isolate_genes.keys()):
            genes = isolate_genes[disease]
            if genes:
                f.write(f"\n# {disease}-isolates (capped at {cap_isolate}): {len(genes)} genes\n")
                for gene in sorted(genes):
                    f.write(f"{gene}\n")
    
    # Also write a clean version (just gene IDs, no comments)
    clean_output = output_file.replace('.txt', '_clean.txt')
    with open(clean_output, 'w') as f:
        for gene in sorted(all_features):
            f.write(f"{gene}\n")
    
    print(f"  Also wrote clean version to {clean_output}", file=sys.stderr)
    
    # Print summary
    print("\n" + "="*60, file=sys.stderr)
    print("PAN-CANCER FEATURE SET SUMMARY", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print(f"PAN_CORE (≥{min_cross_disease} diseases): {len(pan_core)} genes", file=sys.stderr)
    for disease in sorted(disease_specific.keys()):
        print(f"{disease}-specific: {len(disease_specific[disease])} genes", file=sys.stderr)
        print(f"{disease}-isolates: {len(isolate_genes[disease])} genes", file=sys.stderr)
    print(f"TOTAL: {len(all_features)} genes", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    return all_features, pan_core, disease_specific, isolate_genes


def main():
    parser = argparse.ArgumentParser(
        description='Build pan-cancer feature set from disease-specific markers'
    )
    parser.add_argument(
        '--brca-markers',
        required=True,
        help='Path to BRCA DESeq2 markers directory'
    )
    parser.add_argument(
        '--nbl-markers',
        required=True,
        help='Path to NBL DESeq2 markers directory'
    )
    parser.add_argument(
        '--rbl-markers',
        required=True,
        help='Path to RBL DESeq2 markers directory'
    )
    parser.add_argument(
        '--output',
        default='pan_cancer_features.txt',
        help='Output file path (default: pan_cancer_features.txt)'
    )
    parser.add_argument(
        '--min-cross-disease',
        type=int,
        default=DEFAULT_MIN_CROSS_DISEASE,
        help=f'Minimum diseases for PAN_CORE (default: {DEFAULT_MIN_CROSS_DISEASE})'
    )
    parser.add_argument(
        '--cap-specific',
        type=int,
        default=DEFAULT_CAP_SPECIFIC,
        help=f'Cap for disease-specific genes per disease (default: {DEFAULT_CAP_SPECIFIC})'
    )
    parser.add_argument(
        '--cap-isolate',
        type=int,
        default=DEFAULT_CAP_ISOLATE,
        help=f'Cap for isolate genes per disease (default: {DEFAULT_CAP_ISOLATE})'
    )
    parser.add_argument(
        '--keep-ribo-mt',
        action='store_true',
        help='Keep ribosomal/mitochondrial genes (default: will attempt to filter)'
    )
    
    args = parser.parse_args()
    
    build_pan_cancer_features(
        args.brca_markers,
        args.nbl_markers,
        args.rbl_markers,
        args.output,
        min_cross_disease=args.min_cross_disease,
        cap_specific=args.cap_specific,
        cap_isolate=args.cap_isolate,
        remove_ribo_mt=not args.keep_ribo_mt
    )


if __name__ == '__main__':
    main()
