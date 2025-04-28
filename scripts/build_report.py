#!/usr/bin/env python3
import anndata
import yaml
import pandas as pd
import scanpy as sc
import numpy as np
import os
import subprocess
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

def load_data_and_config():
    """Load AnnData objects and configuration"""
    adata_h = anndata.read_h5ad("data/processed/GSE130646.h5ad")
    adata_m = anndata.read_h5ad("data/processed/GSE138707.h5ad")
    
    with open("config.yaml", 'r') as f:
        config = yaml.safe_load(f)
    
    return adata_h, adata_m, config

def calculate_cluster_sizes(adata):
    """Calculate cluster sizes"""
    return adata.obs['leiden'].value_counts().sort_index()

def find_top_markers(adata, n_genes=5):
    """Find top marker genes for each cluster"""
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    markers = {}
    
    for cluster in adata.obs['leiden'].unique():
        genes = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        top_genes = genes.head(n_genes)['names'].tolist()
        markers[cluster] = ', '.join(top_genes)
    
    return markers

def calculate_hdr_scores(adata, hdr_genes):
    """Calculate HDR scores for each cluster"""
    # Filter for HDR genes present in dataset
    present_genes = [g for g in hdr_genes if g in adata.var_names]
    
    if present_genes:
        sc.tl.score_genes(adata, present_genes, score_name='HDR_score')
        
        # Calculate mean HDR score per cluster
        cluster_scores = adata.obs.groupby('leiden')['HDR_score'].agg(['mean', 'std'])
        return cluster_scores
    
    return None

def create_hdr_heatmap(adata, hdr_genes, output_path):
    """Create HDR gene heatmap"""
    present_genes = [g for g in hdr_genes if g in adata.var_names]
    
    if present_genes:
        # Select top 15 HDR genes by variance
        var_genes = adata[:, present_genes].var.copy()
        var_genes['variance'] = adata[:, present_genes].X.var(axis=0)
        top_genes = var_genes.nlargest(15, 'variance').index.tolist()
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.matrixplot(adata, var_names=top_genes, groupby='leiden',
                        standard_scale='var', cmap='viridis', 
                        show=False, ax=ax)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

def get_git_hash():
    """Get current git commit hash"""
    try:
        return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    except:
        return "unknown"

def build_report(adata_h, adata_m, config):
    """Build the final report"""
    # Load HDR genes
    with open("docs/HDR_477_genes.txt", 'r') as f:
        hdr_genes = [line.strip() for line in f if line.strip()]
    
    # Calculate metrics
    cluster_sizes_h = calculate_cluster_sizes(adata_h)
    cluster_sizes_m = calculate_cluster_sizes(adata_m)
    
    markers_h = find_top_markers(adata_h)
    markers_m = find_top_markers(adata_m)
    
    hdr_scores_h = calculate_hdr_scores(adata_h, hdr_genes)
    
    # Create figures
    os.makedirs("results/figures", exist_ok=True)
    create_hdr_heatmap(adata_h, hdr_genes, "results/figures/hdr_heatmap_human.png")
    
    # Generate report content
    report = f"""# Homology-Directed-Repair Gene Landscape in Human and Mouse Skeletal Muscle at Single-Cell Resolution

Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 1. Data Provenance

| Species | GEO Accession | Library prep | Cells (post-QC) | Genes kept | UMAP included |
|---------|---------------|-------------|-----------------|------------|---------------|
| Homo sapiens | GSE130646 | SMART-seq2 | {adata_h.n_obs} | {adata_h.n_vars} | Yes |
| Mus musculus | GSE138707 | 10x v3 | {adata_m.n_obs} | {adata_m.n_vars} | Yes |

## 2. Processing Parameters

| Step | Parameter | Value |
|------|-----------|-------|
| Gene-filter (cells) | min_genes | {config['qc']['min_genes']} |
| Cell-filter (genes) | min_cells | {config['qc']['min_cells']} |
| Mito cutoff | pct_mt_max | {config['qc']['pct_mt_max']}% |
| HVGs | n_top_genes | {config['processing']['n_top_genes']} |
| Neighbours | k | {config['clustering']['n_neighbors']} |
| Leiden resolution | res | {config['clustering']['leiden_resolution']} |

## 3. Cluster Summary - Human

| Leiden | Size | Key markers |
|--------|------|-------------|
"""
    
    for cluster in sorted(cluster_sizes_h.index):
        report += f"| {cluster} | {cluster_sizes_h[cluster]} | {markers_h[cluster]} |\n"
    
    report += f"""
## 4. HDR Gene Expression (Human)

![HDR Gene Heatmap](figures/hdr_heatmap_human.png)

### Key Observations

"""
    
    if hdr_scores_h is not None:
        max_cluster = hdr_scores_h['mean'].idxmax()
        min_cluster = hdr_scores_h['mean'].idxmin()
        report += f"- Cluster {max_cluster} shows highest HDR gene expression (mean score: {hdr_scores_h.loc[max_cluster, 'mean']:.2f})\n"
        report += f"- Cluster {min_cluster} shows lowest HDR gene expression (mean score: {hdr_scores_h.loc[min_cluster, 'mean']:.2f})\n"
    
    report += f"""

## 5. Reproducibility

- Git commit: {get_git_hash()}
- Configuration: See config.yaml
- Environment: See env/requirements.txt

## 6. Methods

Single-cell RNA sequencing data was processed using Scanpy {sc.__version__}. 
Quality control filtered cells with <{config['qc']['min_genes']} genes and 
>{config['qc']['pct_mt_max']}% mitochondrial content. 
Data was normalized to {config['processing']['normalize_target_sum']} counts per cell, 
followed by log transformation. Highly variable genes were identified ({config['processing']['n_top_genes']} genes) 
and used for PCA. Nearest neighbors were computed (k={config['clustering']['n_neighbors']}) 
followed by Leiden clustering (resolution={config['clustering']['leiden_resolution']}) and UMAP visualization.
"""
    
    return report

def main():
    """Main function"""
    # Load data
    adata_h, adata_m, config = load_data_and_config()
    
    # Build report
    report = build_report(adata_h, adata_m, config)
    
    # Save report
    with open("results/final_report.md", 'w') as f:
        f.write(report)
    
    print("Report generated successfully!")

if __name__ == "__main__":
    main()