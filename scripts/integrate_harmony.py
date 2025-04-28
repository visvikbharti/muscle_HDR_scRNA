#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import yaml
try:
    import harmonypy as hm
except ImportError:
    import harmony as hm  # Some installations use different import name

def load_config():
    """Load configuration from YAML file"""
    config_path = 'config.yaml'
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    else:
        return {
            'integration': {
                'method': 'harmony',
                'n_neighbors': 15
            }
        }

def harmonize_genes(adata_human, adata_mouse, ortholog_table):
    """Harmonize gene names between human and mouse using ortholog table"""
    # Load ortholog table
    orthologs = pd.read_csv(ortholog_table, sep='\t')
    
    # Create mapping dictionary
    human_to_mouse = dict(zip(orthologs['human_symbol'], orthologs['mouse_symbol']))
    mouse_to_human = dict(zip(orthologs['mouse_symbol'], orthologs['human_symbol']))
    
    # Convert mouse genes to human nomenclature where possible
    mouse_genes = adata_mouse.var_names
    converted_genes = []
    
    for gene in mouse_genes:
        # Try exact match
        if gene in mouse_to_human:
            converted_genes.append(mouse_to_human[gene])
        # Try case-insensitive match
        elif gene.upper() in [g.upper() for g in mouse_to_human.keys()]:
            matched_gene = [g for g in mouse_to_human.keys() if g.upper() == gene.upper()][0]
            converted_genes.append(mouse_to_human[matched_gene])
        else:
            # Keep original if no match found
            converted_genes.append(gene)
    
    # Update mouse gene names
    adata_mouse.var_names = converted_genes
    adata_mouse.var_names_make_unique()
    
    # Find common genes
    common_genes = list(set(adata_human.var_names) & set(adata_mouse.var_names))
    print(f"Found {len(common_genes)} common genes between human and mouse")
    
    # Subset to common genes
    adata_human = adata_human[:, common_genes].copy()
    adata_mouse = adata_mouse[:, common_genes].copy()
    
    return adata_human, adata_mouse

def integrate_datasets(adata_human, adata_mouse, config):
    """Integrate human and mouse datasets using Harmony"""
    # Add batch information
    adata_human.obs['species'] = 'human'
    adata_mouse.obs['species'] = 'mouse'
    
    # Add technology information
    adata_human.obs['technology'] = 'SMART-seq2'
    adata_mouse.obs['technology'] = '10x'
    
    # Concatenate datasets
    adata_combined = anndata.concat([adata_human, adata_mouse], 
                                  label='dataset', 
                                  keys=['human', 'mouse'],
                                  join='outer',
                                  fill_value=0)
    
    # Create batch key combining species and technology
    adata_combined.obs['batch'] = adata_combined.obs['species'] + '_' + adata_combined.obs['technology']
    
    # Normalize and log transform
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000)
    adata_combined = adata_combined[:, adata_combined.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata_combined, max_value=10)
    
    # PCA
    sc.tl.pca(adata_combined, svd_solver='arpack')
    
    # Run Harmony
    ho = hm.run_harmony(adata_combined.obsm['X_pca'], 
                       adata_combined.obs, 
                       ['batch', 'species', 'technology'])
    
    # Store corrected PCA embeddings
    adata_combined.obsm['X_pca_harmony'] = ho.Z_corr.T
    
    # Compute neighbors and UMAP on integrated space
    sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony', 
                    n_neighbors=config['integration'].get('n_neighbors', 15))
    sc.tl.umap(adata_combined)
    
    # For leiden clustering, ensure we have the right package
    try:
        sc.tl.leiden(adata_combined, resolution=0.5)
    except Exception as e:
        print(f"Leiden clustering failed: {e}")
        # Use louvain as fallback
        sc.tl.louvain(adata_combined, resolution=0.5)
        adata_combined.obs['leiden'] = adata_combined.obs['louvain']
    
    return adata_combined

def main():
    parser = argparse.ArgumentParser(description='Integrate human and mouse datasets using Harmony')
    parser.add_argument('human_h5ad', help='Path to human h5ad file')
    parser.add_argument('mouse_h5ad', help='Path to mouse h5ad file')
    parser.add_argument('ortholog_table', help='Path to ortholog table')
    parser.add_argument('output_h5ad', help='Path to output integrated h5ad file')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config()
    
    # Load datasets
    print("Loading datasets...")
    adata_human = sc.read_h5ad(args.human_h5ad)
    adata_mouse = sc.read_h5ad(args.mouse_h5ad)
    
    # Harmonize gene names
    print("Harmonizing gene names...")
    adata_human, adata_mouse = harmonize_genes(adata_human, adata_mouse, args.ortholog_table)
    
    # Integrate datasets
    print("Integrating datasets with Harmony...")
    adata_integrated = integrate_datasets(adata_human, adata_mouse, config)
    
    # Save integrated dataset
    print(f"Saving integrated dataset to {args.output_h5ad}")
    adata_integrated.write_h5ad(args.output_h5ad)
    
    # Print summary
    print("\nIntegration summary:")
    print(f"Total cells: {adata_integrated.n_obs}")
    print(f"Total genes: {adata_integrated.n_vars}")
    print(f"Batches: {adata_integrated.obs['batch'].unique()}")
    print(f"Clusters: {len(adata_integrated.obs['leiden'].unique())}")

if __name__ == "__main__":
    main()