#!/usr/bin/env python3
import sys
import pandas as pd
import scanpy as sc
import anndata
import yaml
import os

def load_config():
    """Load configuration from YAML file"""
    config_path = 'config.yaml'
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    else:
        # Default values if config file doesn't exist
        return {
            'qc': {
                'min_genes': 200,
                'min_cells': 3,
                'pct_mt_max': 15
            },
            'processing': {
                'n_top_genes': 2000,
                'normalize_target_sum': 10000,
                'scale_max_value': 10
            },
            'clustering': {
                'n_neighbors': 15,
                'n_pcs': 20,
                'leiden_resolution': 0.5,
                'method': 'gauss'
            }
        }

def main(csv_file, output_h5ad):
    # Load configuration
    config = load_config()
    
    print(f"Reading counts from {csv_file}")
    counts = pd.read_csv(csv_file, index_col=0, compression='gzip')
    
    # Create AnnData object (genes x cells -> cells x genes)
    adata = anndata.AnnData(counts.T)
    
    # Fix duplicate variable names
    adata.var_names_make_unique()
    
    # Basic filtering with config parameters
    sc.pp.filter_cells(adata, min_genes=config['qc']['min_genes'])
    sc.pp.filter_genes(adata, min_cells=config['qc']['min_cells'])
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Normalize with config parameters
    sc.pp.normalize_total(adata, target_sum=config['processing']['normalize_target_sum'])
    sc.pp.log1p(adata)
    
    # Store raw counts for later use
    adata.raw = adata
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=config['processing']['n_top_genes'])
    
    # Scale data
    adata_hvg = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata_hvg, max_value=config['processing']['scale_max_value'])
    
    # PCA
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    
    # Compute neighbors
    sc.pp.neighbors(adata_hvg, 
                    n_neighbors=config['clustering']['n_neighbors'], 
                    n_pcs=config['clustering']['n_pcs'],
                    method=config['clustering']['method'])
    
    # UMAP and clustering
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=config['clustering']['leiden_resolution'])
    
    # Transfer results back to the full adata
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsp['distances'] = adata_hvg.obsp['distances']
    adata.obsp['connectivities'] = adata_hvg.obsp['connectivities']
    adata.uns['neighbors'] = adata_hvg.uns['neighbors']
    adata.uns['umap'] = adata_hvg.uns['umap']
    adata.uns['leiden'] = adata_hvg.uns['leiden']
    adata.uns['pca'] = adata_hvg.uns['pca']
    
    # Store processing parameters in uns
    adata.uns['processing_params'] = config
    
    # Save
    print(f"Saving to {output_h5ad}")
    adata.write_h5ad(output_h5ad)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])