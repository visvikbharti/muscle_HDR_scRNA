#!/usr/bin/env python3
import argparse
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import anndata
import os
import yaml

def load_config():
    """Load configuration from YAML file"""
    config_path = 'config.yaml'
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    else:
        return {
            'integration': {
                'method': 'scvi',
                'n_latent': 30,
                'n_layers': 2,
                'n_epochs': 400
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
    """Integrate human and mouse datasets using scVI"""
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
    
    # Store raw counts for scVI
    adata_combined.layers['counts'] = adata_combined.X.copy()
    
    # Normalize and log transform for visualization
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    
    # Set up scVI model
    scvi.model.SCVI.setup_anndata(
        adata_combined,
        layer='counts',
        batch_key='batch',
        categorical_covariate_keys=['species', 'technology']
    )
    
    # Create and train model
    model = scvi.model.SCVI(
        adata_combined,
        n_layers=config['integration']['n_layers'],
        n_latent=config['integration']['n_latent'],
        gene_likelihood='nb'  # negative binomial for count data
    )
    
    # Train the model
    model.train(
        max_epochs=config['integration']['n_epochs'],
        early_stopping=True,
        early_stopping_patience=20
    )
    
    # Get integrated latent representation
    adata_combined.obsm['X_scVI'] = model.get_latent_representation()
    
    # Get normalized expression values
    adata_combined.layers['scvi_normalized'] = model.get_normalized_expression()
    
    # Compute neighbors and UMAP on integrated space
    sc.pp.neighbors(adata_combined, use_rep='X_scVI', n_neighbors=15)
    sc.tl.umap(adata_combined)
    sc.tl.leiden(adata_combined, resolution=0.5)
    
    # Save the model for future use
    model_dir = 'models/scvi_integration'
    os.makedirs(model_dir, exist_ok=True)
    model.save(model_dir, overwrite=True)
    
    # Store model path in uns
    adata_combined.uns['scvi_model_path'] = model_dir
    
    return adata_combined

def main():
    parser = argparse.ArgumentParser(description='Integrate human and mouse datasets using scVI')
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
    print("Integrating datasets with scVI...")
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