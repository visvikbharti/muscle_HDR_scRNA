# Enhanced Streamlit app with batch correction
import streamlit as st
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import warnings
import yaml
import json
warnings.filterwarnings('ignore')


def plot_hdr_by_cluster(adata, hdr_genes):
    """Create detailed HDR gene expression plots by cluster"""
    st.header("HDR Gene Expression Analysis")
    
    # Filter for present HDR genes
    present_genes = [gene for gene in hdr_genes if gene in adata.var_names]
    
    # Calculate HDR score if not present
    if 'HDR_score' not in adata.obs:
        sc.tl.score_genes(adata, gene_list=present_genes, score_name='HDR_score')
    
    # 1. Bar plot of HDR score by cluster
    st.subheader("HDR Score by Cluster")
    fig, ax = plt.subplots(figsize=(10, 6))
    cluster_means = adata.obs.groupby('leiden')['HDR_score'].mean().sort_values(ascending=False)
    cluster_std = adata.obs.groupby('leiden')['HDR_score'].std()
    
    bars = ax.bar(range(len(cluster_means)), cluster_means.values, yerr=cluster_std.values, 
                  capsize=5, alpha=0.8)
    ax.set_xticks(range(len(cluster_means)))
    ax.set_xticklabels(cluster_means.index)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('HDR Score')
    ax.set_title('Mean HDR Score by Cluster')
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + cluster_std.values[i],
                f'{height:.3f}', ha='center', va='bottom')
    
    st.pyplot(fig)
    
    # 2. Violin plot for HDR score distribution
    st.subheader("HDR Score Distribution by Cluster")
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.violin(adata, keys='HDR_score', groupby='leiden', 
                 rotation=90, ax=ax, show=False)
    st.pyplot(fig)
    
    # 3. Individual gene expression
    st.subheader("Individual HDR Gene Expression")
    selected_gene = st.selectbox("Select HDR gene", present_genes)
    
    if selected_gene:
        # Violin plot for selected gene
        fig, ax = plt.subplots(figsize=(12, 6))
        sc.pl.violin(adata, keys=selected_gene, groupby='leiden', 
                     rotation=90, ax=ax, show=False)
        st.pyplot(fig)
        
        # UMAP colored by selected gene
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color=selected_gene, ax=ax, show=False)
        st.pyplot(fig)
    
    # 4. Top expressed HDR genes by cluster
    st.subheader("Top HDR Genes by Cluster")
    selected_cluster = st.selectbox("Select cluster", sorted(adata.obs['leiden'].unique()))
    
    if selected_cluster:
        cluster_cells = adata[adata.obs['leiden'] == selected_cluster]
        gene_means = pd.Series(
            cluster_cells[:, present_genes].X.mean(axis=0),
            index=present_genes
        ).sort_values(ascending=False)
        
        # Bar plot of top 10 HDR genes
        fig, ax = plt.subplots(figsize=(10, 6))
        top_genes = gene_means.head(10)
        bars = ax.bar(range(len(top_genes)), top_genes.values)
        ax.set_xticks(range(len(top_genes)))
        ax.set_xticklabels(top_genes.index, rotation=45, ha='right')
        ax.set_ylabel('Mean Expression')
        ax.set_title(f'Top 10 HDR Genes in Cluster {selected_cluster}')
        plt.tight_layout()
        st.pyplot(fig)
    
    # 5. Correlation matrix
    st.subheader("HDR Gene Correlation Matrix")
    if len(present_genes) > 1:
        # Select subset of genes for correlation
        n_genes = min(20, len(present_genes))
        top_var_genes = adata[:, present_genes].var.copy()
        top_var_genes['variance'] = adata[:, present_genes].X.var(axis=0)
        selected_genes = top_var_genes.nlargest(n_genes, 'variance').index.tolist()
        
        # Calculate correlation
        expr_matrix = adata[:, selected_genes].X
        if hasattr(expr_matrix, 'toarray'):
            expr_matrix = expr_matrix.toarray()
        corr_matrix = pd.DataFrame(expr_matrix, columns=selected_genes).corr()
        
        # Plot correlation matrix
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(corr_matrix, cmap='coolwarm', center=0, ax=ax)
        plt.title('HDR Gene Correlation Matrix')
        plt.tight_layout()
        st.pyplot(fig)
    
    # 6. Differential expression table
    st.subheader("HDR Gene Differential Expression")
    # Calculate differential expression for HDR genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Create table for each cluster
    results_list = []
    for cluster in adata.obs['leiden'].unique():
        degs = sc.get.rank_genes_groups_df(adata, group=cluster)
        degs_hdr = degs[degs['names'].isin(present_genes)].head(10)
        degs_hdr['cluster'] = cluster
        results_list.append(degs_hdr)
    
    if results_list:
        results_df = pd.concat(results_list)
        st.dataframe(results_df.style.background_gradient(subset=['logfoldchanges'], cmap='RdBu_r'))


st.title("Muscle HDR-scRNA Dashboard")

# Load configuration
@st.cache_data
def load_config():
    with open('config.yaml', 'r') as f:
        return yaml.safe_load(f)

config = load_config()

# Check if processed data exists
required_files = [
    "data/processed/GSE130646.h5ad",
    "data/processed/GSE138707.h5ad"
]

if not all(os.path.exists(f) for f in required_files):
    st.error("Please run the Snakemake pipeline first to generate processed data!")
    st.stop()

# Load data
@st.cache_data
def load_data():
    adata_h = sc.read_h5ad("data/processed/GSE130646.h5ad")
    adata_m = sc.read_h5ad("data/processed/GSE138707.h5ad")
    
    # Load integrated data if available
    integrated_path = "data/processed/integrated_scvi.h5ad"
    adata_integrated = None
    if os.path.exists(integrated_path):
        adata_integrated = sc.read_h5ad(integrated_path)
    
    # Fix variable names duplicates warning
    adata_h.var_names_make_unique()
    adata_m.var_names_make_unique()
    if adata_integrated is not None:
        adata_integrated.var_names_make_unique()
    
    return adata_h, adata_m, adata_integrated

adata_h, adata_m, adata_integrated = load_data()

# Load HDR genes
@st.cache_data
def load_hdr_genes():
    if os.path.exists("docs/HDR_477_genes.txt"):
        with open("docs/HDR_477_genes.txt", "r") as f:
            hdr_genes = [line.strip() for line in f if line.strip()]
        return hdr_genes
    return []

hdr_genes = load_hdr_genes()

# Load ortholog mapping
@st.cache_data
def load_ortholog_mapping():
    mapping_path = "docs/orthologs_mapping.json"
    if os.path.exists(mapping_path):
        with open(mapping_path, 'r') as f:
            return json.load(f)
    return {}

ortholog_map = load_ortholog_mapping()

# Sidebar
st.sidebar.header("Controls")

# Dataset selection with batch correction option
dataset_options = ["Human (GSE130646)", "Mouse (GSE138707)"]
if adata_integrated is not None:
    dataset_options.append("Integrated (Batch Corrected)")

sel = st.sidebar.radio("Dataset", dataset_options)

# Use batch corrected representation if available
use_batch_corrected = st.sidebar.checkbox("Use batch-corrected embeddings", 
                                        value=True, 
                                        disabled=sel == "Integrated (Batch Corrected)")

# Select appropriate data
if sel == "Integrated (Batch Corrected)":
    adata = adata_integrated
    representation = "X_scVI"
else:
    adata = adata_h if sel.startswith("Human") else adata_m
    representation = "X_pca"

# Display basic info
st.write(f"Dataset shape: {adata.shape}")
st.write(f"Number of cells: {adata.n_obs}")
st.write(f"Number of genes: {adata.n_vars}")

# Load enhanced markers
@st.cache_data
def load_markers():
    markers_path = "docs/markers.tsv"
    if os.path.exists(markers_path):
        return pd.read_csv(markers_path, sep='\t')
    return None

markers_df = load_markers()

# Cell type annotation function
def annotate_cell_types(adata, markers_df):
    """Annotate cell types using marker genes"""
    if markers_df is None:
        return None
    
    scores = {}
    for _, row in markers_df.iterrows():
        genes = [g.strip() for g in row['gene_list'].split(',')]
        genes_present = [g for g in genes if g in adata.var_names]
        
        if genes_present:
            sc.tl.score_genes(adata, gene_list=genes_present, score_name=row['cell_type'])
            scores[row['cell_type']] = adata.obs[row['cell_type']]
    
    if scores:
        # Assign cell type based on highest score
        score_df = pd.DataFrame(scores)
        adata.obs['cell_type_auto'] = score_df.idxmax(axis=1)
        adata.obs['confidence_score'] = score_df.max(axis=1)
        
        return score_df
    
    return None

# Main visualization
st.header("Cell Clustering")

# Show UMAP/visualization
if 'X_umap' in adata.obsm:
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Color by different attributes based on dataset type
    if sel == "Integrated (Batch Corrected)":
        color_options = ["leiden", "species", "technology", "batch"]
        selected_color = st.selectbox("Color by:", color_options)
    else:
        color_options = ["leiden"]
        if 'cell_type_auto' in adata.obs:
            color_options.append("cell_type_auto")
        selected_color = st.selectbox("Color by:", color_options)
    
    sc.pl.umap(adata, color=selected_color, legend_loc='on data', 
               title=f'UMAP colored by {selected_color}', show=False, ax=ax)
    st.pyplot(fig)
    
    # If integrated, show species-specific clusters
    if sel == "Integrated (Batch Corrected)" and selected_color == "leiden":
        # Cross-tabulation of clusters vs species
        cross_tab = pd.crosstab(adata.obs['leiden'], adata.obs['species'])
        st.write("Cluster composition by species:")
        st.dataframe(cross_tab)
        
        # Visualize as heatmap
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(cross_tab, annot=True, fmt='d', cmap='YlOrRd', ax=ax)
        ax.set_xlabel('Species')
        ax.set_ylabel('Cluster')
        ax.set_title('Cells per Cluster by Species')
        st.pyplot(fig)
else:
    st.info("UMAP not pre-computed for this dataset")

# HDR Gene Analysis
st.header("HDR Gene Analysis")

if hdr_genes:
    # Handle mouse dataset with ortholog mapping
    if sel.startswith("Mouse") and ortholog_map:
        st.write("Mapping human HDR genes to mouse orthologs...")
        mouse_hdr_genes = []
        for gene in hdr_genes:
            if gene in ortholog_map:
                mouse_hdr_genes.append(ortholog_map[gene])
            else:
                # Try simple case conversion
                mouse_gene = gene.capitalize()
                if mouse_gene in adata.var_names:
                    mouse_hdr_genes.append(mouse_gene)
        
        analysis_genes = mouse_hdr_genes
        st.write(f"Found {len(analysis_genes)} mouse orthologs for HDR genes")
    else:
        analysis_genes = hdr_genes
    
    # Filter for present genes
    genes_present = [gene for gene in analysis_genes if gene in adata.var_names]
    st.write(f"HDR genes found in dataset: {len(genes_present)} out of {len(analysis_genes)}")
    
    if genes_present:
        selected_genes = st.multiselect(
            "Select HDR genes to visualize",
            genes_present,
            default=genes_present[:5] if len(genes_present) > 5 else genes_present
        )
        
        if st.button("Analyze HDR Expression"):
            # Calculate HDR score
            sc.tl.score_genes(adata, gene_list=genes_present, score_name='HDR_score')
            
            # Show HDR score on UMAP
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.umap(adata, color='HDR_score', cmap='viridis', show=False, ax=ax)
            st.pyplot(fig)
            
            # HDR genes by cluster
            if 'leiden' in adata.obs:
                fig, ax = plt.subplots(figsize=(12, 8))
                sc.pl.matrixplot(adata, var_names=selected_genes, groupby='leiden',
                               standard_scale='var', cmap='viridis', show=False, ax=ax)
                st.pyplot(fig)
                
                # If integrated, compare HDR expression between species
                if sel == "Integrated (Batch Corrected)":
                    fig, ax = plt.subplots(figsize=(8, 6))
                    human_mask = adata.obs['species'] == 'human'
                    mouse_mask = adata.obs['species'] == 'mouse'
                    
                    data = [
                        adata.obs.loc[human_mask, 'HDR_score'],
                        adata.obs.loc[mouse_mask, 'HDR_score']
                    ]
                    
                    ax.boxplot(data, labels=['Human', 'Mouse'])
                    ax.set_ylabel('HDR Score')
                    ax.set_title('HDR Expression by Species')
                    st.pyplot(fig)
    else:
        st.warning("No HDR genes found in this dataset")
else:
    st.error("HDR gene list not found")

# Cell Type Annotation
st.header("Cell Type Annotation")

if markers_df is not None and st.button("Run Cell Type Annotation"):
    with st.spinner("Annotating cell types..."):
        scores_df = annotate_cell_types(adata, markers_df)
        
        if scores_df is not None:
            # Show UMAP with cell types
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.umap(adata, color='cell_type_auto', legend_loc='on data',
                      title='Auto-assigned Cell Types', show=False, ax=ax)
            st.pyplot(fig)
            
            # Show confidence scores
            fig, ax = plt.subplots(figsize=(8, 6))
            adata.obs['confidence_score'].hist(bins=50, ax=ax)
            ax.set_xlabel('Confidence Score')
            ax.set_ylabel('Number of Cells')
            ax.set_title('Cell Type Assignment Confidence')
            st.pyplot(fig)
            
            # Cell type composition
            if sel == "Integrated (Batch Corrected)":
                cross_tab = pd.crosstab(adata.obs['cell_type_auto'], adata.obs['species'])
                st.write("Cell type composition by species:")
                st.dataframe(cross_tab)

if st.sidebar.checkbox("Advanced HDR Analysis"):
    plot_hdr_by_cluster(adata, hdr_genes)

# Quality Control
st.header("Quality Control")
if st.button("Show QC Plots"):
    if 'n_genes_by_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
    axes[0].set_title('Genes per Cell')
    
    sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
    axes[1].set_title('Total Counts per Cell')
    
    if 'pct_counts_mt' in adata.obs.columns:
        sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
        axes[2].set_title('Mitochondrial %')
    
    plt.tight_layout()
    st.pyplot(fig)

# Export options
st.header("Export Results")
if st.button("Export Current Analysis"):
    # Save the current state
    output_file = f"results/analysis_{sel.replace(' ', '_').replace('(', '').replace(')', '')}.h5ad"
    adata.write_h5ad(output_file)
    st.success(f"Analysis saved to {output_file}")
