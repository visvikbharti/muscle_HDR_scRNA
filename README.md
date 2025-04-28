# muscle_HDR_scRNA

Single-cell RNA-seq analysis pipeline for examining Homology-Directed Repair (HDR) gene expression in muscle cells

## Overview

This repository contains a computational pipeline for analyzing single-cell RNA sequencing data from muscle tissue, with a focus on Homology-Directed Repair (HDR) gene expression patterns. The pipeline integrates data from human and mouse datasets, performs quality control, clustering, cell type annotation, and specialized analysis of HDR gene expression across different cell populations.

## Features

- Automated download of datasets from GEO
- Quality control and preprocessing of scRNA-seq data
- Cross-species integration using scVI
- Cell type annotation based on marker genes
- Specialized analysis of HDR gene expression
- Interactive visualization with Streamlit dashboard
- Automated report generation

## Requirements

### Environment

The pipeline requires Python 3.8+ and the following key packages:
- scanpy
- anndata
- scvi-tools
- pandas
- numpy
- matplotlib
- seaborn
- streamlit
- snakemake

A full environment can be created using the provided configuration:

```bash
conda env create -f env/environment.yml
conda activate muscle_scrna
```

## Pipeline Structure

The pipeline is implemented using Snakemake and consists of the following main components:

1. **Data Acquisition**: Downloads raw data from GEO
2. **Quality Control**: Filters cells and genes based on quality metrics
3. **Preprocessing**: Normalization and feature selection
4. **Integration**: Cross-species integration using scVI
5. **Clustering**: Identification of cell populations
6. **Cell Type Annotation**: Automatic annotation based on marker genes
7. **HDR Analysis**: Focused analysis of HDR genes across cell types
8. **Reporting**: Generation of comprehensive analysis reports

## Getting Started

### Quick Start

1. Clone the repository:
```bash
git clone https://github.com/yourusername/muscle_HDR_scRNA.git
cd muscle_HDR_scRNA
```

2. Create the conda environment:
```bash
conda env create -f env/environment.yml
conda activate muscle_scrna
```

3. Run the full pipeline:
```bash
snakemake --cores all
```

4. Launch the interactive dashboard:
```bash
streamlit run app.py
```

### Configuration

The pipeline is configured via the `config.yaml` file, which contains parameters for:
- Quality control thresholds
- Processing parameters
- Clustering settings
- Integration methods
- File paths
- Species information
- Report settings

Example configuration:
```yaml
qc:
  min_genes: 300
  min_cells: 3
  pct_mt_max: 15

processing:
  n_top_genes: 2000
  normalize_target_sum: 10000
  scale_max_value: 10

clustering:
  n_neighbors: 15
  n_pcs: 20
  leiden_resolution: 0.5
  method: 'gauss'
```

## Key Components

### Scripts

- `csv_to_h5ad.py`: Converts CSV matrices to AnnData format
- `qc_to_h5ad.py`: Processes 10x data and performs QC
- `get_orthologs.py`: Retrieves ortholog mappings between human and mouse genes
- `integrate_scvi.py`: Integrates datasets using scVI
- `build_report.py`: Generates comprehensive analysis reports

### Data Files

- `HDR_477_genes.txt`: List of HDR-related genes for focused analysis
- `markers.tsv`: Cell type marker genes used for annotation

### Rules

- `download.smk`: Rules for downloading raw data
- `qc_scanpy.smk`: Rules for quality control and processing

## Interactive Dashboard

The repository includes a Streamlit app (`app.py`) that provides interactive visualization of the analysis results, including:
- UMAP visualization of cell clusters
- Cell type annotation
- HDR gene expression analysis
- Quality control metrics
- Export options for analysis results

## Report Generation

The pipeline automatically generates a comprehensive report of the analysis using the `build_report.py` script. The report includes:
- Data provenance
- Processing parameters
- Cluster summaries
- HDR gene expression analysis
- Cross-species comparison
- Visualization of key findings

## Future Enhancements

Based on our development roadmap, the following improvements are planned:
- Dynamic Leiden resolution optimization
- Enhanced batch correction
- RNA velocity analysis
- Expanded marker panel for cell type annotation
- Improved UX for the Streamlit dashboard
- Advanced integration of cross-species data

## Contributing

Contributions to this pipeline are welcome. Please feel free to submit issues or pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite:

Rubenstein, A. B., Smith, G. R., Raue, U., Begue, G., Minchev, K., Ruf-Zamojski, F., Nair, V. D., Wang, X., Zhou, L., Zaslavsky, E., Trappe, T. A., Trappe, S., & Sealfon, S. C. (2020). Single-cell transcriptional profiles in human skeletal muscle. Scientific Reports, 10, 229. https://doi.org/10.1038/s41598-019-57110-6

## Acknowledgments

This project builds upon single-cell RNA sequencing datasets from human and mouse skeletal muscle tissues, as described in Rubenstein et al. (2020). We thank the authors for making their data publicly available through the GEO repository (GSE130646, GSE130977, and GSE138707).
