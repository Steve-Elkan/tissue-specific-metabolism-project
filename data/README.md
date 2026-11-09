
# Data Sources

This project uses two publicly available datasets.

## 1. Human Protein Atlas

Dataset: RNA tissue consensus

Source:
https://www.proteinatlas.org/download/

File:
rna_tissue_consensus_tissues.tsv.zip

Description:
Provides gene-expression information across human tissues.

The project uses this dataset to obtain expression profiles for liver
and skeletal muscle.

## 2. Recon3D

Dataset: Recon3D human genome-scale metabolic model

Source:
https://bigg.ucsd.edu/static/models/Recon3D.xml

Description:
A genome-scale reconstruction of human metabolism containing metabolic
reactions, metabolites, genes, and gene-protein-reaction associations.

## Downloading

Run:

python scripts/download_data.py

The datasets are downloaded into this directory automatically.

Raw datasets are not committed to the GitHub repository.
