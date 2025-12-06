# -----------------------------------------------------------
# Script: Sepsis_4.OV_1.Preprocessing.py
# Module: 01_preprocessing
#
# Description:
#   Preprocessing pipeline for OmicVerse-based single-cell analysis.
#   Loads AnnData object, performs quality control, normalization,
#   highly variable gene selection, dimensionality reduction, and
#   neighborhood graph construction.
#
# Inputs:
#   - Raw or partially processed AnnData (.h5ad) file
#
# Outputs:
#   - Preprocessed AnnData object saved for downstream analyses
#
# Dependencies:
#   - Python >= 3.8
#   - omicverse
#   - scanpy
#   - anndata
#   - numpy / pandas / matplotlib
#
# Notes:
#   - All file paths should be adapted to project-relative paths.
# -----------------------------------------------------------

import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------

# Replace with your relative project path
# Example:
# adata = sc.read_h5ad("data/raw/sepsis_raw.h5ad")

# Placeholder if loaded outside
# adata = ...

print("Initial AnnData:", adata)

# -----------------------------------------------------------
# 2. Basic Filtering (Optional)
# -----------------------------------------------------------

# You can adjust min_genes / min_cells depending on your dataset.
# Uncomment if needed:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# -----------------------------------------------------------
# 3. Normalize and Log Transform
# -----------------------------------------------------------

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# -----------------------------------------------------------
# 4. Highly Variable Genes
# -----------------------------------------------------------

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat",
    subset=True
)

print("Number of HVGs:", np.sum(adata.var["highly_variable"]))

# -----------------------------------------------------------
# 5. Scaling and PCA
# -----------------------------------------------------------

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# -----------------------------------------------------------
# 6. Neighborhood Graph
# -----------------------------------------------------------

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

# -----------------------------------------------------------
# 7. UMAP Embedding
# -----------------------------------------------------------

sc.tl.umap(adata)
sc.pl.umap(adata, color=["n_genes", "total_counts"], show=False)
plt.title("UMAP of Preprocessed Cells")

# -----------------------------------------------------------
# 8. Save Output
# -----------------------------------------------------------

# save_path = "results/preprocessing/sepsis_preprocessed.h5ad"
# adata.write(save_path)
# print("Saved preprocessed AnnData to:", save_path)
