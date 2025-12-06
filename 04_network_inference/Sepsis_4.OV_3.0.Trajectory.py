# -----------------------------------------------------------
# Script: Sepsis_4.OV_3.0.Trajectory.py
# Module: 03_crosscell_analysis
#
# Description:
#   OmicVerse-based trajectory inference pipeline for sepsis
#   single-cell data. This script performs pseudotime ordering,
#   trajectory visualization, trend analysis of dynamic genes,
#   and density scatter plots along pseudotime.
#
# Inputs:
#   - AnnData (.h5ad) with:
#       * UMAP / PCA embeddings
#       * Cell-type annotation
#       * Preprocessed expression values
#
# Outputs:
#   - Pseudotime values (adata.obs["pseudotime"])
#   - Gene trend matrices
#   - UMAP, heatmap, and trajectory visualizations
#
# Dependencies:
#   scanpy, omicverse, numpy, pandas, matplotlib, seaborn
# -----------------------------------------------------------

import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------

# Example:
# adata = sc.read_h5ad("results/preprocessing/sepsis_preprocessed.h5ad")
print("Loaded AnnData:", adata)

# -----------------------------------------------------------
# 2. Select HVGs (if needed)
# -----------------------------------------------------------

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat",
    subset=True
)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# -----------------------------------------------------------
# 3. Build kNN Graph and Compute Diffusion Map
# -----------------------------------------------------------

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.diffmap(adata)

# Set diffusion map as OV input
adata.obsm["X_diffmap"] = adata.obsm["X_diffmap"]

# -----------------------------------------------------------
# 4. Trajectory Inference Using OmicVerse
# -----------------------------------------------------------

ov.pp.neighbors(adata, use_rep="X_diffmap")
ov.tl.pseudotime(adata)

# Pseudotime values stored in:
#   adata.obs["pseudotime"]

print("Pseudotime calculation completed.")

# -----------------------------------------------------------
# 5. UMAP + Pseudotime Visualization
# -----------------------------------------------------------

if "X_umap" not in adata.obsm:
    sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color=["pseudotime"],
    cmap="viridis",
    show=False
)
plt.title("UMAP Colored by Pseudotime")
plt.close()

# -----------------------------------------------------------
# 6. Identify Pseudotime-Dependent Genes
# -----------------------------------------------------------

# Example: differential trend genes
pt_genes = ov.tl.trend_genes(adata, n_genes=200)

# pt_genes: list of genes with dynamic expression

print("Top pseudotime-dependent genes:", pt_genes[:10])

# -----------------------------------------------------------
# 7. Plot Gene Expression Trends Along Pseudotime
# -----------------------------------------------------------

for gene in pt_genes[:20]:  # visualize first 20
    ov.pl.trend(
        adata,
        gene=gene,
        show=False
    )
    plt.title(f"Expression Trend: {gene}")
    plt.close()

# -----------------------------------------------------------
# 8. Heatmap of Top Dynamic Genes
# -----------------------------------------------------------

ov.pl.trend_heatmap(
    adata,
    genes=pt_genes[:300],
    show=False
)
plt.title("Pseudotime Gene Trend Heatmap")
plt.close()

# -----------------------------------------------------------
# 9. Density Scatter Along Pseudotime (Example)
# -----------------------------------------------------------

def plot_gene_density(gene):
    plt.figure(figsize=(6,4))
    x = adata.obs["pseudotime"]
    y = adata[:, gene].X.toarray().flatten()
    sns.kdeplot(x=x, y=y, cmap="viridis", fill=True, thresh=0.05)
    plt.xlabel("Pseudotime")
    plt.ylabel(gene)
    plt.title(f"Density Scatter of {gene}")
    plt.close()

# Example:
# plot_gene_density("Il1b")

# -----------------------------------------------------------
# 10. Save Outputs
# -----------------------------------------------------------

# adata.write("results/trajectory/adata_pseudotime.h5ad")
# pd.DataFrame({"gene": pt_genes}).to_csv(
#     "results/trajectory/pseudotime_genes.csv",
#     index=False
# )

print("Trajectory analysis completed.")
