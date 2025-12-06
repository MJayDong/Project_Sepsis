# -----------------------------------------------------------
# Script: Sepsis_4.OV_3.1.Trajectory_Macro-Mrc1.py
# Module: 03_crosscell_analysis
#
# Description:
#   Pseudotime trajectory analysis focused on Mrc1+ macrophages.
#   This script extracts Mrc1+ macrophage subsets, performs
#   OmicVerse-based pseudotime reconstruction, identifies
#   pseudotime-dependent genes, and visualizes transcriptional
#   dynamics unique to this macrophage state.
#
# Inputs:
#   - AnnData (.h5ad) containing SubCellType / macrophage annotation
#
# Outputs:
#   - Pseudotime ordering for Mrc1+ macrophages
#   - Dynamic gene list
#   - Plots (UMAP, trend heatmaps, trajectory curves)
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
# 1. Load Global AnnData and Subset Mrc1+ Macrophages
# -----------------------------------------------------------

# Example:
# adata = sc.read_h5ad("results/preprocessing/sepsis_preprocessed.h5ad")

macro = adata[adata.obs["SubCellType"] == "Macrophages-Mrc1"].copy()
print("Selected Mrc1+ macrophages:", macro)

# -----------------------------------------------------------
# 2. HVG + PCA on Subset
# -----------------------------------------------------------

sc.pp.highly_variable_genes(
    macro,
    n_top_genes=3000,
    flavor="seurat",
    subset=True
)

sc.pp.scale(macro, max_value=10)
sc.tl.pca(macro, n_comps=50)

# -----------------------------------------------------------
# 3. Diffusion Map + OV Neighbors
# -----------------------------------------------------------

sc.pp.neighbors(macro, n_neighbors=30, n_pcs=30)
sc.tl.diffmap(macro)

macro.obsm["X_diffmap"] = macro.obsm["X_diffmap"]

ov.pp.neighbors(macro, use_rep="X_diffmap")
ov.tl.pseudotime(macro)

print("Pseudotime range:", macro.obs["pseudotime"].min(), "→", macro.obs["pseudotime"].max())

# -----------------------------------------------------------
# 4. Visualization: UMAP + Pseudotime
# -----------------------------------------------------------

if "X_umap" not in macro.obsm:
    sc.tl.umap(macro)

sc.pl.umap(
    macro,
    color=["pseudotime"],
    cmap="viridis",
    show=False
)
plt.title("UMAP Colored by Pseudotime (Mrc1+ Macrophages)")
plt.close()

# -----------------------------------------------------------
# 5. Identify Dynamic Genes Along Pseudotime
# -----------------------------------------------------------

pt_genes = ov.tl.trend_genes(macro, n_genes=300)
print("Top dynamic genes:", pt_genes[:10])

# -----------------------------------------------------------
# 6. Plot Gene Trends Along Pseudotime
# -----------------------------------------------------------

for gene in pt_genes[:15]:
    ov.pl.trend(
        macro,
        gene=gene,
        show=False
    )
    plt.title(f"Expression Trend: {gene}")
    plt.close()

# -----------------------------------------------------------
# 7. Trend Heatmap
# -----------------------------------------------------------

ov.pl.trend_heatmap(
    macro,
    genes=pt_genes[:200],
    show=False
)
plt.title("Macrophage Pseudotime Gene Trend Heatmap")
plt.close()

# -----------------------------------------------------------
# 8. Optional: Gene Density KDE Along Pseudotime
# -----------------------------------------------------------

def density_along_pseudotime(adata, gene):
    plt.figure(figsize=(6,4))
    x = adata.obs["pseudotime"]
    y = adata[:, gene].X.toarray().flatten()
    sns.kdeplot(
        x=x, y=y,
        cmap="viridis", fill=True, thresh=0.05
    )
    plt.xlabel("Pseudotime")
    plt.ylabel(gene)
    plt.title(f"Pseudotime Density: {gene}")
    plt.close()

# Example:
# density_along_pseudotime(macro, "Il1b")

# -----------------------------------------------------------
# 9. Save Output
# -----------------------------------------------------------

# macro.write("results/trajectory/Mrc1/adata_Mrc1_pseudotime.h5ad")
# pd.DataFrame({"gene": pt_genes}).to_csv(
#     "results/trajectory/Mrc1/Mrc1_pseudotime_genes.csv",
#     index=False
# )

print("Mrc1+ macrophage trajectory analysis completed.")
