# -----------------------------------------------------------
# Script: Sepsis_4.OV_3.2.Trajectory_.py
# Module: 03_crosscell_analysis
#
# Description:
#   Supplemental pseudotime trajectory analysis for sepsis single-cell
#   data using OmicVerse. This script provides additional trajectory
#   visualization, condition-wise pseudotime mapping, gene dynamics,
#   and trend-based comparisons across sepsis stages.
#
# Inputs:
#   - AnnData containing:
#       * pseudotime (adata.obs["pseudotime"])
#       * cell type or cluster annotations
#       * normalized expression values
#
# Outputs:
#   - Additional pseudotime plots
#   - Condition-separated pseudotime distributions
#   - Trend curves and heatmaps for key genes
#   - UMAP/pseudotime overlays
#
# Dependencies:
#   omicverse, scanpy, numpy, pandas, seaborn, matplotlib
# -----------------------------------------------------------

import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# 1. Ensure Required Fields Exist
# -----------------------------------------------------------

if "pseudotime" not in adata.obs:
    raise ValueError("Pseudotime not found. Run ov.tl.pseudotime() first.")

if "condition" not in adata.obs:
    print("Warning: condition label not found, using 'seurat_clusters' as group.")
    adata.obs["condition"] = adata.obs.get("seurat_clusters", "NA")

print("Pseudotime analysis on:", adata)

# -----------------------------------------------------------
# 2. UMAP Colored by Pseudotime
# -----------------------------------------------------------

sc.pl.umap(
    adata,
    color=["pseudotime"],
    cmap="viridis",
    show=False
)
plt.title("UMAP Colored by Pseudotime")
plt.close()

# -----------------------------------------------------------
# 3. Pseudotime Distribution by Condition
# -----------------------------------------------------------

plt.figure(figsize=(7,4))
sns.boxplot(
    data=adata.obs,
    x="condition",
    y="pseudotime",
    palette="viridis"
)
plt.title("Pseudotime Distribution Across Conditions")
plt.xticks(rotation=30)
# plt.savefig("results/trajectory/pseudotime_by_condition.pdf", dpi=300)
plt.close()

# -----------------------------------------------------------
# 4. Density Plot by Condition
# -----------------------------------------------------------

plt.figure(figsize=(7,4))
sns.kdeplot(
    data=adata.obs,
    x="pseudotime",
    hue="condition",
    common_norm=False
)
plt.title("Pseudotime Density by Condition")
# plt.savefig("results/trajectory/pseudotime_density_condition.pdf", dpi=300)
plt.close()

# -----------------------------------------------------------
# 5. Dynamic Gene Detection (Advanced Supplemental)
# -----------------------------------------------------------

dynamic_genes = ov.tl.trend_genes(adata, n_genes=400)
print("Identified dynamic genes:", dynamic_genes[:10])

# -----------------------------------------------------------
# 6. Gene Trend Visualization (Top genes)
# -----------------------------------------------------------

top_genes_to_plot = dynamic_genes[:20]

for gene in top_genes_to_plot:
    ov.pl.trend(
        adata,
        gene=gene,
        show=False
    )
    plt.title(f"Expression Trend Along Pseudotime: {gene}")
    # plt.savefig(f"results/trajectory/gene_trend_{gene}.pdf", dpi=300)
    plt.close()

# -----------------------------------------------------------
# 7. Heatmap of Dynamic Genes
# -----------------------------------------------------------

ov.pl.trend_heatmap(
    adata,
    genes=dynamic_genes[:300],
    show=False
)
plt.title("Heatmap of Top Dynamic Genes Along Pseudotime")
plt.close()

# -----------------------------------------------------------
# 8. Branch / Local Pattern Exploration (If branch exists)
# -----------------------------------------------------------

if "branch_id" in adata.obs:
    plt.figure(figsize=(6,4))
    sns.scatterplot(
        x=adata.obs["pseudotime"],
        y=adata.obs["branch_id"],
        hue=adata.obs["condition"],
        palette="viridis"
    )
    plt.title("Branch Assignment Across Pseudotime")
    plt.close()

# -----------------------------------------------------------
# 9. Save Results
# -----------------------------------------------------------

# adata.write("results/trajectory/adata_trajectory_supplemental.h5ad")
# pd.DataFrame({"dynamic_gene": dynamic_genes}).to_csv(
#     "results/trajectory/dynamic_genes_supplemental.csv",
#     index=False
# )

print("Supplemental trajectory analysis completed.")
