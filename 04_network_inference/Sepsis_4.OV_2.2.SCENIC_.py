# -----------------------------------------------------------
# Script: Sepsis_4.OV_2.2.SCENIC_.py
# Module: 04_network_inference
#
# Description:
#   Supplemental SCENIC analysis and visualizations.
#   This script loads SCENIC outputs (regulons, AUC matrices),
#   integrates them into AnnData, and produces dotplots, heatmaps,
#   and UMAP overlays of regulon activity for downstream interpretation.
#
# Inputs:
#   - AnnData object with SCENIC_AUC or regulon data
#   - Regulon .csv files (if stored externally)
#
# Outputs:
#   - Updated AnnData with regulon annotations
#   - Dotplots / heatmaps of regulon activities
#   - Optional filtered regulon sets per cell type
#
# Dependencies:
#   scanpy, omicverse, matplotlib, seaborn, pandas
# -----------------------------------------------------------

import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# 1. Load Data (AnnData or AUC Matrix)
# -----------------------------------------------------------

# Example:
# adata = sc.read_h5ad("results/scenic/adata_scenic.h5ad")
# auc_matrix = pd.read_csv("results/scenic/regulon_auc.csv", index_col=0)

print("Loaded object for SCENIC visualization:", adata)

# If AUC is stored in a .csv, integrate manually:
# for col in auc_matrix.columns:
#     adata.obs[col] = auc_matrix.loc[adata.obs_names, col]

# Regulon names (example extraction)
if "SCENIC_AUC" in adata.obsm_keys():
    auc_df = pd.DataFrame(
        adata.obsm["SCENIC_AUC"],
        index=adata.obs_names
    )
else:
    print("Warning: SCENIC AUC matrix not found in AnnData.")

# -----------------------------------------------------------
# 2. Select Regulons for Visualization
# -----------------------------------------------------------

# Example: Select top variable regulons
# top_regulons = auc_df.var(axis=0).sort_values(ascending=False).head(20).index.tolist()

# Alternatively: use curated regulons based on your study
# curated_regulons = ["Regulon1", "Regulon2", "Regulon3"]

reg_list = auc_df.columns[:10].tolist()  # example: take first 10
print("Selected regulons:", reg_list)

# -----------------------------------------------------------
# 3. Dotplot of Regulon Activity by Cell Type
# -----------------------------------------------------------

if "majorclass" in adata.obs.columns:
    adata.obs["group"] = adata.obs["majorclass"]
else:
    adata.obs["group"] = adata.obs["seurat_clusters"]

mean_auc = (
    auc_df[reg_list]
    .join(adata.obs["group"])
    .groupby("group")
    .mean()
)

plt.figure(figsize=(10, 6))
sns.heatmap(mean_auc, cmap="viridis")
plt.title("Mean Regulon Activity Across Cell Types")
plt.xlabel("Regulon")
plt.ylabel("Cell Type")
# plt.savefig("results/scenic/regulon_activity_heatmap.pdf", dpi=300, bbox_inches="tight")
plt.close()

# -----------------------------------------------------------
# 4. UMAP Overlay of Key Regulons
# -----------------------------------------------------------

if "X_umap" not in adata.obsm:
    sc.pp.neighbors(adata, n_neighbors=20)
    sc.tl.umap(adata)

for reg in reg_list[:5]:  # show a few
    adata.obs[reg] = auc_df[reg]
    sc.pl.umap(adata, color=reg, cmap="viridis", show=False)
    plt.title(f"Regulon: {reg}")
    # plt.savefig(f"results/scenic/UMAP_regulon_{reg}.pdf", dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------------------------------------
# 5. Optional: Regulon Ranking / Filtering
# -----------------------------------------------------------

# Example: identify cell-type enriched regulons
group_means = mean_auc
enriched = {
    cell: group_means.loc[cell].sort_values(ascending=False).head(5).index.tolist()
    for cell in group_means.index
}

print("Top enriched regulons per group:")
for k, v in enriched.items():
    print(k, ":", v)

# -----------------------------------------------------------
# 6. Save Outputs
# -----------------------------------------------------------

# adata.write("results/scenic/adata_scenic_visualized.h5ad")
# mean_auc.to_csv("results/scenic/regulon_activity_summary.csv")

print("SCENIC supplemental analysis completed.")
