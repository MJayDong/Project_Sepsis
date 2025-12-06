# -----------------------------------------------------------
# Script: Sepsis_4.OV_2.1.SCENIC_Macro-Mrc1.py
# Module: 04_network_inference
#
# Description:
#   SCENIC analysis focused on Mrc1+ macrophage subpopulation.
#   This script subsets Mrc1+ macrophages from the global AnnData
#   object, performs GRN inference and regulon activity scoring,
#   and generates macrophage-specific regulon visualizations.
#
# Inputs:
#   - Preprocessed AnnData (.h5ad) with:
#       * cell type / subcelltype annotation (including Mrc1+ macrophages)
#       * normalized expression values
#
# Outputs:
#   - Regulon definitions for Mrc1+ macrophages
#   - Regulon activity (AUC) matrix
#   - Updated AnnData subset with regulon scores
#
# Dependencies:
#   - scanpy
#   - omicverse
#   - pyscenic
#   - numpy, pandas, matplotlib
# -----------------------------------------------------------

import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyscenic.aucell import aucell
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.prune import prune2df, df2regulons
from arboreto.algo import grnboost2

# -----------------------------------------------------------
# 1. Load Global AnnData and Subset Mrc1+ Macrophages
# -----------------------------------------------------------

# Example: global object with subcelltype or SubCellType column
# adata = sc.read_h5ad("data/preprocessed/sepsis_preprocessed.h5ad")

# 这里假设:
#   - adata.obs["SubCellType"] 包含 "Macrophages-Mrc1" 或类似标签
#   - 你可以根据实际标签调整筛选条件
macro_mrc1 = adata[adata.obs["SubCellType"] == "Macrophages-Mrc1"].copy()
print("Mrc1+ macrophage subset:", macro_mrc1)

# -----------------------------------------------------------
# 2. Basic Preprocessing for Subset
# -----------------------------------------------------------

# 可选：重复 HVG 选择与 PCA，仅在子集上工作
sc.pp.highly_variable_genes(
    macro_mrc1,
    n_top_genes=3000,
    flavor="seurat",
    subset=True
)

sc.pp.scale(macro_mrc1, max_value=10)
sc.tl.pca(macro_mrc1, n_comps=50)

# -----------------------------------------------------------
# 3. Expression Matrix and TF List
# -----------------------------------------------------------

expr_matrix = macro_mrc1.X
gene_names = macro_mrc1.var_names
cell_names = macro_mrc1.obs_names

print("Subset expression matrix shape:", expr_matrix.shape)

# Example TF list:
# tfs = pd.read_csv("data/TFs_mouse.txt", header=None)[0].tolist()
# tfs = [...]

# -----------------------------------------------------------
# 4. GRN Inference (GRNBoost2)
# -----------------------------------------------------------

# network = grnboost2(
#     expr_matrix,
#     gene_names=gene_names,
#     tf_names=tfs,
#     verbose=True
# )

# -----------------------------------------------------------
# 5. Motif Databases and Regulon Construction
# -----------------------------------------------------------

# ranking_db_paths = [
#     "resources/motifs/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
#     "resources/motifs/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
# ]
# dbs = [RankingDatabase(fname, "mm10") for fname in ranking_db_paths]

# regulon_df = prune2df(network, dbs)
# regulons = df2regulons(regulon_df)

# -----------------------------------------------------------
# 6. AUC Calculation (Regulon Activity)
# -----------------------------------------------------------

# auc_matrix = aucell(expr_matrix, regulons, num_workers=4)
# auc_matrix = pd.DataFrame(
#     auc_matrix,
#     index=cell_names,
#     columns=[r.name for r in regulons]
# )

# macro_mrc1.obsm["SCENIC_AUC_Mrc1"] = auc_matrix

# -----------------------------------------------------------
# 7. Visualization: UMAP + Selected Regulons
# -----------------------------------------------------------

# sc.pp.neighbors(macro_mrc1, n_neighbors=20, n_pcs=30)
# sc.tl.umap(macro_mrc1)

# 示例：绘制前若干 regulon 的 AUC
# top_regulons = auc_matrix.var(axis=0).sort_values(ascending=False).head(5).index.tolist()
# for reg in top_regulons:
#     macro_mrc1.obs[reg] = auc_matrix[reg]
#     sc.pl.umap(macro_mrc1, color=reg, cmap="viridis", show=False)
#     plt.title(f"Regulon activity: {reg}")
#     # plt.savefig(f"results/scenic/Mrc1/UMAP_{reg}.pdf", dpi=300, bbox_inches="tight")
#     plt.close()

# -----------------------------------------------------------
# 8. Save Subset + Regulon Results
# -----------------------------------------------------------

# macro_mrc1.write("results/scenic/Mrc1/adata_Mrc1_SCENIC.h5ad")
# auc_matrix.to_csv("results/scenic/Mrc1/Mrc1_regulon_auc.csv")
# print("Mrc1+ macrophage SCENIC results saved.")
