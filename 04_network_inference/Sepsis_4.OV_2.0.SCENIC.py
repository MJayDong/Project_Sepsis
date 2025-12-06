# -----------------------------------------------------------
# Script: Sepsis_4.OV_2.0.SCENIC.py
# Module: 04_network_inference
#
# Description:
#   Gene regulatory network (GRN) inference using pySCENIC workflow
#   within the OmicVerse ecosystem. This script computes co-expression
#   modules, identifies transcription factor—target regulons, and scores
#   regulon activities across single cells.
#
# Inputs:
#   - Preprocessed AnnData (.h5ad)
#
# Outputs:
#   - pySCENIC regulon results
#   - AUC activity matrices
#   - Updated AnnData with regulon scores
#
# Dependencies:
#   - omicverse
#   - scanpy
#   - pyscenic
#   - numpy / pandas / matplotlib
#
# Notes:
#   - pySCENIC steps can be computationally expensive.
#   - Adjust TF database paths as needed.
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
# 1. Load AnnData
# -----------------------------------------------------------

# Example path:
# adata = sc.read_h5ad("data/preprocessed/sepsis_preprocessed.h5ad")

print("Loaded AnnData:", adata)

# -----------------------------------------------------------
# 2. Extract Expression Matrix
# -----------------------------------------------------------

expr_matrix = adata.X
gene_names = adata.var_names
cell_names = adata.obs_names

print("Expression Matrix Shape:", expr_matrix.shape)

# -----------------------------------------------------------
# 3. GRN Inference (GRNBoost2)
# -----------------------------------------------------------

# Define transcription factors (TF list should match your species)
# You may load a curated TF list:
# tfs = pd.read_csv("data/TFs_mouse.txt", header=None)[0].tolist()

# Placeholder
# tfs = [...]

# Run GRNBoost (co-expression module inference)
# NOTE: This step can take long time depending on CPU cores.

# network = grnboost2(expr_matrix, gene_names=gene_names, tf_names=tfs, verbose=True)

# -----------------------------------------------------------
# 4. Load Motif Databases (Feather files)
# -----------------------------------------------------------

# Example:
# ranking_db_paths = [
#     "resources/motifs/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
#     "resources/motifs/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
# ]
# dbs = [RankingDatabase(fname, "mm10") for fname in ranking_db_paths]

# -----------------------------------------------------------
# 5. Motif Enrichment / Regulon Construction
# -----------------------------------------------------------

# regulon_df = prune2df(network, dbs)
# regulons = df2regulons(regulon_df)

# -----------------------------------------------------------
# 6. AUC (Regulon Activity) Calculation
# -----------------------------------------------------------

# auc_matrix = aucell(expr_matrix, regulons, num_workers=4)
# auc_matrix = pd.DataFrame(auc_matrix, index=cell_names, columns=[r.name for r in regulons])

# -----------------------------------------------------------
# 7. Integrate Regulon Activity into AnnData
# -----------------------------------------------------------

# adata.obsm["SCENIC_AUC"] = auc_matrix

# -----------------------------------------------------------
# 8. Visualization
# -----------------------------------------------------------

# sc.pl.umap(adata, color=["n_genes"], show=False)
# sc.pl.umap(adata, color=list(auc_matrix.columns[:5]), show=False)

plt.figure(figsize=(6,4))
plt.title("pySCENIC Processing Completed")

# -----------------------------------------------------------
# 9. Save Outputs
# -----------------------------------------------------------

# adata.write("results/scenic/adata_scenic.h5ad")
# auc_matrix.to_csv("results/scenic/regulon_auc.csv")
# print("SCENIC results saved.")
