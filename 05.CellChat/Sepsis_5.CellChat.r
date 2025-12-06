# -----------------------------------------------------------
# Script: Sepsis_5.CellChat.R
# Module: 04_network_inference
#
# Description:
#   Cell–cell communication inference using CellChat for
#   sepsis-induced cardiac single-cell transcriptomes. This script
#   constructs ligand–receptor networks, identifies enriched
#   communication pathways, computes outgoing/incoming signaling
#   strengths, and produces visualizations of intercellular signaling.
#
# Inputs:
#   - Seurat object with:
#       * normalized expression data
#       * cell-type annotations (majorclass or subcelltype)
#       * condition labels (optional)
#
# Outputs:
#   - CellChat object with inferred communication networks
#   - Network visualizations (circle, heatmap, bubble)
#   - Signaling pathway enrichment charts
#
# Dependencies:
#   CellChat, Seurat, patchwork, tidyverse
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(tidyverse)
})

# -----------------------------------------------------------
# 1. Prepare Input Data
# -----------------------------------------------------------

# Example:
# seurat_obj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

data_input <- GetAssayData(seurat_obj, slot = "data")
meta_df <- seurat_obj@meta.data

# CellChat requires:
#   - "labels": cell type annotation
#   - rownames(meta_df) == colnames(data_input)

meta_df$labels <- meta_df$majorclass %||% meta_df$subcelltype

cellchat <- createCellChat(
  object = data_input,
  meta = meta_df,
  group.by = "labels"
)

# -----------------------------------------------------------
# 2. Load Mouse Ligand–Receptor Database
# -----------------------------------------------------------

cellchat@DB <- CellChatDB.mouse
# To use secreted signaling only:
# cellchat@DB <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")

# -----------------------------------------------------------
# 3. Preprocessing for CellChat
# -----------------------------------------------------------

cellchat <- subsetData(cellchat)     # subset LR pairs
future::plan("multisession")         # parallel computation

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project expression data into PPI structure
cellchat <- projectData(cellchat, PPI.mouse)

# -----------------------------------------------------------
# 4. Infer Communication Networks
# -----------------------------------------------------------

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Pathway-level communication probabilities
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)

# -----------------------------------------------------------
# 5. Visualizations
# -----------------------------------------------------------

# Circle plot for total signaling
net <- cellchat@net$count
group_size <- as.numeric(table(cellchat@idents))

# circle plot
# netVisual_circle(cellchat@net$count, vertex.weight = group_size)

# Heatmap of interactions
# netVisual_heatmap(cellchat)

# Bubble plot for ligand–receptor pairs
# netVisual_bubble(cellchat, sources.use = 1, targets.use = 2)

# Pathway activity
# netAnalysis_signalingRole_network(cellchat, slot.name = "netP")

# -----------------------------------------------------------
# 6. Centrality Analysis
# -----------------------------------------------------------

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Identify dominant senders/receivers
# netAnalysis_signalingRole_scatter(cellchat)
# netAnalysis_signalingRole_heatmap(cellchat)

# -----------------------------------------------------------
# 7. Save Outputs
# -----------------------------------------------------------

# saveRDS(cellchat, file = "results/cellchat/CellChat_object.rds")

print("CellChat analysis completed.")
