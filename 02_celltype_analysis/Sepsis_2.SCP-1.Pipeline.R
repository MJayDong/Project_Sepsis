# -----------------------------------------------------------
# Script: Sepsis_2.SCP-1.Pipeline.R
# Module: 02_celltype_analysis
#
# Description:
#   Global cell-type annotation pipeline for the sepsis scRNA-seq dataset.
#   This script performs secondary QC, dimensionality reduction, clustering,
#   marker-based cell type assignment, and the generation of core plots used
#   to guide downstream immune/cardiomyocyte/endothelial/etc. analyses.
#
# Inputs:
#   - Integrated Seurat object from preprocessing pipeline
#
# Outputs:
#   - Annotated Seurat object with major cell types
#   - Marker gene lists
#   - UMAP plots and feature plots for cell-type confirmation
#
# Dependencies:
#   Seurat (>=4.0), dplyr, ggplot2, patchwork
# -----------------------------------------------------------

## -----------------------------------------------------------
## 1. Load Libraries
## -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## -----------------------------------------------------------
## 2. Load Preprocessed Object
## -----------------------------------------------------------
# Replace with relative path for GitHub compatibility
# SeuratObj <- readRDS("results/preprocessing/SeuratObj_preprocessed.rds")

# Placeholder if object loaded externally:
# SeuratObj <- <your_preprocessed_object>

## -----------------------------------------------------------
## 3. Re-run PCA / Neighbors / Clustering (if needed)
## -----------------------------------------------------------

# This ensures cell-type clustering is stable after preprocessing
SeuratObj <- RunPCA(SeuratObj, verbose = FALSE)
SeuratObj <- FindNeighbors(SeuratObj, dims = 1:30)
SeuratObj <- FindClusters(SeuratObj, resolution = 0.6)
SeuratObj <- RunUMAP(SeuratObj, dims = 1:30)

## -----------------------------------------------------------
## 4. Identify Marker Genes for Clusters
## -----------------------------------------------------------

DefaultAssay(SeuratObj) <- "RNA"

markers <- FindAllMarkers(
  SeuratObj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save marker table
# write.csv(markers, file = "results/celltype/cluster_markers.csv", row.names = FALSE)

## -----------------------------------------------------------
## 5. Assign Major Cell Types
## -----------------------------------------------------------

# Define reference marker sets (example for mouse heart)
marker_list <- list(
  T_cells = c("Cd3d", "Cd3e", "Cd4", "Cd8a"),
  B_cells = c("Cd79a", "Cd74"),
  NK_cells = c("Nkg7", "Gzmb"),
  Macrophages = c("Lyz2", "Cd68", "Lpl", "Mrc1"),
  Neutrophils = c("S100a8", "S100a9", "Ly6g"),
  Endothelium = c("Pecam1", "Kdr", "Cdh5"),
  Fibroblasts = c("Col1a1", "Dcn", "Pdgfra"),
  Cardiomyocytes = c("Tnnt2", "Actc1", "Myh7")
)

# Assign cell types based on dominant markers
SeuratObj$majorclass <- NA

for (celltype in names(marker_list)) {
  hits <- rownames(SeuratObj)[rownames(SeuratObj) %in% marker_list[[celltype]]]
  if (length(hits) > 0) {
    idx <- WhichCells(SeuratObj, expression = rownames(SeuratObj) %in% marker_list[[celltype]])
    SeuratObj$majorclass[idx] <- celltype
  }
}

## -----------------------------------------------------------
## 6. Visualization
## -----------------------------------------------------------

p_umap <- DimPlot(
  SeuratObj,
  reduction = "umap",
  group.by = "majorclass",
  label = TRUE
) + ggtitle("Major Cell Type Annotation")

print(p_umap)

# FeaturePlot examples (adjust markers as needed)
# FeaturePlot(SeuratObj, features = c("Tnnt2", "Mrc1", "Col1a1", "Pecam1"))

## -----------------------------------------------------------
## 7. Save Outputs
## -----------------------------------------------------------

# saveRDS(SeuratObj, file = "results/celltype/SeuratObj_celltyped.rds")
