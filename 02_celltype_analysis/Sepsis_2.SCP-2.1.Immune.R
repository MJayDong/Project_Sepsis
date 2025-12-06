# -----------------------------------------------------------
# Script: Sepsis_2.SCP-2.1.Immune.R
# Module: 02_celltype_analysis
#
# Description:
#   Immune cell subclustering pipeline for the sepsis scRNA-seq dataset.
#   This script extracts immune populations from the global Seurat object,
#   performs reclustering, identifies subpopulation markers, and generates
#   core UMAP and marker-based visualizations used in downstream analyses.
#
# Inputs:
#   - Global annotated Seurat object (SeuratObj) with $majorclass assigned
#
# Outputs:
#   - Immune-only Seurat object
#   - Subcluster-level UMAPs and markers
#   - Marker table for all immune subclusters
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
## 2. Load Global Object and Subset Immune Population
## -----------------------------------------------------------
# Example:
# SeuratObj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

immune_cells <- subset(SeuratObj, subset = majorclass %in% c("T_cells", "B_cells", "NK_cells", "Macrophages", "Neutrophils"))

## -----------------------------------------------------------
## 3. Re-normalization for Immune-Specific Analysis
## -----------------------------------------------------------

immune_cells <- NormalizeData(immune_cells)
immune_cells <- FindVariableFeatures(immune_cells, selection.method = "vst", nfeatures = 3000)
immune_cells <- ScaleData(immune_cells)

## -----------------------------------------------------------
## 4. Dimensionality Reduction and Re-clustering
## -----------------------------------------------------------

immune_cells <- RunPCA(immune_cells)
immune_cells <- FindNeighbors(immune_cells, dims = 1:30)
immune_cells <- FindClusters(immune_cells, resolution = 0.6)
immune_cells <- RunUMAP(immune_cells, dims = 1:30)

p_umap <- DimPlot(immune_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
          ggtitle("Immune Subclusters")
print(p_umap)

## -----------------------------------------------------------
## 5. Marker Genes for Immune Subclusters
## -----------------------------------------------------------

immune_markers <- FindAllMarkers(
  immune_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# write.csv(immune_markers, file = "results/immune/immune_subcluster_markers.csv", row.names = FALSE)

## -----------------------------------------------------------
## 6. Annotate Immune Subtypes (Example Framework)
## -----------------------------------------------------------

# Reference signatures (example for mouse inflammation)
immune_reference <- list(
  T_naive = c("Lef1", "Tcf7"),
  T_effectors = c("Gzmb", "Ifng"),
  Treg = c("Foxp3", "Il2ra"),
  B_cells = c("Cd79a", "Ms4a1"),
  NK_cells = c("Nkg7", "Gzma"),
  Macrophages_M1 = c("Il1b", "Tnf", "Nos2"),
  Macrophages_M2 = c("Mrc1", "Folr2", "C1qa"),
  Neutrophils = c("S100a8", "S100a9")
)

immune_cells$subcelltype <- "Unassigned"

for (subtype in names(immune_reference)) {
  genes <- immune_reference[[subtype]]
  idx <- WhichCells(immune_cells, expression = rownames(immune_cells) %in% genes)
  immune_cells$subcelltype[idx] <- subtype
}

## -----------------------------------------------------------
## 7. Visualization of Immune Subtypes
## -----------------------------------------------------------

p_immune <- DimPlot(immune_cells, reduction = "umap", label = TRUE, group.by = "subcelltype") +
            ggtitle("Immune Subtype Annotation")
print(p_immune)

# Example feature visualization
# FeaturePlot(immune_cells, features = c("Cd3d", "Gzmb", "Mrc1", "S100a8"))

## -----------------------------------------------------------
## 8. Save Immune Object
## -----------------------------------------------------------

# saveRDS(immune_cells, file = "results/immune/SeuratObj_immune.rds")
