# -----------------------------------------------------------
# Script: Sepsis_2.SCP-2.3.Mesenchymal_Vascular.R
# Module: 02_celltype_analysis
#
# Description:
#   Subclustering and annotation pipeline for mesenchymal and
#   vascular-related populations in sepsis scRNA-seq data.
#   Includes fibroblasts, pericytes, smooth muscle cells (SMCs),
#   and other stromal cell states found in cardiac tissue.
#
# Inputs:
#   - Global Seurat object with $majorclass assigned
#
# Outputs:
#   - Mesenchymal + vascular subset Seurat object
#   - Subcluster markers
#   - Functional subtype annotations (e.g., ECM fibroblasts)
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
## 2. Subset Mesenchymal + Vascular Populations
## -----------------------------------------------------------

# Example:
# SeuratObj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

mv <- subset(
  SeuratObj,
  subset = majorclass %in% c("Fibroblasts", "Mesenchymal", "Pericytes", "Vascular_related")
)

## -----------------------------------------------------------
## 3. Re-normalization and Feature Selection
## -----------------------------------------------------------

mv <- NormalizeData(mv)
mv <- FindVariableFeatures(mv, selection.method = "vst", nfeatures = 3000)
mv <- ScaleData(mv, features = rownames(mv))

## -----------------------------------------------------------
## 4. PCA, Neighbors, Reclustering
## -----------------------------------------------------------

mv <- RunPCA(mv)
mv <- FindNeighbors(mv, dims = 1:20)
mv <- FindClusters(mv, resolution = 0.5)
mv <- RunUMAP(mv, dims = 1:20)

p_mv <- DimPlot(mv, reduction = "umap", label = TRUE) +
        ggtitle("Mesenchymal + Vascular Subclusters")
print(p_mv)

## -----------------------------------------------------------
## 5. Identify Subcluster Marker Genes
## -----------------------------------------------------------

mv_markers <- FindAllMarkers(
  mv,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# write.csv(mv_markers, "results/mesenchymal_vascular/mv_subcluster_markers.csv", row.names = FALSE)

## -----------------------------------------------------------
## 6. Functional Subtype Annotation
## -----------------------------------------------------------

mv_reference <- list(
  Fibroblast_ECM = c("Col1a1", "Col3a1", "Dcn"),
  Fibroblast_progenitor = c("Pdgfra", "Ly6a"),
  Fibroblast_stress = c("Mt1", "Mt2"),
  Pericytes = c("Rgs5", "Pdgfrb"),
  Smooth_muscle_cells = c("Acta2", "Tagln", "Myh11"),
  EndMT_like = c("Klf4", "Cdh2", "Fn1")
)

mv$subcelltype <- "Unassigned"

for (type in names(mv_reference)) {
  genes <- mv_reference[[type]]
  idx <- WhichCells(mv, expression = rownames(mv) %in% genes)
  mv$subcelltype[idx] <- type
}

p_mv_type <- DimPlot(mv, reduction = "umap", group.by = "subcelltype", label = TRUE) +
             ggtitle("Mesenchymal / Vascular Functional States")
print(p_mv_type)

## -----------------------------------------------------------
## 7. Visualization (Optional)
## -----------------------------------------------------------

# FeaturePlot(mv, features = c("Pdgfra", "Rgs5", "Myh11"))
# VlnPlot(mv, features = c("Col1a1", "Acta2"), group.by = "subcelltype")

## -----------------------------------------------------------
## 8. Save Output
## -----------------------------------------------------------

# saveRDS(mv, file = "results/mesenchymal_vascular/SeuratObj_mesenchymal_vascular.rds")
