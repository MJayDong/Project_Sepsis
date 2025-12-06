# -----------------------------------------------------------
# Script: Sepsis_2.SCP-2.2.Cardiomyocytes.R
# Module: 02_celltype_analysis
#
# Description:
#   Cardiomyocyte-specific subclustering and marker analysis pipeline.
#   This script extracts cardiomyocytes from the global dataset,
#   performs reclustering, identifies functional states/subtypes,
#   and generates UMAP and marker visualizations for downstream
#   trajectory and regulatory network analyses.
#
# Inputs:
#   - Cell-typed Seurat object with $majorclass annotations
#
# Outputs:
#   - Cardiomyocyte-only Seurat object
#   - Subcluster markers
#   - Subtype annotation (e.g., Myh6-high, Myh7-high, stressed CM)
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
## 2. Subset Cardiomyocytes
## -----------------------------------------------------------

# Example:
# SeuratObj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

cm <- subset(SeuratObj, subset = majorclass == "Cardiomyocytes")

## -----------------------------------------------------------
## 3. Re-normalization and Variable Feature Selection
## -----------------------------------------------------------

cm <- NormalizeData(cm)
cm <- FindVariableFeatures(cm, selection.method = "vst", nfeatures = 3000)
cm <- ScaleData(cm, features = rownames(cm))

## -----------------------------------------------------------
## 4. PCA, Neighbors, Reclustering
## -----------------------------------------------------------

cm <- RunPCA(cm)
cm <- FindNeighbors(cm, dims = 1:20)
cm <- FindClusters(cm, resolution = 0.5)
cm <- RunUMAP(cm, dims = 1:20)

p_cm <- DimPlot(cm, reduction = "umap", label = TRUE) +
        ggtitle("Cardiomyocyte Subclusters")
print(p_cm)

## -----------------------------------------------------------
## 5. Marker Genes for CM Subclusters
## -----------------------------------------------------------

cm_markers <- FindAllMarkers(
  cm,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# write.csv(cm_markers, "results/cardiomyocytes/cm_subcluster_markers.csv", row.names = FALSE)

## -----------------------------------------------------------
## 6. Cardiomyocyte Functional Subtype Annotation
## -----------------------------------------------------------

cm_reference <- list(
  CM_contraction = c("Tnnt2", "Actc1", "Myh6"),
  CM_stress = c("Nppb", "Gadd45g", "Atf3"),
  CM_Myh7_high = c("Myh7"),
  CM_inflammatory = c("Cxcl2", "Il1b")
)

cm$subcelltype <- "Unassigned"

for (type in names(cm_reference)) {
  idx <- WhichCells(cm, expression = rownames(cm) %in% cm_reference[[type]])
  cm$subcelltype[idx] <- type
}

p_cm_type <- DimPlot(cm, reduction = "umap", group.by = "subcelltype", label = TRUE) +
             ggtitle("Cardiomyocyte Functional States")
print(p_cm_type)

## -----------------------------------------------------------
## 7. Feature Visualizations (Optional)
## -----------------------------------------------------------

# FeaturePlot(cm, features = c("Myh6", "Myh7", "Nppb", "Tnnt2"))
# VlnPlot(cm, features = c("Myh7"), group.by = "subcelltype")

## -----------------------------------------------------------
## 8. Save Output
## -----------------------------------------------------------

# saveRDS(cm, file = "results/cardiomyocytes/SeuratObj_cardiomyocytes.rds")
