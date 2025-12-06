# -----------------------------------------------------------
# Script: Sepsis_2.SCP-2.4.Endothelial_System.R
# Module: 02_celltype_analysis
#
# Description:
#   Endothelial-cell subclustering and functional-state annotation
#   in sepsis-induced cardiac tissue. Includes identification of
#   microvascular, arterial-like, venous-like, and lymphatic EC states,
#   with marker-based functional characterization.
#
# Inputs:
#   - Global Seurat object with $majorclass annotations
#
# Outputs:
#   - Endothelial-only Seurat object
#   - Endothelial subclusters and markers
#   - Functional subtype assignments (capillary / artery / vein / lymphatic EC)
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
## 2. Subset Endothelial Cells
## -----------------------------------------------------------

# Example:
# SeuratObj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

ec <- subset(SeuratObj, subset = majorclass %in% c("Endothelial", "Vascular_related"))

## -----------------------------------------------------------
## 3. Re-normalization and Variable Feature Selection
## -----------------------------------------------------------

ec <- NormalizeData(ec)
ec <- FindVariableFeatures(ec, selection.method = "vst", nfeatures = 3000)
ec <- ScaleData(ec, features = rownames(ec))

## -----------------------------------------------------------
## 4. PCA, Neighbors, Reclustering
## -----------------------------------------------------------

ec <- RunPCA(ec)
ec <- FindNeighbors(ec, dims = 1:20)
ec <- FindClusters(ec, resolution = 0.6)
ec <- RunUMAP(ec, dims = 1:20)

p_ec <- DimPlot(ec, reduction = "umap", label = TRUE) +
        ggtitle("Endothelial System Subclusters")
print(p_ec)

## -----------------------------------------------------------
## 5. Marker Genes for Endothelial Subclusters
## -----------------------------------------------------------

ec_markers <- FindAllMarkers(
  ec,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# write.csv(ec_markers, "results/endothelial/endothelial_markers.csv", row.names = FALSE)

## -----------------------------------------------------------
## 6. Endothelial Functional Subtype Annotation
## -----------------------------------------------------------

# Canonical endothelial lineage markers (mouse)
ec_reference <- list(
  Capillary_EC      = c("Kdr", "Esam", "Eng", "Icam2"),
  Arterial_EC       = c("Dll4", "Efnb2", "Gja4"),
  Venous_EC         = c("Nrp2", "Nr2f2", "Vwf"),
  Lymphatic_EC      = c("Pdpn", "Prox1", "Flt4"),
  Inflammatory_EC   = c("Sele", "Vcam1", "Icam1"),
  EndMT_like        = c("Klf4", "Col1a1", "Fn1")
)

ec$subcelltype <- "Unassigned"

for (type in names(ec_reference)) {
  genes <- ec_reference[[type]]
  idx <- WhichCells(ec, expression = rownames(ec) %in% genes)
  ec$subcelltype[idx] <- type
}

p_ec_type <- DimPlot(ec, reduction = "umap", group.by = "subcelltype", label = TRUE) +
             ggtitle("Endothelial Functional States")
print(p_ec_type)

## -----------------------------------------------------------
## 7. Optional Visualization
## -----------------------------------------------------------

# FeaturePlot(ec, features = c("Kdr", "Dll4", "Vwf", "Prox1"))
# VlnPlot(ec, features = c("Klf4", "Vcam1"), group.by = "subcelltype")

## -----------------------------------------------------------
## 8. Save Output
## -----------------------------------------------------------

# saveRDS(ec, file = "results/endothelial/SeuratObj_endothelial.rds")
