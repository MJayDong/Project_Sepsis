# -----------------------------------------------------------
# Script: Sepsis_2.SCP-2.5.VascularNiche.R
# Module: 02_celltype_analysis
#
# Description:
#   Vascular niche–focused analysis in sepsis-induced cardiac tissue.
#   This script refines mesenchymal / perivascular / endothelial-related
#   populations to characterize vascular niche states, including
#   pericytes, contractile fibroblasts, progenitor-like fibroblasts,
#   and niche-associated endothelial cells.
#
# Inputs:
#   - Global Seurat object or mesenchymal/vascular subset with $majorclass
#     or $subcelltype annotations
#
# Outputs:
#   - Vascular niche–focused Seurat object
#   - Niche subcluster markers
#   - Functional subtype labels (e.g., contractile fibroblast,
#     progenitor-like fibroblast, pericytes, vascular niche EC)
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
## 2. Define Input Object and Subset Vascular Niche
## -----------------------------------------------------------
# 通常你会从 mesenchymal_vascular 或 global object 里抽取：
# mv <- readRDS("results/mesenchymal_vascular/SeuratObj_mesenchymal_vascular.rds")
# 这里假设对象名为 mv

# 示例：根据前一步脚本中的 subcelltype/majorclass 筛 vascular niche 相关群
# 可根据你实际标签调整
vn <- subset(
  mv,
  subset = subcelltype %in% c(
    "Fibroblast_ECM",
    "Fibroblast_progenitor",
    "Pericytes",
    "Smooth_muscle_cells",
    "EndMT_like"
  )
)

## -----------------------------------------------------------
## 3. Re-normalization and Feature Selection (Niche-focused)
## -----------------------------------------------------------

vn <- NormalizeData(vn)
vn <- FindVariableFeatures(vn, selection.method = "vst", nfeatures = 3000)
vn <- ScaleData(vn, features = rownames(vn))

## -----------------------------------------------------------
## 4. PCA, Neighbors, Reclustering
## -----------------------------------------------------------

vn <- RunPCA(vn)
vn <- FindNeighbors(vn, dims = 1:20)
vn <- FindClusters(vn, resolution = 0.5)
vn <- RunUMAP(vn, dims = 1:20)

p_vn <- DimPlot(vn, reduction = "umap", label = TRUE) +
        ggtitle("Vascular Niche Subclusters")
print(p_vn)

## -----------------------------------------------------------
## 5. Identify Vascular Niche Subcluster Markers
## -----------------------------------------------------------

vn_markers <- FindAllMarkers(
  vn,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# write.csv(vn_markers,
#           "results/vascular_niche/vascular_niche_subcluster_markers.csv",
#           row.names = FALSE)

## -----------------------------------------------------------
## 6. Functional Annotation of Niche States
## -----------------------------------------------------------

# 结合你原脚本里的标记（你之前有 Myh11 / Pdgfra / Col1a1 等）
vn_reference <- list(
  Contractile_fibroblast   = c("Myh11", "Acta2", "Tagln"),
  Progenitor_like_fibro    = c("Pdgfra", "Ly6a"),
  ECM_fibroblast           = c("Col1a1", "Col3a1", "Dcn"),
  Pericytes                = c("Rgs5", "Pdgfrb"),
  Vascular_niche_EC_like   = c("Kdr", "Pecam1")
)

vn$subcelltype_niche <- "Unassigned"

for (type in names(vn_reference)) {
  genes <- vn_reference[[type]]
  idx <- WhichCells(vn, expression = rownames(vn) %in% genes)
  vn$subcelltype_niche[idx] <- type
}

p_vn_type <- DimPlot(vn, reduction = "umap",
                     group.by = "subcelltype_niche", label = TRUE) +
             ggtitle("Vascular Niche Functional States")
print(p_vn_type)

## -----------------------------------------------------------
## 7. Optional Feature and Violin Plots
## -----------------------------------------------------------

# FeaturePlot(vn, features = c("Myh11", "Pdgfra", "Col1a1", "Rgs5"))
# VlnPlot(vn, features = c("Myh11", "Pdgfra"), group.by = "subcelltype_niche")

## -----------------------------------------------------------
## 8. Save Output
## -----------------------------------------------------------

# saveRDS(vn,
#         file = "results/vascular_niche/SeuratObj_vascular_niche.rds")
