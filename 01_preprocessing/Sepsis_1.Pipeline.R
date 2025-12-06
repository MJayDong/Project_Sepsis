# -----------------------------------------------------------
# Script: Sepsis_1.Pipeline.R
# Module: 01_preprocessing
#
# Description:
#   Global preprocessing pipeline for sepsis single-cell RNA-seq data.
#   This script performs QC, filtering, normalization, variable gene
#   selection, PCA, UMAP embedding, and initial clustering.
#
# Inputs:
#   - Raw gene expression matrices or pre-constructed Seurat objects.
#
# Outputs:
#   - Preprocessed Seurat object
#   - QC metrics (nFeature_RNA, nCount_RNA, percent.mt/ribo)
#   - PCA and UMAP embeddings
#   - Initial clustering results
#
# Dependencies:
#   Seurat (>=4.0), dplyr, ggplot2, patchwork
#
# Notes:
#   - All file paths should be relative for GitHub reproducibility.
#   - Raw data loading should be adapted by the user.
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
## 2. Load Raw Data
## -----------------------------------------------------------
# Replace this section with your actual data loading step.
# For GitHub, avoid hard-coded HPC or absolute paths.
#
# Example:
# raw_data <- Read10X(data.dir = "data/raw/Sample1/")
# SeuratObj <- CreateSeuratObject(raw_data)

# Placeholder object name:
# SeuratObj <- <your_data_here>

## -----------------------------------------------------------
## 3. Calculate QC Metrics
## -----------------------------------------------------------

# Mitochondrial genes (update pattern if mouse/human differs)
mt_genes <- grep("^mt-", rownames(SeuratObj), value = TRUE)
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, features = mt_genes)

# Ribosomal genes
ribo_genes <- grep("^Rps|^Rpl", rownames(SeuratObj), value = TRUE)
SeuratObj[["percent.ribo"]] <- PercentageFeatureSet(SeuratObj, features = ribo_genes)

# Basic QC filtering (adjust thresholds as needed)
SeuratObj <- subset(
  SeuratObj,
  subset = nFeature_RNA > 500 &
           nCount_RNA > 1000 &
           percent.mt < 20
)

## -----------------------------------------------------------
## 4. Normalization & Variable Feature Selection
## -----------------------------------------------------------

SeuratObj <- NormalizeData(SeuratObj)
SeuratObj <- FindVariableFeatures(
  SeuratObj,
  selection.method = "vst",
  nfeatures = 3000
)

## -----------------------------------------------------------
## 5. Scaling
## -----------------------------------------------------------
SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj))

## -----------------------------------------------------------
## 6. Dimensionality Reduction (PCA)
## -----------------------------------------------------------

SeuratObj <- RunPCA(SeuratObj)

# (Optional) Visualize variance explained
# ElbowPlot(SeuratObj, ndims = 50)

## -----------------------------------------------------------
## 7. Neighbors, Clustering, and UMAP
## -----------------------------------------------------------

SeuratObj <- FindNeighbors(SeuratObj, dims = 1:30)
SeuratObj <- FindClusters(SeuratObj, resolution = 0.5)
SeuratObj <- RunUMAP(SeuratObj, dims = 1:30)

# Example visualization
# p1 <- DimPlot(SeuratObj, reduction = "umap", label = TRUE)
# print(p1)

## -----------------------------------------------------------
## 8. Save Output
## -----------------------------------------------------------
# saveRDS(SeuratObj, file = "results/preprocessing/SeuratObj_preprocessed.rds")
