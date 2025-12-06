# -----------------------------------------------------------
# Script: Sepsis_3.Analysis_miloR&OR.R
# Module: 03_crosscell_analysis
#
# Description:
#   Cross-condition differential abundance analysis using MiloR and
#   odds ratio (OR) statistics in sepsis-induced cardiomyopathy.
#   This script constructs a kNN graph on the low-dimensional space,
#   defines Milo neighborhoods, tests differential abundance between
#   conditions (e.g., Sham vs CLP-24h vs CLP-48h), and summarizes
#   results at the level of cell states / clusters.
#
# Inputs:
#   - Single-cell object with:
#       * Low-dimensional embeddings (e.g., PCA/UMAP)
#       * Cluster labels (e.g., seurat_clusters or subcelltype)
#       * Condition labels (e.g., Sham / CLP_24h / CLP_48h)
#
#   Typically derived from a Seurat object converted to a
#   SingleCellExperiment.
#
# Outputs:
#   - MiloR neighborhood differential abundance statistics
#   - Odds ratios and confidence intervals per state/cluster
#   - Plots summarizing DA neighborhoods and enriched states
#
# Dependencies:
#   Seurat (>=4.0), SingleCellExperiment, miloR, dplyr, ggplot2
# -----------------------------------------------------------


## -----------------------------------------------------------
## 1. Load Libraries
## -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(miloR)
  library(dplyr)
  library(ggplot2)
})


## -----------------------------------------------------------
## 2. Prepare Input Object
## -----------------------------------------------------------
# 假设已有一个 Seurat 对象，带有：
#   - UMAP/PCA
#   - meta.data 中的：
#       condition: Sham / CLP_24h / CLP_48h
#       cluster:   seurat_clusters 或细胞状态标签

# 示例：
# seurat_obj <- readRDS("results/celltype/SeuratObj_celltyped.rds")

# 转换为 SingleCellExperiment 以供 miloR 使用
sce <- as.SingleCellExperiment(seurat_obj)

# 将 UMAP 或 PCA 作为 miloR 的 reducedDim
# 这里以 PCA 为例，可以根据你实际分析替换为 "umap"
reducedDim(sce, "MILO") <- t(seurat_obj@reductions$pca@cell.embeddings)

# 确保 condition / cluster 信息在 colData 中
colData(sce)$condition <- seurat_obj$condition
colData(sce)$cluster   <- seurat_obj$subcelltype %||% seurat_obj$seurat_clusters


## -----------------------------------------------------------
## 3. Initialize Milo Object and Build Graph
## -----------------------------------------------------------

milo_obj <- Milo(sce)

# 构建 kNN 图（参数可根据你的数据规模与噪音水平调整）
milo_obj <- buildGraph(
  milo_obj,
  k = 30,       # number of nearest neighbors
  d = 30,       # dimensions in reduced space
  reduced.dim = "MILO"
)

# 定义 neighborhoods
milo_obj <- makeNhoods(
  milo_obj,
  prop = 0.1,    # proportion of cells used to define neighborhoods
  k = 30,
  d = 30,
  refined = TRUE
)

# 计算每个样本在各 neighborhood 中的细胞数
milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)))


## -----------------------------------------------------------
## 4. Construct Design Matrix and Run Milo DA Test
## -----------------------------------------------------------

# 设计矩阵：condition 作为自变量（可进一步加入 batch 等协变量）
design_df <- data.frame(
  sample_id = milo_obj$nhoodCounts$sample_id,
  condition = milo_obj$nhoodCounts$condition
)

design <- model.matrix(~ 0 + condition, data = design_df)
rownames(design) <- design_df$sample_id

# 运行 Milo 差异丰度检测
da_res <- testNhoods(
  milo_obj,
  design = design,
  design.df = design_df
)

# 添加 FDR 校正
da_res$FDR <- p.adjust(da_res$PValue, method = "BH")

# 保存结果
# write.csv(da_res, "results/milo/milo_neighborhood_DA.csv", row.names = FALSE)


## -----------------------------------------------------------
## 5. Map Neighborhood-Level DA Back to Cell States / Clusters
## -----------------------------------------------------------

# 计算每个 neighborhood 的主导 cluster/state
nhood_annot <- nhoodResults(milo_obj)
nhood_membership <- nhoodGraph(milo_obj)$Nhoods

# 简单示例：按细胞 cluster 的多数票 assigned 给 neighborhood
# 注意：实际可根据你原始脚本逻辑替换 / 加权
cluster_by_nhood <- apply(nhood_membership, 2, function(cells_idx) {
  cl <- colData(milo_obj)$cluster[cells_idx]
  names(sort(table(cl), decreasing = TRUE))[1]
})

da_res$cluster <- cluster_by_nhood[match(da_res$nhood, names(cluster_by_nhood))]

# 以 cluster 为单位计算 OR：示意写法（你可以用你原本更精细的统计）
or_by_cluster <- da_res %>%
  group_by(cluster) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    significant_nhoods = sum(FDR < 0.05)
  )

# write.csv(or_by_cluster, "results/milo/milo_OR_by_cluster.csv", row.names = FALSE)


## -----------------------------------------------------------
## 6. Example Visualization: DA on UMAP
## -----------------------------------------------------------

# 将 neighborhood DA 映射到 UMAP (使用 milo 内置函数或自定义)
# 这里仅示意：用 -log10(P) 作为 DA 强度
da_res$negLog10P <- -log10(da_res$PValue + 1e-300)

# 你可以使用 plotNhoodGraphDA(milo_obj, da_res, ...)，
# 此处给出自定义绘图示意
# plotNhoodGraphDA(milo_obj, da_res, layout = "UMAP")

## -----------------------------------------------------------
## 7. Odds Ratio Analysis (示意框架)
## -----------------------------------------------------------

# 典型 OR 分析：按 cluster × condition 构建 2×2 表
# 这里仅给出一个模板供你替换为你实际的分组方案与统计函数

compute_or_for_cluster <- function(meta_df, cluster_name, cond_case, cond_ctrl) {
  in_cluster  <- meta_df$cluster == cluster_name
  in_case     <- meta_df$condition == cond_case
  in_ctrl     <- meta_df$condition == cond_ctrl

  # 2×2 表
  a <- sum(in_cluster & in_case)
  b <- sum(!in_cluster & in_case)
  c <- sum(in_cluster & in_ctrl)
  d <- sum(!in_cluster & in_ctrl)

  tab <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(
                  cluster = c("in", "out"),
                  condition = c("case", "ctrl")
                ))

  # 使用 fisher.test 估计 OR
  ft <- fisher.test(tab)
  data.frame(
    cluster = cluster_name,
    case = cond_case,
    ctrl = cond_ctrl,
    OR = ft$estimate,
    lower = ft$conf.int[1],
    upper = ft$conf.int[2],
    p.value = ft$p.value
  )
}

meta_df <- data.frame(
  cluster   = colData(sce)$cluster,
  condition = colData(sce)$condition
)

clusters <- unique(meta_df$cluster)

or_results <- do.call(
  rbind,
  lapply(
    clusters,
    function(cl) compute_or_for_cluster(
      meta_df = meta_df,
      cluster_name = cl,
      cond_case = "CLP_48h",
      cond_ctrl = "Sham"
    )
  )
)

or_results$FDR <- p.adjust(or_results$p.value, method = "BH")

# write.csv(or_results, "results/milo/milo_cluster_OR.csv", row.names = FALSE)


## -----------------------------------------------------------
## 8. Save Key Outputs
## -----------------------------------------------------------

# saveRDS(milo_obj, file = "results/milo/milo_object.rds")
# saveRDS(da_res,    file = "results/milo/milo_DA_results.rds")
# saveRDS(or_results,file = "results/milo/milo_OR_results.rds")
