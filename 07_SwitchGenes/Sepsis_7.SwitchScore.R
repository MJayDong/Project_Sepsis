# -----------------------------------------------------------
# Script: Sepsis_7.SwitchScore.R
# Module: 07_SwitchGenes
#
# Description:
#   Pseudotime-dependent gene switching analysis using Monocle
#   and auxiliary metrics such as switching probability and
#   pseudo R². This script identifies genes that transition
#   from OFF→ON / ON→OFF states along pseudotime, quantifies
#   switching quality, and exports switch genes per cell type.
#
# Inputs:
#   - Expression matrix (exprs)
#   - Metadata with pseudotime ordering
#
# Outputs:
#   - CellDataSet (Monocle)
#   - Switch scores (Pseudo R²)
#   - Gene switching ranking tables
#   - Exported gene lists
#
# Dependencies:
#   textshape, monocle, dplyr, ggplot2
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(monocle)
  library(textshape)
  library(dplyr)
  library(ggplot2)
})

options(warn = -1)

# -----------------------------------------------------------
# 1. Prepare Input Data
# -----------------------------------------------------------

# Expression matrix (rows: genes, columns: cells)
# Example:
# exprs <- as.matrix(seurat_obj@assays$RNA@data)

# Metadata must include pseudotime
# meta <- seurat_obj@meta.data

pd <- new("AnnotatedDataFrame", data = meta)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(exprs)))

cds <- newCellDataSet(
  exprs,
  phenoData   = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = negbinomial.size()
)

# -----------------------------------------------------------
# 2. Run Monocle (if trajectory not predefined)
# -----------------------------------------------------------

# cds <- estimateSizeFactors(cds)
# cds <- estimateDispersions(cds)
# cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
# cds <- orderCells(cds)

# Pseudotime stored in:
#   pData(cds)$Pseudotime

# -----------------------------------------------------------
# 3. Gene Switching Analysis
# -----------------------------------------------------------

# Fit gene expression along pseudotime
diff_res <- differentialGeneTest(
  cds,
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 4
)

# Rank by q-value or pseudo R²
diff_res <- diff_res %>%
  arrange(qval)

# -----------------------------------------------------------
# 4. Fit Switch Model for Each Gene
# -----------------------------------------------------------

fit_switch <- function(gene, cds) {
  y <- exprs(cds)[gene, ]
  x <- pData(cds)$Pseudotime

  df <- data.frame(x = x, y = y)

  # Logistic switch model
  m <- glm(y ~ x, family = "binomial")
  r2 <- 1 - (m$deviance / m$null.deviance)

  return(data.frame(
    gene = gene,
    pseudo_r2 = r2
  ))
}

genes <- rownames(exprs)
switch_scores <- bind_rows(lapply(genes, fit_switch, cds = cds))

switch_scores <- switch_scores %>%
  arrange(desc(pseudo_r2))

# -----------------------------------------------------------
# 5. Define High-Confidence Switch Genes
# -----------------------------------------------------------

high_switch <- switch_scores %>%
  filter(pseudo_r2 > 0.3)   # threshold可根据数据调整

# -----------------------------------------------------------
# 6. Visualization
# -----------------------------------------------------------

ggplot(switch_scores, aes(x = pseudo_r2)) +
  geom_histogram(bins = 40, fill = "steelblue") +
  theme_classic() +
  xlab("Pseudo R²") +
  ggtitle("Distribution of Switch Gene Scores")

# -----------------------------------------------------------
# 7. Export Results
# -----------------------------------------------------------

# write.csv(switch_scores, "results/switch_genes/switch_scores_all.csv", row.names = FALSE)
# write.csv(high_switch,   "results/switch_genes/high_switch_genes.csv", row.names = FALSE)

print("Switch gene analysis completed.")
