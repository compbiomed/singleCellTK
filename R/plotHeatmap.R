library(shiny)
library(scater)
library(ComplexHeatmap)
library(biomaRt)
library(circlize)
library(limma)
library(d3heatmap)

plotHeatmap <- function(data) {
  counts(data) <- counts(data) + 1
  y <- log(counts(data))
  
  # Empirical Bayes Moderated T-Test
  design <- model.matrix(~factor(data$Tissue))
  fit <- lmFit(y, design)
  ebayes <- eBayes(fit)
  topGenes <- topTable(ebayes, coef=2, adjust="fdr", n=50)
  
  d3heatmap::d3heatmap(y[rownames(topGenes),])
}