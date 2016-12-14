library(limma)
library(d3heatmap)

#' Plot an interactive (D3) Heatmap
#'   
#' Use this function to run limma differential expression and load an interactive D3 heatmap
#'
#' @param data a SCEset object
#' 
#' @return D3 heatmap will load
#' @export plotHeatmap
#'
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