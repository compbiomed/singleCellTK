#' @title Run Cluster Summary Metrics
#' @description Calculates the mean expression of percent of cells that express the
#' given genes for each cluster
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param gene A string or vector of strings with each gene to aggregate.
#' @param displayName A string that is the name of the column used for genes.
#' @param clusters The name of a colData entry that can be used as groups.
#' @return A dataframe with mean expression and percent of cells in cluster that 
#' express for each cluster.
#' @examples
#' runClusterSummaryMetrics(inSCE=sce, useAssay="counts", gene=c("B2M"), 
#' displayName="feature_name", clusters="cluster")
#' @export

runClusterSummaryMetrics <- function(inSCE, useAssay="logcounts", gene, displayName, clusters="cluster"){
  if (!clusters %in% names(SingleCellExperiment::colData(inSCE))) {
    stop("Specified variable '", clusters, "' not found in colData(inSCE)")
  }
  falseGenes <- setdiff(gene, rowData(inSCE)[[displayName]])
  gene <- intersect(gene, rowData(inSCE)[[displayName]])
  if (length(gene) == 0) {
    stop("All genes in '", toString(falseGenes), "' not found in rowData(inSCE)$", displayName)
  }
  if (length(falseGenes) > 0) {
    warning("Specified genes '", toString(falseGenes), "' not found in rowData(inSCE)$", displayName)
  }
  
  dimnames(inSCE)[[1]] <- rowData(inSCE)[[displayName]]
  dimnames(assay(inSCE, i=useAssay, withDimnames=FALSE))[[1]] <- rowData(inSCE)[[displayName]]
  
  avgExpr <- data.frame(assay(scuttle::aggregateAcrossCells(inSCE, ids=SingleCellExperiment::colData(inSCE)[,clusters], 
                                statistics="mean", use.assay.type=useAssay, 
                                subset.row=gene)), check.names=FALSE)
  avgExpr$Gene <- row.names(avgExpr)
  avgExpr <- tidyr::gather(avgExpr, key="cluster", value="clusterAveExpr", -Gene)
  
  
  percExpr <- data.frame(assay(scuttle::aggregateAcrossCells(inSCE, ids=SingleCellExperiment::colData(inSCE)[,clusters], 
                                statistics="prop.detected", use.assay.type=useAssay, 
                                subset.row=gene)), check.names=FALSE)
  percExpr$Gene <- row.names(percExpr)
  percExpr <- tidyr::gather(percExpr, key="cluster", value="clusterExprPerc", -Gene)
  
  summaryMetrics <- merge(percExpr, avgExpr)
  
  return(summaryMetrics)
}