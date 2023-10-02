#' @title Run Cluster Summary Metrics
#' @description Calculates the mean expression of percent of cells that express the
#' given genes for each cluster
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param featureNames A string or vector of strings with each gene to aggregate.
#' @param displayName A string that is the name of the column used for genes.
#' @param groupNames The name of a colData entry that can be used as groupNames.
#' @param scale Option to scale the data
#' @return A dataframe with mean expression and percent of cells in cluster that 
#' express for each cluster.
#' @examples
#' data("scExample")
#' runClusterSummaryMetrics(inSCE=sce, useAssay="counts", featureNames=c("B2M", "MALAT1"), 
#' displayName="feature_name", groupNames="type")
#' @export

runClusterSummaryMetrics <- function(inSCE, useAssay="logcounts", featureNames, displayName=NULL, groupNames="cluster", scale = FALSE){
  if(isTRUE(scale)){
    runNormalization(inSCE=inSCE, useAssay=useAssay, scale = TRUE, normalizationMethod = NULL, transformation = NULL, 
                     pseudocountsBeforeNorm = NULL, pseudocountsBeforeTransform = NULL)
  }
  if (!groupNames %in% names(SingleCellExperiment::colData(inSCE))) {
    stop("Specified variable '", groupNames, "' not found in colData(inSCE)")
  }
  if (!is.null(displayName)) {
    dimnames(inSCE)[[1]] <- rowData(inSCE)[[displayName]]
    dimnames(assay(inSCE, i=useAssay, withDimnames=FALSE))[[1]] <- rowData(inSCE)[[displayName]]
    falseGenes <- setdiff(featureNames, rowData(inSCE)[[displayName]])
    featureNames <- intersect(featureNames, rowData(inSCE)[[displayName]])
    warning <- paste0("rowData(inSCE)$", displayName)
  }
  else {
    falseGenes <- setdiff(featureNames, dimnames(inSCE)[[1]])
    featureNames <- intersect(featureNames, dimnames(inSCE)[[1]])
    warning <- "dimnames"
  }
  
  if (length(featureNames) == 0) {
    stop("All genes in '", toString(falseGenes), "' not found in ", warning)
  }
  if (length(falseGenes) > 0) {
    warning("Specified genes '", toString(falseGenes), "' not found in ", warning)
  }
  
  tempSCE <- inSCE[featureNames, ]

  
  if(isTRUE(scale)){
    runNormalization(inSCE=tempSCE, useAssay=useAssay,scale = TRUE, normalizationMethod = NULL, transformation = NULL,
                     pseudocountsBeforeNorm = NULL, pseudocountsBeforeTransform = NULL)
    useAssay <- "counts"
  }
  
  avgExpr <- assay(scuttle::aggregateAcrossCells(tempSCE, ids=SingleCellExperiment::colData(inSCE)[,groupNames], 
                                                            statistics="mean", use.assay.type=useAssay, 
                                                 subset.row=NULL))

  
  
  percExpr <- assay(scuttle::aggregateAcrossCells(tempSCE, ids=SingleCellExperiment::colData(inSCE)[,groupNames], 
                                                             statistics="prop.detected", use.assay.type=useAssay, 
                                                  subset.row=NULL))
  

  df <- data.frame(featureNames = featureNames)
  df$avgExpr <- avgExpr
  df$percExpr <- percExpr
  return(df)

}