#' @describeIn getPCA Get t-SNE components for a SCtkE object
#'
#' @return getTSNE(): A SCtkE object with the specified reduecedDim and
#' pcaVariances updated
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #add a CPM assay
#' assay(mouseBrainSubsetSCE, "cpm") <- apply(
#'   assay(mouseBrainSubsetSCE, "counts"), 2, function(x) {
#'     x / (sum(x) / 1000000)
#'   })
#' mouseBrainSubsetSCE <- getTSNE(mouseBrainSubsetSCE, useAssay = "cpm",
#'                                reducedDimName = "TSNE_cpm")
#' reducedDims(mouseBrainSubsetSCE)
#'
getTSNE <- function(inSCE, useAssay="logcounts", reducedDimName="TSNE"){
  if (nrow(inSCE) < 500){
    ntop <- nrow(inSCE)
  } else{
    ntop <- 500
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(useAssay, " not in the assay list")
  }
  exprsMat <- SummarizedExperiment::assay(inSCE, useAssay)
  rv <- matrixStats::rowVars(exprsMat)
  featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprsToPlot <- exprsMat[featureSet, ]
  keepFeature <- (matrixStats::rowVars(exprsToPlot) > 0.001)
  keepFeature[is.na(keepFeature)] <- FALSE
  exprsToPlot <- exprsToPlot[keepFeature, ]
  exprsToPlot <- t(scale(t(exprsToPlot)))
  perplexity <- floor(ncol(inSCE) / 5)
  tsneOut <- Rtsne::Rtsne(t(exprsToPlot), perplexity = perplexity,
                           initial_dims = max(50, ncol(inSCE)))
  tsneOut <- tsneOut$Y[, c(1, 2)]
  rownames(tsneOut) <- colnames(inSCE)
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- tsneOut
  return(inSCE)
}
