#' Run t-SNE dimensionality reduction method on the assay data.
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. The default is "logcounts".
#' @param reducedDimName a name to store the results of the dimension reductions
#' @param n_iterations maximum iterations. Default is 1000
#' @param perplexity perplexity parameter. Default is 5
#'
#' @return A SCtkE object with the specified reducedDim and
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
getTSNE <- function(inSCE, useAssay = "logcounts", reducedDimName = "TSNE",
                    n_iterations = 1000, perplexity = NULL){
  if (nrow(inSCE) < 500){
    ntop <- nrow(inSCE)
  } else{
    ntop <- 500
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(useAssay, " not in the assay list")
  }
  exprsMat <- SummarizedExperiment::assay(inSCE, useAssay)
  if (!is.matrix(exprsMat)){
    stop("Input matrix ", useAssay, " is not a matrix")
  }
  rv <- matrixStats::rowVars(exprsMat)
  featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprsToPlot <- exprsMat[featureSet, ]
  keepFeature <- (matrixStats::rowVars(exprsToPlot) > 0.001)
  keepFeature[is.na(keepFeature)] <- FALSE
  exprsToPlot <- exprsToPlot[keepFeature, ]
  exprsToPlot <- t(scale(t(exprsToPlot)))
  if (is.null(perplexity)){
    perplexity <- floor(ncol(inSCE) / 5)
  }
  tsneOut <- Rtsne::Rtsne(t(exprsToPlot), perplexity = perplexity,
                           initial_dims = max(50, ncol(inSCE)), max_iter = n_iterations)
  tsneOut <- tsneOut$Y[, c(1, 2)]
  rownames(tsneOut) <- colnames(inSCE)
  colnames(tsneOut) <- c("tSNE1", "tSNE2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- tsneOut
  return(inSCE)
}
