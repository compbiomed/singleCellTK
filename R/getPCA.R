#' Get and plot PCA components for a SingleCellExperiment object
#'
#' Selects the 500 most variable genes in the SCE, performs
#' PCA based on them and stores the values in the reducedDims slot of
#' the SCE object.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the PCA data with this name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A \linkS4class{SingleCellExperiment} object with the specified
#' reducedDim
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #add a CPM assay
#' assay(mouseBrainSubsetSCE, "cpm") <- apply(assay(mouseBrainSubsetSCE,
#'                                                  "counts"),
#'                                            2, function(x) {
#'                                              x / (sum(x) / 1000000)
#'                                            })
#' mouseBrainSubsetSCE <- getPCA(mouseBrainSubsetSCE,
#'                               useAssay = "cpm",
#'                               reducedDimName = "PCA_cpm")
#' reducedDims(mouseBrainSubsetSCE)
#'
getPCA <- function(inSCE, useAssay="logcounts", reducedDimName="PCA"){
  if (nrow(inSCE) < 500){
    ntop <- nrow(inSCE)
  } else{
    ntop <- 500
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(useAssay, " not in the assay list")
  }
  exprsMat <- SummarizedExperiment::assay(inSCE, useAssay)
  if (!inherit(exprsMat, 'matrix')){
    exprsMat <- as.matrix(exprsMat)
  }
  rv <- matrixStats::rowVars(exprsMat)
  featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprsToPlot <- exprsMat[featureSet, , drop = FALSE]
  exprsToPlot <- scale(t(exprsToPlot))
  keepFeature <- (matrixStats::colVars(exprsToPlot) > 0.001)
  keepFeature[is.na(keepFeature)] <- FALSE
  exprsToPlot <- exprsToPlot[, keepFeature]
  pca <- stats::prcomp(exprsToPlot)
  #colnames(pc) <- paste("PC", seq_along(1:ncol(inSCE)), sep = "")
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  pca <- pca$x
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- pca
  return(inSCE)
}

