#' Get PCA components for a SCtkE object
#'
#' Selects the 500 most variable genes in the SCE, performs
#' PCA based on them and stores the PCA values in the reducedDims slot of the
#' SCE object.
#'
#' @param countData SCtkE object
#' @param useAssay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the PCA data with this name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCtkE object with the specified reduecedDim and pcaVariances
#' updated
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
getPCA <- function(countData, useAssay="logcounts", reducedDimName="PCA"){
  if (nrow(countData) < 500){
    ntop <- nrow(countData)
  }
  else{
    ntop <- 500
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(countData)))){
    stop(useAssay, " not in the assay list")
  }
  exprsMat <- log2(SummarizedExperiment::assay(countData, useAssay) + 1)
  rv <- matrixStats::rowVars(exprsMat)
  featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprsToPlot <- exprsMat[featureSet, , drop = FALSE]
  exprsToPlot <- scale(t(exprsToPlot))
  keepFeature <- (matrixStats::colVars(exprsToPlot) > 0.001)
  keepFeature[is.na(keepFeature)] <- FALSE
  exprsToPlot <- exprsToPlot[, keepFeature]
  pca <- stats::prcomp(exprsToPlot)
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  pca <- pca$x
  SingleCellExperiment::reducedDim(countData, reducedDimName) <- pca
  if (class(countData) == "SCtkExperiment"){
    pcaVariances(countData) <- S4Vectors::DataFrame(percentVar)
    rownames(pcaVariances(countData)) <- paste0(
      "PC", seq_len(nrow(pcaVariances(countData))))
  }
  return(countData)
}
