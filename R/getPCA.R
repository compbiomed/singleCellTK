#' Get and plot PCA components for a SCtkE object
#'
#' Selects the 500 most variable genes in the SCE, performs
#' PCA based on them and stores the values in the reducedDims slot of
#' the SCE object.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the PCA data with this name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCtkE object with the specified reducedDim and
#' pcaVariances updated
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
getPCA <- function(inSCE, useAssay="logcounts", reducedDimName="PCA", ntop = 500){
  if (nrow(inSCE) < ntop){
    ntop <- nrow(inSCE)
  } else{
    ntop <- ntop
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(useAssay, " not in the assay list")
  }
   # exprsMat <- SummarizedExperiment::assay(inSCE, useAssay)
   # if (!is.matrix(exprsMat)){
   #   #stop("Input matrix ", useAssay, " is not a matrix")
   #   exprsMat <- as.matrix(exprsMat)
   # }
   # rv <- matrixStats::rowVars(exprsMat)
   # featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
   # exprsToPlot <- exprsMat[featureSet, , drop = FALSE]
   # exprsToPlot <- scale(t(exprsToPlot))
   # keepFeature <- (matrixStats::colVars(exprsToPlot) > 0.001)
   # keepFeature[is.na(keepFeature)] <- FALSE
   # exprsToPlot <- exprsToPlot[, keepFeature]
   # pca <- stats::prcomp(exprsToPlot)
   # #colnames(pc) <- paste("PC", seq_along(1:ncol(inSCE)), sep = "")
   # percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
   # pca <- pca$x
   # print("pca")
   # print(pca)
   # print(class(pca))
   # print(type(pca))
   # SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- pca
   # if (class(inSCE) == "SCtkExperiment"){
   #   print("pca vairnaces")
   #   pcaVariances(inSCE) <- S4Vectors::DataFrame(percentVar)
   #   print(pcaVariances(inSCE))
   #   rownames(pcaVariances(inSCE)) <- paste0(
   #     "PC", seq_len(nrow(pcaVariances(inSCE))))
   # }

 # SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- reducedDim(scater::runPCA(inSCE, name = reducedDimName, exprs_values = useAssay), reducedDimName)
  inSCE <- scater::runPCA(inSCE, name = reducedDimName, exprs_values = useAssay, ntop = ntop, scale = TRUE)

  #pcaVariances(inSCE) <- S4Vectors::DataFrame(attr(reducedDims(inSCE)[[reducedDimName]], "percentVar"))
  #rownames(pcaVariances(inSCE)) <- paste0("PC", seq_len(nrow(pcaVariances(inSCE))))
  return(inSCE)
}

