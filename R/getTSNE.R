#' Get t-SNE components for a SCE object
#'
#' Selects the 500 most variable genes in the feature count, performs
#' t-SNE based on them and stores the t-SNE values in the reducedDims slot of the
#' SCE object.
#'
#' @param countData SCE object
#' @param useAssay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the t-SNE data with this name. The default is
#' TSNE. The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCE object with the specified reducedDim updated
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
getTSNE <- function(countData, useAssay="logcounts", reducedDimName="TSNE"){
  if (nrow(countData) < 500){
    ntop <- nrow(countData)
  } else{
    ntop <- 500
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(countData)))){
    stop(useAssay, " not in the assay list")
  }
  exprsMat <- log2(SummarizedExperiment::assay(countData, useAssay) + 1)
  rv <- matrixStats::rowVars(exprsMat)
  featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprsToPlot <- exprsMat[featureSet, ]
  keepFeature <- (matrixStats::rowVars(exprsToPlot) > 0.001)
  keepFeature[is.na(keepFeature)] <- FALSE
  exprsToPlot <- exprsToPlot[keepFeature, ]
  exprsToPlot <- t(scale(t(exprsToPlot)))
  perplexity <- floor(ncol(countData) / 5)
  tsneOut <- Rtsne::Rtsne(t(exprsToPlot), perplexity = perplexity,
                           initial_dims = max(50, ncol(countData)))
  tsneOut <- tsneOut$Y[, c(1, 2)]
  rownames(tsneOut) <- colnames(countData)
  SingleCellExperiment::reducedDim(countData, reducedDimName) <- tsneOut
  return(countData)
}
