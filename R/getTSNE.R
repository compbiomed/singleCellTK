#' Run t-SNE dimensionality reduction method on a SingleCellExperiment Object
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for tSNE computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"logcounts"}
#' @param useAltExp The subset to use for tSNE computation, usually for the
#' selected.variable features. Default \code{NULL}.
#' @param reducedDimName a name to store the results of the dimension
#' reductions. Default \code{"TSNE"}.
#' @param n_iterations maximum iterations. Default \code{1000}.
#' @param perplexity perplexity parameter. Default \code{NULL}.
#' @param run_pca run tSNE on PCA components? Default \code{TRUE}.
#' @param ntop Number of top features to use as a further variable feature
#' selection. Default \code{NULL}.
#' @return A \linkS4class{SingleCellExperiment} object with tSNE computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
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
getTSNE <- function(inSCE, useAssay = "logcounts", useAltExp = NULL,
                    reducedDimName = "TSNE", n_iterations = 1000,
                    perplexity = NULL, run_pca = TRUE, ntop = NULL){
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  if (!is.null(useAltExp)) {
    if (!(useAltExp %in% SingleCellExperiment::altExpNames(inSCE))) {
      stop("Specified altExp '", useAltExp, "' not found. ")
    }
    sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
    if (!(useAssay %in% SummarizedExperiment::assayNames(sce))) {
      stop("Specified assay '", useAssay, "' not found in the ",
           "specified altExp. ")
    }
  } else {
    if (!(useAssay %in% SummarizedExperiment::assayNames(inSCE))) {
      stop("Specified assay '", useAssay, "' not found. ")
    }
    sce <- inSCE
  }
  exprsMat <- as.matrix(SummarizedExperiment::assay(sce, useAssay))
  if (!is.null(ntop) && ntop < nrow(inSCE)) {
    rv <- matrixStats::rowVars(exprsMat)
    featureSet <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    exprsToPlot <- exprsMat[featureSet, ]
    keepFeature <- (matrixStats::rowVars(exprsToPlot) > 0.001)
    keepFeature[is.na(keepFeature)] <- FALSE
    exprsToPlot <- exprsToPlot[keepFeature, ]
    exprsToPlot <- scale(t(exprsToPlot))
  } else {
    exprsToPlot <- t(exprsMat)
  }

  if (is.null(perplexity)){
    perplexity <- floor(ncol(inSCE) / 5)
  }
  tsneOut <- Rtsne::Rtsne(exprsToPlot, perplexity = perplexity,
                          initial_dims = max(50, ncol(inSCE)),
                          max_iter = n_iterations, pca = run_pca)
  tsneOut <- tsneOut$Y[, c(1, 2)]
  rownames(tsneOut) <- colnames(inSCE)
  colnames(tsneOut) <- c("tSNE1", "tSNE2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- tsneOut
  return(inSCE)
}
