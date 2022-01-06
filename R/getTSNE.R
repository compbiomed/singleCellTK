#' Run t-SNE dimensionality reduction method on a SingleCellExperiment Object
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for tSNE computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"logcounts"}
#' @param useAltExp The subset to use for tSNE computation, usually for the
#' selected.variable features. Default \code{NULL}.
#' @param useReducedDim The low dimension representation to use for UMAP
#' computation. Default \code{NULL}.
#' @param reducedDimName a name to store the results of the dimension
#' reductions. Default \code{"TSNE"}.
#' @param nIterations maximum iterations. Default \code{1000}.
#' @param perplexity perplexity parameter. Default \code{30}.
#' @param run_pca run tSNE on PCA components? Default \code{TRUE}.
#' @param ntop Number of top features to use as a further variable feature
#' selection. Default \code{NULL}.
#' @param seed Random seed for reproducibility of tSNE results. 
#' Default \code{NULL} will use global seed in use by the R environment.
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
#'                                reducedDimName = "TSNE_cpm",
#'                                perplexity = NULL)
getTSNE <- function(inSCE, useAssay = "logcounts", useAltExp = NULL,
                    useReducedDim = NULL, reducedDimName = "TSNE",
                    nIterations = 1000, perplexity = 30, run_pca = TRUE,
                    ntop = NULL, seed = NULL){
  
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  
  if (is.null(useAssay) && is.null(useReducedDim)) {
    stop("`useAssay` and `useReducedDim` cannot be NULL at the same time.")
  } else if (!is.null(useAssay) && !is.null(useReducedDim)) {
    stop("`useAssay` and `useReducedDim` cannot be specified at the same time.")
  } else {
    if (!is.null(useReducedDim)) {
      if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
        stop("Specified `useReducedDim` not found.")
      }
      if (!is.null(useAltExp)) {
        warning("`useAltExp` will be ignored when using `useReducedDim`.")
      }
      if (isTRUE(run_pca)) {
        warning("using `useReducedDim`, `run_pca` and `ntop` forced to be FALSE/NULL")
        run_pca <- FALSE
        ntop <- NULL
      }
      sce <- inSCE
    } else {
      if (!is.null(useAltExp)) {
        if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
          stop("Specified `useAltExp` not found.")
        }
        sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
        if (!useAssay %in% SummarizedExperiment::assayNames(sce)) {
          stop("Specified `useAssay` not found in `useAltExp`.")
        }
      } else {
        if (!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
          stop("Specified `useAssay` not found.")
        }
        sce <- inSCE
      }
    }
  }
  
  if (!is.null(useAssay)) {
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
  } else {
    exprsToPlot <- SingleCellExperiment::reducedDim(sce, useReducedDim)
  }
  
  if (is.null(perplexity)){
    perplexity <- floor(ncol(inSCE) / 5)
  }
  
  if(!is.null(seed)){
    withr::with_seed(seed = seed,
                     code = tsneOut <- Rtsne::Rtsne(exprsToPlot, perplexity = perplexity,
                                                    initial_dims = max(50, ncol(inSCE)),
                                                    max_iter = nIterations, pca = run_pca)
                     )
  }
  else{
    tsneOut <- Rtsne::Rtsne(exprsToPlot, perplexity = perplexity,
                            initial_dims = max(50, ncol(inSCE)),
                            max_iter = nIterations, pca = run_pca) 
  }
  
  tsneOut <- tsneOut$Y[, c(1, 2)]
  rownames(tsneOut) <- colnames(inSCE)
  colnames(tsneOut) <- c("tSNE1", "tSNE2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- tsneOut
  return(inSCE)
}