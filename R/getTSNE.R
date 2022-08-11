#' Run t-SNE embedding with Rtsne method
#' @description T-Stochastic Neighbour Embedding (t-SNE) algorithm is commonly
#' for 2D visualization of single-cell data. This function wraps the
#' Rtsne \code{\link[Rtsne]{Rtsne}} function.
#'
#' With this funciton, users can create tSNE embedding directly from raw count
#' matrix, with necessary preprocessing including normalization, scaling,
#' dimension reduction all automated. Yet we still recommend having the PCA as
#' input, so that the result can match with the clustering based on the same
#' input PCA, and will be much faster.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for tSNE computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"logcounts"}.
#' @param useReducedDim The low dimension representation to use for UMAP
#' computation. Default \code{NULL}.
#' @param useAltExp The subset to use for tSNE computation, usually for the
#' selected.variable features. Default \code{NULL}.
#' @param reducedDimName a name to store the results of the dimension
#' reductions. Default \code{"TSNE"}.
#' @param logNorm Whether the counts will need to be log-normalized prior to
#' generating the tSNE via \code{\link{scaterlogNormCounts}}. Ignored when using
#' \code{useReducedDim}. Default \code{FALSE}.
#' @param useFeatureSubset Subset of feature to use for dimension reduction. A
#' character string indicating a \code{rowData} variable that stores the logical
#' vector of HVG selection, or a vector that can subset the rows of
#' \code{inSCE}. Default \code{NULL}.
#' @param nTop Automatically detect this number of variable features to use for
#' dimension reduction. Ignored when using \code{useReducedDim} or using
#' \code{useFeatureSubset}. Default \code{2000}.
#' @param center Whether data should be centered before PCA is applied. Ignored
#' when using \code{useReducedDim}. Default \code{TRUE}.
#' @param scale Whether data should be scaled before PCA is applied. Ignored
#' when using \code{useReducedDim}. Default \code{TRUE}.
#' @param pca Whether an initial PCA step should be performed. Ignored when
#' using \code{useReducedDim}. Default \code{TRUE}.
#' @param partialPCA Whether truncated PCA should be used to calculate principal
#' components (requires the irlba package). This is faster for large input
#' matrices. Ignored when using \code{useReducedDim}. Default \code{FALSE}.
#' @param initialDims Number of dimensions from PCA to use as input in tSNE.
#' Default \code{25}.
#' @param theta Numeric value for speed/accuracy trade-off (increase for less
#' accuracy), set to \code{0.0} for exact TSNE. Default \code{0.5}.
#' @param perplexity perplexity parameter. Should not be bigger than
#' \code{3 * perplexity < ncol(inSCE) - 1}. Default \code{30}. See
#' \code{\link[Rtsne]{Rtsne}} details for interpretation.
#' @param nIterations maximum iterations. Default \code{1000}.
#' @param numThreads Integer, number of threads to use using OpenMP, Default
#' \code{1}. \code{0} corresponds to using all available cores.
#' @param seed Random seed for reproducibility of tSNE results.
#' Default \code{NULL} will use global seed in use by the R environment.
#' @return A \linkS4class{SingleCellExperiment} object with tSNE computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' # Run from raw counts
#' sce <- getTSNE(inSCE = sce, useAssay = "counts", logNorm = TRUE, nTop = 2000,
#'                scale = TRUE, pca = TRUE)
#' \dontrun{
#' # Run from PCA
#' sce <- scaterlogNormCounts(sce, "logcounts")
#' sce <- runModelGeneVar(sce)
#' sce <- scaterPCA(sce, useAssay = "logcounts",
#'                  useFeatureSubset = "HVG_modelGeneVar2000", scale = TRUE)
#' sce <- getTSNE(sce, useReducedDim = "PCA")
#' }
#' @importFrom S4Vectors metadata<-
getTSNE <- function(inSCE, useAssay = "logcounts", useReducedDim = NULL,
                    useAltExp = NULL, reducedDimName = "TSNE", logNorm = FALSE,
                    useFeatureSubset = NULL, nTop = 2000, center = TRUE,
                    scale = TRUE, pca = TRUE, partialPCA = FALSE,
                    initialDims = 25, theta = 0.5, perplexity = 30,
                    nIterations = 1000, numThreads = 1, seed = NULL){
  params <- as.list(environment())
  params$inSCE <- NULL
  # Note: useMat = list(useAssay = useAssay, ...)
  # `.selectSCEMatrix()` does the check and corrects the conditions
  # When using useReducedDim, useAssay will be ignored
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                             useReducedDim = useReducedDim,
                             useAltExp = useAltExp,
                             returnMatrix = FALSE)$names
  params$useAssay <- useMat$useAssay
  useAssay <- useMat$useAssay

  if (!is.null(useAltExp)) {
    sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
  } else {
    sce <- inSCE
  }

  if (!is.null(useAssay)) {
    if (isTRUE(logNorm)) {
      sce <- scaterlogNormCounts(sce, assayName = "logcounts", useAssay = useAssay)
      useAssay <- "logcounts"
    }
    if (!is.null(useFeatureSubset)) {
      subset_row <- .parseUseFeatureSubset(inSCE, useFeatureSubset,
                                           altExpObj = sce,
                                           returnType = "logical")
      sce <- sce[subset_row,]
    } else {
      if (!is.null(nTop) && nTop < nrow(sce)) {
        suppressMessages({
          sce <- runFeatureSelection(sce, useAssay, method = "modelGeneVar")
          sce <- setTopHVG(sce, method = "modelGeneVar", hvgNumber = nTop,
                           featureSubsetName = "tsneHVG")
        })
        sce <- subsetSCERows(sce, rowData = "tsneHVG", returnAsAltExp = FALSE)
      }
    }
    mat <- t(SummarizedExperiment::assay(sce, useAssay))
    if ((!isTRUE(pca) & !isTRUE(partialPCA))) {
      initialDims <- NULL
    }
  } else {
    mat <- SingleCellExperiment::reducedDim(sce, useReducedDim)
    pca <- FALSE
    partialPCA <- FALSE
    scale <- FALSE
    center <- FALSE
    if (initialDims < ncol(mat)) {
      mat <- mat[,seq(initialDims)]
    }
  }

  if (is.null(perplexity)){
    perplexity <- floor(ncol(inSCE) / 5)
  }
  # Rtsne requires a matrix input
  mat <- as.matrix(mat)
  message(paste0(date(), " ... Computing Rtsne."))
  .withSeed(seed, {
    tsneOut <- Rtsne::Rtsne(mat, pca_scale = scale, pca_center = center,
                            pca = pca, partial_pca = partialPCA,
                            perplexity = perplexity,
                            initial_dims = initialDims,
                            max_iter = nIterations, theta = theta,
                            num_threads = numThreads)
  })
  tsneOut <- tsneOut$Y
  rownames(tsneOut) <- colnames(inSCE)
  colnames(tsneOut) <- c("tSNE1", "tSNE2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- tsneOut
  metadata(inSCE)$sctk$runDimReduce$embedding[[reducedDimName]] <- params
  return(inSCE)
}
