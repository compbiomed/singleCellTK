#' Perform scater PCA on a SingleCellExperiment Object
#' @description A wrapper to \link[scater]{runPCA} function to compute principal
#' component analysis (PCA) from a given \linkS4class{SingleCellExperiment}
#' object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for PCA computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"logcounts"}
#' @param useFeatureSubset Subset of feature to use for dimension reduction. A
#' character string indicating a \code{rowData} variable that stores the logical
#' vector of HVG selection, or a vector that can subset the rows of
#' \code{inSCE}. Default \code{NULL}.
#' @param scale Logical scalar, whether to standardize the expression values.
#' Default \code{TRUE}.
#' @param reducedDimName Name to use for the reduced output assay. Default
#' \code{"PCA"}.
#' @param nComponents Number of principal components to obtain from the PCA
#' computation. Default \code{50}.
#' @param ntop Automatically detect this number of variable features to use for
#' dimension reduction. Ignored when using \code{useReducedDim} or using
#' \code{useFeatureSubset}. Default \code{2000}.
#' @param useAltExp The subset to use for PCA computation, usually for the
#' selected.variable features. Default \code{NULL}.
#' @param seed Integer, random seed for reproducibility of PCA results.
#' Default \code{NULL}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' the PCA should be parallelized.
#' @return A \linkS4class{SingleCellExperiment} object with PCA computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, "logcounts")
#' sce <- scaterPCA(sce, "logcounts", scale = TRUE)
#' @importFrom SingleCellExperiment reducedDim altExp rowSubset
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata<-
scaterPCA <- function(inSCE, useAssay = "logcounts", useFeatureSubset = NULL,
                      scale = TRUE, reducedDimName = "PCA", nComponents = 50,
                      ntop = 2000, useAltExp = NULL, seed = NULL,
                      BPPARAM = BiocParallel::SerialParam()) {
  params <- as.list(environment())
  params$inSCE <- NULL
  params$BPPARAM <- NULL
  # Note: useMat = list(names = list(useAssay = useAssay, ...),
  #                     mat = <matrix>)
  # `.selectSCEMatrix()` does the check
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                             useReducedDim = NULL, useAltExp = useAltExp,
                             returnMatrix = TRUE)
  if (!is.null(useAltExp)) {
    sce <- altExp(inSCE, useAltExp)
  } else {
    sce <- inSCE
  }
  subset_row <- .parseUseFeatureSubset(inSCE, useFeatureSubset,
                                       altExpObj = sce, returnType = "logical")
  message(paste0(date(), " ... Computing Scater PCA."))
  withr::with_seed(seed, {
    pca <- scater::calculatePCA(useMat$mat, ncomponents = nComponents,
                                scale = scale, ntop = ntop,
                                subset_row = subset_row, BPPARAM = BPPARAM)
  })
  reducedDim(inSCE, reducedDimName) <- pca
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  return(inSCE)
}



