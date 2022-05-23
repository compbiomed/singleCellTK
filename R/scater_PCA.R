#' Perform scater PCA on a SingleCellExperiment Object
#' @description A wrapper to \link[scater]{runPCA} function to compute principal
#' component analysis (PCA) from a given \linkS4class{SingleCellExperiment} 
#' object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for PCA computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"logcounts"}
#' @param useHVGList A character string indicating a \code{rowData} variable 
#' that stores the logical vector of HVG selection. Default \code{NULL}.
#' @param scale Logical scalar, whether to standardize the expression values.
#' Default \code{TRUE}.
#' @param reducedDimName Name to use for the reduced output assay. Default
#' \code{"PCA"}.
#' @param nComponents Number of principal components to obtain from the PCA
#' computation. Default \code{50}.
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
#' @importFrom S4vectors metadata<-
scaterPCA <- function(inSCE, useAssay = "logcounts", useHVGList = NULL, 
                      scale = TRUE, reducedDimName = "PCA", nComponents = 50, 
                      useAltExp = NULL, seed = NULL, 
                      BPPARAM = BiocParallel::SerialParam()) {
  params <- as.list(environment())
  params$inSCE <- NULL
  params$BPPARAM <- NULL
  if (!is.null(useAltExp)) {
    if (!(useAltExp %in% SingleCellExperiment::altExpNames(inSCE))) {
      stop("Specified altExp '", useAltExp, "' not found. ")
    }
    sce <- altExp(inSCE, useAltExp)
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
  subset_row <- NULL
  if (!is.null(useHVGList)) {
    if (!useHVGList %in% colnames(rowData(inSCE))) {
      stop("Specified HVG list not found")
    }
    hvgs <- rownames(inSCE)[rowSubset(inSCE, useHVGList)]
    subset_row <- rownames(sce) %in% hvgs
  }
  message(paste0(date(), " ... Computing Scater PCA."))
  if (!is.null(seed)) {
    withr::with_seed(seed = seed,
                     code = sce <- scater::runPCA(sce, 
                                                  name = reducedDimName, 
                                                  exprs_values = useAssay,
                                                  ncomponents = nComponents, 
                                                  scale = scale, 
                                                  subset_row = subset_row,
                                                  BPPARAM = BPPARAM))
  } else {
    sce <- scater::runPCA(sce, name = reducedDimName, exprs_values = useAssay,
                          ncomponents = nComponents, scale = scale, 
                          subset_row = subset_row, BPPARAM = BPPARAM)
  }
  reducedDim(inSCE, reducedDimName) <- reducedDim(sce, reducedDimName)
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  return(inSCE)
}
