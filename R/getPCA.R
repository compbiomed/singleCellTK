#' getPCA
#' A wrapper to \link[scater]{runPCA} function to compute principal component analysis (PCA) from a given \code{SingleCellExperiment} object.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Assay to use for PCA computation.
#' @param reducedDimName Name to use for the reduced output assay.
#' @param ndim Number of principal components to obtain from the PCA computation. Default \code{50}.
#' @param ntop Number of top features to use with PCA. Default \code{500}.
#'
#' @return A \code{SingleCellExperiment} object with stored PCA computation
#' @export
#'
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- getPCA(sce_chcl, useAssay = "logcounts")
getPCA <- function(inSCE, useAssay="logcounts", reducedDimName="PCA", ndim = 50, ntop = 500){
  if (nrow(inSCE) < ntop){
    ntop <- nrow(inSCE)
  } else{
    ntop <- ntop
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(paste(useAssay, " not in the assay list"))
  }
  
  inSCE <- scater::runPCA(inSCE, name = reducedDimName, exprs_values = useAssay, ncomponents = ndim,  ntop = ntop, scale = TRUE)

  return(inSCE)
}