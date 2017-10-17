#' A lightweight S4 extension to the SingleCellExperiment class to store
#' additional information.
#'
#' @slot pca_variances The percent variation contained in each PCA dimension
#'
#' @exportClass SingleCelltkExperiment
#'
setClass("SingleCelltkExperiment",
         slots = c(pca_variances = "DataFrame"),
         contains = "SingleCellExperiment")

#' Create a SingleCelltkExperiment
#'
#' @param ... SingleCellExperiment and SummarizedExperiment components
#' @param pca_variances The percent variation contained in each PCA dimension
#'
#' @export SingleCelltkExperiment
#'
SingleCelltkExperiment <- function(..., pca_variances = DataFrame()) {
  sce <- SingleCellExperiment(...)
  out <- new("SingleCelltkExperiment", sce, pca_variances = DataFrame())
  return(out)
}
