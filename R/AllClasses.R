#' A lightweight S4 extension to the SingleCellExperiment class to store
#' additional information.
#'
#' @slot pca_variances The percent variation contained in each PCA dimension
#'
#' @param value The DataFrame of pca_variances
#'
#' @return A SingleCellExperiment like object with an addition pca_variances
#' slot.
#' @exportClass SCtkExperiment
#'
setClass("SCtkExperiment",
         slots = c(pca_variances = "DataFrame"),
         contains = "SingleCellExperiment")

#' Create a SCtkExperiment
#'
#' @param ... SingleCellExperiment and SummarizedExperiment components
#' @param pca_variances The percent variation contained in each PCA dimension
#'
#' @return A SingleCellExperiment like object with an addition pca_variances
#' slot.
#' @export
#'
SCtkExperiment <- function(..., pca_variances = DataFrame()) {
  sce <- SingleCellExperiment(...)
  out <- new("SCtkExperiment", sce, pca_variances = DataFrame())
  return(out)
}
