#' Get PCA variances
#' 
#' @rdname pca_variances
setMethod("pca_variances", "SingleCelltkExperiment", function(x) x@pca_variances)

#' Set PCA variances
#' 
#' @param value The DataFrame of pca_variances
#' 
#' @rdname pca_variances
setReplaceMethod("pca_variances", "SingleCelltkExperiment", function(x, value) {
  x@pca_variances <- value
  return(x)
})