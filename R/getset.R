#' Get PCA variances
#'
#' @return A data frame of percent variation explained by each PC.
#'
#' @rdname pcaVariances
setMethod("pcaVariances", "SCtkExperiment", function(x) x@pcaVariances)

#' Set PCA variances
#'
#' @param value The DataFrame of pcaVariances
#'
#' @return A SCtkExperiment object with the pcaVariances object set.
#'
#' @rdname pcaVariances
setReplaceMethod("pcaVariances", "SCtkExperiment", function(x, value) {
  x@pcaVariances <- value
  return(x)
})
