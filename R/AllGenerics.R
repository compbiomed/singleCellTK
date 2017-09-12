#' Get PCA variances
#' 
#' @param x SCtkE object
#' @param ... other parameters
#'
#' @exportMethod pca_variances
setGeneric("pca_variances", function(x, ...) standardGeneric("pca_variances"))

#' Set PCA variances
#' 
#' @param x SCtkE object
#' @param ... other parameters
#' @param value PCA variances DataFrame()
#'
#' @exportMethod pca_variances<-
setGeneric("pca_variances<-", function(x, ..., value) standardGeneric("pca_variances<-"))