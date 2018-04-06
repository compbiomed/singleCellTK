#' Get PCA variances
#'
#' @param x SCtkE object
#' @param ... other parameters
#'
#' @exportMethod pcaVariances
#' @examples
#' data("mouseBrainSubsetSCE")
#' pcaVariances(mouseBrainSubsetSCE)
setGeneric("pcaVariances", function(x, ...) standardGeneric("pcaVariances"))

#' Set PCA variances
#'
#' @param x SCtkE object
#' @param ... other parameters
#' @param value PCA variances DataFrame()
#'
#' @return A SCtkExperiment object with the pcaVariances slot set.
#'
#' @exportMethod pcaVariances<-
#' @examples
#' data("mouseBrainSubsetSCE")
#' pcaVariances(mouseBrainSubsetSCE)
#' #getPCA() sets the pcaVariances
#' newSCE <- getPCA(mouseBrainSubsetSCE, useAssay = "counts")
#'
#' #alternatively, set the pcaVariances directly
#' pca <- prcomp(assay(mouseBrainSubsetSCE, "logcounts"))
#' percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
#' pcaVariances(mouseBrainSubsetSCE) <- DataFrame(percentVar)
setGeneric("pcaVariances<-",
           function(x, ..., value) standardGeneric("pcaVariances<-"))
