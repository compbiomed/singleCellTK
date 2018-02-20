#' Get PCA variances
#'
#' @param x SCtkE object
#' @param ... other parameters
#'
#' @exportMethod pca_variances
#' @examples
#' data("GSE60361_subset_sce")
#' pca_variances(GSE60361_subset_sce)
#' #getPCA() sets the pca_variances
#' newSCE <- getPCA(GSE60361_subset_sce, use_assay = "counts")
#'
#' #alternatively, set the pca_variances directly
#' pca <- prcomp(assay(GSE60361_subset_sce, "logcounts"))
#' percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
#' pca_variances(GSE60361_subset_sce) <- DataFrame(percentVar)
setGeneric("pca_variances", function(x, ...) standardGeneric("pca_variances"))

#' Set PCA variances
#'
#' @param x SCtkE object
#' @param ... other parameters
#' @param value PCA variances DataFrame()
#'
#' @exportMethod pca_variances<-
setGeneric("pca_variances<-", function(x, ..., value) standardGeneric("pca_variances<-"))
