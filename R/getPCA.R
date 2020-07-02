#' getPCA
#'
#' A wrapper to the scater::runPCA function that inputs a sce object
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the PCA data with this name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCtkE object with the specified reducedDim and
#' pcaVariances updated
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #add a CPM assay
#' assay(mouseBrainSubsetSCE, "cpm") <- apply(assay(mouseBrainSubsetSCE,
#'                                                  "counts"),
#'                                            2, function(x) {
#'                                              x / (sum(x) / 1000000)
#'                                            })
#' mouseBrainSubsetSCE <- getPCA(mouseBrainSubsetSCE,
#'                               useAssay = "cpm",
#'                               reducedDimName = "PCA_cpm")
#' reducedDims(mouseBrainSubsetSCE)
#'
getPCA <- function(inSCE, useAssay="logcounts", reducedDimName="PCA", ntop = 500){
  if (nrow(inSCE) < ntop){
    ntop <- nrow(inSCE)
  } else{
    ntop <- ntop
  }
  if (!(useAssay %in% names(SummarizedExperiment::assays(inSCE)))){
    stop(useAssay, " not in the assay list")
  }
   
  inSCE <- scater::runPCA(inSCE, name = reducedDimName, exprs_values = useAssay, ntop = ntop, scale = TRUE)

  return(inSCE)
}

