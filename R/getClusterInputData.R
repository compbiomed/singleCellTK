#' Get data to use as input clustering algorithms
#'
#' @param countData A SCE object
#' @param inputData A string ("Raw Data", "PCA Components", "tSNE Components")
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"
#' @param reducedDimName If clustering on PCA or t-SNE data, dimension name.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return Cluster input data
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' getClusterInputData(mouseBrainSubsetSCE, "PCA Components",
#'                     useAssay = "logcounts", reducedDimName = "PCA_logcounts")
#'
getClusterInputData <- function(countData, inputData, useAssay="logcounts",
                                reducedDimName=NULL){
  if (inputData == "Raw Data"){
    e <- SummarizedExperiment::assay(countData, useAssay)
  } else if (inputData == "PCA Components") {
    if (is.null(reducedDimName)){
      stop("You must supply a PCA dim name")
    }
    if (is.null(SingleCellExperiment::reducedDim(countData, reducedDimName))) {
      countData <- getPCA(countData, useAssay = useAssay,
                          reducedDimName = reducedDimName)
    }
    e <- SingleCellExperiment::reducedDim(countData, reducedDimName)
  } else if (inputData == "tSNE Components") {
    if (is.null(reducedDimName)){
      stop("You must supply a tSNE dim name")
    }
    if (is.null(SingleCellExperiment::reducedDim(countData, "TSNE"))) {
      countData <- getTSNE(countData, useAssay = useAssay,
                           reducedDimName = reducedDimName)
    }
    e <- SingleCellExperiment::reducedDim(countData, reducedDimName)
  }
  return(e)
}
