#' Get data to use as input clustering algorithms
#'
#' @param count_data A SCE object
#' @param inputData A string ("Raw Data", "PCA Components", "tSNE Components")
#' @param use_assay Indicate which assay to use for PCA. Default is "logcounts"
#' @param reducedDimName If clustering on PCA or t-SNE data, dimension name.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return Cluster input data
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' getClusterInputData(GSE60361_subset_sce, "PCA Components",
#'                     use_assay = "logcounts", reducedDimName = "PCA_logcounts")
#'
getClusterInputData <- function(count_data, inputData, use_assay="logcounts",
                                reducedDimName=NULL){
  if (inputData == "Raw Data"){
    e <- SummarizedExperiment::assay(count_data, use_assay)
  } else if (inputData == "PCA Components") {
    if (is.null(reducedDimName)){
      stop("You must supply a PCA dim name")
    }
    if (is.null(SingleCellExperiment::reducedDim(count_data, reducedDimName))) {
      count_data <- getPCA(count_data, use_assay = use_assay,
                           reducedDimName = reducedDimName)
    }
    e <- SingleCellExperiment::reducedDim(count_data, reducedDimName)
  } else if (inputData == "tSNE Components") {
    if (is.null(reducedDimName)){
      stop("You must supply a tSNE dim name")
    }
    if (is.null(SingleCellExperiment::reducedDim(count_data, "TSNE"))) {
      count_data <- getTSNE(count_data, use_assay = use_assay,
                            reducedDimName = reducedDimName)
    }
    e <- SingleCellExperiment::reducedDim(count_data, reducedDimName)
  }
  return(e)
}
