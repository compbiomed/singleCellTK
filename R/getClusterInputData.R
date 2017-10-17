#' Get data to use as input clustering algorithms
#'
#' @param count_data A SCE object
#' @param inputData A string ("Raw Data", "PCA Components", "tSNE Components")
#'
#' @return Cluster input data
#' @export getClusterInputData
#'
getClusterInputData <- function(count_data, inputData){
  if (inputData == "Raw Data"){
    e <- t(log2(assay(count_data, "counts") + 1))
  } else if (inputData == "PCA Components") {
    if (is.null(reducedDim(count_data, "PCA"))) {
      count_data <- getPCA(count_data)
    }
    e <- reducedDim(count_data, "PCA")
  } else if (inputData == "tSNE Components") {
    if (is.null(reducedDim(count_data, "TSNE"))) {
      count_data <- getTSNE(count_data)
    }
    e <- reducedDim(count_data, "TSNE")
  }
  return(e)
}
