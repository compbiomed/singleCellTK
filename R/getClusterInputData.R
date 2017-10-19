#' Get data to use as input clustering algorithms
#'
#' @param count_data A SCE object
#' @param inputData A string ("Raw Data", "PCA Components", "tSNE Components")
#' @param use_assay Indicate which assay to use for PCA. Default is "logcounts"
#' @param reducedDimName If clustering on PCA or TSNE data, dimension name.
#' The toolkit will store data with the pattern <ASSSAY>_<ALGORITHM>.
#' 
#' @return Cluster input data
#' @export getClusterInputData
#'
getClusterInputData <- function(count_data, inputData, use_assay="logcounts",
                                reducedDimName=NULL){
  if (inputData == "Raw Data"){
    e <- assay(count_data, use_assay)
  } else if (inputData == "PCA Components") {
    if(is.null(reducedDimName)){
      stop("You must supply a PCA dim name")
    }
    if (is.null(reducedDim(count_data, reducedDimName))) {
      count_data <- getPCA(count_data, use_assay = use_assay,
                           reducedDimName = reducedDimName)
    }
    e <- reducedDim(count_data, reducedDimName)
  } else if (inputData == "tSNE Components") {
    if(is.null(reducedDimName)){
      stop("You must supply a tSNE dim name")
    }
    if (is.null(reducedDim(count_data, "TSNE"))) {
      count_data <- getTSNE(count_data, use_assay = use_assay,
                            reducedDimName = reducedDimName)
    }
    e <- reducedDim(count_data, reducedDimName)
  }
  return(e)
}
