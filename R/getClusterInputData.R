#' Get data to use as input clustering algorithms
#'
#' @param count_data A SCE object
#' @param inputData A string ("Raw Data", "PCA Components", "tSNE Components")
#' @param vals Reactive Dataframe
#' 
#' @return Cluster input data
#' @export getClusterInputData
#'

getClusterInputData <- function(count_data, inputData,vals){
  if (inputData == "Raw Data"){
    e <- t(exprs(count_data))
  } else if (inputData == "PCA Components") {
    if (is.null(vals$PCA)) {
      vals$PCA <- getPCA(vals$counts)
    }
    e <- vals$PCA
  } else if (inputData == "tSNE Components") {
    if (is.null(vals$TSNE)) {
      vals$TSNE <- getTSNE(vals$counts)
    }
    e <- vals$TSNE
  }
  return(e)
}
