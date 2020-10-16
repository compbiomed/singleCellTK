#' scran_modelGeneVar
#' Generates and stores variability data from scran::modelGeneVar in the input singleCellExperiment object
#' @param inSCE a singleCellExperiment object
#' @param assayName selected assay to compute variable features from
#' @return inSCE updated singleCellExperiment object that contains variable feature metrics in rowData
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scran_modelGeneVar(sce_chcl, "counts")
#' @importFrom SummarizedExperiment assay rowData rowData<-
scran_modelGeneVar <- function(inSCE, assayName) {
    tempDataFrame <- data.frame(scran::modelGeneVar(assay(inSCE, assayName)))
    rowData(inSCE)$scran_modelGeneVar_mean <- tempDataFrame$mean
    rowData(inSCE)$scran_modelGeneVar_totalVariance <- tempDataFrame$total
    rowData(inSCE)$scran_modelGeneVar_bio <- tempDataFrame$bio
    return(inSCE)
}
