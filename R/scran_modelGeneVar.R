#' Calculate Variable Genes with Scran modelGeneVar
#' 
#' @description Generates and stores variability data in the input 
#' \linkS4class{SingleCellExperiment} object, using 
#' \code{\link[scran]{modelGeneVar}} method. 
#' 
#' Also selects a specified number of top HVGs and store the logical selection 
#' in \code{rowData}. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object
#' @param useAssay A character string to specify an assay to compute variable 
#' features from. Default \code{"logcounts"}.
#' @return \code{inSCE} updated with variable feature metrics in \code{rowData}
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, "logcounts")
#' sce <- runModelGeneVar(sce)
#' hvf <- getTopHVG(sce, method = "modelGeneVar", hvgNumber = 10,
#'           useFeatureSubset = NULL)
#' @seealso \code{\link{runFeatureSelection}}, \code{\link{runSeuratFindHVG}},
#' \code{\link{getTopHVG}}, \code{\link{plotTopHVG}}
#' @importFrom SummarizedExperiment assay rowData rowData<-
#' @importFrom SingleCellExperiment rowSubset
#' @importFrom S4Vectors metadata<-
runModelGeneVar <- function(inSCE,
                            useAssay = "logcounts") {
    tempDataFrame <- data.frame(scran::modelGeneVar(assay(inSCE, useAssay)))
    rowData(inSCE)$scran_modelGeneVar_mean <- tempDataFrame$mean
    rowData(inSCE)$scran_modelGeneVar_totalVariance <- tempDataFrame$total
    rowData(inSCE)$scran_modelGeneVar_bio <- tempDataFrame$bio
    metadata(inSCE)$sctk$runFeatureSelection$modelGeneVar <- 
        list(useAssay = useAssay,
             rowData = c("scran_modelGeneVar_mean", 
                         "scran_modelGeneVar_totalVariance",
                         "scran_modelGeneVar_bio"))
    return(inSCE)
}
