#' Calculate Variable Genes with modelGeneVar
#' 
#' @description Generates and stores variability data in the input 
#' \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment} object, using 
#' \code{\link[scrapper]{modelGeneVariances}}, which replaces the deprecated
#' \code{scran::modelGeneVar}. The mean, total variance, and biological
#' component (the residual from the fitted mean-variance trend) of each feature
#' are stored in \code{rowData}. The method name and \code{rowData} column
#' names are kept from the scran implementation.
#' 
#' Also selects a specified number of top HVGs and store the logical selection 
#' in \code{rowData}. 
#' @param inSCE A \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment} object
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
    fit <- scrapper::modelGeneVariances(assay(inSCE, useAssay))$statistics
    rowData(inSCE)$scran_modelGeneVar_mean <- fit$means
    rowData(inSCE)$scran_modelGeneVar_totalVariance <- fit$variances
    rowData(inSCE)$scran_modelGeneVar_bio <- fit$residuals
    metadata(inSCE)$sctk$runFeatureSelection$modelGeneVar <- 
        list(useAssay = useAssay,
             rowData = c("scran_modelGeneVar_mean", 
                         "scran_modelGeneVar_totalVariance",
                         "scran_modelGeneVar_bio"))
    return(inSCE)
}
