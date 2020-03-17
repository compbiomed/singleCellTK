#' Apply Limma's batch effect correction method to SingleCellExperiment object
#' 
#' Limma's batch effect removal function fits a linear model to the data, then 
#' removes the component due to the batch effects.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"LIMMA"`. The name for the corrected 
#' full-sized expression matrix.
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated 
#' with corrected full-sized expression matrix.
#' @export
#' @references Gordon K Smyth, et al., 2003
#' @examples  
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runLimmaBC(sceBatches)
runLimmaBC <- function(inSCE, useAssay = "logcounts", assayName = "LIMMA", 
                       batch = "batch") {
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found."))
    }
    
    assayName <- gsub(' ', '_', assayName)
    
    ## Run algorithm
    ## One more check for the batch names
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    mat <- SummarizedExperiment::assay(inSCE, useAssay)
    newMat <- limma::removeBatchEffect(mat, batch = batchCol)
    SummarizedExperiment::assay(inSCE, assayName) <- newMat
    return(inSCE)
}
