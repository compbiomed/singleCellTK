#' Apply Limma's batch effect correction method to SingleCellExperiment object
#' 
#' Limma's batch effect removal function fits a linear model to the data, then 
#' removes the component due to the batch effects.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batchKey character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"LIMMA"`. The name for the corrected 
#' full-sized expression matrix.
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated 
#' with corrected full-sized expression matrix.
#' @export
#' @references Gordon K Smyth, et al., 2003
#' @examples  
#' data('sceBatches', package = 'singleCellTK')
#' sceBatches
#' ## class: SingleCellExperiment 
#' ## dim: 27610 1820 
#' ## metadata(0):
#' ## assays(3): normcounts logcounts
#' ## rownames(27610): GCG MALAT1 ... LOC102724004 LOC102724238
#' ## rowData names(0):
#' ## colnames(1820): reads.12732 reads.12733 ... Sample_1598 Sample_1600
#' ## colData names(2): cell_type1 batch
#' ## reducedDimNames(5): PCA
#' ## spikeNames(0):
#' sceCorr <- runLimmaBC(sceBatches)
runLimmaBC <- function(inSCE, exprs = "logcounts", assayName = "LIMMA", 
                       batchKey = "batch") {
    ## Input check
    if(!class(inSCE) == "SingleCellExperiment" && 
       !class(inSCE) == "SCtkExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found."))
    }
    if(!batchKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batchKey\" name:", batchKey, "not found."))
    }
    
    assayName <- gsub(' ', '_', assayName)
    
    ## Run algorithm
    ## One more check for the batch names
    batchCol <- SummarizedExperiment::colData(inSCE)[[batchKey]]
    mat <- SummarizedExperiment::assay(inSCE, exprs)
    newMat <- limma::removeBatchEffect(mat, batch = batchCol)
    SummarizedExperiment::assay(inSCE, assayName) <- newMat
    return(inSCE)
}
