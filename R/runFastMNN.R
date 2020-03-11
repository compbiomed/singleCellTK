#' Apply a fast version of the mutual nearest neighbors (MNN) batch effect 
#' correction method to SingleCellExperiment object
#' 
#' `fastMNN` is a variant of the `mnnCorrect` function, modified for speed and 
#' more robust performance. For introduction of MNN, see `runMNNCorrect`
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`. Alternatively, see `pcInput` parameter.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"fastMNN"`. The name for the 
#' corrected low-dimensional representation.
#' @param pcInput bool, default `FALSE`. Whether to use a reduceDim matrix for 
#' batch effect correction. If TRUE, `useAssay` should exist in 
#' `reduceDimNames(inSCE)`. Note that more dimensions in precomputed reducedDim
#' have shown better correction results.
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation.
#' @export
#' @references Lun ATL, et al., 2016
#' @examples  
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runFastMNN(sceBatches, useAssay = 'PCA', pcInput = TRUE)
runFastMNN <- function(inSCE, useAssay = "logcounts", reducedDimName = "MNN", 
                       batch = 'batch', pcInput = FALSE){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(pcInput){
        if(!useAssay %in% SingleCellExperiment::reducedDimNames(inSCE)) {
            stop(paste("\"useAssay\" (reducedDim) name: ", useAssay, " not found."))
        }
    } else {
        if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
            stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
        }
    }
    
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch name:", batch, "not found."))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    
    ## Run algorithm
    batches <- list()
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    batchFactor <- as.factor(batchCol)
    
    if(pcInput){
        mat <- SingleCellExperiment::reducedDim(inSCE, useAssay)
        redMNN <- batchelor::reducedMNN(mat, batch = batchFactor)
        newRedDim <- redMNN$corrected
    } else {
        mat <- SummarizedExperiment::assay(inSCE, useAssay)
        mnnSCE <- batchelor::fastMNN(mat, batch = batchFactor)
        newRedDim <- SingleCellExperiment::reducedDim(mnnSCE, 'corrected')
    }
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- newRedDim
    return(inSCE)
}
