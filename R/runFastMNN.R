#' Apply a fast version of the mutual nearest neighbors (MNN) batch effect 
#' correction method to SingleCellExperiment object
#' 
#' `fastMNN` is a variant of the `mnnCorrect` function, modified for speed and 
#' more robust performance. For introduction of MNN, see `runMNNCorrect`
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`. Alternatively, see `pcInput` parameter.
#' @param batchKey character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"fastMNN"`. The name for the 
#' corrected low-dimensional representation.
#' @param pcInput bool, default `FALSE`. Whether to use a reduceDim matrix for 
#' batch effect correction. If TRUE, `exprs` should exist in 
#' `reduceDimNames(inSCE)`. Note that more dimensions in precomputed reducedDim
#' have shown better correction results.
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation.
#' @export
#' @references Lun ATL, et al., 2016
#' @examples  
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runFastMNN(sceBatches, exprs = 'PCA', pcInput = TRUE)
runFastMNN <- function(inSCE, exprs = "logcounts", reducedDimName = "MNN", 
                       batchKey = 'batch', pcInput = FALSE){
    ## Input check
    if(!class(inSCE) == "SingleCellExperiment" && 
       !class(inSCE) == "SCtkExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(pcInput){
        if(!exprs %in% SingleCellExperiment::reducedDimNames(inSCE)) {
            stop(paste("\"exprs\" (assay) name: ", exprs, " not found."))
        }
    } else {
        if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
            stop(paste("\"exprs\" (assay) name: ", exprs, " not found."))
        }
    }
    
    if(!batchKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batchKey name:", batchKey, "not found."))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    
    ## Run algorithm
    # Extract each batch
    if(pcInput){
        mat <- SingleCellExperiment::reducedDim(inSCE, exprs)
    } else {
        mat <- SummarizedExperiment::assay(inSCE, exprs)
    }
    batches <- list()
    batchCol <- SummarizedExperiment::colData(inSCE)[[batchKey]]
    batchIndicator <- unique(batchCol)
    for(i in batchIndicator){
        if(pcInput){
            batches[[i]] <- mat[batchCol == i,]
        } else {
            batches[[i]] <- mat[,batchCol == i]
        }
        
    }
    if(pcInput){
        mnn <- batchelor::fastMNN(batches[[1]], 
                              batches[[2]], 
                              pc.input = TRUE)    
    } else {
        mnn <- batchelor::fastMNN(batches[[1]], 
                              batches[[2]])
    }
    rn <- rownames(mat)
    newMat <- mnn$corrected[rn,]
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- newMat
    return(inSCE)
}
