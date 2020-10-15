#' Apply a fast version of the mutual nearest neighbors (MNN) batch effect
#' correction method to SingleCellExperiment object
#'
#' fastMNN is a variant of the classic MNN method, modified for speed and more
#' robust performance. For introduction of MNN, see \code{\link{runMNNCorrect}}.
#' @param inSCE  inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}. Alternatively, see
#' \code{pcInput} parameter.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"fastMNN"}.
#' @param pcInput A logical scalar. Whether to use a low-dimension matrix for
#' batch effect correction. If \code{TRUE}, \code{useAssay} will be searched
#' from \code{reducedDimNames(inSCE)}. Default \code{FALSE}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Lun ATL, et al., 2016
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#' sceCorr <- runFastMNN(sceBatches, useAssay = 'logcounts', pcInput = FALSE)
#' }
runFastMNN <- function(inSCE, useAssay = "logcounts",
                       reducedDimName = "fastMNN", batch = 'batch',
                       pcInput = FALSE){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(isTRUE(pcInput)){
        if(!(useAssay %in% SingleCellExperiment::reducedDimNames(inSCE))) {
            stop(paste("\"useAssay\" (reducedDim) name: ", useAssay, " not found."))
        }
    } else {
        if(!(useAssay %in% SummarizedExperiment::assayNames(inSCE))) {
            stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
        }
    }

    if(!(batch %in% names(SummarizedExperiment::colData(inSCE)))){
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
