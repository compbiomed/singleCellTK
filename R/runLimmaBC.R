#' Apply Limma's batch effect correction method to SingleCellExperiment object
#'
#' Limma's batch effect removal function fits a linear model to the data, then
#' removes the component due to the batch effects.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"LIMMA"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
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
