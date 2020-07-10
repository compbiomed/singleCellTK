#' Apply the mutual nearest neighbors (MNN) batch effect correction method to
#' SingleCellExperiment object
#'
#' MNN is designed for batch correction of single-cell RNA-seq data where the
#' batches are partially confounded with biological conditions of interest. It
#' does so by identifying pairs of MNN in the high-dimensional log-expression
#' space. For each MNN pair, a pairwise correction vector is computed by
#' applying a Gaussian smoothing kernel with bandwidth `sigma`.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param k An integer. Specifies the number of nearest neighbours to
#' consider when defining MNN pairs. This should be interpreted as the minimum
#' frequency of each cell type or state in each batch. Larger values will
#' improve the precision of the correction by increasing the number of MNN
#' pairs, at the cost of reducing accuracy by allowing MNN pairs to form between
#' cells of different type. Default \code{20L}.
#' @param sigma A Numeric scalar. Specifies how much information is
#' shared between MNN pairs when computing the batch effect. Larger values will
#' share more information, approaching a global correction for all cells in the
#' same batch. Smaller values allow the correction to vary across cell types,
#' which may be more accurate but comes at the cost of precision. Default
#' \code{0.1}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"MNN"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Lun ATL, et al., 2016 & 2018
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runMNNCorrect(sceBatches)
runMNNCorrect <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                          assayName = 'MNN', k = 20L, sigma = 0.1){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch name:", batch, "not found."))
    }
    assayName <- gsub(' ', '_', assayName)
    k <- as.integer(k)

    ## Run algorithm
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    batchFactor <- as.factor(batchCol)
    mnnSCE <- batchelor::mnnCorrect(inSCE, batch = batchFactor,
                                    k = k, sigma = sigma)
    corrected <- SummarizedExperiment::assay(mnnSCE, 'corrected')
    SummarizedExperiment::assay(inSCE, assayName) <- corrected
    return(inSCE)
}
