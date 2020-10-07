#' Apply Harmony batch effect correction method to SingleCellExperiment object
#'
#' Harmony is an algorithm that projects cells into a shared embedding in which
#' cells group by cell type rather than dataset-specific conditions.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"HARMONY"}.
#' @param pcInput A logical scalar. Whether to use a low-dimension matrix for
#' batch effect correction. If \code{TRUE}, \code{useAssay} will be searched
#' from \code{reducedDimNames(inSCE)}. Default \code{FALSE}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. If \code{pcInput} is set
#' to \code{TRUE}, the output dimension will follow the low-dimension matrix,
#' so this argument will be ignored. Default \code{50L}.
#' @param nIter An integer. The max number of iterations to perform. Default
#' \code{10L}.
#' @param theta A Numeric scalar. Diversity clustering penalty parameter,
#' Larger value results in more diverse clusters. Default \code{5}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Ilya Korsunsky, et al., 2019
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runHarmony(sceBatches, nComponents = 10L)
runHarmony <- function(inSCE, useAssay = "logcounts", pcInput = FALSE,
                       batch = "batch", reducedDimName = "HARMONY",
                       nComponents = 50L, theta = 5, nIter = 10L){
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
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    nComponents <- as.integer(nComponents)
    nIter <- as.integer(nIter)
    ## Run algorithm
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    if(pcInput){
        mat <- SingleCellExperiment::reducedDim(inSCE, useAssay)
    } else{
        sceTmp <- scater::runPCA(inSCE, exprs_values = useAssay,
                                 ncomponents = nComponents)
        mat <- SingleCellExperiment::reducedDim(sceTmp, 'PCA')
    }
    h <- harmony::HarmonyMatrix(mat, batchCol, do_pca = FALSE,
                                theta = theta, max.iter.harmony = nIter)
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- h
    return(inSCE)
}
