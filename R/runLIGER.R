#' Apply LIGER batch effect correction method to SingleCellExperiment object
#'
#' LIGER relies on integrative non-negative matrix factorization to identify
#' shared and dataset-specific factors.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"LIGER"}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. Default \code{20L}.
#' @param lambda A numeric scalar. Algorithmic parameter, the penalty
#' parameter which limits the dataset-specific component of the factorization.
#' Default \code{5.0}.
#' @param resolution A numeric scalar. Algorithmic paramter, the clustering
#' resolution, increasing this increases the number of communities detected.
#' Default \code{1.0}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Joshua Welch, et al., 2018
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runLIGER(sceBatches)
#' }
runLIGER <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                     reducedDimName = 'LIGER', nComponents = 20L, lambda = 5.0,
                     resolution = 1.0){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)

    ## Run algorithm
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    batches <- unique(batchCol)
    nBatch <- length(batches)
    batchMatrices <- list()
    for(i in 1:nBatch){
        b <- batches[i]
        batchMatrices[[b]] <-
            SummarizedExperiment::assay(inSCE, useAssay)[,batchCol == b]
    }
    ligerex <- liger::createLiger(batchMatrices)
    ligerex <- liger::normalize(ligerex)
    ligerex <- liger::selectGenes(ligerex)
    ligerex <- liger::scaleNotCenter(ligerex)
    ligerex <- liger::optimizeALS(ligerex, k = nComponents, lambda = lambda,
                                  resolution = resolution)
    ligerex <- liger::quantile_norm(ligerex)
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- ligerex@H.norm
    return(inSCE)
}
