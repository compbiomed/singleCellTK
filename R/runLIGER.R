#' Apply LIGER batch effect correction method to SingleCellExperiment object
#' 
#' LIGER relies on integrative non-negative matrix factorization to identify 
#' shared and dataset-specific factors.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"LIGER"`. The name for the 
#' corrected low-dimensional representation.
#' @param nComponents integer, default `20L`. Number of principle components or 
#' dimensionality (factors, for this algorithm) to generate in the resulting 
#' reducedDim.
#' @param lambda numeric, default `5.0`. Algorithmic parameter, the penalty 
#' parameter which limits the dataset-specific component of the factorization.
#' @param resolution numeric, default `1.0`. Algorithmic paramter, the 
#' clustering resolution, increasing this increases the number of communities 
#' detected.
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation.
#' @export
#' @references Joshua Welch, et al., 2018
#' @examples  
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runLIGER(sceBatches)
#' }
runLIGER <- function(inSCE, exprs = 'logcounts', batch = 'batch', 
                     reducedDimName = 'LIGER', nComponents = 20L, lambda = 5.0, 
                     resolution = 1.0){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found"))
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
            SummarizedExperiment::assay(inSCE, exprs)[,batchCol == b]
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
