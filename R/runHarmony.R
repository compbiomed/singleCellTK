#' Apply Harmony batch effect correction method to SingleCellExperiment object
#' 
#' Harmony is an algorithm that projects cells into a shared embedding in which 
#' cells group by cell type rather than dataset-specific conditions. 
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"HARMONY"`. The name for the 
#' corrected low-dimensional representation.
#' @param pcInput bool, default `FALSE`. Whether to do correction on a low-dim
#' matrix. If `TRUE`, will look for `useAssay` in `reducedDim(inSCE)`.
#' @param nComponents integer, default `20`. Number of principle components or 
#' dimensionality to generate in the resulting reducedDim. If `pcInput`, will
#' directly follow the dimensionality of the specified reducedDim.
#' @param nIter integer, default `10`. The max number of iterations to perform. 
#' @param theta Numeric, default 5. Diversity clustering penalty parameter, 
#' Larger value results in more diverse clusters.
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation.
#' @export
#' @references Ilya Korsunsky, et al., 2019
#' @examples  
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runHarmony(sceBatches, nComponents = 10)
runHarmony <- function(inSCE, useAssay = "logcounts", pcInput = FALSE, 
                       batch = "batch", reducedDimName = "HARMONY", 
                       nComponents = 50, theta = 5, nIter = 10){
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
    
    ## Run algorithm
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    if(pcInput){
        mat <- SingleCellExperiment::reducedDim(inSCE, useAssay)
    } else{
        sceTmp <- scater::runPCA(inSCE, useAssay_values = useAssay, 
                                 ncomponents = nComponents)
        mat <- SingleCellExperiment::reducedDim(sceTmp, 'PCA')
    }
    h <- harmony::HarmonyMatrix(mat, batchCol, do_pca = FALSE, 
                                theta = theta, max.iter.harmony = nIter)
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- h
    return(inSCE)
}
