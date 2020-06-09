#' Apply ZINBWaVE Batch effect correction mathod to SingleCellExperiment object
#' 
#' A general and flexible zero-inflated negative binomial model that can be 
#' used to provide a low-dimensional representations of scRNAseq data. The 
#' model accounts for zero inflation (dropouts), over-dispersion, and the count 
#' nature of the data. The model also accounts for the difference in library 
#' sizes and optionally for batch effects and/or other covariates.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"zinbwave"`. The name for the 
#' corrected low-dimensional representation.
#' @param nHVG integer, default `1000`. Number of highly variable genes to use 
#' when fitting the model
#' @param nComponents integer, default `50L`. Number of principle components or 
#' dimensionality to generate in the resulting reducedDim.
#' @param nIter integer, default `10`. The max number of iterations to perform. 
#' @param epsilon integer, default `1000`. Algorithmic parameter, by default, the 
#' epsilon parameter is set to the number of genes. We empirically found that a 
#' high epsilon is often required to obtained a good low-level representation.
#' @export
#' @references Pollen, Alex A et al., 2014
#' @examples 
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runZINBWaVE(sceBatches, nIter=5, filterParams=c(3, 3))
#' }
runZINBWaVE <- function(inSCE, useAssay = 'logcounts', batch = 'batch', 
                        reducedDimName = 'zinbwave', nHVG = 1000, 
                        nComponents = 50, epsilon = 1000, nIter = 10){
    #filterParams = NULL  <<< something told in tutorial but might be ignored
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch name:", batch, "not found."))
    }
    #if(!length(filterParams) == 2 || !class(filterParams) == "numeric" ||
    #   filterParams[1] < 0 || filterParams[2] < 0){
    #    stop("Invalid fillterParams, please follow the form c(a, b) with both a and b not negative. ")
    #}
    reducedDimName <- gsub(' ', '_', reducedDimName)
    
    # Run algorithm
    tmpMatrix <- round(SummarizedExperiment::assay(inSCE, useAssay))
    tmpSCE <- inSCE
    SummarizedExperiment::assay(tmpSCE, useAssay) <- tmpMatrix
    
    #if(!is.null(filterParams)){
    #    a <- filterParams[1]
    #    b <- filterParams[2]
    #    Filter <- rowSums(SummarizedExperiment::assay(tmpSCE) > a) > b
    #    tmpSCE <- tmpSCE[Filter,]
    #}
    
    ##ZINBWaVE tutorial style of HVG selection
    if(nHVG < nrow(inSCE)){
        logAssay <- log1p(SummarizedExperiment::assay(tmpSCE, useAssay))
        vars <- matrixStats::rowVars(logAssay)
        names(vars) <- rownames(tmpSCE)
        vars <- sort(vars, decreasing = TRUE)
        tmpSCE <- tmpSCE[names(vars)[1:nHVG],]
    }
    
    tmpSCE <- zinbwave::zinbwave(tmpSCE, K = nComponents, epsilon = epsilon, 
                                 which_assay = useAssay, 
                                  X = paste('~', batch, sep = ''), 
                                  maxiter.optimize=nIter, verbose = TRUE)
    reducedDim(inSCE, reducedDimName) <- reducedDim(tmpSCE, 'zinbwave')
    return(inSCE)
}
