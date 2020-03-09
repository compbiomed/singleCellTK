#' Apply ZINBWaVE Batch effect correction mathod to SingleCellExperiment object
#' 
#' A general and flexible zero-inflated negative binomial model that can be 
#' used to provide a low-dimensional representations of scRNAseq data. The 
#' model accounts for zero inflation (dropouts), over-dispersion, and the count 
#' nature of the data. The model also accounts for the difference in library 
#' sizes and optionally for batch effects and/or other covariates.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"zinbwave"`. The name for the 
#' corrected low-dimensional representation.
#' @param nHVG integer, default `1000`. Number of highly variable genes to use 
#' when fitting the model
#' @param filterParams vector of two number, default `c(5, 5)`. In the form of 
#' c(a, b), and a filter like `rowSums(assay > a) > b` will be applied, meaning
#' that we only want the genes that has at least `a` counts in at least `b` 
#' cells. 
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
#' sceCorr <- runZINBWaVE(sceBatches, nIter=5)
#' }
runZINBWaVE <- function(inSCE, exprs = 'logcounts', batch = 'batch', 
                        reducedDimName = 'zinbwave', nHVG = 3000, 
                        filterParams = NULL, nComponents = 50, 
                        epsilon = 1000, nIter = 10){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found"))
    }
    if(!batch %in% names(colData(inSCE))){
        stop(paste("\"batch name:", batch, "not found."))
    }
    if(!length(filterParams) == 2 || !class(filterParams) == "numeric" ||
       filterParams[1] < 0 || filterParams[2] < 0){
        stop("Invalid fillterParams, please follow the form c(a, b) with both a and b not negative. ")
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    
    # Run algorithm
    tmpMatrix <- round(SummarizedExperiment::assay(inSCE, exprs))
    tmpSCE <- inSCE
    SummarizedExperiment::assay(tmpSCE, exprs) <- tmpMatrix
    if(!is.null(filterParams)){
        a <- filterParams[1]
        b <- filterParams[2]
        Filter <- rowSums(SummarizedExperiment::assay(tmpSCE) > a) > b
        tmpSCE <- tmpSCE[Filter,]
    }
    # Split the batches
    batches <- list()
    batchCol <- SummarizedExperiment::colData(tmpSCE)[[batch]]
    uniqBatch <- unique(batchCol)
    for(i in uniqBatch){
        batches[[i]] <- tmpSCE[, batchCol == i]
    }
    # Select HVG
    topVarGenesPerBatch <- list()
    for(i in uniqBatch){
        if(nrow(batches[[i]]) <= nHVG){
            topVarGenesPerBatch[[i]] <- 1:nrow(batches[[i]])
        } else {
          sce.var <- scran::modelGeneVar(batches[[i]])
          topVarGenesPerBatch[[i]] <- order(sce.var$bio,decreasing = TRUE)[seq(nHVG)]
#            mvTrend <- scran::trendVar(batches[[i]], use.spikes=FALSE, 
#                                       assay.type=exprs)
#            decomposeTrend <- scran::decomposeVar(batches[[i]], mvTrend, 
#                                                  assay.type=exprs)
#            topVarGenesPerBatch[[i]] <- order(decomposeTrend$bio, 
#                                              decreasing = TRUE)[1:nHVG]
        }    
    }
    selectedHVG <- BiocGenerics::Reduce(intersect, topVarGenesPerBatch)
    tmpSCE <- zinbwave::zinbwave(tmpSCE[selectedHVG,], K = nComponents, 
                                  epsilon = epsilon, which_assay=exprs, 
                                  X = paste('~', batch, sep = ''), 
                                  maxiter.optimize=nIter, verbose = TRUE)
    reducedDim(inSCE, reducedDimName) <- reducedDim(tmpSCE, 'zinbwave')
    return(inSCE)
}
