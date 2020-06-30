#' Apply ZINBWaVE Batch effect correction method to SingleCellExperiment object
#'
#' A general and flexible zero-inflated negative binomial model that can be
#' used to provide a low-dimensional representations of scRNAseq data. The
#' model accounts for zero inflation (dropouts), over-dispersion, and the count
#' nature of the data. The model also accounts for the difference in library
#' sizes and optionally for batch effects and/or other covariates.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Note that ZINBWaVE works for counts (integer) input rather
#' than logcounts that other methods prefer. Default \code{"counts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param nHVG An integer. Number of highly variable genes to use when fitting
#' the model. Default \code{1000L}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. Default \code{50L}.
#' @param nIter An integer, The max number of iterations to perform. Default
#' \code{10L}.
#' @param epsilon An integer. Algorithmic parameter. Empirically, a high epsilon
#' is often required to obtained a good low-level representation. Default
#' \code{1000L}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"zinbwave"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Pollen, Alex A et al., 2014
#' @examples
#' \dontrun{
#'     data('sceBatches', package = 'singleCellTK')
#'     sceCorr <- runZINBWaVE(sceBatches, nIter = 5)
#' }
runZINBWaVE <- function(inSCE, useAssay = 'counts', batch = 'batch',
                        nHVG = 1000L, nComponents = 50L, epsilon = 1000,
                        nIter = 10L, reducedDimName = 'zinbwave'){
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
    reducedDimName <- gsub(' ', '_', reducedDimName)
    nHVG <- as.integer(nHVG)
    nComponents <- as.integer(nComponents)
    epsilon <- as.integer(epsilon)
    nIter <- as.integer(nIter)
    # Run algorithm
    ##ZINBWaVE tutorial style of HVG selection
    if(nHVG < nrow(inSCE)){
        logAssay <- log1p(SummarizedExperiment::assay(inSCE, useAssay))
        vars <- matrixStats::rowVars(logAssay)
        names(vars) <- rownames(inSCE)
        vars <- sort(vars, decreasing = TRUE)
        tmpSCE <- inSCE[names(vars)[1:nHVG],]
    }
    epsilon <- min(nrow(inSCE), epsilon)
    print('start!')
    tmpSCE <- zinbwave::zinbwave(tmpSCE, K = nComponents, epsilon = epsilon,
                                 which_assay = useAssay,
                                 X = paste('~', batch, sep = ''),
                                 maxiter.optimize = nIter, verbose = TRUE)
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <-
        SingleCellExperiment::reducedDim(tmpSCE, 'zinbwave')
    return(inSCE)
}
