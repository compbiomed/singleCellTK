#' Apply BBKNN batch effect correction method to SingleCellExperiment object
#'
#' BBKNN, an extremely fast graph-based data integration algorithm. It modifies
#' the neighbourhood construction step to produce a graph that is balanced
#' across all batches of the data.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"BBKNN"}.
#' @param nComponents An integer. Number of principle components or the
#' dimensionality, adopted in the pre-PCA-computation step, the BBKNN step (for
#' how many PCs the algorithm takes into account), and the final UMAP
#' combination step where the value represent the dimensionality of the updated
#' reducedDim. Default \code{50L}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Krzysztof Polanski et al., 2020
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#' sceCorr <- runBBKNN(sceBatches)
#' }
runBBKNN <-function(inSCE, useAssay = 'logcounts', batch = 'batch',
                    reducedDimName = 'BBKNN', nComponents = 50L){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!reticulate::py_module_available(module = "bbknn")){
        warning("Cannot find python module 'bbknn', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, bbknn can be installed on the local machine",
            "with pip (e.g. pip install bbknn) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
        return(inSCE)
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    nComponents <- as.integer(nComponents)
    ## Run algorithm
    adata <- .sce2adata(inSCE, useAssay = useAssay)
    sc$tl$pca(adata, n_comps = nComponents)
    bbknn$bbknn(adata, batch_key = batch, n_pcs = nComponents)
    sc$tl$umap(adata, n_components = nComponents)
    bbknnUmap <- adata$obsm[["X_umap"]]
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- bbknnUmap
    return(inSCE)
}
