#' Apply BBKNN batch effect correction method to SingleCellExperiment object
#' 
#' BBKNN, an extremely fast graph-based data integration algorithm. It modifies 
#' the neighbourhood construction step to produce a graph that is balanced 
#' across all batches of the data. 
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"BBKNN"`. The name for the 
#' corrected low-dimensional representation.
#' @param nComponents integer, default `50L`. Number of principle components or 
#' dimensionality, adopted in the pre-PCA-computation step, the BBKNN step (for 
#' how many PCs the algorithm takes into account), and the final UMAP 
#' combination step where the value represent the dimensionality of the updated 
#' reducedDim.
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation
#' @export
#' @references Krzysztof Polański et al., 2020
#' @examples  
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
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
    
    ## Run algorithm
    adata <- .sce2adata(inSCE, mainAssay = useAssay)
    sc$tl$pca(adata, n_comps = nComponents)
    bbknn$bbknn(adata, batch_key = batch, n_pcs = nComponents)
    sc$tl$umap(adata, n_components = nComponents)
    bbknnUmap <- adata$obsm[["X_umap"]]
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- bbknnUmap
    return(inSCE)
}
