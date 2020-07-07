#' Apply scGen batch effect correction method to SingleCellExperiment object
#'
#' scGen is a generative model to predict single-cell perturbation response
#' across cell types, studies and species. It works by combining variational
#' autoencoders and latent space vector arithmetics for high-dimensional single-
#' cell gene expression data.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param cellType A single character. A string indicating a field in
#' \code{colData(inSCE)} that defines different cell types. Default
#' \code{'cell_type'}.
#' @param nEpochs An integer. Algorithmic parameter, the number of epochs to
#' iterate and optimize network weights. Default \code{50L}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"SCGEN"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Lotfollahi, Mohammad et al., 2019
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCGEN(sceBatches)
#' }
runSCGEN <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                     cellType = "cell_type", nEpochs = 50L,
                     assayName = 'SCGEN'){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if (!reticulate::py_module_available(module = "scgen")) {
        warning("Cannot find python module 'scgen', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, scgen can be installed on the local machine",
            "with pip (e.g. pip install scgen) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
        return(inSCE)
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    if(!cellType %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"cellType\" name:", cellType, "not found"))
    }
    assayName <- gsub(' ', '_', assayName)
    nEpochs <- as.integer(nEpochs)

    ## Run algorithm
    adata <- .sce2adata(inSCE, useAssay = useAssay)
    network = scgen$VAEArith(x_dimension = adata$n_vars)
    network$train(train_data = adata, n_epochs = nEpochs)
    corrAdata <- scgen$batch_removal(network, adata, batch_key = batch,
                                     cell_label_key = cellType)
    corrMat <- t(corrAdata$X)
    SummarizedExperiment::assay(inSCE, assayName) <- corrMat
    return(inSCE)
}
