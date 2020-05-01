#' Apply scGen batch effect correction method to SingleCellExperiment object
#' 
#' scGen is a generative model to predict single-cell perturbation response 
#' across cell types, studies and species. It works by combining variational 
#' autoencoders and latent space vector arithmetics for high-dimensional single-
#' cell gene expression data.
#' 
#' Result does not look fine for now. Time consuming also even it allocates 32 
#' cores.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the field 
#' of `colData(inSCE)` that defines different batches.
#' @param cellType character, default `"cell_type"`. A string indicating the 
#' field of `colData(inSCE)` that defines different cell types.
#' @param assayName character, default `"SCGEN"`. The name for the corrected 
#' full-sized expression matrix.
#' @param nEpochs integer, default `100L`. Algorithmic parameter, number of 
#' epochs to iterate and optimize network weights. 
#' @export
#' @references Lotfollahi, Mohammad et al., 2019
#' @examples 
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCGEN(sceBatches)
#' }
runSCGEN <- function(inSCE, useAssay = 'logcounts', batch = 'batch', 
                     cellType = "cell_type", assayName = 'SCGEN', 
                     nEpochs = 50L){
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
