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
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batchKey character, default `"batch"`. A string indicating the field 
#' of `colData(inSCE)` that defines different batches.
#' @param cellTypeKey character, default `"cell_type"`. A string indicating the 
#' field of `colData(inSCE)` that defines different cell types.
#' @param assayName character, default `"SCGEN"`. The name for the corrected 
#' full-sized expression matrix.
#' @param nEpochs integer, default `100L`. Algorithmic parameter, number of 
#' epochs to iterate and optimize network weights. 
#' @export
#' @references Lotfollahi, Mohammad et al., 2019
#' @examples 
#' data('sceBatches', package = 'singleCellTK')
#' sceBatches
#' ## class: SingleCellExperiment 
#' ## dim: 27610 1820 
#' ## metadata(0):
#' ## assays(3): normcounts logcounts
#' ## rownames(27610): GCG MALAT1 ... LOC102724004 LOC102724238
#' ## rowData names(0):
#' ## colnames(1820): reads.12732 reads.12733 ... Sample_1598 Sample_1600
#' ## colData names(2): cell_type1 batch
#' ## reducedDimNames(5): PCA
#' ## spikeNames(0):
#' sceCorr <- runSCGEN(sceBatches)
runSCGEN <- function(inSCE, exprs = 'logcounts', batchKey = 'batch', 
                     cellTypeKey = "cell_type", assayName = 'SCGEN', 
                     nEpochs = 50L){
    ## Input check
    if(!class(inSCE) == "SingleCellExperiment" 
       && !class(inSCE) == "SCtkExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!batchKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batchKey\" name:", batchKey, "not found"))
    }
    if(!cellTypeKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"cellTypeKey\" name:", batchKey, "not found"))
    }
    assayName <- gsub(' ', '_', assayName)
    nEpochs <- as.integer(nEpochs)

    ## Run algorithm
    adata <- sce2adata(inSCE, mainAssay = exprs)
    network = scgen$VAEArith(x_dimension = adata$n_vars)
    network$train(train_data = adata, n_epochs = nEpochs)
    corrAdata <- scgen$batch_removal(network, adata, batch_key = batchKey, 
                                     cell_label_key = cellTypeKey)
    corrMat <- t(corrAdata$X)
    SummarizedExperiment::assay(inSCE, assayName) <- corrMat
    return(inSCE)
}
