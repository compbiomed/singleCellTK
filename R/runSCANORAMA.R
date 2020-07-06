#' Apply the mutual nearest neighbors (MNN) batch effect correction method to
#' SingleCellExperiment object
#'
#' SCANORAMA is analogous to computer vision algorithms for panorama stitching
#' that identify images with overlapping content and merge these into a larger
#' panorama.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name
#' of the assay requiring batch correction in "inSCE", should exist in
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the
#' field of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"SCANORAMA"`. The name for the
#' corrected full-sized expression matrix.
#' @param SIGMA numeric, default `15`. Algorithmic parameter, correction
#' smoothing parameter on Gaussian kernel.
#' @param ALPHA numeric, default `0.1`. Algorithmic parameter, alignment score
#' minimum cutoff.
#' @param KNN integer, default `20L`. Algorithmic parameter, number of nearest
#' neighbors to use for matching.
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated
#' with corrected full-sized expression matrix.
#' @export
#' @references Brian Hie et al, 2019
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCANORAMA(sceBatches)
#' }
runSCANORAMA <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                         assayName = 'SCANORAMA', SIGMA = 15, ALPHA = 0.1,
                         KNN = 20L){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!reticulate::py_module_available(module = "scanorama")){
        warning("Cannot find python module 'scanorama', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, scanorama can be installed on the local machine",
            "with pip (e.g. pip install scanorama) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
        return(inSCE)
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    assayName <- gsub(' ', '_', assayName)
    adata <- .sce2adata(inSCE, useAssay)
    py <- reticulate::py
    py$adata <- adata
    py$batch <- batch
    py$sigma <- SIGMA
    py$alpha <- ALPHA
    py$KNN <- as.integer(KNN)
    reticulate::py_run_string('
import numpy as np
import scanorama
batches = list(set(adata.obs[batch]))
adatas = [adata[adata.obs[batch] == b,] for b in batches]
datasets_full = [a.X for a in adatas]
genes_list = [a.var_names.to_list() for a in adatas]
corrected, genes = scanorama.correct(datasets_full, genes_list, sigma=sigma,
                                     alpha=alpha, knn=KNN)
corrected = [m.toarray() for m in corrected]
cellOrders = [adata.obs[batch] == b for b in batches]
integrated = np.zeros(adata.shape)
integrated[:,:] = np.NAN
for i in range(len(batches)):
    integrated[cellOrders[i],] = corrected[i]
geneidx = {gene: idx for idx, gene in enumerate(genes)}
orderIdx = []
for gene in adata.var_names:
    orderIdx.append(geneidx[gene])
integrated = integrated[:, orderIdx]
', convert = FALSE)
    mat <- t(py$integrated)
    rownames(mat) <- rownames(inSCE)
    colnames(mat) <- colnames(inSCE)
    SummarizedExperiment::assay(inSCE, assayName) <- mat
    return(inSCE)
}
