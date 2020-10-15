#' Apply Integration batch effect correction method from Seurat v3 to
#' SingleCellExperiment object
#'
#' Can get either a full-sized corrected assay or a dimension reduced corrected
#' matrix.
#'
#' This method aims to first identify anchors between pairs of datasets, that
#' represent pairwise correspondences between individual cells in each
#' dataset, that is hypothesized to originate from the same cell state. These
#' anchors are then used to harmonize the datasets, or transfer information
#' from one dataset to another.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param nAnchors An integer. The number of features to anchor. The final
#' number of the corrected features depends on this value. Default
#' \code{nrow(inSCE)}.
#' @param verbose A logical scalar. Whether to show detail information of
#' the process. Default \code{TRUE}.
#' @param altExpName A single character. The name for the
#' \code{\link[SingleCellExperiment]{altExp}} that stores the corrected assay.
#' The name of this assay has the same name. Default \code{"Seurat3Int"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{altExp(inSCE, altExpName)} updated.
#' @export
#' @references Stuart et al. 2019
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#' sceCorr <- runSeurat3Integration(sceBatches, nAnchors = 100)
#' }
runSeurat3Integration <- function(inSCE, useAssay = 'logcounts',
                                  batch = 'batch',
                                  altExpName = "Seurat3Int",
                                  nAnchors = nrow(inSCE), verbose = TRUE){

    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
  altExpName <- gsub(' ', '_', altExpName)

    if(nAnchors > nrow(inSCE)){
        stop(paste("Specified nAnchors =", nAnchors,
                   "exceeded total number of available features"))
    }

    ## Run algorithm
    srtObj <- Seurat::as.Seurat(inSCE, counts = useAssay)
    batchSplit <- Seurat::SplitObject(srtObj, split.by = batch)
    nHVG <- nAnchors
    for (i in 1:length(batchSplit)){
        batchSplit[[i]] <- Seurat::NormalizeData(batchSplit[[i]],
                                                 verbose = verbose)
        batchSplit[[i]] <- Seurat::FindVariableFeatures(batchSplit[[i]],
                                                      selection.method = "vst",
                                                      nfeatures = nHVG,
                                                      verbose = verbose)
    }
    anchors <- Seurat::FindIntegrationAnchors(object.list = batchSplit,
                                              anchor.features = nAnchors,
                                              verbose = verbose)
    srtInt <- Seurat::IntegrateData(anchorset = anchors, verbose = verbose)
    IntMat <- as.matrix(Seurat::GetAssayData(srtInt, assay = 'integrated'))
    IntMat <- IntMat[,colnames(inSCE)]
    assayList <- list()
    assayList[[altExpName]] <- IntMat
    AE <- SingleCellExperiment::SingleCellExperiment(assay = assayList)
    SingleCellExperiment::altExp(inSCE, altExpName) <- AE
    return(inSCE)
}
