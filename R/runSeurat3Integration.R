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
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name
#' of the assay requiring batch correction in "inSCE", should exist in
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the
#' field of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"Seurat3Int"`. The name for the
#' corrected full-sized expression matrix. If the number of features returned
#' is smaller the number of total feature, the returned matrix will be saved in
#' `reducedDim(inSCE, assayName)`; if equal, `assay(inSCE, assayName)`.
#' @param nAnchors integer, default `nrow(inSCE)`. The number of features to
#' anchor, and also the final dimensionality of the integrated matrix. Thus
#' default value turns to produce full-sized assay.
#' @param verbose bool, default `TRUE`. Whether to show detail information of
#' the process.
#' @export
#' @references Stuart et al. 2019
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSeurat3Integration(sceBatches, nAnchors = 100)
#' }
runSeurat3Integration <- function(inSCE, useAssay = 'logcounts',
                                  batch = 'batch',
                                  assayName = "Seurat3Int",
                                  nAnchors = nrow(inSCE), verbose = TRUE){

    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    assayName <- gsub(' ', '_', assayName)

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
    assayList[[assayName]] <- IntMat
    AE <- SingleCellExperiment::SingleCellExperiment(assay = assayList)
    SingleCellExperiment::altExp(inSCE, assayName) <- AE
    return(inSCE)
}
