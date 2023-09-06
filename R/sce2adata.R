#' Coverts SingleCellExperiment object from R to anndata.AnnData object in
#' Python
#'
#' The AnnData object here can be saved to .h5ad file and read into Python
#' interactive console. Mostly used senario is when you want to apply
#' reticulated Python function, which only works with an anndata.AnnData object.
#' @param SCE A SingleCellExperiment object.
#' @param useAssay Character, default `"counts"`. The name of assay of
#' interests that will be set as the primary matrix of the output AnnData.
#' Available options can be listed by `assayNames(SCE)`. Thee primary matrix
#' will be saved in `adata$X`, Other assays will be stored in `adata$obsm`
#' together with the low-dimension representations (for now).
#' @return A Python anndata.AnnData object
#' @noRd
.sce2adata <- function(SCE, useAssay = 'counts') {

    # TODO: use zellkonverter in the future, temporary fix for now
    # this is how we used to do it until we started running into problems with the getters and setters
    # in the future, this function might be depreciated altogether since it is only called internally
    # and we will use zellkonverter::writeH5AD directly from an SCE object

    # # Transfer SCE object back to AnnData
    # # Argument check first
    # stopifnot(inherits(SCE, "SingleCellExperiment"))

    # # Extract information that correspond to AnnData structure
    # #X <- as.matrix(t(SummarizedExperiment::assay(SCE, useAssay)))
    # # Sparse matrix conversion supported now, commenting the line above.
    # X <- t(SummarizedExperiment::assay(SCE, useAssay))
    # AnnData <- sc$AnnData(X = X)
    # obs <- as.data.frame(SummarizedExperiment::colData(SCE))
    # if(length(obs) > 0){
    #     AnnData$obs = obs
    # } else {
    #     AnnData$obs_names <- colnames(SCE)
    # }
    # var <- as.data.frame(SummarizedExperiment::rowData(SCE))
    # if(length(var) > 0){
    #     AnnData$var = var
    # } else {
    #     AnnData$var_names <- rownames(SCE)
    # }
    # # uns  <- S4Vectors::metadata(SCE)
    # # if(length(uns) > 0){ AnnData$uns <- uns }
    # obsmNames <- SingleCellExperiment::reducedDimNames(SCE)
    # if(length(obsmNames) > 0){
    #     for (i in seq_along(obsmNames)) {
    #         AnnData$obsm$'__setitem__'(obsmNames[i], SingleCellExperiment::reducedDim(SCE, obsmNames[i]))
    #     }
    # }

    # # Furthermore, the other assays will for now also be saved to .layers
    # allAssayNames <- SummarizedExperiment::assayNames(SCE)
    # for (i in seq_along(allAssayNames)) {
    #     oneName <- allAssayNames[i]
    #     if (!oneName == useAssay) {
    #         AnnData$layers$'__setitem__'(oneName, as.matrix(t(SummarizedExperiment::assay(SCE, oneName))))
    #     }
    # }

    X <- t(SummarizedExperiment::assay(SCE, useAssay))
    AnnData <- sc$AnnData(X = X)
    obs <- as.data.frame(SummarizedExperiment::colData(SCE))
    if(length(obs) > 0){
        AnnData$obs = obs
    } else {
        AnnData$obs_names <- colnames(SCE)
    }
    var <- as.data.frame(SummarizedExperiment::rowData(SCE))
    if(length(var) > 0){
        AnnData$var = var
    } else {
        AnnData$var_names <- rownames(SCE)
    }
    # previously commented out by someone else
    # uns  <- S4Vectors::metadata(SCE)
    # if(length(uns) > 0){ AnnData$uns <- uns }
    obsmNames <- SingleCellExperiment::reducedDimNames(SCE)
    # new method: make a list of dataframes, which Python can process for H5AD construction
    # initialize empty list for obsm and layers, which are the things we need to do
    obsm <- list()

    if(length(obsmNames) > 0){
        for (i in seq_along(obsmNames)) {
            obsm[[obsmNames[i]]] <- data.frame(SingleCellExperiment::reducedDim(SCE, obsmNames[i]))
            #previously using the dunder setter, which doesn't always work
            #AnnData$obsm$'__setitem__'(obsmNames[i], SingleCellExperiment::reducedDim(SCE, obsmNames[i]))
        }
    }
    
    AnnData$obsm <- obsm

    # Furthermore, the other assays will for now also be saved to .layers
    allAssayNames <- SummarizedExperiment::assayNames(SCE)
    layers <- list()
    for (i in seq_along(allAssayNames)) {
        oneName <- allAssayNames[i]
        if (!oneName == useAssay) {
            layers[[allAssayNames[i]]] <- data.frame(as.matrix(t(SummarizedExperiment::assay(SCE, oneName))))
            # this is the way we used to do it, through a Pythonic dunder setter
            #AnnData$layers$'__setitem__'(oneName, as.matrix(t(SummarizedExperiment::assay(SCE, oneName))))
        }
    }
    AnnData$layers <- layers
    return(AnnData)
}
