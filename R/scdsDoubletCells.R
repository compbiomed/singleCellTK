.runDoubletScds <- function(sce,
    sampleColname,
    ...,
    method) {

    doubletScore <- rep(NA, ncol(sce))
    samples <- unique(SummarizedExperiment::colData(sce)[[sampleColname]])

    if ("DelayedMatrix" %in% class(SummarizedExperiment::assay(sce,
        i = "counts"))) {
        SummarizedExperiment::assay(sce, i = "counts") <-
            as.matrix(SummarizedExperiment::assay(sce, i = "counts"))
    }

    for (sample in samples) {
        sceSampleInd <- which(SummarizedExperiment::colData(sce)
            [[sampleColname]] == sample)
        
        if (!method %in% c("cxds", "bcds", "hybrid")){
            warning("Doublet detection methods must be one of the following: ",
                "'cxds', 'bcds' or 'hybrid'")
        } else if (method == "cxds") {
            sceSubset <- scds::cxds(sce[, sceSampleInd], 
                ..., 
                retRes=TRUE)
            score <- sceSubset$cxds_score
        } else if (method == "bcds") {
            sceSubset <- scds::bcds(sce[, sceSampleInd], 
                ..., 
                verb=TRUE)
            score <- sceSubset$bcds_score
        } else if (method == "hybrid") {
            sceSubset <- scds::cxds_bcds_hybrid(sce[, sceSampleInd], ...)
            score <- sceSubset$hybrid_score
        }
        
        doubletScore[sceSampleInd] <- score
    }

    SummarizedExperiment::colData(sce)$scds_doublet_score <- doubletScore
    return(sce)
}


#' @title Detect doublet cells using scds.
#' @description A wrapper function for \link[scds]{cxds}, \link[scds]{bcds} and 
#'  \link[scds]{hybrid}. Identify potential doublet cells based on gene 
#'  co-expression or simulations of putative doublet expression profiles. 
#'  Generate a doublet score for each cell.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object. Must
#'  contain a count matrix.
#' @param sampleColname Character. The column name which specifies the sample
#'  origin in the \link[SummarizedExperiment]{colData} of the provided
#'  \code{sce} object. Default "sample".
#' @param ... Additional arguments to pass to scds doublet detection functions.
#' @details This function is a wrapper function for scds doublet detection
#'  functions. It assumes the provided \code{sce} was generated using the import
#'  functions in this package so there
#'  will be a "sample" column in the \link[SummarizedExperiment]{colData}
#'  specifying the sample origin of each cell. The column name can be specified
#'  using the \code{sampleColname} argument.
#'  \link{runDoubletScds} runs \link[scds]{cxds}, \link[scds]{bcds} or 
#'  \link[scds]{hybrid} for each sample within \code{sce} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link[SummarizedExperiment]{colData} of \code{sce}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  'scds_doublet_score' column added to the
#'  \link[SummarizedExperiment]{colData} slot.
#' @references Bais AS, Kostka D (2019). “scds: Computational Annotation of 
#'  Doublets in Single Cell RNA Sequending Data.” bioRxiv. doi: 10.1101/564021.
#'  \url{https://github.com/kostkalab/scds}
 

#' @seealso \link[scds]{cxds}, \link[scds]{bcds}, \link[scds]{hybrid}
#' @examples
#' data(doubletSceExample, package = "singleCellTK")
#' sce <- runDoubletScds(doubletSceExample)
#' @export
#' @import scds
runDoubletScds <- function(sce,
    sampleColname = "sample",
    ...,
    method = "hybrid") {

    sce <- .runDoubletScds(sce = sce,
        sampleColname = sampleColname,
        ...,
        method = method)
    return(sce)
}