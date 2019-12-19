
.runDoubletCells <- function(sce,
    sampleColname,
    ...,
    assayType) {

    doubletScore <- rep(NA, ncol(sce))
    samples <- unique(SummarizedExperiment::colData(sce)[[sampleColname]])

    if ("DelayedMatrix" %in% class(SummarizedExperiment::assay(sce,
        i = assayType))) {
        SummarizedExperiment::assay(sce, i = assayType) <-
            as.matrix(SummarizedExperiment::assay(sce, i = assayType))
    }

    for (sample in samples) {
        sceSampleInd <- which(SummarizedExperiment::colData(sce)
            [[sampleColname]] == sample)
        scores <- scran::doubletCells(sce[, sceSampleInd],
            ...,
            assay.type = assayType)
        doubletScore[sceSampleInd] <- scores
    }
    SummarizedExperiment::colData(sce)$scran_doublet_score <- doubletScore
    return(sce)
}


#' @title Detect doublet cells using \link[scran]{doubletCells}.
#' @description A wrapper function for \link[scran]{doubletCells}. Identify
#'  potential doublet cells based on simulations of putative doublet expression
#'  profiles. Generate a doublet score for each cell.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object. Must
#'  contain a count matrix.
#' @param sampleColname Character. The column name which specifies the sample
#'  origin in the \link[SummarizedExperiment]{colData} of the provided
#'  \code{sce} object. Default "sample".
#' @param ... Additional arguments to pass to \link[scran]{doubletCells}.
#' @param assayType  A string specifying which assay values contain the count
#'  matrix.
#' @details This function is a wrapper function for \link[scran]{doubletCells}.
#'  It assumes the provided \code{sce} was generated using the import
#'  functions in this package so there
#'  will be a "sample" column in the \link[SummarizedExperiment]{colData}
#'  specifying the sample origin of each cell. The column name can be specified
#'  using the \code{sampleColname} argument.
#'  \link{runDoubletCells} runs \link[scran]{doubletCells} for each
#'  sample within \code{sce} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link[SummarizedExperiment]{colData} of \code{sce}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  'scran_doublet_score' column added to the
#'  \link[SummarizedExperiment]{colData} slot.
#' @references Lun ATL (2018). Detecting doublet cells with scran.
#'  \url{https://ltla.github.io/SingleCellThoughts/software/
#' doublet_detection/bycell.html}
#' @seealso \link[scran]{doubletCells}
#' @examples
#' data(emptyDropsSceExample, package = "scruff")
#' sce <- runDoubletCells(emptyDropsSceExample)
#' @export
#' @import scran
runDoubletCells <- function(sce,
    sampleColname = "sample",
    ...,
    assayType = "counts") {

    sce <- .runDoubletCells(sce = sce,
        sampleColname = sampleColname,
        ...,
        assayType = assayType)
    return(sce)
}
