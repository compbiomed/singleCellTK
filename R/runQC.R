#' @title Generates a SingleCellExperiment object containing the output of
#'  specified QC functions
#' @description A wrapper function for the individual QC algorithms.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object. Must
#'  contain a count matrix.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Possible options are: "emptyDrops", "doubletCells".
#' @param sampleColname Character. The column name which specifies the sample
#'  origin in the \link[SummarizedExperiment]{colData} of the provided
#'  \code{sce} object. Default "sample".
#' @param assayType  A string specifying which assay values contain the count
#'  matrix.
#' @param ... Additional arguments to pass to \code{algorithms}.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms.
#' @examples
#' data(emptyDropsSceExample, package = "scruff")
#' sce <- runQC(emptyDropsSceExample,
#'     algorithms = c("emptyDrops", "doubletCells"))
#' @export
runQC <- function(sce,
    algorithms = c("emptyDrops", "doubletCells"),
    sampleColname = "sample",
    assayType = "counts",
    ...) {

    if ("emptyDrops" %in% algorithms) {
        sce <- runEmptyDrops(sce = sce,
            sampleColname = sampleColname,
            ...,
            assayType = assayType)
    }

    if ("doubletCells" %in% algorithms) {
        sce <- runDoubletCells(sce = sce,
            sampleColname = sampleColname,
            ...,
            assayType = assayType)
    }
    return(sce)
}
