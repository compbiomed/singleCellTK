#' @title Generates a SingleCellExperiment object containing the output of
#'  specified QC functions
#' @description A wrapper function for the individual QC algorithms.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object. Must
#'  contain a count matrix.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are: \link[DropletUtils]{emptyDrops},
#'  \link[scran]{doubletCells}.
#' @param sampleColname Character. The column name which specifies the sample
#'  origin in the \link[SummarizedExperiment]{colData} of the provided
#'  \code{sce} object. Default "sample".
#' @param assayType  A string specifying which assay values contain the count
#'  matrix.
#' @param ... Additional arguments to pass to \code{algorithms}.
#' @details
#'  \code{runQC} by default runs all available QC algorithms
#'  (\link[DropletUtils]{emptyDrops}, \link[scran]{doubletCells}).
#'  \code{runQCFilteredCells} runs QC
#'  algorithm (\link[scran]{doubletCells}) on \code{sce} containing
#'  filtered cells.
#'  \code{runQCAllDroplets}, runs QC
#'  algorithm (\link[DropletUtils]{emptyDrops}) on \code{sce} containing
#'  unfiltered cells.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms.
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
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


#' @rdname runQC
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runQCFilteredCells(emptyDropsSceExample)
#' @export
runQCFilteredCells <- function(sce,
    algorithms = c("doubletCells"),
    sampleColname = "sample",
    assayType = "counts",
    ...) {

    if ("doubletCells" %in% algorithms) {
        sce <- runDoubletCells(sce = sce,
            sampleColname = sampleColname,
            ...,
            assayType = assayType)
    }
    return(sce)
}


#' @rdname runQC
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runQCAllDroplets(emptyDropsSceExample)
#' @export
runQCAllDroplets <- function(sce,
    algorithms = c("emptyDrops"),
    sampleColname = "sample",
    assayType = "counts",
    ...) {

    if ("emptyDrops" %in% algorithms) {
        sce <- runEmptyDrops(sce = sce,
            sampleColname = sampleColname,
            ...,
            assayType = assayType)
    }

    return(sce)
}



