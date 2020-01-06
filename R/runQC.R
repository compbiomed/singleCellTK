#' @title Perform comprehensive single cell QC
#' @description A wrapper function to run several QC algorithms on a SingleCellExperiment
#' object containing cells after empty droplets have been removed.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "decontX" and "doubletCells".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param assayName  A string specifying which assay contains the count
#'  matrix for cells.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{sce}.
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runCellQC(emptyDropsSceExample,
#'   sample = colData(emptyDropsSceExample)$sample)
#' @export
runCellQC <- function(sce,
  algorithms = c("doubletCells", "DecontX"),
  sample = NULL,
  assayName = "counts",
  ...) {

  nonmatch <- setdiff(algorithms, c("doubletCells", "decontX"))
  if(length(nonmatch) > 0) {
    stop(paste0("'", paste(nonmatch, collapse=","), "' are not supported algorithms."))
  }

  if ("doubletCells" %in% algorithms) {
    sce <- runDoubletCells(sce = sce,
      sample = sample,
      ...,
      assayName = assayName)
  }

  if ("decontX" %in% algorithms) {
    sce <- celda::decontX(sce = sce,
      batch = sample,
      ...,
      assayName = assayName)
  }

  return(sce)
}


#' @title Perform comprehensive droplet QC
#' @description A wrapper function to run several QC algorithms for determining
#' empty droplets in single cell RNA-seq data
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the full droplet count matrix
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "emptyDrops" and "barcodeRanks".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param assayName  A string specifying which assay contains the count
#'  matrix for droplets.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{sce}.
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runDropletQC(emptyDropsSceExample,
#'   sample = colData(emptyDropsSceExample)$sample)
#' @export
runDropletQC <- function(sce,
  algorithms = c("emptyDrops", "barcodeRanks"),
  sample = NULL,
  assayName = "counts") {

  nonmatch <- setdiff(algorithms, c("emptyDrops", "barcodeRanks"))
  if(length(nonmatch) > 0) {
    stop(paste0("'", paste(nonmatch, collapse=","), "' are not supported algorithms."))
  }

  ## emptyDrops and barcodeRanks need dgCMatrix objects as input
  ## Convert once for both functions
  counts.class <- class(SummarizedExperiment::assay(sce, i = assayName))
  if (class(counts.class) != "dgCMatrix") {
    SummarizedExperiment::assay(sce, i = assayName) <-
      as(SummarizedExperiment::assay(sce, i = assayName), "dgCMatrix")
  }

  if (any("emptyDrops" %in% algorithms)) {
    sce <- runEmptyDrops(sce = sce,
      sample = sample,
      assayName = assayName)
  }

  if (any("barcodeRanks" %in% algorithms)) {
    sce <- runBarcodeRankDrops(sce = sce,
      sample = sample,
      assayName = assayName)
  }

  ## Convert back to original class
  SummarizedExperiment::assay(sce, i = assayName) <-
    as(SummarizedExperiment::assay(sce, i = assayName), counts.class)

  return(sce)
}



