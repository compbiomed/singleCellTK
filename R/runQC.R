#' @title Perform comprehensive single cell QC
#' @description A wrapper function to run several QC algorithms on a
#'  SingleCellExperiment
#'  object containing cells after empty droplets have been removed.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "QCMetrics", "scrublet", "doubletCells", "cxds", "bcds", "cxds_bcds_hybrid", and "decontX".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param geneSetList See \code{runPerCellQC}. Default NULL.
#' @param geneSetListLocation See \code{runPerCellQC}. Default NULL.
#' @param geneSetCollection See \code{runPerCellQC}. Default NULL.
#' @param assayName  A string specifying which assay contains the count
#'  matrix for cells.
#' @param seed Seed for the random number generator. Default 12345.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{sce}.
#' @examples
#' \dontrun{
#' data(sce_chcl, package = "scds")
#' sce <- runCellQC(sce_chcl)
#' }
#' @export
runCellQC <- function(sce,
  algorithms = c("QCMetrics", "doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "scrublet", "decontX"),
  sample = NULL,
  geneSetList = NULL,
  geneSetListLocation = "rownames",
  geneSetCollection = NULL,
  assayName = "counts",
  seed = 12345) {

  nonmatch <- setdiff(algorithms, c("doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder"))
  if (length(nonmatch) > 0) {
    stop("'", paste(nonmatch, collapse=","), "' are not supported algorithms.")
  }

  if ("QCMetrics" %in% algorithms) {
    sce <- runPerCellQC(sce = sce, assayName = assayName,
                                 geneSetList = geneSetList,
                                 geneSetListLocation = geneSetListLocation,
                                 geneSetCollection = geneSetCollection)
  }

  if ("scrublet" %in% algorithms) {
    sce <- runScrublet(sce = sce,
      sample = sample,
      assayName = assayName,
      seed = seed)
  }

  if ("doubletCells" %in% algorithms) {
    sce <- runDoubletCells(sce = sce,
      sample = sample,
      assayName = assayName,
      seed = seed)
  }
  
  if ("doubletFinder" %in% algorithms) {
    sce <- runDoubletFinder(sce = sce,
      sample = sample,
      seed = seed)
  }

  if ("cxds" %in% algorithms) {
    sce <- runCxds(sce = sce,
      sample = sample,
      seed = seed)
  }

  if ("bcds" %in% algorithms) {
    sce <- runBcds(sce = sce,
      sample = sample,
      seed = seed)
  }

  if ("cxds_bcds_hybrid" %in% algorithms) {
    sce <- runCxdsBcdsHybrid(sce = sce,
      sample = sample,
      seed = seed)
  }

  if ("decontX" %in% algorithms) {
    sce <- runDecontX(sce = sce,
      sample = sample,
      assayName = assayName,
      seed = seed)
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
  algorithms = c("QCMetrics", "emptyDrops", "barcodeRanks"),
  sample = NULL,
  assayName = "counts") {

  nonmatch <- setdiff(algorithms, c("QCMetrics", "emptyDrops", "barcodeRanks"))
  if(length(nonmatch) > 0) {
    stop(paste0("'", paste(nonmatch, collapse=","), "' are not supported algorithms."))
  }

  ## emptyDrops and barcodeRanks need dgCMatrix objects as input
  ## Convert once for both functions
  counts.class <- class(SummarizedExperiment::assay(sce, i = assayName))
  if (class(counts.class) != "dgCMatrix") {
    SummarizedExperiment::assay(sce, i = assayName) <-
      methods::as(SummarizedExperiment::assay(sce, i = assayName), "dgCMatrix")
  }

  if ("QCMetrics" %in% algorithms) {
    sce <- runPerCellQC(sce = sce, assayName = assayName)
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
    methods::as(SummarizedExperiment::assay(sce, i = assayName), counts.class)

  return(sce)
}



