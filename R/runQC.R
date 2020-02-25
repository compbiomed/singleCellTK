#' @title Perform comprehensive single cell QC
#' @description A wrapper function to run several QC algorithms on a
#'  SingleCellExperiment
#'  object containing cells after empty droplets have been removed.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "QCMetrics", "scrublet", "doubletCells", "cxds", "bcds", "cxds_bcds_hybrid", and "decontX".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param geneSetList See \code{runPerCellQC}. Default NULL.
#' @param geneSetListLocation See \code{runPerCellQC}. Default NULL.
#' @param geneSetCollection See \code{runPerCellQC}. Default NULL.
#' @param useAssay  A string specifying which assay contains the count
#'  matrix for cells.
#' @param seed Seed for the random number generator. Default 12345.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{inSCE}.
#' @examples
#' \dontrun{
#' data(sce_chcl, package = "scds")
#' sce <- runCellQC(sce_chcl)
#' }
#' @export
runCellQC <- function(inSCE,
  algorithms = c("QCMetrics", "doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "scrublet", "doubletFinder", "decontX"),
  sample = NULL,
  geneSetList = NULL,
  geneSetListLocation = "rownames",
  geneSetCollection = NULL,
  useAssay = "counts",
  seed = 12345) {

  nonmatch <- setdiff(algorithms, c("doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder"))
  if (length(nonmatch) > 0) {
    stop("'", paste(nonmatch, collapse=","), "' are not supported algorithms.")
  }

  if ("QCMetrics" %in% algorithms) {
    inSCE <- runPerCellQC(inSCE = inSCE, useAssay = useAssay,
                                 geneSetList = geneSetList,
                                 geneSetListLocation = geneSetListLocation,
                                 geneSetCollection = geneSetCollection)
  }

  if ("scrublet" %in% algorithms) {
    inSCE <- runScrublet(inSCE = inSCE,
      sample = sample,
      useAssay = useAssay,
      seed = seed)
  }

  if ("doubletCells" %in% algorithms) {
    inSCE <- runDoubletCells(inSCE = inSCE,
      sample = sample,
      useAssay = useAssay,
      seed = seed)
  }

  if ("doubletFinder" %in% algorithms) {
    inSCE <- runDoubletFinder(inSCE = inSCE,
      sample = sample,
      seed = seed)
  }

  if ("cxds" %in% algorithms) {
    inSCE <- runCxds(inSCE = inSCE,
      sample = sample,
      seed = seed)
  }

  if ("bcds" %in% algorithms) {
    inSCE <- runBcds(inSCE = inSCE,
      sample = sample,
      seed = seed)
  }

  if ("cxds_bcds_hybrid" %in% algorithms) {
    inSCE <- runCxdsBcdsHybrid(inSCE = inSCE,
      sample = sample,
      seed = seed)
  }

  if ("decontX" %in% algorithms) {
    inSCE <- runDecontX(inSCE = inSCE,
      sample = sample,
      useAssay = useAssay,
      seed = seed)
  }

  return(inSCE)
}


#' @title Perform comprehensive droplet QC
#' @description A wrapper function to run several QC algorithms for determining
#' empty droplets in single cell RNA-seq data
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the full droplet count matrix
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "emptyDrops" and "barcodeRanks".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param useAssay  A string specifying which assay contains the count
#'  matrix for droplets.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{inSCE}.
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runDropletQC(emptyDropsSceExample,
#'   sample = colData(emptyDropsSceExample)$sample)
#' @export
runDropletQC <- function(inSCE,
  algorithms = c("QCMetrics", "emptyDrops", "barcodeRanks"),
  sample = NULL,
  useAssay = "counts") {

  nonmatch <- setdiff(algorithms, c("QCMetrics", "emptyDrops", "barcodeRanks"))
  if(length(nonmatch) > 0) {
    stop(paste0("'", paste(nonmatch, collapse=","), "' are not supported algorithms."))
  }

  if ("QCMetrics" %in% algorithms) {
    inSCE <- runPerCellQC(inSCE = inSCE, useAssay = useAssay)
  }

  if (any("emptyDrops" %in% algorithms)) {
    inSCE <- runEmptyDrops(inSCE = inSCE,
      sample = sample,
      useAssay = useAssay)
  }

  if (any("barcodeRanks" %in% algorithms)) {
    inSCE <- runBarcodeRankDrops(inSCE = inSCE,
      sample = sample,
      useAssay = useAssay)
  }

  return(inSCE)
}



