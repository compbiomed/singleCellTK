#' @title Perform comprehensive single cell QC
#' @description A wrapper function to run several QC algorithms on a
#'  SingleCellExperiment
#'  object containing cells after empty droplets have been removed.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "QCMetrics", "scrublet", "doubletCells", "cxds", "bcds", "cxds_bcds_hybrid", and "decontX".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param collectionName Character. Name of a \code{GeneSetCollection} obtained by using one of the importGeneSet* functions. Default \code{NULL}.
#' @param geneSetList See \code{runPerCellQC}. Default NULL.
#' @param geneSetListLocation See \code{runPerCellQC}. Default NULL.
#' @param geneSetCollection See \code{runPerCellQC}. Default NULL.
#' @param useAssay  A string specifying which assay contains the count
#'  matrix for cells.
#' @param seed Seed for the random number generator. Default 12345.
#' @param paramsList A list containing parameters for QC functions. Default NULL.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{inSCE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' sce <- runCellQC(sce)
#' }
#' @export
runCellQC <- function(inSCE,
  algorithms = c("QCMetrics", "doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "scrublet", "doubletFinder", "decontX"),
  sample = NULL,
  collectionName = NULL,
  geneSetList = NULL,
  geneSetListLocation = "rownames",
  geneSetCollection = NULL,
  useAssay = "counts",
  seed = 12345,
  paramsList = NULL) {

  nonmatch <- setdiff(algorithms, c("doubletCells", "cxds", "bcds",
    "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder"))
  if (length(nonmatch) > 0) {
    stop("'", paste(nonmatch, collapse=","), "' are not supported algorithms.")
  }

  if ("QCMetrics" %in% algorithms) {
    inSCE <- do.call(runPerCellQC,
      c(list(inSCE = quote(inSCE),
        useAssay = useAssay,
        collectionName = collectionName,
        geneSetList = geneSetList,
        geneSetListLocation = geneSetListLocation,
        geneSetCollection = geneSetCollection),
        paramsList[["QCMetrics"]]))
  }

  if ("scrublet" %in% algorithms) {

    inSCE <- do.call(runScrublet,
      c(list(inSCE = quote(inSCE),
        sample = sample,
        useAssay = useAssay,
        seed = seed),
        paramsList[["scrublet"]]))
  }

  if ("doubletCells" %in% algorithms) {
    inSCE <- do.call(runDoubletCells,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay,
      seed = seed),
      paramsList[["doubletCells"]]))
  }

  if ("doubletFinder" %in% algorithms) {
    inSCE <- do.call(runDoubletFinder,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      seed = seed),
      paramsList[["doubletFinder"]]))
  }

  if ("cxds" %in% algorithms) {
    inSCE <- do.call(runCxds,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      seed = seed,
      useAssay = useAssay,
      estNdbl = TRUE),
      paramsList[["cxds"]]))
  }

  if ("bcds" %in% algorithms) {
    inSCE <- do.call(runBcds,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      seed = seed,
      useAssay = useAssay,
      estNdbl = TRUE),
      paramsList[["bcds"]]))
  }

  if ("cxds_bcds_hybrid" %in% algorithms) {
    inSCE <- do.call(runCxdsBcdsHybrid,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      seed = seed,
      useAssay = useAssay,
      estNdbl = TRUE),
      paramsList[["cxds_bcds_hybrid"]]))
  }

  if ("decontX" %in% algorithms) {
    inSCE <- do.call(runDecontX,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay,
      seed = seed),
      paramsList[["decontX"]]))
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
#' @param paramsList A list containing parameters for QC functions. Default NULL.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link[SummarizedExperiment]{colData}
#' of \code{inSCE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runDropletQC(sce)
#' }
#' @export
runDropletQC <- function(inSCE,
  algorithms = c("QCMetrics", "emptyDrops", "barcodeRanks"),
  sample = NULL,
  useAssay = "counts",
  paramsList = NULL) {

  nonmatch <- setdiff(algorithms, c("QCMetrics", "emptyDrops", "barcodeRanks"))
  if(length(nonmatch) > 0) {
    stop(paste0("'", paste(nonmatch, collapse=","), "' are not supported algorithms."))
  }

  if ("QCMetrics" %in% algorithms) {
    inSCE <- do.call(runPerCellQC,
      c(list(inSCE = quote(inSCE),
        useAssay = useAssay),
        paramsList[["QCMetrics"]]))
  }

  if (any("emptyDrops" %in% algorithms)) {
    inSCE <- do.call(runEmptyDrops,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay),
      paramsList[["emptyDrops"]]))
  }

  if (any("barcodeRanks" %in% algorithms)) {
    inSCE <- do.call(runBarcodeRankDrops,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay),
      paramsList[["barcodeRanks"]]))
  }

  return(inSCE)
}



