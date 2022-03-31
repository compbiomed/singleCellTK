#' @title Perform comprehensive single cell QC
#' @description A wrapper function to run several QC algorithms on a
#'  SingleCellExperiment
#'  object containing cells after empty droplets have been removed.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param algorithms Character vector. Specify which QC algorithms to run.
#'  Available options are "QCMetrics", "scrublet", "doubletFinder", "scDblFinder", 
#'  "cxds", "bcds", "cxds_bcds_hybrid", "decontX" and "soupX".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Algorithms will be run on cells from each sample separately.
#' @param collectionName Character. Name of a \code{GeneSetCollection} obtained by 
#'  using one of the importGeneSet* functions. Default \code{NULL}.
#' @param geneSetList See \code{runPerCellQC}. Default NULL.
#' @param geneSetListLocation See \code{runPerCellQC}. Default NULL.
#' @param geneSetCollection See \code{runPerCellQC}. Default NULL.
#' @param mitoRef,mitoIDType,mitoPrefix,mitoID,mitoGeneLocation Arguments used to 
#'  import mitochondrial genes and quantify their expression. Please see 
#'  \link[singleCellTK]{runPerCellQC} for detailed information.  
#' @param useAssay  A string specifying which assay contains the count
#'  matrix for cells.
#' @param background A \link[SingleCellExperiment]{SingleCellExperiment}
#' with the matrix located in the assay slot under \code{bgAssayName}. It should have 
#' the same structure as inSCE except it contains the matrix of empty droplets instead 
#' of cells. When supplied, empirical distribution of transcripts from these 
#' empty droplets will be used as the contamination distribution. It is only used in 
#' algorithms "decontX" and "soupX". Default NULL.
#' @param bgAssayName Character. Name of the assay to use if background is a 
#' \link[SingleCellExperiment]{SingleCellExperiment}. If NULL, the function
#' will use the same value as \code{useAssay}. It is only used in algorithms 
#' "decontX" and "soupX". Default is NULL. 
#' @param bgBatch Batch labels for \code{background}. If \code{background} is a 
#' \link[SingleCellExperiment]{SingleCellExperiment} object, this can be a single 
#' character specifying a name that can be found in \code{colData(background)} 
#' to directly use the barcode annotation Its unique values should be the same
#' as those in \code{sample}, such that each batch of cells have their corresponding 
#' batch of empty droplets as background, pointed by this parameter. It is only used in
#' algorithms "decontX" and "soupX". Default to NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param paramsList A list containing parameters for QC functions. Default NULL.
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link{colData}
#' of \code{inSCE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' sce <- runCellQC(sce)
#' }
#' @export

runCellQC <- function(inSCE,
  algorithms = c("QCMetrics", "scDblFinder", "cxds", "bcds",
    "cxds_bcds_hybrid", "decontX", "decontX_bg", "soupX", "soupX_bg"), #"scrublet", "doubletFinder",
  sample = NULL,
  collectionName = NULL,
  geneSetList = NULL,
  geneSetListLocation = "rownames",
  geneSetCollection = NULL,
  mitoRef = NULL,
  mitoIDType = NULL,
  mitoPrefix = NULL,
  mitoID = NULL,
  mitoGeneLocation = NULL,
  useAssay = "counts",
  background = NULL,
  bgAssayName= NULL,
  bgBatch = NULL,
  seed = 12345,
  paramsList = NULL) {

  nonmatch <- setdiff(algorithms, c("scDblFinder", "cxds", "bcds",
    "cxds_bcds_hybrid", "decontX", "decontX_bg", "QCMetrics", "scrublet", 
    "doubletFinder", "soupX", "soupX_bg"))
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
        mitoRef = mitoRef,
        mitoIDType = mitoIDType,
        mitoPrefix = mitoPrefix,
        mitoID = mitoID,
        mitoGeneLocation = mitoGeneLocation,
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

  if ("scDblFinder" %in% algorithms) {
    inSCE <- do.call(runScDblFinder,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay,
      seed = seed),
      paramsList[["scDblFinder"]]))
  }

  if ("doubletFinder" %in% algorithms) {
    inSCE <- do.call(runDoubletFinder,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay,
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
    if (!is.null(paramsList[["decontX"]][['background']])) {
      warning("'decontX' algorithm will decontaminate ambient RNA without the background count ",
        "matrix. Ignore 'background' parameter within 'paramsList[['decontX']]' will be ignored. If you want ",
        "to adjust decontamination with background matrix, please run 'decontX_bg' algorithm.")
      paramsList[["decontX"]][['background']] <- NULL
    }
    inSCE <- do.call(runDecontX,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay,
      seed = seed),
      paramsList[["decontX"]]))
  }

  if ("decontX_bg" %in% algorithms) {
    if (!is.null(background)) {
      inSCE <- do.call(runDecontX,
        c(list(inSCE = quote(inSCE),
        sample = sample,
        useAssay = useAssay,
        seed = seed, 
        background = background,
        bgAssayName = bgAssayName,
        bgBatch = bgBatch),
        paramsList[["decontX_bg"]]))      
    } else {
      warning("'background' is NULL. Skip 'decontX_bg' algorithm.")
    }
  }

  if ("soupX" %in% algorithms) {
    if (!is.null(paramsList[["soupX"]][['background']])) {
      warning("'soupX' algorithm step will decontaminate ambient RNA without the background count ",
        "matrix. 'background' parameter within 'paramsList[['soupX']]' will be ignored. If you want ",
        "to adjust decontamination with background matrix, please run 'soupX_bg' algorithm.")
      paramsList[["soupX"]][['background']] <- NULL
    }
    inSCE <- do.call(runSoupX,
      c(list(inSCE = quote(inSCE),
      sample = sample,
      useAssay = useAssay),
      paramsList[["soupX"]]))
  }

  if ("soupX_bg" %in% algorithms) {
    if (!is.null(background)) {
      inSCE <- do.call(runSoupX,
        c(list(inSCE = quote(inSCE),
        sample = sample,
        useAssay = useAssay,
        background = background,
        bgAssayName = bgAssayName,
        bgBatch = bgBatch),
        paramsList[["soupX_bg"]]))      
    } else {
      warning("'background' is NULL. Skip 'decontX_bg' algorithm. ")
    }
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
#'  specified algorithms in the \link{colData}
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



