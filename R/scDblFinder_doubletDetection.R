#' @title Detect doublet cells using \link[scDblFinder]{scDblFinder}.
#' @description A wrapper function for \link[scDblFinder]{scDblFinder}. Identify
#'  potential doublet cells based on simulations of putative doublet expression
#'  profiles. Generate a doublet score for each cell.
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param sample Character vector or colData variable name. Indicates which
#' sample each cell belongs to. Default \code{NULL}.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#' \code{"counts"}.
#' @param nNeighbors Number of nearest neighbors used to calculate density for
#' doublet detection. Default \code{50}.
#' @param simDoublets Number of simulated doublets created for doublet
#' detection. Default \code{10000}.
#' @param seed Seed for the random number generator, can be set to \code{NULL}.
#' Default \code{12345}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam-class}} object
#' specifying whether the neighbour searches should be parallelized. Default
#' \code{BiocParallel::SerialParam(RNGseed = seed)}.
#' @details This function is a wrapper function for
#' \link[scDblFinder]{scDblFinder}. \code{runScDblFinder} runs
#' \link[scDblFinder]{scDblFinder} for each sample within \code{inSCE}
#' iteratively. The resulting doublet scores for all cells will be appended to
#' the \code{\link{colData}} of \code{inSCE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#' scDblFinder QC outputs added to the \link{colData} slot.
#' @references Lun ATL (2018). Detecting doublet cells with scran.
#'  \url{https://ltla.github.io/SingleCellThoughts/software/doublet_detection/bycell.html}
#' @seealso \code{\link[scDblFinder]{scDblFinder}},
#' \code{\link{plotScDblFinderResults}}, \code{\link{runCellQC}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runScDblFinder(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<- assayNames assayNames<-
#' @importFrom S4Vectors metadata<- metadata
runScDblFinder <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    nNeighbors = 50,
    simDoublets = max(10000, ncol(inSCE)),
    seed = 12345,
    BPPARAM = BiocParallel::SerialParam(RNGseed = seed)
) {
  tempSCE <- inSCE
  #assayNames(inSCE)[which(useAssay %in% assayNames(inSCE))] <- "counts"
  #useAssay <- "counts"

  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  argsList <- argsList[!names(argsList) %in% c("inSCE", "BPPARAM")]
  argsList$packageVersion <- utils::packageDescription("scDblFinder")$Version

  sample <- .manageCellVar(inSCE, var = sample)
  if (is.null(sample)) {
    sample = rep(1, ncol(inSCE))
  }


  message(date(), " ... Running 'scDblFinder'")

  ## Loop through each sample and run barcodeRank

  rm.ix <- which(colSums(assay(inSCE, useAssay)) == 0)
  if (length(rm.ix) > 0) {
    inSCEOrig <- inSCE
    inSCE <- inSCE[,-rm.ix]
    sample <- sample[-rm.ix]
  }
  withr::with_seed(seed, {
    inSCE <- scDblFinder::scDblFinder(sce = inSCE,
                                      samples = sample,
                                      artificialDoublets = simDoublets,
                                      k = nNeighbors,
                                      verbose = FALSE,
                                      BPPARAM = BPPARAM
    )
  })
  if (length(rm.ix) > 0) {
    inSCE <- mergeSCEColData(inSCE1 = inSCEOrig, inSCE2 = inSCE)
  }

  names(colData(inSCE)) <- gsub(pattern = "scDblFinder\\.",
                                "scDblFinder_",
                                names(colData(inSCE)))

  names(colData(inSCE)) <- gsub(pattern = "scDblFinder_score",
                                "scDblFinder_doublet_score",
                                names(colData(inSCE)))
  names(colData(inSCE)) <- gsub(pattern = "scDblFinder_class",
                                "scDblFinder_doublet_call",
                                names(colData(inSCE)))

  levels(inSCE$scDblFinder_doublet_call) <- list(Singlet = "singlet",
                                                 Doublet = "doublet")

  if (all(sample == 1)) {
    metadata(inSCE)$sctk$runScDblFinder$all_cells <- argsList
  } else {
    metadata(inSCE)$sctk$runScDblFinder <- sapply(unique(sample),
                                                  function(x) return(argsList),
                                                  simplify = FALSE,
                                                  USE.NAMES = TRUE)
  }

  colData(tempSCE) <- colData(inSCE)
  metadata(tempSCE) <- metadata(inSCE)

  return(tempSCE)
}
