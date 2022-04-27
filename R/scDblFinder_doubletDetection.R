#' @title Detect doublet cells using \link[scDblFinder]{scDblFinder}.
#' @description A wrapper function for \link[scDblFinder]{scDblFinder}. Identify
#'  potential doublet cells based on simulations of putative doublet expression
#'  profiles. Generate a doublet score for each cell.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[scDblFinder]{scDblFinder} will be run on cells from each sample separately.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @param nNeighbors Number of nearest neighbors used to calculate density for
#'  doublet detection. Default 50.
#' @param simDoublets Number of simulated doublets created for doublet
#'  detection. Default 10000.
#' @param seed Seed for the random number generator. Default 12345.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether the
#'  neighbour searches should be parallelized.
#' @details This function is a wrapper function for \link[scDblFinder]{scDblFinder}.
#'  \code{runScDblFinder} runs \link[scDblFinder]{scDblFinder} for each
#'  \code{sample} within \code{inSCE} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link{colData} of \code{inSCE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  scDblFinder QC outputs added to the
#'  \link{colData} slot.
#' @references Lun ATL (2018). Detecting doublet cells with scran.
#'  \url{https://ltla.github.io/SingleCellThoughts/software/doublet_detection/bycell.html}
#' @seealso \link[scDblFinder]{scDblFinder}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runScDblFinder(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
runScDblFinder <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    nNeighbors = 50,
    simDoublets = max(10000, ncol(inSCE)),
    seed = 12345,
    BPPARAM=BiocParallel::SerialParam()
) {

  tempSCE <- inSCE
  SummarizedExperiment::assayNames(inSCE)[which(useAssay %in% SummarizedExperiment::assayNames(inSCE))] <- "counts"
  useAssay <- "counts"

  argsList <- mget(names(formals()),sys.frame(sys.nframe()))

  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }

  message(paste0(date(), " ... Running 'scDblFinder'"))

  ## Loop through each sample and run barcodeRank

  rm.ix <- which(colSums(assay(inSCE, useAssay)) == 0)
  if(length(rm.ix) > 0){
    inSCEOrig <- inSCE
    inSCE <- inSCE[,-rm.ix]
    sample <- sample[-rm.ix]
  }
  inSCE <- withr::with_seed(seed,
                            scDblFinder::scDblFinder(sce = inSCE,
                            samples = sample,
                            artificialDoublets = simDoublets,
                            k = nNeighbors,
                            verbose = FALSE,
                            BPPARAM = BPPARAM
                            ))
  if(length(rm.ix) > 0){
    inSCE <- mergeSCEColData(inSCE1 = inSCEOrig, inSCE2 = inSCE)
  }

  names(SummarizedExperiment::colData(inSCE)) <- gsub(pattern = "scDblFinder\\.",
                                                      "scDblFinder_",
                                                      names(SummarizedExperiment::colData(inSCE)))

  names(SummarizedExperiment::colData(inSCE)) <- gsub(pattern = "scDblFinder_score",
                                                      "scDblFinder_doublet_score",
                                                      names(SummarizedExperiment::colData(inSCE)))
  names(SummarizedExperiment::colData(inSCE)) <- gsub(pattern = "scDblFinder_class",
                                                      "scDblFinder_doublet_call",
                                                      names(SummarizedExperiment::colData(inSCE)))

  levels(inSCE$scDblFinder_doublet_call) <- list(Singlet = "singlet", Doublet = "doublet")

  argsList <- argsList[!names(argsList) %in% c("BPPARAM")]

  inSCE@metadata$runScDblFinder <- argsList[-1]
  inSCE@metadata$runScDblFinder$packageVersion <- utils::packageDescription("scDblFinder")$Version

  tempSCE@colData <- inSCE@colData
  tempSCE@metadata <- inSCE@metadata

  return(tempSCE)
}
