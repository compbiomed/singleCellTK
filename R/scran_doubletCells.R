# .runDoubletCells <- function(cell.matrix = cell.matrix,
#                               k = k,
#                               nIters = nIters,
#                               size.factors.norm = NULL,
#                               size.factors.content = NULL,
#                               subset.row = NULL,
#                               block = 10000,
#                               d = 50,
#                               force.match=FALSE,
#                               force.k=20,
#                               force.ndist=3,
#                               BNPARAM=BNPARAM,
#                               BSPARAM=BSPARAM,
#                               BPPARAM=BPPARAM
#                               ) {
# 
#   cell.matrix <- .convertToMatrix(cell.matrix)
# 
#   scores <- matrix(scran::doubletCells(cell.matrix, k = k,
#                                        niters = nIters,
#                                        size.factors.norm = NULL,
#                                        size.factors.content = NULL,
#                                        subset.row = NULL,
#                                        block = 10000,
#                                        d = 50,
#                                        force.match=FALSE,
#                                        force.k=20,
#                                        force.ndist=3,
#                                        BNPARAM=BNPARAM,
#                                        BSPARAM=BSPARAM,
#                                        BPPARAM=BPPARAM
#                                        ), ncol=1)
#   scores <- cbind(scores,log10(scores[,1]+1))
#   colnames(scores) <- c("scran_doubletCells_score", "scran_doubletCells_score_log10")
# 
# 
#   return(scores)
# }
# 

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
#'  \code{runDoubletCells} runs \link[scDblFinder]{scDblFinder} for each
#'  \code{sample} within \code{inSCE} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link{colData} of \code{inSCE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  'scran_doubletCells_score' column added to the
#'  \link{colData} slot.
#' @references Lun ATL (2018). Detecting doublet cells with scran.
#'  \url{https://ltla.github.io/SingleCellThoughts/software/
#' doublet_detection/bycell.html}
#' @seealso \link[scran]{doubletCells}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDoubletCells(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
runDoubletCells <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    nNeighbors = 50,
    simDoublets = max(10000, ncol(inSCE)),
    seed = 12345,
    # size.factors.norm = NULL,
    # size.factors.content = NULL,
    # subset.row = NULL,
    # block = 10000,
    # d = 50,
    # force.match=FALSE,
    # force.k=20,
    # force.ndist=3,
    # BNPARAM=BiocNeighbors::KmknnParam(),
    # BSPARAM=BiocSingular::bsparam(),
    BPPARAM=BiocParallel::SerialParam()
) {
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))

  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }

  message(paste0(date(), " ... Running 'doubletCells'"))

  ## Define result matrix for all samples
  # output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
  #           scran_doubletCells_score = numeric(ncol(inSCE)),
  #           scran_doubletCells_score_log10 = numeric(ncol(inSCE)))

  ## Loop through each sample and run barcodeRank
  #samples <- unique(sample)
  
  inSCE <- withr::with_seed(seed,
                            scDblFinder::scDblFinder(sce = inSCE,
                            samples = sample,
                            artificialDoublets = simDoublets,
                            k = nNeighbors,
                            verbose = FALSE
                            ))
  names(SummarizedExperiment::colData(inSCE)) <- gsub(pattern = "scDblFinder\\.",
                                                      "scran_doubletCells_",
                                                      names(SummarizedExperiment::colData(inSCE)))
  
  # for (i in seq_len(length(samples))) {
  #   sceSampleInd <- sample == samples[i]
  #   sceSample <- inSCE[, sceSampleInd]
  # 
  #   mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
  # 
  #   result <- withr::with_seed(seed,
  #             .runDoubletCells(cell.matrix = mat,
  #                              k = nNeighbors,
  #                              nIters = simDoublets,
  #                              size.factors.norm = NULL,
  #                              size.factors.content = NULL,
  #                              subset.row = NULL,
  #                              block = 10000,
  #                              d = 50,
  #                              force.match=FALSE,
  #                              force.k=20,
  #                              force.ndist=3,
  #                              BNPARAM=BNPARAM,
  #                              BSPARAM=BSPARAM,
  #                              BPPARAM=BPPARAM
  #                              ))
  # 
  #   output[sceSampleInd, ] <- result
  # }

  argsList <- argsList[!names(argsList) %in% c("BPPARAM")]

  inSCE@metadata$runDoubletCells <- argsList[-1]
  inSCE@metadata$runDoubletCells$packageVersion <- utils::packageDescription("scDblFinder")$Version

  return(inSCE)
}

