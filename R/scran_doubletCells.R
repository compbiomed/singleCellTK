.runDoubletCells <- function(cell.matrix = cell.matrix,
                              k = k,
                              nIters = nIters,
                              size.factors.norm = NULL,
                              size.factors.content = NULL,
                              subset.row = NULL,
                              block = 10000,
                              d = 50,
                              force.match=FALSE,
                              force.k=20,
                              force.ndist=3,
                              BNPARAM=BNPARAM,
                              BSPARAM=BSPARAM,
                              BPPARAM=BPPARAM
                              ) {

  cell.matrix <- .convertToMatrix(cell.matrix)

  scores <- matrix(scran::doubletCells(cell.matrix, k = k,
                                       niters = nIters,
                                       size.factors.norm = NULL,
                                       size.factors.content = NULL,
                                       subset.row = NULL,
                                       block = 10000,
                                       d = 50,
                                       force.match=FALSE,
                                       force.k=20,
                                       force.ndist=3,
                                       BNPARAM=BNPARAM,
                                       BSPARAM=BSPARAM,
                                       BPPARAM=BPPARAM
                                       ), ncol=1)
  scores <- cbind(scores,log10(scores[,1]+1))
  colnames(scores) <- c("scran_doubletCells_score", "scran_doubletCells_score_log10")


  return(scores)
}


#' @title Detect doublet cells using \link[scran]{doubletCells}.
#' @description A wrapper function for \link[scran]{doubletCells}. Identify
#'  potential doublet cells based on simulations of putative doublet expression
#'  profiles. Generate a doublet score for each cell.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[scran]{doubletCells} will be run on cells from each sample separately.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @param nNeighbors Number of nearest neighbors used to calculate density for
#'  doublet detection. Default 50.
#' @param simDoublets Number of simulated doublets created for doublet
#'  detection. Default 10000.
#' @param seed Seed for the random number generator. Default 12345.
#' @param size.factors.norm A numeric vector of size factors for normalization
#'  of \code{x} prior to PCA and distance calculations. If \code{NULL}, defaults
#'  to size factors derived from the library sizes of \code{x}. For the SingleCellExperiment
#'  method, the default values are taken from \code{\link{sizeFactors}(x)}, if they are available.
#' @param size.factors.content A numeric vector of size factors for RNA content
#'  normalization of \code{x} prior to simulating doublets. #' This is orthogonal to
#'  the values in \code{size.factors.norm}
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param block An integer scalar controlling the rate of doublet generation,
#'  to keep memory usage low.
#' @param d An integer scalar specifying the number of components to retain after the PCA.
#' @param force.match A logical scalar indicating whether remapping of simulated
#'  doublets to original cells should be performed.
#' @param force.k An integer scalar specifying the number of neighbours to use for
#'  remapping if \code{force.match=TRUE}.
#' @param force.ndist A numeric scalar specifying the bandwidth for remapping
#'  if \code{force.match=TRUE}.
#' @param BNPARAM A \code{\link[BiocNeighbors]{BiocNeighborParam}} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link[BiocNeighbors]{findNeighbors}}.
#' @param BSPARAM A \code{\link[BiocSingular]{BiocSingularParam}} object specifying the algorithm to
#'  use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object specifying whether the
#'  neighbour searches should be parallelized.
#' @details This function is a wrapper function for \link[scran]{doubletCells}.
#'  \code{runDoubletCells} runs \link[scran]{doubletCells} for each
#'  \code{sample} within \code{inSCE} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link[SummarizedExperiment]{colData} of \code{inSCE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  'scran_doubletCells_score' column added to the
#'  \link[SummarizedExperiment]{colData} slot.
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
    size.factors.norm = NULL,
    size.factors.content = NULL,
    subset.row = NULL,
    block = 10000,
    d = 50,
    force.match=FALSE,
    force.k=20,
    force.ndist=3,
    BNPARAM=BiocNeighbors::KmknnParam(),
    BSPARAM=BiocSingular::bsparam(),
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
  output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            scran_doubletCells_score = numeric(ncol(inSCE)),
            scran_doubletCells_score_log10 = numeric(ncol(inSCE)))

  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
    sceSample <- inSCE[, sceSampleInd]

    mat <- SummarizedExperiment::assay(sceSample, i = useAssay)

    result <- withr::with_seed(seed,
              .runDoubletCells(cell.matrix = mat,
                               k = nNeighbors,
                               nIters = simDoublets,
                               size.factors.norm = NULL,
                               size.factors.content = NULL,
                               subset.row = NULL,
                               block = 10000,
                               d = 50,
                               force.match=FALSE,
                               force.k=20,
                               force.ndist=3,
                               BNPARAM=BNPARAM,
                               BSPARAM=BSPARAM,
                               BPPARAM=BPPARAM
                               ))

    output[sceSampleInd, ] <- result
  }

  argsList <- argsList[!names(argsList) %in% c("BNPARAM","BSPARAM","BPPARAM")]
  #dotList <- list(...)
  #dotList <- dotList[!names(dotList) %in% c("BNPARAM","BSPARAM","BPPARAM")]
  #argsList <- c(argsList, dotList)
  inSCE@metadata$runDoubletCells <- argsList[-1]
  inSCE@metadata$runDoubletCells$packageVersion <- utils::packageDescription("scran")$Version
  colData(inSCE) = cbind(colData(inSCE), output)

  return(inSCE)
}

