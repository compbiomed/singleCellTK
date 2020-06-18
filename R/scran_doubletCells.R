
.runDoubletCells <- function(cell.matrix, k = k, nIters = nIters, ...) {

  cell.matrix <- .convertToMatrix(cell.matrix)

  scores <- matrix(scran::doubletCells(cell.matrix, k = k,
                                       niters = nIters, ...), ncol=1)
  colnames(scores) <- "scran_doubletCells_score"

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
#' @param ... Additional arguments to pass to \link[scran]{doubletCells}.
#' @details This function is a wrapper function for \link[scran]{doubletCells}.
#'  \code{runDoubletCells} runs \link[scran]{doubletCells} for each
#'  \code{sample} within \code{inSCE} iteratively. The
#'  resulting doublet scores for all cells will be appended to the
#'  \link[SummarizedExperiment]{colData} of \code{inSCE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  'scran_doubletCell_Score' column added to the
#'  \link[SummarizedExperiment]{colData} slot.
#' @references Lun ATL (2018). Detecting doublet cells with scran.
#'  \url{https://ltla.github.io/SingleCellThoughts/software/
#' doublet_detection/bycell.html}
#' @seealso \link[scran]{doubletCells}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- runDoubletCells(sce)
#' @export
runDoubletCells <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    nNeighbors = 50,
    simDoublets = 10000,
    seed = 12345,
    ...
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
            scran_doubletCells_Score = numeric(ncol(inSCE)))

  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
    sceSample <- inSCE[, sceSampleInd]

    mat <- SummarizedExperiment::assay(sceSample, i = useAssay)

    result <- withr::with_seed(seed,
              .runDoubletCells(cell.matrix = mat, k = nNeighbors,
                               nIters = simDoublets, ...))

    output[sceSampleInd, ] <- result
  }

  argsList = argsList[!names(argsList) %in% ("...")]
  inSCE@metadata$runDoubletCells <- argsList[-1]
  inSCE@metadata$runDoubletCells$packageVersion <- utils::packageDescription("scran")$Version
  colData(inSCE) = cbind(colData(inSCE), output)

  return(inSCE)
}

