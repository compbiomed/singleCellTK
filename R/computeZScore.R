#' Compute Z-Score
#'
#' Computes Z-Score from an input count matrix using the formula
#' ((x-mean(x))/sd(x)) for each gene across all cells. The input count matrix
#' can either be a base matrix, dgCMatrix or a DelayedMatrix. Computations are
#' performed using DelayedMatrixStats package to efficiently compute the
#' Z-Score matrix.
#' @param counts matrix (base matrix, dgCMatrix or DelayedMatrix)
#' @return z-score computed counts matrix (DelayedMatrix)
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' assay(sce_chcl, "countsZScore") <- computeZScore(assay(sce_chcl, "counts"))
computeZScore <- function(counts) {
    if (!methods::is(counts, "DelayedArray")) {
        counts <- DelayedArray::DelayedArray(counts)
    }
    counts <- (counts - DelayedMatrixStats::rowMeans2(counts)) / DelayedMatrixStats::rowSds(counts)
    counts[base::is.nan(counts)] <- 0
    return(counts)
}
