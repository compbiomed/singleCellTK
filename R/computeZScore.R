#' Compute Z-Score
#'
#' Computes Z-Score from an input count matrix using the formula ((x-mean(x))/sd(x))
#' for each gene across all cells. The input count matrix can either be a base matrix,
#' dgCMatrix or a DelayedMatrix. Computations are performed using DelayedMatrixStats
#' package to efficiently compute the Z-Score matrix. 
#'
#' @param counts matrix (base matrix, dgCMatrix or DelaymedMatrix)
#'
#' @return z-score computed counts matrix (DelayedMatrix)
#' @export
#'
#' @examples
#' 
#' library(TENxPBMCData)
#' sce <- TENxPBMCData("pbmc3k")
#' rownames(sce) <- rowData(sce)$Symbol_TENx
#' colnames(sce) <- colData(sce)$Barcode
#' assay(sce, "countsZScore") <- computeZScore(assay(sce, "counts"))
#'
computeZScore <- function(counts) {
    if (!methods::is(counts, "DelayedArray")) {
        counts <- DelayedArray(counts)
    }
    counts <- (counts - DelayedMatrixStats::rowMeans2(counts)) / DelayedMatrixStats::rowSds(counts)
    counts[base::is.nan(counts)] <- 0
    return(counts)
}
