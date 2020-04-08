#' Trim Counts
#'
#' Trims an input count matrix such that each value greater than a 
# 'provided upper threshold value (positive) is trimmed to the upper 
#' threshold value (positive) and each value less than a provided
#' lower threshold value (negative) is trimmed to the lower treshold 
#' value (negative).
#'
#' @param counts matrix
#' @param trimValue a value set as +trimValue for upper threshold and -trimValue as lower threshold
#'        (default is 10)	
#'
#' @return trimmed counts matrix
#' @export
#'
#' @examples
#' 
#' library(TENxPBMCData)
#' sce <- TENxPBMCData("pbmc3k")
#' rownames(sce) <- rowData(sce)$Symbol_TENx
#' colnames(sce) <- colData(sce)$Barcode
#' assay(sce, "countsTrimmed") <- trimCounts(assay(sce, "counts"), 10)
#'
trimCounts <- function(counts, trimValue = 10) {
    counts[counts > trimValue] <- trimValue
    counts[counts < -(trimValue)] <- -(trimValue)
    return(counts)
}