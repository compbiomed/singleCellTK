#' Trim Counts
#'
#' Trims an input count matrix such that each value greater than a
# 'provided upper threshold value is trimmed to the upper
#' threshold value and each value less than a provided
#' lower threshold value is trimmed to the lower treshold value.
#' @param counts matrix
#' @param trimValue where trimValue[1] for upper threshold and trimValue[2] as
#' lower threshold. Default is \code{c(10,-10)}
#' @return trimmed counts matrix
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' assay(sce_chcl, "countsTrimmed") <- trimCounts(assay(sce_chcl, "counts"),
#'                                                c(10, -10))
trimCounts <- function(counts, trimValue = c(10,-10)) {
    trimValue <- c(max(trimValue), min(trimValue))
    counts[counts > trimValue[1]] <- trimValue[1]
    counts[counts < trimValue[2]] <- trimValue[2]
    return(counts)
}
