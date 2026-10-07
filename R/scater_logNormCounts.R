#' scaterlogNormCounts
#' Log-normalizes counts by library size, as \code{scater::logNormCounts}
#' did, using \code{\link[scrapper]{normalizeRnaCounts.se}}. Existing
#' \code{sizeFactors(inSCE)} are used (after centering) when present.
#' @param inSCE Input SingleCellExperiment object
#' @param assayName New assay name for log normalized data
#' @param useAssay Input assay 
#' @return inSCE Updated SingleCellExperiment object that contains the new log normalized data
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scaterlogNormCounts(sce_chcl,"logcounts", "counts")
scaterlogNormCounts <- function(inSCE, 
                                 assayName = "ScaterLogNormCounts", 
                                 useAssay = "counts"){
  sizeFactors <- SingleCellExperiment::sizeFactors(inSCE)
  if (is.null(sizeFactors)) {
    sizeFactors <- colSums(assay(inSCE, useAssay))
  }
  if (any(sizeFactors <= 0)) {
    stop("size factors should be positive")
  }
  inSCE <- scrapper::normalizeRnaCounts.se(
    inSCE,
    size.factors = sizeFactors,
    assay.type = useAssay,
    output.name = assayName,
    more.norm.args = list(delayed = FALSE)
  )
  newSizeFactors <- SingleCellExperiment::sizeFactors(inSCE)
  SingleCellExperiment::sizeFactors(inSCE) <- stats::setNames(
    newSizeFactors, names(sizeFactors)
  )
  
  inSCE <- expSetDataTag(inSCE = inSCE, 
                         assayType = "normalized", 
                         assays = assayName)
  return(inSCE)
}
