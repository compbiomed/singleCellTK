#' scater_logNormCounts
#' Uses \link[scater]{logNormCounts} to log normalize input data
#' @param inSCE Input SingleCellExperiment object
#' @param logAssayName New assay name for log normalized data
#' @param useAssay Input assay 
#' @return inSCE Updated SingleCellExperiment object that contains the new log normalized data
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scater_logNormCounts(sce_chcl,"logcounts", "counts")
scater_logNormCounts <- function(inSCE, logAssayName = "ScaterLogNormCounts", useAssay = "counts"){
  inSCE <- scater::logNormCounts(
    x = inSCE, 
    name = logAssayName,
    exprs_values = useAssay)
  return(inSCE)
}
