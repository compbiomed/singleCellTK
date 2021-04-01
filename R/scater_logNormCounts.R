#' scaterlogNormCounts
#' Uses \link{logNormCounts} to log normalize input data
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
  inSCE <- scater::logNormCounts(
    x = inSCE, 
    name = assayName,
    exprs_values = useAssay)
  
  inSCE <- expSetDataTag(inSCE = inSCE, 
                         assayType = "normalized", 
                         assays = assayName)
  return(inSCE)
}
