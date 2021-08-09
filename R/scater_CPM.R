#' scaterCPM
#' Uses CPM from scater library to compute counts-per-million.
#' @param inSCE Input SingleCellExperiment object
#' @param assayName New assay name for cpm data.
#' @param useAssay Input assay 
#' @return inSCE Updated SingleCellExperiment object
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scaterCPM(sce_chcl,"countsCPM", "counts")
scaterCPM <- function(inSCE, 
                       assayName = "ScaterCPMCounts", 
                       useAssay = "counts"){
  assay(inSCE, assayName) <- scater::calculateCPM(
    x = assay(inSCE, useAssay))
  
  inSCE <- expSetDataTag(inSCE = inSCE, 
                         assayType = "normalized", 
                         assays = assayName)
  return(inSCE)
}
