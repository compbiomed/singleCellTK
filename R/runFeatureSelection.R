#' Wrapper function to run all of the feature selection methods integrated
#' within the singleCellTK package including three methods from 
#' Seurat (`vst`, `mean.var.plot` or `dispersion`) and the Scran `modelGeneVar`
#' method.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used.
#' @param hvgMethod Specify the method to use for variable gene selection.
#'  Options include `vst`, `mean.var.plot` or `dispersion` from Seurat
#'  and `modelGeneVar` from Scran.
#' @return A \code{SingleCellExperiment} object that contains the computed
#'  statistics in the \code{rowData} slot of the output object. This function 
#'  does not return the names of the variable features but only computes the 
#'  statistics that are stored in the \code{rowData} slot of the. To get the 
#'  names of the variable features \code{getTopHVG} function should be used 
#'  after computing these statistics.
#' @export
#'
#' @examples
runFeatureSelection <- function(inSCE,
                                useAssay,
                                hvgMethod
                                ){
  
  seuratMethods <- c("vst", "mean.var.plot", "dispersion")
  scranMethods <- c("modelGeneVar")
  
  params <- list(
    inSCE = inSCE
  )
  
  if(hvgMethod %in% seuratMethods){
    params$useAssay <- useAssay
    params$hvgMethod <- hvgMethod
    inSCE <- do.call("seuratFindHVG", args = params)
  }
  else if(hvgMethod %in% scranMethods){
    params$assayName <- useAssay
    inSCE <- do.call("scranModelGeneVar", args = params)
  }
  
  return(inSCE)
}