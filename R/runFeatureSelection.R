#' Run Variable Feature Detection Methods
#' @description Wrapper function to run all of the feature selection methods 
#' integrated within the singleCellTK package including three methods from
#' Seurat (\code{"vst"}, \code{"mean.var.plot"} or \code{dispersion}) and the 
#' Scran \code{modelGeneVar} method.
#' 
#' This function does not return the names of the variable features but only 
#' computes the statistics and the selection which will be stored in the 
#' \code{rowData} slot. To get the names of the variable features, users should 
#' call \code{\link{getTopHVG}} function after computing the statistics.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used. Should use
#' raw counts for \code{"vst"} method, or a normalized assay for other methods.
#' @param method Specify the method to use for variable gene selection.
#' Options include \code{"vst"}, \code{"mean.var.plot"} or \code{"dispersion"}
#' from Seurat and \code{"modelGeneVar"} from Scran.
#' @param hvgNumber An integer for the number of top HVG to select. Default 
#' \code{2000}.
#' @param rowSubsetName A character string for the \code{rowData} variable name
#' to store a logical index of selected HVGs. Default 
#' \code{paste0("HVG_", method, hvgNumber)}
#' @return The input \linkS4class{SingleCellExperiment} object that contains 
#' \item{}{The computed statistics in the \code{rowData} slot} 
#' \item{}{A logical vector indicating the selected HVG in the 
#' \code{rowData} slot, named by \code{rowSubsetName}}  
#' @seealso \code{\link{runModelGeneVar}}, \code{\link{runSeuratFindHVG}},
#' \code{\link{getTopHVG}}, \code{\link{plotTopHVG}}
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runFeatureSelection(mouseBrainSubsetSCE,
#'                                            "logcounts",
#'                                            "modelGeneVar")
runFeatureSelection <- function(inSCE,
                                useAssay,
                                method = c("vst", "mean.var.plot",
                                           "dispersion", "modelGeneVar"),
                                hvgNumber = 2000,
                                rowSubsetName = paste0("HVG_", method, 
                                                       hvgNumber)
                                ){
  method <- match.arg(method)
  seuratMethods <- c("vst", "mean.var.plot", "dispersion")
  scranMethods <- c("modelGeneVar")

  params <- list(
    inSCE = inSCE,
    useAssay = useAssay,
    hvgNumber = hvgNumber,
    rowSubsetName = rowSubsetName
  )

  if(method %in% seuratMethods){
    params$method <- method
    inSCE <- do.call("runSeuratFindHVG", args = params)
  }
  else if(method %in% scranMethods){
    inSCE <- do.call("runModelGeneVar", args = params)
  }

  return(inSCE)
}
