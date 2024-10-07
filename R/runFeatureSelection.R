#' Run Variable Feature Detection Methods
#' @description Wrapper function to run all of the feature selection methods 
#' integrated within the singleCellTK package including three methods from
#' Seurat (\code{"vst"}, \code{"mean.var.plot"} or \code{dispersion}) and the 
#' Scran \code{modelGeneVar} method.
#' 
#' This function does not return the names of the variable features but only 
#' computes the metrics, which will be stored in the \code{rowData} slot. To set
#' a HVG list for downstream use, users should call \code{\link{setTopHVG}} 
#' after computing the metrics. To get the names of the variable features, users
#' should call \code{\link{getTopHVG}} function after computing the metrics.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used. Should use
#' raw counts for \code{"vst"} method, or a normalized assay for other methods.
#' @param method Specify the method to use for variable gene selection.
#' Options include \code{"vst"}, \code{"mean.var.plot"} or \code{"dispersion"}
#' from Seurat and \code{"modelGeneVar"} from Scran. Default \code{"vst"}
#' @return The input \linkS4class{SingleCellExperiment} object that contains 
#' the computed statistics in the \code{rowData} slot
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
                                method = "vst")
  {
  method <- match.arg(method, choices = c("vst", "mean.var.plot", "dispersion", "modelGeneVar", "cell_ranger"))
  seuratMethods <- c("vst", "mean.var.plot", "dispersion")
  scranMethods <- c("modelGeneVar")
  scanpyMethods <- c("cell_ranger")

  params <- list(
    inSCE = inSCE,
    useAssay = useAssay
  )

  if(method %in% seuratMethods){
    params$method <- method
    inSCE <- do.call("runSeuratFindHVG", args = params)
  }
  else if(method %in% scranMethods){
    inSCE <- do.call("runModelGeneVar", args = params)
  }
  else if(method %in% scanpyMethods){
    params$method <- method
    inSCE <- do.call("runScanpyFindHVG", args = params)
  }

  return(inSCE)
}
