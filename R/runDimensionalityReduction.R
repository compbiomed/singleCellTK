#' Wrapper function to run one of the available dimensionality
#' reduction algorithms integrated within the toolkit from `scaterPCA`,
#' `seuratPCA`, `seuratICA`, `rTSNE`, `seuratTSNE`, `uwotUMAP` and `seuratUMAP`.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used.
#' @param reducedDimName Specify the name of the output reducedDim.
#' @param method Specify a method from `scaterPCA`, `seuratPCA`, `seuratICA`, 
#'  `rTSNE`, `seuratTSNE`, `uwotUMAP` and `seuratUMAP`.
#' @param nComponents Specify the number of dimensions to compute with the
#'  selected method.
#' @param ... Additional parameters for the selected method. For `rTSNE`, must
#'  specify `perplexity` (default \code{30}) and  `nIterations` 
#'  (default \code{1000}). For `seuratTSNE`, must specify `perplexity` (default
#'  \code{30}). For `uwotUMAP`, must specify `nNeighbors` (default \code{30}), 
#'  `nIterations` (default \code{200}), `minDist` (default \code{0.01}) and 
#'  `alpha` (default \code{1}). For `seuratUMAP`, must specify `minDist` 
#'  (default \code{0.3}), `nNeighbors` (default \code{30}) and `spread` 
#'  (default \code{1}).
#' @return A \linkS4class{SingleCellExperiment} object with PCA computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
runDimensionalityReduction <- function(inSCE,
                                       useAssay,
                                       reducedDimName,
                                       method = 
                                         c("scaterPCA", 
                                           "seuratPCA", 
                                           "seuratICA",
                                           "rTSNE",
                                           "seuratTSNE",
                                           "uwotUMAP",
                                           "seuratUMAP"),
                                       nComponents = 10,
                                       ...){
  
  tempSCE <- inSCE
  
  params <- list(
    inSCE = tempSCE,
    useAssay = useAssay,
    reducedDimName = reducedDimName
  )
  
  if(useAssay %in% altExpNames(inSCE)){
    if(method %in% c("seuratPCA", "seuratICA")){
      tempSCE <- altExps(tempSCE)[[useAssay]]
      params$inSCE <- tempSCE
    }
    else{
      params$useAltExp = useAssay
    }
  }
  
  if(method %in% c("seuratPCA", "seuratICA")){
    tempSCE <- seuratFindHVG(
      inSCE = tempSCE,
      useAssay = useAssay
    )
    params$inSCE <- tempSCE
    if(method == "seuratPCA") params$nPCs <- nComponents
    if(method == "seuratICA") params$nics <- nComponents
  }
  else{
    params$ndim <- nComponents
  }

  tempSCE <- do.call(method, args = params)
  
  if(useAssay %in% altExpNames(inSCE)){
    if(method %in% c("seuratPCA", "seuratICA")){
      altExps(inSCE)[[useAssay]] <- tempSCE
      reducedDim(inSCE, reducedDimName) <- reducedDim(tempSCE, reducedDimName)
    }
    else{
      inSCE <- tempSCE
    }
  }
  else{
    inSCE <- tempSCE
  }
  
  return(inSCE)
}