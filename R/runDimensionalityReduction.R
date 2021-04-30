#' Wrapper function to run one of the three available dimensionality
#' reduction algorithms integrated within the singleCellTK from 'scaterPCA',
#' 'seuratPCA' and 'seuratICA'.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used.
#' @param reducedDimName Specify the name of the output reducedDim.
#' @param method Specify a method from 'scaterPCA', 'seuratPCA' or 'seuratICA'.
#' @param nComponents Specify the number of dimensions to use with 'seuratPCA' or 
#'  'seuratICA' methods. Not applicable for 'scaterPCA' method. Default
#'  is \code{10}. 
#' @return A \linkS4class{SingleCellExperiment} object with PCA computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
runDimensionalityReduction <- function(inSCE,
                                       useAssay,
                                       reducedDimName,
                                       method = 
                                         c("scaterPCA", 
                                           "seuratPCA", 
                                           "seuratICA"),
                                       nComponents = 10){
  
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