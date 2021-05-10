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
#'  selected method. Only applicable with `scaterPCA`, `seuratPCA`, `seuratICA`,
#'  `seuratTSNE` and `seuratTSNE` methods.
#' @param ... Additional parameters for the selected method. For `rTSNE`, must
#'  specify `perplexity` (default \code{30}) and  `nIterations`
#'  (default \code{1000}). For `seuratTSNE`, must specify `useReduction`
#'  (either `pca` or `ica`) and `perplexity` (default \code{30}).
#'  For `uwotUMAP`, must specify `nNeighbors` (default \code{30}),
#'  `nIterations` (default \code{200}), `minDist` (default \code{0.01}) and
#'  `alpha` (default \code{1}). For `seuratUMAP`, must specify `useReduction`
#'  (either `pca` or `ica`), `minDist` (default \code{0.3}), `nNeighbors`
#'  (default \code{30}) and `spread` (default \code{1}).
#' @return A \linkS4class{SingleCellExperiment} object with PCA computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runDimensionalityReduction(mouseBrainSubsetSCE,
#'                                                   "logcounts",
#'                                                   reducedDimName = "PCA")
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
  method <- match.arg(method)
  args <- list(...)
  tempSCE <- inSCE

  params <- list(
    inSCE = tempSCE,
    useAssay = useAssay,
    reducedDimName = reducedDimName
  )

  params <- c(params, args)

  if(useAssay %in% altExpNames(inSCE)){
    if(method %in% c("seuratPCA", "seuratICA", "seuratTSNE", "seuratUMAP")){
      tempSCE <- altExps(tempSCE)[[useAssay]]
      params$inSCE <- tempSCE
    }
    else{
      params$useAltExp = useAssay
    }
  }

  if(method %in% c("seuratPCA", "seuratICA", "seuratTSNE", "seuratUMAP")){
    if(useAssay %in% altExpNames(inSCE)){
      tempSCE <- seuratFindHVG(
        inSCE = tempSCE,
        useAssay = useAssay,
        altExp = TRUE
      )
    }
    else{
      tempSCE <- seuratFindHVG(
        inSCE = tempSCE,
        useAssay = useAssay
      )
    }
    params$inSCE <- tempSCE
    if(method == "seuratPCA") params$nPCs <- nComponents
    if(method == "seuratICA") params$nics <- nComponents
  }
  if(method == "scaterPCA"){
    params$ndim <- nComponents
  }
  if(method == "rTSNE"){
    method <- "getTSNE"
  }
  if(method == "seuratTSNE"){
    method <- "seuratRunTSNE"
    if(params$useReduction == "pca"){
      params$inSCE <- seuratPCA(inSCE = params$inSCE,
                               useAssay = params$useAssay,
                               reducedDimName = paste0(params$reducedDimName, "_PCA"))
    }
    else{
      params$inSCE <- seuratICA(inSCE = params$inSCE,
                               useAssay = params$useAssay,
                               reducedDimName = paste0(params$reducedDimName, "_ICA"))
    }
    params$useAssay <- NULL
  }
  if(method == "uwotUMAP"){
    method <- "getUMAP"
  }
  if(method == "seuratUMAP"){
    method <- "seuratRunUMAP"
    if(params$useReduction == "pca"){
      params$inSCE <- seuratPCA(inSCE = params$inSCE,
                         useAssay = params$useAssay,
                         reducedDimName = paste0(params$reducedDimName, "_PCA"))
    }
    else{
      params$inSCE <- seuratICA(inSCE = params$inSCE,
                               useAssay = params$useAssay,
                               reducedDimName = paste0(params$reducedDimName, "_ICA"))
    }
    params$useAssay <- NULL
  }

  tempSCE <- do.call(method, args = params)

  if(useAssay %in% altExpNames(inSCE)){
    if(method %in% c("seuratPCA", "seuratICA", "seuratRunTSNE", "seuratRunUMAP")){
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
