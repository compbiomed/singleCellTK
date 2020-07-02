#' Plot PCA run data from its components.
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param pcX User choice for the first principal component
#' @param pcY User choice for the second principal component
#' @param runPCA Run PCA if the reducedDimName does not exist. the Default is
#' FALSE.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "logcounts".
#' @param reducedDimName a name to store the results of the dimension reduction
#' coordinates obtained from this method. This is stored in the SingleCellExperiment
#' object in the reducedDims slot. Required.
#'
#' @return A PCA plot
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotPCA(mouseBrainSubsetSCE, colorBy = "level1class",
#'         reducedDimName = "PCA_counts")
#'
plotPCA <- function(inSCE, colorBy="No Color", shape="No Shape", pcX="PC1",
                    pcY="PC2", reducedDimName="PCA", runPCA=FALSE,
                    useAssay="logcounts"){
  if(!(reducedDimName %in% names(SingleCellExperiment::reducedDims(inSCE)))){
    if (runPCA){
      inSCE <- getPCA(inSCE, useAssay = useAssay,
                      reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName,
           " dimension not found. Run getPCA() or set runPCA to TRUE.")
    }
  }
  pcaDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                       reducedDimName))
  if (!(pcX %in% colnames(pcaDf))){
    stop("pcX dimension ", pcX, " is not in the reducedDim data")
  }
  if (!(pcY %in% colnames(pcaDf))){
    stop("pcY dimension ", pcY, " is not in the reducedDim data")
  }

  
  # Need to add back in variances in the plot axis labels
  pcXlab <- pcX
  pcYlab <- pcY

  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    pcaDf$color <- SingleCellExperiment::colData(inSCE)[, colorBy]
  }
  if (!is.null(shape)){
    pcaDf$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
  }
  pcaDf$Sample <- colnames(inSCE)
  g <- ggplot2::ggplot(pcaDf, ggplot2::aes_string(pcX, pcY, label = "Sample")) +
    ggplot2::geom_point() +
    ggplot2::labs(x = pcXlab, y = pcYlab)
  if (!is.null(colorBy)){
    g <- g + ggplot2::aes_string(color = "color") +
      ggplot2::labs(color = colorBy)
  }
  if (!is.null(shape)){
    g <- g + ggplot2::aes_string(shape = "shape") +
      ggplot2::labs(shape = shape)
  }
  return(g)
}
