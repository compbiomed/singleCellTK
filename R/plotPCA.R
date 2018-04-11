#' Plot PCA
#'
#' Use this function to plot PCA results
#'
#' @param countData A SCE object
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param pcX User choice for the first principal component
#' @param pcY User choice for the second principal component
#' @param reducedDimName PCA dimension name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#' @param runPCA Run PCA if the reducedDimName does not exist. the Default is
#' FALSE.
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"

#'
#' @return A PCA plot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotPCA(mouseBrainSubsetSCE, colorBy = "level1class",
#'         reducedDimName = "PCA_counts")
#'
plotPCA <- function(countData, colorBy="No Color", shape="No Shape", pcX="PC1",
                    pcY="PC2", reducedDimName="PCA", runPCA=FALSE,
                    useAssay="logcounts"){
  if (is.null(SingleCellExperiment::reducedDim(countData, reducedDimName))){
    if (runPCA){
      countData <- getPCA(countData, useAssay = useAssay,
                          reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName,
           " dimension not found. Run getPCA() or set runPCA to TRUE.")
    }
  }
  pcaDf <- data.frame(SingleCellExperiment::reducedDim(countData,
                                                       reducedDimName))
  if (!(pcX %in% colnames(pcaDf))){
    stop("pcX dimension ", pcX, " is not in the reducedDim data")
  }
  if (!(pcY %in% colnames(pcaDf))){
    stop("pcY dimension ", pcY, " is not in the reducedDim data")
  }

  if (class(countData) == "SCtkExperiment"){
    if (all(c(pcX, pcY) %in% rownames(pcaVariances(countData)))){
      #use the variances in pcaVariances
      variances <- pcaVariances(countData)
      pcXlab <- paste0(
        pcX, " ", toString(round(variances[pcX, ] * 100, 2)), "%")
      pcYlab <- paste0(
        pcY, " ", toString(round(variances[pcY, ] * 100, 2)), "%")
    } else {
      #do not use variances in the plot
      pcXlab <- pcX
      pcYlab <- pcY
    }
  } else {
    #do not use variances in the plot
    pcXlab <- pcX
    pcYlab <- pcY
  }

  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    pcaDf$color <- SingleCellExperiment::colData(countData)[, colorBy]
  }
  if (!is.null(shape)){
    pcaDf$shape <- factor(SingleCellExperiment::colData(countData)[, shape])
  }
  pcaDf$Sample <- colnames(countData)
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
