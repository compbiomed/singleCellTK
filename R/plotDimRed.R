#' Plot results either on already run results of reduced dimensions data.
#'
#' @param inSCE Input SCtkExperiment object with saved dimension reduction components
#'  or a variable with saved results. Required
#' @param colorBy color by a condition(any column of the annotation data).
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment object. Required.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param comp1 label for x-axis
#' @param comp2 label for y-axis
#' @param pcX PCA component to be used for plotting(if applicable).
#' Default is first PCA component for PCA data and NULL otherwise.
#' @param pcY PCA component to be used for plotting(if applicable).
#' Default is second PCA component for PCA data and NULL otherwise.
#'
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotDimRed(inSCE = mouseBrainSubsetSCE, colorBy = "No Color", shape = "No Shape",
#'            reducedDimName = "TSNE_counts", useAssay = "counts",
#'            comp1 = "tSNE1", comp2 = "tSNE2")
#'
plotDimRed <- function(inSCE, colorBy = "No Color", shape = "No Shape",
                     reducedDimName = NULL,
                     useAssay = "logcounts", comp1 = NULL, comp2 = NULL,
                     pcX = NULL, pcY = NULL) {
  Df <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                        reducedDimName))
  if (ncol(Df) > 2){
    warning("More than two dimensions. Using the first two.")
  }
  if (!is.null(pcX) & !is.null(pcY)){
    if (!(pcX %in% colnames(Df))){
      stop("X dimension ", pcX, " is not in the reducedDim data")
    }
    if (!(pcY %in% colnames(Df))){
      stop("Y dimension ", pcY, " is not in the reducedDim data")
    }
    xdim <- pcX
    ydim <- pcY
  } else if (!is.null(comp1) & !is.null(comp2)){
    colnames(Df)[1] <- comp1
    colnames(Df)[2] <- comp2
    xdim <- colnames(Df)[1]
    ydim <- colnames(Df)[2]
  } else {
    xdim <- colnames(Df)[1]
    ydim <- colnames(Df)[2]
  }

  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    Df$color <- SingleCellExperiment::colData(inSCE)[, colorBy]
  }
  if (!is.null(shape)){
    Df$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
  }
  Df$Sample <- colnames(inSCE)
  g <- ggplot2::ggplot(Df, ggplot2::aes_string(xdim, ydim,
                                                   label = "Sample")) +
    ggplot2::geom_point()
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
