#' Plot t-SNE plot on dimensionality reduction data run from t-SNE method.
#'
#' @param runTSNE Run t-SNE if the reducedDimName does not exist. the Default is
#' FALSE.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "logcounts".
#' @param reducedDimName a name to store the results of the dimension reduction
#' coordinates obtained from this method. This is stored in the SingleCellExperiment
#' object in the reducedDims slot. Required.
#' @param colorBy color by condition.
#' @param shape add shape to each distinct label.
#'
#' @return A t-SNE plot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotTSNE(mouseBrainSubsetSCE, colorBy = "level1class",
#'          reducedDimName = "TSNE_counts")
#'
plotTSNE <- function(inSCE, colorBy="No Color", shape="No Shape",
                     reducedDimName="TSNE", runTSNE=FALSE,
                     useAssay="logcounts"){
  if(!(reducedDimName %in% names(SingleCellExperiment::reducedDims(inSCE)))){
    if (runTSNE){
      inSCE <- getTSNE(inSCE, useAssay = useAssay,
                       reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName,
           " dimension not found. Run getTSNE() or set runTSNE to TRUE.")
    }
  }
  tsneDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                        reducedDimName))
  if (ncol(tsneDf) > 2){
    warning("More than two t-SNE dimensions. Using the first two.")
  }
  colnames(tsneDf)[1] <- "tSNE1"
  colnames(tsneDf)[2] <- "tSNE2"
  xdim <- colnames(tsneDf)[1]
  ydim <- colnames(tsneDf)[2]
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    tsneDf$color <- SingleCellExperiment::colData(inSCE)[, colorBy]
  }
  if (!is.null(shape)){
    tsneDf$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
  }
  tsneDf$Sample <- colnames(inSCE)
  g <- ggplot2::ggplot(tsneDf, ggplot2::aes_string(xdim, ydim,
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
