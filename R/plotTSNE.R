#' @describeIn getPCA plot t-SNE results
#'
#' @param runTSNE Run t-SNE if the reducedDimName does not exist. the Default is
#' FALSE.
#'
#' @return plotTSNE(): A t-SNE plot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotTSNE(mouseBrainSubsetSCE, colorBy = "level1class",
#'          reducedDimName = "TSNE_counts")
#'
plotTSNE <- function(inSCE, colorBy="No Color", shape="No Shape",
                     reducedDimName="TSNE", runTSNE=FALSE,
                     useAssay="logcounts"){
  if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))){
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
