#' Plot dimensionality reduction
#'
#' @param inSCE Input SCE object
#' @param useReduction Reduction to plot
#' @param showLegend If legends should be plotted or not
#'
#' @return plot object
#' @export
plotDimRed <- function(inSCE, useReduction, showLegend = FALSE) {
  dimRed <- reducedDim(inSCE, type = useReduction)
  x <- dimRed[,1]
  y <- dimRed[,2]
  dimRedPlot <- ggplot() +
    geom_point(aes(x = x, y = y))
  return(dimRedPlot)
}