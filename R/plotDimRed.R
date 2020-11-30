#' Plot dimensionality reduction from computed metrics including PCA, ICA, tSNE
#' and UMAP
#'
#' @param inSCE Input SCE object
#' @param useReduction Reduction to plot
#' @param showLegend If legends should be plotted or not
#' @param xAxisLabel Specify the label for x-axis. Default is \code{NULL} which
#' will specify the label as 'x'.
#' @param yAxisLabel Specify the label for y-axis. Default is \code{NULL} which
#' will specify the label as 'y'.
#'
#' @return plot object
#' @export
plotDimRed <- function(inSCE, useReduction, 
                       showLegend = FALSE,
                       xAxisLabel = NULL,
                       yAxisLabel = NULL) {
  dimRed <- reducedDim(inSCE, type = useReduction)
  
  x <- dimRed[,1]
  y <- dimRed[,2]
  
  dimRedPlot <- ggplot() +
    geom_point(aes(x = x, y = y))
  
  if(!is.null(xAxisLabel)){
    dimRedPlot <- dimRedPlot + xlab(xAxisLabel)
  }
  
  if(!is.null(yAxisLabel)){
    dimRedPlot <- dimRedPlot + ylab(yAxisLabel)
  }
  
  return(dimRedPlot)
}