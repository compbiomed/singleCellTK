#' Plot dimensionality reduction from computed metrics including PCA, ICA, tSNE
#' and UMAP
#'
#' @param inSCE Input SCE object
#' @param useReduction Reduction to plot
#' @param showLegend If legends should be plotted or not
#' @param xDim Numeric value indicating the dimension to use for X-axis.
#'  Default is 1 (refers to PC1).
#' @param yDim Numeric value indicating the dimension to use for Y-axis.
#'  Default is 2 (refers to PC2).
#' @param xAxisLabel Specify the label for x-axis. Default is \code{NULL} which
#' will specify the label as 'x'.
#' @param yAxisLabel Specify the label for y-axis. Default is \code{NULL} which
#' will specify the label as 'y'.
#' @return plot object
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' plotDimRed(mouseBrainSubsetSCE, "PCA_logcounts")
plotDimRed <- function(inSCE, useReduction,
                       showLegend = FALSE,
                       xDim = 1,
                       yDim = 2,
                       xAxisLabel = NULL,
                       yAxisLabel = NULL) {
  dimRed <- reducedDim(inSCE, type = useReduction)

  x <- dimRed[, xDim]
  y <- dimRed[, yDim]

  dimRedPlot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y))
  dimRedPlot <- .ggSCTKTheme(dimRedPlot)

  if(!is.null(xAxisLabel)){
    dimRedPlot <- dimRedPlot + ggplot2::xlab(xAxisLabel)
  }
  else{
    dimRedPlot <- dimRedPlot + ggplot2::xlab(colnames(dimRed)[xDim])
  }

  if(!is.null(yAxisLabel)){
    dimRedPlot <- dimRedPlot + ggplot2::ylab(yAxisLabel)
  }
  else{
    dimRedPlot <- dimRedPlot + ggplot2::ylab(colnames(dimRed)[yDim])
  }

  return(dimRedPlot)
}
