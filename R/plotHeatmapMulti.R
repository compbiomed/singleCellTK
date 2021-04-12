#' plotHeatmapMulti
#' 
#' @param plots Input list of plot objects.
#' @param components Specify the components against which the heatmaps should
#' be plotted. Default is \code{NULL} which will include all available
#' in input list of plots.
#' @param nCol Specify the number of columns in the output plot. Default
#' is \code{NULL} which will auto-compute the number of columns.
#' @return Heatmap plot object.
#' @export
plotHeatmapMulti <- function(plots,
                             components = NULL,
                             nCol = NULL){
  plot <- NULL
  
  if(is.null(nCol)){
    nCol <- floor(sqrt(length(plots)))
  }
  
  if(is.null(components)){
    plot <- cowplot::plot_grid(plotlist = plots, ncol = nCol) 
  }
  else{
    if(!is.numeric(components)){
      components <- as.numeric(gsub("[^0-9.-]", "", components))
    }
    plot <- cowplot::plot_grid(plotlist = plots[components], ncol = nCol) 
  }
  
  return(plot)
}