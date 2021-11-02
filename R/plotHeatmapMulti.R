.plotHeatmapMulti <- function(plots,
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
