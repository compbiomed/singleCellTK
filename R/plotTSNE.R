
#' Plot TSNE
#'   
#' Use this function to plot PCA or tSNE results
#'
#' @param count_data A SCE object
#' @param tsne_df tSNE data frame
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' 
#' @return A TSNE plot
#' @export plotTSNE
#'

plotTSNE <- function(count_data, tsne_df=NULL, colorBy=NULL, shape=NULL){
  if (is.null(tsne_df)){
    tsne_df <- getTSNE(count_data)
  }
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  l <- tsne_df
  if (!is.null(colorBy)){
    l$color <- eval(parse(text = paste("pData(count_data)$",colorBy,sep = "")))
  }
  if (!is.null(shape)){
    l$shape <- factor(eval(parse(text = paste("pData(count_data)$",shape, sep = ""))))
  }
  l$Sample <- rownames(pData(count_data))
  g <- ggplot(l, aes_string("X1", "X2", label = "Sample")) +
    geom_point()
  if (!is.null(colorBy)){
    g <- g + aes_string(color = "color") +
      labs(color = colorBy)
  }
  if (!is.null(shape)){
    g <- g + aes_string(shape = "shape") +
      labs(shape = shape)
  }
  return(g)
}
