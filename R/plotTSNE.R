#' Plot TSNE
#'
#' Use this function to plot PCA or tSNE results
#'
#' @param count_data A SCE object
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#'
#' @return A TSNE plot
#' @export plotTSNE
#'
plotTSNE <- function(count_data, colorBy="No Color", shape="No Shape"){
  if (is.null(reducedDim(count_data, "TSNE"))){
    count_data <- getTSNE(count_data)
  }
  tsne_df <- data.frame(reducedDim(count_data, "TSNE"))
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    tsne_df$color <- eval(parse(text = paste("colData(count_data)$", colorBy, sep = "")))
  }
  if (!is.null(shape)){
    tsne_df$shape <- factor(eval(parse(text = paste("colData(count_data)$", shape, sep = ""))))
  }
  tsne_df$Sample <- colnames(count_data)
  g <- ggplot2::ggplot(tsne_df, ggplot2::aes_string("X1", "X2", label = "Sample")) +
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
