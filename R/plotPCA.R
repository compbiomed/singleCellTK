#' Plot PCA
#'
#' Use this function to plot PCA or tSNE results
#'
#' @param count_data A SCE object
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param pcX User choice for the first principal component
#' @param pcY User choice for the second prinicipal component
#' @param use_assay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName PCA dimension name. The default is PCA.
#' The toolkit will store data with the pattern <ASSSAY>_<ALGORITHM>.
#'
#' @return A PCA plot
#' @export plotPCA
#'
plotPCA <- function(count_data, colorBy="No Color", shape="No Shape", pcX="PC1",
                    pcY="PC2", use_assay="logcounts", reducedDimName="PCA"){
  if (is.null(reducedDim(count_data, reducedDimName))){
    count_data <- getPCA(count_data, use_assay = use_assay, reducedDimName = reducedDimName)
  }
  pca_df <- data.frame(reducedDim(count_data, reducedDimName))
  variances <- pca_variances(count_data)$percentVar
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    pca_df$color <- eval(parse(text = paste("colData(count_data)$", colorBy, sep = "")))
  }
  if (!is.null(shape)){
    pca_df$shape <- factor(eval(parse(text = paste("colData(count_data)$", shape, sep = ""))))
  }
  pca_df$Sample <- colnames(count_data)
  g <- ggplot2::ggplot(pca_df, ggplot2::aes_string(pcX, pcY, label = "Sample")) +
    ggplot2::geom_point() +
    ggplot2::labs(x = paste(pcX, toString(round(variances[strtoi(strsplit(pcX, "PC")[[1]][-1])] * 100, 2)), "%"),
                  y = paste(pcY, toString(round(variances[strtoi(strsplit(pcY, "PC")[[1]][-1])] * 100, 2)), "%"))
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
