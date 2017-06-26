
#' Plot PCA
#'   
#' Use this function to plot PCA or tSNE results
#'
#' @param count_data A SCE object
#' @param pca_df PCA data frame
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param pcX User choice for the first principal component
#' @param pcY User choice for the second prinicipal component
#' 
#' @return A PCA plot
#' @export plotPCA
#'

plotPCA <- function(count_data, pca_df=NULL, colorBy=NULL, shape=NULL, pcX="PC1", pcY="PC2"){
  if(is.null(pca_df)){
    pca_df <- getPCA(count_data)
  }
  if(colorBy == "No Color"){
    colorBy <- NULL
  }
  if(shape == "No Shape"){
    shape <- NULL
  }
  l <- pca_df
  if(!is.null(colorBy)){
    l$color <- eval(parse(text = paste("pData(count_data)$",colorBy,sep="")))
  }
  if(!is.null(shape)){
    l$shape <- factor(eval(parse(text = paste("pData(count_data)$",shape,sep=""))))
  }
  l$Sample <- rownames(pData(count_data))
  variances <- attr(pca_df,"percentVar")
  g <- ggplot(l, aes_string(pcX, pcY, label="Sample")) +
       geom_point()+
       labs(x = paste(pcX,toString(round(variances[strtoi(strsplit(pcX,"PC")[[1]][-1])]*100,2)),"%"),y = paste(pcY,toString(round(variances[strtoi(strsplit(pcY,"PC")[[1]][-1])]*100,2)),"%"))
  if(!is.null(colorBy)){
    g <- g + aes_string(color="color") +
      labs(color = colorBy)
  }
  if(!is.null(shape)){
    g <- g + aes_string(shape="shape") +
      labs(shape = shape)
  }
  return(g)
}
