
#' Plot PCA or tSNE
#'   
#' Use this function to plot PCA or tSNE results
#'
#' @param method Either PCA or tSNE
#' @param vals_method PCA or tSNE data frame
#' @param count_data A SCE object
#' @param colorClusters The variable to color clusters by
#' @param pc1 User choice for the first principal component
#' @param pc2 User choice for the second prinicipal component
#' 
#' @return A reduced dimension object
#' @export plotDimRed
#'

plotDimRed <- function(method, vals_method, count_data, colorClusters, pc1, pc2){
  l <- vals_method
  w <- colorClusters
  l$Treatment <- eval(parse(text = paste("pData(count_data)$",w,sep="")))
  l$Sample <- rownames(pData(count_data))
  if(method == "PCA"){
    variances <- attr(vals_method,"percentVar")
    g <- ggplot(l, aes_string(pc1, pc2, label="Sample", color="Treatment")) + 
      geom_point()+labs(x = paste(pc1,toString(round(variances[strtoi(strsplit(pc1,"PC")[[1]][-1])]*100,2)),"%"),y = paste(pc2,toString(round(variances[strtoi(strsplit(pc2,"PC")[[1]][-1])]*100,2)),"%"))
  } else if(method == "tSNE"){
        g <- ggplot(l, aes(X1, X2, label=Sample, color=Treatment))+geom_point()
  }
  return(g)
}
