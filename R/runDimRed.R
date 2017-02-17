#' Run PCA or tSNE
#'   
#' Use this function to run limma differential expression and load an interactive D3 heatmap
#'
#' @param method Either PCA or tSNE
#' @param count_data A SCE object
#' @param colorClusters The variable to color clusters by
#' @param pc1 User choice for the first principal component
#' @param pc2 User choice for the second prinicipal component
#' 
#' @return A reduced dimension object
#' @export runDimRed
#'

runDimRed <- function(method, count_data, colorClusters, pc1, pc2){
  if(method == "PCA"){
    reducedDimension(count_data) <- reducedDimension(scater::plotPCA(count_data, return_SCESet=TRUE, draw_plot=FALSE))
    variances <- attr(reducedDimension(count_data),"percentVar")
    l <- data.frame(reducedDimension(count_data))
    w <- colorClusters
    d <- c("Treatment")
    l$Treatment <- eval(parse(text = paste("pData(count_data)$",w,sep="")))
    l$Sample <- rownames(pData(count_data))
    g <- ggplot(l, aes_string(pc1, pc2, color=d))+geom_point()+labs(x = paste(pc1,toString(round(variances[strtoi(strsplit(pc1,"PC")[[1]][-1])]*100,2)),"%"),y = paste(pc2,toString(round(variances[strtoi(strsplit(pc2,"PC")[[1]][-1])]*100,2)),"%"))
  } else if(method == "tSNE"){
    tsne <- scater::plotTSNE(count_data, return_SCESet=TRUE)
    g <- reducedDimension(tsne)
    l <- data.frame(g)
    w <- colorClusters
    l$Treatment <- eval(parse(text = paste("pData(count_data)$",w,sep="")))
    l$Sample <- rownames(pData(count_data))
    g <- ggplot(l, aes(X1, X2, label=Sample, color=Treatment))+geom_point()
  }
}
