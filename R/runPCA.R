#' Run multiple PCA approach
#'   
#' Use this function to run Principle component Analysis using different approach and load plot
#'
#' @param plot.type Either 'Single Plot' or 'Paired Plot'
#' @param method regular PCA/randomized PCA/...
#' @param countm A SCE object
#' @param annotm A SCE object
#' @param featurem A SCE object
#' @param involving.variables User choice for the prinicipal components
#' @param additional.variables User choice for feature
#' @param colorClusters The variable to color clusters by
#' 
#' @return A reduced dimension object
#' @export runPCA
#'

runPCA = function(plot.type, method, countm,annotm,featurem, involving.variables, additional.variables, colorClusters){
  dr0col  = function(M) M[, colSums(abs(M)) != 0]
  scaRAW = FromMatrix(exprsArray = as.matrix(countm), cData = annotm, fData = featurem)
  if(plot.type == "Single Plot"){
    variables = c(involving.variables, additional.variables)
    variables = variables[1:2]
    if(method == "regular PCA"){projection = prcomp(dr0col(t(assay(scaRAW))),scale = TRUE)$x}
    else if(method == "randomized PCA"){projection = rpca(dr0col(t(assay(scaRAW))), retx=TRUE, k=2)$x}
    #else if(method = "robust PCA"){}
    pca = data.table(projection,  as.data.frame(colData(scaRAW)))
    text.tmp = paste("ggplot(pca) + geom_point(aes(x=",variables[1],",y=",variables[2],",color=as.factor(",colorClusters,")))",sep = "")
  }
  else {#plot.type == "Paired Plot"
    variables = c(involving.variables, additional.variables)
    variables = variables[1:7]
    if(method == "regular PCA"){projection = prcomp(dr0col(t(assay(scaRAW))),scale = TRUE)$x}
    else if(method == "randomized PCA"){projection = rpca(dr0col(t(assay(scaRAW))), retx=TRUE, k=3)$x}
    pca = data.table(projection,  as.data.frame(colData(scaRAW)))
    text.tmp = paste("ggpairs(pca, columns=variables,mapping=aes(color=",colorClusters,"), upper=list(continuous='blank'))",sep = "")
  }
  g = eval(parse(text = text.tmp))
  return(g)
}
