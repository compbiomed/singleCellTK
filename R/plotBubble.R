#' @title Plot Bubble plot
#' @description Plot a bubble plot with the color of the plot being the mean expression
#' and the size of the dot being the percent of cells in the cluster expressing the gene.
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param features A string or vector of strings with each gene to aggregate.
#' @param displayName A string that is the name of the column used for genes.
#' @param clusters The name of a colData entry that can be used as groups.
#' @param title The title of the bubble plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param colorLow The color to be used for lowest value of mean expression
#' @param colorHigh The color to be used for highest value of mean expression
#' @return A ggplot of the bubble plot.
#' @importFrom rlang .data
#' @examples
#' data("scExample")
#' plotBubble(inSCE=sce, useAssay="counts", features=c("B2M", "MALAT1"), displayName="feature_name", 
#' clusters="type", title="cell type test", xlab="gene", ylab="cluster", 
#' colorLow="white", colorHigh="blue")
#' @export
plotBubble <- function(inSCE, useAssay="logcounts", features, displayName=NULL, clusters="cluster", title="", xlab=NULL, ylab=NULL, colorLow="white", colorHigh="blue"){
  metrics <- runClusterSummaryMetrics(inSCE, useAssay=useAssay, features=features, 
                                      displayName=displayName, clusters=clusters)
  .ggBubble(metrics$features, metrics$avgExpr, metrics$percExpr, colorLow, colorHigh, title)
}

.ggBubble <- function(features, avgExpr, percExpr, colorLow="white", colorHigh="blue", title=""){
  clusters <- colnames(avgExpr)
  clusters <- rep.int(clusters, length(features))
  features = rep(features, each = length(avgExpr[1,]))
  df = data.frame(features = features, clusters = clusters)
  df$avgExpr = unlist(as.list(t(avgExpr)))
  df$percExpr = unlist(as.list(t(percExpr)))
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data[['features']], y = .data[['clusters']])) +
    ggplot2::geom_point(ggplot2::aes(color=.data[['avgExpr']], size=.data[['percExpr']])) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_color_gradient2(low=colorLow, high=colorHigh)
  .ggSCTKTheme(gg)
}



