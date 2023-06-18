#' @title Plot Bubble plot
#' @description Plot a bubble plot with the color of the plot being the mean expression
#' and the size of the dot being the percent of cells in the cluster expressing the gene.
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param feature A string or vector of strings with each gene to aggregate.
#' @param displayName A string that is the name of the column used for genes.
#' @param clusters The name of a colData entry that can be used as groups.
#' @param title The title of the bubble plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param colorLow The color to be used for lowest value of mean expression
#' @param colorHigh The color to be used for highest value of mean expression
#' @return A ggplot of the bubble plot.
#' @examples
#' data("scExample")
#' plotBubble(inSCE=sce, useAssay="counts", feature=c("B2M", "MALAT1"), displayName="feature_name", 
#' clusters="type", title="cell type test", xlab="gene", ylab="cluster", colorLow="white", colorHigh="blue")
#' @export
plotBubble <- function(inSCE, useAssay="logcounts", feature, displayName=NULL, clusters="cluster", title="", xlab=NULL, ylab=NULL, colorLow="white", colorHigh="blue"){
  metrics <- runClusterSummaryMetrics(inSCE, useAssay=useAssay, feature=feature, displayName=displayName, clusters=clusters)
  .ggBubble(metrics$avgExpr, metrics$percExpr, colorLow, colorHigh, title)
}

.ggBubble <- function(avgExpr, percExpr, colorLow="white", colorHigh="blue", title=""){
  x <- NULL
  y <- NULL
  df <- data.frame(x=avgExpr, y=percExpr)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x=x.Gene, y=y.cluster)) +
    ggplot2::geom_point(ggplot2::aes(color=x.clusterAveExpr, size=y.clusterExprPerc)) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_color_gradient2(low=colorLow, high=colorHigh)
  .ggSCTKTheme(gg)
}



