#' @title Plot Bubble plot
#' @description Plot a bubble plot with the color of the plot being the mean expression
#' and the size of the dot being the percent of cells in the cluster expressing the gene.
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param gene A string or vector of strings with each gene to aggregate.
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
#' plotBubble(inSCE=sce, useAssay="counts", gene=c("B2M", "MALAT1"), displayName="feature_name", 
#' clusters="type", title="cell type test", xlab="gene", ylab="cluster", colorLow="white", colorHigh="blue")
#' @export
plotBubble <- function(inSCE, useAssay="logcounts", gene, displayName=NULL, clusters="cluster", title="", xlab=NULL, ylab=NULL, colorLow="white", colorHigh="blue"){
  summaryMetrics <- runClusterSummaryMetrics(inSCE, useAssay=useAssay, gene=gene, displayName=displayName, clusters=clusters)
  .ggBubble(summaryMetrics, colorLow, colorHigh, title, xlab, ylab)
}

.ggBubble <- function(metrics, colorLow="white", colorHigh="blue", title="", xlab=NULL, ylab=NULL){
  gg <- ggplot2::ggplot(metrics, ggplot2::aes(x=Gene, y=cluster)) +
    ggplot2::geom_point(ggplot2::aes(colour=clusterAveExpr, size=clusterExprPerc)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_color_gradient(low=colorLow, high=colorHigh)
  .ggSCTKTheme(gg)
}