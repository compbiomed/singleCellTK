#' @title Plot Bubble plot
#' @description Plot a bubble plot with the color of the plot being the mean expression
#' and the size of the dot being the percent of cells in the cluster expressing the gene.
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param gene A string or vector of strings with each gene to aggregate.
#' @param clusters The name of a colData entry that can be used as groups.
#' @param title The title of the bubble plot
#' @param colorLow The color to be used for lowest value of mean expression
#' @param colorHigh The color to be used for highest value of mean expression
#' @return A ggplot of the bubble plot.
#' @examples
#' plotBubble(inSCE=pbmc3k_2.7.1_sce, useAssay="logcounts", gene=c("IL7R", "CD3E"), 
#' title="cell type test", colorLow="white", colorHigh="blue", clusters="cluster")
#' @export
plotBubble <- function(inSCE, useAssay="logcounts", gene, clusters="cluster", title="", colorLow="white", colorHigh="blue"){
  summaryMetrics <- runClusterSummaryMetrics(inSCE, useAssay=useAssay, gene=gene, clusters=clusters)
  .ggBubble(summaryMetrics, colorLow, colorHigh, title)
}

.ggBubble <- function(metrics, colorLow="white", colorHigh="blue", title=""){
  gg <- ggplot2::ggplot(metrics, ggplot2::aes(x=Gene, y=cluster)) +
    ggplot2::geom_point(ggplot2::aes(color=clusterAveExpr, size=clusterPercExpr)) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_color_gradient2(low=colorLow, high=colorHigh)
  .ggSCTKTheme(gg)
}