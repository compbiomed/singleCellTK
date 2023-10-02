#' @title Plot Bubble plot
#' @description Plot a bubble plot with the color of the plot being the mean expression
#' and the size of the dot being the percent of cells in the cluster expressing the gene.
#'
#' @param inSCE The single cell experiment to use.
#' @param useAssay The assay to use.
#' @param featureNames A string or vector of strings with each gene to aggregate.
#' @param displayName A string that is the name of the column used for genes.
#' @param groupNames The name of a colData entry that can be used as groupNames.
#' @param title The title of the bubble plot
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param colorLow The color to be used for lowest value of mean expression
#' @param colorHigh The color to be used for highest value of mean expression
#' @param scale Option to scale the data
#' @return A ggplot of the bubble plot.
#' @importFrom rlang .data
#' @importFrom reshape2 melt
#' @examples
#' data("scExample")
#' plotBubble(inSCE=sce, useAssay="counts", featureNames=c("B2M", "MALAT1"), displayName="feature_name", 
#' groupNames="type", title="cell type test", xlab="gene", ylab="cluster", 
#' colorLow="white", colorHigh="blue")
#' @export
plotBubble <- function(inSCE, useAssay="logcounts", featureNames, displayName=NULL, groupNames="cluster", title="", xlab=NULL, ylab=NULL, colorLow="white", colorHigh="blue", scale = FALSE){
  metrics <- runClusterSummaryMetrics(inSCE, useAssay=useAssay, featureNames=featureNames, 
                                      displayName=displayName, groupNames=groupNames, scale = scale)
  .ggBubble(avgExpr = metrics$avgExpr, percExpr = metrics$percExpr, colorLow = colorLow, 
            colorHigh = colorHigh, title = title, xlab=xlab, ylab=ylab)
}

.ggBubble <- function(avgExpr, percExpr, groupNames=NULL, featureNames=NULL, colorLow="white", colorHigh="blue", title="", xlab="Features", ylab="Clusters"){
  if(is.null(featureNames)) {
    if(is.null(rownames(avgExpr))) {
      stop("'featureNames' must be supplied or the 'rownames' of the average expression matrix must be set.")
    }
    featureNames <- rownames(avgExpr)
  } else {
    ## Error checking for feature length
    if(length(featureNames) != nrow(avgExpr)) {
      stop("'featureNames' must be the same length as the number of rows in the average expression matrix.")
  }
  }
  if(is.null(groupNames)) {
    if(is.null(colnames(avgExpr))) {
      stop("'featureNames' must be supplied or the 'rownames' of the average expression matrix must be set.")
    }
    groupNames <- colnames(avgExpr)
  } else {
    ## Error checking for feature length
    if(length(groupNames) != ncol(avgExpr)) {
      stop("'featureNames' must be the same length as the number of rows in the average expression matrix.")
    }
  }
  
  avgExpr <- data.frame(avgExpr)
  avgExpr$featureNames <- featureNames
  df <- reshape2::melt(avgExpr, id="featureNames", value.name = "avgExpr", variable.name = "groupNames")
  
  percExpr <- data.frame(percExpr)
  percExpr$featureNames <- featureNames
  dfP <- reshape2::melt(percExpr, id="featureNames", value.name = "percExpr")
  
  df$percExpr = dfP$percExpr
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data[['featureNames']], y = .data[['groupNames']])) +
    ggplot2::geom_point(ggplot2::aes(color=.data[['avgExpr']], size=.data[['percExpr']])) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(ylab) + 
    ggplot2::scale_color_gradient2(low=colorLow, high=colorHigh)
  .ggSCTKTheme(gg)
}

