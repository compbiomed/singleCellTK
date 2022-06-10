#' Plot highly variable genes
#'
#' @param inSCE Input \code{SingleCellExperiment} object containing the
#' computations.
#' @param method Select either \code{"vst"}, \code{"mean.var.plot"}, 
#' \code{"dispersion"} or \code{"modelGeneVar"}.
#' @param useFeatureSubset A character string for the \code{rowData} variable 
#' name to store a logical index of selected features. Default \code{NULL}. See 
#' details.
#' @param hvgNumber Specify the number of top genes to highlight in red. Default
#' \code{NULL}. See details.
#' @param labelsCount Specify the number of data points/genes to label. Should 
#' be less than \code{hvgNumber}. Default \code{20}. See details.
#' @param featureDisplay A character string for the \code{rowData} variable name
#' to indicate what type of feature ID should be displayed. If set by 
#' \code{\link{setSCTKDisplayRow}}, will by default use it. If \code{NULL}, will
#' use \code{rownames(inSCE)}.
#' @return ggplot of HVG metrics and top HVG labels
#' @details When \code{hvgNumber = NULL} and \code{useFeature = NULL}, only plot
#' the mean VS variance/dispersion scatter plot. When only \code{hvgNumber} set,
#' label the top \code{hvgNumber} HVGs ranked by the metrics calculated by 
#' \code{method}. When \code{useFeatureSubset} set, label the features in 
#' the subset on the scatter plot created with \code{method} and ignore 
#' \code{hvgNumber}. 
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runModelGeneVar(mouseBrainSubsetSCE)
#' plotTopHVG(mouseBrainSubsetSCE, method = "modelGeneVar")
#' @seealso \code{\link{runFeatureSelection}}, \code{\link{runSeuratFindHVG}},
#' \code{\link{runModelGeneVar}}, \code{\link{getTopHVG}}
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
plotTopHVG <- function(inSCE,
                       method = c("vst", "mean.var.plot", "dispersion",
                                  "modelGeneVar"),
                       hvgNumber = NULL,
                       useFeatureSubset = NULL,
                       labelsCount = 20,
                       featureDisplay = metadata(inSCE)$featureDisplay
                       )
{
  method <- match.arg(method)
  metric <- .dfFromHVGMetric(inSCE, method)
  yLabelChoice <- list(vst = "Standardized Variance", 
                       mean.var.plot = "Dispersion", dispersion = "Dispersion", 
                       modelGeneVar = "Variance")
  x <- metric[,"mean"]
  y <- metric[,"v_plot"]
  if (method == "vst") x <- log(x)
  yAxisLabel <- yLabelChoice[[method]]
  hvgList <- character()
  if (!is.null(useFeatureSubset)) {
    hvgList <- .parseUseFeatureSubset(inSCE, useFeatureSubset, 
                                      returnType = "character")
    hvgNumber <- length(hvgList)
  } else if (!is.null(hvgNumber)) {
    hvgList <- getTopHVG(inSCE = inSCE, method = method, hvgNumber = hvgNumber,
                         featureDisplay = NULL)
  }
  if (is.null(hvgNumber) || hvgNumber == 0) {
    redIdx <- logical()
  } else {
    redIdx <- rownames(inSCE) %in% hvgList[seq(hvgNumber)]
  }
  if (is.null(hvgNumber) || labelsCount == 0) {
    labelIdx <- logical()
  } else {
    labelIdx <- rownames(inSCE) %in% hvgList[seq(labelsCount)]
  }
  if (!is.null(featureDisplay)) {
    labelTxt <- rowData(inSCE)[[featureDisplay]][labelIdx]
  } else {
    labelTxt <- rownames(inSCE)[labelIdx]
  }
  vfplot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(x = subset(x, redIdx),
                                     y = subset(y, redIdx)),
                        colour = "red") +
    ggplot2::geom_label(ggplot2::aes(x = subset(x, labelIdx),
                                     y = subset(y, labelIdx),
                                     label = labelTxt),
                        colour = "red",
                        size = 2) +
    ggplot2::labs(x = "Mean", y = yAxisLabel)
  vfplot <- .ggSCTKTheme(vfplot)
  return(vfplot)
}
