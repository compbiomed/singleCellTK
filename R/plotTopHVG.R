#' Plot highly variable genes
#'
#' @param inSCE Input \code{SingleCellExperiment} object containing the
#' computations.
#' @param method Select either \code{"vst"}, \code{"mean.var.plot"}, 
#' \code{"dispersion"} or \code{"modelGeneVar"}.
#' @param hvgNumber Specify the number of top genes to highlight in red. Default
#' \code{2000}.
#' @param labelsCount Specify the number of data points/genes to label. Should 
#' be less than \code{hvgNumber}. Default \code{20}.
#' @param featureDisplay A character string for the \code{rowData} variable name
#' to indicate what type of feature ID should be displayed. If set by 
#' \code{\link{setSCTKDisplayRow}}, will by default use it. If \code{NULL}, will
#' use \code{rownames(inSCE)}.
#' @return ggplot of HVG metrics
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
                       hvgNumber = 2000,
                       labelsCount = 20,
                       featureDisplay = metadata(inSCE)$featureDisplay
                       )
{
  method <- match.arg(method)

  hvgList <- getTopHVG(inSCE = inSCE,
                       method = method, 
                       hvgNumber = hvgNumber,
                       featureDisplay = NULL)

  if (method == "vst") {
    x <- log(rowData(inSCE)$seurat_variableFeatures_vst_mean)
    y <- rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized
    labeling <- "Standardized Variance"
  } else if (method == "mean.var.plot") {
    x <- rowData(inSCE)$seurat_variableFeatures_mvp_mean
    y <- rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled
    labeling <- "Dispersion"
  } else if (method == "dispersion") {
    x <- rowData(inSCE)$seurat_variableFeatures_dispersion_mean
    y <- rowData(inSCE)$seurat_variableFeatures_dispersion_dispersionScaled
    labeling <- "Dispersion"
  } else if (method == "modelGeneVar") {
    x <- rowData(inSCE)$scran_modelGeneVar_mean
    y <- rowData(inSCE)$scran_modelGeneVar_totalVariance
    labeling <- "Variance"
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
    ggplot2::labs(x = "Mean", y = labeling)
  
  vfplot <- .ggSCTKTheme(vfplot)
  return(vfplot)
}
