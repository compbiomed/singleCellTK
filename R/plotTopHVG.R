#' Plot highly variable genes
#'
#' @param inSCE Input \code{SingleCellExperiment} object containing the
#' computations.
#' @param method Select either "vst", "mean.var.plot", "dispersion" or
#' "modelGeneVar".
#' @param hvgList Character vector indicating the labels of highly variable
#' genes.
#' @param n Specify the number of top genes to highlight in red. If
#' \code{hvgList}
#'  parameter is not provided, this parameter can be used simply to specify
#'  the number of top genes to highlight in red.
#' @param labelsCount Specify the number of data points/genes to label.
#'  By default, all top genes will be labeled.
#' @return plot object
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- scranModelGeneVar(mouseBrainSubsetSCE, "logcounts")
#' plotTopHVG(mouseBrainSubsetSCE, method = "modelGeneVar",
#'            n = 1000, labelsCount = 0)
plotTopHVG <- function(inSCE,
                       method = c("vst", "mean.var.plot", "dispersion",
                                  "modelGeneVar"),
                       hvgList = NULL,
                       n = NULL,
                       labelsCount = NULL)
{
  method <- match.arg(method)
  if(is.null(n)){
    n = length(hvgList)
  }
  else{
    hvgList <- getTopHVG(
      inSCE = inSCE,
      method = method,
      n = n)
  }

  if(is.null(labelsCount)){
    labelsCount = n
  }

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
  vfplot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(x = subset(x, rownames(inSCE) %in% hvgList[seq(n)]),
                   y = subset(y, rownames(inSCE) %in% hvgList[seq(n)])),
                   colour = "red") +
    ggplot2::geom_label(ggplot2::aes(x = subset(x, rownames(inSCE) %in% hvgList[seq(labelsCount)]),
                   y = subset(y, rownames(inSCE) %in% hvgList[seq(labelsCount)]),
                   label = subset(rownames(inSCE),
                                  rownames(inSCE) %in% hvgList[seq(labelsCount)])),
               colour = "red",
               size = 2) +
    ggplot2::labs(x = "Mean", y = labeling)
  
  vfplot <- .ggSCTKTheme(vfplot)
  return(vfplot)
}
