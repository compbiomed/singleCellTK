#' Plot highly variable genes
#'
#' @param inSCE Input \code{SingleCellExperiment} object containing the computations.
#' @param method Select either "vst", "mean.var.plot", "dispersion" or "modelGeneVar".
#' @param hvgList Character vector indicating the labels of highly variable genes.
#'
#' @return plot object
#' @export
plotTopHVG <- function(inSCE,
                       method = "vst",
                       hvgList)
{
  if (method == "vst") {
    x <- rowData(inSCE)$seurat_variableFeatures_vst_mean
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
    ggplot2::geom_point(ggplot2::aes(x = subset(x, rownames(inSCE) %in% hvgList),
                   y = subset(y, rownames(inSCE) %in% hvgList)),
               colour = "red") +
    ggplot2::geom_label(ggplot2::aes(x = subset(x, rownames(inSCE) %in% hvgList),
                   y = subset(y, rownames(inSCE) %in% hvgList),
                   label = subset(rownames(inSCE),
                                  rownames(inSCE) %in% hvgList)),
               colour = "red",
               size = 2) +
    ggplot2::labs(x = "Mean", y = labeling)
  
  return(vfplot)
}