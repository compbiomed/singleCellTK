#' Plot UMAP results either on already run results or run first and then plot.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components. Required
#' @param colorBy color by a condition(any column of the annotation data).
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param runUMAP If the dimension reduction components are already available
#' set this to FALSE, otherwise set to TRUE. Default is False.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#'
#' @return a UMAP plot of the reduced dimensions.
#' @export
#'
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' plotUMAP(sce, shape = "No Shape", reducedDimName = "UMAP",
#'          runUMAP = TRUE, useAssay = "counts")
#'
plotUMAP <- function(inSCE, colorBy = "No Color", shape = "No Shape",
                     reducedDimName = "UMAP", runUMAP = FALSE,
                     useAssay = "logcounts"){
  if(!(reducedDimName %in% names(SingleCellExperiment::reducedDims(inSCE)))){
    if (runUMAP){
      inSCE <- getUMAP(inSCE, useAssay = useAssay,
                       reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName,
           " dimension not found. Run getUMAP() or set runUMAP to TRUE.")
    }
  }
  UMAPDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                        reducedDimName))
  if (ncol(UMAPDf) > 2){
    warning("More than two UMAP dimensions. Using the first two.")
  }
  colnames(UMAPDf)[1] <- "UMAP1"
  colnames(UMAPDf)[2] <- "UMAP2"
  xdim <- colnames(UMAPDf)[1]
  ydim <- colnames(UMAPDf)[2]
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    UMAPDf$color <- SingleCellExperiment::colData(inSCE)[, colorBy]
  }
  if (!is.null(shape)){
    UMAPDf$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
  }
  UMAPDf$Sample <- colnames(inSCE)
  g <- ggplot2::ggplot(UMAPDf, ggplot2::aes_string(xdim, ydim,
                                                   label = "Sample")) +
    ggplot2::geom_point()
  if (!is.null(colorBy)){
    g <- g + ggplot2::aes_string(color = "color") +
      ggplot2::labs(color = colorBy)
  }
  if (!is.null(shape)){
    g <- g + ggplot2::aes_string(shape = "shape") +
      ggplot2::labs(shape = shape)
  }
  return(g)
}
