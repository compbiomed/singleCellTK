#' Plot t-SNE
#'
#' Use this function to plot t-SNE results
#'
#' @param count_data A SCE object
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param reducedDimName t-SNE dimension name. The default is TSNE
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#' @param runTSNE Run t-SNE if the reducedDimName does not exist. the Default is
#' FALSE.
#' @param use_assay Indicate which assay to use for t-SNE if you are running it.
#' Default is "counts"
#'
#' @return A t-SNE plot
#' @export
#' @examples
#' data("mouse_brain_subset_sce")
#' plotTSNE(mouse_brain_subset_sce, colorBy = "level1class",
#'          reducedDimName = "TSNE_counts")
#'
plotTSNE <- function(count_data, colorBy="No Color", shape="No Shape",
                     reducedDimName="TSNE", runTSNE=FALSE,
                     use_assay="logcounts"){
  if (is.null(SingleCellExperiment::reducedDim(count_data, reducedDimName))){
    if (runTSNE){
      count_data <- getTSNE(count_data, use_assay = use_assay,
                            reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName, " dimension not found. Run getTSNE() or set runTSNE to TRUE.")
    }
  }
  tsne_df <- data.frame(SingleCellExperiment::reducedDim(count_data, reducedDimName))
  if (ncol(tsne_df) > 2){
    warning("More than two t-SNE dimensions. Using the first two.")
  }
  xdim <- colnames(tsne_df)[1]
  ydim <- colnames(tsne_df)[2]
  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    tsne_df$color <- SingleCellExperiment::colData(count_data)[, colorBy]
  }
  if (!is.null(shape)){
    tsne_df$shape <- factor(SingleCellExperiment::colData(count_data)[, shape])
  }
  tsne_df$Sample <- colnames(count_data)
  g <- ggplot2::ggplot(tsne_df, ggplot2::aes_string(xdim, ydim, label = "Sample")) +
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
