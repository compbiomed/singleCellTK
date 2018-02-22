#' Plot PCA
#'
#' Use this function to plot PCA results
#'
#' @param count_data A SCE object
#' @param colorBy The variable to color clusters by
#' @param shape Shape of the points
#' @param pcX User choice for the first principal component
#' @param pcY User choice for the second principal component
#' @param reducedDimName PCA dimension name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#' @param runPCA Run PCA if the reducedDimName does not exist. the Default is
#' FALSE.
#' @param use_assay Indicate which assay to use for PCA. Default is "logcounts"

#'
#' @return A PCA plot
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' plotPCA(GSE60361_subset_sce, colorBy = "level1class",
#'         reducedDimName = "PCA_counts")
#'
plotPCA <- function(count_data, colorBy="No Color", shape="No Shape", pcX="PC1",
                    pcY="PC2", reducedDimName="PCA", runPCA=FALSE,
                    use_assay="logcounts"){
  if (is.null(SingleCellExperiment::reducedDim(count_data, reducedDimName))){
    if (runPCA){
      count_data <- getPCA(count_data, use_assay = use_assay,
                           reducedDimName = reducedDimName)
    } else {
      stop(reducedDimName, " dimension not found. Run getPCA() or set runPCA to TRUE.")
    }
  }
  pca_df <- data.frame(SingleCellExperiment::reducedDim(count_data, reducedDimName))
  if(!(pcX %in% colnames(pca_df))){
    stop("pcX dimension ", pcX, " is not in the reducedDim data")
  }
  if(!(pcY %in% colnames(pca_df))){
    stop("pcY dimension ", pcY, " is not in the reducedDim data")
  }

  if(class(count_data) == "SCtkExperiment"){
    if (all(c(pcX, pcY) %in% rownames(pca_variances(count_data)))){
      #use the variances in pca_variances
      variances <- pca_variances(count_data)
      pcXlab <- paste0(pcX, " ", toString(round(variances[pcX, ] * 100, 2)), "%")
      pcYlab <- paste0(pcY, " ", toString(round(variances[pcY, ] * 100, 2)), "%")
    } else {
      #do not use variances in the plot
      pcXlab <- pcX
      pcYlab <- pcY
    }
  } else {
    #do not use variances in the plot
    pcXlab <- pcX
    pcYlab <- pcY
  }

  if (colorBy == "No Color"){
    colorBy <- NULL
  }
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (!is.null(colorBy)){
    pca_df$color <- SingleCellExperiment::colData(count_data)[, colorBy]
  }
  if (!is.null(shape)){
    pca_df$shape <- factor(SingleCellExperiment::colData(count_data)[, shape])
  }
  pca_df$Sample <- colnames(count_data)
  g <- ggplot2::ggplot(pca_df, ggplot2::aes_string(pcX, pcY, label = "Sample")) +
    ggplot2::geom_point() +
    ggplot2::labs(x = pcXlab, y = pcYlab)
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
