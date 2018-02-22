#' plotBiomarker
#'
#' Given a set of genes, return a ggplot of expression values.
#'
#' @param count_data A SCtkExperiment object
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param visual Type of visualization (PCA or tSNE). Default: "PCA"
#' @param shape visualization shape
#' @param x x coordinate for PCA
#' @param y y coordinate for PCA
#' @param use_assay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName PCA dimension name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A Biomarker plot
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' plotBiomarker(GSE60361_subset_sce, gene="C1qa", shape="level1class")
#'
plotBiomarker <- function(count_data, gene, binary="Binary", visual="PCA",
                          shape="No Shape", x="PC1", y="PC2",
                          use_assay="counts", reducedDimName="PCA"){
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (visual == "PCA"){
    if (is.null(SingleCellExperiment::reducedDim(count_data, reducedDimName))) {
      count_data <- getPCA(count_data, use_assay = use_assay, reducedDimName = reducedDimName)
    }
    axis_df <- data.frame(SingleCellExperiment::reducedDim(count_data, reducedDimName))
    variances <- NULL
    if (class(count_data) == "SCtkExperiment"){
      variances <- pca_variances(count_data)
    }
  }
  if (visual == "tSNE"){
    if (is.null(SingleCellExperiment::reducedDim(count_data, reducedDimName))) {
      count_data <- getTSNE(count_data, use_assay = use_assay, reducedDimName = reducedDimName)
    }
    axis_df <- data.frame(SingleCellExperiment::reducedDim(count_data, reducedDimName))
  }
  if (length(gene) > 9) {
    gene <- gene[1:9]
  }
  for (i in 1:length(gene)){
    bio_df <- getBiomarker(count_data = count_data,
                           gene = gene[i],
                           binary = binary,
                           use_assay = use_assay)
    l <- axis_df
    if (!is.null(shape)){
      l$shape <- factor(SingleCellExperiment::colData(count_data)[, shape])
    }
    gene_name <- colnames(bio_df)[2]
    colnames(bio_df)[2] <- "expression"
    l$Sample <- as.character(bio_df$sample)
    l$expression <- bio_df$expression
    c <- SummarizedExperiment::assay(count_data, use_assay)[c(gene_name), ]
    percent <- round(100 * sum(c > 0) / length(c), 2)
    if (visual == "PCA"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample", color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"), values = c("Blue", "Grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))){
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample", color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression), max(l$expression)), low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(gene_name, " - ", percent, "%", " cells", sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      if(is.null(variances)){
        g <- g + ggplot2::labs(x = x, y = y)
      } else {
        g <- g + ggplot2::labs(x = paste0(x, " ", toString(round(variances[x,] * 100, 2)), "%"),
                               y = paste0(y, " ", toString(round(variances[y,] * 100, 2)), "%"))
      }
    } else if (visual == "tSNE"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes(X1, X2, label = Sample, color = expression)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"), values = c("blue", "grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))) {
          g <- ggplot2::ggplot(l, ggplot2::aes(X1, X2, label = Sample)) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes(X1, X2, label = Sample, color = expression)) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression), max(l$expression)), low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(gene_name, " - ", percent, "%", " cells", sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    if (!is.null(shape)){
      g <- g + ggplot2::aes_string(shape = "shape") +
        ggplot2::labs(shape = shape)
    }
    if (i == 1) {
      plist <- list(g)
    } else{
      plist <- cbind(plist, list(g))
    }
  }
  return(grid::grid.draw(gridExtra::arrangeGrob(grobs = plist, ncol = ceiling(sqrt(length(gene))))))
}
