#' plotBiomarker
#'
#' Given a set of genes, return a ggplot of expression values.
#'
#' @param count_data A SingleCelltkExperiment object
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param visual Type of visualization (PCA or tSNE). Default: "PCA"
#' @param shape visualization shape
#' @param x x coordinate for PCA
#' @param y y coordinate for PCA
#'
#' @return A Biomarker plot
#' @export plotBiomarker
#'
plotBiomarker <- function(count_data, gene, binary="Binary", visual="PCA",
                          shape="No Shape", x="PC1", y="PC2"){
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (visual == "PCA"){
    if (is.null(reducedDim(count_data, "PCA"))) {
      count_data <- getPCA(count_data)
    }
    axis_df <- data.frame(reducedDim(count_data, "PCA"))
    variances <- pca_variances(count_data)$percentVar
  }
  if (visual == "tSNE"){
    if (is.null(reducedDim(count_data, "TSNE"))) {
      count_data <- getTSNE(count_data)
    }
    axis_df <- data.frame(reducedDim(count_data, "TSNE"))
  }
  if (length(gene) > 9) {
    gene <- gene[1:9]
  }
  for (i in 1:length(gene)){
    bio_df <- getBiomarker(count_data, gene[i], binary)
    l <- axis_df
    if (!is.null(shape)){
      l$shape <- factor(eval(parse(text = paste0("colData(count_data)$", shape))))
    }
    gene_name <- colnames(bio_df)[2]
    colnames(bio_df)[2] <- "expression"
    l$Sample <- as.character(bio_df$sample)
    l$expression <- bio_df$expression
    c <- assay(count_data, "counts")[c(gene_name), ]
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
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(x = paste(x, toString(round(variances[strtoi(strsplit(x, "PC")[[1]][-1])] * 100, 2)), "%"), y = paste(y, toString(round(variances[strtoi(strsplit(y, "PC")[[1]][-1])] * 100, 2)), "%"))
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
