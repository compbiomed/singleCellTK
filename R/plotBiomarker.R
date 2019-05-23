#' Given a set of genes, return a ggplot of expression
#' values.
#'
#' @param visual Type of visualization (PCA, tSNE or UMAP). Default: "PCA"
#' @param inSCE Input SCtkExperiment object. Required
#' @param gene genelist to run the method on.
#' @param binary binary/continuous color for the expression.
#' @param shape shape parameter for the ggplot.
#' @param useAssay Indicate which assay to use. The default is "logcounts".
#' @param reducedDimName a name to store the results of the dimension reduction
#' coordinates obtained from this method. This is stored in the SingleCellExperiment
#' object in the reducedDims slot. Required.
#'
#' @return A Biomarker plot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotBiomarker(mouseBrainSubsetSCE, gene="C1qa", shape="level1class")
#'
plotBiomarker <- function(inSCE, gene, binary="Binary", visual="PCA",
                          shape="No Shape",
                          useAssay="counts", reducedDimName="PCA"){
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (visual == "PCA"){
    if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))) {
      inSCE <- getPCA(inSCE, useAssay = useAssay,
                      reducedDimName = reducedDimName)
    }
    variances <- NULL
    if (class(inSCE) == "SCtkExperiment"){
      variances <- pcaVariances(inSCE)
    }
  }
  if (visual == "tSNE"){
    if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))) {
      inSCE <- getTSNE(inSCE, useAssay = useAssay,
                       reducedDimName = reducedDimName)
    }
  }
  if (visual == "UMAP"){
    if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))) {
      inSCE <- getUMAP(inSCE, useAssay = useAssay,
                       reducedDimName = reducedDimName)
    }
  }
  axisDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                        reducedDimName))
  if (!is.null(colnames(axisDf))) {
    x <- colnames(axisDf)[1]
    y <- colnames(axisDf)[2]
  } else {
    colnames(axisDf) <- c("Comp1", "Comp2")
    x <- colnames(axisDf)[1]
    y <- colnames(axisDf)[2]
  }
  print(axisDf)
  if (length(gene) > 9) {
    gene <- gene[seq_len(9)]
  }
  for (i in seq_along(gene)){
    bioDf <- getBiomarker(inSCE = inSCE, gene = gene[i], binary = binary,
                          useAssay = useAssay)
    l <- axisDf
    if (!is.null(shape)){
      l$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
    }
    geneName <- colnames(bioDf)[2]
    colnames(bioDf)[2] <- "expression"
    l$Sample <- as.character(bioDf$sample)
    l$expression <- bioDf$expression
    c <- SummarizedExperiment::assay(inSCE, useAssay)[c(geneName), ]
    percent <- round(100 * sum(c > 0) / length(c), 2)
    if (visual == "PCA"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample",
                                                    color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"),
                                      values = c("Blue", "Grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))){
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample",
                                                      color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression),
                                                      max(l$expression)),
                                           low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(geneName, " - ", percent, "%", " cells",
                               sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      if (is.null(variances)){
        g <- g + ggplot2::labs(x = x, y = y)
      } else {
        g <- g + ggplot2::labs(
          x = paste0(x, " ", toString(round(variances[x, ] * 100, 2)), "%"),
          y = paste0(y, " ", toString(round(variances[y, ] * 100, 2)), "%"))
      }
    } else if (visual == "tSNE"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                    label = "Sample",
                                                    color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"),
                                      values = c("blue", "grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))) {
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                      label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                      label = "Sample",
                                                      color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression),
                                                      max(l$expression)),
                                           low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(geneName, " - ", percent, "%", " cells",
                               sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (visual == "UMAP"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                    label = "Sample",
                                                    color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"),
                                      values = c("blue", "grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))) {
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                      label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y,
                                                      label = "Sample",
                                                      color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression),
                                                      max(l$expression)),
                                           low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(geneName, " - ", percent, "%", " cells",
                               sep = "")) +
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
  return(grid::grid.draw(gridExtra::arrangeGrob(
    grobs = plist, ncol = ceiling(sqrt(length(gene))))))
}
