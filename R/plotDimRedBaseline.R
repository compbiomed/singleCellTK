#' Plot results of reduced dimensions data.
#'
#' @param inSCE Input SCtkExperiment object with saved dimension reduction components
#'  or a variable with saved results. Required
#' @param colorBy color by a condition(any column of the annotation data).
#' @param conditionClass class of the annotation data used in colorBy. Options are
#'  NULL, "factor" or "numeric". If NULL, class will default to the original class.
#'  Default NULL.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment object. Required.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param comp1 label for x-axis
#' @param comp2 label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#' Default is second PCA component for PCA data and NULL otherwise.
#' @param background adds grid to plot when TRUE. Default TRUE.
#' @param size size of dots. Default 2.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#'
#' @return a ggplot of the reduced dimensions.
#' @export
#' @examples
#' plotDimRedBaseline(inSCE = mouseBrainSubsetSCE, colorBy = "No Color",
#'            shape = "No Shape", reducedDimName = "TSNE_counts",
#'            useAssay = "counts", comp1 = "tSNE1", comp2 = "tSNE2")
#'
plotDimRedBaseline <- function(inSCE,
                               colorBy,
                               shape,
                               reducedDimName,
                               conditionClass = NULL,
                               useAssay,
                               labelClusters = FALSE,
                               comp1,
                               comp2,
                               dim1 = NULL,
                               dim2 = NULL,
                               background = FALSE,
                               size = 2,
                               title = NULL,
                               titleSize = 15) {
    Df <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                      reducedDimName))
    if (ncol(Df) > 2){
        warning("More than two dimensions. Using the first two.")
    }
    if (!is.null(dim1) & !is.null(dim2)){
        if (!(dim1 %in% colnames(Df))){
            stop("X dimension ", dim1, " is not in the reducedDim data")
        }
        if (!(dim2 %in% colnames(Df))){
            stop("Y dimension ", dim2, " is not in the reducedDim data")
        }
        xdim <- dim1
        ydim <- dim2
    } else if (!is.null(comp1) & !is.null(comp2)){
        colnames(Df)[1] <- comp1
        colnames(Df)[2] <- comp2
        xdim <- colnames(Df)[1]
        ydim <- colnames(Df)[2]
    } else {
        xdim <- colnames(Df)[1]
        ydim <- colnames(Df)[2]
    }
    if (length(colorBy) == 1){
        if (colorBy == "No Color"){
            colorBy <- NULL
        }
    }
    if (shape == "No Shape"){
        shape <- NULL
    }
    if (!is.null(colorBy)){
        Df$color <- colorBy
    }
    if (!is.null(shape)){
        Df$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
    }
    Df$Sample <- colnames(inSCE)
    g <- ggplot2::ggplot(Df, ggplot2::aes_string(xdim, ydim,
                                                 label = "Sample")) +
        ggplot2::geom_point(size = size)
    if (!is.null(colorBy)){
        g <- g + ggplot2::aes_string(color = "color") +
            ggplot2::labs(color = colorBy)
    }
    if (!is.null(shape)){
        g <- g + ggplot2::aes_string(shape = "shape") +
            ggplot2::labs(shape = shape)
    }
    if (!background){
        g <- g + ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(color = "black"))
    }
    if(!is.null(title)){
        g <- g + ggplot2::ggtitle(label = title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = titleSize))
    }

    if (isTRUE(labelClusters)) {
        centroidList <- lapply(unique(colorBy), function(x) {
            df.sub <- Df[Df$color == x, ]
            median.1 <- stats::median(df.sub[, 1])
            median.2 <- stats::median(df.sub[, 2])
            cbind(median.1, median.2, as.character(x))
        })
        centroid <- do.call(rbind, centroidList)
        centroid <- data.frame(
            Dimension_1 = as.numeric(centroid[, 1]),
            Dimension_2 = as.numeric(centroid[, 2]),
            color = centroid[, 3],
            Sample = seq(1,length(unique(colorBy)))
        )

        colnames(centroid)[seq(2)] <- c(xdim, ydim)
        g <- g + ggplot2::geom_point(
            data = centroid,
            mapping = ggplot2::aes_string(
                x = xdim,
                y = ydim
            ),
            size = 0,
            alpha = 0
        ) +
            ggrepel::geom_text_repel(
                data = centroid,
                mapping = ggplot2::aes(label = color),
                show.legend = F,
                color = "black"
            )
    }

    return(g)
}


#' Plot results of reduced dimensions data of annotations in ColData.
#'
#' @param inSCE Input SCtkExperiment object with saved dimension reduction components
#'  or a variable with saved results. Required
#' @param colorBy color by a condition(any column of the annotation data).
#' @param conditionClass class of the annotation data used in colorBy. Options are
#'  NULL, "factor" or "numeric". If NULL, class will default to the original class.
#'  Default NULL.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment object. Required.
#' @param comp1 label for x-axis
#' @param comp2 label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#' @param background adds grid to plot when TRUE. Default TRUE.
#' @param size size of dots. Default 2.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#'
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotDimRedColData(inSCE = mouseBrainSubsetSCE, colorBy = "tissue",
#'            shape = "No Shape", conditionClass = "factor",
#'            reducedDimName = "TSNE_counts",
#'            comp1 = "tSNE1", comp2 = "tSNE2", labelClusters = TRUE)
plotDimRedColData <- function(inSCE,
                              colorBy = "No Color",
                              shape = "No Shape",
                              reducedDimName = NULL,
                              conditionClass = NULL,
                              comp1 = NULL,
                              comp2 = NULL,
                              dim1 = NULL,
                              dim2 = NULL,
                              background = TRUE,
                              size = 2,
                              title = NULL,
                              titleSize = 15,
                              labelClusters = TRUE) {

    colorPlot <- SingleCellExperiment::colData(inSCE)[, colorBy]
    if(!is.null(conditionClass)){
        colorPlot <- as(colorPlot, Class = conditionClass)
    }

    g <- plotDimRedBaseline(inSCE = inSCE,
                            colorBy = colorPlot,
                            shape = shape,
                            reducedDimName = reducedDimName,
                            conditionClass = conditionClass,
                            useAssay = useAssay,
                            comp1 = comp1,
                            comp2 = comp2,
                            dim1 = dim1,
                            dim2 = dim2,
                            background = background,
                            size = size,
                            title = title,
                            titleSize = titleSize,
                            labelClusters = labelClusters
                            )

    return(g)
}

#' Plot results of reduced dimensions data of counts stored in assay.
#'
#' @param inSCE Input SCtkExperiment object with saved dimension reduction components
#'  or a variable with saved results. Required
#' @param feature name of feature stored in assay of singleCellExperiment object.
#'  Plot will be colored based on feature value.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment object. Required.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param comp1 label for x-axis
#' @param comp2 label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#' Default is second PCA component for PCA data and NULL otherwise.
#' @param background adds grid to plot when TRUE. Default TRUE.
#' @param size size of dots. Default 2.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotDimRedAssay(inSCE = mouseBrainSubsetSCE, feature = "Sox2",
#'            shape = "No Shape", reducedDimName = "TSNE_counts",
#'            useAssay = "counts", comp1 = "tSNE1", comp2 = "tSNE2")

plotDimRedAssay <- function(inSCE,
                              feature,
                              shape = "No Shape",
                              reducedDimName = NULL,
                              useAssay = "logcounts",
                              comp1 = NULL,
                              comp2 = NULL,
                              dim1 = NULL,
                              dim2 = NULL,
                              background = TRUE,
                              size = 2,
                              title = NULL,
                              titleSize = 15) {

    mat <- getBiomarker(inSCE = inSCE,
                        useAssay = useAssay,
                        gene = feature,
                        binary = "Continuous")
    counts <- mat[,2]

    g <- plotDimRedBaseline(inSCE = inSCE,
                            conditionClass = "numeric",
                            colorBy = counts,
                            shape = shape,
                            reducedDimName = reducedDimName,
                            useAssay = useAssay,
                            comp1 = comp1,
                            comp2 = comp2,
                            dim1 = dim1,
                            dim2 = dim2,
                            background = background,
                            size = size,
                            title = title,
                            titleSize = titleSize
    )

    return(g)
}

