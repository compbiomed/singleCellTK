#' @title Plot results of reduced dimensions data.
#' @description Plot results of reduced dimensions data and colors the plots by
#'  the input vector.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param colorBy If provided, colors dots in the scatterplot based on value.
#' @param groupBy If provided, facet wrap the scatterplot based on value.
#' @param conditionClass class of the annotation data used in colorBy. Options
#'  are NULL, "factor" or "numeric". If NULL, class will default to the original
#'  class. Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dim1 1st dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 2.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param colorScale Vector. Needs to be same length as the
#'  number of unique levels of `colorBy`. Will be used only if
#'  conditionClass = "factor" or "character". Default NULL.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'white'. Will be used only if conditionClass = "numeric".
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#'  Default 'gray'. Will be used only if conditionClass = "numeric".
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'. Will be used only if conditionClass = "numeric".
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default 10.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimensions.
.ggScatter <- function(inSCE,
                       sample = NULL,
                       colorBy = NULL,
                       groupBy = NULL,
                       shape = NULL,
                       reducedDimName,
                       conditionClass = NULL,
                       labelClusters = FALSE,
                       xlab = NULL,
                       ylab = NULL,
                       axisSize = 10,
                       axisLabelSize = 10,
                       dim1 = NULL,
                       dim2 = NULL,
                       bin = NULL,
                       binLabel = NULL,
                       dotSize = 2,
                       transparency = 1,
                       colorScale = NULL,
                       colorLow = "white",
                       colorMid = "gray",
                       colorHigh = "blue",
                       defaultTheme = TRUE,
                       title = NULL,
                       titleSize = 15,
                       legendTitle = NULL,
                       legendTitleSize = 12,
                       legendSize = 10,
                       combinePlot = "none",
                       plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop(
        "'sample' must be the same length as the number",
        " of columns in 'inSCE'"
      )
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x){
    sceSampleInd <- which(sample == x)
    inSCESub <- inSCE[, sceSampleInd]
    colorBySub <- colorBy[sceSampleInd]

    dataframe <- data.frame(SingleCellExperiment::reducedDim(
      inSCESub,
      reducedDimName
    ))
    # If dim1 and dim2 are specified
    if (!is.null(dim1) & !is.null(dim2)) {
      if (is.character(dim1) & is.character(dim2)) {
        if (!(dim1 %in% colnames(dataframe))) {
          stop("X dimension ", dim1, " is not in the reducedDim data")
        }
        if (!(dim2 %in% colnames(dataframe))) {
          stop("Y dimension ", dim2, " is not in the reducedDim data")
        }
        dataframe <- dataframe[, c(dim1, dim2)]
      } else if (is.numeric(dim1) && is.numeric(dim2)) {
        dataframe <- dataframe[, c(dim1, dim2)]
      }
    } else if (ncol(dataframe) > 2) {
      warning("More than two dimensions supplied in reducedDims.
              Using the first two.")
    }

    # If xlab and ylab are specified
    if (!is.null(xlab) & !is.null(ylab)) {
      colnames(dataframe) <- c(xlab, ylab)
      # If reduced dimension matrix didnt have colnames
    } else {
      colnames(dataframe) <- c(paste0(reducedDimName, "_1"),
                               paste0(reducedDimName, "_2"))
    }

    xdim <- colnames(dataframe)[1]
    ydim <- colnames(dataframe)[2]

    if (!is.null(conditionClass) & !is.null(colorBySub)) {
      if (conditionClass %in% c("character", "factor", "numeric")) {
        if (conditionClass == "character") {
          colorBySub <- as.character(colorBySub)
        } else if (conditionClass == "factor") {
          colorBySub <- as.factor(colorBySub)
        } else if (conditionClass == "numeric") {
          colorBySub <- as.numeric(colorBySub)
        }
      }
    }

    if (!is.null(bin) & !is.null(colorBySub)) {
      colorBySub <- .binSCTK(
        value = colorBySub,
        bin = bin,
        binLabel = binLabel
      )
    }

    if (!is.null(colorBySub)) {
      dataframe$color <- colorBySub
    }

    if (!is.null(groupBy)){
      dataframe$groups <- factor(SingleCellExperiment::colData(inSCE)@listData[[groupBy]])
    }
    if (!is.null(shape)) {
      dataframe$shape <- factor(SingleCellExperiment::colData(inSCESub)[, shape])
    }
    dataframe$Sample <- colnames(inSCESub)
    g <- ggplot2::ggplot(dataframe, ggplot2::aes_string(xdim, ydim,
                                                        label = "Sample")) + ggplot2::geom_point(size = dotSize,
                                                                                                 alpha = transparency)
    if (!is.null(colorBySub)) {
      g <- g + ggplot2::aes_string(color = "color")
    }
    if (inherits(colorBySub, "numeric")){
      g <- g + ggplot2::scale_color_gradient2(
        low = colorLow,
        mid = colorMid,
        high = colorHigh,
        aesthetics = "colour",
        midpoint = mean(colorBySub))
    }else if (inherits(colorBySub, "character") | inherits(colorBySub, "factor")){
      g <- g +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))
      if(all(!is.null(colorScale))){
        g <- g+ ggplot2::scale_color_manual(values=c(colorScale))
      }
    }
    if (!is.null(shape)) {
      g <- g + ggplot2::aes_string(shape = "shape") +
        ggplot2::labs(shape = shape)
    }
    if (defaultTheme == TRUE) {
      g <- .ggSCTKTheme(g)
    }
    g <- g + ggplot2::theme(axis.title =
                              ggplot2::element_text(size = axisLabelSize),
                            axis.text =
                              ggplot2::element_text(size = axisSize))
    if (!is.null(title)) {
      if (length(samples) > 1) {
        title <- paste(title, x, sep = "_")
      }
      g <- g + ggplot2::ggtitle(label = title) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = titleSize
        ))
    }
    if (!is.null(legendTitle)) {
      g <- g + ggplot2::labs(color = legendTitle) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=legendTitleSize),
                       legend.text=ggplot2::element_text(size=legendSize))
    } else {
      g <- g + ggplot2::labs(color = "") +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legendSize))
    }

    if (!is.null(groupBy)){
      g <- g + ggplot2::facet_wrap(~groups)
    }


    if (isTRUE(labelClusters) && class(colorBySub) %in% c("character", "factor")) {
      centroidList <- lapply(unique(colorBySub), function(x) {
        dataframe.sub <- dataframe[dataframe$color == x, ]
        median.1 <- stats::median(dataframe.sub[, 1])
        median.2 <- stats::median(dataframe.sub[, 2])
        cbind(median.1, median.2, as.character(x))
      })
      centroid <- do.call(rbind, centroidList)
      centroid <- data.frame(
        Dimension_1 = as.numeric(centroid[, 1]),
        Dimension_2 = as.numeric(centroid[, 2]),
        color = centroid[, 3],
        Sample = rep(1, length(unique(colorBySub)))
      )

      if (!is.null(shape)) {
        centroid$shape <- dataframe$shape[1]
      }

      if (!is.null(groupBy)){
        g <- g + ggplot2::facet_wrap(~groups)
      }

      colnames(centroid)[seq_len(2)] <- c(xdim, ydim)
      g <- g + ggplot2::geom_point(
        data = centroid,
        mapping = ggplot2::aes_string(x = xdim, y = ydim),
        size = 0,
        alpha = 0
      ) +
        ggrepel::geom_text_repel(
          data = centroid,
          mapping = ggplot2::aes_string(label = "color"),
          show.legend = FALSE,
          color = "black"
        )
    }
    return(g)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- list(Sample = plotlist)
  }

  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }else if(combinePlot == "none" && length(plotlist) == 1){
    plotlist <- plotlist[[1]]
  }

  return(plotlist)
}
#' @title Dimension reduction plot tool for colData
#' @description Plot results of reduced dimensions data and
#'  colors by annotation data stored in the colData slot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param colorBy Color by a condition(any column of the annotation data).
#'  Required.
#' @param groupBy Group by a condition(any column of the annotation data).
#'  Default NULL.
#' @param conditionClass Class of the annotation data used in colorBy.
#'  Options are NULL, "factor" or "numeric". If NULL, class will default to the
#'  original class. Default NULL.
#' @param shape Add shapes to each condition.
#' @param reducedDimName Saved dimension reduction matrix name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param dim1 1st dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 2.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param colorScale Vector. Needs to be same length as the
#'  number of unique levels of colorBy. Will be used only if
#'  conditionClass = "factor" or "character". Default NULL.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'white'.
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#'  Default 'gray'.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default 10.
#'  Default FALSE.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimension plot of coldata.
#' @examples
#' plotSCEDimReduceColData(
#'   inSCE = mouseBrainSubsetSCE, colorBy = "tissue",
#'   shape = NULL, conditionClass = "factor",
#'   reducedDimName = "TSNE_counts",
#'   xlab = "tSNE1", ylab = "tSNE2", labelClusters = TRUE
#' )
#'
#' plotSCEDimReduceColData(
#'   inSCE = mouseBrainSubsetSCE, colorBy = "age",
#'   shape = NULL, conditionClass = "numeric",
#'   reducedDimName = "TSNE_counts", bin = c(-Inf, 20, 25, +Inf),
#'   xlab = "tSNE1", ylab = "tSNE2", labelClusters = FALSE
#' )
#' @export
plotSCEDimReduceColData <- function(inSCE,
                                    sample = NULL,
                                    colorBy,
                                    groupBy = NULL,
                                    conditionClass = NULL,
                                    shape = NULL,
                                    reducedDimName = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    axisSize = 10,
                                    axisLabelSize = 10,
                                    dim1 = NULL,
                                    dim2 = NULL,
                                    bin = NULL,
                                    binLabel = NULL,
                                    dotSize = 2,
                                    transparency = 1,
                                    colorScale = NULL,
                                    colorLow = "white",
                                    colorMid = "gray",
                                    colorHigh = "blue",
                                    defaultTheme = TRUE,
                                    title = NULL,
                                    titleSize = 15,
                                    labelClusters = TRUE,
                                    legendTitle = NULL,
                                    legendTitleSize = 12,
                                    legendSize = 10,
                                    combinePlot = "none",
                                    plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  colorPlot <- SingleCellExperiment::colData(inSCE)[, colorBy]

  g <- .ggScatter(
    inSCE = inSCE,
    sample = sample,
    colorBy = colorPlot,
    groupBy = groupBy,
    conditionClass = conditionClass,
    shape = shape,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    dim1 = dim1,
    dim2 = dim2,
    axisSize = axisSize,
    axisLabelSize = axisLabelSize,
    bin = bin,
    binLabel = binLabel,
    dotSize = dotSize,
    transparency = transparency,
    colorScale = colorScale,
    colorLow = colorLow,
    colorMid = colorMid,
    colorHigh = colorHigh,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize,
    labelClusters = labelClusters,
    legendTitle = legendTitle,
    legendTitleSize = legendTitleSize,
    legendSize = legendSize,
    combinePlot = combinePlot,
    plotLabels = plotLabels
  )
  return(g)
}


#' @title Dimension reduction plot tool for assay data
#' @description Plot results of reduced dimensions data and
#'  colors by feature data stored in the assays slot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param feature name of feature stored in assay of singleCellExperiment
#'  object. Plot will be colored based on feature value.
#' @param shape add shapes to each condition. Default NULL.
#' @param reducedDimName saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dim1 1st dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#'  Default is second PCA component for PCA data and NULL otherwise.
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 2.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'white'.
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#'  Default 'gray'.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default 10.
#' @param groupBy Facet wrap the scatterplot based on value.
#' Default \code{NULL}.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimension plot of feature data.
#' @examples
#' plotSCEDimReduceFeatures(
#'   inSCE = mouseBrainSubsetSCE, feature = "Apoe",
#'   shape = NULL, reducedDimName = "TSNE_counts",
#'   useAssay = "counts", xlab = "tSNE1", ylab = "tSNE2"
#' )
#' @export
plotSCEDimReduceFeatures <- function(inSCE,
                                     sample = NULL,
                                     feature,
                                     shape = NULL,
                                     reducedDimName,
                                     useAssay = "logcounts",
                                     xlab = NULL,
                                     ylab = NULL,
                                     axisSize = 10,
                                     axisLabelSize = 10,
                                     dim1 = NULL,
                                     dim2 = NULL,
                                     bin = NULL,
                                     binLabel = NULL,
                                     dotSize = 2,
                                     transparency = 1,
                                     colorLow = "white",
                                     colorMid = "gray",
                                     colorHigh = "blue",
                                     defaultTheme = TRUE,
                                     title = NULL,
                                     titleSize = 15,
                                     legendTitle = NULL,
                                     legendSize = 10,
                                     legendTitleSize = 12,
                                     groupBy = NULL,
                                     combinePlot = "none",
                                     plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]

  g <- .ggScatter(
    inSCE = inSCE,
    sample = sample,
    conditionClass = "numeric",
    colorBy = counts,
    shape = shape,
    transparency = 1,
    colorLow = colorLow,
    colorMid = colorMid,
    colorHigh = colorHigh,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    axisLabelSize = axisLabelSize,
    dim1 = dim1,
    dim2 = dim2,
    bin = bin,
    binLabel = binLabel,
    defaultTheme = defaultTheme,
    dotSize = dotSize,
    title = title,
    titleSize = titleSize,
    legendTitle = legendTitle,
    legendTitleSize = legendTitleSize,
    legendSize = legendSize,
    groupBy = groupBy,
    combinePlot = combinePlot,
    plotLabels = plotLabels
  )

  return(g)
}

#' @title Dimension reduction plot tool for all types of data
#' @description Plot results of reduced dimensions data of counts stored in any
#' slot in the SingleCellExperiment object.
#' @param inSCE Input SingleCellExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata", "reducedDims". Default NULL.
#' @param annotation Desired vector within the slot used for plotting. Default NULL.
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object. Will be used only if "assays" slot is chosen. Default NULL.
#' @param groupBy Group by a condition(any column of the annotation data).
#'  Default NULL.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param conditionClass class of the annotation data used in colorBy. Options
#'  are NULL, "factor" or "numeric". If NULL, class will default to the original
#'  class. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dim1 1st dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' Default is second PCA component for PCA data and NULL otherwise.
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 2.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'white'.
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#'  Default 'gray'.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default 10.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimensions.
#' @examples
#' plotSCEScatter(
#'   inSCE = mouseBrainSubsetSCE, legendTitle = NULL,
#'   slot = "assays", annotation = "counts", feature = "Apoe",
#'   reducedDimName = "TSNE_counts", labelClusters = FALSE
#' )
#' @export
#' @import SingleCellExperiment
plotSCEScatter <- function(inSCE,
                           slot = NULL,
                           sample = NULL,
                           annotation,
                           feature = NULL,
                           groupBy = NULL,
                           shape = NULL,
                           reducedDimName = NULL,
                           conditionClass = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           axisSize = 10,
                           axisLabelSize = 10,
                           dim1 = NULL,
                           dim2 = NULL,
                           bin = NULL,
                           binLabel = NULL,
                           dotSize = 2,
                           transparency = 1,
                           colorLow = "white",
                           colorMid = "gray",
                           colorHigh = "blue",
                           defaultTheme = TRUE,
                           title = NULL,
                           titleSize = 15,
                           labelClusters = TRUE,
                           legendTitle = NULL,
                           legendTitleSize = 12,
                           legendSize = 10,
                           combinePlot = "none",
                           plotLabels = NULL){
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!is.null(slot)){
    if (slot == "reducedDims"){
      annotation_clm <- substr(annotation, stringr::str_length(annotation), stringr::str_length(annotation))
      annotation <- substr(annotation, 1, stringr::str_length(annotation) - 2)
    }else if (!slot %in% methods::slotNames(inSCE)) {
      stop("'slot' must be a slot within the SingleCellExperiment object.
             Please run 'methods::slotNames' if you are unsure the
	     specified slot exists.")
    }

    sceSubset <- do.call(slot, args = list(inSCE))

    if (!annotation %in% names(sceSubset)) {
      stop("'annotation' must be an annotation stored within the specified
             slot of the SingleCellExperiment object.")
    }
    annotation.ix <- match(annotation, c(names(sceSubset)))
  }

  if (is.null(slot)){
    colorPlot <- NULL
  }else if (slot == "assays" && !is.null(feature)) {
    counts <- sceSubset[[annotation.ix]]
    if (feature %in% rownames(counts)) {
      colorPlot <- counts[feature, ]
    }
  } else if (slot == "colData") {
    colorPlot <- sceSubset[, annotation.ix]
  } else if (slot == "metadata") {
    colorPlot <- sceSubset[[annotation.ix]]
  } else if (slot == "reducedDims") {
    colorPlot <- sceSubset[[annotation.ix]][, as.numeric(annotation_clm)]
  }

  g <- .ggScatter(
    inSCE = inSCE,
    sample = sample,
    colorBy = colorPlot,
    groupBy = groupBy,
    conditionClass = conditionClass,
    shape = shape,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    axisLabelSize = axisLabelSize,
    dim1 = dim1,
    dim2 = dim2,
    bin = bin,
    binLabel = binLabel,
    dotSize = dotSize,
    transparency = transparency,
    colorLow = colorLow,
    colorMid = colorMid,
    colorHigh = colorHigh,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize,
    labelClusters = labelClusters,
    legendTitle = legendTitle,
    legendTitleSize = legendTitleSize,
    legendSize = legendSize,
    combinePlot = combinePlot
  )
  return(g)
}

#' @title Violin plot plotting tool.
#' @description Visualizes specified values via a violin plot.
#' @param y Numeric values to be plotted on y-axis.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @return a ggplot of the reduced dimensions.
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
.ggViolin <- function(y,
                      groupBy = NULL,
                      violin = TRUE,
                      boxplot = TRUE,
                      dots = TRUE,
                      xlab = NULL,
                      ylab = NULL,
                      axisSize = 10,
                      axisLabelSize = 10,
                      dotSize = 1,
                      transparency = 1,
                      defaultTheme = TRUE,
                      gridLine = FALSE,
                      summary = NULL,
                      title = NULL,
                      titleSize = 15) {
  if (is.null(groupBy)) {
    groupBy <- rep("Sample", length(y))
  }
  df <- data.frame(groupBy = groupBy, y = y)

  p <- ggplot2::ggplot(df) +
    ggplot2::aes_string(
      x = "groupBy",
      y = "y"
    )
  if (dots == TRUE) {
    p <- p + ggplot2::geom_jitter(
      color = "blue",
      width = 0.2,
      height = 0,
      size = dotSize,
      alpha = transparency
    )
  }
  if (boxplot == TRUE) {
    p <- p + ggplot2::geom_boxplot(width = 0.5,
                                   alpha = 0)
  }
  if (violin == TRUE) {
    p <- p + ggplot2::geom_violin(trim = TRUE,
                                  scale = "width",
                                  size = 1,
                                  fill = "grey",
                                  alpha = 0.75)
  }
  if (defaultTheme == TRUE) {
    p <- .ggSCTKTheme(p)
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(label = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = titleSize
      ))
  }

  ###
  p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = axisSize))
  ###

  if(length(unique(df$groupBy)) > 1){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                hjust = 1,
                                                                size = axisSize))
  }else{
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank())
  }

  if (gridLine == TRUE){
    p <- p + ggplot2::theme(panel.grid.major.y = ggplot2::element_line("grey"))
  }
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize))
  }
  if (!is.null(ylab)) {
    p <- p + ggplot2::ylab(ylab) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize))
  }
  if (!is.null(summary)){
    if(summary == "mean"){
      summ <- df %>% dplyr::group_by(groupBy) %>% dplyr::summarize(value = base::mean(y))
      fun <- base::mean
    }else if(summary == "median"){
      summ <- df %>% dplyr::group_by(groupBy) %>% dplyr::summarize(value = stats::median(y))
      fun <- stats::median
    }else{
      stop("`summary`` must be either `mean` or `median`.")
    }
    summ$statY <-  max(df$y) + (max(df$y) - min(df$y)) * 0.1
    summary <- paste(toupper(substr(summary, 1, 1)),
                     substr(summary, 2, nchar(summary)), sep="")
    summ$label <- paste0(summary,": ", round(summ$value, 5))

    p <- p + ggplot2::geom_text(data = summ,
                                ggplot2::aes_string(x = "groupBy",
                                                    y = "statY",
                                                    label = "label"),
                                size = 5)
    p <- p + ggplot2::stat_summary(fun = fun, fun.min = fun,
                                   fun.max = fun,
                                   geom = "crossbar",
                                   color = "red",
                                   linetype = "dashed")
  }

  return(p)
}


#' @title Violin plot of colData.
#' @description Visualizes values stored in the colData slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param coldata colData value that will be plotted.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot of coldata.
#' @examples
#' plotSCEViolinColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEViolinColData <- function(inSCE,
                                 sample = NULL,
                                 coldata,
                                 groupBy = NULL,
                                 violin = TRUE,
                                 boxplot = TRUE,
                                 dots = TRUE,
                                 xlab = NULL,
                                 ylab = NULL,
                                 axisSize = 10,
                                 axisLabelSize = 10,
                                 dotSize = 1,
                                 transparency = 1,
                                 defaultTheme = TRUE,
                                 gridLine = FALSE,
                                 summary = NULL,
                                 title = NULL,
                                 titleSize = NULL,
                                 combinePlot = "none",
                                 plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!is.null(coldata)) {
    if (!coldata %in% names(SummarizedExperiment::colData(inSCE))) {
      stop("'", paste(coldata), "' is not found in ColData.")
    }
    coldata <- SummarizedExperiment::colData(inSCE)[, coldata]
  } else {
    stop("You must define the desired colData to plot.")
  }

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(coldata)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number",
           " of columns in 'inSCE'")
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    coldataSub <- coldata[sampleInd]
    if(!is.null(groupBy)){
      groupbySub <- groupBy[sampleInd]
    }else{
      groupbySub <- NULL
    }

    if(!is.null(title) && length(samples) > 1){
      title = paste(title, x, sep = "_")
    }

    p <- .ggViolin(
      y = coldataSub,
      groupBy = groupbySub,
      violin = violin,
      boxplot = boxplot,
      dots = dots,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      dotSize = dotSize,
      transparency = transparency,
      defaultTheme = defaultTheme,
      gridLine = gridLine,
      summary = summary,
      title = title,
      titleSize = titleSize
    )
    return(p)
  })

  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }else if(combinePlot == "none" && length(plotlist) == 1){
    plotlist <- plotlist[[1]]
  }

  return(plotlist)
}


#' @title Violin plot of assay data.
#' @description Visualizes values stored in the assay slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param feature Name of feature stored in assay of SingleCellExperiment
#'  object.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot of assay data.
#' @examples
#' plotSCEViolinAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEViolinAssayData <- function(inSCE,
                                   sample = NULL,
                                   useAssay = "counts",
                                   feature,
                                   groupBy = NULL,
                                   violin = TRUE,
                                   boxplot = TRUE,
                                   dots = TRUE,
                                   xlab = NULL,
                                   ylab = NULL,
                                   axisSize = 10,
                                   axisLabelSize = 10,
                                   dotSize = 1,
                                   transparency = 1,
                                   defaultTheme = TRUE,
                                   gridLine = FALSE,
                                   summary = NULL,
                                   title = NULL,
                                   titleSize = NULL,
                                   combinePlot = "none",
                                   plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]
  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(counts)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }
  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number",
           " of columns in 'inSCE'")
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    countSub <- counts[sampleInd]
    if(!is.null(groupBy)){
      groupbySub <- groupBy[sampleInd]
    }else{
      groupbySub <- NULL
    }

    p <- .ggViolin(
      y = countSub,
      groupBy = groupbySub,
      violin = violin,
      boxplot = boxplot,
      dots = dots,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      dotSize = dotSize,
      transparency = transparency,
      defaultTheme = defaultTheme,
      gridLine = gridLine,
      summary = summary,
      title = title,
      titleSize = titleSize
    )
    return(p)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    if(combinePlot == "sample"){
      plotlist <- c(list(Sample = plotlist))
    }
  } else {
    plotlist <- plotlist[[1]]
  }
  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }else if(combinePlot == "none" && length(plotlist) == 1){
    plotlist <- plotlist[[1]]
  }

  return(plotlist)
}

#' @title Violin plot of any data stored in the SingleCellExperiment object.
#' @description Visualizes values stored in any slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata"
#' @param annotation Desired vector within the slot used for plotting.
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object.
#'  Will be used only if "assays" slot is chosen. Default NULL.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot.
#' @examples
#' plotSCEViolin(
#'   inSCE = mouseBrainSubsetSCE, slot = "assays",
#'   annotation = "counts", feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEViolin <- function(inSCE,
                          sample = NULL,
                          slot,
                          annotation,
                          feature,
                          groupBy = NULL,
                          violin = TRUE,
                          boxplot = TRUE,
                          dots = TRUE,
                          xlab = NULL,
                          ylab = NULL,
                          axisSize = 10,
                          axisLabelSize = 10,
                          dotSize = 1,
                          transparency = 1,
                          defaultTheme = TRUE,
                          gridLine = FALSE,
                          summary = NULL,
                          title = NULL,
                          titleSize = NULL,
                          combinePlot = "none",
                          plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!slot %in% methods::slotNames(inSCE)) {
    stop("'slot' must be a slot within the SingleCellExperiment object.
             Please run 'methods::slotNames' if you are unsure the
	 specified slot exists.")
  }

  sceSubset <- do.call(slot, args = list(inSCE))

  if (!annotation %in% names(sceSubset)) {
    stop("'annotation' must be an annotation stored within the specified
             slot of the SingleCellExperiment object.")
  }

  annotation.ix <- match(annotation, names(sceSubset))

  if (slot == "assays" && !is.null(feature)) {
    counts <- sceSubset[[annotation.ix]]
    if (feature %in% rownames(counts)) {
      counts <- counts[feature, ]
    }
  } else if (slot == "colData") {
    counts <- sceSubset[, annotation.ix]
  } else if (slot == "metadata") {
    counts <- sceSubset[[annotation.ix]]
  }

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(counts)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number",
           " of columns in 'inSCE'")
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }
  samples <- unique(sample)
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    countSub <- counts[sampleInd]
    if(!is.null(groupBy)){
      groupbySub <- groupBy[sampleInd]
    }else{
      groupbySub <- NULL
    }
    p <- .ggViolin(
      y = countSub,
      groupBy = groupbySub,
      violin = violin,
      boxplot = boxplot,
      dots = dots,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      dotSize = dotSize,
      transparency = transparency,
      defaultTheme = defaultTheme,
      gridLine = gridLine,
      summary = summary,
      title = title,
      titleSize = titleSize
    )

    return(p)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    if(!is.null(combinePlot)){
      if(combinePlot == "sample"){
        plotlist <- c(list(Sample = plotlist))
      }
    }
  } else {
    plotlist <- plotlist[[1]]
    # plotlist <- unlist(plotlist, recursive=F)
  }

  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample") &&
     length(unique(samples)) > 1){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }else if(combinePlot == "none" && length(plotlist) == 1){
    plotlist <- plotlist[[1]]
  }
  return(plotlist)
}

#' @title Density plot plotting tool.
#' @description Visualizes values stored in the specified slot of a
#'  SingleCellExperiment object via a density plot.
#' @param value Numeric value that will be plotted via density plot.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param cutoff Numeric value. The plot will be annotated with a vertical line
#'  if set. Default NULL.
#' @return density plot, in .ggplot.
.ggDensity <- function(value,
                       groupBy = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       axisSize = 10,
                       axisLabelSize = 10,
                       defaultTheme = TRUE,
                       title = NULL,
                       titleSize = 18,
                       cutoff = NULL) {
  if (is.null(groupBy)) {
    groupBy <- rep("Sample", length(value))
  }
  df <- data.frame(x = groupBy, y = value)

  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = value)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(. ~ x)

  if (defaultTheme == TRUE) {
    p <- .ggSCTKTheme(p) +
      ggplot2::theme(strip.background = ggplot2::element_blank())
  }

  if (all(unique(groupBy) == "Sample")) {
    p <- p + ggplot2::theme(strip.text.x = ggplot2::element_blank())
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(label = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = titleSize
      ))
  }

  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize))
  }

  if (!is.null(ylab)) {
    p <- p + ggplot2::ylab(ylab) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize))
  }
  p <- p + ggplot2::theme(axis.text = ggplot2::element_text(size = axisSize))

  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = cutoff, color = "red")
  }

  return(p)
}

#' @title Density plot of colData.
#' @description Visualizes values stored in the colData slot of a
#'  SingleCellExperiment object via a density plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param coldata colData value that will be plotted.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param cutoff Numeric value. The plot will be annotated with a vertical line
#'  if set. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the density plot of colData.
#' @examples
#' plotSCEDensityColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEDensityColData <- function(inSCE,
                                  sample = NULL,
                                  coldata,
                                  groupBy = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  axisSize = 10,
                                  axisLabelSize = 10,
                                  defaultTheme = TRUE,
                                  title = NULL,
                                  titleSize = 18,
                                  cutoff = NULL,
                                  combinePlot = "none",
                                  plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!is.null(coldata)) {
    if (!coldata %in% names(SummarizedExperiment::colData(inSCE))) {
      stop("'", paste(coldata), "' is not found in ColData.")
    }
    coldata <- SummarizedExperiment::colData(inSCE)[, coldata]
  } else {
    stop("You must define the desired colData to plot.")
  }

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(coldata)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop(
        "'sample' must be the same length as the number",
        " of columns in 'inSCE'"
      )
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    coldataSub <- coldata[sampleInd]
    if (!is.null(groupBy)) {
      groupbySub <- groupBy[sampleInd]
    } else {
      groupbySub <- NULL
    }

    if (!is.null(title) && length(samples) > 1) {
      title <- paste(title, x, sep = "_")
    }
    p <- .ggDensity(
      value = coldataSub,
      groupBy = groupbySub,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      defaultTheme = defaultTheme,
      title = title,
      titleSize = titleSize,
      cutoff = cutoff
    )
    return(p)
  })
  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }

  return(plotlist)
}

#' @title Density plot of assay data.
#' @description Visualizes values stored in the assay slot of a
#'  SingleCellExperiment object via a density plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param feature Name of feature stored in assay of SingleCellExperiment
#'  object.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param cutoff Numeric value. The plot will be annotated with a vertical line
#'  if set. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the density plot of assay data.
#' @examples
#' plotSCEDensityAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe"
#' )
#' @export
plotSCEDensityAssayData <- function(inSCE,
                                    sample = NULL,
                                    useAssay = "counts",
                                    feature,
                                    groupBy = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    axisSize = 10,
                                    axisLabelSize = 10,
                                    defaultTheme = TRUE,
                                    cutoff = NULL,
                                    title = NULL,
                                    titleSize = 18,
                                    combinePlot = "none",
                                    plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(counts)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop(
        "'sample' must be the same length as the number",
        " of columns in 'inSCE'"
      )
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    countsSub <- counts[sampleInd]
    if (!is.null(groupBy)) {
      groupbySub <- groupBy[sampleInd]
    } else {
      groupbySub <- NULL
    }

    if (!is.null(title) && length(samples) > 1) {
      title <- paste(title, x, sep = "_")
    }

    p <- .ggDensity(
      value = countsSub,
      groupBy = groupbySub,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      defaultTheme = defaultTheme,
      title = title,
      titleSize = titleSize,
      cutoff = cutoff
    )
    return(p)
  })

  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }
  return(plotlist)
}

#' @title Density plot of any data stored in the SingleCellExperiment object.
#' @description Visualizes values stored in any slot of a
#'  SingleCellExperiment object via a densityn plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata"
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param annotation Desired vector within the slot used for plotting.
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object. Will be used only if "assays" slot is chosen. Default NULL.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param cutoff Numeric value. The plot will be annotated with a vertical line
#'  if set. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot object of the density plot.
#' @examples
#' plotSCEDensity(
#'   inSCE = mouseBrainSubsetSCE, slot = "assays",
#'   annotation = "counts", feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEDensity <- function(inSCE,
                           slot,
                           annotation,
                           sample = NULL,
                           useAssay = "counts",
                           feature = NULL,
                           groupBy = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           axisSize = 10,
                           axisLabelSize = 10,
                           defaultTheme = TRUE,
                           title = NULL,
                           titleSize = 18,
                           cutoff = NULL,
                           combinePlot = "none",
                           plotLabels = NULL) {
  if(!combinePlot %in% c("all", "sample", "none")){
    stop("'combinePlot' must be set to either 'all', 'sample', or 'none'.")
  }

  if (!slot %in% methods::slotNames(inSCE)) {
    stop("'slot' must be a slot within the SingleCellExperiment object.
             Please run 'methods::slotNames' if you are unsure the
             specified slot exists.")
  }

  sceSubset <- do.call(slot, args = list(inSCE))

  if (!annotation %in% names(sceSubset)) {
    stop("'annotation' must be an annotation stored within the specified
             slot of the SingleCellExperiment object.")
  }

  annotation.ix <- match(annotation, names(sceSubset))

  if (slot == "assays" && !is.null(feature)) {
    counts <- sceSubset[[annotation.ix]]
    if (feature %in% rownames(counts)) {
      counts <- counts[feature, ]
    }
  } else if (slot == "colData") {
    counts <- sceSubset[, annotation.ix]
  } else if (slot == "metadata") {
    counts <- sceSubset[[annotation.ix]]
  }

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(counts)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop(
        "'sample' must be the same length as the number",
        " of columns in 'inSCE'"
      )
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    countsSub <- counts[sampleInd]
    if (!is.null(groupBy)) {
      groupbySub <- groupBy[sampleInd]
    } else {
      groupbySub <- NULL
    }

    if (!is.null(title) && length(samples) > 1) {
      title <- paste(title, x, sep = "_")
    }

    p <- .ggDensity(
      value = countsSub,
      groupBy = groupbySub,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      defaultTheme = defaultTheme,
      title = title,
      titleSize = titleSize
    )
    return(p)
  })
  if(!is.null(feature)){
    names(plotlist) <- feature
  }

  ##Needs to be turned off for Shiny User Interface
  if(combinePlot %in% c("all", "sample")){
    figNcol = NULL
    if(!is.null(groupBy)){
      if(length(unique(groupBy)) > 1){
        figNcol = 1
      }
    }
    plotlist <- .ggSCTKCombinePlots(plotlist,
                                    combinePlot = combinePlot,
                                    ncols = figNcol,
                                    labels = plotLabels)
  }

  return(plotlist)
}

#' @title Plots for runEmptyDrops outputs.
#' @description A plotting function which visualizes outputs from the
#'  runEmptyDrops function stored in the colData slot of the SingleCellExperiment
#'  object via scatterplot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runEmptyDrops}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param fdrCutoff Numeric. Thresholds barcodes based on the FDR values from
#'  runEmptyDrops as "Empty Droplet" or "Putative Cell". Default 0.01.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 18.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 12.
#' @param axisLabelSize Size of x/y-axis labels. Default 15.
#' @param legendTitle Title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default 10.
#' @param combinePlot Boolean. If multiple plots are generated (multiple
#'  samples, etc.), will combined plots using `cowplot::plot_grid`.
#'  Default TRUE.
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return a ggplot object of the scatter plot.
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- runEmptyDrops(inSCE=sce)
#' plotEmptyDropsScatter(inSCE=sce)
#' @export
plotEmptyDropsScatter <- function(inSCE,
                                  sample = NULL,
                                  fdrCutoff = 0.01,
                                  defaultTheme = TRUE,
                                  dotSize = 1,
                                  title = NULL,
                                  titleSize = 18,
                                  xlab = NULL,
                                  ylab = NULL,
                                  axisSize = 12,
                                  axisLabelSize = 15,
                                  legendTitle = NULL,
                                  legendTitleSize = 12,
                                  legendSize = 10,
                                  combinePlot = "none",
                                  relHeights=1,
                                  relWidths=1,
                                  samplePerColumn = TRUE,
                                  sampleRelHeights = 1,
                                  sampleRelWidths = 1
){
  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop(
        "'sample' must be the same length as the number",
        " of columns in 'inSCE'"
      )
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  samples <- unique(sample)

  plotlist <- lapply(samples, function(x){
    sceSampleInd <- which(sample == x)
    inSCESub <- inSCE[, sceSampleInd]
    inSCESub = inSCESub[,!is.na(inSCESub$dropletUtils_emptyDrops_fdr)]
    isCell <- unlist(lapply(inSCESub$dropletUtils_emptyDrops_fdr, function(x){
      if(!is.na(x)){
        if(x <= fdrCutoff){
          return("Putative Cell")
        }else{
          return("Empty Droplet")
        }
      }

    }))

    df <- data.frame(x = inSCESub$dropletUtils_emptyDrops_total,
                     y = -(inSCESub$dropletUtils_emptyDrops_logprob),
                     isCell = isCell)

    p <- ggplot2::ggplot(df, ggplot2::aes_string("x",
                                                 "y", color = "isCell")) +
      ggplot2::geom_point(size = dotSize) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2))) +
      ggplot2::scale_color_manual(values = c("gray", "red"))

    if (defaultTheme == TRUE) {
      p <- .ggSCTKTheme(p)
    }

    if (!is.null(title)) {
      if(length(samples) > 1){
        title = paste(title, x, sep = "_")
      }
      p <- p + ggplot2::ggtitle(label = title) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = titleSize
        ))
    }
    if (!is.null(xlab)) {
      p <- p + ggplot2::xlab(xlab) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize),
                       axis.text.x = ggplot2::element_text(size = axisSize))
    }
    if (!is.null(ylab)) {
      p <- p + ggplot2::ylab(ylab) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize),
                       axis.text.y = ggplot2::element_text(size = axisSize))
    }
    if (!is.null(legendTitle)) {
      p <- p + ggplot2::labs(color = legendTitle) +
        ggplot2::theme(legend.title=ggplot2::element_text(size=legendTitleSize),
                       legend.text=ggplot2::element_text(size=legendSize))
    } else {
      p <- p + ggplot2::labs(color = "") +
        ggplot2::theme(legend.text=ggplot2::element_text(size=legendSize))
    }
    return(p)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- list(Sample = plotlist)
  } else {
    plotlist <- plotlist[[1]]
  }

  ##Needs to be turned off for Shiny User Interface
  if(!combinePlot == "none"){
    if(combinePlot == "all" && length(unique(samples)) > 1){
      return(cowplot::plot_grid(plotlist = unlist(plotlist,
                                                  recursive = FALSE),
                                align = "h",
                                vjust = 0,
                                rel_heights = sampleRelHeights,
                                rel_widths = sampleRelWidths))

    }else{
      return(plotlist)
    }
  }
  return(plotlist)
}


#' @title Plots for runBarcodeRankDrops outputs.
#' @description A plotting function which visualizes outputs from the
#'  runBarcodeRankDrops function stored in the colData slot of the SingleCellExperiment
#'  object via scatterplot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runBarcodeRankDrops}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 18.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 12.
#' @param axisLabelSize Size of x/y-axis labels. Default 15.
#' @param legendSize size of legend. Default 10.
#' @param combinePlot Boolean. If multiple plots are generated (multiple
#'  samples, etc.), will combined plots using `cowplot::plot_grid`.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#'  Default TRUE.
#' @return a ggplot object of the scatter plot.
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- runBarcodeRankDrops(inSCE=sce)
#' plotBarcodeRankScatter(inSCE=sce)
#' @export
plotBarcodeRankScatter <- function(inSCE,
                                   sample = NULL,
                                   defaultTheme = TRUE,
                                   dotSize = 1,
                                   title = NULL,
                                   titleSize = 18,
                                   xlab = NULL,
                                   ylab = NULL,
                                   axisSize = 12,
                                   axisLabelSize = 15,
                                   legendSize = 10,
                                   combinePlot = "none",
                                   sampleRelHeights = 1,
                                   sampleRelWidths = 1){
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop(
                "'sample' must be the same length as the number",
                " of columns in 'inSCE'"
            )
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)
    meta <- S4Vectors::metadata(inSCE)$runBarcodeRanksMetaOutput
  plotlist <- lapply(samples, function(x){

    sampleMeta <- meta[[x]]
    knee <- sampleMeta$dropletUtils_barcodeRank_knee
    inflection <- sampleMeta$dropletUtils_barcodeRank_inflection
    df <- data.frame(rank = sampleMeta$dropletUtils_barcodeRank_rank,
                     umi = sampleMeta$dropletUtils_barcodeRank_total)


    p <- ggplot2::ggplot(df, ggplot2::aes_string(x="rank", y="umi")) +
      ggplot2::geom_point(size=dotSize, shape=20) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10()

    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept=knee, linetype = "Knee"), colour = 'red') +
      ggplot2::geom_hline(ggplot2::aes(yintercept=inflection, linetype = "Inflection"), colour= 'blue') +
      ggplot2::scale_linetype_manual(name = "", values = c(2, 2),
                                     guide = ggplot2::guide_legend(label.theme = ggplot2::element_text(size = legendSize),
                                                                   override.aes = list(color = c("blue", "red"))))

    if (defaultTheme == TRUE) {
      p <- .ggSCTKTheme(p)
    }
    if (!is.null(title)) {
      p <- p + ggplot2::ggtitle(label = title) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = titleSize
        ))
    }
    if (!is.null(xlab)) {
      p <- p + ggplot2::xlab(xlab) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize),
                       axis.text.x = ggplot2::element_text(size = axisSize))
    }else{
      p <- p + ggplot2::xlab("Rank") +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize),
                       axis.text.x = ggplot2::element_text(size = axisSize))
    }

    if (!is.null(ylab)) {
      p <- p + ggplot2::ylab(ylab) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize),
                       axis.text.y = ggplot2::element_text(size = axisSize))
    }else{
      p <- p + ggplot2::ylab("Total UMI Counts") +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize),
                       axis.text.y = ggplot2::element_text(size = axisSize))
    }
    return(p)
  })
  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- list(Sample = plotlist)
  } else {
      plotlist <- plotlist[[1]]
  }

  ##Needs to be turned off for Shiny User Interface
  if(!combinePlot == "none"){
      if(combinePlot %in% c("all") && length(unique(sample)) > 1){
          return(cowplot::plot_grid(plotlist = unlist(plotlist,
          recursive = FALSE),
          align = "h",
          vjust = 0,
          rel_heights = sampleRelHeights,
          rel_widths = sampleRelWidths))
      }else if(combinePlot == "sample"){
          return(plotlist)
      }
  }
  return(plotlist)

}

#' @title Bar plot plotting tool.
#' @description Visualizes specified values via a violin plot.
#' @param y Numeric values to be plotted on y-axis.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @return a ggplot of the reduced dimensions.
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
.ggBar <- function(y,
                   groupBy = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   axisSize = 10,
                   axisLabelSize = 10,
                   dotSize = 1,
                   transparency = 1,
                   defaultTheme = TRUE,
                   gridLine = FALSE,
                   summary = NULL,
                   title = NULL,
                   titleSize = 15) {
  if (is.null(groupBy)) {
    groupBy <- rep("Sample", length(y))
  }

  df <- data.frame(x = groupBy, y = y)

  p <- ggplot2::ggplot(df) +
    ggplot2::aes_string(
      x = "groupBy",
      y = "y"
    )

  p <- p + ggplot2::geom_bar(stat = "identity")

  if (defaultTheme == TRUE) {
    p <- .ggSCTKTheme(p)
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(label = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = titleSize
      ))
  }

  ###
  p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = axisSize))
  ###

  if(length(unique(df$groupBy)) > 1){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                hjust = 1,
                                                                size = axisSize))
  }else{
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank())
  }

  if (gridLine == TRUE){
    p <- p + ggplot2::theme(panel.grid.major.y = ggplot2::element_line("grey"))
  }
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisLabelSize))
  }
  if (!is.null(ylab)) {
    p <- p + ggplot2::ylab(ylab) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisLabelSize))
  }
  if (!is.null(summary)){
    if(summary == "mean"){
      summ <- df %>% dplyr::group_by(groupBy) %>% dplyr::summarize(value = base::mean(y))
      fun <- base::mean
    }else if(summary == "median"){
      summ <- df %>% dplyr::group_by(groupBy) %>% dplyr::summarize(value = stats::median(y))
      fun <- stats::median
    }else{
      stop("`summary`` must be either `mean` or `median`.")
    }
    summ$statY <-  max(df$y) + (max(df$y) - min(df$y)) * 0.1
    summary <- paste(toupper(substr(summary, 1, 1)),
                     substr(summary, 2, nchar(summary)), sep="")
    summ$label <- paste0(summary,": ", round(summ$value, 5))

    p <- p + ggplot2::geom_text(data = summ,
      ggplot2::aes_string(x = "groupBy",
        y = "statY",
        label = "label"),
      size = 5)
    p <- p + ggplot2::stat_summary(fun = fun, fun.min = fun,
                                   fun.max = fun,
                                   geom = "crossbar",
                                   color = "red",
                                   linetype = "dashed")
  }

  return(p)
}



#' @title Bar plot of colData.
#' @description Visualizes values stored in the colData slot of a
#'  SingleCellExperiment object via a bar plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param coldata colData value that will be plotted.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param combinePlot Boolean. If multiple plots are generated (multiple
#'  samples, etc.), will combined plots using `cowplot::plot_grid`.
#'  Default TRUE.
#' @return a ggplot of the barplot of coldata.
#' @examples
#' plotSCEBarColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEBarColData <- function(inSCE,
                              sample = NULL,
                              coldata,
                              groupBy = NULL,
                              dots = TRUE,
                              xlab = NULL,
                              ylab = NULL,
                              axisSize = 10,
                              axisLabelSize = 10,
                              dotSize = 1,
                              transparency = 1,
                              defaultTheme = TRUE,
                              gridLine = FALSE,
                              summary = NULL,
                              title = NULL,
                              titleSize = NULL,
                              combinePlot = TRUE) {
  if (!is.null(coldata)) {
    if (!coldata %in% names(SummarizedExperiment::colData(inSCE))) {
      stop("'", paste(coldata), "' is not found in ColData.")
    }
    coldata <- SummarizedExperiment::colData(inSCE)[, coldata]
  } else {
    stop("You must define the desired colData to plot.")
  }

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(coldata)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number",
           " of columns in 'inSCE'")
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  p <- .ggBar(
    y = coldata,
    groupBy = groupBy,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    axisLabelSize = axisLabelSize,
    dotSize = dotSize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize
  )

  return(p)
}

#' @title Bar plot of assay data.
#' @description Visualizes values stored in the assay slot of a
#'  SingleCellExperiment object via a bar plot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param feature Name of feature stored in assay of SingleCellExperiment
#'  object.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param combinePlot Boolean. If multiple plots are generated (multiple
#'  samples, etc.), will combined plots using `cowplot::plot_grid`.
#'  Default TRUE.
#' @return a ggplot of the barplot of assay data.
#' @examples
#' plotSCEBarAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEBarAssayData <- function(inSCE,
                                sample = NULL,
                                useAssay = "counts",
                                feature,
                                groupBy = NULL,
                                xlab = NULL,
                                ylab = NULL,
                                axisSize = 10,
                                axisLabelSize = 10,
                                dotSize = 1,
                                transparency = 1,
                                defaultTheme = TRUE,
                                gridLine = FALSE,
                                summary = NULL,
                                title = NULL,
                                titleSize = NULL,
                                combinePlot = TRUE) {
  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]

  if (!is.null(groupBy)) {
    if (length(groupBy) > 1) {
      if (length(groupBy) != length(counts)) {
        stop("The input vector for 'groupBy' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupBy %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupBy), "' is not found in ColData.")
      }
      groupBy <- as.character(SummarizedExperiment::colData(inSCE)[, groupBy])
    }
  }

  p <- .ggBar(
    y = counts,
    groupBy = groupBy,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    axisLabelSize = axisLabelSize,
    dotSize = dotSize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize
  )

  return(p)
}

.binSCTK <- function(value, bin, binLabel = NULL) {
  if (!is.null(binLabel)) {
    if (length(bin) == 1) {
      if (bin != length(binLabel)) {
        stop("'binLabel' must be equal to the bin length")
      }
    } else if (length(bin) > 1) {
      if (bin != length(binLabel) + 1) {
        stop("'binLabel' must be equal to the bin length")
      }
    }
  }
  value.bin <- cut(x = value, breaks = bin, labels = binLabel)
  return(value.bin)
}

.ggSCTKCombinePlots <- function(plotlist,
                                ncols = NULL,
                                nrows = NULL,
                                combinePlot = "all",
                                relHeights = 1,
                                relWidths = 1,
                                labels = "default",
                                labelPositionX = NULL,
                                labelPositionY = NULL,
                                labelSize = 20,
                                samplePerColumn = TRUE,
                                sampleRelHeights = 1,
                                sampleRelWidths = 1) {
  if ("Violin" %in% names(plotlist)) {
    plotlistViolin <- plotlist$Violin
  } else {
    plotlistViolin <- NULL
  }

  if ("Sample" %in% names(plotlist)) {
    plotlistSample <- plotlist$Sample
    if (samplePerColumn) {
      ncols <- 1
      nrowSub <- 1
      sampleRelHeights <- 1
    }else{
      nrowSub = NULL
    }
    plotlistSample <- lapply(plotlistSample, function(x) {
      if(all(class(x) %in% c("gg","ggplot"))){
        return(x)
      }else if (inherits(x, "list")){
        return(cowplot::plot_grid(
          plotlist = x,
          align = "h",
          nrow = nrowSub,
          vjust = 0,
          rel_heights = sampleRelHeights,
          rel_widths = sampleRelWidths
        ))
      }
    })
  }else{
    plotlistSample <- NULL
  }

  if(!is.null(plotlistViolin) | !is.null(plotlistSample)){
    plotlist <- c(plotlistViolin, plotlistSample)
  }
  # To make the resulting plot close to a square as possible
  if (is.null(ncols) && is.null(nrows)) {
    ncols <- round(sqrt(length(plotlist)))
  }

  if (combinePlot == "all") {
    if (!is.null(labels) && labels != "none") {
      # If default, sample name will be used as labels
      if(length(labels) == 1){
        if (labels == "default") {
          if(!is.null(names(plotlist))){
            labels <- names(plotlist)
          }
        }
      }

      listNamePlot <- list()

      if(is.null(labelPositionX)){
        labelPositionX = rep(0, length(plotlist))
      }
      if(is.null(labelPositionY)){
        labelPositionY = rep(1, length(plotlist))
      }

      for(x in seq_along(plotlist)){
        labeled <- plotlist[[x]] + cowplot::draw_plot_label(
          labels[x],
          x = labelPositionX[x],
          y = labelPositionY[x],
          size = labelSize)
        listNamePlot[[x]] <- labeled
      }

      plotlist <- listNamePlot
    }

    plotRes <- cowplot::plot_grid(
      plotlist = plotlist,
      ncol = ncols,
      nrow = nrows,
      rel_heights = relHeights,
      rel_widths = relWidths
    )

    return(plotRes)
  } else if (combinePlot == "sample") {
    #Will happen if "sample" is chosen and multiple samples exist,
    #whcih means there will be a plotlistViolin object
    if (!is.null(plotlistViolin)) {
      return(list(Violin = plotlistViolin, Sample = plotlistSample))
      #Will happen when calling the non-QC plot fxns (ie, plotSCEScatter, plotSCEViolin)
      #for multiple samples, meaning no merging has occurred across plots
    }else if (!is.null(plotlistSample)){
      return(plotlistSample)
      # #Will happen when?
      # }else if(length(plotlist) == 1) {
      #     return(plotlist[[1]])
      #Will happen when sample = NULL, combinePlot = "sample", meaning up to this point
      #the "plotlist" should be a list of individual plots for only one sample
    } else{
      return(cowplot::plot_grid(
        plotlist = plotlist,
        ncol = ncols,
        nrow = nrows,
        rel_heights = relHeights,
        rel_widths = relWidths
      ))
    }
  }
}
.ggSCTKTheme <- function(gg) {
  return(gg + ggplot2::theme_bw() +
           ggplot2::theme(
             panel.grid.major = ggplot2::element_blank(),
             panel.grid.minor = ggplot2::element_blank(),
             axis.text = ggplot2::element_text(size = 10),
             axis.title = ggplot2::element_text(size = 10)
           ))
}
