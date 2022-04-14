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
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
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
#' @param dotSize Size of dots. Default 0.5.
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
#' @param clusterLabelSize Numeric. Determines the size of cluster label
#'  when `labelClusters` is set to TRUE. Default 3.5.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimensions.
.ggScatter <- function(inSCE,
                       reducedDimName,
                       sample = NULL,
                       colorBy = NULL,
                       groupBy = NULL,
                       shape = NULL,
                       conditionClass = NULL,
                       labelClusters = FALSE,
                       clusterLabelSize = 3.5,
                       xlab = NULL,
                       ylab = NULL,
                       baseSize = 12,
                       axisSize = NULL,
                       axisLabelSize = NULL,
                       dim1 = NULL,
                       dim2 = NULL,
                       bin = NULL,
                       binLabel = NULL,
                       dotSize = 0.5,
                       transparency = 1,
                       colorScale = NULL,
                       colorLow = "white",
                       colorMid = "gray",
                       colorHigh = "blue",
                       defaultTheme = TRUE,
                       title = NULL,
                       titleSize = NULL,
                       legendTitle = NULL,
                       legendTitleSize = NULL,
                       legendSize = NULL,
                       combinePlot = "none",
                       plotLabels = NULL) {
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

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
      g <- .ggSCTKTheme(g, baseSize, groupBy = factor(sample),
                        combinePlot)
    }else{
      g <- g + ggplot2::theme_gray(base_size = baseSize)
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
          color = "black",
          size = clusterLabelSize
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
#' dimension reduction components or a variable with saved results. Required.
#' @param reducedDimName Saved dimension reduction matrix name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param colorBy Color by a condition(any column of the annotation data).
#'  Required.
#' @param groupBy Group by a condition(any column of the annotation data).
#'  Default NULL.
#' @param conditionClass Class of the annotation data used in colorBy.
#'  Options are NULL, "factor" or "numeric". If NULL, class will default to the
#'  original class. Default NULL.
#' @param shape Add shapes to each condition.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param dim1 1st dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Can either be a string which specifies
#'  the name of the dimension to be plotted from reducedDims, or a numeric value which specifies
#'  the index of the dimension to be plotted. Default is NULL.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 0.5.
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
#' @param clusterLabelSize Numeric. Determines the size of cluster label
#'  when `labelClusters` is set to TRUE. Default 3.5.
#' @param legendTitle title of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default 12.
#' @param legendSize size of legend. Default NULL.
#'  Default FALSE.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the reduced dimension plot of coldata.
#' @examples
#' data("mouseBrainSubsetSCE")
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
                                    colorBy,
                                    reducedDimName,
                                    sample = NULL,
                                    groupBy = NULL,
                                    conditionClass = NULL,
                                    shape = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    baseSize = 12,
                                    axisSize = NULL,
                                    axisLabelSize = NULL,
                                    dim1 = NULL,
                                    dim2 = NULL,
                                    bin = NULL,
                                    binLabel = NULL,
                                    dotSize = 0.5,
                                    transparency = 1,
                                    colorScale = NULL,
                                    colorLow = "white",
                                    colorMid = "gray",
                                    colorHigh = "blue",
                                    defaultTheme = TRUE,
                                    title = NULL,
                                    titleSize = 15,
                                    labelClusters = TRUE,
                                    clusterLabelSize = 3.5,
                                    legendTitle = NULL,
                                    legendTitleSize = NULL,
                                    legendSize = NULL,
                                    combinePlot = "none",
                                    plotLabels = NULL) {
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

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
    baseSize = baseSize,
    title = title,
    titleSize = titleSize,
    labelClusters = labelClusters,
    clusterLabelSize = clusterLabelSize,
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
#' dimension reduction components or a variable with saved results. Required.
#' @param reducedDimName saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param feature Name of feature stored in assay of SingleCellExperiment
#'  object.
#' @param featureLocation Indicates which column name of rowData to query gene.
#' @param featureDisplay Indicates which column name of rowData to use
#' to display feature for visualization.
#' @param shape add shapes to each condition. Default NULL.
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
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 0.5.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEDimReduceFeatures(
#'   inSCE = mouseBrainSubsetSCE, feature = "Apoe",
#'   shape = NULL, reducedDimName = "TSNE_counts",
#'   useAssay = "counts", xlab = "tSNE1", ylab = "tSNE2"
#' )
#' @export
plotSCEDimReduceFeatures <- function(inSCE,
                                     feature,
                                     reducedDimName,
                                     sample = NULL,
                                     featureLocation = NULL,
                                     featureDisplay = NULL,
                                     shape = NULL,
                                     useAssay = "logcounts",
                                     xlab = NULL,
                                     ylab = NULL,
                                     axisSize = 10,
                                     axisLabelSize = 10,
                                     dim1 = NULL,
                                     dim2 = NULL,
                                     bin = NULL,
                                     binLabel = NULL,
                                     dotSize = 0.5,
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
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

  if(!is.null(featureDisplay)){
    featureDisplay <- match.arg(featureDisplay,
                                colnames(SummarizedExperiment::rowData(inSCE)))
  }else{
    if(exists(x = "featureDisplay", inSCE@metadata)){
      featureDisplay <- inSCE@metadata$featureDisplay
    }
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous",
    featureLocation = featureLocation,
    featureDisplay = featureDisplay
  )
  counts <- mat[, 2]

  if(!is.null(featureDisplay)){
    title = utils::tail(colnames(mat),1)
  }

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
#'  components or a variable with saved results. Required.
#' @param reducedDimName saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata", "reducedDims". Default NULL.
#' @param annotation Desired vector within the slot used for plotting. Default NULL.
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object. Will be used only if "assays" slot is chosen. Default NULL.
#' @param groupBy Group by a condition(any column of the annotation data).
#'  Default NULL.
#' @param shape add shapes to each condition.
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
#' @param bin Numeric vector. If single value, will divide the numeric values into the `bin` groups.
#'  If more than one value, will bin numeric values using values as a cut point.
#' @param binLabel Character vector. Labels for the bins created by the `bin` parameter.
#'  Default NULL.
#' @param dotSize Size of dots. Default 0.5.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEScatter(
#'   inSCE = mouseBrainSubsetSCE, legendTitle = NULL,
#'   slot = "assays", annotation = "counts", feature = "Apoe",
#'   reducedDimName = "TSNE_counts", labelClusters = FALSE
#' )
#' @export
#' @import SingleCellExperiment
plotSCEScatter <- function(inSCE,
                           annotation,
                           reducedDimName = NULL,
                           slot = NULL,
                           sample = NULL,
                           feature = NULL,
                           groupBy = NULL,
                           shape = NULL,
                           conditionClass = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           axisSize = 10,
                           axisLabelSize = 10,
                           dim1 = NULL,
                           dim2 = NULL,
                           bin = NULL,
                           binLabel = NULL,
                           dotSize = 0.5,
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
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

  if (!is.null(slot)){
    if (slot == "reducedDims"){
      annotation_clm <- substr(annotation, stringr::str_length(annotation), stringr::str_length(annotation))
      annotation <- substr(annotation, 1, stringr::str_length(annotation) - 2)
    }else if (!slot %in% methods::slotNames(inSCE)) {
      stop("'slot' must be a slot within the SingleCellExperiment object.",
           "Please run 'methods::slotNames' if you are unsure the",
           "specified slot exists.")
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
#' @param plotOrder Character vector. If set, reorders the violin plots
#'  in the order of the character vector when `groupBy` is set.
#'  Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param dotSize Size of dots. Default 0.5.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param hcutoff Adds a horizontal line with the y-intercept at given value. Default NULL.
#' @param hcolor Character. A color available from `colors()`.
#'  Controls the color of the horizontal cutoff line, if drawn.
#'  Default 'black'.
#' @param hsize Size of horizontal line, if drawn. Default 0.5.
#' @param hlinetype Type of horizontal line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param vcutoff Adds a vertical line with the x-intercept at given value. Default NULL.
#' @param vcolor Character. A color available from `colors()`.
#'  Controls the color of the vertical cutoff line, if drawn.
#'  Default 'black'.
#' @param vsize Size of vertical line, if drawn. Default 0.5.
#' @param vlinetype Type of vertical line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @return a ggplot of the reduced dimensions.
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
.ggViolin <- function(y,
                      groupBy = NULL,
                      violin = TRUE,
                      boxplot = TRUE,
                      dots = TRUE,
                      plotOrder = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      baseSize = 12,
                      axisSize = NULL,
                      axisLabelSize = NULL,
                      dotSize = 0.5,
                      transparency = 1,
                      defaultTheme = TRUE,
                      gridLine = FALSE,
                      summary = NULL,
                      summaryTextSize = 3,
                      combinePlot = "none",
                      title = NULL,
                      titleSize = NULL,
                      hcutoff = NULL,
                      hcolor = "red",
                      hsize = 1,
                      hlinetype = 1,
                      vcutoff = NULL,
                      vcolor = "red",
                      vsize = 1,
                      vlinetype = 1) {
  if (is.null(groupBy)) {
    groupBy <- rep("Sample", length(y))
  }


  if(!is.factor(groupBy)){
    if(is.null(plotOrder)){
      plotOrder = unique(groupBy)
    }
    groupBy <- factor(groupBy, levels = plotOrder)
  }else{
    if(!is.null(plotOrder)){
      groupBy <- factor(groupBy, levels = plotOrder)
    }
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
    p <- .ggSCTKTheme(p, baseSize, groupBy, combinePlot)
  }else{
    p <- p + ggplot2::theme_gray(base_size = baseSize)
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(label = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = titleSize
      ))
  }

  p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = axisSize))

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
    summ$statY <-  max(df$y) + (max(df$y) - min(df$y)) * 0.05
    summary <- paste(toupper(substr(summary, 1, 1)),
                     substr(summary, 2, nchar(summary)), sep="")

    ##Truncate label of mean/median if too many sample types
    if(length(levels(groupBy)) > 5){
      if(all(summ$value>1)){
        summ$label <- round(summ$value, 1)
      }else{
        summ$label <- signif(summ$value, 1)
      }
      p <- p + ggplot2::labs(subtitle = paste0(summary," values shown"))
    }else{
      if(all(summ$value>1)){
        summ$label <- paste0(summary,": ", round(summ$value, 2))
      }else{
        summ$label <- paste0(summary,": ", signif(summ$value, 2))
      }
    }

    if(!is.null(groupBy)){
      summaryTextSize = summaryTextSize/length(levels(groupBy)) + 2
    }

    p <- p + ggplot2::geom_text(data = summ, size = summaryTextSize,
                                ggplot2::aes_string(x = "groupBy",
                                                    y = "statY",
                                                    label = "label"))
    p <- p + ggplot2::stat_summary(fun = fun, fun.min = fun,
                                   fun.max = fun,
                                   geom = "crossbar",
                                   color = "red",
                                   linetype = "dashed")
  }
  if(!is.null(hcutoff)){
    p <- .ggAddLine(p, hcutoff = hcutoff, hcolor = hcolor,
                    hsize = hsize, hlinetype = hlinetype)
  }
  if(!is.null(vcutoff)){
    p <- .ggAddLine(p, vcutoff = vcutoff, vcolor = vcolor,
                    vsize = vsize, vlinetype = vlinetype)
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
#' @param plotOrder Character vector. If set, reorders the violin plots
#'  in the order of the character vector when `groupBy` is set.
#'  Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param dotSize Size of dots. Default 0.5.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param hcutoff Adds a horizontal line with the y-intercept at given value. Default NULL.
#' @param hcolor Character. A color available from `colors()`.
#'  Controls the color of the horizontal cutoff line, if drawn.
#'  Default 'black'.
#' @param hsize Size of horizontal line, if drawn. Default 0.5.
#' @param hlinetype Type of horizontal line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param vcutoff Adds a vertical line with the x-intercept at given value. Default NULL.
#' @param vcolor Character. A color available from `colors()`.
#'  Controls the color of the vertical cutoff line, if drawn.
#'  Default 'black'.
#' @param vsize Size of vertical line, if drawn. Default 0.5.
#' @param vlinetype Type of vertical line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot of coldata.
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotSCEViolinColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEViolinColData <- function(inSCE,
                                 coldata,
                                 sample = NULL,
                                 groupBy = NULL,
                                 violin = TRUE,
                                 boxplot = TRUE,
                                 dots = TRUE,
                                 plotOrder = NULL,
                                 xlab = NULL,
                                 ylab = NULL,
                                 baseSize = 12,
                                 axisSize = NULL,
                                 axisLabelSize = NULL,
                                 dotSize = 0.5,
                                 transparency = 1,
                                 defaultTheme = TRUE,
                                 gridLine = FALSE,
                                 summary = NULL,
                                 summaryTextSize = 3,
                                 title = NULL,
                                 titleSize = NULL,
                                 hcutoff = NULL,
                                 hcolor = "red",
                                 hsize = 1,
                                 hlinetype = 1,
                                 vcutoff = NULL,
                                 vcolor = "red",
                                 vsize = 1,
                                 vlinetype = 1,
                                 combinePlot = "none",
                                 plotLabels = NULL) {
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

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
      title = paste(title, x, sep = ", ")
    }

    p <- .ggViolin(
      y = coldataSub,
      groupBy = groupbySub,
      violin = violin,
      boxplot = boxplot,
      dots = dots,
      plotOrder = plotOrder,
      xlab = xlab,
      ylab = ylab,
      baseSize=baseSize,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      dotSize = dotSize,
      transparency = transparency,
      defaultTheme = defaultTheme,
      gridLine = gridLine,
      summary = summary,
      summaryTextSize=summaryTextSize,
      combinePlot = combinePlot,
      title = title,
      titleSize = titleSize
    )
    if(!is.null(hcutoff)){
      p <- .ggAddLine(p, hcutoff = hcutoff, hcolor = hcolor,
                      hsize = hsize, hlinetype = hlinetype)
    }
    if(!is.null(vcutoff)){
      p <- .ggAddLine(p, vcutoff = vcutoff, vcolor = vcolor,
                      vsize = vsize, vlinetype = vlinetype)
    }
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
#' @param featureLocation Indicates which column name of rowData to query gene.
#' @param featureDisplay Indicates which column name of rowData to use
#' to display feature for visualization.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param plotOrder Character vector. If set, reorders the violin plots
#'  in the order of the character vector when `groupBy` is set.
#'  Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 0.5.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param hcutoff Adds a horizontal line with the y-intercept at given value. Default NULL.
#' @param hcolor Character. A color available from `colors()`.
#'  Controls the color of the horizontal cutoff line, if drawn.
#'  Default 'black'.
#' @param hsize Size of horizontal line, if drawn. Default 0.5.
#' @param hlinetype Type of horizontal line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param vcutoff Adds a vertical line with the x-intercept at given value. Default NULL.
#' @param vcolor Character. A color available from `colors()`.
#'  Controls the color of the vertical cutoff line, if drawn.
#'  Default 'black'.
#' @param vsize Size of vertical line, if drawn. Default 0.5.
#' @param vlinetype Type of vertical line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot of assay data.
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotSCEViolinAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEViolinAssayData <- function(inSCE,
                                   feature,
                                   sample = NULL,
                                   useAssay = "counts",
                                   featureLocation = NULL,
                                   featureDisplay = NULL,
                                   groupBy = NULL,
                                   violin = TRUE,
                                   boxplot = TRUE,
                                   dots = TRUE,
                                   plotOrder = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   axisSize = 10,
                                   axisLabelSize = 10,
                                   dotSize = 0.5,
                                   transparency = 1,
                                   defaultTheme = TRUE,
                                   gridLine = FALSE,
                                   summary = NULL,
                                   title = NULL,
                                   titleSize = NULL,
                                   hcutoff = NULL,
                                   hcolor = "red",
                                   hsize = 1,
                                   hlinetype = 1,
                                   vcutoff = NULL,
                                   vcolor = "red",
                                   vsize = 1,
                                   vlinetype = 1,
                                   combinePlot = "none",
                                   plotLabels = NULL) {
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

  if(!is.null(featureDisplay)){
    featureDisplay <- match.arg(featureDisplay,
                                colnames(SummarizedExperiment::rowData(inSCE)))
  }else{
    if(exists(x = "featureDisplay", inSCE@metadata)){
      featureDisplay <- inSCE@metadata$featureDisplay
    }
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    featureLocation = featureLocation,
    featureDisplay = featureDisplay,
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
  if(!is.null(featureDisplay) && is.null(title)){
    title = utils::tail(colnames(mat),1)
  }
  if(is.null(xlab)){
    ylab = "Expression"
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
      plotOrder = plotOrder,
      xlab = xlab,
      ylab = ylab,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      dotSize = dotSize,
      transparency = transparency,
      defaultTheme = defaultTheme,
      gridLine = gridLine,
      summary = summary,
      combinePlot = combinePlot,
      title = title,
      titleSize = titleSize
    )
    if(!is.null(hcutoff)){
      p <- .ggAddLine(p, hcutoff = hcutoff, hcolor = hcolor,
                      hsize = hsize, hlinetype = hlinetype)
    }
    if(!is.null(vcutoff)){
      p <- .ggAddLine(p, vcutoff = vcutoff, vcolor = vcolor,
                      vsize = vsize, vlinetype = vlinetype)
    }
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
#' dimension reduction components or a variable with saved results. Required.
#' @param slotName Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata", "reducedDims". Required.
#' @param itemName Desired vector within the slot used for plotting. Required.
#' @param feature Desired name of feature stored in assay of SingleCellExperiment
#'  object. Only used when "assays" slotName is selected. Default NULL.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param dimension Desired dimension stored in the specified reducedDims.
#'  Either an integer which indicates the column or a character vector specifies
#'  column name. By default, the 1st dimension/column will be used.
#'  Only used when "reducedDims" slotName is selected. Default NULL.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param plotOrder Character vector. If set, reorders the violin plots
#'  in the order of the character vector when `groupBy` is set.
#'  Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 0.5.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param gridLine Adds a horizontal grid line if TRUE. Will still be
#'  drawn even if defaultTheme is TRUE. Default FALSE.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param hcutoff Adds a horizontal line with the y-intercept at given value. Default NULL.
#' @param hcolor Character. A color available from `colors()`.
#'  Controls the color of the horizontal cutoff line, if drawn.
#'  Default 'black'.
#' @param hsize Size of horizontal line, if drawn. Default 0.5.
#' @param hlinetype Type of horizontal line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param vcutoff Adds a vertical line with the x-intercept at given value. Default NULL.
#' @param vcolor Character. A color available from `colors()`.
#'  Controls the color of the vertical cutoff line, if drawn.
#'  Default 'black'.
#' @param vsize Size of vertical line, if drawn. Default 0.5.
#' @param vlinetype Type of vertical line, if drawn. can be specified with either an integer or
#'  a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#'  Default 1.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @return a ggplot of the violin plot.
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotSCEViolin(
#'   inSCE = mouseBrainSubsetSCE, slotName = "assays",
#'   itemName = "counts", feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEViolin <- function(inSCE,
                          slotName,
                          itemName,
                          feature = NULL,
                          sample = NULL,
                          dimension = NULL,
                          groupBy = NULL,
                          violin = TRUE,
                          boxplot = TRUE,
                          dots = TRUE,
                          plotOrder = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          axisSize = 10,
                          axisLabelSize = 10,
                          dotSize = 0.5,
                          transparency = 1,
                          defaultTheme = TRUE,
                          gridLine = FALSE,
                          summary = NULL,
                          title = NULL,
                          titleSize = NULL,
                          hcutoff = NULL,
                          hcolor = "red",
                          hsize = 1,
                          hlinetype = 1,
                          vcutoff = NULL,
                          vcolor = "red",
                          vsize = 1,
                          vlinetype = 1,
                          combinePlot = "none",
                          plotLabels = NULL) {
    combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

    if (!slotName %in% c("rowData", "colData", "assays", "metadata", "reducedDims")) {
        stop("'slotName' must be a slotName within the SingleCellExperiment object.",
             "Please run 'methods::slot' if you are unsure the",
             "specified slotName exists.")
    }

    sceSubset <- do.call(slotName, args = list(inSCE))

    if (!itemName %in% names(sceSubset)) {
        stop("'itemName' must be an itemName stored within the specified
             slotName of the SingleCellExperiment object.")
    }

    itemName.ix <- match(itemName, names(sceSubset))

    if (slotName == "assays" && !is.null(feature)) {
        counts <- sceSubset[[itemName.ix]]
        if (feature %in% rownames(counts)) {
            counts <- counts[feature, ]
        }
    } else if (slotName == "colData") {
        counts <- sceSubset[, itemName.ix]
    } else if (slotName == "metadata") {
        counts <- sceSubset[[itemName.ix]]
    } else if (slotName == "reducedDims") {
        if(is.null(dimension)){
            dimension <- 1
        }else if(is.character(dimension)){
            dimension <- match(dimension, colnames(sceSubset[[itemName.ix]]))
        }
        counts <- sceSubset[[itemName.ix]][,dimension]
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
            plotOrder = plotOrder,
            xlab = xlab,
            ylab = ylab,
            axisSize = axisSize,
            axisLabelSize = axisLabelSize,
            dotSize = dotSize,
            transparency = transparency,
            defaultTheme = defaultTheme,
            gridLine = gridLine,
            summary = summary,
            combinePlot = combinePlot,
            title = title,
            titleSize = titleSize
        )
        if(!is.null(hcutoff)){
          p <- .ggAddLine(p, hcutoff = hcutoff, hcolor = hcolor,
                          hsize = hsize, hlinetype = hlinetype)
        }
        if(!is.null(vcutoff)){
          p <- .ggAddLine(p, vcutoff = vcutoff, vcolor = vcolor,
                          vsize = vsize, vlinetype = vlinetype)
        }
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
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title Title of plot. Default NULL.
#' @param titleSize Size of title of plot. Default 15.
#' @param cutoff Numeric value. The plot will be annotated with a vertical line
#'  if set. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "none".
#' @return density plot, in .ggplot.
.ggDensity <- function(value,
                       groupBy = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       baseSize = 12,
                       axisSize = NULL,
                       axisLabelSize = NULL,
                       defaultTheme = TRUE,
                       title = NULL,
                       titleSize = NULL,
                       combinePlot = "none",
                       cutoff = NULL) {
  if (is.null(groupBy)) {
    groupBy <- rep("Sample", length(value))
  }
  groupBy <- factor(groupBy, levels = unique(groupBy))
  df <- data.frame(x = groupBy, y = value)

  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = value)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(. ~ x)

  if (defaultTheme == TRUE) {
    p <- .ggSCTKTheme(p, baseSize, groupBy, combinePlot) +
      ggplot2::theme(strip.background = ggplot2::element_blank())
  }else{
    p <- p + ggplot2::theme_gray(base_size = baseSize)
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
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEDensityColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEDensityColData <- function(inSCE,
                                  coldata,
                                  sample = NULL,
                                  groupBy = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  baseSize = 12,
                                  axisSize = NULL,
                                  axisLabelSize = NULL,
                                  defaultTheme = TRUE,
                                  title = NULL,
                                  titleSize = 18,
                                  cutoff = NULL,
                                  combinePlot = "none",
                                  plotLabels = NULL) {
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

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
      title <- paste(title, x, sep = ", ")
    }
    p <- .ggDensity(
      value = coldataSub,
      groupBy = groupbySub,
      xlab = xlab,
      ylab = ylab,
      baseSize = baseSize,
      axisSize = axisSize,
      axisLabelSize = axisLabelSize,
      defaultTheme = defaultTheme,
      title = title,
      titleSize = titleSize,
      combinePlot = combinePlot,
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
#' @param featureLocation Indicates which column name of rowData to query gene.
#' @param featureDisplay Indicates which column name of rowData to use
#' to display feature for visualization.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEDensityAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe"
#' )
#' @export
plotSCEDensityAssayData <- function(inSCE,
                                    feature,
                                    sample = NULL,
                                    useAssay = "counts",
                                    featureLocation = NULL,
                                    featureDisplay = NULL,
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
  combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

  if(!is.null(featureDisplay)){
    featureDisplay <- match.arg(featureDisplay,
                                colnames(SummarizedExperiment::rowData(inSCE)))
  }else{
    if(exists(x = "featureDisplay", inSCE@metadata)){
      featureDisplay <- inSCE@metadata$featureDisplay
    }
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous",
    featureLocation = featureLocation,
    featureDisplay = featureDisplay
  )
  counts <- mat[, 2]

  if(!is.null(featureDisplay)){
    title = utils::tail(colnames(mat),1)
  }
  if(is.null(xlab)){
    xlab = "Expression"
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
#' dimension reduction components or a variable with saved results. Required.
#' @param slotName Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata", "reducedDims". Required.
#' @param itemName Desired vector within the slot used for plotting. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param feature Desired name of feature stored in assay of SingleCellExperiment
#'  object. Only used when "assays" slotName is selected. Default NULL.
#' @param dimension Desired dimension stored in the specified reducedDims.
#'  Either an integer which indicates the column or a character vector specifies
#'  column name. By default, the 1st dimension/column will be used.
#'  Only used when "reducedDims" slotName is selected. Default NULL.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEDensity(
#'   inSCE = mouseBrainSubsetSCE, slotName = "assays",
#'   itemName = "counts", feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEDensity <- function(inSCE,
                           slotName,
                           itemName,
                           sample = NULL,
                           feature = NULL,
                           dimension = NULL,
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
    combinePlot <- match.arg(combinePlot,c("all", "sample", "none"))

    if (!slotName %in% c("rowData", "colData", "assays", "metadata", "reducedDims")) {
        stop("'slotName' must be a slotName within the SingleCellExperiment object.",
             "Please run 'methods::slotNames' if you are unsure the",
             "specified slot exists.")
    }

    sceSubset <- do.call(slotName, args = list(inSCE))

    if (!itemName %in% names(sceSubset)) {
        stop("'itemName' must be an itemName stored within the specified
             slot of the SingleCellExperiment object.")
    }

    itemName.ix <- match(itemName, names(sceSubset))

    if (slotName == "assays" && !is.null(feature)) {
        counts <- sceSubset[[itemName.ix]]
        if (feature %in% rownames(counts)) {
            counts <- counts[feature, ]
        }
    } else if (slotName == "colData") {
        counts <- sceSubset[, itemName.ix]
    } else if (slotName == "metadata") {
        counts <- sceSubset[[itemName.ix]]
    } else if (slotName == "reducedDims") {
        if(is.null(dimension)){
            dimension <- 1
        }else if(is.character(dimension)){
            dimension <- match(dimension, colnames(sceSubset[[itemName.ix]]))
        }
        counts <- sceSubset[[itemName.ix]][,dimension]
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
    }else if(combinePlot == "none" && length(plotlist) == 1){
        plotlist <- plotlist[[1]]
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
#' @param dotSize Size of dots. Default 0.5.
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
                                  dotSize = 0.5,
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
#' @param dotSize Size of dots. Default 0.5.
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
                                   dotSize = 0.5,
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
#' @param dotSize Size of dots. Default 0.5.
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
                   dotSize = 0.5,
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

    p <- p + ggrepel::geom_text_repel(data = summ,
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
#' @param dotSize Size of dots. Default 0.5.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEBarColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupBy = "sex"
#' )
#' @export
plotSCEBarColData <- function(inSCE,
                              coldata,
                              sample = NULL,
                              groupBy = NULL,
                              dots = TRUE,
                              xlab = NULL,
                              ylab = NULL,
                              axisSize = 10,
                              axisLabelSize = 10,
                              dotSize = 0.5,
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
#' @param featureLocation Indicates which column name of rowData to query gene.
#' @param featureDisplay Indicates which column name of rowData to use
#' to display feature for visualization.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default 10.
#' @param axisLabelSize Size of x/y-axis labels. Default 10.
#' @param dotSize Size of dots. Default 0.5.
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
#' data("mouseBrainSubsetSCE")
#' plotSCEBarAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Apoe", groupBy = "sex"
#' )
#' @export
plotSCEBarAssayData <- function(inSCE,
                                feature,
                                sample = NULL,
                                useAssay = "counts",
                                featureLocation = NULL,
                                featureDisplay = NULL,
                                groupBy = NULL,
                                xlab = NULL,
                                ylab = NULL,
                                axisSize = 10,
                                axisLabelSize = 10,
                                dotSize = 0.5,
                                transparency = 1,
                                defaultTheme = TRUE,
                                gridLine = FALSE,
                                summary = NULL,
                                title = NULL,
                                titleSize = NULL,
                                combinePlot = TRUE) {
  if(!is.null(featureDisplay)){
    featureDisplay <- match.arg(featureDisplay,
                                colnames(SummarizedExperiment::rowData(inSCE)))
  }else{
    if(exists(x = "featureDisplay", inSCE@metadata)){
      featureDisplay <- inSCE@metadata$featureDisplay
    }
  }

  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous",
    featureLocation = featureLocation,
    featureDisplay = featureDisplay
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

#' @title Indicates which rowData to use for visualization
#' @description This function is to be used to specify which
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required.
#' @param featureDisplayRow Indicates which column name of rowData to be used for plots.
#' @return A SingleCellExperiment object with the specific column name of rowData
#'  to be used for plotting stored in metadata.
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- setSCTKDisplayRow(inSCE = sce, featureDisplayRow = "feature_name")
#' plotSCEViolinAssayData(inSCE = sce, feature = "ENSG00000019582")
#' @export
setSCTKDisplayRow <- function(inSCE,
                              featureDisplayRow) {
  inSCE@metadata$featureDisplay <- featureDisplayRow
  return(inSCE)
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
.ggSCTKTheme <- function(gg, baseSize = 12,
                         groupBy = NULL, combinePlot = "none") {

  scaleFactor <- .ggSetScaleFactor(groupBy = groupBy,
                    combinePlot = combinePlot)
  return(gg + ggplot2::theme_bw(base_size = baseSize * scaleFactor) +
           ggplot2::theme(
             panel.grid.major = ggplot2::element_blank(),
             panel.grid.minor = ggplot2::element_blank(),
             axis.text = ggplot2::element_text(),
             axis.title = ggplot2::element_text()
           ))
}

.ggSetScaleFactor <- function(groupBy = NULL,
                              combinePlot = "none"){
  if(!is.null(groupBy)){
    scaleFactor = 1/length(levels(groupBy)) + 0.5
  }else{
    scaleFactor = 1
  }
  if(combinePlot == "all"){
    scaleFactor = scaleFactor * 0.75
  }
  return(scaleFactor)
}

.ggAddLine <- function(plot,
                       hcutoff = NULL,
                       vcutoff = NULL,
                       hcolor = "red",
                       hsize = 1,
                       hlinetype = 1,
                       vcolor = "red",
                       vsize = 1,
                       vlinetype = 1){
  if(!is.null(hcutoff)){
    plot <- plot + ggplot2::geom_hline(yintercept = hcutoff, colour = hcolor,
                                       size = hsize, linetype = hlinetype)
  }
  if(!is.null(vcutoff)){
    plot <- plot + ggplot2::geom_vline(xintercept = vcutoff, colour = vcolor,
                                       size = vsize, linetype = vlinetype)
  }
  return(plot)
}
