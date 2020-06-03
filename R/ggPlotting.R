### Functions for plotting a scatter plot based on values stored in the
### SingleCellExperiment object (Ex. clinical data, QC info)

## Functions:
### plotSCEDimReduceColData - wrapper function; takes values from the colData
### slot in SingleCellExperiment object and passes it to .ggScatter fxn

### .ggScatter - Baseline scatter plot function.

### .ggSCTKTheme - Converts the ggplot to a specific theme (No background,
### no grid, etc.)

### .binSCTK - Bins input values based on the bin parameter.

## Example data can be found here: /restricted/projectnb/camplab/home/ykoga07/codereview/

#' @title Dimension reduction plot tool for colData
#' @description Plot results of reduced dimensions data and
#'  colors by annotation data stored in the colData slot.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param colorBy color by a condition(any column of the annotation data).
#' @param conditionClass class of the annotation data used in colorBy.
#'  Options are NULL, "factor" or "numeric". If NULL, class will default to the
#'  original class. Default NULL.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment
#'  object. Required.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#' @param dotsize size of dots. Default 2.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#' @param legendTitle title of legend. Default NULL.
#'  Default FALSE.
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotSCEDimReduceColData(
#'   inSCE = mouseBrainSubsetSCE, colorBy = "tissue",
#'   shape = "No Shape", conditionClass = "factor",
#'   reducedDimName = "TSNE_counts",
#'   xlab = "tSNE1", ylab = "tSNE2", labelClusters = TRUE
#' )
#'
#' plotSCEDimReduceColData(
#'   inSCE = mouseBrainSubsetSCE, colorBy = "age", labelClusters = FALSE,
#'   shape = "No Shape", bin = 3, binLabel = c("Young", "Medium", "Old")
#'   reducedDimName = "TSNE_counts",
#'   xlab = "tSNE1", ylab = "tSNE2"
#' )

plotSCEDimReduceColData <- function(inSCE,
  colorBy = "No Color",
  shape = "No Shape",
  reducedDimName = NULL,
  conditionClass = NULL,
  xlab = NULL,
  ylab = NULL,
  dim1 = NULL,
  dim2 = NULL,
  bin = NULL,
  binLabel = NULL,
  dotsize = 2,
  transparency = 1,
  defaultTheme = TRUE,
  title = NULL,
  titleSize = 15,
  labelClusters = TRUE,
  legendTitle = NULL,
  groupBy = NULL) {
  if (colorBy != "No Color"){
    colorPlot <- SingleCellExperiment::colData(inSCE)[, colorBy]
  }else{
    colorPlot <- "No Color"
  }

  g <- .ggScatter(
    inSCE = inSCE,
    colorBy = colorPlot,
    conditionClass = conditionClass,
    shape = shape,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    dim1 = dim1,
    dim2 = dim2,
    bin = bin,
    binLabel = binLabel,
    dotsize = dotsize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize,
    labelClusters = labelClusters,
    legendTitle = legendTitle,
    groupBy = groupBy
  )

  return(g)
}


.ggScatter <- function(inSCE,
  colorBy = "No Color",
  shape,
  reducedDimName,
  conditionClass = NULL,
  labelClusters = FALSE,
  xlab,
  ylab,
  dim1 = NULL,
  dim2 = NULL,
  bin = NULL,
  binLabel = NULL,
  dotsize = 2,
  transparency = 1,
  defaultTheme = TRUE,
  title = NULL,
  titleSize = 15,
  legendTitle = NULL,
  groupBy = NULL) {
  Df <- data.frame(SingleCellExperiment::reducedDim(
    inSCE,
    reducedDimName
  ))
  if (ncol(Df) > 2) {
    warning("More than two dimensions. Using the first two.")
  }
  if (!is.null(dim1) & !is.null(dim2)) {
    if (!(dim1 %in% colnames(Df))) {
      stop("X dimension ", dim1, " is not in the reducedDim data")
    }
    if (!(dim2 %in% colnames(Df))) {
      stop("Y dimension ", dim2, " is not in the reducedDim data")
    }
    xdim <- dim1
    ydim <- dim2
  } else if (!is.null(xlab) & !is.null(ylab)) {
    colnames(Df)[1] <- xlab
    colnames(Df)[2] <- ylab
    xdim <- colnames(Df)[1]
    ydim <- colnames(Df)[2]
  } else {
    colnames(Df)[1] <- paste0(reducedDimName, "_1")
    colnames(Df)[2] <- paste0(reducedDimName, "_2")
    xdim <- colnames(Df)[1]
    ydim <- colnames(Df)[2]
  }

  if (!is.null(conditionClass)){
    if (conditionClass %in% c("character", "factor", "numeric")) {
      if (conditionClass == "character") {
        colorBy <- as.character(colorBy)
      } else if (conditionClass == "factor") {
        colorBy <- as.factor(colorBy)
      } else if (conditionClass == "numeric") {
        colorBy <- as.numeric(colorBy)
      }
    }
  }

  if (!is.null(bin)){
    colorBy <- .binSCTK(value = colorBy,
      bin = bin,
      binLabel = binLabel)
  }

  if (length(colorBy) == 1) {
    if (colorBy == "No Color") {
      colorBy <- NULL
    }
  }
  if (shape == "No Shape") {
    shape <- NULL
  }
  if (!is.null(colorBy)) {
    Df$color <- colorBy
  }
  if (!is.null(shape)) {
    Df$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
  }
  Df$Sample <- colnames(inSCE)
  g <- ggplot2::ggplot(Df, ggplot2::aes_string(xdim, ydim,
    label = "Sample"
  )) +
    ggplot2::geom_point(size = dotsize, alpha = transparency)
  if (!is.null(colorBy)) {
    g <- g + ggplot2::aes_string(color = "color")
  }
  if (!is.null(shape)) {
    g <- g + ggplot2::aes_string(shape = "shape") +
      ggplot2::labs(shape = shape)
  }
  if (defaultTheme == TRUE) {
    g <- .ggSCTKTheme(g)
  }
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(label = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = titleSize
      ))
  }
  if (!is.null(legendTitle)) {
    g <- g + ggplot2::labs(color = legendTitle)
  } else {
    g <- g + ggplot2::labs(color = "")
  }

  if (!is.null(groupBy)){
    df$groups <- as.factor(colData(inSCE)@listData[[groupBy]] %>% data.frame())
    g <- g + facet_wrap(~groups)
  }

  if (isTRUE(labelClusters) && class(colorBy) %in% c("character", "factor")) {
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
      Sample = seq(1, length(unique(colorBy)))
    )

    if (!is.null(shape)) {
      centroid$shape <- Df$shape[1]
    }

    colnames(centroid)[seq(2)] <- c(xdim, ydim)
    g <- g + ggplot2::geom_point(
      data = centroid,
      mapping = ggplot2::aes_string(x = xdim, y = ydim),
      size = 0,
      alpha = 0
    ) +
      ggrepel::geom_text_repel(
        data = centroid,
        mapping = ggplot2::aes_string(label = "color"),
        show.legend = F,
        color = "black"
      )
  }

  return(g)
}

.ggSCTKTheme <- function(gg) {
  gg <- gg + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10)
    )
  return(gg)
}

.binSCTK <- function(value, bin, binLabel = NULL) {
  if(!is.null(binLabel)){
    if(length(bin) == 1){
      if(bin != length(binLabel)){
        stop("'binLabel' must be equal to the bin length")
      }
    }else if(length(bin) > 1){
      if(bin != length(binLabel)+1){
        stop("'binLabel' must be equal to the bin length")
      }
    }
  }

  # Since binning will only be done on numeric values, coerces
  #all values into numeric values.
  # Alternatively, could force users to only be able to use fxn
  #on numeric values
  value <- as.numeric(as.character(value))

  value.bin <- cut(x = value, breaks = bin, labels = binLabel)
  return(value.bin)
}


#' @title Dimension reduction plot tool for assay data
#' @description Plot results of reduced dimensions data and
#'  colors by feature data stored in the assays slot.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param feature name of feature stored in assay of singleCellExperiment
#'  object. Plot will be colored based on feature value.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment
#'  object. Required.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#'  Default is second PCA component for PCA data and NULL otherwise.
#' @param dotsize size of dots. Default 2.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @param legendTitle title of legend. Default NULL.
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotSCEDimReduceFeatures(
#'   inSCE = mouseBrainSubsetSCE, feature = "Sox2",
#'   shape = "No Shape", reducedDimName = "TSNE_counts",
#'   useAssay = "counts", xlab = "tSNE1", ylab = "tSNE2"
#' )
plotSCEDimReduceFeatures <- function(inSCE,
                                     feature,
                                     shape = "No Shape",
                                     reducedDimName,
                                     useAssay = "logcounts",
                                     xlab = NULL,
                                     ylab = NULL,
                                     dim1 = NULL,
                                     dim2 = NULL,
                                     dotsize = 2,
                                     transparency = 1,
                                     defaultTheme = TRUE,
                                     title = NULL,
                                     titleSize = 15,
                                     legendTitle = NULL,
                                     groupBy = NULL) {
  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]

  g <- .ggScatter(
    inSCE = inSCE,
    conditionClass = "numeric",
    colorBy = counts,
    shape = shape,
    transparency = 1,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    dim1 = dim1,
    dim2 = dim2,
    defaultTheme = defaultTheme,
    dotsize = dotsize,
    title = title,
    titleSize = titleSize,
    legendTitle = legendTitle,
    groupBy = groupBy
  )

  return(g)
}

#' @title Dimension reduction plot tool for all types of data
#' @description Plot results of reduced dimensions data of counts stored in any
#' slot in the SingleCellExperiment object.
#'
#' @param inSCE Input SingleCellExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata"
#' @param annotation Desired vector within the slot used for plotting.
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object. Will be used only if "assays" slot is chosen. Deafult NULL.
#' @param shape add shapes to each condition.
#' @param reducedDimName saved dimension reduction name in the SCtkExperiment
#'  object. Required.
#' @param conditionClass class of the annotation data used in colorBy. Options
#'  are NULL, "factor" or "numeric". If NULL, class will default to the original
#'  class. Default NULL.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param dim1 1st dimension to be used for plotting. Default is NULL.
#' @param dim2 2nd dimension to be used for plotting. Default is NULL.
#' Default is second PCA component for PCA data and NULL otherwise.
#' @param dotsize size of dots. Default 2.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme adds grid to plot when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#' @param legendTitle title of legend. Default NULL.
#' @return a ggplot of the reduced dimensions.
#' @export
#'
#' @examples
#' plotSCEScatter(
#'   inSCE = mouseBrainSubsetSCE, legendTitle = NULL,
#'   slot = "assays", annotation = "counts", feature = "Tspan12",
#'   reducedDimName = "TSNE_counts", labelClusters = FALSE
#' )
plotSCEScatter <- function(inSCE,
                           slot,
                           annotation,
                           feature = NULL,
                           shape = "No Shape",
                           reducedDimName = NULL,
                           conditionClass = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           dim1 = NULL,
                           dim2 = NULL,
                           dotsize = 2,
                           transparency = 1,
                           defaultTheme = TRUE,
                           title = NULL,
                           titleSize = 15,
                           labelClusters = TRUE,
                           legendTitle = NULL) {
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
      colorPlot <- counts[feature, ]
    }
  } else if (slot == "colData") {
    colorPlot <- sceSubset[, annotation.ix]
  } else if (slot == "metadata") {
    colorPlot <- sceSubset[[annotation.ix]]
  }

  if (!is.null(conditionClass)) {
    if (conditionClass %in% c("character", "factor", "numeric")) {
      if (conditionClass == "character") {
        colorPlot <- as.character(colorPlot)
      } else if (conditionClass == "factor") {
        colorPlot <- as.factor(colorPlot)
      } else if (conditionClass == "numeric") {
        colorPlot <- as.numeric(colorPlot)
      }
    }
  }
  g <- .ggScatter(
    inSCE = inSCE,
    colorBy = colorPlot,
    conditionClass = conditionClass,
    shape = shape,
    reducedDimName = reducedDimName,
    xlab = xlab,
    ylab = ylab,
    dim1 = dim1,
    dim2 = dim2,
    dotsize = dotsize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize,
    labelClusters = labelClusters,
    legendTitle = legendTitle
  )

  return(g)
}



#' @title Violin plot plotting tool.
#' @description Visualizes specified values via a violin plot.
#' @param y Numeric values to be plotted on y-axis.
#' @param groupby Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab label for x-axis. Default NULL.
#' @param ylab label for y-axis. Default NULL.
#' @param axisSize size of x/y-axis labels. Default 10.
#' @param dotSize size of dots. Default 1.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.
#' @return a ggplot of the reduced dimensions.
.ggViolin <- function(y,
                      groupby = NULL,
                      violin = TRUE,
                      boxplot = TRUE,
                      dots = TRUE,
                      xlab = NULL,
                      ylab = NULL,
                      axisSize = 10,
                      dotSize = 1,
                      transparency = 1,
                      defaultTheme = TRUE,
                      title = NULL,
                      titleSize = 15) {
  if (is.null(groupby)) {
    groupby <- rep("Sample", length(y))
  }
  df <- data.frame(x = groupby, y = y)

  p <- ggplot2::ggplot(df) +
    ggplot2::aes_string(
      x = "x",
      y = "y"
    )
  if (violin == TRUE) {
    p <- p + ggplot2::geom_violin(trim = TRUE, scale = "width")
  }
  if (boxplot == TRUE) {
    p <- p + ggplot2::geom_boxplot(width = 0.1)
  }
  if (dots == TRUE) {
    p <- p + ggplot2::geom_jitter(
      height = 0,
      size = dotSize,
      alpha = transparency
    )
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
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = axisSize))
  }
  if (!is.null(ylab)) {
    p <- p + ggplot2::ylab(ylab) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = axisSize))
  }

  return(p)
}

#' @title Violin plot of colData.
#' @description Visualizes values stored in the colData slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param coldata colData value that will be plotted.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab label for x-axis. Default NULL.
#' @param ylab label for y-axis. Default NULL.
#' @param axisSize size of x/y-axis labels. Default 10.
#' @param dotSize size of dots. Default 1.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.

#' @examples
#' plotSCEViolinColData(
#'   inSCE = mouseBrainSubsetSCE,
#'   coldata = "age", groupby = "sex"
#' )
#' @export
plotSCEViolinColData <- function(inSCE,
                                 coldata,
                                 groupby = NULL,
                                 violin = TRUE,
                                 boxplot = TRUE,
                                 dots = TRUE,
                                 xlab = NULL,
                                 ylab = NULL,
                                 axisSize = 10,
                                 dotSize = 1,
                                 transparency = 1,
                                 defaultTheme = TRUE,
                                 title = NULL,
                                 titleSize = NULL) {
  if (!is.null(coldata)) {
    if (!coldata %in% names(SummarizedExperiment::colData(inSCE))) {
      stop("'", paste(coldata), "' is not found in ColData.")
    }
    coldata <- SummarizedExperiment::colData(inSCE)[, coldata]
  } else {
    stop("You must define the desired colData to plot.")
  }

  if (!is.null(groupby)) {
    if (length(groupby) > 1) {
      if (length(groupby) != length(coldata)) {
        stop("The input vector for 'groupby' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupby %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupby), "' is not found in ColData.")
      }
      groupby <- as.character(SummarizedExperiment::colData(inSCE)[, groupby])
    }
  }


  p <- .ggViolin(
    y = coldata,
    groupby = groupby,
    violin = violin,
    boxplot = boxplot,
    dots = dots,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    dotSize = dotSize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize
  )

  return(p)
}

#' @title Violin plot of assay data.
#' @description Visualizes values stored in the assay slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param feature Name of feature stored in assay of SingleCellExperiment
#'  object.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab label for x-axis. Default NULL.
#' @param ylab label for y-axis. Default NULL.
#' @param axisSize size of x/y-axis labels. Default 10.
#' @param dotSize size of dots. Default 1.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.

#' @examples
#' plotSCEViolinAssayData(
#'   inSCE = mouseBrainSubsetSCE,
#'   feature = "Sox2", groupby = "sex"
#' )
#' @export
plotSCEViolinAssayData <- function(inSCE,
                                   useAssay = "counts",
                                   feature,
                                   groupby = NULL,
                                   violin = TRUE,
                                   boxplot = TRUE,
                                   dots = TRUE,
                                   xlab = NULL,
                                   ylab = NULL,
                                   axisSize = 10,
                                   dotSize = 1,
                                   transparency = 1,
                                   defaultTheme = TRUE,
                                   title = NULL,
                                   titleSize = NULL) {
  mat <- getBiomarker(
    inSCE = inSCE,
    useAssay = useAssay,
    gene = feature,
    binary = "Continuous"
  )
  counts <- mat[, 2]

  if (!is.null(groupby)) {
    if (length(groupby) > 1) {
      if (length(groupby) != length(counts)) {
        stop("The input vector for 'groupby' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupby %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupby), "' is not found in ColData.")
      }
      groupby <- as.character(SummarizedExperiment::colData(inSCE)[, groupby])
    }
  }


  p <- .ggViolin(
    y = counts,
    groupby = groupby,
    violin = violin,
    boxplot = boxplot,
    dots = dots,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    dotSize = dotSize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize
  )

  return(p)
}

#' @title Violin plot of any data stored in the SingleCellExperiment object.
#' @description Visualizes values stored in any slot of a
#'  SingleCellExperiment object via a violin plot.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results. Required
#' @param slot Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata"
#' @param annotation Desired vector within the slot used for plotting.
#' @param feature name of feature stored in assay of SingleCellExperiment
#'  object.
#'  Will be used only if "assays" slot is chosen. Deafult NULL.
#' @param groupby Groupings for each numeric value. A user may input a vector
#' equal length to the number of the samples in the SingleCellExperiment
#' object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param xlab label for x-axis. Default NULL.
#' @param ylab label for y-axis. Default NULL.
#' @param axisSize size of x/y-axis labels. Default 10.
#' @param dotSize size of dots. Default 1.
#' @param transparency transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param title title of plot. Default NULL.
#' @param titleSize size of title of plot. Default 15.

#' @examples
#' plotSCEViolin(
#'   inSCE = mouseBrainSubsetSCE, slot = "assays",
#'   annotation = "counts", feature = "Sox2", groupby = "sex"
#' )
#' @export
plotSCEViolin <- function(inSCE,
                          slot,
                          annotation,
                          feature,
                          groupby = NULL,
                          violin = TRUE,
                          boxplot = TRUE,
                          dots = TRUE,
                          xlab = NULL,
                          ylab = NULL,
                          axisSize = 10,
                          dotSize = 1,
                          transparency = 1,
                          defaultTheme = TRUE,
                          title = NULL,
                          titleSize = NULL) {
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

  if (!is.null(groupby)) {
    if (length(groupby) > 1) {
      if (length(groupby) != length(counts)) {
        stop("The input vector for 'groupby' needs to be the same
                     length as the number of samples in your
                     SingleCellExperiment object.")
      }
    } else {
      if (!groupby %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("'", paste(groupby), "' is not found in ColData.")
      }
      groupby <- as.character(SummarizedExperiment::colData(inSCE)[, groupby])
    }
  }


  p <- .ggViolin(
    y = counts,
    groupby = groupby,
    violin = violin,
    boxplot = boxplot,
    dots = dots,
    xlab = xlab,
    ylab = ylab,
    axisSize = axisSize,
    dotSize = dotSize,
    transparency = transparency,
    defaultTheme = defaultTheme,
    title = title,
    titleSize = titleSize
  )

  return(p)
}