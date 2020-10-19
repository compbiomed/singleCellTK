#' @title Plots for runPerCellQC outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runPerCellQC function stored in the colData slot of the SingleCellExperiment
#'  object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' runPerCellQC. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default FALSE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' \dontrun{
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runPerCellQC(sce)
#' plotRunPerCellQCResults(inSCE=sce)
#' }
#' @export
plotRunPerCellQCResults <- function(inSCE,
                                    sample=NULL,
                                    groupBy=NULL,
                                    combinePlot="all",
                                    violin=TRUE,
                                    boxplot=FALSE,
                                    dots=TRUE,
                                    dotSize=1,
                                    axisSize=15,
                                    axisLabelSize=18,
                                    transparency=1,
                                    defaultTheme=TRUE,
                                    titleSize=18,
                                    relHeights=c(1.5, 1.5, 1, 1),
                                    relWidths=c(1, 1, 1, 1),
                                    plotNCols = NULL,
                                    plotNRows = NULL,
                                    plotLabels = "none",
                                    plotLabelSize = 20,
                                    plotLabelPositionX = NULL,
                                    plotLabelPositionY = NULL,
                                    samplePerColumn = TRUE,
                                    sampleRelHeights = 1,
                                    sampleRelWidths = 1) {
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
  sampleVector <- sample

  samples <- unique(sample)

  if(combinePlot == "sample" && length(samples) == 1){
    warning("'combinePlot' was set to 'sample' but the sample was not set,
            or there is only one type of sample specified.")
    combinePlot = "all"
  }

  if (length(samples) > 1) {
    combined.sum <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="sum",
      groupBy=sampleVector,
      xlab="",
      ylab="Counts",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Total counts per cell",
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      titleSize=titleSize,
      combinePlot = "all",
      plotLabels = "none"
    )

    combined.detected <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="detected",
      groupBy=sampleVector,
      xlab="",
      ylab="Features",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Total features detected per cell",
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      titleSize=titleSize,
      combinePlot = "all",
      plotLabels = "none"
    )
    combined.toppercent <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="percent_top_50",
      groupBy=sampleVector,
      xlab="",
      ylab="Gene expression percentage (%)",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Top 50 gene expression percentage",
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      titleSize=titleSize,
      combinePlot = "all",
      plotLabels = "none"
    )

    merged.plots <- list(combined.sum, combined.detected, combined.toppercent)
    names(merged.plots) <- c("Sum", "Detected", "TopPercent")

    if (any(grepl(
      pattern="subsets_",
      names(colData(inSCE))
    ))) {
      subsets <- grep(
        pattern="subsets_",
        names(colData(inSCE)), value=TRUE
      )

      combined.subset <- lapply(subsets, function(x) {
        plotSCEViolinColData(
          inSCE=inSCE,
          coldata=x,
          groupBy = sampleVector,
          xlab="",
          ylab=x,
          violin=violin,
          boxplot=boxplot,
          dots=dots,
          transparency=transparency,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          title=paste0(x, " per cell"),
          dotSize=dotSize,
          titleSize=titleSize,
          summary="median",
          combinePlot = "all",
          plotLabels = "none"
        )
      })
      names(combined.subset) <- c("Gene_Subset_Sum",
                                  "Gene_Subset_Features",
                                  "Gene_Subset_Top50_Percent")
      merged.plots <- c(merged.plots, combined.subset)
    } else {
      combined.subset <- NULL
    }

    merged.plots <- list(Violin = merged.plots)
  }

  res.list <- c()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    if(combinePlot == "sample" | combinePlot == "none" | length(samples) == 1){
      violin.sum <- list(sum = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata="sum",
        sample=sampleSub,
        xlab="",
        ylab="Counts",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title="Total counts per cell",
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary="median",
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.sum)

      violin.detected <- list(detected = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata="detected",
        sample=sampleSub,
        xlab="",
        ylab="Features",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title="Total features detected per cell",
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary="median",
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.detected)

      violin.toppercent <- list(toppercent = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata="percent_top_50",
        sample=sampleSub,
        xlab="",
        ylab="Gene expression percentage (%)",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title="Top 50 gene expression percentage",
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary="median",
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.toppercent)
      names(res.list) <- c("Sum", "Detected", "TopPercent")

      if (any(grepl(
        pattern="subsets_",
        names(colData(inSCESub))
      ))) {
        subsets <- grep(
          pattern="subsets_",
          names(colData(inSCESub)), value=TRUE
        )

        violin.subset <- lapply(subsets, function(x) {
          plotSCEViolinColData(
            inSCE=inSCESub,
            coldata=x,
            sample=sampleSub,
            xlab="",
            ylab=x,
            groupBy=groupBy,
            violin=violin,
            boxplot=boxplot,
            dots=dots,
            transparency=transparency,
            axisSize=axisSize,
            axisLabelSize=axisLabelSize,
            title=paste0(x, " per cell"),
            dotSize=dotSize,
            titleSize=titleSize,
            summary="median",
            combinePlot="all"
          )
        })
        names(violin.subset) <- subsets
      } else {
        violin.subset <- NULL
      }

      if (!is.null(violin.subset)) {
        res.list <- c(res.list, violin.subset)
      }
    }

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    if(combinePlot == "all"){
      plotlist <- c(merged.plots)
    }else if(combinePlot == "sample"){
      plotlist <- c(merged.plots, list(Sample = plotlist))
    }else if(combinePlot == "none"){
      plotlist <- c(merged.plots, list(Sample = plotlist))
    }
  } else {
    plotlist <- unlist(plotlist, recursive=FALSE)
  }

  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      nrows = plotNRows,
                                      ncols = plotNCols,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runEmptyDrops outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runEmptyDrops function stored in the colData slot of the SingleCellExperiment
#'  object via plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runScrublet}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param fdrCutoff Numeric. Thresholds barcodes based on the FDR values from
#'  runEmptyDrops as "Empty Droplet" or "Putative Cell". Default 0.01.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- runEmptyDrops(inSCE=sce)
#' plotEmptyDropsResults(inSCE=sce)
#' @export
plotEmptyDropsResults <- function(inSCE,
                                  sample=NULL,
                                  combinePlot="all",
                                  fdrCutoff=0.01,
                                  defaultTheme=TRUE,
                                  dotSize=1,
                                  titleSize=18,
                                  axisLabelSize=18,
                                  axisSize=15,
                                  legendSize=15,
                                  legendTitleSize=16,
                                  relHeights=1,
                                  relWidths=1,
                                  samplePerColumn = TRUE,
                                  sampleRelHeights = 1,
                                  sampleRelWidths = 1) {
  scatterEmptyDrops <- plotEmptyDropsScatter(inSCE,
    sample=sample,
    combinePlot=combinePlot,
    fdrCutoff=fdrCutoff,
    dotSize=dotSize,
    title="EmptyDrops, Total UMI counts vs Log-Probability",
    titleSize=titleSize,
    defaultTheme=TRUE,
    xlab="Total UMI count",
    ylab="-Log Probability",
    axisLabelSize=axisLabelSize,
    axisSize=axisSize,
    legendTitle=paste0("Cutoff:\nFDR < ", fdrCutoff),
    legendTitleSize=legendTitleSize,
    legendSize=legendSize,
    relHeights = relHeights,
    relWidths = relWidths,
    samplePerColumn = samplePerColumn,
    sampleRelHeights = sampleRelHeights,
    sampleRelWidths = sampleRelWidths
  )

  res.list <- list(scatterEmptyDrops)
  names(res.list) <- c("scatterEmptyDrops")
  return(res.list)
}

#' @title Plots for runEmptyDrops outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runEmptyDrops function stored in the colData slot of the SingleCellExperiment
#'  object via plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runBarcodeRankDrops}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- runBarcodeRankDrops(inSCE=sce)
#' plotBarcodeRankDropsResults(inSCE=sce)
#' @export
plotBarcodeRankDropsResults <- function(inSCE,
                                  sample=NULL,
                                  defaultTheme=TRUE,
                                  dotSize=1,
                                  titleSize=18,
                                  axisLabelSize=18,
                                  axisSize=15,
                                  legendSize=15) {
    scatterBarcodeRank <- plotBarcodeRankScatter(inSCE,
                                               sample=sample,
                                               dotSize=dotSize,
                                               title="BarcodeRanks Rank Plot",
                                               titleSize=titleSize,
                                               defaultTheme=TRUE,
                                               axisLabelSize=axisLabelSize,
                                               axisSize=axisSize,
                                               legendSize=legendSize
    )

    res.list <- list(scatterBarcodeRank)
    names(res.list) <- c("scatterBarcodeRank")
    return(res.list)
}

#' @title Plots for runScrublet outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runScrublet function stored in the colData slot of the SingleCellExperiment
#'  object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runScrublet}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runScrublet(sce)
#' plotScrubletResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotScrubletResults <- function(inSCE,
                                sample=NULL,
                                shape=NULL,
                                groupBy=NULL,
                                combinePlot="all",
                                violin=TRUE,
                                boxplot=FALSE,
                                dots=TRUE,
                                reducedDimName,
                                xlab=NULL,
                                ylab=NULL,
                                dim1=NULL,
                                dim2=NULL,
                                bin=NULL,
                                binLabel=NULL,
                                defaultTheme=TRUE,
                                dotSize=1,
                                transparency=1,
                                titleSize=18,
                                axisLabelSize=18,
                                axisSize=15,
                                legendSize=15,
                                legendTitleSize=16,
                                relHeights=c(1.5, 1, 1),
                                relWidths=c(1, 1, 1),
                                plotNCols = NULL,
                                plotNRows = NULL,
                                plotLabels = "default",
                                plotLabelSize = 20,
                                plotLabelPositionX = NULL,
                                plotLabelPositionY = NULL,
                                samplePerColumn = TRUE,
                                sampleRelHeights = 1,
                                sampleRelWidths = 1) {
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
  sampleVector <- sample
  samples <- unique(sample)

  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scrublet_score",
      groupBy=sampleVector,
      xlab="",
      ylab="Doublet Score",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Scrublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scrublet_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        axisSize=15,
        axisLabelSize=18,
        defaultTheme=defaultTheme,
        cutoff=0.5,
        title="Density, Scrublet Score",
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list <- c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scrublet_score",
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="Scrublet Doublet Score",
      titleSize=titleSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendSize=legendSize,
      legendTitleSize=legendTitleSize,
      combinePlot="all"
    ))
    res.list <- c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
      violinScore <- list(violin_doubletScore = plotSCEViolinColData(
      inSCE=inSCESub, coldata="scrublet_score",
      sample=sampleSub,
      xlab="",
      ylab="Doublet Score",
      groupBy=groupBy,
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="Scrublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize, axisLabelSize=axisLabelSize,
      summary="median",
      combinePlot="all"
    ))
    res.list <- c(res.list, violinScore)
    }

    scatterCall <- list(scatter_doubletCall = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scrublet_call",
      conditionClass="factor",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      transparency=transparency,
      colorScale = c("lightgray","red"),
      defaultTheme=defaultTheme,
      title="Scrublet Doublet Assignment",
      titleSize=titleSize,
      axisSize=axisSize, axisLabelSize=axisLabelSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nAssignment",
      legendTitleSize=legendTitleSize,
      legendSize=legendSize,
      combinePlot="all"
    ))
    res.list <- c(res.list, scatterCall)

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }

  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }

  return(plotlist)
}

#' @title Plots for runDoubletFinder outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDoubletFinder function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runDoubletFinder}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDoubletFinder(sce)
#' plotDoubletFinderResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDoubletFinderResults <- function(inSCE,
                                     sample=NULL,
                                     shape=NULL,
                                     groupBy=NULL,
                                     combinePlot="all",
                                     violin=TRUE,
                                     boxplot=FALSE,
                                     dots=TRUE,
                                     reducedDimName=NULL,
                                     xlab=NULL,
                                     ylab=NULL,
                                     dim1=NULL,
                                     dim2=NULL,
                                     bin=NULL,
                                     binLabel=NULL,
                                     defaultTheme=TRUE,
                                     dotSize=1,
                                     transparency=1,
                                     titleSize=18,
                                     axisLabelSize=18,
                                     axisSize=15,
                                     legendSize=15,
                                     legendTitleSize=16,
                                     relHeights=c(1.5, 1, 1),
                                     relWidths=c(1, 1, 1),
                                     plotNCols = NULL,
                                     plotNRows = NULL,
                                     plotLabels = "default",
                                     plotLabelSize = 20,
                                     plotLabelPositionX = NULL,
                                     plotLabelPositionY = NULL,
                                     samplePerColumn = TRUE,
                                     sampleRelHeights = 1,
                                     sampleRelWidths = 1) {
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
  sampleVector <- sample
  samples <- unique(sample)
  df.scores <- grep(
    pattern="doubletFinder_doublet_score_resolution_",
    names(colData(inSCE)), value=TRUE
  )

  df.labels <- grep(
    pattern="doubletFinder_doublet_label_resolution_",
    names(colData(inSCE)), value=TRUE
  )
  if (length(samples) > 1) {
    merged.plots <- lapply(df.scores, function(x) {
      plotSCEViolinColData(
        inSCE=inSCE,
        coldata=x,
        sample=NULL,
        xlab="",
        ylab="Doublet Score",
        groupBy=sampleVector,
        violin=violin,
        boxplot=boxplot,
        dots=TRUE,
        transparency=transparency,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        titleSize=titleSize,
        title=paste(
          "DoubletFinder Score Resolution",
          gsub(
            pattern="doubletFinder_doublet_score_Resolution_",
            "", x
          )
        ),
        dotSize=dotSize,
        summary="median",
        combinePlot = "all",
        plotLabels = "none"
      )
    })

    names(merged.plots) <- vapply(df.scores, function(x) {
      paste0("Violin_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=x
      ))
    }, character(1))

    # merged.plots <- list(merged.plots)
    merged.plots <- list(Violin = merged.plots)
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    densityScore <- lapply(df.scores, function(x) {
        plotSCEDensityColData(
            inSCE=inSCESub,
            sample=sampleSub,
            coldata=x,
            groupBy=groupBy,
            xlab="Score",
            ylab="Density",
            axisSize=axisSize, axisLabelSize=axisLabelSize,
            defaultTheme=defaultTheme,
            cutoff=0.5,
            combinePlot="all",
            titleSize=titleSize,
            title=paste(
                "Density, Doublet Score Resolution",
                gsub(
                    pattern="doubletFinder_doublet_score_resolution_",
                    "", x
                )
            )
        )
    })
    names(densityScore) <- vapply(df.scores, function(x) {
        paste0("Density_", gsub(
            pattern="doubletFinder_doublet_score_",
            "", x=x
        ))
    }, character(1))
    res.list <- c(res.list, densityScore)

    scatterScore <- lapply(df.scores, function(x) {
      plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        conditionClass="numeric",
        shape=shape,
        colorBy=x,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        defaultTheme=defaultTheme,
        title=paste(
          "Doublet Score Resolution",
          gsub(
            pattern="doubletFinder_doublet_score_resolution_",
            "", x
          )
        ),
        titleSize=titleSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nScore",
        legendSize=legendSize,
        legendTitleSize=legendTitleSize,
        combinePlot="all"
      )
    })

    names(scatterScore) <- vapply(df.scores, function(x) {
      paste0("Scatter_Score_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=x
      ))
    }, character(1))
    res.list <- c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
    violinScore <- lapply(df.scores, function(x) {
      plotSCEViolinColData(
        inSCE=inSCESub,
        coldata=x,
        sample=sampleSub,
        xlab="",
        ylab="Doublet Score",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        defaultTheme=defaultTheme,
        summary="median",
        title=paste(
          "Doublet Score Resolution",
          gsub(
            pattern="doubletFinder_doublet_score_resolution_",
            "", x
          )
        ),
        titleSize=titleSize,
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        combinePlot="all"
      )
    })

    names(violinScore) <- vapply(df.scores, function(x) {
      paste0("violin_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=x
      ))
    }, character(1))
    res.list <- c(res.list, violinScore)
    }


    scatterCall <- lapply(df.labels, function(x) {
      plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        conditionClass="factor",
        shape=shape,
        colorBy=x,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        colorScale = c("red","lightgray"),
        defaultTheme=defaultTheme,
        title=paste(
          "Doublet Call Resolution",
          gsub(
            pattern="doubletFinder_doublet_label_resolution_",
            "", x
          )
        ),
        titleSize=titleSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet Score",
        legendSize=legendSize,
        legendTitleSize=legendTitleSize,
        combinePlot="all"
      )
    })

    names(scatterCall) <- vapply(df.labels, function(x) {
      paste0("Scatter_Call_", gsub(
        pattern="doubletFinder_doublet_label_",
        "", x=x
      ))
    }, character(1))
    res.list <- c(res.list, scatterCall)
    return(res.list)
  })
  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runDoubletCells outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDoubletCells function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runDoubletCells}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param logScore Boolean. If TRUE, the log normalized doublet score will be used.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDoubletCells(sce)
#' plotDoubletCellsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDoubletCellsResults <- function(inSCE,
                                    sample=NULL,
                                    shape=NULL,
                                    groupBy=NULL,
                                    combinePlot="all",
                                    violin=TRUE,
                                    boxplot=FALSE,
                                    dots=TRUE,
                                    logScore=TRUE,
                                    reducedDimName=NULL,
                                    xlab=NULL,
                                    ylab=NULL,
                                    dim1=NULL,
                                    dim2=NULL,
                                    bin=NULL,
                                    binLabel=NULL,
                                    defaultTheme=TRUE,
                                    dotSize=1,
                                    transparency=1,
                                    titleSize=18,
                                    axisLabelSize=18,
                                    axisSize=15,
                                    legendSize=15,
                                    legendTitleSize=16,
                                    relHeights=c(1.5, 1, 1),
                                    relWidths=c(1, 1, 1),
                                    plotNCols = NULL,
                                    plotNRows = NULL,
                                    plotLabels = "default",
                                    plotLabelSize = 20,
                                    plotLabelPositionX = NULL,
                                    plotLabelPositionY = NULL,
                                    samplePerColumn = TRUE,
                                    sampleRelHeights = 1,
                                    sampleRelWidths = 1) {
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
  sampleVector <- sample

  if (logScore) {
    coldata = "scran_doubletCells_score_log10"
    titleDoubletCells <- "DoubletCells Doublet Score, log10"
  } else {
    coldata = "scran_doubletCells_score"
    titleDoubletCells <- "DoubletCells Doublet Score"
  }

  samples <- unique(sample)
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata=coldata,
      groupBy=sampleVector,
      xlab="",
      ylab="Doublet Score",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      title=titleDoubletCells,
      titleSize=titleSize,
      dotSize=dotSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata=coldata,
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=paste0("Density, ", titleDoubletCells),
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy=coldata,
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      colorLow="gray",
      colorHigh="blue",
      transparency=transparency,
      defaultTheme=defaultTheme,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      title=titleDoubletCells,
      titleSize=titleSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendSize=legendSize,
      legendTitleSize=legendTitleSize,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
    violinScore <- list(violin_doubletScore = plotSCEViolinColData(
      inSCE=inSCESub,
      coldata=coldata,
      sample=sampleSub,
      xlab="",
      ylab="Doublet Score",
      groupBy=groupBy,
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title=titleDoubletCells,
      titleSize=titleSize,
      defaultTheme=defaultTheme,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      dotSize=dotSize,
      summary="median",
      combinePlot="all"
    ))
    res.list = c(res.list, violinScore)
    }

    return(res.list)
  })
  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runCxds outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runCxds function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runCxds}.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runCxds(sce)
#' plotCxdsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotCxdsResults <- function(inSCE,
                            sample=NULL,
                            shape=NULL,
                            groupBy=NULL,
                            combinePlot="all",
                            violin=TRUE,
                            boxplot=FALSE,
                            dots=TRUE,
                            reducedDimName=NULL,
                            xlab=NULL,
                            ylab=NULL,
                            dim1=NULL,
                            dim2=NULL,
                            bin=NULL,
                            binLabel=NULL,
                            defaultTheme=TRUE,
                            dotSize=1,
                            transparency=1,
                            titleSize=18,
                            axisLabelSize=18,
                            axisSize=15,
                            legendSize=15,
                            legendTitleSize=16,
                            relHeights=c(1.5, 1, 1),
                            relWidths=c(1, 1, 1),
                            plotNCols = NULL,
                            plotNRows = NULL,
                            plotLabels = "default",
                            plotLabelSize = 20,
                            plotLabelPositionX = NULL,
                            plotLabelPositionY = NULL,
                            samplePerColumn = TRUE,
                            sampleRelHeights = 1,
                            sampleRelWidths = 1) {
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
  sampleVector <- sample
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_cxds_score",
      groupBy=sampleVector,
      xlab="",
      ylab="Doublet Score",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="CXDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scds_cxds_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, CXDS Score",
        titleSize=titleSize,
        plotLabels = NULL,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scds_cxds_score",
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="CXDS Doublet Score",
      titleSize=titleSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendSize=legendSize,
      legendTitleSize=legendTitleSize,
      plotLabels = NULL,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
    violinScore <- list(violin_doubletScore = plotSCEViolinColData(
      inSCE=inSCESub,
      coldata="scds_cxds_score",
      sample=sampleSub,
      xlab="",
      ylab="Doublet Score",
      groupBy=groupBy,
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="CXDS Doublet Score",
      titleSize=titleSize,
      defaultTheme=defaultTheme,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      summary="median",
      combinePlot="all"
    ))
    res.list = c(res.list, violinScore)
  }

    if("scds_cxds_call" %in% names(SingleCellExperiment::colData(inSCE))){
      scatterCall <- list(scatter_doubletCall = plotSCEDimReduceColData(
          inSCE=inSCESub,
          sample=sampleSub,
          colorBy="scds_cxds_call",
          conditionClass="factor",
          shape=shape,
          reducedDimName=reducedDimName,
          xlab=xlab,
          ylab=ylab,
          dim1=dim1,
          dim2=dim2,
          bin=bin,
          binLabel=binLabel,
          dotSize=dotSize,
          transparency=transparency,
          colorScale = c("lightgray","red"),
          defaultTheme=defaultTheme,
          title="CXDS Doublet Assignment",
          titleSize=titleSize,
          axisSize=axisSize, axisLabelSize=axisLabelSize,
          labelClusters=FALSE,
          legendTitle="Doublet \nAssignment",
          legendTitleSize=legendTitleSize,
          legendSize=legendSize,
          combinePlot="all"
      ))
      res.list <- c(res.list, scatterCall)
    }
    return(res.list)
  })
  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runBcds outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runBcds function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runBcds}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runBcds(sce)
#' plotBcdsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotBcdsResults <- function(inSCE,
                            sample=NULL,
                            shape=NULL,
                            groupBy=NULL,
                            combinePlot="all",
                            violin=TRUE,
                            boxplot=FALSE,
                            dots=TRUE,
                            reducedDimName=NULL,
                            xlab=NULL,
                            ylab=NULL,
                            dim1=NULL,
                            dim2=NULL,
                            bin=NULL,
                            binLabel=NULL,
                            defaultTheme=TRUE,
                            dotSize=1,
                            transparency=1,
                            titleSize=18,
                            axisLabelSize=18,
                            axisSize=15,
                            legendSize=15,
                            legendTitleSize=16,
                            relHeights=c(1.5, 1, 1),
                            relWidths=c(1, 1, 1),
                            plotNCols = NULL,
                            plotNRows = NULL,
                            plotLabels = "default",
                            plotLabelSize = 20,
                            plotLabelPositionX = NULL,
                            plotLabelPositionY = NULL,
                            samplePerColumn = TRUE,
                            sampleRelHeights = 1,
                            sampleRelWidths = 1) {
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
  sampleVector <- sample
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_bcds_score",
      groupBy=sampleVector,
      xlab="",
      ylab="Doublet Score",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="BCDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    densityScore <- list(density_doubletScore = plotSCEDensityColData(
      inSCE=inSCESub,
      sample=sampleSub,
      coldata="scds_bcds_score",
      groupBy=groupBy,
      xlab="Score",
      ylab="Density",
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      defaultTheme=defaultTheme,
      title="Density, BCDS Score",
      titleSize=titleSize,
      plotLabels = NULL,
      combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scds_bcds_score",
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="BCDS Doublet Score",
      titleSize=titleSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendSize=legendSize,
      legendTitleSize=legendTitleSize,
      plotLabels = NULL,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
      violinScore <- list(violin_doubletScore = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata="scds_bcds_score",
        sample=sampleSub,
        xlab="",
        ylab="Doublet Score",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title="BCDS Doublet Score",
        titleSize=titleSize,
        defaultTheme=defaultTheme,
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary="median",
        combinePlot="all"
      ))
      res.list = c(res.list, violinScore)
    }

    if("scds_bcds_call" %in% names(SingleCellExperiment::colData(inSCE))){
      scatterCall <- list(scatter_doubletCall = plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        colorBy="scds_bcds_call",
        conditionClass="factor",
        shape=shape,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title="BCDS Doublet Assignment",
        titleSize=titleSize,
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendTitleSize=legendTitleSize,
        legendSize=legendSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, scatterCall)
    }
    return(res.list)
  })
  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
    plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}


#' @title Plots for runCxdsBcdsHybrid outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runCxdsBcdsHybrid function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runCxdsBcdsHybrid}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runCxdsBcdsHybrid(sce)
#' plotScdsHybridResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotScdsHybridResults <- function(inSCE,
                                  sample=NULL,
                                  shape=NULL,
                                  groupBy=NULL,
                                  combinePlot="all",
                                  violin=TRUE,
                                  boxplot=FALSE,
                                  dots=TRUE,
                                  reducedDimName=NULL,
                                  xlab=NULL,
                                  ylab=NULL,
                                  dim1=NULL,
                                  dim2=NULL,
                                  bin=NULL,
                                  binLabel=NULL,
                                  defaultTheme=TRUE,
                                  dotSize=1,
                                  transparency=1,
                                  titleSize=18,
                                  axisLabelSize=18,
                                  axisSize=15,
                                  legendSize=15,
                                  legendTitleSize=16,
                                  relHeights=c(1.5, 1, 1),
                                  relWidths=c(1, 1, 1),
                                  plotNCols = NULL,
                                  plotNRows = NULL,
                                  plotLabels = "default",
                                  plotLabelSize = 20,
                                  plotLabelPositionX = NULL,
                                  plotLabelPositionY = NULL,
                                  samplePerColumn = TRUE,
                                  sampleRelHeights = 1,
                                  sampleRelWidths = 1) {
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
  sampleVector <- sample
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_hybrid_score",
      groupBy=sampleVector,
      xlab="",
      ylab="Doublet Score",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="CXDS BCDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list <- c()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scds_hybrid_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, CXDS BCDS Hybrid Score",
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scds_hybrid_score",
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="CXDS BCDS Doublet Score",
      titleSize=titleSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendTitleSize=legendTitleSize,
      legendSize=legendSize,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){
    violinScore <- list(violin_doubletScore = plotSCEViolinColData(
      inSCE=inSCESub,
      coldata="scds_hybrid_score",
      sample=sampleSub,
      xlab="",
      ylab="Doublet Score",
      groupBy=groupBy,
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="CXDS BCDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      summary="median",
      combinePlot="all"
    ))
    res.list = c(res.list, violinScore)
    }

    if("scds_hybrid_call" %in% names(SingleCellExperiment::colData(inSCE))){
      scatterCall <- list(scatter_doubletCall = plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        colorBy="scds_hybrid_call",
        conditionClass="factor",
        shape=shape,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title="CXDS BCDS Doublet Assignment",
        titleSize=titleSize,
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendTitleSize=legendTitleSize,
        legendSize=legendSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, scatterCall)
    }

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runDecontX outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDecontX function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runDecontX}. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupBy Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the
#' \linkS4class{SingleCellExperiment} object. Required.
#' @param xlab Character vector. Label for x-axis. Default NULL.
#' @param ylab Character vector. Label for y-axis. Default NULL.
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
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param titleSize Size of title of plot. Default 18.
#' @param axisSize Size of x/y-axis ticks. Default 15.
#' @param axisLabelSize Size of x/y-axis labels. Default 18.
#' @param legendSize size of legend. Default 15.
#' @param legendTitleSize size of legend title. Default 16.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param plotLabels labels to each plot. If set to "default", will use the name of the samples
#'  as the labels. If set to "none", no label will be plotted.
#' @param plotLabelSize size of labels
#' @param plotLabelPositionX Numeric vector. The X position of the plot label.
#' @param plotLabelPositionY Numeric vector. The Y position of the plot label.
#' @param samplePerColumn If TRUE, when there are multiple samples and combining by "all",
#'  the output .ggplot will have plots from each sample on a single column. Default TRUE.
#' @param sampleRelHeights If there are multiple samples and combining by "all",
#'  the relative heights for each plot.
#' @param sampleRelWidths If there are multiple samples and combining by "all",
#'  the relative widths for each plot.
#' @return list of .ggplot objects
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDecontX(sce)
#' plotDecontXResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDecontXResults <- function(inSCE,
                               sample=NULL,
                               shape=NULL,
                               groupBy=NULL,
                               combinePlot="all",
                               violin=TRUE,
                               boxplot=FALSE,
                               dots=TRUE,
                               reducedDimName=NULL,
                               xlab=NULL,
                               ylab=NULL,
                               dim1=NULL,
                               dim2=NULL,
                               bin=NULL,
                               binLabel=NULL,
                               defaultTheme=TRUE,
                               dotSize=1,
                               transparency=1,
                               titleSize=18,
                               axisLabelSize=18,
                               axisSize=15,
                               legendSize=15,
                               legendTitleSize=16,
                               relHeights=c(1.5, 1, 1),
                               relWidths=c(1, 1, 1),
                               plotNCols = NULL,
                               plotNRows = NULL,
                               plotLabels = "default",
                               plotLabelSize = 20,
                               plotLabelPositionX = NULL,
                               plotLabelPositionY = NULL,
                               samplePerColumn = TRUE,
                               sampleRelHeights = 1,
                               sampleRelWidths = 1) {
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
  sampleVector <- sample
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata="decontX_contamination",
      groupBy=sampleVector,
      xlab="",
      ylab="DecontX Contamination",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      transparency=transparency,
      title="DecontX Contamination Score",
      titleSize=titleSize,
      dotSize=dotSize,
      gridLine=TRUE,
      summary="median",
      combinePlot = "all",
      plotLabels = "none"
    ))
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Violin"
  }

  res.list = list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]
    densityContamination <- list(density_decontXContamination = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="decontX_contamination",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, DecontX Contamination Score",
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityContamination)

    scatterContamination <- list(scatter_decontXContamination = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="decontX_contamination",
      conditionClass="numeric",
      shape=shape,
      reducedDimName=reducedDimName,
      xlab=xlab,
      ylab=ylab,
      dim1=dim1,
      dim2=dim2,
      bin=bin,
      binLabel=binLabel,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      transparency=transparency,
      defaultTheme=defaultTheme,
      title="DecontX Contamination",
      titleSize=titleSize,
      labelClusters=FALSE,
      legendTitle="Contamination",
      legendTitleSize=legendTitleSize,
      legendSize=legendSize,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterContamination)

    if(combinePlot != "all" | length(samples) == 1){
    violinContamination <- list(violin_decontXContamination = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata="decontX_contamination",
        sample=sampleSub,
        xlab="", ylab="DecontX Contamination",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title="DecontX Contamination Score",
        titleSize=titleSize,
        defaultTheme=defaultTheme,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        dotSize=dotSize,
        summary="median",
        combinePlot="all"
    ))
    res.list = c(res.list, violinContamination)
    }


    scatterCluster <- list(scatter_decontXClusters = plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        colorBy="decontX_clusters",
        conditionClass="factor",
        shape=shape,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        dotSize=dotSize,
        transparency=transparency,
        defaultTheme=defaultTheme,
        title="DecontX Clusters",
        titleSize=titleSize,
        labelClusters=TRUE,
        legendTitle="Clusters",
        legendSize=legendSize,
        legendTitleSize=legendTitleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, scatterCluster)

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
      names(plotlist) <- samples
      plotlist <- c(merged.plots, list(Sample = plotlist))
  } else {
      plotlist <- unlist(plotlist, recursive=FALSE)
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      labels = plotLabels,
                                      labelSize = plotLabelSize,
                                      labelPositionX = plotLabelPositionX,
                                      labelPositionY = plotLabelPositionY,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}
