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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default "median".
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param baseSize The base font size for all text. Default 15.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single
#' .ggplot object, while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
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
                                    dotSize=0.5,
                                    summary="median",
                                    summaryTextSize=3,
                                    baseSize=15,
                                    axisSize=NULL,
                                    axisLabelSize=NULL,
                                    transparency=1,
                                    defaultTheme=TRUE,
                                    titleSize=NULL,
                                    relHeights=1,
                                    relWidths=1,
                                    labelSamples = TRUE,
                                    plotNCols = NULL,
                                    plotNRows = NULL,
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
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
      titleSize=titleSize,
      combinePlot = "all",
      plotLabels = "none"
    )
    combined.toppercent <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="percent.top_50",
      groupBy=sampleVector,
      xlab="",
      ylab="Gene expression percentage (%)",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Top 50 gene expression percentage",
      dotSize=dotSize,
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
      titleSize=titleSize,
      combinePlot = "all",
      plotLabels = "none"
    )

    merged.plots <- list(combined.sum, combined.detected, combined.toppercent)
    names(merged.plots) <- c("Sum", "Detected", "TopPercent")

    if (any(grepl(pattern="subsets_",names(colData(inSCE))
    ) | grepl(pattern="mito_", names(colData(inSCE))))) { 
      subsets <- grep(
        pattern="subsets_",
        names(colData(inSCE)), value=TRUE
      )
      mitos <- grep(
        pattern="mito_",
        names(colData(inSCE)), value=TRUE
      )
      subsets <- c(subsets, mitos)

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
          baseSize=baseSize,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          title=paste0(x, " per cell"),
          dotSize=dotSize,
          titleSize=titleSize,
          summary=summary,
          summaryTextSize=summaryTextSize,
          combinePlot = "all",
          plotLabels = "none"
        )
      })
      names(combined.subset) <- subsets
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
      title = "Total counts per cell"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }

      if(combinePlot == "sample" | combinePlot == "all"){
        baseSize = baseSize * 0.5
      }
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
        title=title,
        dotSize=dotSize,
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary=summary,
        summaryTextSize=summaryTextSize,
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.sum)

      title = "Total features detected per cell"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
        title=title,
        dotSize=dotSize,
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary=summary,
        summaryTextSize=summaryTextSize,
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.detected)

      topPattern <- grep(
        pattern="percent.top_50$",
        names(colData(inSCESub)), value=TRUE
      )
      title = "Top 50 gene expression percentage"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
      violin.toppercent <- list(toppercent = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata=topPattern,
        sample=sampleSub,
        xlab="",
        ylab="Gene expression percentage (%)",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        title=title,
        dotSize=dotSize,
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary=summary,
        summaryTextSize=summaryTextSize,
        titleSize=titleSize,
        combinePlot="all"
      ))
      res.list <- c(res.list, violin.toppercent)
      names(res.list) <- c("Sum", "Detected", "TopPercent")

      if (any(grepl(pattern="subsets_", names(colData(inSCESub))) |
              grepl(pattern="mito_", names(colData(inSCESub))))) {
        subsets <- grep(
          pattern="subsets_",
          names(colData(inSCESub)), value=TRUE
        )
        mitos <- grep(
          pattern="mito_",
          names(colData(inSCESub)), value=TRUE
        )
        subsets <- c(subsets, mitos)

        violin.subset <- lapply(subsets, function(y) {
          title = paste0(y, " per cell")
          if(labelSamples && length(samples) > 1){
            title = paste0(title, ", ", x)
          }
          plotSCEViolinColData(
            inSCE=inSCESub,
            coldata=y,
            sample=sampleSub,
            xlab="",
            ylab=y,
            groupBy=groupBy,
            violin=violin,
            boxplot=boxplot,
            dots=dots,
            transparency=transparency,
            baseSize=baseSize,
            axisSize=axisSize,
            axisLabelSize=axisLabelSize,
            title=title,
            dotSize=dotSize,
            titleSize=titleSize,
            summary=summary,
            summaryTextSize=summaryTextSize,
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
    relHeights=1
  }

  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      nrows = plotNRows,
                                      ncols = plotNCols,
                                      labels = "none",
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
#' @param dotSize Size of dots. Default 0.5.
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
                                  dotSize=0.5,
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
#' @param dotSize Size of dots. Default 0.5.
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
                                  dotSize=0.5,
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
#' sce <- runScrublet(sce)
#' plotScrubletResults(inSCE=sce, reducedDimName="UMAP")
#' }
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
                                dotSize=0.5,
                                summary="median",
                                summaryTextSize=3,
                                transparency=1,
                                baseSize=15,
                                titleSize=NULL,
                                axisLabelSize=NULL,
                                axisSize=NULL,
                                legendSize=NULL,
                                legendTitleSize=NULL,
                                relHeights=1,
                                relWidths=c(1, 1, 1),
                                plotNCols = NULL,
                                plotNRows = NULL,
                                labelSamples = TRUE,
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
      baseSize=baseSize,
      title="Scrublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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

    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }

    title = "Density, Scrublet Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scrublet_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list <- c(res.list, densityScore)

    title = "Scrublet Doublet Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
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
      title = "Scrublet Score"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }

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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize, axisLabelSize=axisLabelSize,
      summary=summary,
      summaryTextSize=summaryTextSize,
      combinePlot="all"
    ))
    res.list <- c(res.list, violinScore)
    }

    title = "Scrublet Doublet Assignment"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
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
      baseSize=baseSize,
      colorScale = c("lightgray","red"),
      defaultTheme=defaultTheme,
      title=title,
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
      plotLabels <- "none"
      relHeights=1
  }

  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
                                     reducedDimName="UMAP",
                                     xlab=NULL,
                                     ylab=NULL,
                                     dim1=NULL,
                                     dim2=NULL,
                                     bin=NULL,
                                     binLabel=NULL,
                                     defaultTheme=TRUE,
                                     dotSize=0.5,
                                     summary="median",
                                     summaryTextSize=3,
                                     transparency=1,
                                     baseSize=15,
                                     titleSize=NULL,
                                     axisLabelSize=NULL,
                                     axisSize=NULL,
                                     legendSize=NULL,
                                     legendTitleSize=NULL,
                                     relHeights=1,
                                     relWidths=c(1, 1, 1),
                                     plotNCols = NULL,
                                     plotNRows = NULL,
                                     labelSamples = TRUE,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
  }

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
        baseSize=baseSize,
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
        summary=summary,
        summaryTextSize=summaryTextSize,
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

    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }

    densityScore <- lapply(df.scores, function(y) {
        title <- paste(
          "Density, Doublet Score Resolution",
          gsub(
            pattern="doubletFinder_doublet_score_resolution_","", y))
        if(labelSamples && length(samples) > 1){
          title = paste0(title, ", ", x)
        }
        plotSCEDensityColData(
            inSCE=inSCESub,
            sample=sampleSub,
            coldata=y,
            groupBy=groupBy,
            xlab="Score",
            ylab="Density",
            baseSize=baseSize,
            axisSize=axisSize,
            axisLabelSize=axisLabelSize,
            defaultTheme=defaultTheme,
            cutoff=0.5,
            combinePlot="all",
            titleSize=titleSize,
            title=title
            )
    })
    names(densityScore) <- vapply(df.scores, function(y) {
        paste0("Density_", gsub(
            pattern="doubletFinder_doublet_score_",
            "", x=y
        ))
    }, character(1))
    res.list <- c(res.list, densityScore)

    scatterScore <- lapply(df.scores, function(y) {
      title <- paste(
        "Doublet Score Resolution",
        gsub(
          pattern="doubletFinder_doublet_score_resolution_","", y))
        if(labelSamples && length(samples) > 1){
          title = paste0(title, ", ", x)
        }

      plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        conditionClass="numeric",
        shape=shape,
        colorBy=y,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        baseSize=baseSize,
        defaultTheme=defaultTheme,
        title=title,
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

    names(scatterScore) <- vapply(df.scores, function(y) {
      paste0("Scatter_Score_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=y
      ))
    }, character(1))
    res.list <- c(res.list, scatterScore)

    if(combinePlot != "all" | length(samples) == 1){

    violinScore <- lapply(df.scores, function(y) {
      title <- paste(
        "Doublet Score Resolution",
        gsub(
          pattern="doubletFinder_doublet_score_resolution_",
          "", y))

        if(labelSamples && length(samples) > 1){
          title = paste0(title, ", ", x)
        }

      plotSCEViolinColData(
        inSCE=inSCESub,
        coldata=y,
        sample=sampleSub,
        xlab="",
        ylab="Doublet Score",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        baseSize=baseSize,
        defaultTheme=defaultTheme,
        summary=summary,
        summaryTextSize=summaryTextSize,
        title=title,
        titleSize=titleSize,
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        combinePlot="all"
      )
    })

    names(violinScore) <- vapply(df.scores, function(y) {
      paste0("violin_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=y
      ))
    }, character(1))
    res.list <- c(res.list, violinScore)
    }

    scatterCall <- lapply(df.labels, function(y) {
      title <- paste(
        "Doublet Call Resolution",
        gsub(
          pattern="doubletFinder_doublet_label_resolution_",
          "", y))

      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }

      plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        conditionClass="factor",
        shape=shape,
        colorBy=y,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        dotSize=dotSize,
        transparency=transparency,
        baseSize=baseSize,
        colorScale = c("lightgray", "red"),
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendSize=legendSize,
        legendTitleSize=legendTitleSize,
        combinePlot="all"
      )
    })

    names(scatterCall) <- vapply(df.labels, function(y) {
      paste0("Scatter_Call_", gsub(
        pattern="doubletFinder_doublet_label_",
        "", x=y
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
      plotLabels <- "none"
      relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}

#' @title Plots for runScDblFinder outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runScDblFinder function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' \link{runScDblFinder}. Required.
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
#' sce <- runScDblFinder(sce)
#' plotScDblFinderResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotScDblFinderResults <- function(inSCE,
                                    sample=NULL,
                                    shape=NULL,
                                    groupBy=NULL,
                                    combinePlot="all",
                                    violin=TRUE,
                                    boxplot=FALSE,
                                    dots=TRUE,
                                    reducedDimName="UMAP",
                                    xlab=NULL,
                                    ylab=NULL,
                                    dim1=NULL,
                                    dim2=NULL,
                                    bin=NULL,
                                    binLabel=NULL,
                                    defaultTheme=TRUE,
                                    dotSize=0.5,
                                    summary="median",
                                    summaryTextSize=3,
                                    transparency=1,
                                    baseSize=15,
                                    titleSize=NULL,
                                    axisLabelSize=NULL,
                                    axisSize=NULL,
                                    legendSize=NULL,
                                    legendTitleSize=NULL,
                                    relHeights=1,
                                    relWidths=c(1, 1, 1),
                                    plotNCols = NULL,
                                    plotNRows = NULL,
                                    labelSamples = TRUE,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
  }
  coldata = "scDblFinder_doublet_score"
  titleScDblFinder <- "ScDblFinder Doublet Score"

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
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      title=titleScDblFinder,
      titleSize=titleSize,
      dotSize=dotSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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

    title = paste0("Density, ", titleScDblFinder)
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata=coldata,
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    title = titleScDblFinder
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      title=title,
      titleSize=titleSize,
      labelClusters=FALSE,
      legendTitle="Doublet \nScore",
      legendSize=legendSize,
      legendTitleSize=legendTitleSize,
      combinePlot="all"
    ))
    res.list = c(res.list, scatterScore)

    if("scDblFinder_doublet_call" %in% names(SingleCellExperiment::colData(inSCE))){
      title = "scDblFinder Doublet Assignment"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
      scatterCall <- list(scatter_doubletCall = plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        colorBy="scDblFinder_doublet_call",
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
        baseSize=baseSize,
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title=title,
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

    if(combinePlot != "all" | length(samples) == 1){
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
      baseSize=baseSize,
      title=title,
      titleSize=titleSize,
      defaultTheme=defaultTheme,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      dotSize=dotSize,
      summary=summary,
      summaryTextSize=summaryTextSize,
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
      plotLabels <- "none"
      relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
                            reducedDimName="UMAP",
                            xlab=NULL,
                            ylab=NULL,
                            dim1=NULL,
                            dim2=NULL,
                            bin=NULL,
                            binLabel=NULL,
                            defaultTheme=TRUE,
                            dotSize=0.5,
                            summary="median",
                            summaryTextSize=3,
                            transparency=1,
                            baseSize=15,
                            titleSize=NULL,
                            axisLabelSize=NULL,
                            axisSize=NULL,
                            legendSize=NULL,
                            legendTitleSize=NULL,
                            relHeights=1,
                            relWidths=c(1, 1, 1),
                            plotNCols = NULL,
                            plotNRows = NULL,
                            labelSamples = TRUE,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
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
      baseSize=baseSize,
      title="CXDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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

    title = "Density, CXDS Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scds_cxds_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        plotLabels = NULL,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    title = "CXDS Doublet Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
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
      title = "CXDS Doublet Score"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
      baseSize=baseSize,
      title=title,
      titleSize=titleSize,
      defaultTheme=defaultTheme,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      summary=summary,
      summaryTextSize=summaryTextSize,
      combinePlot="all"
    ))
    res.list = c(res.list, violinScore)
  }

    if("scds_cxds_call" %in% names(SingleCellExperiment::colData(inSCE))){
      title = "CXDS Doublet Assignment"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
          baseSize=baseSize,
          colorScale = c("lightgray","red"),
          defaultTheme=defaultTheme,
          title=title,
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
      plotLabels <- "none"
      relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 15.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
                            reducedDimName="UMAP",
                            xlab=NULL,
                            ylab=NULL,
                            dim1=NULL,
                            dim2=NULL,
                            bin=NULL,
                            binLabel=NULL,
                            defaultTheme=TRUE,
                            dotSize=0.5,
                            summary="median",
                            summaryTextSize=3,
                            transparency=1,
                            baseSize=15,
                            titleSize=NULL,
                            axisLabelSize=NULL,
                            axisSize=NULL,
                            legendSize=NULL,
                            legendTitleSize=NULL,
                            relHeights=1,
                            relWidths=c(1, 1, 1),
                            plotNCols = NULL,
                            plotNRows = NULL,
                            labelSamples = TRUE,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
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
      baseSize=baseSize,
      title="BCDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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

    title = "Density, BCDS Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
      inSCE=inSCESub,
      sample=sampleSub,
      coldata="scds_bcds_score",
      groupBy=groupBy,
      xlab="Score",
      ylab="Density",
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      defaultTheme=defaultTheme,
      title=title,
      titleSize=titleSize,
      plotLabels = NULL,
      combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    title = "BCDS Doublet Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
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
      title = "BCDS Doublet Score"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
        baseSize=baseSize,
        title=title,
        titleSize=titleSize,
        defaultTheme=defaultTheme,
        dotSize=dotSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        summary=summary,
        summaryTextSize=summaryTextSize,
        combinePlot="all"
      ))
      res.list = c(res.list, violinScore)
    }

    if("scds_bcds_call" %in% names(SingleCellExperiment::colData(inSCE))){
      title = "BCDS Doublet Assignment"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
        baseSize=baseSize,
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title=title,
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
    plotLabels <- "none"
    relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
                                  reducedDimName="UMAP",
                                  xlab=NULL,
                                  ylab=NULL,
                                  dim1=NULL,
                                  dim2=NULL,
                                  bin=NULL,
                                  binLabel=NULL,
                                  defaultTheme=TRUE,
                                  dotSize=0.5,
                                  summary="median",
                                  summaryTextSize=3,
                                  transparency=1,
                                  baseSize=15,
                                  titleSize=NULL,
                                  axisLabelSize=NULL,
                                  axisSize=NULL,
                                  legendSize=NULL,
                                  legendTitleSize=NULL,
                                  relHeights=1,
                                  relWidths=c(1, 1, 1),
                                  plotNCols = NULL,
                                  plotNRows = NULL,
                                  labelSamples = TRUE,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
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
      baseSize=baseSize,
      title="CXDS BCDS Doublet Score",
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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
    title = "Density, CXDS BCDS Hybrid Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
    }
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scds_hybrid_score",
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityScore)

    title = "CXDS BCDS Hybrid Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
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
      title = "CXDS BCDS Hybrid Score"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
      baseSize=baseSize,
      defaultTheme=defaultTheme,
      title=title,
      titleSize=titleSize,
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      summary=summary,
      summaryTextSize=summaryTextSize,
      combinePlot="all"
    ))
    res.list = c(res.list, violinScore)
    }

    if("scds_hybrid_call" %in% names(SingleCellExperiment::colData(inSCE))){
      title = "CXDS BCDS Doublet Assignment"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
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
        baseSize=baseSize,
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title=title,
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
      plotLabels <- "none"
      relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
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
#' @param bgResult Boolean. If TRUE, will plot decontX results generated with
#' raw/droplet matrix Default FALSE. 
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
#' \linkS4class{SingleCellExperiment} object. Required. Default = "UMAP"
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
#' @param dotSize Size of dots. Default 0.5.
#' @param summary Adds a summary statistic, as well as a crossbar to the
#'  violin plot. Options are "mean" or "median". Default NULL.
#' @param summaryTextSize The text size of the summary statistic displayed
#'  above the violin plot. Default 3.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param baseSize The base font size for all text. Default 12.
#'  Can be overwritten by titleSize, axisSize, and axisLabelSize,
#'  legendSize, legendTitleSize.
#' @param titleSize Size of title of plot. Default NULL.
#' @param axisSize Size of x/y-axis ticks. Default NULL.
#' @param axisLabelSize Size of x/y-axis labels. Default NULL.
#' @param legendSize size of legend. Default NULL.
#' @param legendTitleSize size of legend title. Default NULL.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#' @param clusterLabelSize Numeric. Determines the size of cluster label
#'  when `labelClusters` is set to TRUE. Default 3.5.
#' @param combinePlot Must be either "all", "sample", or "none". "all" will combine all plots into a single .ggplot object,
#' while "sample" will output a list of plots separated by sample. Default "all".
#' @param relHeights Relative heights of plots when combine is set.
#' @param relWidths Relative widths of plots when combine is set.
#' @param plotNCols Number of columns when plots are combined in a grid.
#' @param plotNRows Number of rows when plots are combined in a grid.
#' @param labelSamples Will label sample name in title of plot if TRUE. Default TRUE.
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
#' sce <- runDecontX(sce)
#' plotDecontXResults(inSCE=sce, reducedDimName="decontX_UMAP")
#' @export
plotDecontXResults <- function(inSCE,
                               sample=NULL,
                               bgResult = FALSE,
                               shape=NULL,
                               groupBy=NULL,
                               combinePlot="all",
                               violin=TRUE,
                               boxplot=FALSE,
                               dots=TRUE,
                               reducedDimName="UMAP",
                               xlab=NULL,
                               ylab=NULL,
                               dim1=NULL,
                               dim2=NULL,
                               bin=NULL,
                               binLabel=NULL,
                               defaultTheme=TRUE,
                               dotSize=0.5,
                               summary="median",
                               summaryTextSize=3,
                               transparency=1,
                               baseSize=15,
                               titleSize=NULL,
                               axisLabelSize=NULL,
                               axisSize=NULL,
                               legendSize=NULL,
                               legendTitleSize=NULL,
                               relHeights=1,
                               relWidths=c(1, 1, 1),
                               plotNCols = NULL,
                               plotNRows = NULL,
                               labelSamples = TRUE,
                               labelClusters = TRUE,
                               clusterLabelSize = 3.5,
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

  if (!(reducedDimName %in% reducedDimNames(inSCE))){
    stop("Specified `reducedDimName` is not found in input
         SingleCellExperiment object. Please check for spelling errors
         with reducedDimNames().")
  }

  scoreCol <- "decontX_contamination"
  clusterCol <- "decontX_clusters"

  if (!isTRUE(bgResult) & !scoreCol %in% colnames(SummarizedExperiment::colData(inSCE))) {
      stop("The result of running decontX without raw/droplet matrix 
           is not found in the input SingleCellExperiment object. 
           Please check whether runDecontX has been run without
           'background' parameter. ")    
  }

  if (isTRUE(bgResult)) {
    bgColId <- grep('decontX_contamination_bg', colnames(SummarizedExperiment::colData(inSCE)))

    if (length(bgColId) == 0) {
      stop("The result of running decontX with raw/droplet matrix 
           is not found in the input SingleCellExperiment object. 
           Please check whether runDecontX has been run with
           'background' parameter. ")
    } else {
      scoreCol <- "decontX_contamination_bg"
      clusterCol <- "decontX_clusters_bg"
    }
  }

  samples <- unique(sample)
  sampleVector <- sample
  if (length(samples) > 1) {
    merged.plots <- list(Score = plotSCEViolinColData(
      inSCE=inSCE,
      coldata=scoreCol,
      groupBy=sampleVector,
      xlab="",
      ylab="DecontX Contamination",
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      baseSize=baseSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      transparency=transparency,
      title="DecontX Contamination Score",
      titleSize=titleSize,
      dotSize=dotSize,
      gridLine=TRUE,
      summary=summary,
      summaryTextSize=summaryTextSize,
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
    title = "Density, DecontX Contamination Score"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }
    if(combinePlot == "sample" | combinePlot == "all"){
      baseSize = baseSize * 0.5
      clusterLabelSize <- clusterLabelSize * 0.5
    }
    densityContamination <- list(density_decontXContamination = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata=scoreCol,
        groupBy=groupBy,
        xlab="Score",
        ylab="Density",
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        combinePlot="all"
    ))
    res.list = c(res.list, densityContamination)

    scatterContamination <- list(scatter_decontXContamination = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy=scoreCol,
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
      baseSize=baseSize,
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
      title = "DecontX Contamination Score"
      if(labelSamples && length(samples) > 1){
        title = paste0(title, ", ", x)
      }
    violinContamination <- list(violin_decontXContamination = plotSCEViolinColData(
        inSCE=inSCESub,
        coldata=scoreCol,
        sample=sampleSub,
        xlab="", ylab="DecontX Contamination",
        groupBy=groupBy,
        violin=violin,
        boxplot=boxplot,
        dots=dots,
        transparency=transparency,
        baseSize=baseSize,
        title=title,
        titleSize=titleSize,
        defaultTheme=defaultTheme,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        dotSize=dotSize,
        summary=summary,
        summaryTextSize=summaryTextSize,
        combinePlot="all"
    ))
    res.list = c(res.list, violinContamination)
    }

    if(is.null(legendSize) && !is.null(baseSize)){
      legendSizeScatterCluster = baseSize - 1
    }else{
      legendSizeScatterCluster = legendSize
    }
    title = "DecontX Clusters"
    if(labelSamples && length(samples) > 1){
      title = paste0(title, ", ", x)
    }

    scatterCluster <- list(scatter_decontXClusters = plotSCEDimReduceColData(
        inSCE=inSCESub,
        sample=sampleSub,
        colorBy=clusterCol,
        conditionClass="factor",
        shape=shape,
        reducedDimName=reducedDimName,
        xlab=xlab,
        ylab=ylab,
        dim1=dim1,
        dim2=dim2,
        bin=bin,
        binLabel=binLabel,
        baseSize=baseSize,
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        dotSize=dotSize,
        transparency=transparency,
        defaultTheme=defaultTheme,
        title=title,
        titleSize=titleSize,
        labelClusters=labelClusters,
        clusterLabelSize = clusterLabelSize,
        legendTitle="Clusters",
        legendSize=legendSizeScatterCluster,
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
      plotLabels <- "none"
      relHeights=1
  }
  if(!is.null(combinePlot)){
    if(combinePlot %in% c("all", "sample")){
      plotlist <- .ggSCTKCombinePlots(plotlist, combinePlot = combinePlot,
                                      relHeights = relHeights,
                                      relWidths = relWidths,
                                      ncols = plotNCols,
                                      nrows = plotNRows,
                                      samplePerColumn = samplePerColumn,
                                      sampleRelHeights = sampleRelHeights,
                                      sampleRelWidths = sampleRelWidths)
    }
  }
  return(plotlist)
}
