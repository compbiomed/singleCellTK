#' @title Plots for runPerCellQC outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runPerCellQC function stored in the colData slot of the SingleCellExperiment
#'  object via various plots.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results from
#' runPerCellQC. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' \donttest{
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runPerCellQC(sce)
#' plotRunPerCellQCResults(inSCE=sce)
#' }
#' @export
plotRunPerCellQCResults <- function(inSCE,
                                    sample=NULL,
                                    groupby=NULL,
                                    combinePlot=FALSE,
                                    violin=TRUE,
                                    boxplot=FALSE,
                                    dots=TRUE,
                                    dotSize=1,
                                    axisSize=15,
                                    axisLabelSize=18,
                                    transparency=1,
                                    defaultTheme=TRUE,
                                    titleSize=18) {
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

  if (length(samples) > 1) {
    combined.sum <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="sum",
      groupby=sample,
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
      titleSize=titleSize
    )
    combined.detected <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="detected",
      groupby=sample,
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
      titleSize=titleSize
    )
    merged.plots <- list(combined.sum, combined.detected)
    names(merged.plots) <- c("Sum", "Detected")
  }

  res.list <- c()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    if (length(samples) == 1) {
        violin.sum <- list(sum = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="sum",
          sample=sampleSub,
          xlab="",
          ylab="Counts",
          groupby=groupby,
          violin=violin,
          boxplot=boxplot,
          dots=dots,
          transparency=transparency,
          title="Total counts per cell",
          dotSize=dotSize,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          summary="median",
          titleSize=titleSize
        ))
        res.list <- c(res.list, violin.sum)
    }

    if (length(samples) == 1) {
        violin.detected <- list(detected = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="detected",
          sample=sampleSub,
          xlab="",
          ylab="Features",
          groupby=groupby,
          violin=violin,
          boxplot=boxplot,
          dots=dots,
          transparency=transparency,
          title="Total features detected per cell",
          dotSize=dotSize,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          summary="median",
          titleSize=titleSize
        ))
        res.list <- c(res.list, violin.detected)
    }

    violin.toppercent <- list(toppercent = plotSCEViolinColData(
      inSCE=inSCESub,
      coldata="percent_top_50",
      sample=sampleSub,
      xlab="",
      ylab="Gene expression percentage (%)",
      groupby=groupby,
      violin=violin,
      boxplot=boxplot,
      dots=dots,
      transparency=transparency,
      title="Top 50 gene expression percentage",
      dotSize=dotSize,
      axisSize=axisSize,
      axisLabelSize=axisLabelSize,
      summary="median",
      titleSize=titleSize
    ))
    res.list <- c(res.list, violin.toppercent)

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
          groupby=groupby,
          violin=violin,
          boxplot=boxplot,
          summary="median",
          dots=dots,
          transparency=transparency,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          title=paste0(x, " per cell"),
          dotSize=dotSize,
          titleSize=titleSize
        )
      })
      names(violin.subset) <- subsets
    } else {
      violin.subset <- NULL
    }

    if (!is.null(violin.subset)) {
      res.list <- c(res.list, violin.subset)
    }
    return(res.list)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,
                                      relHeights = c(1.5, 1.5, 1, 1),
                                      labels = samples)
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- runEmptyDrops(inSCE=sce)
#' plotEmptyDropsResults(inSCE=sce)
#' @export
plotEmptyDropsResults <- function(inSCE,
                                  sample=NULL,
                                  fdrCutoff=0.01,
                                  defaultTheme=TRUE,
                                  dotSize=1,
                                  titleSize=18,
                                  axisLabelSize=18,
                                  axisSize=15,
                                  legendSize=15,
                                  legendTitleSize=16) {
  scatterEmptyDrops <- plotEmptyDropsScatter(inSCE,
    sample=sample,
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
    legendSize=legendSize
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runScrublet(sce)
#' plotScrubletResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotScrubletResults <- function(inSCE,
                                sample=NULL,
                                shape=NULL,
                                groupby=NULL,
                                combinePlot=FALSE,
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
                                legendTitleSize=16) {
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

  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scrublet_score",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "Scrublet_Score"
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
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=15,
        axisLabelSize=18,
        defaultTheme=defaultTheme,
        cutoff=0.5,
        title="Density, Scrublet Score",
        titleSize=titleSize
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
      legendTitleSize=legendTitleSize
    ))
    res.list <- c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- list(violin_doubletScore = plotSCEViolinColData(
          inSCE=inSCESub, coldata="scrublet_score",
          sample=sampleSub,
          xlab="",
          ylab="Doublet Score",
          groupby=groupby,
          violin=violin,
          boxplot=boxplot,
          dots=dots,
          transparency=transparency,
          defaultTheme=defaultTheme,
          title="Scrublet Score",
          titleSize=titleSize,
          dotSize=dotSize,
          axisSize=axisSize, axisLabelSize=axisLabelSize,
          summary="median"
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
      legendTitleSize=16,
      legendSize=15
    ))
    res.list <- c(res.list, scatterCall)

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }

  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist, relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDoubletFinder(sce)
#' plotDoubletFinderResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDoubletFinderResults <- function(inSCE,
                                     sample=NULL,
                                     shape=NULL,
                                     groupby=NULL,
                                     combinePlot=FALSE,
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
                                     legendTitleSize=16) {
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
        groupby=sample,
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
            pattern="doubletFinder_doublet_score_resolution_",
            "", x
          )
        ),
        dotSize=dotSize,
        summary="median"
      )
    })

    names(merged.plots) <- sapply(df.scores, function(x) {
      paste0("Violin_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=x
      ))
    })
    # merged.plots <- list(merged.plots)
    names(merged.plots) <- "DoubletFinder_Score"
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
            groupby=groupby,
            xlab="Score",
            ylab="Density",
            axisSize=axisSize, axisLabelSize=axisLabelSize,
            defaultTheme=defaultTheme,
            cutoff=0.5,
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
    names(densityScore) <- sapply(df.scores, function(x) {
        paste0("Density_", gsub(
            pattern="doubletFinder_doublet_score_",
            "", x=x
        ))
    })
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
        legendTitleSize=legendTitleSize
      )
    })

    names(scatterScore) <- sapply(df.scores, function(x) {
      paste0("scatter_score_", gsub(
        pattern="doubletFinder_doublet_score_",
        "", x=x
      ))
    })
    res.list <- c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- lapply(df.scores, function(x) {
          plotSCEViolinColData(
            inSCE=inSCESub,
            coldata=x,
            sample=sampleSub,
            xlab="",
            ylab="Doublet Score",
            groupby=groupby,
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
            axisLabelSize=axisLabelSize
          )
        })

        names(violinScore) <- sapply(df.scores, function(x) {
          paste0("violin_", gsub(
            pattern="doubletFinder_doublet_score_",
            "", x=x
          ))
        })
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
        legendTitle="Doublet Label",
        legendSize=legendSize,
        legendTitleSize=legendTitleSize
      )
    })

    names(scatterCall) <- sapply(df.labels, function(x) {
      paste0("Scatter_Call_", gsub(
        pattern="doubletFinder_doublet_label_",
        "", x=x
      ))
    })
    res.list <- c(res.list, scatterCall)
    return(res.list)
  })
  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param logScore Boolean. If TRUE, the doublet score will be log normalized.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDoubletCells(sce)
#' plotDoubletCellsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDoubletCellsResults <- function(inSCE,
                                    sample=NULL,
                                    shape=NULL,
                                    groupby=NULL,
                                    combinePlot=FALSE,
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
                                    legendTitleSize=16) {
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

  if (logScore) {
    colData(inSCE)$scran_doubletCells_score <- log10(colData(inSCE)$scran_doubletCells_score + 1)
    titleDoubletCells <- "DoubletCells Doublet Score, log10"
  } else {
    titleDoubletCells <- "DoubletCells Doublet Score"
  }

  samples <- unique(sample)
  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scran_doubletCells_score",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "DoubletCells_Score"
  }

  res.list <- list()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]
    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scran_doubletCells_score",
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title=paste0("Density, ", titleDoubletCells),
        titleSize=titleSize
    ))
    res.list = c(res.list, densityScore)

    scatterScore <- list(scatter_doubletScore = plotSCEDimReduceColData(
      inSCE=inSCESub,
      sample=sampleSub,
      colorBy="scran_doubletCells_score",
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
      legendTitleSize=legendTitleSize
    ))
    res.list = c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- list(violin_doubletScore = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="scran_doubletCells_score",
          sample=sampleSub,
          xlab="",
          ylab="Doublet Score",
          groupby=groupby,
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
          summary="median"
        ))
        res.list = c(res.list, violinScore)
    }

    return(res.list)
  })
  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runCxds(sce)
#' plotCxdsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotCxdsResults <- function(inSCE,
                            sample=NULL,
                            shape=NULL,
                            groupby=NULL,
                            combinePlot=FALSE,
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
                            legendTitleSize=16) {
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
  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_cxds_score",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "CXDS_Score"
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
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, CXDS Score",
        titleSize=titleSize
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
      legendTitleSize=legendTitleSize
    ))
    res.list = c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- list(violin_doubletScore = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="scds_cxds_score",
          sample=sampleSub,
          xlab="",
          ylab="Doublet Score",
          groupby=groupby,
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
          summary="median"
        ))
        res.list = c(res.list, violinScore)
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
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title="CXDS Doublet Assignment",
        titleSize=titleSize,
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendTitleSize=16,
        legendSize=15
    ))
    res.list <- c(res.list, scatterCall)
    return(res.list)
  })
  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runBcds(sce)
#' plotBcdsResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotBcdsResults <- function(inSCE,
                            sample=NULL,
                            shape=NULL,
                            groupby=NULL,
                            combinePlot=FALSE,
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
                            legendTitleSize=16) {
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
  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_bcds_score",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "BCDS_Score"
  }

  res.list <- c()
  plotlist <- lapply(samples, function(x) {
    sampleInd <- which(sample == x)
    sampleSub <- sample[sampleInd]
    inSCESub <- inSCE[, sampleInd]

    densityScore <- list(density_doubletScore = plotSCEDensityColData(
        inSCE=inSCESub,
        sample=sampleSub,
        coldata="scds_bcds_score",
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, BCDS Score"
    ))
    res.list <- c(res.list, densityScore)

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
      legendTitleSize=legendTitleSize
    ))
    res.list <- c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- list(violin_doubletScore = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="scds_bcds_score",
          sample=sampleSub,
          xlab="",
          ylab="Doublet Score",
          groupby=groupby,
          violin=violin,
          boxplot=boxplot,
          dots=dots,
          transparency=transparency,
          defaultTheme=defaultTheme,
          title="BCDS Doublet Score",
          titleSize=titleSize,
          dotSize=dotSize,
          axisSize=axisSize,
          axisLabelSize=axisLabelSize,
          summary="median"
        ))
        res.list <- c(res.list, violinScore)
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
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title="BCDS Doublet Assignment",
        titleSize=titleSize,
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendTitleSize=16,
        legendSize=15
    ))
    res.list <- c(res.list, scatterCall)

    return(res.list)
  })
  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runCxdsBcdsHybrid(sce)
#' plotScdsHybridResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotScdsHybridResults <- function(inSCE,
                                  sample=NULL,
                                  shape=NULL,
                                  groupby=NULL,
                                  combinePlot=FALSE,
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
                                  legendTitleSize=16) {
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
  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="scds_hybrid_score",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "CXDS_BCDS_Score"
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
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, CXDS BCDS Hybrid Score",
        titleSize=titleSize
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
      legendSize=legendSize
    ))
    res.list = c(res.list, scatterScore)

    if (length(samples) == 1) {
        violinScore <- list(violin_doubletScore = plotSCEViolinColData(
          inSCE=inSCESub,
          coldata="scds_hybrid_score",
          sample=sampleSub,
          xlab="",
          ylab="Doublet Score",
          groupby=groupby,
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
          summary="median"
        ))
        res.list = c(res.list, violinScore)
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
        colorScale = c("lightgray","red"),
        defaultTheme=defaultTheme,
        title="CXDS BCDS Doublet Assignment",
        titleSize=titleSize,
        axisSize=axisSize, axisLabelSize=axisLabelSize,
        labelClusters=FALSE,
        legendTitle="Doublet \nAssignment",
        legendTitleSize=16,
        legendSize=15
    ))
    res.list <- c(res.list, scatterCall)
    return(res.list)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,relHeights = c(1.5, 1, 1),
                                      labels = samples)
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
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param combinePlot Boolean. Will combine plots using `cowplot::plot_grid`.
#'  Default FALSE.
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
#' @examples
#' data(scExample, package="singleCellTK")
#' sce <- sce[, colData(sce)$type != "EmptyDroplet"]
#' sce <- getUMAP(inSCE=sce, useAssay="counts", reducedDimName="UMAP")
#' sce <- runDecontX(sce)
#' plotDecontXResults(inSCE=sce, reducedDimName="UMAP")
#' @export
plotDecontXResults <- function(inSCE,
                               sample=NULL,
                               shape=NULL,
                               groupby=NULL,
                               combinePlot=FALSE,
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
                               legendTitleSize=16) {
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

  if (length(samples) > 1) {
    merged.plots <- plotSCEViolinColData(
      inSCE=inSCE,
      coldata="decontX_contamination",
      groupby=sample,
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
      summary="median"
    )
    merged.plots <- list(merged.plots)
    names(merged.plots) <- "DecontX_Contamination"
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
        groupby=groupby,
        xlab="Score",
        ylab="Density",
        axisSize=axisSize,
        axisLabelSize=axisLabelSize,
        defaultTheme=defaultTheme,
        title="Density, DecontX Contamination Score",
        titleSize=titleSize
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
      legendSize=legendSize
    ))
    res.list = c(res.list, scatterContamination)

    if (length(samples) == 1) {
        violinContamination <- list(violin_decontXContamination = plotSCEViolinColData(
            inSCE=inSCESub,
            coldata="decontX_contamination",
            sample=sampleSub,
            xlab="", ylab="DecontX Contamination",
            groupby=groupby,
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
            summary="median"
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
        legendTitleSize=legendTitleSize
    ))
    res.list = c(res.list, scatterCluster)

    return(res.list)
  })

  if (length(unique(samples)) > 1) {
    names(plotlist) <- samples
    plotlist <- c(merged.plots, plotlist)
  } else {
    plotlist <- unlist(plotlist, recursive=F)
  }
  if(combinePlot){
      plotlist <- .ggSCTKCombinePlots(plotlist,
                                     relHeights = c(1.5, 1, 1),
                                      labels = samples)
  }
  return(plotlist)
}
