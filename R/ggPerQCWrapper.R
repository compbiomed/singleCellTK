#' @title Plots for runPerCellQC outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runPerCellQC function stored in the colData slot of the SingleCellExperiment
#'  object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runPerCellQC. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default FALSE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param dotSize Size of dots. Default 1.
#' @param axisSize Size of x/y-axis labels. Default 10.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @examples
#' \donttest{
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runPerCellQC(sce)
#' plotRunPerCellQCResults(inSCE = sce)
#' }
#' @export
plotRunPerCellQCResults <- function(inSCE,
                                    sample = NULL,
                                    groupby = NULL,
                                    violin = TRUE,
                                    boxplot = FALSE,
                                    dots = TRUE,
                                    dotSize = 0.5,
                                    axisSize = 10,
                                    transparency = 1,
                                    defaultTheme = TRUE) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)

    if(length(samples) > 1){
        combined.sum <- plotSCEViolinColData(inSCE=inSCE,
                                              coldata="sum",
                                              groupby=sample,
                                              xlab="",
                                              ylab="Counts",
                                              violin = violin,
                                              boxplot = boxplot,
                                              dots=dots,
                                              transparency = transparency,
                                              title="Total counts per cell",
                                              dotSize=dotSize,
                                              gridLine = TRUE,
                                              summary = "mean")
        combined.detected <- plotSCEViolinColData(inSCE=inSCE,
                                                  coldata="detected",
                                                  groupby=sample,
                                                  xlab="",
                                                  ylab="Features",
                                                  violin = violin,
                                                  boxplot = boxplot,
                                                  dots=dots,
                                                  transparency = transparency,
                                                  title="Total features detected per cell",
                                                  dotSize=dotSize,
                                                  gridLine = TRUE,
                                                  summary = "mean")
        combined.plots <- list(combined.sum, combined.detected)
        names(combined.plots) <- c("Sum", "Detected")
    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        violin.sum <-plotSCEViolinColData(inSCE=inSCESub,
                                          coldata="sum",
                                          sample = sampleSub,
                                          xlab="",
                                          ylab="Counts",
                                          groupby = groupby,
                                          violin = violin,
                                          boxplot = boxplot,
                                          dots=dots,
                                          transparency = transparency,
                                          title="Total counts per cell",
                                          dotSize=dotSize)

        violin.detected <- plotSCEViolinColData(inSCE=inSCESub,
                                                coldata="detected",
                                                sample = sampleSub,
                                                xlab="",
                                                ylab="Features",
                                                groupby = groupby,
                                                violin = violin,
                                                boxplot = boxplot,
                                                dots=dots,
                                                transparency = transparency,
                                                title="Total features detected per cell",
                                                dotSize=dotSize)

        violin.toppercent <- plotSCEViolinColData(inSCE=inSCESub,
                                                  coldata="percent_top_50",
                                                  sample = sampleSub,
                                                  xlab="",
                                                  ylab="Gene expression percentage (%)",
                                                  groupby = groupby,
                                                  violin = violin,
                                                  boxplot = boxplot,
                                                  dots=dots,
                                                  transparency = transparency,
                                                  title="Top 50 gene expression percentage",
                                                  dotSize=dotSize)

        if(any(grepl(pattern = "subsets_",
                     names(colData(inSCESub))))){

            subsets <- grep(pattern = "subsets_",
                            names(colData(inSCESub)), value = TRUE)

            violin.subset <- lapply(subsets, function(x) plotSCEViolinColData(inSCE=inSCESub,
                                                                              coldata=x,
                                                                              sample = sampleSub,
                                                                              xlab="",
                                                                              ylab=x,
                                                                              groupby = groupby,
                                                                              violin = violin,
                                                                              boxplot = boxplot,
                                                                              dots=dots,
                                                                              transparency = transparency,
                                                                              title=paste0(x," per cell"),
                                                                              dotSize=dotSize))
            names(violin.subset) <- subsets
        }else{
            violin.subset <- NULL
        }

        res.list <- list(violin.sum,
                         violin.detected,
                         violin.toppercent)
        names(res.list) <- c("sum", "detected", "toppercent")
        if(!is.null(violin.subset)){
            res.list <- c(res.list, violin.subset)
        }
        return(res.list)
    })

    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)
}

#' @title Plots for runScrublet outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runScrublet function stored in the colData slot of the SingleCellExperiment
#'  object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runScrublet. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runScrublet(sce)
#' plotScrubletResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotScrubletResults <- function(inSCE,
                                sample = NULL,
                                shape = NULL,
                                groupby = NULL,
                                violin = TRUE,
                                boxplot = FALSE,
                                dots = TRUE,
                                reducedDimName,
                                xlab = NULL,
                                ylab = NULL,
                                dim1 = NULL,
                                dim2 = NULL,
                                bin = NULL,
                                binLabel = NULL,
                                dotSize = 0.5,
                                transparency = 1,
                                defaultTheme = TRUE,
                                titleSize = 15) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)

    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                            coldata="scrublet_score",
                                            groupby=sample,
                                            xlab="",
                                            ylab="Doublet Score",
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            title="Scrublet Score",
                                            dotSize=dotSize,
                                            gridLine = TRUE,
                                            summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "Scrublet_Score"

    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        scatterScore <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "scrublet_score",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "Scrublet Doublet Score",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Doublet Score")

        densityScore <- plotSCEDensityColData(inSCE = inSCESub,
                                              sample = sampleSub,
                                              coldata = "scrublet_score",
                                              groupby = groupby,
                                              xlab = "Score",
                                              ylab = "Density",
                                              axisSize = 10,
                                              defaultTheme = defaultTheme,
                                              title = "Density, Scrublet Score")

        violinScore <- plotSCEViolinColData(inSCE=inSCESub, coldata="scrublet_score",
                                            sample = sampleSub,
                                            xlab="",
                                            ylab="Doublet Score",
                                            groupby = groupby,
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            defaultTheme = defaultTheme,
                                            title="Scrublet Score",
                                            dotSize=dotSize)

        scatterCall <- plotSCEDimReduceColData(inSCE = inSCESub,
                                               sample = sampleSub,
                                               colorBy = "scrublet_call",
                                               conditionClass = "factor",
                                               shape = shape,
                                               reducedDimName = reducedDimName,
                                               xlab = xlab,
                                               ylab = ylab,
                                               dim1 = dim1,
                                               dim2 = dim2,
                                               bin = bin,
                                               binLabel = binLabel,
                                               dotSize = dotSize,
                                               transparency = transparency,
                                               defaultTheme = defaultTheme,
                                               title = "Scrublet Doublet Call",
                                               titleSize = titleSize,
                                               labelClusters = FALSE,
                                               legendTitle = "Doublet Call")

        res.list <- list(scatterScore, densityScore, violinScore, scatterCall)
        names(res.list) <- c("scatterScore", "densityScore", "violinScore", "scatterCall")
        return(res.list)
    })
    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)

}

#' @title Plots for runDoubletFinder outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDoubletFinder function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runDoubletFinder. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runDoubletFinder(sce)
#' plotDoubletFinderResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotDoubletFinderResults <- function(inSCE,
                                     sample = NULL,
                                     shape = NULL,
                                     groupby = NULL,
                                     violin = TRUE,
                                     boxplot = FALSE,
                                     dots = TRUE,
                                     reducedDimName = NULL,
                                     xlab = NULL,
                                     ylab = NULL,
                                     dim1 = NULL,
                                     dim2 = NULL,
                                     bin = NULL,
                                     binLabel = NULL,
                                     dotSize = 0.5,
                                     transparency = 1,
                                     defaultTheme = TRUE,
                                     titleSize = 15) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }
    samples <- unique(sample)
    df.scores <- grep(pattern = "doubletFinder_doublet_score_Resolution_",
                      names(colData(inSCE)), value = TRUE)

    df.labels <- grep(pattern = "doubletFinder_doublet_label_Resolution_",
                      names(colData(inSCE)), value = TRUE)
    if(length(samples) > 1){
        combined.plots <- lapply(df.scores, function(x) plotSCEViolinColData(
            inSCE=inSCE,
            coldata=x,
            sample = NULL,
            xlab="",
            ylab="Doublet Score",
            groupby = sample,
            violin = violin,
            boxplot = boxplot,
            dots=TRUE,
            transparency = transparency,
            defaultTheme = defaultTheme,
            title=paste("DoubletFinder Score Resolution",
                        gsub(pattern = "doubletFinder_doublet_score_Resolution_",
                             "", x)),
            dotSize=dotSize))

        names(combined.plots) <- sapply(df.scores, function(x)
            paste0("Violin_", gsub(pattern = "doubletFinder_doublet_score_",
                                   "", x = x)))
        # combined.plots <- list(combined.plots)
        names(combined.plots) <- "DoubletFinder_Score"
    }


    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]

        scatterScore <- lapply(df.scores, function(x) plotSCEDimReduceColData(
            inSCE = inSCESub,
            sample = sampleSub,
            conditionClass = "numeric",
            shape = shape,
            colorBy = x,
            reducedDimName = reducedDimName,
            xlab = xlab,
            ylab = ylab,
            dim1 = dim1,
            dim2 = dim2,
            bin = bin,
            binLabel = binLabel,
            dotSize = dotSize,
            transparency = transparency,
            defaultTheme = defaultTheme,
            title = paste("DoubletFinder Doublet Score Resolution",
                          gsub(pattern = "doubletFinder_doublet_score_Resolution_",
                               "", x)),
            titleSize = titleSize,
            labelClusters = FALSE,
            legendTitle = "Doublet Score"
        ))

        names(scatterScore) <- sapply(df.scores, function(x)
            paste0("Scatter_Score_", gsub(pattern = "doubletFinder_doublet_score_",
                                          "", x = x)))

        densityScore <- lapply(df.scores, function(x) plotSCEDensityColData(
            inSCE = inSCESub,
            sample = sampleSub,
            coldata = x,
            groupby = groupby,
            xlab = "Score",
            ylab = "Density",
            axisSize = 10,
            defaultTheme = defaultTheme,
            title = paste("Density, DoubletFinder Score Resolution",
                          gsub(pattern = "doubletFinder_doublet_score_Resolution_",
                               "", x))))

        names(densityScore) <- sapply(df.scores, function(x)
            paste0("Density_", gsub(pattern = "doubletFinder_doublet_score_",
                                    "", x = x)))

        violinScore <- lapply(df.scores, function(x) plotSCEViolinColData(
            inSCE=inSCESub,
            coldata=x,
            sample = sampleSub,
            xlab="",
            ylab="Doublet Score",
            groupby = groupby,
            violin = violin,
            boxplot = boxplot,
            dots=dots,
            transparency = transparency,
            defaultTheme = defaultTheme,
            title=paste("DoubletFinder Score Resolution",
                        gsub(pattern = "doubletFinder_doublet_score_Resolution_",
                             "", x)),
            dotSize=dotSize))

        names(violinScore) <- sapply(df.scores, function(x)
            paste0("Violin_", gsub(pattern = "doubletFinder_doublet_score_",
                                    "", x = x)))

        scatterCall <- lapply(df.labels, function(x) plotSCEDimReduceColData(
            inSCE = inSCESub,
            sample = sampleSub,
            conditionClass = "factor",
            shape = shape,
            colorBy = x,
            reducedDimName = reducedDimName,
            xlab = xlab,
            ylab = ylab,
            dim1 = dim1,
            dim2 = dim2,
            bin = bin,
            binLabel = binLabel,
            dotSize = dotSize,
            transparency = transparency,
            defaultTheme = defaultTheme,
            title = paste("DoubletFinder Doublet Call Resolution",
                          gsub(pattern = "doubletFinder_doublet_label_Resolution_",
                               "", x)),
            titleSize = titleSize,
            labelClusters = FALSE,
            legendTitle = "Doublet Score"
        ))

        names(scatterCall) <- sapply(df.labels, function(x)
            paste0("Scatter_Call_",gsub(pattern = "doubletFinder_doublet_label_",
                                        "", x = x)))

        res.list <- c(scatterScore, densityScore, violinScore, scatterCall)
        return(res.list)
    })
    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)
}

#' @title Plots for runDoubletCells outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDoubletCells function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runDoubletCells. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runDoubletCells(sce)
#' plotDoubletCellsResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotDoubletCellsResults <- function(inSCE,
                                    sample = NULL,
                                    shape = NULL,
                                    groupby = NULL,
                                    violin = TRUE,
                                    boxplot = FALSE,
                                    dots = TRUE,
                                    reducedDimName = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    dim1 = NULL,
                                    dim2 = NULL,
                                    bin = NULL,
                                    binLabel = NULL,
                                    dotSize = 0.5,
                                    transparency = 1,
                                    defaultTheme = TRUE,
                                    titleSize = 15) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }
    samples <- unique(sample)
    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                                coldata="scran_doubletCells_Score",
                                                groupby=sample,
                                                xlab="",
                                                ylab="Doublet Score",
                                                violin = violin,
                                                boxplot = boxplot,
                                                dots=dots,
                                                transparency = transparency,
                                                title="DoubletCells Doublet Score",
                                                dotSize=dotSize,
                                                gridLine = TRUE,
                                                summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "DoubletCells_Score"
    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        scatterScore <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "scran_doubletCells_Score",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "DoubletCells Doublet Score",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Doublet Score")

        densityScore <- plotSCEDensityColData(inSCE = inSCESub,
                                              sample = sampleSub,
                                              coldata = "scran_doubletCells_Score",
                                              groupby = groupby,
                                              xlab = "Score",
                                              ylab = "Density",
                                              axisSize = 10,
                                              defaultTheme = defaultTheme,
                                              title = "Density, DoubletCells Score")

        violinScore <- plotSCEViolinColData(inSCE=inSCESub,
                                            coldata="scran_doubletCells_Score",
                                            sample = sampleSub,
                                            xlab="",
                                            ylab="Doublet Score",
                                            groupby = groupby,
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            title="DoubletCells Doublet Score",
                                            defaultTheme = defaultTheme,
                                            dotSize=dotSize)

        res.list <- list(scatterScore, densityScore, violinScore)
        names(res.list) <- c("scatterScore", "densityScore", "violinScore")
        return(res.list)
    })
    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)

}

#' @title Plots for runCxds outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runCxds function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runCxds. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runCxds(sce)
#' plotCxdsResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotCxdsResults <- function(inSCE,
                            sample = NULL,
                            shape = NULL,
                            groupby = NULL,
                            violin = TRUE,
                            boxplot = FALSE,
                            dots = TRUE,
                            reducedDimName = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            dim1 = NULL,
                            dim2 = NULL,
                            bin = NULL,
                            binLabel = NULL,
                            dotSize = 0.5,
                            transparency = 1,
                            defaultTheme = TRUE,
                            titleSize = 15) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)
    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                               coldata="scds_cxds_score",
                                               groupby=sample,
                                               xlab="",
                                               ylab="Doublet Score",
                                               violin = violin,
                                               boxplot = boxplot,
                                               dots=dots,
                                               transparency = transparency,
                                               title="CXDS Doublet Score",
                                               dotSize=dotSize,
                                               gridLine = TRUE,
                                               summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "CXDS_Score"
    }


    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        scatterScore <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "scds_cxds_score",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "CXDS Doublet Score",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Doublet Score")

        densityScore <- plotSCEDensityColData(inSCE = inSCESub,
                                              sample = sampleSub,
                                              coldata = "scds_cxds_score",
                                              groupby = groupby,
                                              xlab = "Score",
                                              ylab = "Density",
                                              axisSize = 10,
                                              defaultTheme = defaultTheme,
                                              title = "Density, CXDS Score")

        violinScore <- plotSCEViolinColData(inSCE=inSCESub,
                                            coldata="scds_cxds_score",
                                            sample = sampleSub,
                                            xlab="",
                                            ylab="Doublet Score",
                                            groupby = groupby,
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            title="CXDS Doublet Score",
                                            defaultTheme = defaultTheme,
                                            dotSize=dotSize)

        res.list <- list(scatterScore, densityScore, violinScore)
        names(res.list) <- c("scatterScore", "densityScore", "violinScore")
        return(res.list)
    })
    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)
}

#' @title Plots for runBcds outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runBcds function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runBcds. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runBcds(sce)
#' plotBcdsResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotBcdsResults <- function(inSCE,
                            sample = NULL,
                            shape = NULL,
                            groupby = NULL,
                            violin = TRUE,
                            boxplot = FALSE,
                            dots = TRUE,
                            reducedDimName = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            dim1 = NULL,
                            dim2 = NULL,
                            bin = NULL,
                            binLabel = NULL,
                            dotSize = 0.5,
                            transparency = 1,
                            defaultTheme = TRUE,
                            titleSize = 15) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)
    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                               coldata="scds_bcds_score",
                                               groupby=sample,
                                               xlab="",
                                               ylab="Doublet Score",
                                               violin = violin,
                                               boxplot = boxplot,
                                               dots=dots,
                                               transparency = transparency,
                                               title="BCDS Doublet Score",
                                               dotSize=dotSize,
                                               gridLine = TRUE,
                                               summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "BCDS_Score"
    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        scatterScore <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "scds_bcds_score",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "BCDS Doublet Score",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Doublet Score")

        densityScore <- plotSCEDensityColData(inSCE = inSCESub,
                                              sample = sampleSub,
                                              coldata = "scds_bcds_score",
                                              groupby = groupby,
                                              xlab = "Score",
                                              ylab = "Density",
                                              axisSize = 10,
                                              defaultTheme = defaultTheme,
                                              title = "Density, BCDS Score")


        violinScore <- plotSCEViolinColData(inSCE=inSCESub,
                                            coldata="scds_bcds_score",
                                            sample = sampleSub,
                                            xlab="",
                                            ylab="Doublet Score",
                                            groupby = groupby,
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            defaultTheme = defaultTheme,
                                            title="BCDS Doublet Score",
                                            dotSize=dotSize)

        res.list <- list(scatterScore, densityScore, violinScore)
        names(res.list) <- c("scatterScore", "densityScore", "violinScore")
        return(res.list)
    })
    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)

}

#' @title Plots for runCxdsBcdsHybrid outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runCxdsBcdsHybrid function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runCxdsBcdsHybrid.
#'  Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runCxdsBcdsHybrid(sce)
#' plotScdsHybridResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotScdsHybridResults <- function(inSCE,
                                  sample = NULL,
                                  shape = NULL,
                                  groupby = NULL,
                                  violin = TRUE,
                                  boxplot = FALSE,
                                  dots = TRUE,
                                  reducedDimName = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  dim1 = NULL,
                                  dim2 = NULL,
                                  bin = NULL,
                                  binLabel = NULL,
                                  dotSize = 0.5,
                                  transparency = 1,
                                  defaultTheme = TRUE,
                                  titleSize = 15) {
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)
    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                               coldata="scds_hybrid_score",
                                               groupby=sample,
                                               xlab="",
                                               ylab="Doublet Score",
                                               violin = violin,
                                               boxplot = boxplot,
                                               dots=dots,
                                               transparency = transparency,
                                               title="CXDS BCDS Doublet Score",
                                               dotSize=dotSize,
                                               gridLine = TRUE,
                                               summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "CXDS_BCDS_Score"
    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub = inSCE[,sampleInd]
        scatterScore <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "scds_hybrid_score",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "CXDS BCDS Doublet Score",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Doublet Score")

        densityScore <- plotSCEDensityColData(inSCE = inSCESub,
                                              sample = sampleSub,
                                              coldata = "scds_hybrid_score",
                                              groupby = groupby,
                                              xlab = "Score",
                                              ylab = "Density",
                                              axisSize = 10,
                                              defaultTheme = defaultTheme,
                                              title = "Density, CXDS BCDS Hybrid Score")

        violinScore <- plotSCEViolinColData(inSCE=inSCESub,
                                            coldata="scds_hybrid_score",
                                            sample = sampleSub,
                                            xlab="",
                                            ylab="Doublet Score",
                                            groupby = groupby,
                                            violin = violin,
                                            boxplot = boxplot,
                                            dots=dots,
                                            transparency = transparency,
                                            defaultTheme = defaultTheme,
                                            title="CXDS BCDS Doublet Score",
                                            dotSize=dotSize)

        res.list <- list(scatterScore, densityScore, violinScore)
        names(res.list) <- c("scatterScore", "densityScore", "violinScore")
        return(res.list)
    })

    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)
    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }
    return(plotlist)
}

#' @title Plots for runDecontX outputs.
#' @description A wrapper function which visualizes outputs from the
#'  runDecontX function stored in the colData slot of the
#'  SingleCellExperiment object via various plots.
#' @param inSCE Input SCtkExperiment object with saved dimension reduction
#'  components or a variable with saved results from runDecontX. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Default NULL.
#' @param shape If provided, add shapes based on the value.
#' @param groupby Groupings for each numeric value. A user may input a vector
#'  equal length to the number of the samples in the SingleCellExperiment
#'  object, or can be retrieved from the colData slot. Default NULL.
#' @param violin Boolean. If TRUE, will plot the violin plot. Default TRUE.
#' @param boxplot Boolean. If TRUE, will plot boxplots for each violin plot.
#'  Default TRUE.
#' @param dots Boolean. If TRUE, will plot dots for each violin plot.
#'  Default TRUE.
#' @param reducedDimName Saved dimension reduction name in the SCtkExperiment
#'  object. Required.
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
#' @param dotSize Size of dots. Default 1.
#' @param transparency Transparency of the dots, values will be 0-1. Default 1.
#' @param defaultTheme Removes grid in plot and sets axis title size to 10
#'  when TRUE. Default TRUE.
#' @param titleSize Size of title of plot. Default 15.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", reducedDimName = "UMAP")
#' sce <- runDecontX(sce)
#' plotDecontXResults(inSCE = sce, reducedDimName = "UMAP")
#' @export
plotDecontXResults <- function(inSCE,
                               sample = NULL,
                               shape = NULL,
                               groupby = NULL,
                               violin = TRUE,
                               boxplot = FALSE,
                               dots = TRUE,
                               reducedDimName = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               dim1 = NULL,
                               dim2 = NULL,
                               bin = NULL,
                               binLabel = NULL,
                               dotSize = 0.5,
                               transparency = 1,
                               defaultTheme = TRUE,
                               titleSize = 15) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    samples <- unique(sample)

    if(length(samples) > 1){
        combined.plots <- plotSCEViolinColData(inSCE=inSCE,
                                                coldata="decontX_contamination",
                                                groupby=sample,
                                                xlab="",
                                                ylab="DecontX Contamination",
                                                violin = violin,
                                                boxplot = boxplot,
                                                dots=dots,
                                                transparency = transparency,
                                                title="DecontX Contamination Score",
                                                dotSize=dotSize,
                                                gridLine = TRUE,
                                                summary = "mean")
        combined.plots <- list(combined.plots)
        names(combined.plots) <- "DecontX_Contamination"
    }

    plotlist <- lapply(samples, function(x) {
        sampleInd <- which(sample == x)
        sampleSub <- sample[sampleInd]
        inSCESub <- inSCE[,sampleInd]
        scatterDecon <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                sample = sampleSub,
                                                colorBy = "decontX_contamination",
                                                conditionClass = "numeric",
                                                shape = shape,
                                                reducedDimName = reducedDimName,
                                                xlab = xlab,
                                                ylab = ylab,
                                                dim1 = dim1,
                                                dim2 = dim2,
                                                bin = bin,
                                                binLabel = binLabel,
                                                dotSize = dotSize,
                                                transparency = transparency,
                                                defaultTheme = defaultTheme,
                                                title = "DecontX Contamination",
                                                titleSize = titleSize,
                                                labelClusters = FALSE,
                                                legendTitle = "Contamination")

        densityContamination <- plotSCEDensityColData(inSCE = inSCESub,
                                                      sample = sampleSub,
                                                      coldata = "decontX_contamination",
                                                      groupby = groupby,
                                                      xlab = "Score",
                                                      ylab = "Density",
                                                      axisSize = 10,
                                                      defaultTheme = defaultTheme,
                                                      title = "Density, DecontX Contamination Score")

        violinContamination <- plotSCEViolinColData(inSCE=inSCESub,
                                                    coldata="decontX_contamination",
                                                    sample = sampleSub,
                                                    xlab="", ylab="DecontX Contamination",
                                                    groupby = groupby, violin = violin, boxplot = boxplot, dots=dots, transparency = transparency,
                                                    title="DecontX Contamination Score",
                                                    defaultTheme = defaultTheme,
                                                    dotSize=dotSize)

        scatterCluster <- plotSCEDimReduceColData(inSCE = inSCESub,
                                                  sample = sampleSub,
                                                  colorBy = "decontX_clusters",
                                                  conditionClass = "factor",
                                                  shape = shape,
                                                  reducedDimName = reducedDimName,
                                                  xlab = xlab,
                                                  ylab = ylab,
                                                  dim1 = dim1,
                                                  dim2 = dim2,
                                                  bin = bin,
                                                  binLabel = binLabel,
                                                  dotSize = dotSize,
                                                  transparency = transparency,
                                                  defaultTheme = defaultTheme,
                                                  title = "DecontX Clusters",
                                                  titleSize = titleSize,
                                                  labelClusters = TRUE,
                                                  legendTitle = "Clusters")

        res.list <- list(scatterDecon, densityContamination,
                         violinContamination, scatterCluster)
        names(res.list) <- c("scatterDecon",
                             "densityContamination",
                             "violinContamination",
                             "scatterCluster")
        return(res.list)
    })

    if(length(unique(samples)) > 1){
        names(plotlist) <- samples
        plotlist <- c(combined.plots, plotlist)

    }else{
        plotlist <- unlist(plotlist, recursive = F)
    }

    return(plotlist)
}

