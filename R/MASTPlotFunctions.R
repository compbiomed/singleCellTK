#' Check if the specified MAST result in SCtkExperiment object is complete.
#' But does not garantee the biological correctness.
#' @param inSCE SCtkExperiment object. `runMAST()` has to be run in advance.
#' @param useResult character. A string specifying the `comparisonName` used.
checkMASTResult <- function(inSCE, useResult){
    if(!inherits(inSCE, 'SingleCellExperiment')){
        stop('Given object is not a valid SingleCellExperiment object.')
    }
    if(!'MAST' %in% names(S4Vectors::metadata(inSCE))){
        stop('"MAST" not in metadata, please run runMAST() first.')
    }
    if(!useResult %in% names(S4Vectors::metadata(inSCE)$MAST)){
        stop(paste0('"', useResult, '"', ' not in metadata(inSCE)$MAST. '),
             'Please check.')
    }
    result <- S4Vectors::metadata(inSCE)$MAST[[useResult]]
    if(!all(c('groupNames', 'select', 'result', 'useAssay') %in% names(result))){
        stop(paste0('"', useResult, '"', ' result is not complete. '),
             'You might need to rerun it.')
    }

    return(TRUE)
}

#' plot the violin plot to show visualize the expression distribution of DEGs
#' identified by MAST
#' @param inSCE SCtkExperiment object. `runMAST()` has to be run in advance.
#' @param useResult character. A string specifying the `comparisonName` used.
#' @param threshP logical, default FALSE. Whether to plot threshold values from
#' adaptive thresholding, instead of using the assay used by `runMAST()`
#' @param nrow integer, default 6. Number of rows in the plot grid.
#' @param ncol integer, default 6. Number of columns in the plot grid.
#' @return A ggplot object of MAST violin plot
#' @export
plotMASTViolin <- function(inSCE, useResult, threshP = FALSE,
                           nrow = 6, ncol = 6){
    #TODO: DO we split the up/down regulation too?
    # Check
    checkMASTResult(inSCE, useResult)
    # Extract
    result <- S4Vectors::metadata(inSCE)$MAST[[useResult]]
    deg <- result$result
    useAssay <- result$useAssay
    geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), Gene]
    groupName1 <- result$groupNames[1]
    ix1 <- result$select$ix1
    cells1 <- colnames(inSCE)[ix1]
    groupName2 <- result$groupNames[2]
    ix2 <- result$select$ix2
    cells2 <- colnames(inSCE)[ix2]
    expres <- SummarizedExperiment::assay(inSCE[geneToPlot, c(cells1, cells2)],
                                          useAssay)
    # Format
    cdat <- data.frame(wellKey = colnames(expres),
                       condition = factor(c(rep(groupName1, length(cells1)),
                                            rep(groupName2, length(cells2))),
                                          levels = result$groupNames),
                       ngeneson = rep("", (length(cells1) + length(cells2))),
                       stringsAsFactors = FALSE)
    sca <- suppressMessages(MAST::FromMatrix(expres, cdat))
    if(threshP){
        #TODO: if nrow*ncol < `min_per_bin`` below, there would be an error.
        invisible(utils::capture.output(thres <-
            MAST::thresholdSCRNACountMatrix(expres, nbins = 20,
                                            min_per_bin = 30)))
        SummarizedExperiment::assay(sca) <- thres$counts_threshold
    }
    flatDat <- methods::as(sca, "data.table")
    flatDat$primerid <- factor(flatDat$primerid, levels = geneToPlot)
    names(flatDat)[5] <- useAssay
    # Plot
    violinplot <- ggplot2::ggplot(flatDat,
        ggplot2::aes_string(x = 'condition', y = useAssay,
                            color = 'condition')) +
        ggplot2::geom_jitter() +
        ggplot2::facet_wrap(~primerid, scale = "free_y",
                            ncol = ncol) +
        ggplot2::geom_violin() +
        ggplot2::ggtitle(paste0("Violin Plot for ", useResult))
    return(violinplot)
}

#' plot the linear regression to show visualize the expression the of DEGs
#' identified by MAST
#' @param inSCE SCtkExperiment object. `runMAST()` has to be run in advance.
#' @param useResult character. A string specifying the `comparisonName` used.
#' @param threshP logical, default FALSE. Whether to plot threshold values from
#' adaptive thresholding, instead of using the assay used by `runMAST()`
#' @param nrow integer, default 6. Number of rows in the plot grid.
#' @param ncol integer, default 6. Number of columns in the plot grid.
#' @return A ggplot object of MAST linear regression
#' @export
plotMASTRegression <- function(inSCE, useResult, threshP = FALSE,
                               nrow = 6, ncol = 6){
    #TODO: DO we split the up/down regulation too?
    # Check
    checkMASTResult(inSCE, useResult)
    # Extract
    result <- S4Vectors::metadata(inSCE)$MAST[[useResult]]
    deg <- result$result
    useAssay <- result$useAssay
    geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), Gene]
    groupName1 <- result$groupNames[1]
    ix1 <- result$select$ix1
    cells1 <- colnames(inSCE)[ix1]
    groupName2 <- result$groupNames[2]
    ix2 <- result$select$ix2
    cells2 <- colnames(inSCE)[ix2]
    expres <- SummarizedExperiment::assay(inSCE[geneToPlot, c(cells1, cells2)],
                                          useAssay)
    # Format
    cdat <- data.frame(wellKey = colnames(expres),
                       condition = factor(c(rep(groupName1, length(cells1)),
                                            rep(groupName2, length(cells2))),
                                          levels = result$groupNames),
                       ngeneson = rep("", (length(cells1) + length(cells2))),
                       stringsAsFactors = FALSE)
    sca <- suppressMessages(MAST::FromMatrix(expres, cData = cdat))
    cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
    SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
    if(threshP){
        #TODO: if nrow*ncol < `min_per_bin`` below, there would be an error.
        invisible(utils::capture.output(thres <-
              MAST::thresholdSCRNACountMatrix(expres, nbins = 20,
                                              min_per_bin = 30)))
        SummarizedExperiment::assay(sca) <- thres$counts_threshold
    }
    flatDat <- methods::as(sca, "data.table")
    flatDat$primerid <- factor(flatDat$primerid, levels = geneToPlot)
    names(flatDat)[6] <- useAssay
    # Calculate
    resData <- NULL
    for (i in unique(flatDat$primerid)){
        resdf <- flatDat[flatDat$primerid == i, ]
        resdf$lmPred <- stats::lm(
            stats::as.formula(paste0(useAssay, "~cngeneson+", 'condition')),
            data = flatDat[flatDat$primerid == i, ])$fitted
        if (is.null(resData)){
            resData <- resdf
        } else {
            resData <- rbind(resData, resdf)
        }
    }
    # Plot
    ggbase <- ggplot2::ggplot(resData, ggplot2::aes_string(
                                  x = 'condition',
                                  y = useAssay,
                                  color = 'condition')) +
                              ggplot2::geom_jitter() +
                              ggplot2::facet_wrap(~primerid, scale = "free_y",
                                                  ncol = ncol)
    regressionplot <- ggbase +
                      ggplot2::aes_string(x = "cngeneson") +
                      ggplot2::geom_line(ggplot2::aes_string(y = "lmPred"),
                                         lty = 1) +
                      ggplot2::xlab("Standardized Cellular Detection Rate") +
                      ggplot2::ggtitle("Linear Model Plot")
    return(regressionplot)
}

#' plotMASTHeatmap helper function
#' Extract subsetted columns in factor class from col/rowData of SCE object
#' @param inSCE
#' @param axis choose from 'col' or 'row'
#' @param columns character. column name(s) of col/rowData
#' @param index The character/numeric/logical indices for subsetting.
#' @return A formatted data.frame object
extractData <- function(inSCE, axis = NULL, columns = NULL, index = NULL){
    if(is.null(columns)){
        return(NULL)
    } else {
        if(is.null(axis) || !axis %in% c('col', 'row')){
            stop("axis should be 'col' or 'row'.")
        } else if(axis == 'col'){
            data <- SummarizedExperiment::colData(inSCE)
        } else if(axis == 'row'){
            data <- SummarizedExperiment::rowData(inSCE)
        }
        if(!is.null(index)){
            data <- data[index,]
        }
        df <- data.frame(data[columns])
        for(i in colnames(df)){
            if(is.character(df[[i]]) || is.logical(df[[i]])){
                df[[i]] <- as.factor(df[[i]])
            }
        }
        return(df)
    }
}

#' plotMASTHeatmap helper function
#' Automatically setting the color list for the col/row annotation
#' @param inSCE
#' @param axis choose from 'col' or 'row'
#' @param columns character. column name(s) of col/rowData
#' @param index The character/numeric/logical indices for subsetting.
#' @return A list object formatted as required by celda::semiHeatmap
dataCol <- function(inSCE, axis = NULL, columns = NULL, index = NULL){
    if(is.null(columns)){
      return(NULL)
    } else {
        if(is.null(axis) || !axis %in% c('col', 'row')){
            stop("axis should be 'col' or 'row'.")
        } else if(axis == 'col'){
            data <- SummarizedExperiment::colData(inSCE)
        } else if(axis == 'row'){
            data <- SummarizedExperiment::rowData(inSCE)
        }
        #if(!is.null(index)){
        #    data <- data[index,]
        #}
        df <- data.frame(data[columns])
        colList <- list()
        for(i in colnames(df)){
            if(is.factor(df[[i]])){
                allLevels <- levels(df[[i]])
            } else {
                allLevels <- unique(df[[i]])
            }
            allLevels <- stringr::str_sort(allLevels, numeric = TRUE)
            allCol <- celda::distinctColors(length(allLevels))
            names(allCol) <- allLevels
            if(!is.null(index)){
                subseted <- df[index,i]
                uniqSubseted <- as.vector(unique(subseted))
                allCol <- allCol[names(allCol) %in% uniqSubseted]
            }
            colList[[i]] <- allCol
        }
        return(colList)
    }
}

#' Heatmap visualization of DEG result called by MAST
#'
#' runMAST() has to be run in advance so that information is stored in the
#' metadata of the input SCE object. This function is an modified version of
#' celda::plotHeatmap(). Additional arguments are annotated in document. For
#' others please refer to celda's document.
#' @param inSCE SingleCellExperiment object. Should be output of runMAST(),
#' with DEG information saved in metadata.
#' @param useResult character or numeric, default NULL. Specify which
#' comparison result to plot. Character refers to metadata(inSCE)$MAST; numeric
#' stands for index.
#' @param onlyPos logical, default FALSE. Whether to only plot DEG with
#' positive log2_FC value.
#' @param log2fcThreshold numeric, default 1. Only plot DEGs with the absolute
#' values of log2FC larger than this value.
#' @param fdrThreshold numeric, default 0.05. Only plot DEGs with FDR value
#' smaller than this value.
#' @param labelColData character vector, default NULL. Column names in
#' colData(inSCE), that needs to be annotated.
#' @param main character, default useResult.
#' @export
#' @note Make it available to adjust the order of legends
plotMASTHeatmap <- function(inSCE, useResult = NULL, onlyPos = FALSE,
    log2fcThreshold = 1, fdrThreshold = 0.05, annotationCell = NULL,
    annotationFeature = NULL, labelColData = NULL, labelRowData = NULL,
    trim = c(-2, 2), colorScheme = c("divergent", "sequential"),
    colorSchemeSymmetric = TRUE, colorSchemeCenter = 0, col = NULL,
    annotationColor = NULL, breaks = NULL, hclustMethod = "ward.D2",
    treeheightFeature = 30, treeheightCell = 30, silent = FALSE,
    main = useResult, ...){
    # Check
    checkMASTResult(inSCE, useResult)
    # Extract
    result <- S4Vectors::metadata(inSCE)$MAST[[useResult]]
    deg <- result$result
    useAssay <- result$useAssay
    expres <- as.matrix(SummarizedExperiment::assay(inSCE, useAssay))
    expres <- featureNameDedup(expres)
    ix1 <- result$select$ix1
    ix2 <- result$select$ix2
    filter <- which(deg$FDR < fdrThreshold & abs(deg$Log2_FC) > log2fcThreshold)
    deg.filtered <- deg[filter,]
    if(onlyPos){
        deg.filtered <- deg.filtered[which(deg.filtered$Log2_FC > 0),]
    }
    if(dim(deg.filtered)[1] <= 1){
        stop('Too few genes that pass filtration, unable to plot')
    }
    gene.ix <- rownames(inSCE) %in% deg.filtered$Gene
    cell.ix <- which(ix1 | ix2)
    expres <- expres[gene.ix, cell.ix]
    # Annotation organization
    ## Cells
    group <- vector()
    group[ix1] <- result$groupNames[1]
    group[ix2] <- result$groupNames[2]
    group <- factor(group[cell.ix], levels = result$groupNames)
    if(!is.null(annotationCell)){
        if(any(rownames(annotationCell) != colnames(expres))){
            stop("Rownames of 'annotationCell' are not matching cells ",
                 "involved in comparison")
        }
        annotationCell <- data.frame(condition = group,
                                     annotationCell)
    } else {
        annotationCell <- data.frame(condition = group,
                                     row.names = colnames(expres))
    }
    if(!is.null(labelColData)){
        if(!all(labelColData %in% names(SummarizedExperiment::colData(inSCE)))){
            stop('Given "labelColData" not matching colData.')
        }
        extractedCol <- extractData(inSCE, 'col', labelColData, cell.ix)
        annotationCell <- cbind(annotationCell, extractedCol)
        colDataColor <- dataCol(inSCE, 'col', labelColData, cell.ix)
        newEntries <- setdiff(names(colDataColor), names(annotationColor))
        colDataColor <- colDataColor[newEntries]
        annotationColor <- c(colDataColor, annotationColor)
    }
    kCol <- celda::distinctColors(2)
    names(kCol) <- result$groupNames
    if(!"comparisonGroup" %in% names(annotationColor)){
        annotationColor <- c(list(condition = kCol), annotationColor)
    }
    ## Genes
    regulation <- vector()
    genes.up <- deg.filtered[deg.filtered$Log2_FC > 0, Gene]
    genes.down <- deg.filtered[deg.filtered$Log2_FC < 0, Gene]
    regulation[rownames(inSCE) %in% genes.up] <- 'up'
    regulation[rownames(inSCE) %in% genes.down] <- 'down'
    regulation <- factor(regulation[gene.ix], levels = c('up', 'down'))
    if(!is.null(annotationFeature)){
        if(any(rownames(annotationFeature) != rownames(expres))){
            stop("Rownames of 'annotationFeature' are not matching genes involved in comparison")
        }
        annotationFeature <- data.frame(regulation = regulation,
                                        annotationFeature)
    } else {
        annotationFeature <- data.frame(regulation = regulation,
                                        row.names = rownames(expres))
    }
    if(!is.null(labelRowData)){
        if(!all(labelRowData %in% names(SummarizedExperiment::rowData(inSCE)))){
            stop('Given "labelRowData" not matching rowData.')
        }
        extractedRow <- extractData(inSCE, 'row', labelRowData, gene.ix)
        annotationFeature <- cbind(annotationFeature, extractedRow)
        rowDataColor <- dataCol(inSCE, 'row', labelRowData, gene.ix)
        newEntries <- setdiff(names(rowDataColor), names(annotationColor))
        rowDataColor <- rowDataColor[newEntries]
        annotationColor <- c(rowDataColor, annotationColor)
    }
    lCol <- celda::distinctColors(2)
    names(lCol) <- c('up', 'down')
    if(!"regulation" %in% names(annotationColor)){
        annotationColor <- c(list(regulation = lCol), annotationColor)
    }
    # celda::plotHeatmap remaining
    colorScheme <- match.arg(colorScheme)

    cn <- colnames(expres)
    expres <- t(base::apply(expres, 1, scale))
    colnames(expres) <- cn

    if (!is.null(trim)) {
        if (length(trim) != 2) {
            stop("'trim' should be a 2 element vector specifying the lower",
                 " and upper boundaries")
        }
        trim <- sort(trim)
        expres[expres < trim[1]] <- trim[1]
        expres[expres > trim[2]] <- trim[2]
    }
    uBoundRange <- max(expres, na.rm = TRUE)
    lboundRange <- min(expres, na.rm = TRUE)

    if (colorScheme == "divergent") {
        if (colorSchemeSymmetric == TRUE) {
            uBoundRange <- max(abs(uBoundRange), abs(lboundRange))
            lboundRange <- -uBoundRange
        }
        if (is.null(col)) {
            col <- colorRampPalette(c("#1E90FF", "#FFFFFF", "#CD2626"),
                                    space = "Lab")(100)
        }
        colLen <- length(col)
        if (is.null(breaks)) {
           breaks <- c(seq(lboundRange, colorSchemeCenter,
                           length.out = round(colLen/2) + 1),
                       seq(colorSchemeCenter + 1e-06, uBoundRange,
                           length.out = colLen - round(colLen/2)))
        }
    }
    else {
        if (is.null(col)) {
            col <- colorRampPalette(c("#FFFFFF",
                                      brewer.pal(n = 9, name = "Blues")))(100)
        }
        colLen <- length(col)
        if (is.null(breaks)) {
            breaks <- seq(lboundRange, uBoundRange, length.out = colLen)
        }
    }
    sp <- celda:::semiPheatmap(mat = expres, color = col, breaks = breaks,
        clusterCols = TRUE, clusterRows = TRUE,
        annotationRow = annotationFeature, annotationCol = annotationCell,
        annotationColors = annotationColor, legend = TRUE,
        annotationLegend = TRUE, annotationNamesRow = TRUE,
        annotationNamesCol = TRUE, showRownames = FALSE,
        showColnames = FALSE, clusteringMethod = hclustMethod,
        treeHeightRow = treeheightFeature, treeHeightCol = treeheightCell,
        rowLabel = regulation, colLabel = group, silent = TRUE, main = main, ...)
    if (!isTRUE(silent)) {
        grid::grid.newpage()
        grid::grid.draw(sp$gtable)
    }
    invisible(sp)
}

#' @describeIn MAST Identify adaptive thresholds
#'
#' @return thresholdGenes(): list of thresholded counts (on natural scale),
#' thresholds, bins, densities estimated on each bin, and the original data from
#' MAST::thresholdSCRNACountMatrix
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- thresholdGenes(mouseBrainSubsetSCE)
thresholdGenes <- function(inSCE, useAssay="logcounts"){
    # data preparation
    expres <- SummarizedExperiment::assay(inSCE, useAssay)
    fdata <- data.frame(Gene = rownames(expres))
    rownames(fdata) <- fdata$Gene
    SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                               fdata)
    SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
    invisible(utils::capture.output(thres <- MAST::thresholdSCRNACountMatrix(
        SummarizedExperiment::assay(SCENew), nbins = 20, min_per_bin = 30)))
    return(thres)
}