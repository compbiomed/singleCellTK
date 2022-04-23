#' Check if the specified MAST result in SingleCellExperiment object is
#' complete. But does not garantee the biological correctness.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. a
#' differential expression analysis function has to be run in advance.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param labelBy A single character for a column of \code{rowData(inSCE)} as
#' where to search for the labeling text. Default \code{NULL}.
#' @return Stop point if found
.checkDiffExpResultExists <- function(inSCE, useResult, labelBy = NULL){
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop('Given object is not a valid SingleCellExperiment object.')
  }
  if(!'diffExp' %in% names(S4Vectors::metadata(inSCE))){
    stop('"diffExp" not in metadata, please run runMAST() first.')
  }
  if(!useResult %in% names(S4Vectors::metadata(inSCE)$diffExp)){
    stop(paste0('"', useResult, '"', ' not in metadata(inSCE)$diffExp. '),
         'Please check.')
  }
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  if(!all(c('groupNames', 'select', 'result', 'useAssay') %in% names(result))){
    stop(paste0('"', useResult, '"', ' result is not complete. '),
         'You might need to rerun it.')
  }
  if(!is.null(labelBy)){
    if(!labelBy %in% names(SummarizedExperiment::rowData(inSCE))){
      stop("labelBy: '", labelBy, "' not found.")
    }
  }
}

#' Generate violin plot to show the expression of top DEGs
#' @details Any of the differential expression analysis method from SCTK should 
#' be performed prior to using this function
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param threshP logical. Whether to plot threshold values from adaptive
#' thresholding, instead of using the assay used by \code{runMAST()}. Default
#' \code{FALSE}.
#' @param labelBy A single character for a column of \code{rowData(inSCE)} as
#' where to search for the labeling text. Default \code{NULL}.
#' @param nrow Integer. Number of rows in the plot grid. Default \code{6}.
#' @param ncol Integer. Number of columns in the plot grid. Default \code{6}.
#' @param defaultTheme Logical scalar. Whether to use default SCTK theme in
#' ggplot. Default \code{TRUE}.
#' @param isLogged Logical scalar. Whether the assay used for the analysis is
#' logged. If not, will do a \code{log(assay + 1)} transformation. Default
#' \code{TRUE}.
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return A ggplot object of violin plot
#' @export
#' @examples
#' data("sceBatches")
#' logcounts(sceBatches) <- log(counts(sceBatches) + 1)
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- runWilcox(sce.w, class = "cell_type", classGroup1 = "alpha",
#'                    groupName1 = "w.alpha", groupName2 = "w.beta",
#'                    analysisName = "w.aVSb")
#' plotDEGViolin(sce.w, "w.aVSb")
plotDEGViolin <- function(inSCE, useResult, threshP = FALSE, labelBy = NULL,
                          nrow = 6, ncol = 6, defaultTheme = TRUE,
                          isLogged = TRUE, check_sanity = TRUE){
  #TODO: DO we split the up/down regulation too?
  # Check
  .checkDiffExpResultExists(inSCE, useResult, labelBy)
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  deg <- result$result
  useAssay <- result$useAssay
  useReducedDim <- result$useReducedDim
  deg <- deg[order(deg$FDR),]
  geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), "Gene"]
  groupName1 <- result$groupNames[1]
  ix1 <- result$select$ix1
  cells1 <- colnames(inSCE)[ix1]
  groupName2 <- result$groupNames[2]
  ix2 <- result$select$ix2
  cells2 <- colnames(inSCE)[ix2]
  if (!is.null(useAssay)) {
    if(!is.null(labelBy)){
      replGeneName <- SummarizedExperiment::rowData(inSCE[geneToPlot,])[[labelBy]]
    } else {
      replGeneName <- geneToPlot
    }
    expres <- expData(inSCE[geneToPlot, c(cells1, cells2)],
                      useAssay)
    useMat <- useAssay
  } else {
    if(!is.null(labelBy)){
      warning("Analysis performed on reducedDim. Cannot use rowData for ", 
              "labelBy. Ignored.")
    }
    replGeneName <- geneToPlot
    expres <- t(expData(inSCE[, c(cells1, cells2)], useReducedDim))[geneToPlot,]
    useMat <- useReducedDim
  }
  
  if(!is.matrix(expres)){
    expres <- as.matrix(expres)
  }
  rownames(expres) <- replGeneName
  # Format
  cdat <- data.frame(wellKey = colnames(expres),
                     condition = factor(c(rep(groupName1, length(cells1)),
                                          rep(groupName2, length(cells2))),
                                        levels = result$groupNames),
                     ngeneson = rep("", (length(cells1) + length(cells2))),
                     stringsAsFactors = FALSE)
  if (!isTRUE(isLogged)) {
    expres <- log(expres + 1)
  }
  sca <- suppressMessages(MAST::FromMatrix(expres, cdat,
                                           check_sanity = check_sanity))
  if(threshP){
    #TODO: if nrow*ncol < `min_per_bin`` below, there would be an error.
    invisible(utils::capture.output(thres <-
                                      MAST::thresholdSCRNACountMatrix(expres, nbins = 20,
                                                                      min_per_bin = 30)))
    SummarizedExperiment::assay(sca) <- thres$counts_threshold
  }
  flatDat <- methods::as(sca, "data.table")
  flatDat$primerid <- factor(flatDat$primerid, levels = replGeneName)
  names(flatDat)[5] <- useMat
  # Plot
  violinplot <- ggplot2::ggplot(flatDat,
                                ggplot2::aes_string(x = 'condition',
                                                    y = useMat,
                                                    color = 'condition')) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y",
                        ncol = ncol) +
    ggplot2::geom_violin() +
    ggplot2::ggtitle(paste0("Violin Plot for ", useResult))
  if(isTRUE(defaultTheme)){
    violinplot <- .ggSCTKTheme(violinplot)
  }
  return(violinplot)
}

#' Create linear regression plot to show the expression the of top DEGs
#' @details Any of the differential expression analysis method from SCTK should 
#' be performed prior to using this function
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param threshP logical. Whether to plot threshold values from adaptive
#' thresholding, instead of using the assay used by when performing DE analysis.
#' Default \code{FALSE}.
#' @param labelBy A single character for a column of \code{rowData(inSCE)} as
#' where to search for the labeling text. Default \code{NULL}.
#' @param nrow Integer. Number of rows in the plot grid. Default \code{6}.
#' @param ncol Integer. Number of columns in the plot grid. Default \code{6}.
#' @param defaultTheme Logical scalar. Whether to use default SCTK theme in
#' ggplot. Default \code{TRUE}.
#' @param isLogged Logical scalar. Whether the assay used for the analysis is
#' logged. If not, will do a \code{log(assay + 1)} transformation. Default
#' \code{TRUE}.
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return A ggplot object of linear regression
#' @export
#' @examples
#' data("sceBatches")
#' logcounts(sceBatches) <- log(counts(sceBatches) + 1)
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- runWilcox(sce.w, class = "cell_type", classGroup1 = "alpha",
#'                    groupName1 = "w.alpha", groupName2 = "w.beta",
#'                    analysisName = "w.aVSb")
#' plotDEGRegression(sce.w, "w.aVSb")
plotDEGRegression <- function(inSCE, useResult, threshP = FALSE, labelBy = NULL,
                              nrow = 6, ncol = 6, defaultTheme = TRUE,
                              isLogged = TRUE, check_sanity = TRUE){
  #TODO: DO we split the up/down regulation too?
  # Check
  .checkDiffExpResultExists(inSCE, useResult, labelBy)
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  deg <- result$result
  useReducedDim <- result$useReducedDim
  useAssay <- result$useAssay
  geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), "Gene"]
  groupName1 <- result$groupNames[1]
  ix1 <- result$select$ix1
  cells1 <- colnames(inSCE)[ix1]
  groupName2 <- result$groupNames[2]
  ix2 <- result$select$ix2
  cells2 <- colnames(inSCE)[ix2]
  if (!is.null(useAssay)) {
    if(!is.null(labelBy)){
      replGeneName <- SummarizedExperiment::rowData(inSCE[geneToPlot,])[[labelBy]]
    } else {
      replGeneName <- geneToPlot
    }
    expres <- expData(inSCE[geneToPlot, c(cells1, cells2)],
                      useAssay)
    useMat <- useAssay
  } else {
    if (!is.null(labelBy)) {
      warning("Analysis performed on reducedDim. Cannot use rowData for ", 
              "labelBy. Ignored.")
    }
    replGeneName <- geneToPlot
    expres <- t(expData(inSCE[,c(cells1, cells2)], useReducedDim))[geneToPlot,]
    useMat <- useReducedDim
  }
  if(!is.matrix(expres)){
    expres <- as.matrix(expres)
  }
  rownames(expres) <- replGeneName
  # Format
  cdat <- data.frame(wellKey = colnames(expres),
                     condition = factor(c(rep(groupName1, length(cells1)),
                                          rep(groupName2, length(cells2))),
                                        levels = result$groupNames),
                     ngeneson = rep("", (length(cells1) + length(cells2))),
                     stringsAsFactors = FALSE)
  if (!isTRUE(isLogged)) {
    expres <- log(expres + 1)
  }
  sca <- suppressMessages(MAST::FromMatrix(expres, cData = cdat,
                                           check_sanity = check_sanity))
  cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
  SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
  if(length(unique(SummarizedExperiment::colData(sca)$cngeneson)) < 2){
    stop("Standardized cellular detection rate not various, unable to plot.")
  }
  if(threshP){
    #TODO: if nrow*ncol < `min_per_bin`` below, there would be an error.
    invisible(utils::capture.output(thres <-
                                      MAST::thresholdSCRNACountMatrix(expres,
                                                                      nbins = 20,
                                                                      min_per_bin = 30)))
    SummarizedExperiment::assay(sca) <- thres$counts_threshold
  }
  flatDat <- methods::as(sca, "data.table")
  flatDat$primerid <- factor(flatDat$primerid, levels = replGeneName)
  names(flatDat)[6] <- useMat
  # Calculate
  resData <- NULL
  for (i in unique(flatDat$primerid)){
    resdf <- flatDat[flatDat$primerid == i, ]
    resdf$lmPred <- stats::lm(
      stats::as.formula(paste0(useMat, "~cngeneson+", 'condition')),
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
    y = useMat,
    color = 'condition')) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y",
                        ncol = ncol)
  regressionplot <- ggbase +
    ggplot2::aes_string(x = "cngeneson") +
    ggplot2::geom_line(ggplot2::aes_string(y = "lmPred"),
                       lty = 1) +
    ggplot2::xlab("Standardized Cellular Detection Rate") +
    ggplot2::ggtitle(paste0("Linear Model Plot for ",
                            useResult))
  if(isTRUE(defaultTheme)){
    regressionplot <- .ggSCTKTheme(regressionplot)
  }
  return(regressionplot)
}

#' Get Top Table of a DEG analysis
#' @description Users have to run \code{runDEAnalysis()} first, any of the 
#' wrapped functions of this generic function. Users can set further filters on
#' the result. A \code{data.frame} object, with variables of \code{Gene}, 
#' \code{Log2_FC}, \code{Pvalue}, and \code{FDR}, will be returned. 
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object, with of the
#' singleCellTK DEG method performed in advance. 
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param labelBy A single character for a column of \code{rowData(inSCE)} as
#' where to search for the labeling text. Default \code{NULL} for the 
#' "rownames".
#' @param onlyPos logical. Whether to only fetch DEG with positive log2_FC
#' value. Default \code{FALSE}.
#' @param log2fcThreshold numeric. Only fetch DEGs with the absolute values of
#' log2FC larger than this value. Default \code{0.25}.
#' @param fdrThreshold numeric. Only fetch DEGs with FDR value smaller than this
#' value. Default \code{0.05}.
#' @param minGroup1MeanExp numeric. Only fetch DEGs with mean expression in 
#' group1 greater then this value. Default \code{NULL}.
#' @param maxGroup2MeanExp numeric. Only fetch DEGs with mean expression in 
#' group2 less then this value. Default \code{NULL}.
#' @param minGroup1ExprPerc numeric. Only fetch DEGs expressed in greater then 
#' this fraction of cells in group1. Default \code{NULL}.
#' @param maxGroup2ExprPerc numeric. Only fetch DEGs expressed in less then this 
#' fraction of cells in group2. Default \code{NULL}.
#' @return A \code{data.frame} object of the top DEGs, with variables of 
#' \code{Gene}, \code{Log2_FC}, \code{Pvalue}, and \code{FDR}.
#' @export
#' @examples
#' data("sceBatches")
#' sceBatches <- scaterlogNormCounts(sceBatches, "logcounts")
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- runWilcox(sce.w, class = "cell_type", classGroup1 = "alpha",
#'                    groupName1 = "w.alpha", groupName2 = "w.beta",
#'                    analysisName = "w.aVSb")
#' getDEGTopTable(sce.w, "w.aVSb")
getDEGTopTable <- function(inSCE, useResult, labelBy = NULL, onlyPos = FALSE,
                        log2fcThreshold = 0.25, fdrThreshold = 0.05,
                        minGroup1MeanExp = NULL, maxGroup2MeanExp = NULL, 
                        minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL){
  # Check
  .checkDiffExpResultExists(inSCE, useResult, labelBy)
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]$result
  if (!is.null(labelBy)) {
    genes <- result$Gene
    result$Gene <- SummarizedExperiment::rowData(inSCE[genes,])[[labelBy]]
  }
  # Filter
  result <- .filterDETable(result, onlyPos, log2fcThreshold, fdrThreshold, 
                           minGroup1MeanExp, maxGroup2MeanExp, 
                           minGroup1ExprPerc, maxGroup2ExprPerc)
  return(result)
}

#' Heatmap visualization of DEG result
#'
#' @details A differential expression analysis function has to be run in advance 
#' so that information is stored in the metadata of the input SCE object. This 
#' function wraps \code{\link{plotSCEHeatmap}}.
#' A feature annotation basing on the log2FC level called \code{"regulation"}
#' will be automatically added. A cell annotation basing on the condition
#' selection while running the analysis called \code{"condition"}, and the
#' annotations used from \code{colData(inSCE)} while setting the condition and
#' covariates will also be added.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param doLog Logical scalar. Whether to do \code{log(assay + 1)}
#' transformation on the assay used for the analysis. Default \code{FALSE}.
#' @param onlyPos logical. Whether to only plot DEG with positive log2_FC
#' value. Default \code{FALSE}.
#' @param log2fcThreshold numeric. Only plot DEGs with the absolute values of
#' log2FC larger than this value. Default \code{0.25}.
#' @param fdrThreshold numeric. Only plot DEGs with FDR value smaller than this
#' value. Default \code{0.05}.
#' @param minGroup1MeanExp numeric. Only plot DEGs with mean expression in 
#' group1 greater then this value. Default \code{NULL}.
#' @param maxGroup2MeanExp numeric. Only plot DEGs with mean expression in 
#' group2 less then this value. Default \code{NULL}.
#' @param minGroup1ExprPerc numeric. Only plot DEGs expressed in greater then 
#' this fraction of cells in group1. Default \code{NULL}.
#' @param maxGroup2ExprPerc numeric. Only plot DEGs expressed in less then this 
#' fraction of cells in group2. Default \code{NULL}.
#' @param useAssay character. A string specifying an assay of expression value
#' to plot. By default the assay used for \code{runMAST()} will be used.
#' Default \code{NULL}.
#' @param featureAnnotations \code{data.frame}, with \code{rownames} containing
#' all the features going to be plotted. Character columns should be factors.
#' Default \code{NULL}.
#' @param cellAnnotations \code{data.frame}, with \code{rownames} containing
#' all the cells going to be plotted. Character columns should be factors.
#' Default \code{NULL}.
#' @param featureAnnotationColor A named list. Customized color settings for
#' feature labeling. Should match the entries in the \code{featureAnnotations}
#' or \code{rowDataName}. For each entry, there should be a list/vector of
#' colors named with categories. Default \code{NULL}.
#' @param cellAnnotationColor A named list. Customized color settings for
#' cell labeling. Should match the entries in the \code{cellAnnotations} or
#' \code{colDataName}. For each entry, there should be a list/vector of colors
#' named with categories. Default \code{NULL}.
#' @param rowDataName character. The column name(s) in \code{rowData} that need
#' to be added to the annotation. Default \code{NULL}.
#' @param colDataName character. The column name(s) in \code{colData} that need
#' to be added to the annotation. Default \code{NULL}.
#' @param rowSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{rowDataName} or
#' \code{names(featureAnnotations)}. Default \code{"regulation"}.
#' @param colSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{colDataName} or
#' \code{names(cellAnnotations)}. Default \code{"condition"}.
#' @param title character. Main title of the heatmap. Default
#' \code{"DE Analysis: <useResult>"}.
#' @param ... Other arguments passed to \code{\link{plotSCEHeatmap}}
#' @examples
#' data("sceBatches")
#' logcounts(sceBatches) <- log(counts(sceBatches) + 1)
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- runWilcox(sce.w, class = "cell_type", classGroup1 = "alpha",
#'                    groupName1 = "w.alpha", groupName2 = "w.beta",
#'                    analysisName = "w.aVSb")
#' plotDEGHeatmap(sce.w, "w.aVSb")
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#' @author Yichen Wang
plotDEGHeatmap <- function(inSCE, useResult, doLog = FALSE, onlyPos = FALSE,
                           log2fcThreshold = 0.25, fdrThreshold = 0.05,
                           minGroup1MeanExp = NULL, maxGroup2MeanExp = NULL, 
                           minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL,
                           useAssay = NULL, featureAnnotations = NULL,
                           cellAnnotations = NULL,
                           featureAnnotationColor = NULL,
                           cellAnnotationColor = NULL,
                           rowDataName = NULL, colDataName = NULL,
                           colSplitBy = 'condition', rowSplitBy = 'regulation',
                           title = paste0("DE Analysis: ", useResult), ...){
  # Check
  .checkDiffExpResultExists(inSCE, useResult)
  extraArgs <- list(...)
  warnArgs <- c('featureIndex', 'cellIndex')
  if(any(warnArgs %in% names(extraArgs))){
    warning('"', paste(warnArgs[warnArgs %in% names(extraArgs)],
                       collapse = ', '), '" are not allowed at this point.')
    extraArgs[c('cellIndex', 'featureIndex')] <- NULL
  }
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  if (!is.null(result$useReducedDim)) {
    # Analysis performed with reducedDim, cannot set useAssay in plotDEGHeatmap
    if (!is.null(useAssay)) {
      warning("Analysis performed on reducedDim, cannot set `useAssay`, ", 
              "ignored.")
    }
    useAssay <- NULL
    useReducedDim <- result$useReducedDim
  } else {
    if(is.null(useAssay)){
      useAssay <- result$useAssay
    }
    useReducedDim <- NULL
  }
  
  ix1 <- result$select$ix1
  ix2 <- result$select$ix2
  deg.filtered <- getDEGTopTable(inSCE, useResult = useResult, labelBy = NULL, 
                                 onlyPos = onlyPos, 
                                 log2fcThreshold = log2fcThreshold, 
                                 fdrThreshold = fdrThreshold,
                                 minGroup1MeanExp, maxGroup2MeanExp, 
                                 minGroup1ExprPerc, maxGroup2ExprPerc)
  if(dim(deg.filtered)[1] <= 1){
    stop('Too few genes that pass filtration, unable to plot')
  }
  # Not directly using deg.filtered$Gene because rownames might be deduplicated
  # to avoid error when performing DEG.
  if (!is.null(useReducedDim)) {
    mat <- t(expData(inSCE, useReducedDim))
    assayList <- list(mat)
    names(assayList) <- useReducedDim
    tmpSCE <- SingleCellExperiment::SingleCellExperiment(assays = assayList)
    SummarizedExperiment::colData(tmpSCE) <- SummarizedExperiment::colData(inSCE)
    assayName <- useReducedDim
  } else {
    tmpSCE <- inSCE
    assayName <- useAssay
  }
  gene.ix <- rownames(tmpSCE) %in% deg.filtered$Gene
  cell.ix <- which(ix1 | ix2)
  allGenes <- rownames(tmpSCE)[gene.ix]
  allCells <- colnames(tmpSCE)[cell.ix]

  # Annotation organization
  ## Cells
  group <- vector()
  group[ix1] <- result$groupNames[1]
  group[ix2] <- result$groupNames[2]
  group <- factor(group[cell.ix], levels = result$groupNames)
  if(!is.null(cellAnnotations)){
    if(!all(allCells %in% rownames(cellAnnotations))){
      stop('Not all cells involved in comparison found in given ',
           '`cellAnnotations`. ')
    }
    cellAnnotations <- cellAnnotations[allCells, , drop = FALSE]
    cellAnnotations <- data.frame(cellAnnotations, condition = group)
  } else {
    cellAnnotations <- data.frame(condition = group,
                                  row.names = allCells)
  }
  kCol <- celda::distinctColors(2)
  names(kCol) <- result$groupNames
  if(!is.null(cellAnnotationColor)){
    if(!"condition" %in% names(cellAnnotationColor)){
      cellAnnotationColor <- c(list(condition = kCol),
                               cellAnnotationColor)
    }
  } else {
    cellAnnotationColor <- list(condition = kCol)
  }
  if(!is.null(colDataName)){
    if (length(which(colDataName %in% result$annotation)) > 0) {
      colDataName <- colDataName[-which(colDataName %in% result$annotation)]
    }
    colDataName <- c(colDataName, result$annotation)
  } else {
    colDataName <- result$annotation
  }
  ## Genes
  regulation <- vector()
  genes.up <- deg.filtered[deg.filtered$Log2_FC > 0, "Gene"]
  genes.down <- deg.filtered[deg.filtered$Log2_FC < 0, "Gene"]
  regulation[rownames(tmpSCE) %in% genes.up] <- 'up'
  regulation[rownames(tmpSCE) %in% genes.down] <- 'down'
  regulation <- factor(regulation[gene.ix], levels = c('up', 'down'))
  if(!is.null(featureAnnotations)){
    if(!all(allGenes %in% rownames(featureAnnotations))){
      stop('Not all genes involved in comparison found in given ',
           '`featureAnnotations`. ')
    }
    featureAnnotations <- featureAnnotations[allGenes, , drop = FALSE]
    featureAnnotations <- data.frame(featureAnnotations,
                                     regulation = regulation)
  } else {
    featureAnnotations <- data.frame(regulation = regulation,
                                     row.names = allGenes)
  }
  lCol <- celda::distinctColors(2)
  names(lCol) <- c('up', 'down')
  if(!is.null(featureAnnotationColor)){
    if(!"regulation" %in% names(featureAnnotationColor)){
      featureAnnotationColor <- c(list(regulation = lCol),
                                  featureAnnotationColor)
    }
  } else {
    featureAnnotationColor <- list(regulation = lCol)
  }
  # Plot
  hm <- plotSCEHeatmap(inSCE = tmpSCE, useAssay = assayName, doLog = doLog,
                       featureIndex = gene.ix, cellIndex = cell.ix,
                       featureAnnotations = featureAnnotations,
                       cellAnnotations = cellAnnotations,
                       rowDataName = rowDataName,
                       colDataName = colDataName,
                       featureAnnotationColor = featureAnnotationColor,
                       cellAnnotationColor = cellAnnotationColor,
                       rowSplitBy = rowSplitBy, colSplitBy = colSplitBy,
                       title = title, ...)
  return(hm)
}

#' Generate volcano plot for DEGs
#' @details Any of the differential expression analysis method from SCTK should 
#' be performed prior to using this function to generate volcano plots. 
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param labelTopN Integer, label this number of top DEGs that pass the 
#' filters.
#' @param log2fcThreshold numeric. Label genes with the absolute values of
#' log2FC greater than this value as regulated. Default \code{0.25}.
#' @param fdrThreshold numeric. Label genes with FDR value less than this
#' value as regulated. Default \code{0.05}.
#' @return A \code{ggplot} object of volcano plot
#' @export
#' @examples 
#' data("sceBatches")
#' sceBatches <- scaterlogNormCounts(sceBatches, "logcounts")
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- runWilcox(sce.w, class = "cell_type", classGroup1 = "alpha",
#'                    groupName1 = "w.alpha", groupName2 = "w.beta",
#'                    analysisName = "w.aVSb")
#' plotDEGVolcano(sce.w, "w.aVSb")
plotDEGVolcano <- function(inSCE,
                           useResult,
                           labelTopN = 10,
                           log2fcThreshold = 0.25, 
                           fdrThreshold = 0.05) {
  .checkDiffExpResultExists(inSCE, useResult)
  deg <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]$result
  deg <- deg[order(deg$FDR),]
  rownames(deg) <- deg$Gene
  groupNames <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]$groupNames
  # Prepare for coloring that shows the filtering
  deg$Regulation <- NA
  deg$Regulation[deg$Log2_FC > 0] <- "Up"
  deg$Regulation[deg$Log2_FC < 0] <- "Down"
  if (!is.null(log2fcThreshold)) {
    deg$Regulation[abs(deg$Log2_FC) < log2fcThreshold] <- "No"
  }
  if (!is.null(fdrThreshold)) {
    deg$Regulation[deg$FDR > fdrThreshold] <- "No"
  }
  # Prepare for Top DEG text labeling
  passIdx <- deg$Regulation != "No"
  deg$label <- NA
  labelTopN <- min(labelTopN, length(which(passIdx)))
  deg.pass <- deg[passIdx,]
  label.origTable.idx <- deg$Gene %in% deg.pass$Gene[seq(labelTopN)]
  deg$label[label.origTable.idx] <- deg$Gene[label.origTable.idx]
  # Prepare for lines that mark the cutoffs
  vlineLab <- data.frame(
    X = c(-log2fcThreshold, log2fcThreshold),
    text = c(paste("lower log2FC cutoff:", -log2fcThreshold),
             paste("upper log2FC cutoff:", log2fcThreshold)),
    h = c(1.01, -0.01)
  )
  hlineLab <- data.frame(
    Y = c(-log10(fdrThreshold)),
    text = paste("FDR cutoff:", fdrThreshold)
  )
  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_point(data = deg, 
                        ggplot2::aes_string(x = "Log2_FC", y = "-log10(FDR)", 
                                            col = "Regulation")) +
    ggplot2::scale_color_manual(values = c("Down" = "#619cff", 
                                           "No" = "light grey", 
                                           "Up" = "#f8766d")) +
    ggrepel::geom_text_repel(data = deg,
                             ggplot2::aes_string(x = "Log2_FC", 
                                                 y = "-log10(FDR)",
                                                 label = "label"),
                             colour = "black", na.rm = TRUE) +
    ggplot2::geom_vline(data = vlineLab,
                        ggplot2::aes_string(xintercept = "X"),
                        linetype = "longdash") +
    ggplot2::geom_text(data = vlineLab, 
                       ggplot2::aes_string(x = "X", y = 0, label = "text", 
                                           hjust = "h"),
                       size = 3, vjust = 1) +
    ggplot2::geom_hline(data = hlineLab,
                        ggplot2::aes_string(yintercept = "Y"),
                        linetype = "longdash") +
    ggplot2::geom_text(data = hlineLab,
                       ggplot2::aes_string(x = -Inf, y = "Y", label = "text"),
                       size = 3, vjust = -0.5, hjust = -.03) +
    ggplot2::xlab("Fold Change (log2)") +
    ggplot2::ylab("FDR (-Log10 q-value)") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::xlim(-max(abs(deg$Log2_FC)), max(abs(deg$Log2_FC))) +
    ggplot2::ggtitle(paste("DEG between", groupNames[1], 
                           "and", groupNames[2]))
}

#' MAST Identify adaptive thresholds
#'
#' Calculate and produce a list of thresholded counts (on natural scale),
#' thresholds, bins, densities estimated on each bin, and the original data from
#' \code{\link[MAST]{thresholdSCRNACountMatrix}}
#' @param inSCE SingleCellExperiment object
#' @param useAssay character, default \code{"logcounts"}
#' @param doPlot Logical scalar. Whether to directly plot in the plotting area.
#' If \code{FALSE}, will return a graphical object which can be visualized with
#' \code{grid.draw()}. Default \code{TRUE}.
#' @param isLogged Logical scalar. Whether the assay used for the analysis is
#' logged. If not, will do a \code{log(assay + 1)} transformation. Default
#' \code{TRUE}.
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return Plot the thresholding onto the plotting region if \code{plot == TRUE}
#' or a graphical object if \code{plot == FALSE}.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotMASTThresholdGenes(mouseBrainSubsetSCE)
plotMASTThresholdGenes <- function(inSCE, useAssay="logcounts", doPlot = TRUE,
                                   isLogged = TRUE, check_sanity = TRUE){
  # data preparation
  expres <- expData(inSCE, useAssay)
  if(!is.matrix(expres)){
    expres <- as.matrix(expres)
  }
  expres <- dedupRowNames(expres)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                             fdata, check_sanity = check_sanity)
  SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
  invisible(utils::capture.output(thres <- MAST::thresholdSCRNACountMatrix(
    SummarizedExperiment::assay(SCENew), nbins = 20, min_per_bin = 30,
    data_log = isLogged)))
  # plotting
  plotNRow <- ceiling(length(thres$valleys) / 4)
  thres.grob <- ggplotify::as.grob(function(){
    graphics::par(mfrow = c(plotNRow, 4), mar = c(3, 3, 2, 1),
        mgp = c(2, 0.7, 0), tck = -0.01, new = TRUE)
    plot(thres)
  })
  if (isTRUE(doPlot)) {
    grid::grid.draw(thres.grob)
  } else {
    return(thres.grob)
  }
}
