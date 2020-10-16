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

#' plot the violin plot to show visualize the expression distribution of DEGs
#' identified by differential expression analysis
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' \code{runMAST()} has to be run in advance.
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
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return A ggplot object of violin plot
#' @export
plotDEGViolin <- function(inSCE, useResult, threshP = FALSE, labelBy = NULL,
                          nrow = 6, ncol = 6, defaultTheme = TRUE,
                          check_sanity = TRUE){
  #TODO: DO we split the up/down regulation too?
  # Check
  .checkDiffExpResultExists(inSCE, useResult, labelBy)
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  deg <- result$result
  useAssay <- result$useAssay
  deg <- deg[order(deg$FDR),]
  geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), "Gene"]
  groupName1 <- result$groupNames[1]
  ix1 <- result$select$ix1
  cells1 <- colnames(inSCE)[ix1]
  groupName2 <- result$groupNames[2]
  ix2 <- result$select$ix2
  cells2 <- colnames(inSCE)[ix2]
  if(!is.null(labelBy)){
    replGeneName <- SummarizedExperiment::rowData(inSCE[geneToPlot,])[[labelBy]]
  } else {
    replGeneName <- geneToPlot
  }
  expres <- SummarizedExperiment::assay(inSCE[geneToPlot, c(cells1, cells2)],
                                        useAssay)
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
  names(flatDat)[5] <- useAssay
  # Plot
  violinplot <- ggplot2::ggplot(flatDat,
                                ggplot2::aes_string(x = 'condition',
                                                    y = useAssay,
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

#' plot the linear regression to show visualize the expression the of DEGs
#' identified by differential expression analysis
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' \code{runMAST()} has to be run in advance.
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
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return A ggplot object of linear regression
#' @export
plotDEGRegression <- function(inSCE, useResult, threshP = FALSE, labelBy = NULL,
                              nrow = 6, ncol = 6, defaultTheme = TRUE,
                              check_sanity = TRUE){
  #TODO: DO we split the up/down regulation too?
  # Check
  .checkDiffExpResultExists(inSCE, useResult, labelBy)
  # Extract
  result <- S4Vectors::metadata(inSCE)$diffExp[[useResult]]
  deg <- result$result
  useAssay <- result$useAssay
  geneToPlot <- deg[seq_len(min(nrow(deg), nrow*ncol)), "Gene"]
  groupName1 <- result$groupNames[1]
  ix1 <- result$select$ix1
  cells1 <- colnames(inSCE)[ix1]
  groupName2 <- result$groupNames[2]
  ix2 <- result$select$ix2
  cells2 <- colnames(inSCE)[ix2]
  if(!is.null(labelBy)){
    replGeneName <- SummarizedExperiment::rowData(inSCE[geneToPlot,])[[labelBy]]
  } else {
    replGeneName <- geneToPlot
  }
  expres <- SummarizedExperiment::assay(inSCE[geneToPlot, c(cells1, cells2)],
                                        useAssay)
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
    ggplot2::ggtitle(paste0("Linear Model Plot for ",
                            useResult))
  if(isTRUE(defaultTheme)){
    regressionplot <- .ggSCTKTheme(regressionplot)
  }
  return(regressionplot)
}

#' Heatmap visualization of DEG result
#'
#' A differential expression analysis function has to be run in advance so that
#' information is stored in the metadata of the input SCE object. This function
#' wraps plotSCEHeatmap.
#' A feature annotation basing on the log2FC level called \code{"regulation"}
#' will be automatically added. A cell annotation basing on the condition
#' selection while running the analysis called \code{"condition"}, and the
#' annotations used from \code{colData(inSCE)} while setting the condition and
#' covariates will also be added.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' \code{runMAST()} has to be run in advance.
#' @param useResult character. A string specifying the \code{analysisName}
#' used when running a differential expression analysis function.
#' @param onlyPos logical. Whether to only plot DEG with positive log2_FC
#' value. Default \code{FALSE}.
#' @param log2fcThreshold numeric. Only plot DEGs with the absolute values of
#' log2FC larger than this value. Default \code{1}.
#' @param fdrThreshold numeric. Only plot DEGs with FDR value smaller than this
#' value. Default \code{0.05}.
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
#' \code{"MAST Result: <useResult>"}.
#' @param ... Other arguments passed to \code{\link{plotSCEHeatmap}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDEAnalysis(inSCE = sce, groupName1 = "Sample1", method = "DESeq2",
#'  groupName2 = "Sample2", index1 = 1:100, index2 = 101:190, analysisName = "DESeq2")
#' plotDEGHeatmap(sce, useResult = "DESeq2", fdrThreshold = 1)
#' }
#'
#' @return A ComplexHeatmap::Heatmap object
#' @export
#' @author Yichen Wang
plotDEGHeatmap <- function(inSCE, useResult, onlyPos = FALSE,
                           log2fcThreshold = 0.25, fdrThreshold = 0.05,
                           useAssay = NULL, featureAnnotations = NULL,
                           cellAnnotations = NULL,
                           featureAnnotationColor = NULL,
                           cellAnnotationColor = NULL,
                           rowDataName = NULL, colDataName = NULL,
                           colSplitBy = 'condition', rowSplitBy = 'regulation',
                           title = paste0("MAST Result: ", useResult), ...){
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
  deg <- result$result
  if(is.null(useAssay)){
    useAssay <- result$useAssay
  } else {
    if(useAssay != result$useAssay){
      warning("`useAssay` is different to the one used for `runMAST()`")
    }
  }
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
  allGenes <- rownames(inSCE)[gene.ix]
  allCells <- colnames(inSCE)[cell.ix]

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
  regulation[rownames(inSCE) %in% genes.up] <- 'up'
  regulation[rownames(inSCE) %in% genes.down] <- 'down'
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
  hm <- plotSCEHeatmap(inSCE = inSCE, useAssay = useAssay,
                       featureIndex = gene.ix, cellIndex = cell.ix,
                       featureAnnotations = featureAnnotations,
                       cellAnnotations = cellAnnotations,
                       rowDataName = rowDataName,
                       colDataName = colDataName,
                       featureAnnotationColor = featureAnnotationColor,
                       cellAnnotationColor = cellAnnotationColor,
                       rowSplitBy = rowSplitBy, colSplitBy = colSplitBy,
                       title = title)
  return(hm)
}

#' MAST Identify adaptive thresholds
#'
#' Calculate and produce a list of thresholded counts (on natural scale),
#' thresholds, bins, densities estimated on each bin, and the original data from
#' \code{\link[MAST]{thresholdSCRNACountMatrix}}
#' @param inSCE SingleCellExperiment object
#' @param useAssay character, default \code{"logcounts"}
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}
#' @return Plot the thresholding onto the plotting region.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotMASTThresholdGenes(mouseBrainSubsetSCE)
plotMASTThresholdGenes <- function(inSCE, useAssay="logcounts",
                                   check_sanity = TRUE){
  # data preparation
  expres <- SummarizedExperiment::assay(inSCE, useAssay)
  if(!is.matrix(expres)){
    expres <- as.matrix(expres)
  }
  expres <- featureNameDedup(expres)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                             fdata, check_sanity = check_sanity)
  SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
  invisible(utils::capture.output(thres <- MAST::thresholdSCRNACountMatrix(
    SummarizedExperiment::assay(SCENew), nbins = 20, min_per_bin = 30)))
  # plotting
  graphics::par(mfrow = c(5, 4))
  graphics::plot(thres)
  graphics::par(mfrow = c(1, 1))
  # return(thres)
}
