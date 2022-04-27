#' Find the marker gene set for each cluster
#' With an input SingleCellExperiment object and specifying the clustering
#' labels, this function iteratively call the differential expression analysis
#' on each cluster against all the others.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param method A single character for specific differential expression
#' analysis method. Choose from \code{'wilcox'}, \code{'MAST'}, \code{'DESeq2'},
#' \code{'Limma'}, and \code{'ANOVA'}. Default \code{"wilcox"}.
#' @param cluster One single character to specify a column in
#' \code{colData(inSCE)} for the clustering label. Alternatively, a vector or
#' a factor is also acceptable. Default \code{"cluster"}.
#' @param covariates A character vector of additional covariates to use when
#' building the model. All covariates must exist in
#' \code{names(colData(inSCE))}. Not applicable when \code{method} is
#' \code{"MAST"} method. Default \code{NULL}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this
#' value. Default \code{1}
#' @param minClustExprPerc A numeric scalar. The minimum cutoff of the
#' percentage of cells in the cluster of interests that expressed the marker
#' gene. Default \code{0.7}.
#' @param maxCtrlExprPerc A numeric scalar. The maximum cutoff of the
#' percentage of cells out of the cluster (control group) that expressed the
#' marker gene. Default \code{0.4}.
#' @param minMeanExpr A numeric scalar. The minimum cutoff of the mean
#' expression value of the marker in the cluster of interests. Default \code{1}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$findMarker} updated with a data.table of the up-
#' regulated DEGs for each cluster.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- findMarkerDiffExp(mouseBrainSubsetSCE,
#'                                          useAssay = "logcounts",
#'                                          cluster = "level1class")
findMarkerDiffExp <- function(inSCE, useAssay = 'logcounts',
                              method = c('wilcox', 'MAST', "DESeq2", "Limma",
                                         "ANOVA"),
                              cluster = 'cluster', covariates = NULL,
                              log2fcThreshold = 0.25, fdrThreshold = 0.05,
                              minClustExprPerc = 0.6, maxCtrlExprPerc = 0.4,
                              minMeanExpr = 0.5){
  # Input checks
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop('"inSCE" should be a SingleCellExperiment inherited Object.')
  }
  if(!useAssay %in% expDataNames(inSCE)){
    stop('"useAssay" name: ', useAssay, ' not found.')
  }
  method <- match.arg(method)
  message("Running with ", method)
  ## Check whether use colData or customized vector/factor
  if(is.character(cluster) && length(cluster) == 1){
    if(!cluster %in% names(SummarizedExperiment::colData(inSCE))){
      stop('"cluster": ', cluster, ' not found in colData(inSCE).')
    }
    clusterName <- cluster
    cluster <- SummarizedExperiment::colData(inSCE)[[cluster]]
  } else if(!length(cluster) == ncol(inSCE)){
    stop('The number of cluster labels does not match the number of cells.')
  } else {
    clusterName <- 'findMarker_cluster'
    SummarizedExperiment::colData(inSCE)[[clusterName]] <- cluster
  }
  # Iterate
  if(is.factor(cluster)){
    # In case inSCE is a subset, when "levels" is a full list of all
    # clusters, want to keep the ordering stored in the factor.
    uniqClust <- as.vector(unique(cluster))
    levelOrder <- levels(cluster)
    uniqClust <- levelOrder[levelOrder %in% uniqClust]
  } else {
    uniqClust <- unique(cluster)
  }
  for(c in uniqClust){
    message('Computing for cluster: ', c)
    clusterIndex <- cluster == c
    inSCE <- runDEAnalysis(method = method, inSCE = inSCE,
                           useAssay = useAssay,
                           index1 = clusterIndex,
                           analysisName = paste0('findMarker', c),
                           onlyPos = TRUE,
                           log2fcThreshold = log2fcThreshold,
                           fdrThreshold = fdrThreshold,
                           covariates = covariates,
                           groupName1 = c, groupName2 = 'others')
  }
  degFull <- NULL
  for(c in uniqClust){
    degTable <-
      S4Vectors::metadata(inSCE)$diffExp[[paste0('findMarker', c)]]$result
    degTable <- degTable[,c(-5,-6,-7,-8)]
    degTable$Gene <- as.character(degTable$Gene)
    if (nrow(degTable) > 0) {
      degTable[[clusterName]] <- c
      if(is.null(degFull)){
        degFull <- degTable
      } else {
        degFull <- rbind(degFull, degTable)
      }
    }
  }
  S4Vectors::metadata(inSCE)$diffExp[paste0('findMarker', uniqClust)] <- NULL
  if(length(names(S4Vectors::metadata(inSCE)$diffExp)) == 0){
    S4Vectors::metadata(inSCE)$diffExp <- NULL
  }
  degFull <- degFull[stats::complete.cases(degFull),]
  attr(degFull, "useAssay") <- useAssay
  degFull <- .calcMarkerExpr(degFull, inSCE, clusterName)
  if (!is.null(minClustExprPerc)) {
    degFull <- degFull[degFull$clusterExprPerc > minClustExprPerc,]
  }
  if (!is.null(maxCtrlExprPerc)) {
    degFull <- degFull[degFull$ControlExprPerc < maxCtrlExprPerc,]
  }
  if (!is.null(minMeanExpr)) {
    degFull <- degFull[degFull$clusterAveExpr > minMeanExpr,]
  }
  attr(degFull, "method") <- method
  attr(degFull, "params") <- list(log2fcThreshold = log2fcThreshold,
                                  fdrThreshold = fdrThreshold,
                                  minClustExprPerc = minClustExprPerc,
                                  maxCtrlExprPerc = maxCtrlExprPerc,
                                  minMeanExpr = minMeanExpr)
  S4Vectors::metadata(inSCE)$findMarker <- degFull
  return(inSCE)
}

.calcMarkerExpr <- function(markerTable, inSCE, clusterName) {
  markerTable <- markerTable[markerTable$Gene %in% rownames(inSCE),]
  genes <- markerTable$Gene
  uniqClust <- unique(markerTable[[clusterName]])
  useAssay <- attr(markerTable, "useAssay")
  cells.in.col <- rep(NA, length(genes))
  cells.out.col <- rep(NA, length(genes))
  exprs.avg <- rep(NA, length(genes))
  for (i in uniqClust) {
    markers <- genes[markerTable[[clusterName]] == i]
    cells.ix <- inSCE[[clusterName]] == i
    assay.in <- expData(inSCE, useAssay)[markers, cells.ix, drop = FALSE]
    assay.out <- expData(inSCE, useAssay)[markers, !cells.ix, drop = FALSE]
    assay.in.expr.perc <- rowSums(assay.in > 0) / ncol(assay.in)
    cells.in.col[markerTable[[clusterName]] == i] <- assay.in.expr.perc
    assay.out.expr.perc <- rowSums(assay.out > 0) / ncol(assay.out)
    cells.out.col[markerTable[[clusterName]] == i] <- assay.out.expr.perc
    assay.in.mean <- rowMeans(assay.in)
    exprs.avg[markerTable[[clusterName]] == i] <- assay.in.mean
  }
  markerTable$clusterExprPerc <- cells.in.col
  markerTable$ControlExprPerc <- cells.out.col
  markerTable$clusterAveExpr <- exprs.avg
  return(markerTable)
}


#' Fetch the table of top markers that pass the filtering
#'
#' @details Users have to run \code{findMarkerDiffExp()} prior to using this
#' function to extract a top marker table.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param log2fcThreshold Only use DEGs with the absolute values of log2FC
#' larger than this value. Default \code{1}
#' @param fdrThreshold Only use DEGs with FDR value smaller than this value.
#' Default \code{0.05}
#' @param minClustExprPerc A numeric scalar. The minimum cutoff of the
#' percentage of cells in the cluster of interests that expressed the marker
#' gene. Default \code{0.7}.
#' @param maxCtrlExprPerc A numeric scalar. The maximum cutoff of the
#' percentage of cells out of the cluster (control group) that expressed the
#' marker gene. Default \code{0.4}.
#' @param minMeanExpr A numeric scalar. The minimum cutoff of the mean
#' expression value of the marker in the cluster of interests. Default \code{1}.
#' @param topN An integer. Only to fetch this number of top markers for each
#' cluster in maximum, in terms of log2FC value. Use \code{NULL} to cancel the
#' top N subscription. Default \code{10}.
#' @return An organized \code{data.frame} object, with the top marker gene
#' information.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- findMarkerDiffExp(mouseBrainSubsetSCE,
#'                                          useAssay = "logcounts",
#'                                          cluster = "level1class")
#' findMarkerTopTable(mouseBrainSubsetSCE)
findMarkerTopTable <- function(inSCE, log2fcThreshold = 1,
                               fdrThreshold = 0.05, minClustExprPerc = 0.7,
                               maxCtrlExprPerc = 0.4, minMeanExpr = 1,
                               topN = 10) {
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop('"inSCE" should be a SingleCellExperiment inherited Object.')
  }
  if(!'findMarker' %in% names(S4Vectors::metadata(inSCE))){
    stop('"findMarker" result not found in metadata. ',
         'Run findMarkerDiffExp() in advance')
  }

  # Extract and basic filter
  degFull <- S4Vectors::metadata(inSCE)$findMarker
  if(!all(c("Gene", "Pvalue", "Log2_FC", "FDR") %in%
          colnames(degFull)[seq_len(4)])){
    stop('"findMarker" result cannot be interpreted properly')
  }
  degFull$Gene <- as.character(degFull$Gene)
  if(length(which(!degFull$Gene %in% rownames(inSCE))) > 0){
    # Remove genes happen in deg table but not in sce. Weird.
    degFull <- degFull[-which(!degFull$Gene %in% rownames(inSCE)),]
  }
  if(!is.null(log2fcThreshold)){
    degFull <- degFull[degFull$Log2_FC > log2fcThreshold,]
  }
  if(!is.null(fdrThreshold)){
    degFull <- degFull[degFull$FDR < fdrThreshold,]
  }
  if (!is.null(minClustExprPerc)) {
    degFull <- degFull[degFull$clusterExprPerc >= minClustExprPerc,]
  }
  if (!is.null(maxCtrlExprPerc)) {
    degFull <- degFull[degFull$ControlExprPerc <= maxCtrlExprPerc,]
  }
  if (!is.null(minMeanExpr)) {
    degFull <- degFull[degFull$clusterAveExpr >= minMeanExpr,]
  }
  # Remove duplicate by assigning the duplicated genes to the cluster where
  # their log2FC is the highest.
  # Done by keeping all unique genes and the rows with highest Log2FC entry
  # for each duplicated gene.
  dup.gene <- unique(degFull$Gene[duplicated(degFull$Gene)])
  for(g in dup.gene){
    deg.gix <- degFull$Gene == g
    deg.gtable <- degFull[deg.gix,]
    toKeep <- which.max(deg.gtable$Log2_FC)
    toRemove <- which(deg.gix)[-toKeep]
    degFull <- degFull[-toRemove,]
  }
  clusterName <- colnames(degFull)[5]
  selected <- character()
  if (!is.null(topN)) {
    for (c in unique(degFull[[clusterName]])) {
      deg.cluster <- degFull[degFull[[clusterName]] == c,]
      deg.cluster <- deg.cluster[order(deg.cluster$Log2_FC,
                                       decreasing = TRUE),]
      if (dim(deg.cluster)[1] > topN) {
        deg.cluster <- deg.cluster[seq_len(topN),]
      }
      selected <- c(selected, deg.cluster$Gene)
    }
  } else {
    selected <- degFull$Gene
  }
  degFull <- degFull[degFull$Gene %in% selected,]
  return(degFull)
}

#' Plot a heatmap to visualize the result of \code{\link{findMarkerDiffExp}}
#' @description This function will first reads the result saved in
#' \code{metadata} slot, named by \code{"findMarker"} and generated by
#' \code{\link{findMarkerDiffExp}}. Then it do the filtering on the statistics
#' based on the input parameters and get unique genes to plot. We choose the
#' genes that are identified as up-regulated only. As for the genes identified
#' as up-regulated for multiple clusters, we only keep the belonging towards the
#' one they have the highest Log2FC value.
#' In the heatmap, there will always be a cell annotation for the cluster
#' labeling used when finding the markers, and a feature annotation for which
#' cluster each gene belongs to. And by default we split the heatmap by these
#' two annotations. Additional legends can be added and the splitting can be
#' canceled.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param log2fcThreshold Only use DEGs with the absolute values of log2FC
#' larger than this value. Default \code{1}
#' @param fdrThreshold Only use DEGs with FDR value smaller than this value.
#' Default \code{0.05}
#' @param minClustExprPerc A numeric scalar. The minimum cutoff of the
#' percentage of cells in the cluster of interests that expressed the marker
#' gene. Default \code{0.7}.
#' @param maxCtrlExprPerc A numeric scalar. The maximum cutoff of the
#' percentage of cells out of the cluster (control group) that expressed the
#' marker gene. Default \code{0.4}.
#' @param minMeanExpr A numeric scalar. The minimum cutoff of the mean
#' expression value of the marker in the cluster of interests. Default \code{1}.
#' @param topN An integer. Only to plot this number of top markers for each
#' cluster in maximum, in terms of log2FC value. Use \code{NULL} to cancel the
#' top N subscription. Default \code{10}.
#' @param orderBy The ordering method of the clusters on the splitted heatmap.
#' Can be chosen from \code{"size"} or \code{"name"}, specified with vector of
#' ordered unique cluster labels, or set as \code{NULL} for unsplitted heatmap.
#' Default \code{"size"}.
#' @param decreasing Order the cluster decreasingly. Default \code{TRUE}.
#' @param rowDataName character. The column name(s) in \code{rowData} that need
#' to be added to the annotation. Default \code{NULL}.
#' @param colDataName character. The column name(s) in \code{colData} that need
#' to be added to the annotation. Default \code{NULL}.
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
#' @param colSplitBy character vector. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{colDataName} or
#' \code{names(cellAnnotations)}. Default is the value of \code{cluster} in
#' \code{\link{findMarkerDiffExp}} when \code{orderBy} is not \code{NULL}, or
#' \code{NULL} otherwise.
#' @param rowSplitBy character vector. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{rowDataName} or
#' \code{names(featureAnnotations)}. Default \code{"marker"}, which indicates an
#' auto generated annotation for this plot.
#' @param ... Other arguments passed to \code{\link{plotSCEHeatmap}}.
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object
#' @author Yichen Wang
#' @export
#' @examples
#' data("sceBatches")
#' logcounts(sceBatches) <- log(counts(sceBatches) + 1)
#' sce.w <- subsetSCECols(sceBatches, colData = "batch == 'w'")
#' sce.w <- findMarkerDiffExp(sce.w, method = "wilcox", cluster = "cell_type")
#' plotMarkerDiffExp(sce.w)
plotMarkerDiffExp <- function(inSCE, orderBy = 'size',
                              log2fcThreshold = 1, fdrThreshold = 0.05, minClustExprPerc = 0.7,
                              maxCtrlExprPerc = 0.4, minMeanExpr = 1, topN = 10, decreasing = TRUE,
                              rowDataName = NULL, colDataName = NULL,
                              featureAnnotationColor = NULL,
                              colSplitBy = ifelse(is.null(orderBy), NULL,
                                                  colnames(inSCE@metadata$findMarker)[5]),
                              rowSplitBy = "marker", ...){
  if(!is.null(orderBy)){
    if(length(orderBy) == 1){
      if(!orderBy %in% c('size', 'name')){
        stop('Single charater setting for "orderBy" should be chosen',
             'from "size" or "name".')
      }
    }# else if(any(!SummarizedExperiment::colData(inSCE)[[cluster]] %in%
    #             orderBy)){
    #   stop('Invalid "orderBy", please input a vector of unique ordered ',
    #        'cluster identifiers that match all clusters in ',
    #        'colData(inSCE) specified by "cluster" to adjust the order ',
    #        'of clusters.')
    #}
  }
  degFull <- findMarkerTopTable(inSCE = inSCE, topN = topN,
                                log2fcThreshold = log2fcThreshold,
                                fdrThreshold = fdrThreshold,
                                minClustExprPerc = minClustExprPerc,
                                maxCtrlExprPerc = maxCtrlExprPerc,
                                minMeanExpr = minMeanExpr)
  useAssay <- attr(degFull, "useAssay")
  inSCE <- inSCE[degFull$Gene,]
  clusterName <- colnames(degFull)[5]
  z <- SummarizedExperiment::colData(inSCE)[[clusterName]]
  if(is.factor(z)){
    z.order <- levels(z)
    # When 'z' is a subset factor, there would be redundant levels.
    z.order <- z.order[z.order %in% as.vector(unique(z))]
    z <- factor(z, levels = z.order)
    SummarizedExperiment::rowData(inSCE)[[clusterName]] <-
      factor(degFull[[clusterName]], levels = z.order)
  } else {
    SummarizedExperiment::rowData(inSCE)[[clusterName]] <-
      degFull[[clusterName]]
  }
  y <- SummarizedExperiment::rowData(inSCE)[[clusterName]]
  if(!is.null(orderBy)){
    if(length(orderBy) == 1 && orderBy == 'size'){
      z.order <- names(table(z)[order(table(z), decreasing = decreasing)])
    } else if(length(orderBy) == 1 && orderBy == 'name'){
      if(is.factor(z)){
        z.order <- sort(levels(z), decreasing = decreasing)
      } else {
        z.order <- sort(unique(z), decreasing = decreasing)
      }
    } else {
      z.order <- orderBy[-which(!orderBy %in% z)]
    }
  }
  SummarizedExperiment::colData(inSCE)[[clusterName]] <-
    factor(z, levels = z.order)
  SummarizedExperiment::rowData(inSCE)[["marker"]] <-
    factor(y, levels = z.order)
  # Organize plotSCEHeatmap arguments
  colDataName <- c(colDataName, clusterName)
  rowDataName <- c(rowDataName, "marker")
  markerConsistColor <-
    list(marker = dataAnnotationColor(inSCE, 'col')[[clusterName]])
  featureAnnotationColor <- c(featureAnnotationColor, markerConsistColor)
  hm <- plotSCEHeatmap(inSCE, useAssay = useAssay, colDataName = colDataName,
                       rowDataName = rowDataName, colSplitBy = colSplitBy,
                       rowSplitBy = rowSplitBy,
                       featureAnnotationColor = featureAnnotationColor,
                       cluster_row_slices = FALSE,
                       cluster_column_slices = FALSE, ...)
  return(hm)
}
