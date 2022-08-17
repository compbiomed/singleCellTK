#' Find the marker gene set for each cluster
#' @rdname runFindMarker
#' @description With an input SingleCellExperiment object and specifying the
#' clustering labels, this function iteratively call the differential expression
#' analysis on each cluster against all the others. \code{findMarkerDiffExp}
#' will be deprecated in the future.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param useReducedDim character. A string specifying which reducedDim to use
#' for MAST calculations. Set \code{useAssay} to \code{NULL} when using.
#' Required.
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
#' value. Default \code{NULL}
#' @param minClustExprPerc A numeric scalar. The minimum cutoff of the
#' percentage of cells in the cluster of interests that expressed the marker
#' gene. From 0 to 1. Default \code{NULL}.
#' @param maxCtrlExprPerc A numeric scalar. The maximum cutoff of the
#' percentage of cells out of the cluster (control group) that expressed the
#' marker gene. From 0 to 1. Default \code{NULL}.
#' @param minMeanExpr A numeric scalar. The minimum cutoff of the mean
#' expression value of the marker in the cluster of interests. Default
#' \code{NULL}.
#' @param detectThresh A numeric scalar, above which a matrix value will be
#' treated as expressed when calculating cluster/control expression percentage.
#' Default \code{0}.
#' @details
#' The returned marker table, in the \code{metadata} slot, consists of 8
#' columns: \code{"Gene"}, \code{"Log2_FC"}, \code{"Pvalue"}, \code{"FDR"},
#' \code{cluster}, \code{"clusterExprPerc"}, \code{"ControlExprPerc"} and
#' \code{"clusterAveExpr"}.
#'
#' \code{"clusterExprPerc"} is the fraction of cells,
#' that has marker value (e.g. gene expression counts) larger than
#' \code{detectThresh}, in the cell population of the cluster. As for each
#' cluster, we set all cells out of this cluster as control. Similarly,
#' \code{"ControlExprPerc"} is the fraction of cells with marker value larger
#' than \code{detectThresh} in the control cell group.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$findMarker} updated with a data.table of the up-
#' regulated DEGs for each cluster.
#' @seealso \code{\link{runDEAnalysis}}
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runFindMarker(mouseBrainSubsetSCE,
#'                                      useAssay = "logcounts",
#'                                      cluster = "level1class")
runFindMarker <- function(inSCE, useAssay = 'logcounts',
                          useReducedDim = NULL,
                          method = c('wilcox', 'MAST', "DESeq2", "Limma",
                                     "ANOVA"),
                          cluster = 'cluster', covariates = NULL,
                          log2fcThreshold = NULL, fdrThreshold = 0.05,
                          minClustExprPerc = NULL, maxCtrlExprPerc = NULL,
                          minMeanExpr = NULL, detectThresh = 0){
  method <- match.arg(method)
  # Input checks will be done in `runDEAnalysis()`
  if (is.character(cluster) && length(cluster) == 1) clusterName <- cluster
  else clusterName <- 'findMarker_cluster'
  clusterVar <- .manageCellVar(inSCE, var = cluster)
  # Iterate
  if(is.factor(clusterVar)){
    # In case inSCE is a subset, when "levels" is a full list of all
    # clusters, want to keep the ordering stored in the factor.
    uniqClust <- as.vector(unique(clusterVar))
    levelOrder <- levels(clusterVar)
    uniqClust <- levelOrder[levelOrder %in% uniqClust]
  } else {
    uniqClust <- unique(clusterVar)
  }
  for(c in uniqClust){
    message(date(), " ... Identifying markers for cluster '", c,
            "', using DE method '", method, "'")
    clusterIndex <- clusterVar == c
    inSCE <- runDEAnalysis(method = method, inSCE = inSCE,
                           useAssay = useAssay,
                           useReducedDim = useReducedDim,
                           index1 = clusterIndex,
                           analysisName = paste0('findMarker', c),
                           onlyPos = TRUE,
                           log2fcThreshold = log2fcThreshold,
                           fdrThreshold = fdrThreshold,
                           covariates = covariates,
                           groupName1 = c, groupName2 = 'others',
                           verbose = FALSE)
  }
  message(date(), " ... Organizing findMarker result")
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
  if (!is.null(useAssay)) {
    attr(degFull, "useAssay") <- useAssay
  } else {
    attr(degFull, "useReducedDim") <- useReducedDim
  }

  degFull <- .calcMarkerExpr(degFull, inSCE, clusterName, clusterVar,
                             detectThresh)
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
  attr(degFull, "cluster") <- clusterVar
  attr(degFull, "clusterName") <- clusterName
  attr(degFull, "params") <- list(log2fcThreshold = log2fcThreshold,
                                  fdrThreshold = fdrThreshold,
                                  minClustExprPerc = minClustExprPerc,
                                  maxCtrlExprPerc = maxCtrlExprPerc,
                                  minMeanExpr = minMeanExpr)
  S4Vectors::metadata(inSCE)$findMarker <- degFull
  return(inSCE)
}

#' @rdname runFindMarker
#' @export
findMarkerDiffExp <- function(inSCE, useAssay = 'logcounts',
                              useReducedDim = NULL,
                              method = c('wilcox', 'MAST', "DESeq2", "Limma",
                                         "ANOVA"),
                              cluster = 'cluster', covariates = NULL,
                              log2fcThreshold = NULL, fdrThreshold = 0.05,
                              minClustExprPerc = NULL, maxCtrlExprPerc = NULL,
                              minMeanExpr = NULL, detectThresh = 0) {
  .Deprecated("runFindMarker")
  runFindMarker(inSCE = inSCE, useAssay = useAssay,
                useReducedDim = useReducedDim,
                method = method, cluster = cluster, covariates = covariates,
                log2fcThreshold = log2fcThreshold, fdrThreshold = fdrThreshold,
                minClustExprPerc = minClustExprPerc,
                maxCtrlExprPerc = maxCtrlExprPerc,
                minMeanExpr = minMeanExpr, detectThresh = detectThresh)
}

.calcMarkerExpr <- function(markerTable, inSCE, clusterName, clusterVar,
                            detectThresh) {
  #markerTable <- markerTable[markerTable$Gene %in% rownames(inSCE),]
  genes <- markerTable$Gene
  uniqClust <- unique(markerTable[[clusterName]])
  useAssay <- attr(markerTable, "useAssay" )
  useReducedDim <- attr(markerTable, "useReducedDim")

  cells.in.col <- rep(NA, length(genes))
  cells.out.col <- rep(NA, length(genes))
  exprs.avg <- rep(NA, length(genes))
  for (i in uniqClust) {
    markers <- genes[markerTable[[clusterName]] == i]
    cells.ix <- clusterVar == i
    if (!is.null(useReducedDim)) {
      assay.in <-
        t(expData(inSCE, useReducedDim))[markers, cells.ix, drop = FALSE]
      assay.out <-
        t(expData(inSCE, useReducedDim))[markers, !cells.ix, drop = FALSE]

    } else {
      assay.in <- expData(inSCE, useAssay)[markers, cells.ix, drop = FALSE]
      assay.out <- expData(inSCE, useAssay)[markers, !cells.ix, drop = FALSE]

    }
    assay.in.expr.perc <- rowSums(assay.in > detectThresh) / ncol(assay.in)
    cells.in.col[markerTable[[clusterName]] == i] <- assay.in.expr.perc
    assay.out.expr.perc <- rowSums(assay.out > detectThresh) / ncol(assay.out)
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
#' @rdname getFindMarkerTopTable
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
#' mouseBrainSubsetSCE <- runFindMarker(mouseBrainSubsetSCE,
#'                                      useAssay = "logcounts",
#'                                      cluster = "level1class")
#' getFindMarkerTopTable(mouseBrainSubsetSCE)
getFindMarkerTopTable <- function(inSCE, log2fcThreshold = 1,
                                  fdrThreshold = 0.05, minClustExprPerc = 0.7,
                                  maxCtrlExprPerc = 0.4, minMeanExpr = 1,
                                  topN = 10) {
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop('"inSCE" should be a SingleCellExperiment inherited Object.')
  }
  if(!'findMarker' %in% names(S4Vectors::metadata(inSCE))){
    stop('"findMarker" result not found in metadata. ',
         'Run runFindMarker() in advance')
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

#' @rdname getFindMarkerTopTable
#' @export
findMarkerTopTable <- function(inSCE, log2fcThreshold = 1,
                               fdrThreshold = 0.05, minClustExprPerc = 0.7,
                               maxCtrlExprPerc = 0.4, minMeanExpr = 1,
                               topN = 10) {
  .Deprecated("getFindMarkerTopTable")
  getFindMarkerTopTable(inSCE = inSCE, log2fcThreshold = log2fcThreshold,
                        fdrThreshold = fdrThreshold,
                        minClustExprPerc = minClustExprPerc,
                        maxCtrlExprPerc = maxCtrlExprPerc,
                        minMeanExpr = minMeanExpr,
                        topN = topN)
}
