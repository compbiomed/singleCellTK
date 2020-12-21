#' Find the marker gene set for each cluster
#' With an input SingleCellExperiment object and specifying the clustering
#' labels, this function iteratively call the differential expression analysis
#' on each cluster against all the others.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param method A single character for specific differential expression
#' analysis method. Choose from \code{'MAST'}, \code{'DESeq2'}, \code{'Limma'}.
#' Default \code{"MAST"}.
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
#' @author Yichen Wang
findMarkerDiffExp <- function(inSCE, useAssay = 'logcounts',
                              method = c('MAST', "DESeq2", "Limma"),
                              cluster = 'cluster', covariates = NULL,
                              log2fcThreshold = 0.25, fdrThreshold = 0.05,
                              minClustExprPerc = 0.6, maxCtrlExprPerc = 0.4,
                              minMeanExpr = 0.5){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
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
    degFull <- .calcMarkerExpr(degFull, inSCE, clusterName,
                               c(minClustExprPerc, maxCtrlExprPerc,
                                 minMeanExpr))
    S4Vectors::metadata(inSCE)$findMarker <- degFull
    return(inSCE)
}

.calcMarkerExpr <- function(markerTable, inSCE, clusterName, params) {
    genes <- markerTable$Gene
    cellLabels <- SummarizedExperiment::colData(inSCE)[[clusterName]]
    useAssay <- attr(markerTable, "useAssay")
    cells.in.col <- c()
    cells.out.col <- c()
    exprs.avg <- c()
    toRemove <- integer()
    for (i in seq_len(nrow(markerTable))) {
        gene <- genes[i]
        cluster <- markerTable[[clusterName]][i]
        if (gene %in% rownames(inSCE)) {
            # Logical indices
            cells.in <- cellLabels == cluster
            cells.out <- cellLabels != cluster
            # cells.in expressed percentage
            cells.in.expr <- SummarizedExperiment::assay(inSCE, useAssay)[gene, cells.in]
            cells.in.perc <- sum(cells.in.expr > 0) / length(which(cells.in))
            # cells.out expressed percentage
            cells.out.expr <- SummarizedExperiment::assay(inSCE, useAssay)[gene, cells.out]
            cells.out.perc <- sum(cells.out.expr > 0) / length(which(cells.out))
            # Average Expression in the cluster
            cells.in.mean <- mean(cells.in.expr)
            if (cells.in.perc >= params[1] &&
                cells.out.perc <= params[2] &&
                cells.in.mean >= params[3]) {
                cells.in.col <- c(cells.in.col, cells.in.perc)
                cells.out.col <- c(cells.out.col, cells.out.perc)
                exprs.avg <- c(exprs.avg, cells.in.mean)
            } else {
                toRemove <- c(toRemove, i)
            }
        }
    }
    markerTable <- markerTable[-toRemove,]
    markerTable$clusterExprPerc <- cells.in.col
    markerTable$ControlExprPerc <- cells.out.col
    markerTable$clusterAveExpr <- exprs.avg
    return(markerTable)
}
