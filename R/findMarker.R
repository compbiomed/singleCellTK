#' Find the marker gene set for each cluster
#' With an input SingleCellExperiment object and specifying the clustering
#' labels, this function iteratively call the differential expression analysis
#' on each cluster against all the others.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param cluster One single character to specify a column in
#' \code{colData(inSCE)} for the clustering label. Alternatively, a vector or
#' a factor is also acceptable. Default \code{"cluster"}.
#' @param useThresh Whether to use adaptive thresholding to filter genes.
#' Default \code{FALSE}.
#' @param freqExpressed A numeric threshold that the genes expressed in less
#' than this fraction of cells will be removed. Default \code{0.1}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this
#' value. Default \code{1}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$findMarker} updated with a data.table of the up-
#' regulated DEGs for each cluster.
#' @export
#' @author Yichen Wang
findMarkerDiffExp <- function(inSCE, useAssay = 'logcounts',
                              cluster = 'cluster', log2fcThreshold = NULL,
                              fdrThreshold = 1, useThresh = FALSE,
                              freqExpressed = 0.1){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop('"useAssay" name: ', useAssay, ' not found.')
    }
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
        inSCE <- runMAST(inSCE, useAssay = useAssay, index1 = clusterIndex,
                         analysisName = paste0('findMarker', c),
                         onlyPos = TRUE, log2fcThreshold = log2fcThreshold,
                         fdrThreshold = fdrThreshold, useThresh = useThresh,
                         groupName1 = c, groupName2 = 'others')
                         #freqExpressed = freqExpressed)
    }
    degFull <- NULL
    for(c in uniqClust){
        degTable <-
            S4Vectors::metadata(inSCE)$diffExp[[paste0('findMarker', c)]]$result
        degTable[[clusterName]] <- c
        if(is.null(degFull)){
            degFull <- degTable
        } else {
            degFull <- rbind(degFull, degTable)
        }
    }
    S4Vectors::metadata(inSCE)$diffExp[paste0('findMarker', uniqClust)] <- NULL
    if(length(names(S4Vectors::metadata(inSCE)$diffExp)) == 0){
        S4Vectors::metadata(inSCE)$diffExp <- NULL
    }
    degFull <- degFull[stats::complete.cases(degFull),]
    S4Vectors::metadata(inSCE)$findMarker <- degFull
    return(inSCE)
}
