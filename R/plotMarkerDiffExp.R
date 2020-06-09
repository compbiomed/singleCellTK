#' Plot a heatmap to visualize the result of \code{\link{findMarkerDiffExp}}
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' expression values. Default \code{"logcounts"}.
#' @param log2fcThreshold Only use DEGs with the absolute values of log2FC
#' larger than this value. Default \code{1}
#' @param fdrThreshold Only use DEGs with FDR value smaller than this value.
#' Default \code{0.05}
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
#' @param rowSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{rowDataName} or
#' \code{names(featureAnnotations)}. Default is the value of \code{cluster} in
#' \code{\link{findMarkerDiffExp}} when \code{orderBy} is not \code{NULL}, or
#' \code{NULL} otherwise.
#' @param colSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{colDataName} or
#' \code{names(cellAnnotations)}. Default is the same as \code{rowSplitBy}.
#' @param ... Other arguments passed to \code{\link{plotSCEHeatmap}}.
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object
#' @author Yichen Wang
#' @export
plotMarkerDiffExp <- function(inSCE, useAssay = 'logcounts', orderBy = 'size',
    log2fcThreshold = 1, fdrThreshold = 0.05, decreasing = TRUE,
    rowDataName = NULL, colDataName = NULL, featureAnnotations = NULL,
    cellAnnotations = NULL, featureAnnotationColor = NULL,
    cellAnnotationColor = NULL,
    rowSplitBy = ifelse(is.null(orderBy), NULL,
                        colnames(inSCE@metadata$findMarker)[5]),
    colSplitBy = rowSplitBy, ...){
    if(!inherits(inSCE, 'SingleCellExperiment')){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!'findMarker' %in% names(S4Vectors::metadata(inSCE))){
        stop('"findMarker" result not found in metadata. ',
             'Run findMarkerDiffExp() before plotting.')
    }
    if(!is.null(orderBy)){
        if(length(orderBy) == 1){
            if(!orderBy %in% c('size', 'name')){
                stop('Single charater setting for "orderBy" should be chosen',
                     'from "size" or "name".')
            }
        }# else if(any(!SummarizedExperiment::colData(inSCE)[[cluster]] %in%
         #             orderBy)){
         #   stop('Invalid "orderBy", please input a vector of unique ordered ',
         #        'cluster identifiers that match all clusters in colData(inSCE) ',
         #        'specified by "cluster" to adjust the order of clusters.')
        #}
    }
    # Extract and basic filter
    degFull <- S4Vectors::metadata(inSCE)$findMarker
    if(!all(colnames(degFull)[1:4] ==
           c("Gene", "Pvalue", "Log2_FC", "FDR"))){
        stop('"findMarker" result cannot be interpreted properly')
    }
    if(!is.null(log2fcThreshold)){
        degFull <- degFull[degFull$Log2_FC > log2fcThreshold,]
    }
    if(!is.null(fdrThreshold)){
        degFull <- degFull[degFull$FDR < fdrThreshold,]
    }
    # Remove duplicate by assigning the duplicated genes to the cluster where
    # their log2FC is the highest.
    # Done by keeping all unique genes and the rows  with highest Log2FC entry
    # for each duplicated gene.
    dup.gene <- unique(degFull$Gene[duplicated(degFull$Gene)])
    for(g in dup.gene){
        deg.gix <- degFull$Gene == g
        deg.gtable <- degFull[deg.gix,]
        toKeep <- which.max(deg.gtable$Log2_FC)
        toRemove <- which(deg.gix)[-toKeep]
        degFull <- degFull[-toRemove,]
    }
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
        SummarizedExperiment::rowData(inSCE)[[clusterName]] <- degFull[[clusterName]]
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
    SummarizedExperiment::rowData(inSCE)[[clusterName]] <-
        factor(y, levels = z.order)
    # Organize plotSCEHeatmap arguments
    colDataName <- c(colDataName, clusterName)
    rowDataName <- c(rowDataName, clusterName)
    markerConsistColor <-
        list(marker = dataAnnotationColor(inSCE, 'col')[[clusterName]])
    featureAnnotationColor <- c(featureAnnotationColor, markerConsistColor)

    hm <- plotSCEHeatmap(inSCE, useAssay = useAssay, colDataName = colDataName,
        rowDataName = rowDataName, colSplitBy = colSplitBy, rowSplitBy = rowSplitBy,
        featureAnnotations = featureAnnotations,
        cellAnnotations = cellAnnotations,
        featureAnnotationColor = featureAnnotationColor,
        cellAnnotationColor = cellAnnotationColor,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE, ...)
    return(hm)
}
