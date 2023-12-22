###################################################
###  Setting accessor functions
###################################################
#' @title getTSCANResults accessor function
#' @description SCTK allows user to access all TSCAN related results with
#' \code{"getTSCANResults"}. See details.
#' @param x Input \linkS4class{SingleCellExperiment} object.
#' @param analysisName Algorithm name implemented, should be one of
#' \code{"Pseudotime"}, \code{"DEG"}, or \code{"ClusterDEAnalysis"}.
#' @param pathName Sub folder name within the \code{analysisName}. See details.
#' @param value Value to be stored within the \code{pathName} or
#' \code{analysisName}
#' @details
#' When \code{analysisName = "Pseudotime"}, returns the list result from
#' \code{\link{runTSCAN}}, including the MST structure.
#'
#' When \code{analysisName = "DEG"}, returns the list result from
#' \code{\link{runTSCANDEG}}, including \code{DataFrame}s containing genes that
#' increase/decrease along each the pseudotime paths. \code{pathName} indicates
#' the path index, the available options of which can be listed by
#' \code{listTSCANTerminalNodes}.
#'
#' When \code{analysisName = "ClusterDEAnalysis"}, returns the list result from
#' \code{\link{runTSCANClusterDEAnalysis}}. Here \code{pathName} needs to match
#' with the \code{useCluster} argument when running the algorithm.
#' @export
#' @rdname getTSCANResults
#' @return Get or set TSCAN results
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' results <- getTSCANResults(mouseBrainSubsetSCE, "Pseudotime")
setGeneric("getTSCANResults", signature = "x",
           function(x, analysisName = NULL, pathName = NULL) {
               standardGeneric("getTSCANResults")
           }
)

#' @export
#' @rdname getTSCANResults
#'
setMethod("getTSCANResults", signature(x = "SingleCellExperiment"),
          function(x, analysisName = NULL, pathName = NULL){
              result.names <- listTSCANResults(x)
              if(!analysisName %in% result.names) {
                  stop("The analysis was not found for TSCAN results")
              }
              if (analysisName == "Pseudotime")
                  results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN$Pseudotime
              else {
                  if (is.null(pathName)) {
                      results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]]
                  } else {
                      pathName <- as.character(pathName)
                      results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]][[pathName]]
                  }
              }
              return(results)
          })

#' @export
#' @rdname getTSCANResults
#'
setGeneric("getTSCANResults<-",
           function(x, analysisName, pathName = NULL, value)
               standardGeneric("getTSCANResults<-")
)

#' @export
#' @rdname getTSCANResults
setReplaceMethod("getTSCANResults", signature(x = "SingleCellExperiment"),
                 function(x, analysisName, pathName = NULL, value) {
                     if (analysisName == "Pseudotime")
                         S4Vectors::metadata(x)$sctk$Traj$TSCAN$Pseudotime <- value
                     else {
                         if (is.null(pathName))
                             stop("Have to specify `pathName` for TSCAN DE analyses")
                         pathName <- as.character(pathName)
                         S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]][[pathName]] <- value
                     }
                     return(x)
                 })

#' @export
#' @rdname getTSCANResults
setGeneric("listTSCANResults", signature = "x",
           function(x) {standardGeneric("listTSCANResults")}
)

#' @export
#' @rdname getTSCANResults
setMethod("listTSCANResults", "SingleCellExperiment", function(x){
    all.results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN
    if (is.null(all.results) ||
        !"Pseudotime" %in% names(all.results)) {
        stop("No TSCAN result found. Please run `runTSCAN` first.")
    }
    return(names(all.results))
})

#' @export
#' @rdname getTSCANResults
setGeneric("listTSCANTerminalNodes", signature = "x",
           function(x) {
               standardGeneric("listTSCANTerminalNodes")
           }
)

#' @export
#' @rdname getTSCANResults
setMethod("listTSCANTerminalNodes", signature(x = "SingleCellExperiment"),
          function(x){
              if (is.null(S4Vectors::metadata(x)$sctk$Traj$TSCAN$Pseudotime))
                  stop("No TSCAN result found. Please run `runTSCAN()` first.")
              results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN$Pseudotime$pathClusters
              return(names(results))
          })

###################################################
###  STEP 1:: creating cluster and MST
###################################################

#' @title Run TSCAN to obtain pseudotime values for cells
#' @description Wrapper for obtaining a pseudotime ordering of the cells by
#' projecting them onto the minimum spanning tree (MST)
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useReducedDim Character. A low-dimension representation in
#' \code{reducedDims}, will be used for both clustering if \code{cluster} not
#' specified and MST construction. Default \code{"PCA"}.
#' @param cluster Grouping for each cell in \code{inSCE}. A vector with equal
#' length to the number of the cells in \code{inSCE}, or a single character for
#' retriving \code{colData} variable. Default \code{NULL}, will run
#' \code{runScranSNN} to obtain.
#' @param starter Character. Specifies the starting node from which to compute
#' the pseudotime. Default \code{NULL}, will select an arbitrary node.
#' @param seed An integer. Random seed for clustering if \code{cluster} is not
#' specified. Default \code{12345}.
#' @return The input \code{inSCE} object with pseudotime ordering of the cells
#' along the paths and the cluster label stored in \code{colData}, and other
#' unstructured information in \code{metadata}.
#' @export
#' @author Nida Pervaiz
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
runTSCAN <- function(inSCE,
                     useReducedDim = "PCA",
                     cluster = NULL,
                     starter = NULL,
                     seed = 12345) {
    if (is.null(cluster)) {
        # DON'T RETURN TO `inSCE`
        # It really overwrites existing cluster labels
        sce <- runScranSNN(inSCE, useReducedDim = useReducedDim, seed = seed)
        cluster <- colData(sce)$scranSNN_cluster
        rm(sce)
    } else {
        cluster <- .manageCellVar(inSCE, var = cluster)
    }
    if (length(unique(cluster)) == 1) {
        stop("Only one cluster found. Unable to calculate.")
    }
    message(date(), " ... Running TSCAN to estimate pseudotime")
    inSCE <- scran::computeSumFactors(inSCE, clusters = cluster)
    by.cluster <- scuttle::aggregateAcrossCells(inSCE, ids = cluster)
    centroids <- SingleCellExperiment::reducedDim(by.cluster, useReducedDim)
    mst <- TSCAN::createClusterMST(centroids, clusters = NULL)

    # Map each cell to the closest edge on the MST, reporting also the distance to
    # the corresponding vertices.
    map.tscan <- TSCAN::mapCellsToEdges(inSCE , mst = mst,
                                        use.dimred = useReducedDim,
                                        clusters = cluster)

    # Compute a pseudotime for each cell lying on each path through the MST from a
    # given starting node.
    tscan.pseudo <- TSCAN::orderCells(map.tscan, mst, start = starter)
    common.pseudo <- TrajectoryUtils::averagePseudotime(tscan.pseudo)

    colData(inSCE)$TSCAN_clusters <- cluster
    colData(inSCE)$TSCAN_pseudotime <- common.pseudo

    pathClusters <- list()
    pathList <- data.frame()
    for(i in seq(ncol(tscan.pseudo))){
        pathIndex = as.character(colnames(tscan.pseudo)[i])
        colDataVarName <- paste0("Path_",pathIndex,"_pseudotime")
        inSCE[[colDataVarName]] <- TrajectoryUtils::pathStat(tscan.pseudo)[,pathIndex]
        nonNA.idx <- !is.na(colData(inSCE)[[colDataVarName]])
        pathClusters[[pathIndex]] <- unique(colData(inSCE)$TSCAN_clusters[nonNA.idx])
        pathList[i,1] <- pathIndex
        pathList[i,2] <- toString(pathClusters[[pathIndex]])
        message(date(), " ...   Clusters involved in path index ", pathIndex,
                " are: ", paste(pathClusters[[i]], collapse = ", "))
    }

    maxWidth <- max(stringr::str_length(pathList[, 1]))
    choiceList <- paste0(stringr::str_pad(pathList[, 1], width=maxWidth,
                                          side="right"),
                         "|", pathList[, 2])

    branchClusters <- c()
    edgeList <- mst
    for (i in seq_along(unique(colData(inSCE)$TSCAN_clusters))){
        degree <- sum(edgeList[i] != 0)
        if (degree > 1) branchClusters <- c(branchClusters, i)
    }

    ## Save these results in a list and then make S4 accessor that passes the
    ## entire list
    result <- list(pseudo = tscan.pseudo, mst = mst,
                   maptscan = map.tscan,
                   pathClusters = pathClusters,
                   branchClusters = branchClusters,
                   pathIndexList = choiceList)
    getTSCANResults(inSCE, analysisName = "Pseudotime") <- result

    message(date(), " ...   Number of estimated paths is ", ncol(tscan.pseudo))

    return (inSCE)
}

###################################################
###  plot pseudotime values
###################################################

#' @title Plot MST pseudotime values on cell 2D embedding
#' @description A wrapper function which visualizes outputs from the
#' \code{\link{runTSCAN}} function. Plots the pseudotime ordering of the cells
#' and project them onto the MST.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useReducedDim Saved dimension reduction name in \code{inSCE} object.
#' Required.
#' @return A \code{.ggplot} object with the pseudotime ordering of the cells
#' colored on a cell 2D embedding, and the MST path drawn on it.
#' @export
#' @author Nida Pervaiz
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' plotTSCANResults(inSCE = mouseBrainSubsetSCE,
#'                  useReducedDim = "TSNE_logcounts")
plotTSCANResults <- function(inSCE, useReducedDim = "UMAP") {

    results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
    clusters <- colData(inSCE)$TSCAN_clusters
    by.cluster <- scuttle::aggregateAcrossCells(inSCE, ids = clusters)
    line.data <- TSCAN::reportEdges(by.cluster, mst = results$mst,
                                    clusters = NULL, use.dimred = useReducedDim)

    scater::plotReducedDim(inSCE, dimred = useReducedDim,
                           colour_by = I(colData(inSCE)$TSCAN_pseudotime),
                           text_by = "TSCAN_clusters", text_colour = "black") +
        suppressMessages(ggplot2::scale_colour_viridis_c(name = "Pseudotime")) +
        ggplot2::geom_path(data = line.data,
                           ggplot2::aes_string(x = names(line.data)[2],
                                               y = names(line.data)[3],
                                               group = 'edge'))
}

###################################################
###  STEP 2:: identify expressive genes
###################################################
#' @title Test gene expression changes along a TSCAN trajectory path
#' @description Wrapper for identifying genes with significant changes with
#' respect to one of the TSCAN pseudotime paths
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path index for which the pseudotime values should be used.
#' This corresponds to the terminal node of specific path from the root
#' node to the terminal node. Run \code{listTSCANTerminalNodes(inSCE)} for
#' available options.
#' @param useAssay Character. The name of the assay to use for testing the
#' expression change. Should be log-normalized. Default \code{"logcounts"}
#' @param discardCluster Cluster(s) which are not of use or masks other
#' interesting effects can be discarded. Default \code{NULL}.
#' @return The input \code{inSCE} with results updated in \code{metadata}.
#' @export
#' @author Nida Pervaiz
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' terminalNodes <- listTSCANTerminalNodes(mouseBrainSubsetSCE)
#' mouseBrainSubsetSCE <- runTSCANDEG(inSCE = mouseBrainSubsetSCE,
#'                                    pathIndex = terminalNodes[1])
runTSCANDEG <- function(inSCE,
                        pathIndex,
                        useAssay = "logcounts",
                        discardCluster = NULL) {
    nx <- inSCE
    if (!is.null(discardCluster)) {
        if (any(!discardCluster %in% unique(colData(nx)$TSCAN_clusters))) {
            stop("Not all `discardCluster` exist in TSCAN clusters")
        }
        nx <- nx[, !(colData(nx)$TSCAN_clusters %in% discardCluster)]
    }

    pseudo <- TSCAN::testPseudotime(nx,
                                    pseudotime = nx[[paste0("Path_", pathIndex,
                                                            "_pseudotime")]],
                                    assay.type = useAssay)
    pseudo <- pseudo[order(pseudo$FDR),]
    up.left <- pseudo[pseudo$logFC < 0,]
    up.right <- pseudo[pseudo$logFC > 0,]

    ## Save these results in a list and then make S4 accessor that passes the
    ## entire list
    pathresults <- list(discardClusters = discardCluster, upLeft = up.left,
                        upRight = up.right, useAssay = useAssay)
    getTSCANResults(inSCE, analysisName = "DEG",
                    pathName = as.character(pathIndex)) <- pathresults

    return (inSCE)
}

###################################################
###  plot heatmap of top genes
###################################################
#' @title Plot heatmap of genes with expression change along TSCAN pseudotime
#' @description A wrapper function which visualizes outputs from the
#' \code{\link{runTSCANDEG}} function. Plots the top genes that change in
#' expression with increasing pseudotime along the path in the MST.
#' \code{\link{runTSCANDEG}} has to be run in advance with using the same
#' \code{pathIndex} of interest.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path index for which the pseudotime values should be used.
#' Should have being used in \code{\link{runTSCANDEG}}.
#' @param direction Should we show features with expression increasing or
#' decreeasing along the increase in TSCAN pseudotime? Choices are
#' \code{"both"}, \code{"increasing"} or \code{"decreasing"}.
#' @param topN An integer. Only to plot this number of top genes along the path
#' in the MST, in terms of FDR value. Use \code{NULL} to cancel the top N
#' subscription. Default \code{30}.
#' @param log2fcThreshold Only output DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}.
#' @param useAssay A single character to specify a feature expression matrix in
#' \code{assays} slot. The expression of top features from here will be
#' visualized. Default \code{NULL} use the one used for
#' \code{\link{runTSCANDEG}}.
#' @param featureDisplay Whether to display feature ID and what ID type to
#' display. Users can set default ID type by \code{\link{setSCTKDisplayRow}}.
#' \code{NULL} will display when number of features to display is less than 60.
#' \code{FALSE} for no display. Variable name in \code{rowData} to indicate ID
#' type. \code{"rownames"} or \code{TRUE} for using \code{rownames(inSCE)}.
#' @return A ComplexHeatmap in \code{.ggplot} class
#' @export
#' @author Nida Pervaiz
#' @importFrom S4Vectors metadata
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' terminalNodes <- listTSCANTerminalNodes(mouseBrainSubsetSCE)
#' mouseBrainSubsetSCE <- runTSCANDEG(inSCE = mouseBrainSubsetSCE,
#'                                    pathIndex = terminalNodes[1])
#' plotTSCANPseudotimeHeatmap(mouseBrainSubsetSCE,
#'                            pathIndex = terminalNodes[1])
plotTSCANPseudotimeHeatmap <- function(inSCE,
                                       pathIndex,
                                       direction = c("both", "increasing",
                                                     "decreasing"),
                                       topN = 50,
                                       log2fcThreshold = NULL,
                                       useAssay = NULL,
                                       featureDisplay = metadata(inSCE)$featureDisplay){
    results <- getTSCANResults(inSCE, analysisName = "DEG", pathName = pathIndex)
    direction <- match.arg(direction)
    # Organizing cells
    cell.idx <- rep(TRUE, ncol(inSCE))
    if(!is.null(results$discardClusters)){
        cell.idx <- cell.idx &
            (!colData(inSCE)$TSCAN_clusters %in% results$discardClusters)
    }
    colPathPseudo <- paste0("Path_", pathIndex, "_pseudotime")
    cell.idx <- cell.idx & !is.na(colData(inSCE)[[colPathPseudo]])
    cell.order <- order(colData(inSCE)[[colPathPseudo]])
    cell.order <- cell.order[cell.order %in% which(cell.idx)]

    # Organizing genes
    if (direction == "both") {
        degR <- results$upRight
        degL <- results$upLeft
        if (!is.null(log2fcThreshold)) {
            degR <- degR[abs(degR$logFC) > log2fcThreshold,]
            degL <- degL[abs(degL$logFC) > log2fcThreshold,]
        }
        genesR <- rownames(degR)[seq(min(topN, nrow(degR)))]
        genesL <- rownames(degL)[seq(min(topN, nrow(degL)))]
        genesR <- genesR[!is.na(genesR)]
        genesL <- genesL[!is.na(genesL)]
        genes <- c(genesL, genesR)
        if (length(genes) < 1) stop("No feature passed the filter.")
        direction.df <- data.frame(
            Direction = factor(c(rep("decreasing", length(genesL)),
                                 rep("increasing", length(genesR))),
                               levels = c("increasing", "decreasing")),
            row.names = genes)
    } else {
        if (direction == "increasing") deg <- results$upRight
        if (direction == "decreasing") deg <- results$upLeft
        if (!is.null(log2fcThreshold)) {
            deg <- deg[abs(deg$logFC) > log2fcThreshold,]
        }
        topN <- min(topN, nrow(deg))
        if (topN < 1) stop("No feature passed the filter.")
        genes <- rownames(deg)[seq(topN)]
        direction.df <- data.frame(Direction = rep(direction, length(genes)),
                                   row.names = genes)
    }
    rowLabel <- featureDisplay
    if (is.null(featureDisplay)) {
        if (length(genes) <= 60) rowLabel <- TRUE else rowLabel <- FALSE
    } else if (featureDisplay == "rownames") {
        rowLabel <- TRUE
    }
    cellAnnotationColor = list(
        .viridisPseudoTimeColor(colData(inSCE)[[colPathPseudo]])
    )
    names(cellAnnotationColor) <- colPathPseudo
    if (is.null(useAssay)) useAssay <- results$useAssay
    plotSCEHeatmap(inSCE = inSCE,
                   useAssay = useAssay,
                   cellIndex = cell.order,
                   featureIndex = genes,
                   colDend = FALSE,
                   rowDend = FALSE,
                   cluster_columns = FALSE, cluster_rows = TRUE,
                   colDataName = c("TSCAN_clusters", colPathPseudo),
                   rowLabel = rowLabel,
                   featureAnnotations = direction.df,
                   rowSplitBy = "Direction",
                   cellAnnotationColor = cellAnnotationColor
    )
}

#' Basing on pseudotime is presented within range from 0 to 100
#' Use the viridis color scale to match with `plotTSCANResult` embedding color
#' @param x numeric vector of pseudotime
#' @return function object returned by circlize::colorRamp2
#' @noRd
.viridisPseudoTimeColor <- function(x) {
    a <- max(x, na.rm = TRUE)
    b <- min(x, na.rm = TRUE)
    chunk <- (a - b) / 4
    circlize::colorRamp2(seq(0,4) * chunk + b,
                         c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))
}

###################################################
###  plot expressive genes
###################################################
#' @title Plot expression changes of top features along a TSCAN pseudotime path
#' @description A wrapper function which visualizes outputs from the
#' \code{\link{runTSCANDEG}} function. Plots the genes that increase or decrease
#' in expression with increasing pseudotime along the path in the MST.
#' \code{\link{runTSCANDEG}} has to be run in advance with using the same
#' \code{pathIndex} of interest.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path index for which the pseudotime values should be used.
#' Should have being used in \code{\link{runTSCANDEG}}.
#' @param direction Should we show features with expression increasing or
#' decreeasing along the increase in TSCAN pseudotime? Choices are
#' \code{"increasing"} or \code{"decreasing"}.
#' @param topN An integer. Only to plot this number of top genes that are
#' increasing/decreasing in expression with increasing pseudotime along
#' the path in the MST. Default 10
#' @param useAssay A single character to specify a feature expression matrix in
#' \code{assays} slot. The expression of top features from here will be
#' visualized. Default \code{NULL} use the one used for
#' \code{\link{runTSCANDEG}}.
#' @param featureDisplay Specify the feature ID type to display. Users can set
#' default value with \code{\link{setSCTKDisplayRow}}. \code{NULL} or
#' \code{"rownames"} specifies the rownames of \code{inSCE}. Other character
#' values indicates \code{rowData} variable.
#' @return A \code{.ggplot} object with the facets of the top genes. Expression
#' on y-axis, pseudotime on x-axis.
#' @export
#' @author Nida Pervaiz
#' @importFrom utils head
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment rowData
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' terminalNodes <- listTSCANTerminalNodes(mouseBrainSubsetSCE)
#' mouseBrainSubsetSCE <- runTSCANDEG(inSCE = mouseBrainSubsetSCE,
#'                                    pathIndex = terminalNodes[1])
#' plotTSCANPseudotimeGenes(mouseBrainSubsetSCE,
#'                          pathIndex = terminalNodes[1],
#'                          useAssay = "logcounts")
plotTSCANPseudotimeGenes <- function(inSCE,
                                     pathIndex,
                                     direction = c("increasing", "decreasing"),
                                     topN = 10,
                                     useAssay = NULL,
                                     featureDisplay = metadata(inSCE)$featureDisplay){
    results <- getTSCANResults(inSCE, analysisName = "DEG", pathName = pathIndex)
    direction = match.arg(direction)
    if (!is.null(featureDisplay) &&
        featureDisplay == "rownames")
        featureDisplay <- NULL
    if(!is.null(results$discardClusters)){
        inSCE <- inSCE[,!colData(inSCE)$TSCAN_clusters %in% results$discardClusters]
    }
    if (direction == "decreasing") features <- rownames(results$upLeft)
    if (direction == "increasing") features <- rownames(results$upRight)
    features <- head(features, topN)
    if (!is.null(featureDisplay)) {
        features.idx <- featureIndex(features, inSCE)
        features <- rowData(inSCE)[[featureDisplay]][features.idx]
    }
    if (is.null(useAssay)) useAssay <- results$useAssay
    scater::plotExpression(inSCE,
                           features = features,
                           exprs_values = useAssay,
                           swap_rownames = featureDisplay,
                           x = paste0("Path_",pathIndex,"_pseudotime"),
                           colour_by = "TSCAN_clusters")
}

###################################################
###  STEP3:: identify DE genes in branch cluster
###################################################
#' @title Find DE genes between all TSCAN paths rooted from given cluster
#' @description This function finds all paths that root from a given cluster
#' \code{useCluster}, and performs tests to identify significant features for
#' each path, and are not significant and/or changing in the opposite direction
#' in the other paths. Using a branching cluster (i.e. a node with degree > 2)
#' may highlight features which are responsible for the branching event. MST has
#' to be pre-calculated with \code{\link{runTSCAN}}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useCluster The cluster to be regarded as the root, has to existing in
#' \code{colData(inSCE)$TSCAN_clusters}.
#' @param useAssay Character. The name of the assay to use. This assay should
#' contain log normalized counts. Default \code{"logcounts"}.
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this value.
#' Default \code{0.05}.
#' @return The input \code{inSCE} with results updated in \code{metadata}.
#' @export
#' @author Nida Pervaiz
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' mouseBrainSubsetSCE <- runTSCANClusterDEAnalysis(inSCE = mouseBrainSubsetSCE,
#'                                          useCluster = 1)
runTSCANClusterDEAnalysis <- function(inSCE,
                                      useCluster,
                                      useAssay = "logcounts",
                                      fdrThreshold = 0.05) {
    results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
    message(date(), " ... Finding DEG between TSCAN branches")
    tscan.pseudo <- TSCAN::orderCells(results$maptscan, mst = results$mst,
                                      start = useCluster)
    nPaths <- colnames(tscan.pseudo)

    pathClusters <- list()
    pathList <- data.frame()

    x <- inSCE[,colData(inSCE)$TSCAN_clusters == useCluster]
    store <- list()
    genes <- list()

    for(i in seq(nPaths)){
        # Get Pseudotime of a path
        pathIndex = as.character(nPaths[i])
        branchPseudotime <- TrajectoryUtils::pathStat(tscan.pseudo)[,pathIndex]
        colDataVarName <- paste0("TSCAN_Cluster", useCluster,
                                 "_Path_", pathIndex, "_pseudotime")
        colData(inSCE)[[colDataVarName]] <- branchPseudotime
        nonNA.idx <- !is.na(branchPseudotime)
        pathClusters[[pathIndex]] <- unique(colData(inSCE)$TSCAN_clusters[nonNA.idx])
        pathList[i,1] <- pathIndex
        pathList[i,2] <- toString(pathClusters[[pathIndex]])
        message(date(), " ...   Clusters involved in path index ", pathIndex,
                " are: ", paste(pathClusters[[i]], collapse = ", "))
        # Test along the path
        pseudo.df <- TSCAN::testPseudotime(
            x, df = 1,
            pseudotime = branchPseudotime[colData(inSCE)$TSCAN_clusters == useCluster],
            assay.type = useAssay)
        pseudo.df <- pseudo.df[order(pseudo.df$p.value),]
        store[[pathIndex]] <- pseudo.df
    }
    # Filter against all other paths
    genes <- list()
    for(i in seq(nPaths)){
        pseudo.df <- store[[nPaths[i]]]
        paths <- setdiff(seq(nPaths), i)
        idx <- pseudo.df$FDR < fdrThreshold
        for(j in seq_along(paths)){
            pseudo.otherPath <- store[[paths[j]]]
            idx <- idx & (pseudo.otherPath$p.value >= 0.05 |
                              sign(pseudo.otherPath$logFC) != sign(pseudo.df$logFC))
        }
        pseudo.df <- pseudo.df[which(idx), ]
        genes[[nPaths[i]]] <- pseudo.df[order(pseudo.df$FDR),]
    }

    maxWidth <- max(stringr::str_length(pathList[, 1]))
    choiceList <- paste0(stringr::str_pad(pathList[, 1], width=maxWidth,
                                          side="right"),
                         "|", pathList[, 2])
    expGenes <- list(DEgenes = genes, useAssay = useAssay,
                     terminalNodes = tscan.pseudo,
                     pathIndexList = choiceList)
    getTSCANResults(inSCE, analysisName = "ClusterDEAnalysis",
                    pathName = useCluster) <- expGenes
    return(inSCE)
}

###################################################
###  plot branch cluster
###################################################
#' @title Plot TSCAN pseudotime rooted from given cluster
#' @description This function finds all paths that root from a given cluster
#' \code{useCluster}. For each path, this function plots the recomputed
#' pseudotime starting from the root on a scatter plot which contains cells only
#' in this cluster. MST has to be pre-calculated with \code{\link{runTSCAN}}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useCluster The cluster to be regarded as the root, has to existing in
#' \code{colData(inSCE)$TSCAN_clusters}.
#' @param useReducedDim Saved dimension reduction name in the
#' SingleCellExperiment object. Required.
#' @param combinePlot Must be either \code{"all"} or \code{"none"}. \code{"all"}
#' will combine plots of pseudotime along each path into a single \code{.ggplot}
#' object, while \code{"none"} will output a list of plots. Default
#' \code{"all"}.
#' @return
#' \item{combinePlot = "all"}{A \code{.ggplot} object}
#' \item{combinePlot = "none"}{A list of \code{.ggplot}}
#' @export
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' plotTSCANClusterPseudo(mouseBrainSubsetSCE, useCluster = 1,
#'                        useReducedDim = "TSNE_logcounts")
plotTSCANClusterPseudo <- function(inSCE, useCluster, useReducedDim = "UMAP",
                                   combinePlot = c("all", "none")){
    # Get the plotting info of edges
    results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
    combinePlot <- match.arg(combinePlot)
    clusters <- colData(inSCE)$TSCAN_clusters
    by.cluster <- scuttle::aggregateAcrossCells(inSCE, ids = clusters)
    line.data <- TSCAN::reportEdges(by.cluster, mst = results$mst,
                                    clusters = NULL, use.dimred = useReducedDim)
    line.data.sub <- .getClustersLineData(line.data, useCluster)

    # Get Branch pseudotime
    tscan.pseudo <- TSCAN::orderCells(results$maptscan, mst = results$mst,
                                      start = useCluster)
    nPaths <- colnames(tscan.pseudo)
    x <- inSCE[,colData(inSCE)$TSCAN_clusters == useCluster]
    plotList <- list()
    for (i in seq(nPaths)) {
        branchPseudotime <- TrajectoryUtils::pathStat(tscan.pseudo)[,i]
        branchPseudotime <- branchPseudotime[colData(inSCE)$TSCAN_clusters == useCluster]
        x$pseudotime <- branchPseudotime
        plotList[[i]] <- scater::plotReducedDim(x, dimred = useReducedDim,
                                                colour_by = "pseudotime",
                                                text_by = "TSCAN_clusters",
                                                text_colour = "black") +
            ggplot2::geom_line(data = line.data.sub,
                               ggplot2::aes_string(x = colnames(line.data.sub)[2],
                                                   y = colnames(line.data.sub)[3],
                                                   group = "edge"))
    }
    if (combinePlot == "all") {
        plotList <- cowplot::plot_grid(plotlist = plotList)
    }
    return(plotList)
}

###################################################
###  plot gene of interest in branch cluster
###################################################

#' @title Plot features identified by \code{\link{runTSCANClusterDEAnalysis}} on
#' cell 2D embedding with MST overlaid
#' @description A wrapper function which plot the top features expression
#' identified by \code{\link{runTSCANClusterDEAnalysis}} on the 2D embedding of
#' the cells cluster used in the analysis. The related MST edges are overlaid.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useCluster Choose a cluster used for identifying DEG with
#' \code{\link{runTSCANClusterDEAnalysis}}. Required.
#' @param pathIndex Specifies one of the branching paths from \code{useCluster}
#' and plot the top DEGs on this path. Ususally presented by the terminal
#' cluster of a path. By default \code{NULL} plot top DEGs of all paths.
#' @param useReducedDim A single character for the matrix of 2D embedding.
#' Should exist in \code{reducedDims} slot. Default \code{"UMAP"}.
#' @param topN Integer. Use top N genes identified. Default \code{9}.
#' @param useAssay A single character for the feature expression matrix. Should
#' exist in \code{assayNames(inSCE)}. Default \code{NULL} for using the one used
#' in \code{\link{runTSCANClusterDEAnalysis}}.
#' @param featureDisplay Specify the feature ID type to display. Users can set
#' default value with \code{\link{setSCTKDisplayRow}}. \code{NULL} or
#' \code{"rownames"} specifies the rownames of \code{inSCE}. Other character
#' values indicates \code{rowData} variable.
#' @param combinePlot Must be either \code{"all"} or \code{"none"}. \code{"all"}
#' will combine plots of each feature into a single \code{.ggplot} object,
#' while \code{"none"} will output a list of plots. Default \code{"all"}.
#' @return A \code{.ggplot} object of cell scatter plot, colored by the
#' expression of a gene identified by \code{\link{runTSCANClusterDEAnalysis}},
#' with the layer of trajectory.
#' @export
#' @author Yichen Wang
#' @importFrom S4Vectors metadata
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' mouseBrainSubsetSCE <- runTSCANClusterDEAnalysis(inSCE = mouseBrainSubsetSCE,
#'                                                  useCluster = 1)
#' plotTSCANClusterDEG(mouseBrainSubsetSCE, useCluster = 1,
#'                     useReducedDim = "TSNE_logcounts")
plotTSCANClusterDEG <- function(
        inSCE,
        useCluster,
        pathIndex = NULL,
        useReducedDim = "UMAP",
        topN = 9,
        useAssay = NULL,
        featureDisplay = metadata(inSCE)$featureDisplay,
        combinePlot = c("all", "none")
) {
    MSTResults <- getTSCANResults(inSCE, analysisName = "Pseudotime")
    if (length(useCluster) != 1) stop("Can only specify one cluster")
    DEResults <- getTSCANResults(inSCE, analysisName = "ClusterDEAnalysis",
                                 pathName = useCluster)
    if (is.null(DEResults)) {
        stop("Branch DE result of specified cluster not found. Run ",
             "`runTSCANClusterDEAnalysis()` with the same `useCluster` first.")
    }
    combinePlot <- match.arg(combinePlot)
    if (is.null(pathIndex)) {
        deg <- do.call(rbind, DEResults$DEgenes)
        deg <- deg[order(deg$FDR),]
    } else {
        pathIndex <- as.character(pathIndex)
        deg <- DEResults$DEgenes[[pathIndex]]
    }
    features <- rownames(deg)[seq(min(topN, nrow(deg), na.rm = TRUE))]
    if (is.null(useAssay)) useAssay <- DEResults$useAssay
    plotTSCANDimReduceFeatures(inSCE = inSCE, features = features,
                               useReducedDim = useReducedDim,
                               useAssay = useAssay, useCluster = useCluster,
                               featureDisplay = featureDisplay,
                               combinePlot = combinePlot)
}

#' @title Plot feature expression on cell 2D embedding with MST overlaid
#' @description A wrapper function which plots all cells or cells in chosen
#' cluster. Each point is a cell colored by the expression of a feature of
#' interest, the relevant edges of the MST are overlaid on top.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param features Choose the feature of interest to explore the expression
#' level on the trajectory. Required.
#' @param useReducedDim A single character for the matrix of 2D embedding.
#' Should exist in \code{reducedDims} slot. Default \code{"UMAP"}.
#' @param useAssay A single character for the feature expression matrix. Should
#' exist in \code{assayNames(inSCE)}. Default \code{"logcounts"}.
#' @param by Where should \code{features} be found? \code{NULL},
#' \code{"rownames"} for \code{rownames(inSCE)}, otherwise will be regarded as
#' \code{rowData} variable.
#' @param useCluster Choose specific clusters where gene expression needs to be
#' visualized. By default \code{NULL}, all clusters are chosen.
#' @param featureDisplay Specify the feature ID type to display. Users can set
#' default value with \code{\link{setSCTKDisplayRow}}. \code{NULL} or
#' \code{"rownames"} specifies the rownames of \code{inSCE}. Other character
#' values indicates \code{rowData} variable.
#' @param combinePlot Must be either \code{"all"} or \code{"none"}. \code{"all"}
#' will combine plots of each feature into a single \code{.ggplot} object,
#' while \code{"none"} will output a list of plots. Default \code{"all"}.
#' @return A \code{.ggplot} object of cell scatter plot, colored by the
#' expression of a gene of interest, with the layer of trajectory.
#' @export
#' @author Yichen Wang
#' @importFrom S4Vectors metadata
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- runTSCAN(inSCE = mouseBrainSubsetSCE,
#'                                 useReducedDim = "PCA_logcounts")
#' plotTSCANDimReduceFeatures(inSCE = mouseBrainSubsetSCE,
#'                            features = "Tshz1",
#'                            useReducedDim = "TSNE_logcounts")
plotTSCANDimReduceFeatures <- function(
        inSCE,
        features,
        useReducedDim = "UMAP",
        useAssay = "logcounts",
        by = "rownames",
        useCluster = NULL,
        featureDisplay = metadata(inSCE)$featureDisplay,
        combinePlot = c("all", "none"))
{
    results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
    combinePlot <- match.arg(combinePlot)
    clusters <- colData(inSCE)$TSCAN_clusters
    by.cluster <- scuttle::aggregateAcrossCells(inSCE, ids = clusters)
    line.data <- TSCAN::reportEdges(by.cluster, mst = results$mst,
                                    clusters = NULL, use.dimred = useReducedDim)
    if (!is.null(useCluster)) {
        line.data <- .getClustersLineData(line.data, useCluster)
        inSCE <- inSCE[, clusters %in% useCluster]
    }
    if (by == "rownames") by <- NULL
    plotList <- list()
    features <- stats::na.omit(features)
    for (f in features) {
        # Should not enter the loop if features is length zero after NA omit
        g <- plotSCEDimReduceFeatures(inSCE, feature = f,
                                      useAssay = useAssay,
                                      featureLocation = by,dim1 = 1, dim2 = 2,
                                      featureDisplay = featureDisplay,
                                      reducedDimName = useReducedDim,
                                      title = f) +
            ggplot2::geom_line(data = line.data,
                               ggplot2::aes_string(x = colnames(line.data)[2],
                                                   y = colnames(line.data)[3],
                                                   group = "edge"),
                               inherit.aes = FALSE)
        # Text labeling code credit: scater
        reddim <- as.data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                                 useReducedDim))
        text_out <- scater::retrieveCellInfo(inSCE, "TSCAN_clusters",
                                             search = "colData")
        by_text_x <- vapply(split(reddim[,1], text_out$val),
                            stats::median, FUN.VALUE = 0)
        by_text_y <- vapply(split(reddim[,2], text_out$val),
                            stats::median, FUN.VALUE = 0)
        g <- g + ggrepel::geom_text_repel(
            data = data.frame(x = by_text_x,
                              y = by_text_y,
                              label = names(by_text_x)),
            mapping = ggplot2::aes_string(x = "x",
                                          y = "y",
                                          label = "label"),
            inherit.aes = FALSE, na.rm = TRUE,
        )
        plotList[[f]] <- g
    }
    if (combinePlot == "all") {
        if (length(plotList) > 0)
            plotList <- cowplot::plot_grid(plotlist = plotList)
    }
    return(plotList)
}

#' Extract edges that connects to the clusters of interest
#' @param line.data data.frame of three columns. First column should be named
#' "edge", the other two are coordinates of one of the vertices that this
#' edge connects.
#' @param useCluster Vector of clusters, that are going to be the vertices for
#' edge finding.
#' @return data.frame, subset rows of line.data
#' @noRd
.getClustersLineData <- function(line.data, useCluster) {
    idx <- vapply(useCluster, function(x){
        grepl(paste0("^", x, "--"), line.data$edge) |
            grepl(paste0("--", x, "$"), line.data$edge)
    }, logical(nrow(line.data)))
    # `idx` returned was c("matrix", "array") class. Each column is a logical
    # vector for edge selection of each cluster.
    # Next, do "or" for each row.
    idx <- rowSums(idx) > 0
    line.data[idx,]
}
