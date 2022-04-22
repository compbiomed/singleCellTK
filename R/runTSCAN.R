###################################################
###  Setting accessor functions
###################################################
#' @title getTSCANResults accessor function 
#' @param x Input \linkS4class{SingleCellExperiment} object.
#' @param analysisName Algorithm name implemented
#' @param pathName Sub folder name within the \code{analysisName}
#' @param value Value to be stored within the \code{pathName} or 
#' \code{analysisName}
#' @export
#' @rdname getTSCANResults
#' @return Get or set TSCAN results
setGeneric("getTSCANResults", signature = "x",
           function(x, analysisName=NULL, pathName = NULL) {
             standardGeneric("getTSCANResults")
           }
)

#' @export
#' @rdname getTSCANResults
#' 
setMethod("getTSCANResults", signature(x = "SingleCellExperiment"), 
          function(x, analysisName=NULL, pathName = NULL){
  result.names <- listTSCANResults(x)
  if(!analysisName %in% result.names) {
    stop("The analysis was not found in the results for tool 'Trajectory'")  
  }
  if (is.null(pathName)){
    results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]]
    
  }
  else{
    results <- S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]][[pathName]]
    
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
  
  if (is.null(pathName)){
    S4Vectors::metadata(x)$sctk$Traj$TSCAN[[analysisName]] <- value
    
  }
  else{
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
  if(is.null(all.results) || length(names(all.results)) == 0) {
    stop("No results from 'Trajectory' are found. Please run `runTSCAN` first.") 
  }
  return(names(all.results))
})


###################################################
###  STEP 1:: creating cluster and MST
###################################################

#' @title Run runTSCAN function to obtain pseudotime values for cells
#' @description Wrapper for obtaining a pseudotime ordering of the cells by 
#' projecting them onto the MST
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param cluster Grouping for each cell in \code{inSCE}. A user may input a 
#' vector equal length to the number of the samples in \code{inSCE}, or can be 
#' retrieved from the \code{colData} slot. Default \code{NULL}.
#' @param seed An integer. Set the seed for random process that happens only in 
#' "random" generation. Default \code{12345}.
#' @param useReducedDim Character. Saved dimension reduction name in 
#' \code{inSCE} object. Required. Used for specifying which low-dimension 
#' representation to perform the clustering algorithm and building nearest 
#' neighbor graph on. Default \code{"PCA"}
#'
#' @return A \linkS4class{SingleCellExperiment} object with pseudotime ordering 
#' of the cells along the paths
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
runTSCAN <- function(inSCE, 
                     useReducedDim, 
                     cluster = NULL, 
                     seed = 12345) { 
  
  if (is.null(seed)) {
    if(is.null(cluster)){
      inSCE <- runScranSNN(inSCE, useReducedDim = useReducedDim)
      cluster <- colData(inSCE)$"scranSNN_cluster"
    }
  }
  else{
    withr::with_seed(   
      seed, 
      if(is.null(cluster)){
        inSCE <- runScranSNN(inSCE, useReducedDim = useReducedDim)
        cluster <- colData(inSCE)$"scranSNN_cluster"
      })
  }
  inSCE <- scran::computeSumFactors(inSCE, clusters = cluster)
  by.cluster <- scuttle::aggregateAcrossCells(inSCE , ids = cluster)
  centroids <- SingleCellExperiment::reducedDim(by.cluster, useReducedDim)         
  mst <- TSCAN::createClusterMST(centroids, clusters = NULL)      
  
  # Map each cell to the closest edge on the MST, reporting also the distance to 
  # the corresponding vertices.
  map.tscan <- TSCAN::mapCellsToEdges(inSCE , mst = mst, 
                                      use.dimred = useReducedDim , 
                                      clusters = cluster)
  
  # Compute a pseudotime for each cell lying on each path through the MST from a 
  # given starting node.
  tscan.pseudo <- TSCAN::orderCells(map.tscan, mst)
  common.pseudo <- TrajectoryUtils::averagePseudotime(tscan.pseudo) 
  
  colData(inSCE)$"TSCAN_clusters" <- cluster
  colData(inSCE)$"TSCAN_pseudotime" <- common.pseudo 
  
  pathClusters <- list()
  pathList <- data.frame()
  for(i in seq(ncol(tscan.pseudo))){
    pathIndex = as.character(colnames(tscan.pseudo)[i])
    inSCE[[paste0("Path_",pathIndex,"_pseudotime") ]] <- TrajectoryUtils::pathStat(tscan.pseudo)[,pathIndex]
    nonNA <- inSCE[,!is.na(inSCE[[paste0("Path_",pathIndex,"_pseudotime") ]])]
    pathClusters[[pathIndex]] <- (sort(unique(colData(nonNA)$TSCAN_clusters)))
    
    pathList[i,1] <- pathIndex
    pathList[i,2] <- toString (pathClusters[[pathIndex]])
    message("Cluster involved in path ",pathIndex," are: ", pathClusters[i])
  }
  
  maxWidth <- max(stringr::str_length(pathList[, 1]))
  choiceList <- pathList[, 1]
  names(choiceList) <- paste0(stringr::str_pad(pathList[, 1], width=maxWidth, 
                                               side="right"), 
                              "|", pathList[, 2])
  
  branchClustersList <- c()
  edgeList <- mst
  for (i in seq_along(unique(colData(inSCE)$TSCAN_clusters))){
    clustEdges <- sum(edgeList[i] != 0)
    if(clustEdges > 1){
      branchClustersList <- sort(BiocGenerics::append(i,branchClustersList))
    }
    
  }
  
  ## Save these results in a list and then make S4 accessor that passes the 
  ## entire list
  result <- list(pseudo = tscan.pseudo, mst = mst, clusters = cluster, 
                 by.cluster = by.cluster, maptscan = map.tscan, 
                 pathClusters = pathClusters, 
                 branchClustersList = branchClustersList, 
                 pathIndexList = names(choiceList))
  getTSCANResults(inSCE, analysisName = "Pseudotime") <- result 
  
  message("Number of estimated paths is ", ncol(tscan.pseudo) )
  
  return (inSCE)
}  

###################################################
###  plot pseudotime values
###################################################

#' @title Plot MST pseudotime values for cells
#' @description A wrapper function which visualizes outputs from the 
#' \code{\link{runTSCAN}} function. Plots the pseudotime ordering of the cells 
#' by projecting them onto the MST
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useReducedDim Saved dimension reduction name in \code{inSCE} object. 
#' Required.
#'
#' @return A plot with the pseudotime ordering of the cells by projecting them 
#' onto the MST. 
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' plotTSCANResults(inSCE = sce, useReducedDim = "TSNE")

plotTSCANResults <- function(inSCE, 
                             useReducedDim){
  
  results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
  line.data <- TSCAN::reportEdges(results$by.cluster, mst = results$mst, 
                                  clusters = NULL, use.dimred = useReducedDim) 
  
  a <- b <- c <- NULL
  a = unlist(line.data[,2])
  b = unlist(line.data[,3])
  c = unlist(line.data[,1])
  
  scater::plotReducedDim(inSCE, dimred = useReducedDim, 
                         colour_by = I(colData(inSCE)$"TSCAN_pseudotime"), 
                         text_by = "TSCAN_clusters", text_colour = "black") + 
    ggplot2::geom_path(data = line.data, ggplot2::aes(x = a, y = b, group = c)) 
}


###################################################
###  STEP 2:: identify expressive genes
###################################################
#' @title Run runTSCANDEG function to obtain changes along a trajectory
#' @description Wrapper for identifying genes with significant changes with 
#' respect to one 
#' of the TSCAN pseudotimes
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path number for which the pseudotime values should be used. 
#' \code{PathIndex} corresponds to one path from the root node to one of the 
#' terminal nodes.
#' @param useAssay Character. The name of the assay to use. This assay should 
#' contain log normalized counts.
#' @param discardCluster Optional. Clusters which are not of use or masks other 
#' interesting effects can be discarded.
#' @param log2fcThreshold Only output DEGs with the absolute values of log2FC 
#' larger than this value. Default \code{0}
#'
#' @return A \linkS4class{SingleCellExperiment} object with genes that decrease 
#' and increase in expression with increasing pseudotime along the path in the 
#' MST.
#' @export 
#' @author Nida Pervaiz
#'
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)

runTSCANDEG <- function(inSCE, 
                        pathIndex, 
                        useAssay = "logcounts", 
                        discardCluster = NULL, 
                        log2fcThreshold = 0 ) {  
  
  
  if(!is.null(discardCluster)){
    for (i in seq_along(discardCluster)){
      if (i == 1){
        nx <- inSCE[, colData(inSCE)$"TSCAN_clusters" != discardCluster[i]]
      }
      else{
        nx <- nx[, colData(nx)$"TSCAN_clusters" != discardCluster[i]]
        
      }
    }
    
  }
  else{
    nx <- inSCE
  }
  
  pseudo <- TSCAN::testPseudotime(nx, 
                                  pseudotime = nx[[paste0("Path_", pathIndex, 
                                                          "_pseudotime")]], 
                                  assay.type = useAssay)
  pseudo$SYMBOL <- rowData(nx)$Symbol
  sorted <- pseudo[order(pseudo$p.value),]
  up.left <- sorted[sorted$logFC < log2fcThreshold,]
  up.right <- sorted[sorted$logFC > log2fcThreshold,]
  
  ## Save these results in a list and then make S4 accessor that passes the 
  ## entire list
  pathresults <- list(discardClusters = discardCluster, upLeft = up.left, 
                      upRight = up.right) 
  getTSCANResults(inSCE, analysisName = "DEG", 
                  pathName = as.character(pathIndex)) <- pathresults
  
  return (inSCE)
}

###################################################
###  plot heatmap of top genes
###################################################
#' @title Run plotTSCANPseudotimeHeatmap function to plot heatmap for top genes
#' @description A wrapper function which visualizes outputs from the 
#' \code{\link{runTSCANDEG}} function. 
#' Plots the top genes that increase in expression with increasing pseudotime 
#' along the path in the MST 
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path number for which the pseudotime values should be used. 
#' PathIndex corresponds to one path from the root node to one of the terminal 
#' nodes.
#' @param topN An integer. Only to plot this number of top genes along the path 
#' in the MST, in terms of log2FC value. Use \code{NULL} to cancel the top N 
#' subscription. Default \code{50}.
#' @return A plot with the top genes that increase in expression with increasing
#' pseudotime along the path in the MST.
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
#' plotTSCANPseudotimeHeatmap(inSCE = sce, pathIndex = 4,topN = 5)

plotTSCANPseudotimeHeatmap <- function(inSCE, 
                                       pathIndex, 
                                       topN = 50){
  
  pathIndex = as.character(pathIndex)
  
  results <- getTSCANResults(inSCE, analysisName = "DEG", pathName = pathIndex)  
  if(!is.null(results$discardClusters)){
    for (i in seq_along(results$discardClusters)){
      if (i == 1){
        x <- inSCE[, colData(inSCE)$"TSCAN_clusters" != results$discardClusters[i]]
      }
      else{
        x <- x[, colData(x)$"TSCAN_clusters" != results$discardClusters[i]]
        
      }
    }
    
  }
  else{
    x <- inSCE
  }
  
  colPathPseudo <- paste0("Path_",pathIndex,"_pseudotime")
  on.path <- !is.na(colData(x)[names(colData(x)) == colPathPseudo])
  on.path <- as.vector(on.path[,max(ncol(on.path))])
  
  x <- x[,on.path]
  plotSCEHeatmap(inSCE = x,
                 featureIndex = utils::head(results$upRight$SYMBOL, topN),
                 colDend = FALSE, 
                 rowDend = FALSE,
                 cluster_columns = FALSE, 
                 colDataName = c("TSCAN_clusters", colPathPseudo), 
                 rowLabel = TRUE, 
                 colSplitBy = colPathPseudo
  )
}

###################################################
###  plot expressive genes
###################################################
#' @title Run plotTSCANPseudotimeGenes function to plot genes with significant 
#' changes
#' @description A wrapper function which visualizes outputs from the 
#' \code{\link{runTSCANDEG}} function. Plots the genes that increase or decrease
#' in expression with increasing pseudotime along the path in the MST.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param pathIndex Path number for which the pseudotime values should be used. 
#' PathIndex corresponds to one path from the root node to one of the terminal 
#' nodes.
#' @param direction Which direction to use. Choices are increasing or 
#' decreasing.
#' @param n An integer. Only to plot this number of top genes that are 
#' increasing/decreasing in expression with increasing pseudotime along the path
#' in the MST. Default \code{10}.
#'
#' @return A plot with the top genes that increase/decrease in expression with 
#' increasing pseudotime along the path in the MST 
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
#' plotTSCANPseudotimeGenes(inSCE = sce, pathIndex = 4, 
#'                          direction = "increasing")

plotTSCANPseudotimeGenes <- function (inSCE, 
                                      pathIndex, 
                                      direction = c("increasing", "decreasing"), 
                                      n = 10){
  
  pathIndex = as.character(pathIndex)
  results <- getTSCANResults(inSCE, analysisName = "DEG", pathName = pathIndex) 
  if(!is.null(results$discardClusters)){
    for (i in seq_along(results$discardClusters)){
      if (i == 1){
        y <- inSCE[, colData(inSCE)$"TSCAN_clusters" != results$discardClusters[i]]
      }
      else{
        y <- y[, colData(y)$"TSCAN_clusters" != results$discardClusters[i]]
        
      }
    }
    
  }
  else{
    y <- inSCE
  }
  direction = match.arg(direction)
  if (direction == "decreasing"){
    scater::plotExpression(y, 
                           features = utils::head(results$upLeft$SYMBOL, n), 
                           swap_rownames = "Symbol",
                           x = paste0("Path_",pathIndex,"_pseudotime"), 
                           colour_by = "TSCAN_clusters")  
  }
  else{ 
    scater::plotExpression(y, 
                           features = utils::head(results$upRight$SYMBOL, n), 
                           swap_rownames = "Symbol",
                           x = paste0("Path_",pathIndex,"_pseudotime"), 
                           colour_by = "TSCAN_clusters")  
    
  }
}


###################################################
###  STEP3:: identify DE genes in branch cluster 
###################################################
#' @title Run runTSCANClusterDEAnalysis function to observe changes between 
#' paths and to obtain DE genes
#' @description Wrapper for looking for differences in expression between paths 
#' of a branched trajectory. The differential expression analysis may highlight 
#' genes which are responsible for the branching event
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useClusters Choose the cluster containing the branch point in the data
#' in order to recompute the pseudotimes so that the root lies at the cluster 
#' center, allowing us to detect genes that are associated with the divergence 
#' of the branches.
#' @param useAssay Character. The name of the assay to use. This assay should 
#' contain log normalized counts.
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this value.
#' Default \code{0.05}.
#' @return A \linkS4class{SingleCellExperiment} object with DE genes that are 
#' significant in our path of interest and are not significant and/or changing 
#' in the opposite direction in the other paths. 
#' @export 
#' @author Nida Pervaiz
#'
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
#' sce <- runTSCANClusterDEAnalysis(inSCE = sce, useClusters = 5)

runTSCANClusterDEAnalysis <- function(inSCE, 
                                      useClusters , 
                                      useAssay = "logcounts", 
                                      fdrThreshold = 0.05){ 
  
  results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
  starter <- colData(inSCE)$"TSCAN_clusters" == useClusters
  tscan.pseudo <- TSCAN::orderCells(results$maptscan, mst = results$mst, 
                                    start = useClusters)
  nPaths <- colnames(tscan.pseudo)
  
  pathClusters <- list()
  for(i in seq(nPaths)){ 
    pathIndex = as.character(nPaths[i])
    inSCE$"branchPseudotime" <- TrajectoryUtils::pathStat(tscan.pseudo)[,nPaths[i]]
    nonNA <- inSCE[,!is.na(inSCE$"branchPseudotime")]
    pathClusters[[pathIndex]] <- sort(unique(colData(nonNA)$TSCAN_clusters))
    names(colData(inSCE))[names(colData(inSCE)) == 'branchPseudotime'] <- 
      paste0("TSCAN_Cluster",useClusters,"_Path_" ,nPaths[i],"_pseudotime") 
    
    message("Clusters involved in path ",pathIndex," are: ", pathClusters[i])
  }
  
  x <- inSCE[,starter]
  pseudo1 <- {}
  store <- {}
  
  for (i in seq(nPaths)){
    pseudo <- colData(x)[names(colData(x)) == paste0("TSCAN_Cluster",
                                                     useClusters,
                                                     "_Path_" ,
                                                     nPaths[i],
                                                     "_pseudotime")]
    pseudo <- TSCAN::testPseudotime(x, df = 1, pseudotime = pseudo, 
                                    assay.type = useAssay)
    pseudo[[1]]$SYMBOL <- rowData(x)$SYMBOL
    pseudo1 <- pseudo[[1]][order(pseudo[[1]]$p.value),]
    store[[nPaths[i]]] <- list(allGenes = pseudo1)
  }
  
  genes <- {}
  thresh <- {}
  
  for(i in seq(nPaths)){
    
    thresh <- data.frame(store[[nPaths[i]]]$allGenes[which(store[[nPaths[i]]]$allGenes$FDR <= fdrThreshold) ,])
    paths <- setdiff(seq(nPaths), i)
    
    for(j in seq_along(paths)){ 
      
      pvals <- store[[paths[j]]]$allGenes$p.value
      fc <- store[[paths[j]]]$allGenes$logFC
      nonde.ix <- store[[nPaths[i]]]$allGenes[which(pvals >= fdrThreshold | sign(fc) != sign(store[[nPaths[i]]]$allGenes$logFC)),]
      
      if(!is.null(nonde.ix)){
        
        pathgenes <- data.frame(nonde.ix)
        thresh <- generics::intersect(thresh, pathgenes)
      }
    }
    
    genes[[nPaths[i]]] <- thresh[order(thresh$p.value),]
    
  }
  
  expGenes <- list(allgenes = store, DEgenes = genes, 
                   terminalNodes = tscan.pseudo)
  getTSCANResults(inSCE, analysisName = "ClusterDEAnalysis", 
                  pathName =  as.character(useClusters)) <- expGenes
  
  message("Number of estimated paths of cluster ", useClusters, " is ",
          ncol(tscan.pseudo), 
          ". Following are the terminal nodes for each path respectively: ",
          data.frame(colnames(tscan.pseudo)))
  return(inSCE)
}

###################################################
###  plot branch cluster
###################################################
#' @title Run plotClusterPseudo function to plot TSCAN-derived pseudotimes 
#' around cluster in the dataset.
#' @description A wrapper function which visualizes outputs from the 
#' \code{\link{runTSCANClusterDEAnalysis}} function. Each point is a cell in the
#' cluster and is colored by its pseudotime value along the path to which it was
#' assigned. 
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useClusters Choose the cluster containing the branch point in the data
#' in order to recompute the pseudotimes so that the root lies at the cluster 
#' center, allowing us to detect genes that are associated with the divergence 
#' of the branches.
#' @param pathIndex Path number for which the pseudotime values should be used. 
#' PathIndex corresponds to one path from the root node to one of the terminal 
#' nodes.
#' @param useReducedDim Saved dimension reduction name in \code{inSCE}. 
#' Required.
#'
#' @return A plots with the TSCAN-derived pseudotimes of all the cells along the
#' path belonging to the cluster
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
#' sce <- runTSCANClusterDEAnalysis(inSCE = sce, useClusters = 5)
#' plotClusterPseudo(inSCE = sce, useClusters = 5, pathIndex = NULL, 
#'                   useReducedDim = "TSNE")

plotClusterPseudo <- function(inSCE, 
                              useClusters, 
                              pathIndex = NULL, 
                              useReducedDim){
  
  starter <- colData(inSCE)$"TSCAN_clusters" == useClusters
  inSCE <- inSCE[,starter]
  results <- getTSCANResults(inSCE, analysisName = "Pseudotime")
  
  line.data <- TSCAN::reportEdges(results$by.cluster, mst = results$mst, 
                                  clusters = NULL, use.dimred = useReducedDim) 
  line.data.sub <- line.data[grepl(paste0("^", useClusters, "--"), 
                                   line.data$edge) | 
                               grepl(paste0("--", useClusters, "$"), 
                                     line.data$edge),]
  
  a <- b <- c <- NULL
  a = unlist(line.data.sub[,2])
  b = unlist(line.data.sub[,3])
  c = unlist(line.data.sub[,1])
  
  if (is.null(pathIndex)){
    scater::plotReducedDim(inSCE, 
                           dimred = useReducedDim, 
                           colour_by = I(colData(inSCE)$"TSCAN_pseudotime"), 
                           text_by = "TSCAN_clusters", 
                           text_colour = "black")+
      ggplot2::geom_line(data = line.data.sub, 
                         mapping = ggplot2::aes(x = a, y = b, group = c))
  }
  
  else{
    scater::plotReducedDim(inSCE, 
                           dimred = useReducedDim, 
                           colour_by = paste0("TSCAN_Cluster",
                                              useClusters,
                                              "_Path_",
                                              pathIndex,
                                              "_pseudotime"),  
                           text_by = "TSCAN_clusters", 
                           text_colour = "black") +
      ggplot2::geom_line(data = line.data.sub, 
                         mapping = ggplot2::aes(x = a, y = b, group = c)) 
  }
  
}

###################################################
###  plot gene of interest in branch cluster
###################################################
#' @title Run plotTSCANDEgenes function to plot cells colored by the expression 
#' of a gene of interest
#' @description A wrapper function which plots all the cells in the cluster 
#' containing the branch point of the MST in the dataset. Each point is a cell 
#' colored by the expression of a gene of interest and the relevant edges of 
#' the MST are overlaid on top.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param geneSymbol Choose the gene of interest from the DE genes in order to 
#' know the level of expression of gene in clusters.
#' @param useClusters Choose the cluster containing the branch point in the data
#' in order to recompute the pseudotimes so that the root lies at the cluster 
#' center, allowing us to detect genes that are associated with the divergence 
#' of the branches.
#' @param useReducedDim Saved dimension reduction name in \code{inSCE}. 
#' Required.
#'
#' @return A plots with the cells colored by the expression of a gene of 
#' interest.
#' @export 
#' @author Nida Pervaiz
#' @importFrom utils head
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rowData(sce)$Symbol <- rowData(sce)$feature_name
#' rownames(sce) <- rowData(sce)$Symbol
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
#'                     useAssay = "logcounts", reducedDimName = "PCA")
#' sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", 
#'                     reducedDimName = "TSNE")
#' sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
#' sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
#' sce <- runTSCANClusterDEAnalysis(inSCE = sce, useClusters = 5)
#' plotTSCANDEgenes(inSCE = sce, geneSymbol = "CD74", useReducedDim = "TSNE")

plotTSCANDEgenes <- function(inSCE, 
                             geneSymbol, 
                             useClusters=NULL, 
                             useReducedDim){
  
  if(!is.null(useClusters)){
    discardClusters <- setdiff(c(unique(colData(inSCE)$"TSCAN_clusters")), 
                               useClusters)
    for (i in seq_along(discardClusters)){
      if (i == 1){
        nx <- inSCE[, colData(inSCE)$"TSCAN_clusters" != discardClusters[i]]
      }
      else{
        nx <- nx[, colData(nx)$"TSCAN_clusters" != discardClusters[i]]
        
      }
    }
    
  }
  else{
    nx <- inSCE
  }
  scater::plotReducedDim(nx, 
                         dimred = useReducedDim, 
                         colour_by = geneSymbol, 
                         swap_rownames="Symbol", 
                         text_by = "TSCAN_clusters", 
                         text_colour = "black")
  
}
