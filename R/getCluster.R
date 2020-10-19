#' Get clustering with SNN graph
#' @description Perform SNN graph clustering on a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object, with graph
#' construction by \code{\link[scran]{buildSNNGraph}} and graph clustering by
#' "igraph" package.
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param useAssay A single \code{character}, specifying which
#' \code{\link[SummarizedExperiment]{assay}} to perform the clustering algorithm
#' on. Default \code{NULL}.
#' @param useReducedDim A single \code{character}, specifying which
#' low-dimension representation (\code{\link[SingleCellExperiment]{reducedDim}})
#' to perform the clustering algorithm on. Default \code{NULL}.
#' @param useAltExp A single \code{character}, specifying the assay which
#' \code{\link[SingleCellExperiment]{altExp}} to perform the clustering
#' algorithm on. Default \code{NULL}.
#' @param altExpAssay A single \code{character}, specifying which
#' \code{\link[SummarizedExperiment]{assay}} in the chosen
#' \code{\link[SingleCellExperiment]{altExp}} to work on. Only used when
#' \code{useAltExp} is set. Default \code{"counts"}.
#' @param clusterName A single \code{character}, specifying the name to store
#' the cluster label in \code{\link[SummarizedExperiment]{colData}}. Default
#' \code{"scranSNN_cluster"}.
#' @param k An \code{integer}, the number of nearest neighbors used to construct
#' the graph. Smaller value indicates higher resolution and larger number of
#' clusters. Default \code{10}.
#' @param nComp An \code{integer}, the number of components to use when
#' \code{useAssay} or \code{useAltExp} is specified. WON'T work with
#' \code{useReducedDim}. Default \code{50}.
#' @param weightType A single \code{character}, that specifies the edge weighing
#' scheme when constructing the Shared Nearest-Neighbor (SNN) graph. Choose from
#' \code{"rank"}, \code{"number"}, \code{"jaccard"}. Default \code{"rank"}.
#' @param algorithm A single \code{character}, that specifies the community
#' detection algorithm to work on the SNN graph. Choose from \code{"walktrap"},
#' \code{"louvain"}, \code{"infomap"}, \code{"fastGreedy"}, \code{"labelProp"},
#' \code{"leadingEigen"}. Default \code{"walktrap"}.
#' @return The input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object with \code{factor} cluster labeling updated in
#' \code{colData(inSCE)[[clusterName]]}.
#' @references Aaron Lun and et. al., 2016
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- runScranSNN(mouseBrainSubsetSCE,
#'                                    useReducedDim = "PCA_logcounts")
runScranSNN <- function(inSCE, useAssay = NULL, useReducedDim = NULL,
                        useAltExp = NULL, altExpAssay = "counts",
                        clusterName = "scranSNN_cluster",
                        k = 10, nComp = 50,
                        weightType = c("rank", "number", "jaccard"),
                        algorithm = c("walktrap", "louvain", "infomap",
                                      "fastGreedy", "labelProp",
                                      "leadingEigen")) {
  # TODO: parallele parameter
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("SCE")
  }
  if (is.null(useAssay) + is.null(useReducedDim) + is.null(useAltExp) != 2) {
    stop("Scran SNN clustering requires one and only one of 'useAssay', ",
         "'useReducedDim', and 'useAltExp'.")
  }
  graphClustAlgoList = list(walktrap = igraph::cluster_walktrap,
                            louvain = igraph::cluster_louvain,
                            infomap = igraph::cluster_infomap,
                            fastGreedy = igraph::cluster_fast_greedy,
                            labelProp = igraph::cluster_label_prop,
                            leadingEigen = igraph::cluster_leading_eigen)
  if (all(algorithm == c("walktrap", "louvain", "infomap", "fastGreedy",
                         "labelProp", "leadingEigen"))) {
    algorithm = "walktrap"
  } else if (!algorithm %in% names(graphClustAlgoList)){
    stop("Specified algorithm '", algorithm, "' not supported")
  }

  message(paste0(date(), " ... Running 'scran SNN clustering'"))
  if (!is.null(useAssay)){
    if (!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
      stop("Specified assay '", useAssay, "' not found.")
    }
    g <- scran::buildSNNGraph(x = inSCE, k = k, assay.type = useAssay,
                              d = nComp, type = weightType, use.dimred = NULL)
  } else if (!is.null(useReducedDim)) {
    if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
      stop("Specified reducedDim '", useReducedDim, "' not found.")
    }
    g <- scran::buildSNNGraph(x = inSCE, k = k, use.dimred = useReducedDim,
                              type = weightType)
  } else if (!is.null(useAltExp)) {
    if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
      stop("Specified altExp '", useAltExp, "' not found.")
    }
    ae <- SingleCellExperiment::altExp(inSCE, useAltExp)
    if (!altExpAssay %in% SummarizedExperiment::assayNames(ae)) {
      stop("altExpAssay: '", altExpAssay, "' not in specified altExp.")
    }
    g <- scran::buildSNNGraph(x = ae, k = k, assay.type = altExpAssay,
                              d = nComp, type = weightType, use.dimred = NULL)
  }

  clustFunc = graphClustAlgoList[[algorithm]]
  clust <- clustFunc(g)$membership
  SummarizedExperiment::colData(inSCE)[[clusterName]] <- factor(clust)
  return(inSCE)
}

#' Get clustering with KMeans
#' @description Perform KMeans clustering on a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object, with
#' \code{\link[stats]{kmeans}}.
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param useReducedDim A single \code{character}, specifying which
#' low-dimension representation to perform the clustering algorithm on. Default
#' \code{"PCA"}.
#' @param clusterName A single \code{character}, specifying the name to store
#' the cluster label in \code{\link[SummarizedExperiment]{colData}}. Default
#' \code{"scranSNN_cluster"}.
#' @param nCenters An \code{integer}, the number of centroids (clusters).
#' @param nIter An \code{integer}, the maximum number of iterations allowed.
#' Default \code{10}.
#' @param nStart An \code{integer}, the number of random sets to choose. Default
#' \code{1}.
#' @param seed An \code{integer}. The seed for the random number generator.
#' Default \code{12345}.
#' @param algorithm A single \code{character}. Choose from
#' \code{"Hartigan-Wong"}, \code{"Lloyd"}, \code{"MacQueen"}. May be
#' abbreviated. Default \code{"Hartigan-Wong"}.
#' @return The input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object with \code{factor} cluster labeling updated in
#' \code{colData(inSCE)[[clusterName]]}.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- runKMeans(mouseBrainSubsetSCE,
#'                                  useReducedDim = "PCA_logcounts",
#'                                  nCenters = 2)
runKMeans <- function(inSCE, useReducedDim = "PCA",
                      clusterName = "KMeans_cluster", nCenters, nIter = 10,
                      nStart = 1, seed = 12345,
                      algorithm = c("Hartigan-Wong", "Lloyd", "MacQueen")){
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("SCE")
  }
  if (is.null(useReducedDim)) {
    stop("runKMeans requires 'useReducedDim' input.")
  }
  if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
    stop("Specified reducedDim '", useReducedDim, "' not found.")
  }
  if(all(algorithm == c("Hartigan-Wong", "Lloyd", "MacQueen"))){
    algorithm = "Hartigan-Wong"
  }
  message(paste0(date(), " ... Running 'KMeans clustering'"))
  mat <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)
  
  clust.kmeans <- withr::with_seed(seed = seed,
                {stats::kmeans(mat, centers = nCenters, iter.max = nIter,
                        nstart = nStart, algorithm = algorithm)})
  
  clust.kmeans <- factor(clust.kmeans$cluster)
  SummarizedExperiment::colData(inSCE)[[clusterName]] <- clust.kmeans
  return(inSCE)
}
