#' Get clustering with SNN graph
#' @description Perform SNN graph clustering on a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object, with graph
#' construction by \code{\link[scran]{buildSNNGraph}} and graph clustering by
#' "igraph" package.
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param useAssay A single \code{character}, specifying which
#' \code{\link{assay}} to perform the clustering algorithm
#' on. Default \code{NULL}.
#' @param useReducedDim A single \code{character}, specifying which
#' low-dimension representation (\code{\link{reducedDim}})
#' to perform the clustering algorithm on. Default \code{NULL}.
#' @param useAltExp A single \code{character}, specifying the assay which
#' \code{\link{altExp}} to perform the clustering
#' algorithm on. Default \code{NULL}.
#' @param altExpAssay A single \code{character}, specifying which
#' \code{\link{assay}} in the chosen
#' \code{\link{altExp}} to work on. Only used when
#' \code{useAltExp} is set. Default \code{"counts"}.
#' @param altExpRedDim A single \code{character}, specifying which
#' \code{\link{reducedDim}} within the \code{\link{altExp}} specified by
#' \code{useAltExp} to use. Only used when \code{useAltExp} is set. Default
#' \code{NULL}.
#' @param clusterName A single \code{character}, specifying the name to store
#' the cluster label in \code{\link{colData}}. Default
#' \code{"scranSNN_cluster"}.
#' @param k An \code{integer}, the number of nearest neighbors used to construct
#' the graph. Smaller value indicates higher resolution and larger number of
#' clusters. Default \code{8}.
#' @param nComp An \code{integer}, the number of components to use when
#' \code{useAssay} or \code{useAltExp} is specified. WON'T work with
#' \code{useReducedDim}. Default \code{50}.
#' @param weightType A single \code{character}, that specifies the edge weighing
#' scheme when constructing the Shared Nearest-Neighbor (SNN) graph. Choose from
#' \code{"rank"}, \code{"number"}, \code{"jaccard"}. Default \code{"rank"}.
#' @param algorithm A single \code{character}, that specifies the community
#' detection algorithm to work on the SNN graph. Choose from \code{"leiden"}, 
#' \code{"louvain"}, \code{"walktrap"}, \code{"infomap"}, \code{"fastGreedy"}, 
#' \code{"labelProp"}, \code{"leadingEigen"}. Default \code{"louvain"}. See 
#' Detail.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object to use 
#' for processing the SNN graph generation step in parallel.
#' @param seed Random seed for reproducibility of results. Default \code{NULL} 
#' will use global seed in use by the R environment.
#' @param ... Other optional parameters passed to the \code{\link{igraph}} 
#' clustering functions. See Details.
#' @return The input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object with \code{factor} cluster labeling updated in
#' \code{colData(inSCE)[[clusterName]]}.
#' @details Different graph based clustering algorithms have diverse sets of 
#' parameters that users can tweak. The help information can be found here:
#' \itemize{
#'  \item{for \code{"louvain"}, see function help 
#'  \code{\link[igraph]{cluster_louvain}}}
#'  \item{for \code{"leiden"}, see function help 
#'  \code{\link[igraph]{cluster_leiden}}}
#'  \item{for \code{"walktrap"}, see function help 
#'  \code{\link[igraph]{cluster_walktrap}}}
#'  \item{for \code{"infomap"}, see function help 
#'  \code{\link[igraph]{cluster_infomap}}}
#'  \item{for \code{"fastGreedy"}, see function help 
#'  \code{\link[igraph]{cluster_fast_greedy}}}
#'  \item{for \code{"labelProp"}, see function help 
#'  \code{\link[igraph]{cluster_label_prop}}}
#'  \item{for \code{"leadingEigen"}, see function help 
#'  \code{\link[igraph]{cluster_leading_eigen}}}
#' }
#' @references Aaron Lun and et. al., 2016
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- runScranSNN(mouseBrainSubsetSCE,
#'                                    useReducedDim = "PCA_logcounts")
runScranSNN <- function(inSCE, useAssay = NULL, useReducedDim = NULL,
                        useAltExp = NULL, altExpAssay = "counts",
                        altExpRedDim = NULL,
                        clusterName = "scranSNN_cluster",
                        k = 8, nComp = 50, 
                        weightType = c("rank", "number", "jaccard"),
                        algorithm = c("louvain", "leiden", "walktrap", 
                                      "infomap", "fastGreedy", "labelProp",
                                      "leadingEigen"),
                        BPPARAM = BiocParallel::SerialParam(),
                        seed = 12345,
                        ...) {
  if (!is.null(seed)) {
    # If set seed, run the whole function again but enter the next condition
    # Not sure if this is computationally inefficient. 
    withr::with_seed(
      seed = seed,
      code = runScranSNN(inSCE, useAssay = useAssay, 
                         useReducedDim = useReducedDim, useAltExp = useAltExp, 
                         altExpAssay = altExpAssay, altExpRedDim = altExpRedDim,
                         clusterName = clusterName, k = k, nComp = nComp, 
                         weightType = weightType, algorithm = algorithm,
                         BPPARAM = BPPARAM, seed = NULL))
  } else {
    if (!inherits(inSCE, "SingleCellExperiment")) {
      stop("'inSCE' should be a SingleCellExperiment object.")
    }
    if (is.null(useAssay) + is.null(useReducedDim) + is.null(useAltExp) != 2) {
      stop("Scran SNN clustering requires one and only one of 'useAssay', ",
           "'useReducedDim', and 'useAltExp'.")
    }
    weightType <- match.arg(weightType)
    algorithm <- match.arg(algorithm)
    
    graphClustAlgoList = list(leiden = igraph::cluster_leiden,
                              walktrap = igraph::cluster_walktrap,
                              louvain = igraph::cluster_louvain,
                              infomap = igraph::cluster_infomap,
                              fastGreedy = igraph::cluster_fast_greedy,
                              labelProp = igraph::cluster_label_prop,
                              leadingEigen = igraph::cluster_leading_eigen)
    message(paste0(date(), " ... Running 'scran SNN clustering'"))
    if (!is.null(useAssay)){
      if (!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop("Specified assay '", useAssay, "' not found.")
      }
      g <- scran::buildSNNGraph(x = inSCE, k = k, assay.type = useAssay,
                                d = nComp, type = weightType, use.dimred = NULL,
                                BPPARAM = BPPARAM)
    } else if (!is.null(useReducedDim)) {
      if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
        stop("Specified reducedDim '", useReducedDim, "' not found.")
      }
      g <- scran::buildSNNGraph(x = inSCE, k = k, use.dimred = useReducedDim,
                                type = weightType, BPPARAM = BPPARAM)
    } else if (!is.null(useAltExp)) {
      if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
        stop("Specified altExp '", useAltExp, "' not found.")
      }
      ae <- SingleCellExperiment::altExp(inSCE, useAltExp)
      if (!is.null(altExpRedDim)) {
        if (!altExpRedDim %in% SingleCellExperiment::reducedDimNames(ae)) {
          stop("altExpRedDim: '", altExpRedDim, "' not in specified altExp.")
        }
        g <- scran::buildSNNGraph(x = ae, k = k, use.dimred = altExpRedDim,
                                  type = weightType, BPPARAM = BPPARAM)
      } else {
        if (!altExpAssay %in% SummarizedExperiment::assayNames(ae)) {
          stop("altExpAssay: '", altExpAssay, "' not in specified altExp.")
        }
        g <- scran::buildSNNGraph(x = ae, k = k, assay.type = altExpAssay,
                                  d = nComp, type = weightType, 
                                  use.dimred = NULL, BPPARAM = BPPARAM)
      }
    }
    
    clustFunc = graphClustAlgoList[[algorithm]]
    if (!is.null(seed)) {
      withr::with_seed(
        seed = seed,
        code = clust <- clustFunc(g, ...)$membership)
    } else {
      clust <- clustFunc(g, ...)$membership
    }
    clust <- clustFunc(g, ...)$membership
    SummarizedExperiment::colData(inSCE)[[clusterName]] <- factor(clust)
    return(inSCE)
  }
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
#' the cluster label in \code{\link{colData}}. Default
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
    stop("'inSCE' should be a SingleCellExperiment object.")
  }
  if (is.null(useReducedDim)) {
    stop("runKMeans requires 'useReducedDim' input.")
  }
  if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
    stop("Specified reducedDim '", useReducedDim, "' not found.")
  }
  algorithm <- match.arg(algorithm)
  message(paste0(date(), " ... Running 'KMeans clustering'"))
  mat <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)

  clust.kmeans <- withr::with_seed(seed = seed,
                {stats::kmeans(mat, centers = nCenters, iter.max = nIter,
                        nstart = nStart, algorithm = algorithm)})

  clust.kmeans <- factor(clust.kmeans$cluster)
  SummarizedExperiment::colData(inSCE)[[clusterName]] <- clust.kmeans
  return(inSCE)
}
