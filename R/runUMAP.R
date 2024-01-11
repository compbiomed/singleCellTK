#' Run UMAP embedding with scater method
#' @rdname runUMAP
#' @description Uniform Manifold Approximation and Projection (UMAP) algorithm
#' is commonly for 2D visualization of single-cell data. These functions wrap 
#' the scater \code{\link[scater]{calculateUMAP}} function.
#' 
#' Users can use \code{runQuickUMAP} to directly create UMAP embedding from raw
#' count matrix, with necessary preprocessing including normalization, variable
#' feature selection, scaling, dimension reduction all automated. Therefore, 
#' \code{useReducedDim} is disabled for \code{runQuickUMAP}. 
#' 
#' In a complete analysis, we still recommend having dimension reduction such as
#' PCA created beforehand and select proper numbers of dimensions for using
#' \code{runUMAP}, so that the result can match with the clustering based on the
#' same input PCA. 
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useReducedDim The low dimension representation to use for UMAP
#' computation. If \code{useAltExp} is specified, \code{useReducedDim} has to
#' exist in \code{reducedDims(altExp(inSCE, useAltExp))}. Default \code{"PCA"}.
#' @param useAssay Assay to use for UMAP computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Ignored when using
#' \code{useReducedDim}. Default \code{NULL}.
#' @param useAltExp The subset to use for UMAP computation, usually for the
#' selected variable features. Default \code{NULL}.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' If given a single character, will take the annotation from \code{colData}.
#' Default \code{NULL}.
#' @param reducedDimName A name to store the results of the UMAP embedding
#' coordinates obtained from this method. Default \code{"UMAP"}.
#' @param logNorm Whether the counts will need to be log-normalized prior to
#' generating the UMAP via \code{\link{scaterlogNormCounts}}. Ignored when using
#' \code{useReducedDim}. Default \code{TRUE}.
#' @param useFeatureSubset Subset of feature to use for dimension reduction. A
#' character string indicating a \code{rowData} variable that stores the logical
#' vector of HVG selection, or a vector that can subset the rows of
#' \code{inSCE}. Default \code{NULL}.
#' @param nTop Automatically detect this number of variable features to use for
#' dimension reduction. Ignored when using \code{useReducedDim} or using
#' \code{useFeatureSubset}. Default \code{2000}.
#' @param scale Whether \code{useAssay} matrix will need to be standardized.
#' Default \code{TRUE}.
#' @param pca Logical. Whether to perform dimension reduction with PCA before
#' UMAP. Ignored when using \code{useReducedDim}. Default \code{TRUE}.
#' @param initialDims Number of dimensions from PCA to use as input in UMAP.
#' Default \code{25}.
#' @param nNeighbors The size of local neighborhood used for manifold
#' approximation. Larger values result in more global views of the manifold,
#' while smaller values result in more local data being preserved. Default
#' \code{30}. See \code{\link[scater]{calculateUMAP}} for more information.
#' @param nIterations The number of iterations performed during layout
#' optimization. Default is \code{200}.
#' @param alpha The initial value of "learning rate" of layout optimization.
#' Default is \code{1}.
#' @param minDist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where nearby
#' points on the manifold are drawn closer together, while larger values will
#' result on a more even dispersal of points. Default \code{0.01}. See
#' \code{\link[scater]{calculateUMAP}} for more information.
#' @param spread The effective scale of embedded points. In combination with
#' \code{minDist}, this determines how clustered/clumped the embedded points
#' are. Default \code{1}. See \code{\link[scater]{calculateUMAP}} for more
#' information.
#' @param seed Random seed for reproducibility of UMAP results.
#' Default \code{NULL} will use global seed in use by the R environment.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' the PCA should be parallelized.
#' @param verbose Logical. Whether to print log messages. Default \code{TRUE}.
#' @return A \linkS4class{SingleCellExperiment} object with UMAP computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' # Run from raw counts
#' sce <- runQuickUMAP(sce)
#' plotDimRed(sce, "UMAP")
#' 
#' @importFrom S4Vectors metadata<-
#' @importFrom BiocParallel SerialParam
runUMAP <- function(inSCE, useReducedDim = "PCA", useAssay = NULL, 
                    useAltExp = NULL, sample = NULL, reducedDimName = "UMAP",
                    logNorm = TRUE, useFeatureSubset = NULL, nTop = 2000,
                    scale = TRUE, pca = TRUE, initialDims = 25, nNeighbors = 30,
                    nIterations = 200, alpha = 1, minDist = 0.01, spread = 1,
                    seed = 12345, verbose = TRUE, BPPARAM = SerialParam()) {
  params <- as.list(environment())
  params$inSCE <- NULL
  params$BPPARAM <- NULL
  # Note: useMat = list(useAssay = useAssay, ...)
  # `.selectSCEMatrix()` does the check and corrects the conditions
  # When using useReducedDim, useAssay will be ignored
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                             useReducedDim = useReducedDim,
                             useAltExp = useAltExp,
                             returnMatrix = FALSE)$names
  params$useAssay <- useMat$useAssay
  useAssay <- useMat$useAssay
  if (!is.null(useAltExp)) sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
  else sce <- inSCE
  sample <- .manageCellVar(inSCE, sample)
  if (is.null(sample)) sample <- rep(1, ncol(inSCE))
  samples <- unique(sample)
  umapDims <- matrix(nrow = ncol(inSCE), ncol = 2)
  # Start looping for each sample
  for (i in seq_along(samples)){
    sceSampleInd <- sample == samples[i]
    sceSample <- sce[, sceSampleInd]
    # Do logNorm
    useAssayTemp <- useAssay
    if (!is.null(useAssay)) {
      if (isTRUE(logNorm)) {
        sceSample <- scaterlogNormCounts(sceSample, useAssay = useAssay)
        useAssayTemp = "ScaterLogNormCounts"
      }
    }
    subset_row <- .parseUseFeatureSubset(inSCE, useFeatureSubset,
                                         altExpObj = sce, returnType = "logic")
    # initialDims passed to `pca` of scran::calculateUMAP.
    if (!isTRUE(pca) & !is.null(useAssay)) {
      initialDims <- NULL
    }
    
    nNeighbors <- min(ncol(sceSample), nNeighbors)
    if(isTRUE(verbose)) {
      p <- paste0(date(), " ... Computing Scater UMAP for sample '",
                  samples[i], "'.")
      message(p)
    }
    
    umapRes <- withr::with_seed(seed, {
      scater::calculateUMAP(sceSample, exprs_values = useAssayTemp,
                                       dimred = useReducedDim, scale = scale,
                                       n_neighbors = nNeighbors,
                                       learning_rate = alpha,
                                       min_dist = minDist, spread = spread,
                                       n_sgd_threads = 1, pca = initialDims,
                                       n_epochs = nIterations,
                                       subset_row = subset_row, ntop = nTop,
                                       BPPARAM = BPPARAM)
    })
    if (is.null(rownames(sceSample))) {
      rownames(umapRes) <- colnames(sceSample)
    }
    umapDims[sceSampleInd, ] <- umapRes
  }
  colnames(umapDims) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- umapDims
  metadata(inSCE)$sctk$runDimReduce$embedding[[reducedDimName]] <- params
  return(inSCE)
}

#' @rdname runUMAP
#' @param ... Parameters passed to \code{runUMAP}
#' @importFrom BiocParallel SerialParam
#' @export
runQuickUMAP <- function(inSCE, useAssay = "counts", sample = "sample", ...) {
  args <- list(...)
  if (!is.null(args$useReducedDim)) {
    warning("Forcing `useReducedDim` to be `NULL`. Please use `runUMAP` for ",
            "using reducedDim.")
  }
  # Here `useReducedDim` entry is removed from list
  # Later, we add the entry with value `NULL`
  args$useReducedDim <- NULL
  args <- c(list(inSCE = inSCE, useAssay = useAssay, useReducedDim = NULL, 
                 sample = sample), args)
  inSCE <- do.call("runUMAP", args = args)
  return(inSCE)
}

#' @rdname runUMAP
#' @export
#' @importFrom BiocParallel SerialParam
getUMAP <- function(inSCE, useReducedDim = "PCA", useAssay = NULL, 
                    useAltExp = NULL, sample = NULL, reducedDimName = "UMAP",
                    logNorm = TRUE, useFeatureSubset = NULL, nTop = 2000,
                    scale = TRUE, pca = TRUE, initialDims = 25, nNeighbors = 30,
                    nIterations = 200, alpha = 1, minDist = 0.01, spread = 1,
                    seed = 12345, BPPARAM = SerialParam()) {
  .Deprecated("runUMAP")
  runUMAP(inSCE, useReducedDim = useReducedDim, useAssay = useAssay, 
          useAltExp = useAltExp, sample = sample, 
          reducedDimName = reducedDimName,
          logNorm = logNorm, useFeatureSubset = useFeatureSubset, nTop = nTop,
          scale = scale, pca = pca, initialDims = initialDims, 
          nNeighbors = nNeighbors, nIterations = nIterations, alpha = alpha, 
          minDist = minDist, spread = spread, seed = seed, BPPARAM = BPPARAM)
}
