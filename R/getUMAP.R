#' Run UMAP embedding with scater method
#' @description Uniform Manifold Approximation and Projection (UMAP) algorithm 
#' is commonly for 2D visualization of single-cell data. This function wraps the 
#' scater \code{\link[scater]{calculateUMAP}} function. 
#' 
#' With this funciton, users can create UMAP embedding directly from raw count
#' matrix, with necessary preprocessing including normalization, scaling, 
#' dimension reduction all automated. Yet we still recommend having the PCA as 
#' input, so that the result can match with the clustering based on the same 
#' input PCA.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for UMAP computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Ignored when using 
#' \code{useReducedDim}. Default \code{"logcounts"}.
#' @param useReducedDim The low dimension representation to use for UMAP
#' computation. If \code{useAltExp} is specified, \code{useReducedDim} has to 
#' exist in \code{reducedDims(altExp(inSCE, useAltExp))}. Default \code{NULL}.
#' @param useAltExp The subset to use for UMAP computation, usually for the
#' selected variable features. Default \code{NULL}.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' If given a single character, will take the annotation from \code{colData}. 
#' Default \code{NULL}.
#' @param reducedDimName A name to store the results of the UMAP embedding
#' coordinates obtained from this method. Default \code{"UMAP"}.
#' @param logNorm Whether the counts will need to be log-normalized prior to
#' generating the UMAP via \code{\link{scaterlogNormCounts}}. Ignored when using
#' \code{useReducedDim}. Default \code{FALSE}.
#' @param useHVGList A character string indicating a \code{rowData} variable 
#' that stores the logical vector of HVG selection. Ignored when using
#' \code{useReducedDim}. Default \code{NULL}.
#' @param nTop Automatically detect this number of variable features to use for
#' dimension reduction. Ignored when using \code{useReducedDim} or using 
#' \code{useHVGList}. Default \code{2000}.
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
#' @return A \linkS4class{SingleCellExperiment} object with UMAP computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' # Run from raw counts
#' sce <- getUMAP(inSCE = sce, useAssay = "counts", logNorm = TRUE, nTop = 2000,
#'                scale = TRUE, pca = TRUE)
#' \dontrun{
#' # Run from PCA
#' sce <- scaterlogNormCounts(sce, "logcounts")
#' sce <- runModelGeneVar(sce)
#' sce <- scaterPCA(sce, useAssay = "logcounts", 
#'                  useHVGList = "HVG_modelGeneVar2000", scale = TRUE)
#' sce <- getUMAP(sce, useReducedDim = "PCA")
#' }
getUMAP <- function(inSCE, useAssay = "logcounts", useReducedDim = NULL, 
                    useAltExp = NULL, sample = NULL, reducedDimName = "UMAP", 
                    logNorm = FALSE, useHVGList = NULL, nTop = 2000, scale = TRUE, 
                    pca = TRUE, initialDims = 25, nNeighbors = 30, 
                    nIterations = 200, alpha = 1, minDist = 0.01, spread = 1, 
                    seed = NULL, BPPARAM = BiocParallel::SerialParam()) {
  # input class check
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  # matrix specification check
  # When useReducedDim, never useAssay; 
  # when useAltExp, useAssay/reducedDim from there
  if (is.null(useAssay) && is.null(useReducedDim)) {
    stop("`useAssay` and `useReducedDim` cannot be NULL at the same time.")
  } 
  if (!is.null(useAltExp)) {
    if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
      stop("Specified `useAltExp` not found.")
    }
    sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
  } else {
    sce <- inSCE
  }
  if (!is.null(useReducedDim)) {
    useAssay <- NULL
    if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(sce)) {
      stop("Specified `useReducedDim` not found.")
    }
  } else {
    # In this case, there is definitely `useAssay` specified
    if (!useAssay %in% SummarizedExperiment::assayNames(sce)) {
      stop("Specified `useAssay` not found.")
    }
  }
  # Check sample variable
  if(!is.null(sample)) {
    if (is.character(sample) && length(sample) == 1) {
      if (!sample %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("Given sample annotation '", sample, "' not found.")
      }
      sample <- SummarizedExperiment::colData(inSCE)[[sample]]
    } else if (length(sample) > 1) {
      if(length(sample) != ncol(inSCE)) {
        stop("'sample' must be the same length as the number of columns in ", 
             "'inSCE'")
      }
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }
  samples <- unique(sample)
  umapDims = matrix(nrow = ncol(inSCE), ncol = 2)
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
    if (!is.null(useHVGList) & !is.null(useAssay)) {
      if (!useHVGList %in% colnames(SummarizedExperiment::rowData(inSCE))) {
        stop("Specified HVG list not found")
      }
      hvgs <- rownames(inSCE)[rowSubset(inSCE, useHVGList)]
      subset_row <- rownames(sceSample) %in% hvgs
    } else {
      subset_row <- NULL
    }
    # initialDims passed to `pca` of scran::calculateUMAP. 
    if (!isTRUE(pca) & !is.null(useAssay)) {
      initialDims <- NULL
    }
    
    nNeighbors <- min(ncol(sceSample), nNeighbors)
    message(paste0(date(), " ... Computing Scater UMAP for sample '",
                   samples[i], "'."))
    if (!is.null(seed)) {
      withr::with_seed(
        seed = seed,
        code = umapRes <- scater::calculateUMAP(sceSample, scale = scale,
                                                exprs_values = useAssayTemp,
                                                dimred = useReducedDim,
                                                n_neighbors = nNeighbors,
                                                learning_rate = alpha,
                                                min_dist = minDist, 
                                                spread = spread,
                                                n_sgd_threads = 1, 
                                                pca = initialDims,
                                                n_epochs = nIterations, 
                                                subset_row = subset_row, 
                                                ntop = nTop, BPPARAM = BPPARAM)
      )
    } else {
      umapRes <- scater::calculateUMAP(sceSample, exprs_values = useAssayTemp,
                                       dimred = useReducedDim, scale = scale,
                                       n_neighbors = nNeighbors,
                                       learning_rate = alpha,
                                       min_dist = minDist, spread = spread,
                                       n_sgd_threads = 1, pca = initialDims,
                                       n_epochs = nIterations, 
                                       subset_row = subset_row, ntop = nTop, 
                                       BPPARAM = BPPARAM)
    }
    if (is.null(rownames(sceSample))) {
      rownames(umapRes) <- colnames(sceSample)
    }
    umapDims[sceSampleInd, ] <- umapRes
  }
  colnames(umapDims) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- umapDims
  return(inSCE)
}
