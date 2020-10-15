#' Uniform Manifold Approximation and Projection(UMAP) algorithm for
#' dimension reduction.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Assay to use for UMAP computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"counts"}.
#' @param useAltExp The subset to use for UMAP computation, usually for the
#' selected.variable features. Default \code{NULL}.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' If given a single character, will take the annotation from
#' \code{\link[SummarizedExperiment]{colData}}. Default \code{NULL}.
#' @param reducedDimName A name to store the results of the dimension reduction
#' coordinates obtained from this method. Default \code{"UMAP"}.
#' @param logNorm Whether the counts will need to be log-normalized prior to
#' generating the UMAP via \code{\link[scater]{logNormCounts}}. Default
#' \code{TRUE}.
#' @param nNeighbors The size of local neighborhood used for manifold
#' approximation. Larger values result in more global views of the manifold,
#' while smaller values result in more local data being preserved. Default
#' \code{30}. See `?uwot::umap` for more information.
#' @param nIterations The number of iterations performed during layout
#' optimization. Default is \code{200}.
#' @param alpha The initial value of "learning rate" of layout optimization.
#' Default is \code{1}.
#' @param minDist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where nearby
#' points on the manifold are drawn closer together, while larger values will
#' result on a more even dispersal of points. Default \code{0.01}. See
#' `?uwot::umap` for more information.
#' @param spread The effective scale of embedded points. In combination with
#' ‘min_dist’, this determines how clustered/clumped the embedded points are.
#' Default \code{1}. See `?uwot::umap` for more information.
#' @param pca Logical. Whether to perform dimension reduction with PCA before
#' UMAP. Default \code{TRUE}
#' @param initialDims  Number of dimensions from PCA to use as input in UMAP.
#' Default \code{50}.
#' @return A \linkS4class{SingleCellExperiment} object with UMAP computation
#' updated in \code{reducedDim(inSCE, reducedDimName)}.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' umap_res <- getUMAP(inSCE = sce, useAssay = "counts",
#'                     reducedDimName = "UMAP", logNorm = TRUE,
#'                     nNeighbors = 30, alpha = 1,
#'                     nIterations = 200, spread = 1, pca = TRUE,
#'                     initialDims = 50)
#' reducedDims(umap_res)
getUMAP <- function(inSCE, useAssay = "counts", useAltExp = NULL,
                    sample = NULL, reducedDimName = "UMAP", logNorm = TRUE,
                    nNeighbors = 30, nIterations = 200, alpha = 1,
                    minDist = 0.01, spread = 1, pca = TRUE,
                    initialDims = 50) {
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  if (!is.null(useAltExp)) {
    if (!(useAltExp %in% SingleCellExperiment::altExpNames(inSCE))) {
      stop("Specified altExp '", useAltExp, "' not found. ")
    }
    sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
    if (!(useAssay %in% SummarizedExperiment::assayNames(sce))) {
      stop("Specified assay '", useAssay, "' not found in the ",
           "specified altExp. ")
    }
  } else {
    if (!(useAssay %in% SummarizedExperiment::assayNames(inSCE))) {
      stop("Specified assay '", useAssay, "' not found. ")
    }
    sce <- inSCE
  }

  if(!is.null(sample)) {
    if (is.character(sample) && length(sample) == 1) {
      if (!sample %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("Given sample annotation '", sample, "' not found.")
      }
      sample <- SummarizedExperiment::colData(inSCE)[[sample]]
    } else if (length(sample) > 1){
      if(length(sample) != ncol(inSCE)) {
        stop("'sample' must be the same length as the number of columns in 'inSCE'")
      }
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }
  samples <- unique(sample)
  umapDims = matrix(nrow = ncol(inSCE), ncol = 2)
  for (i in seq_along(samples)){
    useAssayTemp <- useAssay
    sceSampleInd <- sample == samples[i]
    sceSample <- sce[, sceSampleInd]
    if(logNorm){
      sceSample <- scater_logNormCounts(sceSample, useAssay = useAssay)
      useAssayTemp = "ScaterLogNormCounts"
    }

    matColData <- SummarizedExperiment::assay(sceSample, useAssayTemp)
    matColData <- as.matrix(matColData)

    if (isTRUE(pca)) {
      if(initialDims > ncol(matColData)){
        doPCA <- ncol(matColData)
      }else{
        doPCA <- initialDims
      }
    } else {
      doPCA <- NULL
    }
    if(nNeighbors > ncol(matColData)){
      nNeighbors <- ncol(matColData)
    }

    umapRes <- uwot::umap(t(matColData), n_neighbors = nNeighbors,
                          learning_rate = alpha,
                          min_dist = minDist, spread = spread,
                          n_sgd_threads = 1, pca = doPCA,
                          n_epochs = nIterations)
    if (is.null(rownames(sceSample))) {
      rownames(umapRes) <- colnames(sceSample)
    }
    umapDims[sceSampleInd, ] <- umapRes
  }
  colnames(umapDims) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- umapDims
  return(inSCE)
}
