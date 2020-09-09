#' Uniform Manifold Approximation and Projection(UMAP) algorithm for
#' dimension reduction.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "counts".
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param reducedDimName a name to store the results of the dimension reduction
#' coordinates obtained from this method. This is stored in the SingleCellExperiment
#' object in the reducedDims slot. Default "UMAP".
#' @param logNorm Whether the counts will need to be log-normalized prior to
#' generating the UMAP via scater::logNormCounts. Default TRUE.
#' @param nNeighbors The size of local neighborhood used for
#'   manifold approximation. Larger values result in more global
#'   views of the manifold, while smaller values result in more
#'   local data being preserved. Default 30.
#'    See `?uwot::umap` for more information.
#' @param nIterations number of iterations performed during layout optimization.
#' Default is 200.
#' @param alpha initial value of "learning rate" of layout optimization. Default is 1.
#' @param minDist The effective minimum distance between embedded points.
#'    Smaller values will result in a more clustered/clumped
#'    embedding where nearby points on the manifold are drawn
#'    closer together, while larger values will result on a more
#'    even dispersal of points. Default 0.01.
#'    See `?uwot::umap` for more information.
#' @param spread The effective scale of embedded points. In combination with
#'    ‘min_dist’, this determines how clustered/clumped the
#'    embedded points are. Default 1.
#'    See `?uwot::umap` for more information.
#' @param pca Logical. Whether to perform dimensionality reduction with PCA
#' before UMAP.
#' @param initialDims  Number of dimensions from PCA to use as
#' input in UMAP. Default 50.
#'
#' @return A \linkS4class{SingleCellExperiment} object with the reduced
#' dimensions updated under reducedDimName specified.
#' @export
#'
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' umap_res <- getUMAP(inSCE = sce, useAssay = "counts",
#'                     reducedDimName = "UMAP", logNorm = TRUE,
#'                     nNeighbors = 30, alpha = 1,
#'                     nIterations = 200, spread = 1, pca = TRUE,
#'                     initialDims = 50)
#' reducedDims(umap_res)

getUMAP <- function(inSCE, useAssay = "counts",
                    sample = NULL,
                    reducedDimName = "UMAP",
                    logNorm = TRUE,
                    nNeighbors = 30,
                    nIterations = 200,
                    alpha = 1,
                    minDist = 0.01,
                    spread = 1,
                    pca = TRUE,
                    initialDims = 50) {
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  #test for assay existing
    if (!all(useAssay %in% names(assays(inSCE)))){
        stop("assay '", useAssay, "' does not exist.")
    }

    if(!is.null(sample)) {
        if(length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number of columns in 'inSCE'")
        }
    } else {
        sample = rep(1, ncol(inSCE))
    }
    samples <- unique(sample)
    umapDims = matrix(nrow = ncol(inSCE), ncol = 2)
    for (i in seq_len(length(samples))){
        useAssayTemp = useAssay
	sceSampleInd <- sample == samples[i]
        sceSample <- inSCE[, sceSampleInd]
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
        umapDims[sceSampleInd, ] = umapRes
    }
    colnames(umapDims) <- c("UMAP1", "UMAP2")
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- umapDims
    return(inSCE)
}
