#' Apply BBKNN batch effect correction method to SingleCellExperiment object
#'
#' BBKNN, an extremely fast graph-based data integration algorithm. It modifies
#' the neighbourhood construction step to produce a graph that is balanced
#' across all batches of the data.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"BBKNN"}.
#' @param nComponents An integer. Number of principle components or the
#' dimensionality, adopted in the pre-PCA-computation step, the BBKNN step (for
#' how many PCs the algorithm takes into account), and the final UMAP
#' combination step where the value represent the dimensionality of the updated
#' reducedDim. Default \code{50L}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Krzysztof Polanski et al., 2020
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceBatches <- runBBKNN(sceBatches, useAssay = "logcounts",
#'                        nComponents = 10)
#' }
runBBKNN <-function(inSCE, useAssay = 'logcounts', batch = 'batch',
                    reducedDimName = 'BBKNN', nComponents = 50L){
  if(!reticulate::py_module_available(module = "bbknn")){
    warning("Cannot find python module 'bbknn', please install Conda and",
            " run sctkPythonInstallConda() or run ",
            "sctkPythonInstallVirtualEnv(). If one of these have been ",
            "previously run to install the modules, make sure to run ",
            "selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module ",
            "installation. Alternatively, bbknn can be installed on the local ",
            "machine with pip (e.g. pip install bbknn) and then the ",
            "'use_python()' function from the 'reticulate' package can be used",
            " to select the correct Python environment.")
    return(inSCE)
  }
  useAssay <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                               returnMatrix = FALSE)$names$useAssay
  batchVec <- .manageCellVar(inSCE, batch, as.factor = FALSE)
  reducedDimName <- gsub(' ', '_', reducedDimName)
  nComponents <- as.integer(nComponents)
  ## Run algorithm
  adata <- .sce2adata(inSCE, useAssay = useAssay)
  adata$obs[["batch"]] <- batchVec
  sc$tl$pca(adata, n_comps = nComponents)
  bbknn$bbknn(adata, batch_key = "batch", n_pcs = nComponents)
  sc$tl$umap(adata, n_components = nComponents)
  bbknnUmap <- adata$obsm[["X_umap"]]
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- bbknnUmap
  S4Vectors::metadata(inSCE)$batchCorr[[reducedDimName]] <- list(useAssay = useAssay,
                                                            origLogged = TRUE,
                                                            method = "BBKNN",
                                                            matType = "reducedDim",
                                                            batch = batch)
  return(inSCE)
}

#' Apply ComBat-Seq batch effect correction method to SingleCellExperiment
#' object
#'
#' The ComBat-Seq batch adjustment approach assumes that batch effects represent
#' non-biological but systematic shifts in the mean or variability of genomic
#' features for all samples within a processing batch. It uses either parametric
#' or non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects.
#'
#' For the parameters \code{covariates} and \code{useSVA}, when the cell type
#' information is known, it is recommended to specify the cell type annotation
#' to the argument \code{covariates}; if the cell types are unknown but
#' expected to be balanced, it is recommended to run with default settings, yet
#' informative covariates could still be useful. If the cell types are unknown
#' and are expected to be unbalanced, it is recommended to set \code{useSVA}
#' to \code{TRUE}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"counts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param covariates A character vector indicating the fields in
#' \code{\link[SummarizedExperiment]{colData}} that annotates other covariates,
#' such as the cell types. Default \code{NULL}.
#' @param bioCond A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the biological
#' conditions. Default \code{NULL}.
#' @param useSVA A logical scalar. Whether to estimate surrogate variables and
#' use them as an empirical control. Default \code{FALSE}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"ComBat"}.
#' @param shrink A logical scalar. Whether to apply shrinkage on parameter
#' estimation. Default \code{FALSE}.
#' @param shrinkDisp A logical scalar. Whether to apply shrinkage on dispersion.
#' Default \code{FALSE}.
#' @param nGene An integer. Number of random genes to use in empirical Bayes
#' estimation, only useful when \code{shrink} is set to \code{TRUE}. Default
#' \code{NULL}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceBatches <- sample(sceBatches, 40)
#' # Cell type known
#' sceBatches <- runComBatSeq(sceBatches, "counts", "batch",
#'                            covariates = "cell_type",
#'                            assayName = "ComBat_cell_seq")
#' # Cell type unknown but balanced
#' #sceBatches <- runComBatSeq(sceBatches, "counts", "batch",
#' #                           assayName = "ComBat_seq")
#' # Cell type unknown and unbalanced
#' #sceBatches <- runComBatSeq(sceBatches, "counts", "batch",
#' #                           useSVA = TRUE,
#' #                           assayName = "ComBat_sva_seq")
#' @export
runComBatSeq <- function(inSCE, useAssay = "counts", batch = 'batch',
                         covariates = NULL, bioCond = NULL, useSVA = FALSE,
                         assayName = "ComBatSeq", shrink = FALSE,
                         shrinkDisp = FALSE, nGene = NULL) {
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(any(!c(batch, covariates, bioCond) %in%
         names(SummarizedExperiment::colData(inSCE)))){
    anns <- c(batch, covariates, bioCond)
    notFound <- which(!anns %in% names(SummarizedExperiment::colData(inSCE)))
    notFound <- anns[notFound]
    stop("\"annotation\" name:", paste(notFound, collapse = ', '), "not found")
  }
  if (!is.null(covariates) &&
      isTRUE(useSVA)) {
    stop("Only one of 'covariates' and 'useSVA' can be specified.")
  }
  # Modeling
  annot <- data.frame(SummarizedExperiment::colData(inSCE))
  if (!is.null(covariates)) {
    # Do ComBat-Cell-Seq
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = annot)
  } else if (isTRUE(useSVA)) {
    # Do ComBat-SVA-Seq
    mod.batch = stats::model.matrix(
      stats::as.formula(paste0("~", batch)),
      data = annot)
    mod0 = stats::model.matrix(~1, data = annot)
    svobj_filt = sva::svaseq(as.matrix(assay(inSCE, useAssay) + 2.5),
                             mod.batch, mod0)
    nSV <- dim(svobj_filt$sv)[2]
    for (i in seq_len(nSV)) {
      annot[[paste0("SV", i)]] <- svobj_filt$sv[,i]
    }
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(paste0("SV", seq_len(nSV)),
                                           collapse = "+"))),
      data = annot)
  } else {
    mod <- NULL
  }
  # Running
  mat <- as.matrix(SummarizedExperiment::assay(inSCE, useAssay))
  batchCol <- annot[[batch]]
  if (!is.null(bioCond)) {
    groupCol <- annot[[bioCond]]
  } else {
    groupCol <- NULL
  }
  mat.corr <-  sva::ComBat_seq(counts = mat, batch = batchCol, group = groupCol,
                               covar_mod = mod, shrink = shrink,
                               shrink.disp = shrinkDisp, gene.subset.n = nGene)
  expData(inSCE, assayName, tag = "batchCorrected", altExp = FALSE) <- mat.corr
  S4Vectors::metadata(inSCE)$batchCorr[[assayName]] <- list(useAssay = useAssay,
                                                            origLogged = FALSE,
                                                            method = "ComBatSeq",
                                                            matType = "assay",
                                                            batch = batch,
                                                            condition = covariates)
  return(inSCE)
}

#' Apply a fast version of the mutual nearest neighbors (MNN) batch effect
#' correction method to SingleCellExperiment object
#'
#' fastMNN is a variant of the classic MNN method, modified for speed and more
#' robust performance. For introduction of MNN, see \code{\link{runMNNCorrect}}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param useReducedDim A single character indicating the dimension reduction
#' used for batch correction. Will ignore \code{useAssay} when using.
#' Default \code{NULL}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Default \code{"fastMNN"}.
#' @param k An integer scalar specifying the number of nearest neighbors to
#' consider when identifying MNNs. See "See Also". Default \code{20}.
#' @param propK A numeric scalar in (0, 1) specifying the proportion of cells in
#' each dataset to use for mutual nearest neighbor searching. See "See Also".
#' Default \code{NULL}.
#' @param ndist A numeric scalar specifying the threshold beyond which
#' neighbours are to be ignored when computing correction vectors. See "See
#' Also". Default \code{3}.
#' @param minBatchSkip Numeric scalar specifying the minimum relative magnitude
#' of the batch effect, below which no correction will be performed at a given
#' merge step. See "See Also". Default \code{0}.
#' @param cosNorm A logical scalar indicating whether cosine normalization
#' should be performed on \code{useAssay} prior to PCA. See "See Also". Default
#' \code{TRUE}.
#' @param nComponents An integer scalar specifying the number of dimensions to
#' produce. See "See Also". Default \code{50}.
#' @param weights The weighting scheme to use. Passed to
#' \code{\link[batchelor]{multiBatchPCA}}. Default \code{NULL}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' the SVD should be parallelized.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @seealso \code{\link[batchelor]{fastMNN}} for using \code{useAssay}, and
#' \code{\link[batchelor]{reducedMNN}} for using \code{useReducedDim}
#' @export
#' @references Lun ATL, et al., 2016
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceCorr <- runFastMNN(sceBatches, useAssay = 'logcounts')
runFastMNN <- function(inSCE, useAssay = "logcounts", useReducedDim = NULL,
                       batch = 'batch', reducedDimName = "fastMNN", k = 20,
                       propK = NULL, ndist = 3, minBatchSkip = 0,
                       cosNorm = TRUE, nComponents = 50, weights = NULL,
                       BPPARAM = BiocParallel::SerialParam()){
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                             useReducedDim = useReducedDim,
                             returnMatrix = TRUE)
  mat <- useMat$mat
  batchVec <- .manageCellVar(inSCE, batch, as.factor = TRUE)
  reducedDimName <- gsub(' ', '_', reducedDimName)

  if (!is.null(useReducedDim)) {
    redMNN <- batchelor::reducedMNN(mat, batch = batchVec, BPPARAM = BPPARAM,
                                    k = k, prop.k = propK, ndist = ndist,
                                    min.batch.skip = minBatchSkip)
    newRedDim <- redMNN$corrected
  } else {
    mnnSCE <- batchelor::fastMNN(mat, batch = batchVec, BPPARAM = BPPARAM, k = k,
                                 prop.k = propK, ndist = ndist,
                                 min.batch.skip = minBatchSkip,
                                 cos.norm = cosNorm, d = nComponents,
                                 weights = weights)
    newRedDim <- SingleCellExperiment::reducedDim(mnnSCE, 'corrected')
  }
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- newRedDim
  S4Vectors::metadata(inSCE)$batchCorr[[reducedDimName]] <-
    list(useAssay = useAssay, origLogged = TRUE, method = "fastMNN",
         matType = "reducedDim", batch = batch)
  return(inSCE)
}

#' Apply Harmony batch effect correction method to SingleCellExperiment object
#' @description Harmony is an algorithm that projects cells into a shared
#' embedding in which cells group by cell type rather than dataset-specific
#' conditions.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param useReducedDim A single character indicating the name of the reducedDim
#' used to be corrected. Specifying this will ignore \code{useAssay}. Default
#' \code{NULL}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"HARMONY"}.
#' @param nComponents An integer. The number of PCs to use and generate.
#' Default \code{50L}.
#' @param lambda A Numeric scalar. Ridge regression penalty parameter. Must be
#' strictly positive. Smaller values result in more aggressive correction.
#' Default \code{0.1}.
#' @param theta A Numeric scalar. Diversity clustering penalty parameter. Larger
#' values of theta result in more diverse clusters. theta=0 does not encourage
#' any diversity. Default \code{5}.
#' @param sigma A Numeric scalar. Width of soft kmeans clusters. Larger values
#' of sigma result in cells assigned to more clusters. Smaller values of sigma
#' make soft kmeans cluster approach hard clustering. Default \code{0.1}.
#' @param nIter An integer. The max number of iterations to perform. Default
#' \code{10L}.
#' @param verbose Whether to print progress messages. Default \code{TRUE}.
#' @param ... Other arguments passed to \code{\link[harmony]{HarmonyMatrix}}.
#' See details.
#' @details Since some of the arguments of \code{\link[harmony]{HarmonyMatrix}}
#' is controlled by this wrapper function. The additional arguments users can
#' work with only include: \code{nclust}, \code{tau}, \code{block.size},
#' \code{max.iter.cluster}, \code{epsilon.cluster}, \code{epsilon.harmony},
#' \code{plot_convergence}, \code{reference_values} and \code{cluster_prior}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Ilya Korsunsky, et al., 2019
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' \dontrun{
#' if (require("harmony"))
#'     sceCorr <- runHarmony(sceBatches)
#' }
runHarmony <- function(inSCE, useAssay = "logcounts", useReducedDim = NULL,
                       batch = "batch", reducedDimName = "HARMONY",
                       nComponents = 50, lambda = 0.1, theta = 5,
                       sigma = 0.1, nIter = 10, verbose = TRUE, ...) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("The Harmony package is required to run this function. ",
         "Install harmony with: ",
         "install.packages('harmony')",
         call. = FALSE)
  }
  ## Input check
  useMat <- .selectSCEMatrix(inSCE, useAssay, useReducedDim,
                             useAltExp = NULL, returnMatrix = TRUE)
  batchVec <- .manageCellVar(inSCE, batch, as.factor = TRUE)
  reducedDimName <- gsub(' ', '_', reducedDimName)
  ## Run algorithm
  message(Sys.Date(), " ... Running harmony batch correction")
  mat <- useMat$mat
  do_pca <- ifelse(is.null(useMat$names$useReducedDim), TRUE, FALSE)
  if (!is.null(useMat$names$useReducedDim)) {
    if (nComponents > ncol(mat)) {
      warning("Specified number of components more than available, ",
              "using all (", ncol(mat), ") components.")
      nComponents <- ncol(mat)
    }
    mat <- mat[,seq(nComponents)]
  }
  h <- harmony::HarmonyMatrix(data_mat = mat, meta_data = batchVec,
                              do_pca = do_pca, npcs = nComponents,
                              lambda = lambda, theta = theta,
                              sigma = sigma, max.iter.harmony = nIter,
                              verbose = verbose, return_object = FALSE, ...)
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- h
  S4Vectors::metadata(inSCE)$batchCorr[[reducedDimName]] <-
    list(useAssay = useAssay, origLogged = TRUE, method = "harmony",
         matType = "reducedDim", batch = batch)
  return(inSCE)
}

# #' Apply LIGER batch effect correction method to SingleCellExperiment object
# #'
# #' LIGER relies on integrative non-negative matrix factorization to identify
# #' shared and dataset-specific factors.
# #' @param Input \linkS4class{SingleCellExperiment} object
# #' @param useAssay A single character indicating the name of the assay requiring
# #' batch correction. Default \code{"logcounts"}.
# #' @param batch A single character indicating a field in
# #' \code{\link{colData}} that annotates the batches.
# #' Default \code{"batch"}.
# #' @param reducedDimName A single character. The name for the corrected
# #' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
# #' Default \code{"LIGER"}.
# #' @param nComponents An integer. The number of principle components or
# #' dimensionality to generate in the resulting matrix. Default \code{20L}.
# #' @param lambda A numeric scalar. Algorithmic parameter, the penalty
# #' parameter which limits the dataset-specific component of the factorization.
# #' Default \code{5.0}.
# #' @param resolution A numeric scalar. Algorithmic paramter, the clustering
# #' resolution, increasing this increases the number of communities detected.
# #' Default \code{1.0}
# #' @return The input \linkS4class{SingleCellExperiment} object with
# #' \code{reducedDim(inSCE, reducedDimName)} updated.
# #' @export
# #' @references Joshua Welch, et al., 2018
# #' @examples
# #' \dontrun{
# #' data('sceBatches', package = 'singleCellTK')
# #' sceCorr <- runLIGER(sceBatches)
# #' }
# runLIGER <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
#                      reducedDimName = 'LIGER', nComponents = 20L, lambda = 5.0,
#                      resolution = 1.0){
#   if (!requireNamespace("liger", quietly = TRUE)) {
#     stop("The Liger package is required to run this function. ",
#          "Install liger with: ",
#          "devtools::install_github('joshua-d-campbell/liger') \n",
#          "NOTICE that the one on cran/BiocManager is not what we need.",
#          call. = FALSE)
#   }
#   if (!exists('createLiger', where=asNamespace('liger'), mode='function')) {
#     stop("You are using the wrong source of Liger, please reinstall with: ",
#          "devtools::install_github('joshua-d-campbell/liger') \n",
#          "NOTICE that the one on cran/BiocManager is not what we need.",
#          call. = FALSE)
#   }
#   ## Input check
#   if(!inherits(inSCE, "SingleCellExperiment")){
#     stop("\"inSCE\" should be a SingleCellExperiment Object.")
#   }
#   if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
#     stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
#   }
#   if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
#     stop(paste("\"batch\" name:", batch, "not found"))
#   }
#   reducedDimName <- gsub(' ', '_', reducedDimName)

#   ## Run algorithm
#   batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
#   batches <- unique(batchCol)
#   batchMatrices <- list()
#   for(i in seq_along(batches)){
#     b <- batches[i]
#     batchMatrices[[b]] <-
#       SummarizedExperiment::assay(inSCE, useAssay)[,batchCol == b]
#   }
#   ligerex <- liger::createLiger(batchMatrices)
#   ligerex <- liger::normalize(ligerex)
#   ligerex <- liger::selectGenes(ligerex)
#   ligerex <- liger::scaleNotCenter(ligerex)
#   ligerex <- liger::optimizeALS(ligerex, k = nComponents, lambda = lambda,
#                                 resolution = resolution)
#   ligerex <- liger::quantile_norm(ligerex)
#   SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- ligerex@H.norm
#   return(inSCE)
# }

#' Apply Limma's batch effect correction method to SingleCellExperiment object
#'
#' Limma's batch effect removal function fits a linear model to the data, then
#' removes the component due to the batch effects.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link{assay}}. Default
#' \code{"LIMMA"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Gordon K Smyth, et al., 2003
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceCorr <- runLimmaBC(sceBatches)
runLimmaBC <- function(inSCE, useAssay = "logcounts", assayName = "LIMMA",
                       batch = "batch") {
  ## Input check
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay, returnMatrix = TRUE)
  mat <- useMat$mat
  batchVec <- .manageCellVar(inSCE, batch, as.factor = TRUE)
  assayName <- gsub(' ', '_', assayName)
  ## Run algorithm
  ## One more check for the batch names
  newMat <- limma::removeBatchEffect(mat, batch = batchVec)
  expData(inSCE, assayName, tag = "batchCorrected", altExp = FALSE) <- newMat
  S4Vectors::metadata(inSCE)$batchCorr[[assayName]] <- list(useAssay = useAssay,
                                                            origLogged = TRUE,
                                                            method = "Limma",
                                                            matType = "assay",
                                                            batch = batch)
  return(inSCE)
}

#' Apply the mutual nearest neighbors (MNN) batch effect correction method to
#' SingleCellExperiment object
#'
#' MNN is designed for batch correction of single-cell RNA-seq data where the
#' batches are partially confounded with biological conditions of interest. It
#' does so by identifying pairs of MNN in the high-dimensional log-expression
#' space. For each MNN pair, a pairwise correction vector is computed by
#' applying a Gaussian smoothing kernel with bandwidth `sigma`.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link{assay}}. Default
#' \code{"MNN"}.
#' @param k An integer scalar specifying the number of nearest neighbors to
#' consider when identifying MNNs. See "See Also". Default \code{20}.
#' @param propK A numeric scalar in (0, 1) specifying the proportion of cells in
#' each dataset to use for mutual nearest neighbor searching. See "See Also".
#' Default \code{NULL}.
#' @param sigma A numeric scalar specifying the bandwidth of the Gaussian
#' smoothing kernel used to compute the correction vector for each cell. See
#' "See Also". Default \code{0.1}.
#' @param cosNormIn A logical scalar indicating whether cosine normalization
#' should be performed on the input data prior to calculating distances between
#' cells. See "See Also". Default \code{TRUE}.
#' @param cosNormOut A logical scalar indicating whether cosine normalization
#' should be performed prior to computing corrected expression values. See "See
#' Also". Default \code{TRUE}.
#' @param varAdj A logical scalar indicating whether variance adjustment should
#' be performed on the correction vectors. See "See Also". Default \code{TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' the PCA and nearest-neighbor searches should be parallelized.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @seealso \code{\link[batchelor]{mnnCorrect}}
#' @export
#' @references Haghverdi L, Lun ATL, et. al., 2018
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceCorr <- runMNNCorrect(sceBatches)
runMNNCorrect <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                          assayName = 'MNN', k = 20L, propK = NULL,
                          sigma = 0.1, cosNormIn = TRUE, cosNormOut = TRUE,
                          varAdj = TRUE, BPPARAM = BiocParallel::SerialParam()){
  ## Input check
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay, returnMatrix = TRUE)
  mat <- useMat$mat
  batchVec <- .manageCellVar(inSCE, batch, as.factor = TRUE)
  assayName <- gsub(' ', '_', assayName)
  k <- as.integer(k)

  ## Run algorithm
  corr.sce <- batchelor::mnnCorrect(mat, batch = batchVec, k = k,
                                    prop.k = propK, sigma = sigma,
                                    cos.norm.in = cosNormIn,
                                    cos.norm.out = cosNormOut, var.adj = varAdj,
                                    BPPARAM = BPPARAM)
  expData(inSCE, assayName, tag = "batchCorrected", altExp = FALSE) <-
    SummarizedExperiment::assay(corr.sce, "corrected")
  S4Vectors::metadata(inSCE)$batchCorr[[assayName]] <- list(useAssay = useAssay,
                                                            origLogged = TRUE,
                                                            method = "MNN",
                                                            matType = "assay",
                                                            batch = batch)
  return(inSCE)
}

#' Apply the mutual nearest neighbors (MNN) batch effect correction method to
#' SingleCellExperiment object
#'
#' SCANORAMA is analogous to computer vision algorithms for panorama stitching
#' that identify images with overlapping content and merge these into a larger
#' panorama.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Scanorama requires a transformed normalized expression
#' assay. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in \code{colData} that
#' annotates the batches of each cell; or a vector/factor with the same length
#' as the number of cells. Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link{assay}}. Default
#' \code{"SCANORAMA"}.
#' @param SIGMA A numeric scalar. Algorithmic parameter, correction smoothing
#' parameter on Gaussian kernel. Default \code{15}.
#' @param ALPHA A numeric scalar. Algorithmic parameter, alignment score
#' minimum cutoff. Default \code{0.1}.
#' @param KNN An integer. Algorithmic parameter, number of nearest neighbors to
#' use for matching. Default \code{20}.
#' @param approx Boolean. Use approximate nearest neighbors, greatly speeds up
#' matching runtime. Default \code{TRUE}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Brian Hie et al, 2019
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceCorr <- runSCANORAMA(sceBatches, "ScaterLogNormCounts")
#' }
runSCANORAMA <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                         assayName = 'SCANORAMA', SIGMA = 15, ALPHA = 0.1,
                         KNN = 20, approx = TRUE){
  ## Input check
  if(!reticulate::py_module_available(module = "scanorama")){
    warning("Cannot find python module 'scanorama', please install Conda and",
            " run sctkPythonInstallConda() or run ",
            "sctkPythonInstallVirtualEnv(). If one of these have been ",
            "previously run to install the modules, make sure to run ",
            "selectSCTKConda() or selectSCTKVirtualEnvironment(), ",
            "respectively, if R has been restarted since the module ",
            "installation. Alternatively, scanorama can be installed on the ",
            "local machine with pip (e.g. pip install scanorama) and then the ",
            "'use_python()' function from the 'reticulate' package can be used",
            " to select the correct Python environment.")
    return(inSCE)
  }
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay)
  batchVec <- .manageCellVar(inSCE, batch)
  assayName <- gsub(' ', '_', assayName)
  adata <- .sce2adata(inSCE, useAssay)
  adata$obs[["batch"]] <- batchVec
  py <- reticulate::py
  py$adata <- adata
  py$sigma <- SIGMA
  py$alpha <- ALPHA
  py$KNN <- as.integer(KNN)
  reticulate::py_run_string('
import numpy as np
import scanorama
batches = list(set(adata.obs["batch"]))
adatas = [adata[adata.obs["batch"] == b,] for b in batches]
datasets_full = [a.X for a in adatas]
genes_list = [a.var_names.to_list() for a in adatas]
corrected, genes = scanorama.correct(datasets_full, genes_list, sigma=sigma,
                                     alpha=alpha, knn=KNN)
corrected = [m.toarray() for m in corrected]
cellOrders = [adata.obs["batch"] == b for b in batches]
integrated = np.zeros(adata.shape)
integrated[:,:] = np.NAN
for i in range(len(batches)):
    integrated[cellOrders[i],] = corrected[i]
geneidx = {gene: idx for idx, gene in enumerate(genes)}
orderIdx = []
for gene in adata.var_names:
    orderIdx.append(geneidx[gene])
integrated = integrated[:, orderIdx]
', convert = FALSE)
  mat <- t(py$integrated)
  rownames(mat) <- rownames(inSCE)
  colnames(mat) <- colnames(inSCE)
  expData(inSCE, assayName, tag = "batchCorrected", altExp = FALSE) <- mat
  S4Vectors::metadata(inSCE)$batchCorr[[assayName]] <-
    list(useAssay = useAssay,  origLogged = TRUE, method = "Scanorama",
         matType = "assay", batch = batch)
  return(inSCE)
}

#' Apply scMerge batch effect correction method to SingleCellExperiment object
#'
#' The scMerge method leverages factor analysis, stably expressed genes (SEGs)
#' and (pseudo-) replicates to remove unwanted variations and merge multiple
#' scRNA-Seq data.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link{assay}}. Default \code{"scMerge"}.
#' @param hvgExprs A single characeter. The assay that to be used for highly
#' variable genes identification. Default \code{"counts"}.
#' @param kmeansK An integer vector. Indicating the kmeans' K-value for each
#' batch (i.e. how many subclusters in each batch should exist), in order to
#' construct pseudo-replicates. The length of code{kmeansK} needs to be the same
#' as the number of batches. Default \code{NULL}, and this value will be
#' auto-detected by default, depending on \code{cellType}.
#' @param cellType A single character. A string indicating a field in
#' \code{colData(inSCE)} that defines different cell types. Default
#' \code{'cell_type'}.
#' @param seg A vector of gene names or indices that specifies SEG (Stably
#' Expressed Genes) set as negative control. Pre-defined dataset with human and
#' mouse SEG lists is available with \code{\link[scMerge]{segList}} or
#' \code{\link[scMerge]{segList_ensemblGeneID}}. Default
#' \code{NULL}, and this value will be auto-detected by default with
#' \code{\link[scMerge]{scSEGIndex}}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' should be parallelized. Default \code{BiocParallel::SerialParam()}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Hoa, et al., 2020
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#' logcounts(sceBatches) <- log1p(counts(sceBatches))
#' sceCorr <- runSCMerge(sceBatches)
#' }
runSCMerge <- function(inSCE, useAssay = "logcounts", batch = 'batch',
                       assayName = "scMerge", hvgExprs = "counts", seg = NULL,
                       kmeansK = NULL, cellType = NULL,
                       BPPARAM = BiocParallel::SerialParam()){
  ## Input check
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay, returnMatrix = TRUE)
  .selectSCEMatrix(inSCE, useAssay = hvgExprs, returnMatrix = FALSE)
  mat <- useMat$mat
  batchVec <- .manageCellVar(inSCE, batch)
  cellTypeCol <- .manageCellVar(inSCE, cellType)
  #if(is.null(cellType) & is.null(kmeansK)){
  #  stop("\"cellType\" and \"kmeansK\" cannot be NULL at the same time")
  #}
  assayName <- gsub(' ', '_', assayName)

  ## Run algorithm
  uniqBatch <- unique(batch)

  ## kmeansK
  if(!is.null(cellType) && is.null(kmeansK)){
    # If kmeansK not given, detect by cell type.
    kmeansK <- c()
    for (i in seq_along(uniqBatch)){
      cellTypePerBatch <- cellType[batchVec == uniqBatch[i]]
      kmeansK <- c(kmeansK, length(unique(cellTypePerBatch)))
    }
    cat("Detected kmeansK:\n")
    print(t(data.frame(K = kmeansK, row.names = uniqBatch)))
  }
  ## SEG
  if(is.null(seg)){
    seg <- scMerge::scSEGIndex(mat, cell_type = cellTypeCol, BPPARAM = BPPARAM)
    ctl <- rownames(seg[order(seg$segIdx, decreasing = TRUE)[seq(1000)],])
  } else {
    ctl <- seg
  }

  inSCE <- scMerge::scMerge(sce_combine = inSCE, exprs = useAssay,
                            hvg_exprs = useAssay, batch_name = batch,
                            assay_name = assayName,
                            ctl = ctl, kmeansK = kmeansK,
                            #marker_list = topVarGenesPerBatch,
                            cell_type = cellTypeCol,
                            BPPARAM = BPPARAM)
  # scMerge's function automatically returns the SCE object with information
  # completed, thus using this helper function to simply add the tag.
  inSCE <- expSetDataTag(inSCE, "batchCorrected", assayName)
  S4Vectors::metadata(inSCE)$batchCorr[[assayName]] <-
    list(useAssay = useAssay, origLogged = TRUE, method = "scMerge",
         matType = "assay", batch = batch, condition = cellType)
  return(inSCE)
}

#' Apply ZINBWaVE Batch effect correction method to SingleCellExperiment object
#'
#' A general and flexible zero-inflated negative binomial model that can be
#' used to provide a low-dimensional representations of scRNAseq data. The
#' model accounts for zero inflation (dropouts), over-dispersion, and the count
#' nature of the data. The model also accounts for the difference in library
#' sizes and optionally for batch effects and/or other covariates.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Note that ZINBWaVE works for counts (integer) input rather
#' than logcounts that other methods prefer. Default \code{"counts"}.
#' @param batch A single character indicating a field in
#' \code{\link{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param nHVG An integer. Number of highly variable genes to use when fitting
#' the model. Default \code{1000L}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. Default \code{50L}.
#' @param nIter An integer, The max number of iterations to perform. Default
#' \code{10L}.
#' @param epsilon An integer. Algorithmic parameter. Empirically, a high epsilon
#' is often required to obtained a good low-level representation. Default
#' \code{1000L}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"zinbwave"}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#' should be parallelized. Default \code{BiocParallel::SerialParam()}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Pollen, Alex A et al., 2014
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#'     sceCorr <- runZINBWaVE(sceBatches, nIter = 5)
#' }
runZINBWaVE <- function(inSCE, useAssay = 'counts', batch = 'batch',
                        nHVG = 1000L, nComponents = 50L, epsilon = 1000,
                        nIter = 10L, reducedDimName = 'zinbwave',
                        BPPARAM = BiocParallel::SerialParam()){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch name:", batch, "not found."))
  }
  reducedDimName <- gsub(' ', '_', reducedDimName)
  nHVG <- as.integer(nHVG)
  nComponents <- as.integer(nComponents)
  epsilon <- as.integer(epsilon)
  nIter <- as.integer(nIter)
  # Run algorithm
  ##ZINBWaVE tutorial style of HVG selection
  if(nHVG < nrow(inSCE)){
    logAssay <- log1p(SummarizedExperiment::assay(inSCE, useAssay))
    vars <- matrixStats::rowVars(logAssay)
    names(vars) <- rownames(inSCE)
    vars <- sort(vars, decreasing = TRUE)
    hvgs <- names(vars)[seq_len(nHVG)]
    #tmpSCE <- inSCE[names(vars)[seq_len(nHVG)],]
  } else {
    #tmpSCE <- inSCE
    hvgs <- rownames(inSCE)
  }
  epsilon <- min(nrow(inSCE), epsilon)
  inSCE <- zinbwave::zinbwave(inSCE, K = nComponents, epsilon = epsilon,
                               which_assay = useAssay, which_genes = hvgs,
                               X = paste('~', batch, sep = ''),
                               maxiter.optimize = nIter, verbose = TRUE,
                              BPPARAM = BPPARAM)
  reddimName <- reducedDimNames(inSCE)
  reddimName[reddimName == "zinbwave"] <- reducedDimName
  #SingleCellExperiment::reducedDim(inSCE, reducedDimName) <-
  #  SingleCellExperiment::reducedDim(inSCE, 'zinbwave')
  S4Vectors::metadata(inSCE)$batchCorr[[reducedDimName]] <-
    list(useAssay = useAssay, origLogged = FALSE, method = "ZINBWaVE",
         matType = "reducedDim", batch = batch)
  return(inSCE)
}
