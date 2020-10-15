#' Apply BBKNN batch effect correction method to SingleCellExperiment object
#'
#' BBKNN, an extremely fast graph-based data integration algorithm. It modifies
#' the neighbourhood construction step to produce a graph that is balanced
#' across all batches of the data.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
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
#' sceCorr <- runBBKNN(sceBatches)
#' }
runBBKNN <-function(inSCE, useAssay = 'logcounts', batch = 'batch',
                    reducedDimName = 'BBKNN', nComponents = 50L){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
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
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found"))
  }
  reducedDimName <- gsub(' ', '_', reducedDimName)
  nComponents <- as.integer(nComponents)
  ## Run algorithm
  adata <- .sce2adata(inSCE, useAssay = useAssay)
  sc$tl$pca(adata, n_comps = nComponents)
  bbknn$bbknn(adata, batch_key = batch, n_pcs = nComponents)
  sc$tl$umap(adata, n_components = nComponents)
  bbknnUmap <- adata$obsm[["X_umap"]]
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- bbknnUmap
  return(inSCE)
}

#' Apply ComBat batch effect correction method to SingleCellExperiment object
#'
#' The ComBat batch adjustment approach assumes that batch effects represent
#' non-biological but systematic shifts in the mean or variability of genomic
#' features for all samples within a processing batch. It uses either parametric
#' or non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param par.prior A logical scalar. TRUE indicates parametric adjustments
#' will be used, FALSE indicates non-parametric adjustments will be used.
#' Default \code{TRUE}.
#' @param covariates List of other column names in colData to be added to the
#' ComBat model as covariates. Default \code{NULL}.
#' @param mean.only If TRUE ComBat only corrects the mean of the batch effect.
#' Default \code{FALSE}.
#' @param ref.batch If given, will use the selected batch as a reference for
#' batch adjustment. Default \code{NULL}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"ComBat"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' # parametric adjustment
#' sceCorr <- runComBat(sceBatches)
#' # non-parametric adjustment, mean-only version
#' sceCorr <- runComBat(sceBatches, par.prior=FALSE, mean.only=TRUE)
#' # reference-batch version, with covariates
#' sceCorr <- runComBat(sceBatches, covariates = "cell_type", ref.batch = 'w')
#' }
#' @export
runComBat <- function(inSCE, useAssay = "logcounts", batch = 'batch',
                      par.prior = TRUE, covariates = NULL,
                      mean.only = FALSE, ref.batch = NULL,
                      assayName = "ComBat") {
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(any(!c(batch, covariates) %in%
         names(SummarizedExperiment::colData(inSCE)))){
    anns <- c(batch, covariates)
    notFound <- which(!anns %in% names(SummarizedExperiment::colData(inSCE)))
    notFound <- anns[notFound]
    stop("\"annotation\" name:", paste(notFound, collapse = ', '), "not found")
  }
  #prepare model matrix
  mod <- NULL
  if (!is.null(covariates)){
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SummarizedExperiment::colData(inSCE)))
  }

  resassay <-
    sva::ComBat(dat = SummarizedExperiment::assay(inSCE, useAssay),
                batch = SummarizedExperiment::colData(inSCE)[[batch]],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)

  SummarizedExperiment::assay(inSCE, assayName) <- resassay
  return(inSCE)
}

#' Apply a fast version of the mutual nearest neighbors (MNN) batch effect
#' correction method to SingleCellExperiment object
#'
#' fastMNN is a variant of the classic MNN method, modified for speed and more
#' robust performance. For introduction of MNN, see \code{\link{runMNNCorrect}}.
#' @param inSCE  inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}. Alternatively, see
#' \code{pcInput} parameter.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"fastMNN"}.
#' @param pcInput A logical scalar. Whether to use a low-dimension matrix for
#' batch effect correction. If \code{TRUE}, \code{useAssay} will be searched
#' from \code{reducedDimNames(inSCE)}. Default \code{FALSE}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Lun ATL, et al., 2016
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runFastMNN(sceBatches, useAssay = 'logcounts', pcInput = FALSE)
#' }
runFastMNN <- function(inSCE, useAssay = "logcounts",
                       reducedDimName = "fastMNN", batch = 'batch',
                       pcInput = FALSE){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(isTRUE(pcInput)){
    if(!(useAssay %in% SingleCellExperiment::reducedDimNames(inSCE))) {
      stop(paste("\"useAssay\" (reducedDim) name: ", useAssay, " not found."))
    }
  } else {
    if(!(useAssay %in% SummarizedExperiment::assayNames(inSCE))) {
      stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
    }
  }

  if(!(batch %in% names(SummarizedExperiment::colData(inSCE)))){
    stop(paste("\"batch name:", batch, "not found."))
  }
  reducedDimName <- gsub(' ', '_', reducedDimName)

  ## Run algorithm
  batches <- list()
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  batchFactor <- as.factor(batchCol)

  if(pcInput){
    mat <- SingleCellExperiment::reducedDim(inSCE, useAssay)
    redMNN <- batchelor::reducedMNN(mat, batch = batchFactor)
    newRedDim <- redMNN$corrected
  } else {
    mat <- SummarizedExperiment::assay(inSCE, useAssay)
    mnnSCE <- batchelor::fastMNN(mat, batch = batchFactor)
    newRedDim <- SingleCellExperiment::reducedDim(mnnSCE, 'corrected')
  }
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- newRedDim
  return(inSCE)
}

#' Apply Harmony batch effect correction method to SingleCellExperiment object
#'
#' Harmony is an algorithm that projects cells into a shared embedding in which
#' cells group by cell type rather than dataset-specific conditions.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"HARMONY"}.
#' @param pcInput A logical scalar. Whether to use a low-dimension matrix for
#' batch effect correction. If \code{TRUE}, \code{useAssay} will be searched
#' from \code{reducedDimNames(inSCE)}. Default \code{FALSE}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. If \code{pcInput} is set
#' to \code{TRUE}, the output dimension will follow the low-dimension matrix,
#' so this argument will be ignored. Default \code{50L}.
#' @param nIter An integer. The max number of iterations to perform. Default
#' \code{10L}.
#' @param theta A Numeric scalar. Diversity clustering penalty parameter,
#' Larger value results in more diverse clusters. Default \code{5}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Ilya Korsunsky, et al., 2019
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runHarmony(sceBatches, nComponents = 10L)
runHarmony <- function(inSCE, useAssay = "logcounts", pcInput = FALSE,
                       batch = "batch", reducedDimName = "HARMONY",
                       nComponents = 50L, theta = 5, nIter = 10L){
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("The Harmony package is required to run this function. ",
         "Install harmony with: ",
         "devtools::install_github('joshua-d-campbell/harmony') ",
         call. = FALSE)
  }
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(pcInput){
    if(!useAssay %in% SingleCellExperiment::reducedDimNames(inSCE)) {
      stop(paste("\"useAssay\" (reducedDim) name: ", useAssay, " not found."))
    }
  } else {
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
      stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
    }
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found"))
  }
  reducedDimName <- gsub(' ', '_', reducedDimName)
  nComponents <- as.integer(nComponents)
  nIter <- as.integer(nIter)
  ## Run algorithm
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  if(pcInput){
    mat <- SingleCellExperiment::reducedDim(inSCE, useAssay)
  } else{
    sceTmp <- scater::runPCA(inSCE, exprs_values = useAssay,
                             ncomponents = nComponents)
    mat <- SingleCellExperiment::reducedDim(sceTmp, 'PCA')
  }
  h <- harmony::HarmonyMatrix(mat, batchCol, do_pca = FALSE,
                              theta = theta, max.iter.harmony = nIter)
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- h
  return(inSCE)
}

#' Apply LIGER batch effect correction method to SingleCellExperiment object
#'
#' LIGER relies on integrative non-negative matrix factorization to identify
#' shared and dataset-specific factors.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param reducedDimName A single character. The name for the corrected
#' low-dimensional representation. Will be saved to \code{reducedDim(inSCE)}.
#' Default \code{"LIGER"}.
#' @param nComponents An integer. The number of principle components or
#' dimensionality to generate in the resulting matrix. Default \code{20L}.
#' @param lambda A numeric scalar. Algorithmic parameter, the penalty
#' parameter which limits the dataset-specific component of the factorization.
#' Default \code{5.0}.
#' @param resolution A numeric scalar. Algorithmic paramter, the clustering
#' resolution, increasing this increases the number of communities detected.
#' Default \code{1.0}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Joshua Welch, et al., 2018
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runLIGER(sceBatches)
#' }
runLIGER <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                     reducedDimName = 'LIGER', nComponents = 20L, lambda = 5.0,
                     resolution = 1.0){
  if (!requireNamespace("liger", quietly = TRUE)) {
    stop("The Liger package is required to run this function. ",
         "Install liger with: ",
         "devtools::install_github('joshua-d-campbell/liger') \n",
         "NOTICE that the one on cran/BiocManager is not what we need.",
         call. = FALSE)
  }
  if (!exists('createLiger', where=asNamespace('liger'), mode='function')) {
    stop("You are using the wrong source of Liger, please reinstall with: ",
         "devtools::install_github('joshua-d-campbell/liger') \n",
         "NOTICE that the one on cran/BiocManager is not what we need.",
         call. = FALSE)
  }
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found"))
  }
  reducedDimName <- gsub(' ', '_', reducedDimName)

  ## Run algorithm
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  batches <- unique(batchCol)
  batchMatrices <- list()
  for(i in seq_along(batches)){
    b <- batches[i]
    batchMatrices[[b]] <-
      SummarizedExperiment::assay(inSCE, useAssay)[,batchCol == b]
  }
  ligerex <- liger::createLiger(batchMatrices)
  ligerex <- liger::normalize(ligerex)
  ligerex <- liger::selectGenes(ligerex)
  ligerex <- liger::scaleNotCenter(ligerex)
  ligerex <- liger::optimizeALS(ligerex, k = nComponents, lambda = lambda,
                                resolution = resolution)
  ligerex <- liger::quantile_norm(ligerex)
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- ligerex@H.norm
  return(inSCE)
}

#' Apply Limma's batch effect correction method to SingleCellExperiment object
#'
#' Limma's batch effect removal function fits a linear model to the data, then
#' removes the component due to the batch effects.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"LIMMA"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Gordon K Smyth, et al., 2003
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runLimmaBC(sceBatches)
runLimmaBC <- function(inSCE, useAssay = "logcounts", assayName = "LIMMA",
                       batch = "batch") {
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found."))
  }

  assayName <- gsub(' ', '_', assayName)

  ## Run algorithm
  ## One more check for the batch names
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  mat <- SummarizedExperiment::assay(inSCE, useAssay)
  newMat <- limma::removeBatchEffect(mat, batch = batchCol)
  SummarizedExperiment::assay(inSCE, assayName) <- newMat
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
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param k An integer. Specifies the number of nearest neighbours to
#' consider when defining MNN pairs. This should be interpreted as the minimum
#' frequency of each cell type or state in each batch. Larger values will
#' improve the precision of the correction by increasing the number of MNN
#' pairs, at the cost of reducing accuracy by allowing MNN pairs to form between
#' cells of different type. Default \code{20L}.
#' @param sigma A Numeric scalar. Specifies how much information is
#' shared between MNN pairs when computing the batch effect. Larger values will
#' share more information, approaching a global correction for all cells in the
#' same batch. Smaller values allow the correction to vary across cell types,
#' which may be more accurate but comes at the cost of precision. Default
#' \code{0.1}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"MNN"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Lun ATL, et al., 2016 & 2018
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runMNNCorrect(sceBatches)
runMNNCorrect <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                          assayName = 'MNN', k = 20L, sigma = 0.1){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch name:", batch, "not found."))
  }
  assayName <- gsub(' ', '_', assayName)
  k <- as.integer(k)

  ## Run algorithm
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  batchFactor <- as.factor(batchCol)
  mnnSCE <- batchelor::mnnCorrect(inSCE, assay.type = useAssay,
                                  batch = batchFactor, k = k, sigma = sigma)
  corrected <- SummarizedExperiment::assay(mnnSCE, 'corrected')
  SummarizedExperiment::assay(inSCE, assayName) <- corrected
  return(inSCE)
}

#' Apply the mutual nearest neighbors (MNN) batch effect correction method to
#' SingleCellExperiment object
#'
#' SCANORAMA is analogous to computer vision algorithms for panorama stitching
#' that identify images with overlapping content and merge these into a larger
#' panorama.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param SIGMA A numeric scalar. Algorithmic parameter, correction smoothing
#' parameter on Gaussian kernel. Default \code{15}.
#' @param ALPHA A numeric scalar. Algorithmic parameter, alignment score
#' minimum cutoff. Default \code{0.1}.
#' @param KNN An integer. Algorithmic parameter, number of nearest neighbors to
#' use for matching. Default \code{20L}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"SCANORAMA"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Brian Hie et al, 2019
#' @examples
#' \dontrun{
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCANORAMA(sceBatches)
#' }
runSCANORAMA <- function(inSCE, useAssay = 'logcounts', batch = 'batch',
                         SIGMA = 15, ALPHA = 0.1, KNN = 20L,
                         assayName = 'SCANORAMA'){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
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
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found"))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found"))
  }
  assayName <- gsub(' ', '_', assayName)
  adata <- .sce2adata(inSCE, useAssay)
  py <- reticulate::py
  py$adata <- adata
  py$batch <- batch
  py$sigma <- SIGMA
  py$alpha <- ALPHA
  py$KNN <- as.integer(KNN)
  reticulate::py_run_string('
import numpy as np
import scanorama
batches = list(set(adata.obs[batch]))
adatas = [adata[adata.obs[batch] == b,] for b in batches]
datasets_full = [a.X.toarray() for a in adatas]
genes_list = [a.var_names.to_list() for a in adatas]
corrected, genes = scanorama.correct(datasets_full, genes_list, sigma=sigma,
                                     alpha=alpha, knn=KNN)
corrected = [m.toarray() for m in corrected]
cellOrders = [adata.obs[batch] == b for b in batches]
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
  SummarizedExperiment::assay(inSCE, assayName) <- mat
  return(inSCE)
}

#' Apply scMerge batch effect correction method to SingleCellExperiment object
#'
#' The scMerge method leverages factor analysis, stably expressed genes (SEGs)
#' and (pseudo-) replicates to remove unwanted variations and merge multiple
#' scRNA-Seq data.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
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
#' mouse SEG lists is available to user by running \code{data('SEG')}. Default
#' \code{NULL}, and this value will be auto-detected by default with
#' \code{\link[scMerge]{scSEGIndex}}.
#' @param nCores An integer. The number of cores of processors to allocate for
#' the task. Default \code{1L}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"scMerge"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Hoa, et al., 2020
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCMerge(sceBatches)
runSCMerge <- function(inSCE, useAssay = "logcounts", batch = 'batch',
                       assayName = "scMerge", seg = NULL, kmeansK = NULL,
                       cellType = 'cell_type',
                       nCores = 1L){
  ## Input check
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
    stop(paste("\"batch\" name:", batch, "not found"))
  }
  if(is.null(cellType) & is.null(kmeansK)){
    stop("\"cellType\" and \"kmeansK\" cannot be NULL at the same time")
  }
  if(!cellType %in% names(SummarizedExperiment::colData(inSCE))){
    # If NULL, scMerge still works
    stop(paste("\"cellType\" name:", cellType, "not found"))
  }

  nCores <- min(as.integer(nCores), parallel::detectCores())
  assayName <- gsub(' ', '_', assayName)

  ## Run algorithm

  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  uniqBatch <- unique(batchCol)

  # Infer parameters
  if(is.null(cellType)){
    cellTypeCol <- NULL
  } else {
    cellTypeCol <- SummarizedExperiment::colData(inSCE)[[cellType]]
  }
  ## kmeansK
  if(!is.null(cellType) && is.null(kmeansK)){
    # If kmeansK not given, detect by cell type.
    cellTypeCol <- SummarizedExperiment::colData(inSCE)[[cellType]]
    kmeansK <- c()
    for (i in seq_along(uniqBatch)){
      cellTypePerBatch <- cellTypeCol[batchCol == uniqBatch[i]]
      kmeansK <- c(kmeansK, length(unique(cellTypePerBatch)))
    }
    cat("Detected kmeansK:\n")
    print(t(data.frame(K = kmeansK, row.names = uniqBatch)))
  }
  ## SEG
  if(is.null(seg)){
    bpParam <- BiocParallel::MulticoreParam(workers = nCores)
    seg <- scMerge::scSEGIndex(SummarizedExperiment::assay(inSCE, useAssay),
                               cell_type = cellTypeCol,
                               BPPARAM = bpParam)
    ctl <- rownames(seg[order(seg$segIdx, decreasing = TRUE)[seq_len(1000)],])
  } else {
    ctl <- seg
  }

  # scMerge automatically search for the column called "batch"...
  colDataNames <- names(SummarizedExperiment::colData(inSCE))
  names(SummarizedExperiment::colData(inSCE))[colDataNames == batch] <- 'batch'
  bpParam <- BiocParallel::MulticoreParam(workers = nCores)
  inSCE <- scMerge::scMerge(sce_combine = inSCE, exprs = useAssay,
                            hvg_exprs = useAssay,
                            assay_name = assayName,
                            ctl = ctl, kmeansK = kmeansK,
                            #marker_list = topVarGenesPerBatch,
                            cell_type = cellTypeCol,
                            BPPARAM = bpParam)
  colDataNames <- names(SummarizedExperiment::colData(inSCE))
  names(SummarizedExperiment::colData(inSCE))[colDataNames == 'batch'] <- batch
  return(inSCE)
}

#' Apply ZINBWaVE Batch effect correction method to SingleCellExperiment object
#'
#' A general and flexible zero-inflated negative binomial model that can be
#' used to provide a low-dimensional representations of scRNAseq data. The
#' model accounts for zero inflation (dropouts), over-dispersion, and the count
#' nature of the data. The model also accounts for the difference in library
#' sizes and optionally for batch effects and/or other covariates.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Note that ZINBWaVE works for counts (integer) input rather
#' than logcounts that other methods prefer. Default \code{"counts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
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
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim(inSCE, reducedDimName)} updated.
#' @export
#' @references Pollen, Alex A et al., 2014
#' @examples
#' \dontrun{
#'     data('sceBatches', package = 'singleCellTK')
#'     sceCorr <- runZINBWaVE(sceBatches, nIter = 5)
#' }
runZINBWaVE <- function(inSCE, useAssay = 'counts', batch = 'batch',
                        nHVG = 1000L, nComponents = 50L, epsilon = 1000,
                        nIter = 10L, reducedDimName = 'zinbwave'){
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
    tmpSCE <- inSCE[names(vars)[seq_len(nHVG)],]
  }
  epsilon <- min(nrow(inSCE), epsilon)
  tmpSCE <- zinbwave::zinbwave(tmpSCE, K = nComponents, epsilon = epsilon,
                               which_assay = useAssay,
                               X = paste('~', batch, sep = ''),
                               maxiter.optimize = nIter, verbose = TRUE)
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <-
    SingleCellExperiment::reducedDim(tmpSCE, 'zinbwave')
  return(inSCE)
}
