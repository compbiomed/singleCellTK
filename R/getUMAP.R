
#' Uniform Manifold Approximation and Projection(UMAP) algorithm for
#' dimension reduction.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. The default is "logcounts".
#' @param reducedDimName a name to store the results of the dimension reduction
#' coordinates obtained from this method. This is stored in the SingleCellExperiment
#' object in the reducedDims slot. Required.
#' @param n_neighbors specify the number of nearest neighbors. Default is 5.
#' @param n_iterations number of iterations performed during layout optimization.
#' Default is 200.
#' @param alpha initial value of "learning rate" of layout optimization. Default is 1.
#' @param metric distance metric. Default is euclidean, other options are 'manhattan', 'cosine', 'pearson'.
#' @param run_pca run tSNE on PCA components? Default is TRUE.
#'
#' @return a SCtkExperiment object with the reduced dimensions updated under
#' reducedDimName specified.
#' @export
#'
#' @examples
#' umap_res <- getUMAP(inSCE = mouseBrainSubsetSCE, useAssay = "counts",
#'                     reducedDimName = "UMAP", n_neighbors = 3, n_iterations = 200,
#'                     alpha = 1, metric = "euclidean",  run_pca = TRUE)
#' reducedDims(umap_res)
#'
getUMAP <- function(inSCE, useAssay = "logcounts", reducedDimName = "UMAP",
                    n_neighbors = 5, n_iterations = 200, alpha = 1, metric="euclidean",
                    run_pca = TRUE) {
  if (!(class(inSCE) %in% c("SingleCellExperiment", "SCtkExperiment", "SummarizedExperiment"))){
    stop("Please use a SingleCellExperiment or a SCtkExperiment object")
  }
  #test for assay existing
  if (!all(useAssay %in% names(assays(inSCE)))){
    stop("assay '", useAssay, "' does not exist.")
  }
  matColData <- SummarizedExperiment::assay(inSCE, useAssay)
  data <- as.data.frame(t(matColData))
  if(run_pca == TRUE) {
    
    getpca <- getPCA(inSCE = inSCE, useAssay = useAssay, reducedDimName = "PCA_UMAP")
    retpca <- reducedDim(getpca, "PCA_UMAP")
    pc_mat <- as.matrix(t(retpca))[, 1:2]
   
    #running umap with p1
    umap_results <- uwot::umap(data,
                               init = pc_mat,
                               n_neighbors = n_neighbors,
                               metric = metric,
                               learning_rate = alpha,
                               n_epochs = n_iterations
                               )
  } else {
    umap_results <- uwot::umap(data,
                               n_neighbors = n_neighbors,
                               metric = metric,
                               learning_rate = alpha,
                               n_epochs = n_iterations)
  }
  if (is.null(rownames(inSCE))) {
    rownames(umap_results) <- colnames(inSCE)
  }
  colnames(umap_results) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- umap_results
  return(inSCE)
}
