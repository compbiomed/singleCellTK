#' Get PCA components for a SCtkE object
#'
#' Selects the 500 most variable genes in the SCE, performs
#' PCA based on them and stores the PCA values in the reducedDims slot of the
#' SCE object.
#'
#' @param count_data SCtkE object
#' @param use_assay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the PCA data with this name. The default is PCA.
#' The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCtkE object with the specified reduecedDim and pca_variances
#' updated
#' @export
#' @examples
#' data("mouse_brain_subset_sce")
#' #add a CPM assay
#' assay(mouse_brain_subset_sce, "cpm") <- apply(assay(mouse_brain_subset_sce,
#'                                                  "counts"),
#'                                            2, function(x) {
#'                                              x / (sum(x) / 1000000)
#'                                            })
#' mouse_brain_subset_sce <- getPCA(mouse_brain_subset_sce,
#'                               use_assay = "cpm",
#'                               reducedDimName = "PCA_cpm")
#' reducedDims(mouse_brain_subset_sce)
#'
getPCA <- function(count_data, use_assay="logcounts", reducedDimName="PCA"){
  if (nrow(count_data) < 500){
    ntop <- nrow(count_data)
  }
  else{
    ntop <- 500
  }
  if (!(use_assay %in% names(SummarizedExperiment::assays(count_data)))){
    stop(use_assay, " not in the assay list")
  }
  exprs_mat <- log2(SummarizedExperiment::assay(count_data, use_assay) + 1)
  rv <- matrixStats::rowVars(exprs_mat)
  feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprs_to_plot <- exprs_mat[feature_set, , drop = FALSE]
  exprs_to_plot <- scale(t(exprs_to_plot))
  keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_to_plot <- exprs_to_plot[, keep_feature]
  pca <- stats::prcomp(exprs_to_plot)
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  pca <- pca$x
  SingleCellExperiment::reducedDim(count_data, reducedDimName) <- pca
  if (class(count_data) == "SCtkExperiment"){
    pca_variances(count_data) <- S4Vectors::DataFrame(percentVar)
    rownames(pca_variances(count_data)) <- paste0(
      "PC", 1:nrow(pca_variances(count_data)))
  }
  return(count_data)
}
