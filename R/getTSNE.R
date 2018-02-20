#' Get t-SNE components for a SCE object
#'
#' Selects the 500 most variable genes in the feature count, performs
#' t-SNE based on them and stores the t-SNE values in the reducedDims slot of the
#' SCE object.
#'
#' @param count_data SCE object
#' @param use_assay Indicate which assay to use for PCA. Default is "counts"
#' @param reducedDimName Store the t-SNE data with this name. The default is
#' TSNE. The toolkit will store data with the pattern <ASSAY>_<ALGORITHM>.
#'
#' @return A SCE object with the specified reducedDim updated
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' #add a CPM assay
#' assay(GSE60361_subset_sce, "cpm") <- apply(assay(GSE60361_subset_sce,
#'                                                  "counts"),
#'                                            2, function(x) {
#'                                              x / (sum(x) / 1000000)
#'                                            })
#' GSE60361_subset_sce <- getTSNE(GSE60361_subset_sce,
#'                                use_assay = "cpm",
#'                                reducedDimName = "TSNE_cpm")
#' reducedDims(GSE60361_subset_sce)
#'
getTSNE <- function(count_data, use_assay="logcounts", reducedDimName="TSNE"){
  if (nrow(count_data) < 500){
    ntop <- nrow(count_data)
  } else{
    ntop <- 500
  }
  if (!(use_assay %in% names(assays(count_data)))){
    stop(use_assay, " not in the assay list")
  }
  exprs_mat <- log2(assay(count_data, use_assay) + 1)
  rv <- matrixStats::rowVars(exprs_mat)
  feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprs_to_plot <- exprs_mat[feature_set, ]
  keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_to_plot <- exprs_to_plot[keep_feature, ]
  exprs_to_plot <- t(scale(t(exprs_to_plot)))
  perplexity <- floor(ncol(count_data) / 5)
  tsne_out <- Rtsne::Rtsne(t(exprs_to_plot), perplexity = perplexity,
                           initial_dims = max(50, ncol(count_data)))
  tsne_out <- tsne_out$Y[, 1:2]
  rownames(tsne_out) <- colnames(count_data)
  reducedDim(count_data, reducedDimName) <- tsne_out
  return(count_data)
}
