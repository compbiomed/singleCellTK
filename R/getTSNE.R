#' Get t-SNE components for a SCE object
#'
#' Selects the 500 most variable genes in the feature count, performs
#' t-SNE based on them and outputs the principal components in a data frame
#' and their variances as percentVar attribute
#'
#' @param count_data SCE object
#'
#' @return A reduced dimension object
#' @export getTSNE
#'
getTSNE <- function(count_data){
  if (nrow(count_data) < 500){
    ntop <- nrow(count_data)
  }
  else{
    ntop <- 500
  }
  scale_features <- TRUE
  exprs_mat <- log2(assay(count_data, "counts") + 1)
  rv <- matrixStats::rowVars(exprs_mat)
  feature_set <-
    order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprs_to_plot <- exprs_mat[feature_set, ]
  keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_to_plot <- exprs_to_plot[keep_feature, ]
  exprs_to_plot <- t(scale(t(exprs_to_plot), scale = scale_features))
  perplexity <- floor(ncol(count_data) / 5)
  tsne_out <- Rtsne::Rtsne(t(exprs_to_plot), perplexity = perplexity,
                           initial_dims = max(50, ncol(count_data)))
  tsne_out <- data.frame(tsne_out$Y[, 1:2],
                           row.names = colnames(count_data))
  return(tsne_out)
}
