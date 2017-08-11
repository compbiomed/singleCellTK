#' Get PCA components for a SCE object
#'
#' Selects the 500 most variable genes in the SCE, performs
#' PCA based on them and outputs the principal components in a data frame
#' and their variances as percentVar attribute
#'
#' @param count_data SCE set
#'
#' @return A reduced dimension object
#' @export getPCA
#'
getPCA <- function(count_data){
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
  exprs_to_plot <- exprs_mat[feature_set, , drop = FALSE]
  exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
  keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_to_plot <- exprs_to_plot[, keep_feature]
  pca <- prcomp(exprs_to_plot)
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  pca <- data.frame(pca$x)
  attr(pca, "percentVar") <- percentVar
  return(pca)
}
