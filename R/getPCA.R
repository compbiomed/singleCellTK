#' Get PCA components for a feature counts table
#'   
#' Selects the 500 most variable genes in the feature count, performs
#' PCA based on them and outputs the principal components in a data frame
#' and their variances as percentVar attribute
#'
#' @param count_data Matrix of expresion values
#' 
#' @return A reduced dimension object
#' @export getPCA
#'
getPCA <- function(count_data){
  if (nrow(counts(count_data)) < 500){
    ntop <- nrow(counts(count_data))
  }
  else{
    ntop <- 500
  }
  scale_features <- TRUE
  exprs_mat <- exprs(count_data)
  rv <- matrixStats::rowVars(exprs_mat)
  feature_set <-
    order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  exprs_to_plot <- exprs_mat[feature_set,,drop=FALSE]
  exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
  keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_to_plot <- exprs_to_plot[, keep_feature]
  pca <- prcomp(exprs_to_plot)
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  pca <- data.frame(pca$x)
  attr(pca,"percentVar") <- percentVar
  return(pca)
}
