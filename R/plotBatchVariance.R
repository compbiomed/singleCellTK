#' Plot the percent of the variation that is explained by batch and condition
#' in the data
#'
#' Visualize the percent variation in the data that is explained by batch and
#' condition if it is given.
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param use_assay Indicate which assay to use for PCA. Default is "logcounts"
#' @param batch The column in the annotation data that corresponds to batch.
#' Required
#' @param condition The column in the annotation data that corresponds to
#' condition. Optional
#'
#' @return A boxplot of variation explained by batch, condition, and
#' batch+condition (if applicable).
#' @export
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'   dat <- as(as(bladderEset, "SummarizedExperiment"), "SCtkExperiment")
#'   plotBatchVariance(dat, use_assay="exprs", batch="batch", condition = "cancer")
#' }
#'
plotBatchVariance <- function(inSCESet, use_assay="logcounts", batch,
                              condition=NULL){
  nlb <- nlevels(as.factor(SingleCellExperiment::colData(inSCESet)[, batch]))
  if(nlb <= 1){
    batch_mod <- matrix(rep(1, ncol(inSCESet)), ncol = 1)
  } else {
    batch_mod <- stats::model.matrix(~as.factor(SingleCellExperiment::colData(inSCESet)[, batch]))
  }
  if(is.null(condition)){
    stop("condition required for now")
  } else {
    nlc <- nlevels(as.factor(SingleCellExperiment::colData(inSCESet)[, condition]))
    if(nlc <= 1){
      cond_mod <- matrix(rep(1, ncol(inSCESet)), ncol = 1)
    } else {
      cond_mod <- stats::model.matrix(~as.factor(SingleCellExperiment::colData(inSCESet)[, condition]))
    }
  }

  mod <- cbind(cond_mod, batch_mod[, -1])

  cond_test <- batchqc_f.pvalue(SummarizedExperiment::assay(inSCESet, use_assay), mod, batch_mod)
  batch_test <- batchqc_f.pvalue(SummarizedExperiment::assay(inSCESet, use_assay), mod, cond_mod)

  cond_ps <- cond_test$p
  batch_ps <- batch_test$p

  r2_full <- cond_test$r2_full
  cond_r2 <- batch_test$r2_reduced
  batch_r2 <- cond_test$r2_reduced
  explained_variation <- round(cbind(`Full (Condition+Batch)` = r2_full,
                                     Condition = cond_r2, Batch = batch_r2), 5) * 100
  ex_var_m <- reshape2::melt(explained_variation)
  colnames(ex_var_m) <- c("Gene", "Model", "Value")
  a <- ggplot2::ggplot(ex_var_m, ggplot2::aes(factor(Model), Value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = factor(Model))) +
    ggplot2::geom_boxplot(width = .1) +
    ggplot2::xlab("Model") +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1"), guide = FALSE)
  return(a)
}

batchqc_f.pvalue <- function(dat, mod, mod0) {
  ## F-test (full/reduced model) and returns R2 values
  ## (full/reduced) as well.
  mod00 <- matrix(rep(1, ncol(dat)), ncol = 1)
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)

  resid <- dat - dat %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  rss1 <- rowSums(resid * resid)
  rm(resid)

  resid0 <- dat - dat %*% mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)

  resid00 <- dat - dat %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
  rss00 <- rowSums(resid00 * resid00)
  rm(resid00)

  r2_full <- 1 - rss1 / rss00
  r2_reduced <- 1 - rss0 / rss00

  p <- 1
  if (df1 > df0)  {
    fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
    p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  }
  return(list(p = p, r2_full = r2_full, r2_reduced = r2_reduced))
}
