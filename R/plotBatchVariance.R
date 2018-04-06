#' Plot the percent of the variation that is explained by batch and condition
#' in the data
#'
#' Visualize the percent variation in the data that is explained by batch and
#' condition if it is given.
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"
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
#'   plotBatchVariance(dat, useAssay="exprs", batch="batch", condition = "cancer")
#' }
#'
plotBatchVariance <- function(inSCESet, useAssay="logcounts", batch,
                              condition=NULL){
  nlb <- nlevels(as.factor(SingleCellExperiment::colData(inSCESet)[, batch]))
  if (nlb <= 1){
    batchMod <- matrix(rep(1, ncol(inSCESet)), ncol = 1)
  } else {
    batchMod <- stats::model.matrix(
      ~as.factor(SingleCellExperiment::colData(inSCESet)[, batch]))
  }
  if (is.null(condition)){
    stop("condition required for now")
  } else {
    nlc <- nlevels(as.factor(
      SingleCellExperiment::colData(inSCESet)[, condition]))
    if (nlc <= 1){
      condMod <- matrix(rep(1, ncol(inSCESet)), ncol = 1)
    } else {
      condMod <- stats::model.matrix(
        ~as.factor(SingleCellExperiment::colData(inSCESet)[, condition]))
    }
  }

  mod <- cbind(condMod, batchMod[, -1])

  condTest <- batchqc_f.pvalue(SummarizedExperiment::assay(inSCESet, useAssay),
                                mod, batchMod)
  batchTest <- batchqc_f.pvalue(
    SummarizedExperiment::assay(inSCESet, useAssay), mod, condMod)

  r2Full <- condTest$r2Full
  condR2 <- batchTest$r2Reduced
  batchR2 <- condTest$r2Reduced
  explainedVariation <- round(cbind(`Full (Condition+Batch)` = r2Full,
                                     Condition = condR2,
                                     Batch = batchR2), 5) * 100
  exVarM <- reshape2::melt(explainedVariation)
  colnames(exVarM) <- c("Gene", "Model", "Value")
  exVarM$Model <- factor(exVarM$Model)
  a <- ggplot2::ggplot(exVarM, ggplot2::aes_string("Model", "Value")) +
    ggplot2::geom_violin(ggplot2::aes_string(fill = "Model")) +
    ggplot2::geom_boxplot(width = .1) +
    ggplot2::xlab("Model") +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1"),
                               guide = FALSE)
  return(a)
}

batchqc_f.pvalue <- function(dat, mod, mod0) {
  # F-test (full/reduced model) and returns R2 values
  # (full/reduced) as well.
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

  r2Full <- 1 - rss1 / rss00
  r2Reduced <- 1 - rss0 / rss00

  p <- 1
  if (df1 > df0)  {
    fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
    p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  }
  return(list(p = p, r2Full = r2Full, r2Reduced = r2Reduced))
}
