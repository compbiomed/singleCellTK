#' Plot the percent of the variation that is explained by batch and condition
#' in the data
#'
#' Visualize the percent variation in the data that is explained by batch and
#' condition if it is given.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"
#' @param batch The column in the annotation data that corresponds to batch.
#' Required
#' @param condition The column in the annotation data that corresponds to
#' condition. Optional
#' @param pcInput A logical scalar indicating whether \code{useAssay} is in
#' \code{names(reducedDims(inSCE))}. Default \code{FALSE}.
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
plotBatchVariance <- function(inSCE, useAssay="logcounts", batch='batch',
                              condition=NULL, pcInput = FALSE){
  if(isTRUE(pcInput)){
    mat <- t(SingleCellExperiment::reducedDim(inSCE, useAssay))
  } else {
    mat <- SummarizedExperiment::assay(inSCE, useAssay)
  }
  batchCol <- SummarizedExperiment::colData(inSCE)[, batch]
  nlb <- nlevels(as.factor(batchCol))
  if (nlb <= 1){
    batchMod <- matrix(rep(1, ncol(inSCE)), ncol = 1)
  } else {
    batchMod <- stats::model.matrix(~as.factor(batchCol))
  }
  condCol <- SingleCellExperiment::colData(inSCE)[, condition]
  if (is.null(condition)){
    stop("condition required for now")
  } else {
    nlc <- nlevels(as.factor(condCol))
    if (nlc <= 1){
      condMod <- matrix(rep(1, ncol(inSCE)), ncol = 1)
    } else {
      condMod <- stats::model.matrix(~as.factor(condCol))
    }
  }
  mod <- cbind(condMod, batchMod[, -1])
  condTest <- batchqc_f.pvalue(mat, mod, batchMod)
  batchTest <- batchqc_f.pvalue(mat, mod, condMod)
  r2Full <- condTest$r2Full
  condR2 <- batchTest$r2Reduced
  batchR2 <- condTest$r2Reduced
  explainedVariation <- round(cbind(`Full (Condition+Batch)` = r2Full,
                                     Condition = condR2,
                                     Batch = batchR2), 5) * 100
  exVarM <- reshape2::melt(explainedVariation)
  colnames(exVarM) <- c("Gene", "Model", "Percent.Explained.Variation")
  exVarM$Model <- factor(exVarM$Model)
  a <- ggplot2::ggplot(exVarM, ggplot2::aes_string("Model", "Percent.Explained.Variation")) +
    ggplot2::geom_violin(ggplot2::aes_string(fill = "Model")) +
    ggplot2::geom_boxplot(width = .1) +
    ggplot2::xlab("Model") +
    ggplot2::ylab("Percent Explained Variation") +
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

#' Plot mean feature value in each batch of a SingleCellExperiment object
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay The name of the assay that stores the value to plot. Use
#' \code{useReddim} for dimension reduced matrix instead. Default \code{NULL}.
#' @param useReddim The name of the dimension reduced matrix that stores the
#' value to plot. Default \code{NULL}.
#' @param batch The column name of \code{colData(inSCE)} that indicates the
#' batch annotation. Default \code{"batch"}.
#' @param xlab label for x-axis. Default \code{"batch"}.
#' @param ylab label for y-axis. Default \code{"Feature Mean"}.
#' @param ... Additional arguments passed to \code{\link{.ggViolin}}.
#' @return ggplot
#' @export
plotSCEBatchFeatureMean <- function(inSCE, useAssay = NULL, useReddim = NULL,
  batch = 'batch', xlab='batch', ylab='Feature Mean', ...){
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop("'inSCE' must inherit from 'SingleCellExperiment'.")
  }
  if(is.null(useAssay) & is.null(useReddim)){
    stop("Either `useAssay` or `useReddim` has to be specified.")
  } else if(!is.null(useAssay) & !is.null(useReddim)){
    stop("Only one of `useAssay` and `useReddim` can be specified.")
  }
  if(!is.null(useAssay)){
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
      stop("'useAssay' not found in 'inSCE'.")
    }
    mat <- SummarizedExperiment::assay(inSCE, useAssay)
  }
  if(!is.null(useReddim)){
    if(!useReddim %in% SingleCellExperiment::reducedDimNames(inSCE)){
      stop("'useReddim not found in 'inSCE'.")
    }
    mat <- t(SingleCellExperiment::reducedDim(inSCE, useReddim))
  }
  if(is.null(batch)){
    stop("Batch annotation has to be given.")
  } else{
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
      stop("'batch' not found in 'inSCE'.")
    }
  }
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  uniqBatch <- as.vector(unique(batchCol)) #as.vector in case batchCol is factor
  allMeans <- numeric()
  groupBy <- character()
  for(i in uniqBatch){
    allMeans <- c(allMeans, DelayedArray::rowMeans(mat[,batchCol == i]))
    groupBy <- c(groupBy, rep(i, nrow(inSCE)))
  }
  p <- .ggViolin(allMeans, groupby = groupBy, xlab = xlab, ylab = ylab, ...)
  return(p)
}