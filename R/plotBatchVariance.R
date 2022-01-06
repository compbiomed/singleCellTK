.searchBCDefaultInfo <- function(inSCE, corrMat, origAssay, matType) {
  if (is.null(origAssay)) {
    if ("counts" %in% expDataNames(inSCE)) {
      origAssay <- "counts"
    } else {
      origAssay <- expDataNames(inSCE)[1]
    }
    warning("using '", origAssay, "' for comparison.")
  }

  if (is.null(matType)) {
    if (corrMat %in% SummarizedExperiment::assayNames(inSCE)) {
      matType <- "assay"
    } else if (corrMat %in% SingleCellExperiment::altExpNames(inSCE)) {
      matType <- "altExp"
    } else if (corrMat %in% SingleCellExperiment::reducedDimNames(inSCE)) {
      matType <- "reducedDim"
    } else {
      stop("Corrected Matrix name '", corrMat, "' not found in inSCE")
    }
  }

  return(c(origAssay, matType))
}

.checkBCMeta <- function(inSCE, corrMat, origAssay, origLogged, method, matType,
                         batch, condition) {
  if (!is.null(matType)) {
    if (!matType %in% c("assay", "altExp", "reducedDim")) {
      stop("Wrong matrix type '", matType, "'. Choose from 'assay', 'altExp', ",
           "'reducedDim'.")
    }
  }
  if (!"batchCorr" %in% names(S4Vectors::metadata(inSCE))) {
    warning("Batch correction result from SCTK not found.")
    s <- .searchBCDefaultInfo(inSCE, corrMat, origAssay, matType)
    origAssay <- ifelse(is.null(origAssay), s[1], origAssay)
    method <- ifelse(is.null(method), "Unidentified Method", method)
    matType <- ifelse(is.null(matType), s[2], matType)
  } else {
    if (!corrMat %in% names(S4Vectors::metadata(inSCE)$batchCorr)) {
      warning("'", corrMat, "' not identified as a Batch correction result ",
              "from SCTK")
      s <- .searchBCDefaultInfo(inSCE, corrMat, origAssay, matType)
      origAssay <- ifelse(is.null(origAssay), s[1], origAssay)
      method <- ifelse(is.null(method), "Unidentified Method", method)
      matType <- ifelse(is.null(matType), s[2], matType)
    } else {
      bcInfo <- S4Vectors::metadata(inSCE)$batchCorr[[corrMat]]
      origAssay <- ifelse(is.null(origAssay), bcInfo$useAssay, origAssay)
      origLogged <- ifelse(is.null(origLogged), bcInfo$origLogged, origLogged)
      method <- ifelse(is.null(method), bcInfo$method, method)
      if (!is.null(matType) && matType != bcInfo$matType) {
        warning("User specified matType different from SCTK identified ",
                "matType. Force using user specification.")
      }
      matType <- ifelse(is.null(matType), bcInfo$matType, matType)
      batch <- ifelse(is.null(batch), bcInfo$batch, batch)
      if (is.null(condition)) condition <- bcInfo$condition
      #condition <- ifelse(is.null(condition), bcInfo$condition, condition)
    }
  }
  return(list(origAssay = origAssay,
              origLogged = origLogged,
              method = method,
              matType = matType,
              batch = batch,
              condition = condition))
}

#' Plot comparison of batch corrected result against original assay
#' @details Four plots will be combined. Two of them are violin/box-plots for
#' percent variance explained by the batch variation, and optionally the
#' covariate, for original and corrected. The other two are UMAPs of the
#' original assay and the correction result matrix. If SCTK batch correction
#' methods are performed in advance, this function will automatically detect
#' necessary input. Otherwise, users can also customize the input. Future
#' improvement might include solution to reduce redundant UMAP calculation.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param corrMat A single character indicating the name of the corrected matrix.
#' @param batch A single character. The name of batch annotation column in
#' \code{colData(inSCE)}.
#' @param condition A single character. The name of an additional covariate
#' annotation column in \code{colData(inSCE)}.
#' @param origAssay A single character indicating what the original assay used
#' for batch correction is.
#' @param origLogged Logical scalar indicating whether \code{origAssay} is
#' log-normalized.
#' @param method A single character indicating the name of the batch correction
#' method. Only used for the titles of plots.
#' @param matType A single character indicating the type of the batch correction
#' result matrix, choose from \code{"assay"}, \code{"altExp"},
#' \code{"reducedDim"}.
#' @return An object of class \code{"gtable"}, combining four \code{ggplot}s.
#' @examples
#' data("sceBatches")
#' sceBatches <- scaterlogNormCounts(sceBatches, "logcounts")
#' sceBatches <- runLimmaBC(sceBatches)
#' plotBatchCorrCompare(sceBatches, "LIMMA", condition = "cell_type")
#' @export
#' @author Yichen Wang
plotBatchCorrCompare <- function(inSCE, corrMat, batch = NULL, condition = NULL,
                                 origAssay = NULL, origLogged = NULL,
                                 method = NULL, matType = NULL) {
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  m <- .checkBCMeta(inSCE, corrMat, origAssay, origLogged, method, matType,
                    batch, condition)
  origAssay <- m$origAssay
  origLogged <- m$origLogged
  method <- m$method
  matType <- m$matType
  batch <- m$batch
  condition <- m$condition

  if (isFALSE(origLogged)) {
    inSCE <- scaterlogNormCounts(inSCE, origAssay, origAssay)
  }

  # Batch Variance Plot for origAssay
  bv.before <- plotBatchVariance(inSCE, useAssay = origAssay, useReddim = NULL,
                                 useAltExp = NULL, batch = batch,
                                 condition = condition,
                                 title = "Batch Variance before correction") +
    ggplot2::theme(text=ggplot2::element_text(size=10))

  inSCE <- getUMAP(inSCE, useAssay = origAssay, reducedDimName = "umap.before")
  umap.before <- plotSCEDimReduceColData(inSCE, batch, "umap.before",
                                         shape = condition, axisLabelSize = 9,
                                         axisSize = 8, dotSize = 1,
                                         titleSize = 12, labelClusters = FALSE,
                                         legendSize = 10, legendTitle = "batch",
                                         legendTitleSize = 10,
                                         title = "UMAP before correction")

  if (matType == "assay") {
    if (isFALSE(origLogged)) {
      inSCE <- scaterlogNormCounts(inSCE, corrMat, corrMat)
    }
    # Batch Variance Plot for CorrMat
    bv.after <- plotBatchVariance(inSCE, useAssay = corrMat, batch = batch,
                                  condition = condition,
                                  title = paste0("Batch Variance corrected with ",
                                                 method)) +
      ggplot2::theme(text=ggplot2::element_text(size=10))

    if (method == "ComBatSeq") {
      inSCE <- getUMAP(inSCE, useAssay = corrMat, reducedDimName = "umap.after")
    } else {
      inSCE <- getUMAP(inSCE, useAssay = corrMat, reducedDimName = "umap.after",
                       logNorm = FALSE)
    }
  } else if (matType == "altExp") {
    # Doing log, because only Seurat returns altExp,
    # and the assay inside is not logged
    ae <- SingleCellExperiment::altExp(inSCE, corrMat)
    ae <- scaterlogNormCounts(ae, corrMat, corrMat)
    SingleCellExperiment::altExp(inSCE, corrMat) <- ae
    bv.after <- plotBatchVariance(inSCE, useAltExp = corrMat, batch = batch,
                                  condition = condition,
                                  title = paste0("Batch Variance corrected with ",
                                                 method)) +
      ggplot2::theme(text=ggplot2::element_text(size=10))
    inSCE <- getUMAP(inSCE, useAltExp = corrMat, useAssay = corrMat,
                     reducedDimName = "umap.after")
  } else if (matType == "reducedDim") {
    bv.after <- plotBatchVariance(inSCE, useReddim = corrMat, batch = batch,
                                  condition = condition,
                                  title = paste0("Batch Variance corrected with ",
                                                 method)) +
      ggplot2::theme(text=ggplot2::element_text(size=10))
    if (method == "BBKNN") {
      SingleCellExperiment::reducedDim(inSCE, "umap.after") <-
        SingleCellExperiment::reducedDim(inSCE, corrMat)
    } else {
      inSCE <- getUMAP(inSCE, useAssay = NULL, useReducedDim = corrMat,
                       reducedDimName = "umap.after")
    }
  } else {
    stop("Cannot identify result matrix type")
  }
  umap.after <- plotSCEDimReduceColData(inSCE, batch, "umap.after", dim1 = 1,
                                        dim2 = 2,
                                        shape = condition, axisLabelSize = 9,
                                        axisSize = 8, dotSize = 1,
                                        titleSize = 12, labelClusters = FALSE,
                                        legendSize = 10, legendTitle = "batch",
                                        legendTitleSize = 10,
                                        title = "UMAP after correction") +
    ggplot2::theme(text=ggplot2::element_text(size=8))
  return(gridExtra::grid.arrange(bv.before, bv.after,
                                 umap.before, umap.after, nrow = 2))
}

#' Plot the percent of the variation that is explained by batch and condition
#' in the data
#'
#' Visualize the percent variation in the data that is explained by batch and
#' condition, individually, and that explained by combining both annotations.
#' Plotting only the variation explained by batch is supported but not
#' recommended, because this can be confounded by potential condition.
#'
#' When condition and batch both are causing some variation, if the difference
#' between full variation and condition variation is close to batch variation,
#' this might imply that batches are causing some effect; if the difference is
#' much less than batch variation, then the batches are likely to be confounded
#' by the conditions.
#'
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay A single character. The name of the assay that stores the
#' value to plot. For \code{useReddim} and \code{useAltExp} also. Default
#' \code{NULL}.
#' @param useReddim A single character. The name of the dimension reduced
#' matrix that stores the value to plot. Default \code{NULL}.
#' @param useAltExp A single character. The name of the alternative experiment
#' that stores an assay of the value to plot. Default \code{NULL}.
#' @param batch A single character. The name of batch annotation column in
#' \code{colData(inSCE)}. Default \code{"batch"}.
#' @param condition A single character. The name of an additional condition
#' annotation column in \code{colData(inSCE)}. Default \code{NULL}.
#' @param title A single character. The title text on the top. Default
#' \code{NULL}.
#' @return A ggplot object of a boxplot of variation explained by batch,
#' condition, and batch+condition.
#' @export
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' plotBatchVariance(sceBatches,
#'                   useAssay="counts",
#'                   batch="batch",
#'                   condition = "cell_type")
plotBatchVariance <- function(inSCE, useAssay = NULL, useReddim = NULL,
                              useAltExp = NULL, batch = 'batch',
                              condition=NULL, title = NULL) {
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop("'inSCE' must inherit from 'SingleCellExperiment'.")
  }
  if(is.null(useAssay) + is.null(useReddim) + is.null(useAltExp) != 2){
    stop("One and only one of `useAssay`, `useReddim`, ",
         "`usAltExp` has to be specified.")
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
  if(!is.null(useAltExp)){
    if(!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)){
      stop("'useAltExp not found in 'inSCE'.")
    }
    ae <- SingleCellExperiment::altExp(inSCE, useAltExp)
    mat <- SummarizedExperiment::assay(ae)
  }
  if(is.null(batch)){
    stop("Batch annotation has to be given.")
  } else{
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
      stop("'batch' not found in 'inSCE'.")
    }
  }
  if(!inherits(mat, 'matrix')){
    mat <- as.matrix(mat)
  }
  batchCol <- SummarizedExperiment::colData(inSCE)[, batch]
  nlb <- nlevels(as.factor(batchCol))
  if (nlb <= 1){
    stop("No more than one batch found in specified annotation")
  } else {
    batchMod <- stats::model.matrix(~as.factor(batchCol))
  }
  if (is.null(condition)){
    condMod <- matrix(rep(1, ncol(mat)), ncol = 1)
  } else {
    condCol <- SingleCellExperiment::colData(inSCE)[, condition]
    nlc <- nlevels(as.factor(condCol))
    if (nlc <= 1){
      condMod <- matrix(rep(1, ncol(mat)), ncol = 1)
    } else {
      condMod <- stats::model.matrix(~as.factor(condCol))
    }
  }
  mod <- cbind(condMod, batchMod[, -1])
  condTest <- .batchqc_f.pvalue(mat, mod, batchMod)
  batchTest <- .batchqc_f.pvalue(mat, mod, condMod)
  r2Full <- condTest$r2Full
  condR2 <- batchTest$r2Reduced
  batchR2 <- condTest$r2Reduced
  explainedVariation <- round(cbind(`Full (Condition+Batch)` = r2Full,
                                     Condition = condR2,
                                     Batch = batchR2), 5) * 100
  colnames(explainedVariation) <- c("Full",
                                    ifelse(is.null(condition), "No Condition", condition),
                                    batch)
  exVarM <- reshape2::melt(explainedVariation)
  colnames(exVarM) <- c("Gene", "Model", "Percent.Explained.Variation")
  exVarM$Model <- factor(exVarM$Model)
  a <- ggplot2::ggplot(exVarM,
                       ggplot2::aes_string("Model",
                                           "Percent.Explained.Variation")) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2),
                        size = 1, alpha = 0.9) +
    ggplot2::geom_violin(ggplot2::aes_string(fill = "Model"), alpha = 0.7, ) +
    ggplot2::geom_boxplot(alpha = 0.4, width = 0.2) +
    ggplot2::ylim(0, 100) +
    ggplot2::xlab("Model") +
    ggplot2::ylab("Explained Variation %") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  return(a)
}

.batchqc_f.pvalue <- function(dat, mod, mod0) {
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
#' @param useAssay A single character. The name of the assay that stores the
#' value to plot. For \code{useReddim} and \code{useAltExp} also. Default
#' \code{NULL}.
#' @param useReddim A single character. The name of the dimension reduced
#' matrix that stores the value to plot. Default \code{NULL}.
#' @param useAltExp A single character. The name of the alternative experiment
#' that stores an assay of the value to plot. Default \code{NULL}.
#' @param batch A single character. The name of batch annotation column in
#' \code{colData(inSCE)}. Default \code{"batch"}.
#' @param xlab label for x-axis. Default \code{"batch"}.
#' @param ylab label for y-axis. Default \code{"Feature Mean"}.
#' @param ... Additional arguments passed to \code{\link{.ggViolin}}.
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' plotSCEBatchFeatureMean(sceBatches, useAssay = "counts")
#' @return ggplot
#' @export
plotSCEBatchFeatureMean <- function(inSCE, useAssay = NULL, useReddim = NULL,
  useAltExp = NULL, batch = 'batch', xlab='batch', ylab='Feature Mean', ...){
  if(!inherits(inSCE, 'SingleCellExperiment')){
    stop("'inSCE' must inherit from 'SingleCellExperiment'.")
  }
  if(is.null(useAssay) + is.null(useReddim) + is.null(useAltExp) != 2){
    stop("One and only one of `useAssay`, `useReddim`, ",
         "`usAltExp` has to be specified.")
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
  if(!is.null(useAltExp)){
    if(!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)){
      stop("'useAltExp not found in 'inSCE'.")
    }
    ae <- SingleCellExperiment::altExp(inSCE, useAltExp)
    mat <- SummarizedExperiment::assay(ae)
  }
  if(is.null(batch)){
    stop("Batch annotation has to be given.")
  } else{
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
      stop("'batch' not found in 'inSCE'.")
    }
  }
  if(!inherits(mat, 'matrix')){
    mat <- as.matrix(mat)
  }
  batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
  uniqBatch <- as.vector(unique(batchCol)) #as.vector in case batchCol is factor
  allMeans <- numeric()
  groupBy <- character()
  for(i in uniqBatch){
    allMeans <- c(allMeans, DelayedArray::rowMeans(mat[,batchCol == i]))
    groupBy <- c(groupBy, rep(i, nrow(mat)))
  }
  p <- .ggViolin(allMeans, groupBy = groupBy, xlab = xlab, ylab = ylab, ...)
  p <- .ggSCTKTheme(p)
  return(p)
}
