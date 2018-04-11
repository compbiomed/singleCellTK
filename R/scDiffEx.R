#' Create a heatmap for differential expression analysis
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. Default is "logcounts"
#' @param condition The name of the condition to use for differential
#' expression. Required
#' @param covariates Additional covariates to add to the model. Currently only
#' supported for ANOVA. Default is NULL
#' @param significance FDR corrected significance cutoff for differentially
#' expressed genes. Required
#' @param ntop Number of top differentially expressed genes to display in the
#' heatmap. Required
#' @param usesig If TRUE, only display genes that meet the significance cutoff,
#' up to ntop genes. Required
#' @param diffexmethod The method for performing differential expression
#' analysis. Available options are DESeq2, limma, and ANOVA. Required
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor for condition. The differential
#' expression results will compare the factor in levelofinterest to all other
#' data.
#' @param analysisType For conditions with more than two levels, limma and
#' DESeq2 can be run using multiple methods. See scDiffExlimma() and
#' scDiffExDESeq2() for details.
#' @param controlLevel If the condition has more than two labels, controlLevel
#' should contain one factor from condition to use as the control.
#' @param adjust Method for p-value correction. See options in p.adjust().
#' The default is fdr.
#'
#' @return A list of differentially expressed genes.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffEx(mouseBrainSubsetSCE,
#'                 useAssay = "logcounts",
#'                 "level1class",
#'                 diffexmethod = "limma")
#'
scDiffEx <- function(inSCESet, useAssay="logcounts", condition,
                     covariates=NULL, significance=0.05, ntop=500, usesig=TRUE,
                     diffexmethod, levelofinterest=NULL, analysisType=NULL,
                     controlLevel=NULL, adjust = "fdr"){
  for (i in c(condition, covariates)){
    if (is.factor(SingleCellExperiment::colData(inSCESet)[, i])){
      SummarizedExperiment::colData(inSCESet)[, i] <- droplevels(
        SummarizedExperiment::colData(inSCESet)[, i])
    }
  }
  if (length(condition) == 1){
    if (is.factor(SingleCellExperiment::colData(inSCESet)[, i])){
      in.condition <- droplevels(as.factor(
        SingleCellExperiment::colData(inSCESet)[, condition]))
    }
  } else if (diffexmethod != "ANOVA"){
    stop("Only submit one condition for this method.")
  }

  if (diffexmethod == "DESeq2"){
    diffex.results <- scDiffExDESeq2(inSCESet = inSCESet,
                                      useAssay = useAssay,
                                      condition = condition,
                                      analysisType = analysisType,
                                      levelofinterest = levelofinterest,
                                      controlLevel = controlLevel,
                                      covariates = covariates,
                                      adjust = adjust)
  }
  else if (diffexmethod == "limma"){
    diffex.results <- scDiffExlimma(inSCESet = inSCESet,
                                     useAssay = useAssay,
                                     condition = condition,
                                     analysisType = analysisType,
                                     levelofinterest = levelofinterest,
                                     covariates = covariates,
                                     adjust = adjust)
  }
  else if (diffexmethod == "ANOVA"){
    diffex.results <- scDiffExANOVA(inSCESet = inSCESet,
                                     useAssay = useAssay,
                                     condition = condition,
                                     covariates = covariates,
                                     adjust = adjust)
  }
  else{
    stop("Unsupported differential expression method, ", diffexmethod)
  }
  ngenes <- nrow(inSCESet)
  if (usesig){
    if (length(which(diffex.results$padj <= significance)) < ntop){
      newgenes <- rownames(diffex.results)[
        which(diffex.results$padj <= significance)]
    }
    else{
      newgenes <- rownames(diffex.results)[
        order(diffex.results$padj)[seq_len(min(ntop, ngenes))]]
    }
  }
  else{
    newgenes <- rownames(diffex.results)[
      order(diffex.results$padj)[seq_len(min(ntop, ngenes))]]
  }

  return(diffex.results[newgenes, ])
}

#' Plot Differential Expression
#'
#' @param inSCESet Input data object that contains the data to be plotted.
#' Required
#' @param useAssay Indicate which assay to use. Default is "logcounts"
#' @param condition The condition used for plotting the heatmap. Required
#' @param geneList The list of genes to put in the heatmap. Required
#' @param clusterRow Cluster the rows. The default is TRUE
#' @param clusterCol Cluster the columns. The default is TRUE
#' @param displayRowLabels Display the row labels on the heatmap. The default
#' is TRUE.
#' @param displayColumnLabels Display the column labels on the heatmap. The
#' default is TRUE
#' @param displayRowDendrograms Display the row dendrograms on the heatmap. The
#' default is TRUE
#' @param displayColumnDendrograms Display the column dendrograms on the
#' heatmap. The default is TRUE.
#' @param annotationColors Set of annotation colors for color bar. If null,
#' no color bar is shown. default is NULL.
#' @param columnTitle Title to be displayed at top of heatmap.
#'
#' @return ComplexHeatmap object for the provided geneList annotated with the
#' condition.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffEx(mouseBrainSubsetSCE,
#'                 useAssay = "logcounts",
#'                 "level1class",
#'                 diffexmethod = "limma")
#' plotDiffEx(mouseBrainSubsetSCE, condition = "level1class",
#'             geneList = rownames(res)[1:50], annotationColors = "auto")
#'
plotDiffEx <- function(inSCESet, useAssay="logcounts", condition, geneList,
                        clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=TRUE,
                        displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
                        displayColumnDendrograms=TRUE, annotationColors=NULL,
                        columnTitle="Differential Expression"){
  if (is.null(annotationColors)){
    topha <- NULL
  } else if (annotationColors == "auto") {
    colors <- RColorBrewer::brewer.pal(9, "Set1")
    condLevels <- unique(SingleCellExperiment::colData(inSCESet)[, condition])
    if (length(condLevels) > 9){
      stop("Too many levels in condition for auto coloring")
    }
    col <- list()
    col[[condition]] <- stats::setNames(colors[seq_along(condLevels)],
                                        condLevels)
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCESet)[, condition, drop = FALSE],
      col = col)
  } else {
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCESet)[, condition, drop = FALSE],
      col = annotationColors)
  }

  heatmap <- ComplexHeatmap::Heatmap(
    t(scale(t(SummarizedExperiment::assay(inSCESet, useAssay)[geneList, ]))),
    name = "Expression", column_title = columnTitle, cluster_rows = clusterRow,
    cluster_columns = clusterCol, top_annotation = topha,
    show_row_names = displayRowLabels, show_column_names = displayColumnLabels,
    show_row_dend = displayRowDendrograms,
    show_column_dend = displayColumnDendrograms)
  return(heatmap)
}

#' Perform differential expression analysis with DESeq2
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. Default is "counts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param analysisType Choose "biomarker" to compare the levelofinterest to all
#' other samples. Choose "contrast" to compare the levelofinterest to a
#' controlLevel (see below). Choose "fullreduced" to perform DESeq2 in LRT mode
#' comparing the model with condition to a model without condition.
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor of interest from condition.
#' @param controlLevel If the condition has more than two labels, controlLevel
#' should contain one factor from condition to use as the control.
#' @param covariates Additional covariates to add to the model.
#' @param adjust Method for p-value correction. See options in p.adjust().
#' The default is fdr.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #sort first 100 expressed genes
#' ord <- rownames(mouseBrainSubsetSCE)[
#'   order(rowSums(assay(mouseBrainSubsetSCE, "counts")), 
#'         decreasing = TRUE)][1:100]
#' #subset to those first 100 genes
#' subset <- mouseBrainSubsetSCE[ord, ]
#' res <- scDiffExDESeq2(subset, condition = "level1class")
#'
scDiffExDESeq2 <- function(inSCESet, useAssay="counts", condition,
                            analysisType="biomarker", levelofinterest=NULL,
                            controlLevel=NULL, covariates=NULL, adjust="fdr"){
  cnts <- SummarizedExperiment::assay(inSCESet, useAssay)
  annotData <-
    SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                            drop = FALSE]

  if (is.factor(annotData[, c(condition)])){
    conditionFactor <- factor(annotData[, c(condition)])
    if (length(levels(conditionFactor)) < 2){
      stop("Problem with deseq2 condition")
    } else if (length(levels(conditionFactor)) == 2){
      analysisType <- "standard"
    } else if (is.null(analysisType)){
      stop("You must supply an analysis type")
    } else if (!(analysisType %in% c("standard", "biomarker",
                                     "contrast", "fullreduced"))){
      stop("Unrecognized analysis type, ", analysisType)
    }
  } else {
    analysisType <- "standard"
  }

  if (analysisType == "standard"){
    levelofinterest <- NULL
    controlLevel <- NULL
  } else if (analysisType == "biomarker"){
    if (is.null(levelofinterest)){
      stop("You must specify a level of interest for biomarker analysis.")
    } else {
      annotData[, condition] <- factor(
        ifelse(annotData[, condition] == levelofinterest, levelofinterest,
               paste0("not_", levelofinterest)),
        levels = c(paste0("not_", levelofinterest), levelofinterest))
    }
    controlLevel <- NULL
  } else if (analysisType == "contrast"){
    if (is.null(levelofinterest) || is.null(controlLevel)){
      stop("You must specify a level of interest and a control level for ",
           "contrast analysis.")
    }
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts, colData = annotData,
                                        design = stats::as.formula(
                                          paste0("~", c(condition, covariates),
                                                 collapse = "+")))

  if (analysisType == "fullreduced"){
    if (is.null(covariates)){
      dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1)
    } else {
      dds <- DESeq2::DESeq(dds, test = "LRT",
                           reduced = stats::as.formula(paste0("~",
                                                              c(covariates),
                                                              collapse = "+")))
    }
  } else {
    dds <- DESeq2::DESeq(dds)
  }

  if (analysisType == "standard"){
    res <- DESeq2::results(dds, pAdjustMethod = adjust)
  } else if (analysisType == "fullreduced"){
    res <- DESeq2::results(dds, pAdjustMethod = adjust)
  } else if (analysisType == "contrast") {
    res <- DESeq2::results(dds, contrast = c(condition,
                                             levelofinterest,
                                             controlLevel),
                           pAdjustMethod = adjust)
  } else {
    res <- DESeq2::results(dds, contrast = c(condition,
                                            levelofinterest,
                                            levels(annotData[, condition])[1]),
                          pAdjustMethod = adjust)
  }

  return(data.frame(res))
}

#' Perform differential expression analysis with limma
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. Default is "logcounts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param analysisType If there are more than two levels in your condition
#' variable, select the analysisType. Choose "biomarker" to compare the
#' levelofinterest to all other samples. Choose "coef" to select a coefficient
#' of interset with levelofinterest (see below). Choose "allcoef" to test if
#' any coefficient is different from zero.
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor of interest from condition.
#' @param covariates Additional covariates to add to the model.
#' @param adjust Method for p-value correction. See options in p.adjust().
#' The default is fdr.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffExlimma(mouseBrainSubsetSCE, condition = "level1class")
#'
scDiffExlimma <- function(inSCESet, useAssay="logcounts", condition,
                           analysisType="biomarker",
                           levelofinterest=NULL, covariates=NULL, adjust="fdr"){
  if (is.factor(SingleCellExperiment::colData(inSCESet)[, c(condition)])){
    conditionFactor <- factor(
      SingleCellExperiment::colData(inSCESet)[, c(condition)])
    if (length(levels(conditionFactor)) < 2){
      stop("Problem with limma condition")
    } else if (length(levels(conditionFactor)) == 2){
      analysisType <- "standard"
    } else if (is.null(analysisType)){
      stop("You must supply an analysis type")
    } else if (!(analysisType %in% c("standard", "biomarker",
                                     "coef", "allcoef"))){
      stop("Unrecognized analysis type, ", analysisType)
    }
  } else {
    analysisType <- "standard"
  }

  annotData <- data.frame(
    SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                            drop = FALSE])
  if (analysisType == "biomarker"){
    if (is.null(levelofinterest)){
      stop("You must supply a level of interest")
    }
    annotData[, condition] <- factor(
      ifelse(annotData[, condition] == levelofinterest, levelofinterest,
             paste0("not_", levelofinterest)),
      levels = c(paste0("not_", levelofinterest), levelofinterest))
  }
  design <- stats::model.matrix(
    stats::as.formula(paste0("~", paste0(c(condition, covariates),
                                         collapse = "+"))),
    data = annotData)

  fit <- limma::lmFit(SummarizedExperiment::assay(inSCESet, useAssay), design)
  ebayes <- limma::eBayes(fit)
  if (analysisType == "standard"){
    topGenes <- limma::topTable(ebayes, adjust = adjust,
                                number = nrow(inSCESet))
  } else if (analysisType == "biomarker") {
    topGenes <- limma::topTable(ebayes, coef = 2, adjust = adjust,
                                number = nrow(inSCESet))
  } else if (analysisType == "coef") {
    topGenes <- limma::topTable(
      ebayes, coef = which(levels(conditionFactor) == levelofinterest),
      adjust = adjust, number = nrow(inSCESet))
  } else if (analysisType == "allcoef") {
    topGenes <- limma::topTable(ebayes, adjust = adjust,
                                number = nrow(inSCESet))
  }

  colnames(topGenes)[which(colnames(topGenes) == "adj.P.Val")] <- "padj"
  return(topGenes)
}

#' Perform ANOVA analysis
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. Default is "logcounts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param covariates Additional covariates to add to the model.
#' @param adjust Method for p-value correction. See options in p.adjust().
#' The default is fdr.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffExANOVA(mouseBrainSubsetSCE, condition = "level1class")
#'
scDiffExANOVA <- function(inSCESet, useAssay="logcounts", condition,
                           covariates=NULL, adjust = "fdr"){
  if (is.null(covariates)) {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition), collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition),
                                                                drop = FALSE]))
    mod0 <- stats::model.matrix(~1,
                                data = SingleCellExperiment::colData(inSCESet))
  } else {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition, covariates),
                                           collapse = "+"))),
      data = data.frame(
        SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                                drop = FALSE]))
    mod0 <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(
        SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                                drop = FALSE]))
  }
  dat <- SummarizedExperiment::assay(inSCESet, useAssay)
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
                      t(mod))
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*%
                       t(mod0))
  rss1 <- resid ^ 2 %*% rep(1, n)
  rss0 <- resid0 ^ 2 %*% rep(1, n)
  fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
  p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  results <- data.frame(row.names = rownames(dat), p.value = p,
                        padj = stats::p.adjust(p, method = adjust))
  return(results)
}
