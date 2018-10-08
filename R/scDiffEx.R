#' Perform differential expression analysis on a SCtkExperiment object
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use. Default is "logcounts" for limma
#' and ANOVA, and "counts" for DESeq2.
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param covariates Additional covariates to add to the model. Default is NULL
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
#' expression results will use levelofinterest depending on the analysisType
#' parameter.
#' @param analysisType For conditions with more than two levels, limma and
#' DESeq2 can be run using multiple methods. For DESeq2, choose "biomarker" to
#' compare the levelofinterest to all other samples. Choose "contrast" to
#' compare the levelofinterest to a controlLevel (see below). Choose
#' "fullreduced" to perform DESeq2 in LRT mode. For limma, Choose "biomarker" to
#' compare the levelofinterest to all other samples. Choose "coef" to select a
#' coefficient of interset with levelofinterest (see below). Choose "allcoef" to
#' test if any coefficient is different from zero.
#' @param controlLevel If the condition has more than two labels, controlLevel
#' should contain one factor from condition to use as the control.
#' @param adjust Method for p-value correction. See options in p.adjust().
#' The default is fdr.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffEx(mouseBrainSubsetSCE,
#'                 useAssay = "logcounts",
#'                 "level1class",
#'                 diffexmethod = "limma")
#'
scDiffEx <- function(inSCE, useAssay="logcounts", condition,
                     covariates=NULL, significance=0.05, ntop=500, usesig=TRUE,
                     diffexmethod, levelofinterest=NULL, analysisType=NULL,
                     controlLevel=NULL, adjust = "fdr"){
  #Check for NAs, if true throw error
  if (any(is.na(SingleCellExperiment::colData(inSCE)[, c(condition, covariates)]))){
     stop("Annotation data has NA values. Filter them to continue.")
  }
  for (i in c(condition, covariates)){
    if (is.factor(SingleCellExperiment::colData(inSCE)[, i])){
      SummarizedExperiment::colData(inSCE)[, i] <- droplevels(
        SummarizedExperiment::colData(inSCE)[, i])
    }
  }
  if (length(condition) != 1 & diffexmethod != "ANOVA"){
    stop("Only submit one condition for this method.")
  }

  if (diffexmethod == "DESeq2"){
    diffex.results <- scDiffExDESeq2(inSCE = inSCE,
                                      useAssay = useAssay,
                                      condition = condition,
                                      analysisType = analysisType,
                                      levelofinterest = levelofinterest,
                                      controlLevel = controlLevel,
                                      covariates = covariates,
                                      adjust = adjust)
  }
  else if (diffexmethod == "limma"){
    diffex.results <- scDiffExlimma(inSCE = inSCE,
                                     useAssay = useAssay,
                                     condition = condition,
                                     analysisType = analysisType,
                                     levelofinterest = levelofinterest,
                                     covariates = covariates,
                                     adjust = adjust)
  }
  else if (diffexmethod == "ANOVA"){
    diffex.results <- scDiffExANOVA(inSCE = inSCE,
                                     useAssay = useAssay,
                                     condition = condition,
                                     covariates = covariates,
                                     adjust = adjust)
  }
  else{
    stop("Unsupported differential expression method, ", diffexmethod)
  }
  ngenes <- nrow(inSCE)
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

#' @describeIn scDiffEx Perform differential expression analysis with DESeq2
#'
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
scDiffExDESeq2 <- function(inSCE, useAssay="counts", condition,
                           analysisType="biomarker", levelofinterest=NULL,
                           controlLevel=NULL, covariates=NULL, adjust="fdr"){
  cnts <- SummarizedExperiment::assay(inSCE, useAssay)
  annotData <-
    SingleCellExperiment::colData(inSCE)[, c(condition, covariates),
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

#' @describeIn scDiffEx Perform differential expression analysis with limma
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffExlimma(mouseBrainSubsetSCE, condition = "level1class")
#'
scDiffExlimma <- function(inSCE, useAssay="logcounts", condition,
                          analysisType="biomarker", levelofinterest=NULL,
                          covariates=NULL, adjust="fdr"){
  if (is.factor(SingleCellExperiment::colData(inSCE)[, c(condition)])){
    conditionFactor <- factor(
      SingleCellExperiment::colData(inSCE)[, c(condition)])
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
    SingleCellExperiment::colData(inSCE)[, c(condition, covariates),
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

  fit <- limma::lmFit(SummarizedExperiment::assay(inSCE, useAssay), design)
  ebayes <- limma::eBayes(fit)
  if (analysisType == "standard"){
    topGenes <- limma::topTable(ebayes, adjust = adjust,
                                number = nrow(inSCE))
  } else if (analysisType == "biomarker") {
    topGenes <- limma::topTable(ebayes, coef = 2, adjust = adjust,
                                number = nrow(inSCE))
  } else if (analysisType == "coef") {
    topGenes <- limma::topTable(
      ebayes, coef = which(levels(conditionFactor) == levelofinterest),
      adjust = adjust, number = nrow(inSCE))
  } else if (analysisType == "allcoef") {
    topGenes <- limma::topTable(ebayes, adjust = adjust,
                                number = nrow(inSCE))
  }

  colnames(topGenes)[which(colnames(topGenes) == "adj.P.Val")] <- "padj"
  return(topGenes)
}

#' @describeIn scDiffEx Perform differential expression analysis with ANOVA
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffExANOVA(mouseBrainSubsetSCE, condition = "level1class")
#'
scDiffExANOVA <- function(inSCE, useAssay="logcounts", condition,
                          covariates=NULL, adjust = "fdr"){
  if (is.null(covariates)) {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition), collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCE)[, c(condition),
                                                                drop = FALSE]))
    mod0 <- stats::model.matrix(~1,
                                data = SingleCellExperiment::colData(inSCE))
  } else {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition, covariates),
                                           collapse = "+"))),
      data = data.frame(
        SingleCellExperiment::colData(inSCE)[, c(condition, covariates),
                                                drop = FALSE]))
    mod0 <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(
        SingleCellExperiment::colData(inSCE)[, c(condition, covariates),
                                                drop = FALSE]))
  }
  dat <- SummarizedExperiment::assay(inSCE, useAssay)
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

#' Plot Differential Expression
#'
#' @param inSCE Input data object that contains the data to be plotted.
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
#' @param scaleExpression Row scale the heatmap values. The default is TRUE.
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
plotDiffEx <- function(inSCE, useAssay="logcounts", condition, geneList,
                       clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=TRUE,
                       displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
                       displayColumnDendrograms=TRUE, annotationColors=NULL,
                       scaleExpression=TRUE,
                       columnTitle="Differential Expression"){
  if (is.null(annotationColors)){
    topha <- NULL
  } else if (annotationColors == "auto") {
    if (length(condition) != 1) {
      stop("Only one condition may be used for auto coloring.")
    }
    condLevels <- unique(SingleCellExperiment::colData(inSCE)[, condition])
    if (length(condLevels) > 9){
      colors <- distinctColors(length(condLevels))
    } else {
      colors <- RColorBrewer::brewer.pal(9, "Set1")
    }
    col <- list()
    col[[condition]] <- stats::setNames(colors[seq_along(condLevels)],
                                        condLevels)
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCE)[, condition, drop = FALSE],
      col = col)
  } else {
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCE)[, condition, drop = FALSE],
      col = annotationColors)
  }
  heatmapkey <- useAssay
  heatdata <- SummarizedExperiment::assay(inSCE, useAssay)[geneList, ]
  if (scaleExpression){
    heatdata <- t(scale(t(heatdata)))
    heatmapkey <- paste("Scaled", heatmapkey, sep = "\n")
  }

  heatmap <- ComplexHeatmap::Heatmap(
    heatdata,
    name = heatmapkey, column_title = columnTitle, cluster_rows = clusterRow,
    cluster_columns = clusterCol, top_annotation = topha,
    show_row_names = displayRowLabels, show_column_names = displayColumnLabels,
    show_row_dend = displayRowDendrograms,
    show_column_dend = displayColumnDendrograms)
  return(heatmap)
}
