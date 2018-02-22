#' Create a heatmap for differential expression analysis
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param use_assay Indicate which assay to use. Default is "logcounts"
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
#' @param clusterRow Cluster the rows. The default is TRUE
#' @param clusterCol Cluster the columns. The default is TRUE
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor for condition. The differential
#' expression results will compare the factor in levelofinterest to all other
#' data.
#' @param analysis_type Choose "biomarker" to compare the levelofinterest to all
#' other samples. Choose "contrast" to compare the levelofinterest to a
#' controlLevel.
#' @param controlLevel If the condition has more than two labels, controlLevel
#' should contain one factor from condition to use as the control.
#'
#' @return A list of differentially expressed genes.
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' res <- scDiffEx(GSE60361_subset_sce,
#'                 use_assay = "logcounts",
#'                 "level1class",
#'                 diffexmethod = "limma")
#'
scDiffEx <- function(inSCESet, use_assay="logcounts", condition,
                     covariates=NULL, significance=0.05, ntop=500, usesig=TRUE,
                     diffexmethod, clusterRow=TRUE, clusterCol=TRUE,
                     levelofinterest=NULL, analysis_type=NULL, controlLevel=NULL){
  if(length(condition) == 1){
    in.condition <- droplevels(as.factor(SingleCellExperiment::colData(inSCESet)[, condition]))
  } else if (diffexmethod != "ANOVA"){
    stop("Only submit one condition for this method.")
  }

  if (diffexmethod == "DESeq2"){
    diffex.results <- scDiffEx_deseq2(inSCESet = inSCESet,
                                      use_assay = use_assay,
                                      condition = condition,
                                      analysis_type = analysis_type,
                                      levelofinterest = levelofinterest,
                                      controlLevel = controlLevel,
                                      covariates = covariates)
  }
  else if (diffexmethod == "limma"){
    diffex.results <- scDiffEx_limma(inSCESet = inSCESet,
                                     use_assay = use_assay,
                                     condition = condition,
                                     levelofinterest = levelofinterest,
                                     covariates = covariates)
  }
  else if (diffexmethod == "ANOVA"){
    diffex.results <- scDiffEx_anova(inSCESet = inSCESet,
                                     use_assay = use_assay,
                                     condition = condition,
                                     covariates = covariates)
  }
  else{
    stop("Unsupported differential expression method, ", diffexmethod)
  }

  if (usesig){
    if (length(which(diffex.results$padj <= significance)) < ntop){
      newgenes <- rownames(diffex.results)[which(diffex.results$padj <= significance)]
    }
    else{
      newgenes <- rownames(diffex.results)[order(diffex.results$padj)[1:ntop]]
    }
  }
  else{
    newgenes <- rownames(diffex.results)[order(diffex.results$padj)[1:ntop]]
  }

  return(diffex.results[newgenes, ])
}

#' Plot Differential Expression
#'
#' @param inSCESet Input data object that contains the data to be plotted.
#' Required
#' @param use_assay Indicate which assay to use. Default is "logcounts"
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
#' data("GSE60361_subset_sce")
#' res <- scDiffEx(GSE60361_subset_sce,
#'                 use_assay = "logcounts",
#'                 "level1class",
#'                 diffexmethod = "limma")
#' plot_DiffEx(GSE60361_subset_sce, condition = "level1class",
#'             geneList = rownames(res)[1:50], annotationColors = "auto")
#'
plot_DiffEx <- function(inSCESet, use_assay="logcounts", condition, geneList,
                        clusterRow=TRUE, clusterCol=TRUE, displayRowLabels=TRUE,
                        displayColumnLabels=TRUE, displayRowDendrograms=TRUE,
                        displayColumnDendrograms=TRUE, annotationColors=NULL,
                        columnTitle="Differential Expression"){
  if (is.null(annotationColors)){
    topha <- NULL
  } else if (annotationColors == "auto") {
    colors <- RColorBrewer::brewer.pal(9, "Set1")
    cond_levels <- unique(SingleCellExperiment::colData(inSCESet)[, condition])
    if (length(cond_levels) > 9){
      stop("Too many levels in condition for auto coloring")
    }
    col <- list()
    col[[condition]] <- stats::setNames(colors[1:length(cond_levels)],
                                        cond_levels)
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCESet)[, condition, drop = FALSE], col = col)
  } else {
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = SingleCellExperiment::colData(inSCESet)[, condition, drop = FALSE], col = annotationColors)
  }

  heatmap <- ComplexHeatmap::Heatmap(t(scale(t(SummarizedExperiment::assay(inSCESet, use_assay)[geneList, ]))),
                                     name = "Expression",
                                     column_title = columnTitle,
                                     cluster_rows = clusterRow,
                                     cluster_columns = clusterCol,
                                     top_annotation = topha,
                                     show_row_names = displayRowLabels,
                                     show_column_names = displayColumnLabels,
                                     show_row_dend = displayRowDendrograms,
                                     show_column_dend = displayColumnDendrograms)
  return(heatmap)
}

#' Perform differential expression analysis with DESeq2
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param use_assay Indicate which assay to use. Default is "counts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param analysis_type Choose "biomarker" to compare the levelofinterest to all
#' other samples. Choose "contrast" to compare the levelofinterest to a
#' controlLevel (see below).
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor of interest from condition.
#' @param controlLevel If the condition has more than two labels, controlLevel
#' should contain one factor from condition to use as the control.
#' @param covariates Additional covariates to add to the model.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' #subset to 100 genes
#' subset <- GSE60361_subset_sce[rownames(GSE60361_subset_sce)[order(rowSums(assay(GSE60361_subset_sce, "counts")), decreasing = TRUE)][1:100], ]
#' res <- scDiffEx_deseq2(subset, condition = "level1class")
#'
scDiffEx_deseq2 <- function(inSCESet, use_assay="counts", condition,
                            analysis_type="biomarker", levelofinterest=NULL,
                            controlLevel=NULL, covariates=NULL){
  cnts <- SummarizedExperiment::assay(inSCESet, use_assay)
  annot_data <- SingleCellExperiment::colData(inSCESet)[, c(condition, covariates), drop = FALSE]

  if (length(levels(annot_data[, condition])) > 2){
    if(analysis_type == "biomarker"){
      if(is.null(levelofinterest)){
        stop("You must specify a level of interest for biomarker analysis.")
      } else {
        annot_data[, condition] <- factor(ifelse(annot_data[,condition] == levelofinterest,
                                                 levelofinterest,
                                                 paste0("not_", levelofinterest)),
                                          levels = c(paste0("not_", levelofinterest),
                                                     levelofinterest))
      }
      controlLevel <- NULL
    } else if (analysis_type == "contrast"){
      if(is.null(levelofinterest) || is.null(controlLevel)){
        stop("You must specify a level of interest and a control level for contrast analysis.")
      }
    } else {
      stop("Invalid analysis type: ", analysis_type)
    }
  } else {
    levelofinterest <- NULL
    controlLevel <- NULL
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts, colData = annot_data,
                                        design = stats::as.formula(
                                          paste0("~", c(condition, covariates),
                                                 collapse = "+")))
  dds <- DESeq2::DESeq(dds)
  if(is.null(levelofinterest) && is.null(controlLevel)){
    res <- DESeq2::results(dds, contrast = c(condition,
                                             levels(annot_data[,condition])[2],
                                             levels(annot_data[,condition])[1]))
    res <- DESeq2::lfcShrink(dds, coef = 2)
  } else if(is.null(controlLevel)){
    res <- DESeq2::results(dds, contrast = c(condition,
                                             levelofinterest,
                                             levels(annot_data[,condition])[1]))
    res <- DESeq2::lfcShrink(dds, coef = which(levels(annot_data[,condition]) == levelofinterest))
  } else {
    res <- DESeq2::results(dds, contrast = c(condition,
                                             levelofinterest,
                                             controlLevel))
    res <- DESeq2::lfcShrink(dds, coef = which(levels(SingleCellExperiment::colData(inSCESet)[, c(condition)]) == levelofinterest))
  }

  return(data.frame(res))
}

#' Perform differential expression analysis with limma
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param use_assay Indicate which assay to use. Default is "logcounts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor of interest from condition.
#' @param covariates Additional covariates to add to the model.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' res <- scDiffEx_limma(GSE60361_subset_sce, condition = "level1class")
#'
scDiffEx_limma <- function(inSCESet, use_assay="logcounts", condition,
                           levelofinterest=NULL, covariates=NULL){
  condition_factor <- factor(SingleCellExperiment::colData(inSCESet)[, c(condition)])
  if (length(levels(condition_factor)) < 2){
    stop("Problem with limma condition")
  }
  if(is.null(covariates)) {
    design <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition), collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition), drop = FALSE]))
  } else {
    design <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition, covariates),
                                           collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                          drop = FALSE]))
  }
  fit <- limma::lmFit(SummarizedExperiment::assay(inSCESet, use_assay), design)
  ebayes <- limma::eBayes(fit)
  if (length(levels(condition_factor)) == 2){
    topGenes <- limma::topTable(ebayes, coef = 2, adjust = "fdr",
                                number = nrow(inSCESet))
  } else {
    if(is.null(levelofinterest)){
      stop("Level of interest required.")
    }
    topGenes <- limma::topTable(
      ebayes, coef = which(levels(condition_factor) == levelofinterest),
      adjust = "fdr", number = nrow(inSCESet))
  }

  colnames(topGenes)[5] <- "padj"
  return(topGenes)
}

#' Perform ANOVA analysis
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param use_assay Indicate which assay to use. Default is "logcounts"
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from colData that contains at least
#' two labels. Required
#' @param covariates Additional covariates to add to the model.
#'
#' @return A data frame of gene names and adjusted p-values
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' res <- scDiffEx_anova(GSE60361_subset_sce, condition = "level1class")
#'
scDiffEx_anova <- function(inSCESet, use_assay="logcounts", condition,
                           covariates=NULL){
  if(is.null(covariates)) {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition), collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition), drop = FALSE]))
    mod0 <- stats::model.matrix(~1, data = SingleCellExperiment::colData(inSCESet))
  } else {
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(c(condition, covariates),
                                           collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                          drop = FALSE]))
    mod0 <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCESet)[, c(condition, covariates),
                                          drop = FALSE]))
  }
  dat <- SummarizedExperiment::assay(inSCESet, use_assay)
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
                        padj = stats::p.adjust(p, method = "fdr"))
  return(results)
}
