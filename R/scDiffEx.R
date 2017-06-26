#' Create a heatmap for differential expression analysis
#'
#' @param inSCESet Input SCESet object. Required
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from pData that contains two labels.
#' Required
#' @param significance FDR corrected significance cutoff for differentially
#' expressed genes. Required
#' @param ntop Number of top differentially expressed genes to display in the
#' heatmap. Required
#' @param usesig If TRUE, only display genes that meet the significance cutoff,
#' up to ntop genes. Required
#' @param diffexmethod The method for performing differential expression
#' analyis. Available options are DESeq, DESeq2, and limma. Required
#' @param clusterRow Cluster the rows. The default is TRUE
#' @param clusterCol Cluster the columns. The default is TRUE
#' @param displayRowLabels Display the row labels on the heatmap. The default
#' is TRUE.
#' @param levelofinterest If the condition has more than two labels,
#' levelofinterest should contain one factor for condition. The differential
#' expression results will compare the factor in levelofinterest to all other
#' data.
#'
#' @return A list of differentially expressed genes.
#' @export scDiffEx
#'
#' @examples
#' \dontrun{
#' scDiffEx(newsceset_david, "age", "0.1")
#' }
#'
scDiffEx <- function(inSCESet, condition, significance=0.05, ntop=500,
                     usesig=TRUE, diffexmethod, clusterRow=TRUE,
                     clusterCol=TRUE, displayRowLabels=TRUE,
                     levelofinterest){
  in.condition <- droplevels(as.factor(Biobase::pData(inSCESet)[, condition]))

  if (length(levels(in.condition)) < 2){
    stop("You must submit a condition with more than 1 labels: ", condition,
         " has ", length(levels(in.condition)), " labels")
  } else if (length(levels(in.condition)) > 2){
    in.condition <- droplevels(as.factor(ifelse(in.condition == levelofinterest,
                           levelofinterest,
                           paste("not", levelofinterest, sep = ""))))
  }

  if (diffexmethod == "DESeq"){
    diffex.results <- scDiffEx_deseq(inSCESet, in.condition)
  }
  else if (diffexmethod == "DESeq2"){
    diffex.results <- scDiffEx_deseq2(inSCESet, in.condition)
  }
  else if (diffexmethod == "limma"){
    diffex.results <- scDiffEx_limma(inSCESet, in.condition)
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
#' @export plot_DiffEx
#'
plot_DiffEx <- function(inSCESet, condition, geneList, clusterRow=TRUE,
                     clusterCol=TRUE, displayRowLabels=TRUE, displayColumnLabels=TRUE,
                     displayRowDendrograms=TRUE, displayColumnDendrograms=TRUE,
                     annotationColors=NULL, columnTitle="Differential Expression"){
  if (is.null(annotationColors)){
    topha <- NULL
  } else {
    topha <- ComplexHeatmap::HeatmapAnnotation(df = Biobase::pData(inSCESet)[, condition, drop = FALSE],
                                               col = annotationColors)
  }

  heatmap <- ComplexHeatmap::Heatmap(t(scale(t(Biobase::exprs(inSCESet)[geneList, ]))),
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

#' Plot Interactive Differential Expression
#'
#' @param inSCESet Input data object that contains the data to be plotted.
#' Required
#' @param condition The condition used for plotting the heatmap. Required
#' @param geneList The list of genes to put in the heatmap. Required
#' @param clusterRow Cluster the rows. The default is TRUE
#' @param clusterCol Cluster the columns. The default is TRUE
#'
#' @return A d3heatmap object is plotted
#' @export plot_d3DiffEx
#'
plot_d3DiffEx <- function(inSCESet, condition, geneList, clusterRow=TRUE,
                          clusterCol=TRUE){
  diffex.annotation <- data.frame(Biobase::pData(inSCESet)[, condition])
  colnames(diffex.annotation) <- condition
  topha <- ComplexHeatmap::HeatmapAnnotation(df = diffex.annotation,
                                             height = unit(0.333, "cm"))

  d3heatmap::d3heatmap(t(scale(t(Biobase::exprs(inSCESet)[geneList, ]))),
                       Rowv = clusterRow,
                       Colv = clusterCol,
                       ColSideColors = RColorBrewer::brewer.pal(8, "Set1")[as.numeric(factor(diffex.annotation[, 1]))])
}

#' Perform differential expression analysis with DESeq2
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCESet object. Required
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from pData that contains two labels.
#' Required
#'
#' @return A data frame of gene names and adjusted p-values
#' @export scDiffEx_deseq2
#'
scDiffEx_deseq2 <- function(inSCESet, condition){
  cnts <- scater::counts(inSCESet)
  annot_data <- data.frame(condition)
  colnames(annot_data) <- "condition"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts,
                                        colData = annot_data,
                                        design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  return(data.frame(res))
}

#' Perform differential expression analysis with DESeq
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCESet object. Required
#' @param condition A factor for the condition to use for differential
#' expression. Must be a two level factor. Required
#'
#' @return A data frame of gene names and adjusted p-values
#' @export scDiffEx_deseq
#'
scDiffEx_deseq <- function(inSCESet, condition){
  countData <- DESeq::newCountDataSet(scater::counts(inSCESet), condition)
  countData <- DESeq::estimateSizeFactors(countData)
  countData <- DESeq::estimateDispersions(countData, method = "pooled",
                                          fitType = "local")
  diff.results <- DESeq::nbinomTest(countData, levels(condition)[1],
                                    levels(condition)[2])
  top.results <- stats::p.adjust(diff.results$pval, method = "fdr")
  diff.results$padj <- top.results
  rownames(diff.results) <- diff.results$id
  diff.results$id <- NULL
  return(diff.results)
}

#' Perform differential expression analysis with limma
#'
#' Returns a data frame of gene names and adjusted p-values
#'
#' @param inSCESet Input SCESet object. Required
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from pData that contains two labels.
#' Required
#'
#' @return A data frame of gene names and adjusted p-values
#' @export scDiffEx_limma
#'
scDiffEx_limma <- function(inSCESet, condition){
  design <- stats::model.matrix(~factor(condition))
  fit <- limma::lmFit(Biobase::exprs(inSCESet), design)
  ebayes <- limma::eBayes(fit)
  topGenes <- limma::topTable(ebayes, coef = 2, adjust = "fdr", number = nrow(inSCESet))
  colnames(topGenes)[5] <- "padj"
  return(topGenes)
}
