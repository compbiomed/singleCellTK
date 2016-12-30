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
                     usesig=TRUE, diffexmethod){
  in.condition <- as.factor(pData(inSCESet)[,condition])
  if (length(levels(in.condition)) != 2)
    stop("only two labels supported, ", condition, " has ",
         length(levels(in.condition)), " labels")
  if(diffexmethod == "DESeq"){
    diffex.results <- scDiffEx_deseq(inSCESet, in.condition)
  }
  else if(diffexmethod == "DESeq2"){
    diffex.results <- scDiffEx_deseq2(inSCESet, condition)
  }
  else if(diffexmethod == "limma"){
    diffex.results <- scDiffEx_limma(inSCESet, condition)
  }
  else{
    stop("Unsupported differential expression method, ", diffexmethod)
  }
  
  if(usesig){
    if(length(which(diffex.results$padj <= significance)) < ntop){
      newgenes <- rownames(diffex.results)[which(diffex.results$padj <= significance)]
    }
    else{
      newgenes <- rownames(diffex.results)[order(diffex.results$padj)[1:ntop]]
    }
  }
  else{
    newgenes <- rownames(diffex.results)[order(diffex.results$padj)[1:ntop]]
  }
  
  diffex.annotation <- data.frame(pData(inSCESet)[,condition])
  colnames(diffex.annotation) <- condition
  topha <- ComplexHeatmap::HeatmapAnnotation(df = diffex.annotation,
                                             height = unit(0.333, "cm"))
  
  heatmap <- ComplexHeatmap::Heatmap(t(scale(t(exprs(inSCESet)[newgenes,]))),
                                     name="Expression",
                                     column_title = "Differential Expression",
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     top_annotation = topha)
  if(FALSE){
    scDiffEx_deseq2(inSCESet, in.condition)
  }
  return(heatmap)
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
  cnts <- counts(inSCESet)
  annot_data <- pData(inSCESet)[,condition,drop=F]
  colnames(annot_data) <- "condition"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts,
                                        colData = annot_data,
                                        design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  return(data.frame(res[,"padj",drop=F]))
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
  countData <- DESeq::newCountDataSet(counts(inSCESet), condition)
  countData <- DESeq::estimateSizeFactors(countData)
  countData <- DESeq::estimateDispersions(countData, method="pooled",
                                          fitType="local")
  diff.results <- DESeq::nbinomTest( countData, levels(condition)[1],
                                     levels(condition)[2])
  top.results <- stats::p.adjust( diff.results$pval, method="fdr" )
  return(data.frame(row.names = diff.results$id, padj=top.results))
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
  design <- model.matrix(~factor(pData(inSCESet)[,condition]))
  fit <- lmFit(exprs(inSCESet), design)
  ebayes <- eBayes(fit)
  topGenes <- topTable(ebayes, coef=2, adjust="fdr", number=nrow(inSCESet))
  topGenes <- topGenes[,"adj.P.Val", drop=F]
  colnames(topGenes) <- "padj"
  return(topGenes)
}
