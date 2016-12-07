#' Perform differential expression on a SCESet object using DESeq
#'
#' @param inSCESet Input SCESet object. Required
#' @param condition The name of the condition to use for differential
#' expression. Must be a name of a column from pData that contains two labels.
#' @param significance FDR corrected significance cutoff for differentially
#' expressed genes
#'
#' @return A list of differentially expressed genes.
#' @export scDiffEx
#'
#' @examples
#' 
#' scDiffEx(newsceset_david, "age", "0.1")
scDiffEx <- function(inSCESet, condition, significance, ntop, usesig){
  in.condition <- as.factor(pData(inSCESet)[,condition])
  if (length(levels(in.condition)) != 2)
    stop("only two labels supported, ", condition, " has ",
         length(levels(in.condition)), " labels")
  countData <- DESeq::newCountDataSet( counts(inSCESet), in.condition )
  countData <- DESeq::estimateSizeFactors( countData )
  countData <- DESeq::estimateDispersions( countData, method="pooled",
                                           fitType="local" )
  diff.results <- DESeq::nbinomTest( countData, levels(in.condition)[1],
                                     levels(in.condition)[2])
  top.results <- stats::p.adjust( diff.results$pval, method="fdr" )
  
  if(usesig){
    if(length(which(top.results <= significance)) < ntop){
      newgenes <- rownames(fData(inSCESet))[which(top.results <= significance)]
    }
    else{
      newgenes <- rownames(fData(inSCESet))[order(top.results)[1:ntop]]
    }
  }
  else{
    newgenes <- rownames(fData(inSCESet))[order(top.results)[1:ntop]]
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
  return(heatmap)
}
