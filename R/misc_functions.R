#' Summarize SCESet
#' 
#' Creates a table of summary metrics from an input SCESet.
#'
#' @param indata Input SCESet
#'
#' @return A data.frame object of summary metrics.
#' @export summarizeTable
summarizeTable <- function(indata){
  return(data.frame("Metric"=c("Number of Samples","Number of Genes",
                               "Samples with <1700 detected genes"),
                    "Value"=c(ncol(indata),
                              nrow(indata),
                              sum(apply(indata, 2, function(x) sum(as.numeric(x)==0)) < 1700))))
}

#' Create a SCESet object
#' 
#' From a file of counts and a file of annotation information, create a SCESet
#' object.
#'
#' @param countfile The path to a text file that contains a header row of sample
#' names, and rows of raw counts per gene for those samples
#' @param annotfile The path to a text file that contains columns of annotation
#' information for each sample in the countfile. This file should have the same
#' number of rows as there are columns in the countfile.
#'
#' @return a SCESet object
#' @export createSCESet
createSCESet <- function(countfile, annotfile){
  countsin <- read.table(countfile, sep="\t", header=T, row.names=1)
  annotin <- read.table(annotfile, sep="\t", header=T, row.names=1)
  pd <- new("AnnotatedDataFrame", data = annotin)
  
  gene_df <- data.frame(Gene = rownames(countsin))
  rownames(gene_df) <- gene_df$Gene
  fd <- new("AnnotatedDataFrame", data = gene_df)
  return(newSCESet(countData = countsin, phenoData = pd,
                   featureData = fd))
}