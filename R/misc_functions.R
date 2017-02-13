#' Summarize SCESet
#' 
#' Creates a table of summary metrics from an input SCESet.
#'
#' @param indata Input SCESet
#'
#' @return A data.frame object of summary metrics.
#' @export summarizeTable
summarizeTable <- function(indata){
  return(data.frame("Metric"=c("Number of Samples",
                               "Number of Genes",
                               "Samples with <1700 detected genes",
                               "Genes with no expression across all samples"),
                    "Value"=c(ncol(indata),
                              nrow(indata),
                              sum(apply(scater::counts(indata), 2, function(x) sum(as.numeric(x)==0)) < 1700),
                              sum(rowSums(scater::counts(indata)) == 0))))
}

#' Create a SCESet object
#' 
#' From a file of counts and a file of annotation information, create a SCESet
#' object.
#'
#' @param countfile The path to a text file that contains a header row of sample
#' names, and rows of raw counts per gene for those samples.
#' @param annotfile The path to a text file that contains columns of annotation
#' information for each sample in the countfile. This file should have the same
#' number of rows as there are columns in the countfile.
#' @param featurefile The path to a text file that contains columns of
#' annotation information for each gene in the count matrix. This file should
#' have the same genes in the same order as countfile. This is optional.
#' @param inputdataframes If TRUE, countfile and annotfile are read as data
#' frames instead of file paths. The default is FALSE.
#' instead of 
#'
#' @return a SCESet object
#' @export createSCESet
createSCESet <- function(countfile=NULL, annotfile=NULL, featurefile=NULL,
                         inputdataframes=FALSE){
  if(is.null(countfile)){
    stop("You must supply a count file.")
  }
  if(inputdataframes){
    countsin <- countfile
    annotin <- annotfile
    featurein <- featurefile
  } else{
    countsin <- utils::read.table(countfile, sep="\t", header=T, row.names=1)
    if(!is.null(annotfile)){
      annotin <- utils::read.table(annotfile, sep="\t", header=T, row.names=1)
    }
    if(!is.null(featurefile)){
      featurein <- utils::read.table(featurefile, sep="\t", header=T, row.names=1)
    }
  }
  if(is.null(annotfile)){
    annotin <- data.frame(row.names=colnames(countsin))
    annotin$Sample <- rownames(annotin)
  }
  if(is.null(featurefile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
  }
  pd <- methods::new("AnnotatedDataFrame", data = annotin)
  fd <- methods::new("AnnotatedDataFrame", data = featurein)
  return(scater::newSCESet(countData = countsin, phenoData = pd,
                           featureData = fd))
}
