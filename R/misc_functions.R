#' Summarize SCESet
#'
#' Creates a table of summary metrics from an input SCESet.
#'
#' @param indata Input SCESet
#'
#' @return A data.frame object of summary metrics.
#' @export summarizeTable
summarizeTable <- function(indata){
  return(data.frame("Metric" = c("Number of Samples",
                               "Number of Genes",
                               "Average number of reads per cell",
                               "Average number of genes per cell",
                               "Samples with <1700 detected genes",
                               "Genes with no expression across all samples"),
                    "Value" = c(ncol(indata),
                              nrow(indata),
                              as.integer(mean(apply(scater::counts(indata), 2, function(x) sum(x)))),
                              as.integer(mean(apply(scater::counts(indata), 2, function(x) sum(x > 0)))),
                              sum(apply(scater::counts(indata), 2, function(x) sum(as.numeric(x) == 0)) < 1700),
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
#'
#' @return a SCESet object
#' @export createSCESet
#' @examples
#' library(scater)
#' data("GSE60361_subset")
#' GSE60361_SCESet <- createSCESet(countfile = GSE60361_subset$counts,
#'                                 annotfile = GSE60361_subset$annot,
#'                                 inputdataframes = TRUE)
createSCESet <- function(countfile=NULL, annotfile=NULL, featurefile=NULL,
                         inputdataframes=FALSE){
  if (is.null(countfile)){
    stop("You must supply a count file.")
  }
  if (inputdataframes){
    countsin <- countfile
    annotin <- annotfile
    featurein <- featurefile
  } else{
    countsin <- utils::read.table(countfile, sep = "\t", header = T, row.names = 1)
    if (!is.null(annotfile)){
      annotin <- utils::read.table(annotfile, sep = "\t", header = T, row.names = 1)
    }
    if (!is.null(featurefile)){
      featurein <- utils::read.table(featurefile, sep = "\t", header = T, row.names = 1)
    }
  }
  if (is.null(annotfile)){
    annotin <- data.frame(row.names = colnames(countsin))
    annotin$Sample <- rownames(annotin)
  }
  if (is.null(featurefile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
  }
  pd <- methods::new("AnnotatedDataFrame", data = annotin)
  fd <- methods::new("AnnotatedDataFrame", data = featurein)
  return(scater::newSCESet(countData = countsin, phenoData = pd,
                           featureData = fd))
}

#' Filter Genes and Samples from a Single Cell Object
#'
#' @param insceset Input single cell object, required
#' @param deletesamples List of samples to delete from the object.
#' @param remove_noexpress Remove genes that have no expression across all
#' samples. The default is true
#' @param remove_bottom Fraction of low expression genes to remove from the
#' single cell object. This occurs after remove_noexpress. The default is 0.50.
#' @param minimum_detect_genes Minimum number of genes with at least 1
#' count to include a sample in the single cell object. The default is 1700.
#'
#' @return The filtered single cell object.
#' @export filterSCData
#'
#' @examples
#' library(scater)
#' data("GSE60361_subset")
#' GSE60361_SCESet <- createSCESet(countfile = GSE60361_subset$counts,
#'                                 annotfile = GSE60361_subset$annot,
#'                                 inputdataframes = TRUE)
#' GSE60361_SCESet_filtered <- filterSCData(GSE60361_SCESet,
#'                                          deletesamples="X1772063061_G11")
filterSCData <- function(insceset, deletesamples=NULL, remove_noexpress=TRUE,
                         remove_bottom=0.5, minimum_detect_genes=1700){
  insceset <- insceset[, !(colnames(insceset) %in% deletesamples)]
  if (remove_noexpress){
    insceset <- insceset[rowSums(counts(insceset)) != 0, ]
  }
  nkeeprows <- ceiling((1 - remove_bottom) * as.numeric(nrow(insceset)))
  tokeeprow <- order(rowSums(counts(insceset)), decreasing = TRUE)[1:nkeeprows]
  tokeepcol <- apply(counts(insceset), 2, function(x) sum(as.numeric(x) == 0)) >= minimum_detect_genes
  insceset <- insceset[tokeeprow, tokeepcol]
  return(insceset)
}
