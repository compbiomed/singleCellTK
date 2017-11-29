#' Summarize SCtkExperiment
#'
#' Creates a table of summary metrics from an input SCtkExperiment.
#'
#' @param indata Input SCtkExperiment
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
                              as.integer(mean(colSums(assay(indata, "counts")))),
                              as.integer(mean(colSums(assay(indata, "counts") > 0))),
                              sum(colSums(assay(indata, "counts") != 0) < 1700),
                              sum(rowSums(assay(indata, "counts")) == 0))))
}

#' Create a SCtkExperiment object
#'
#' From a file of counts and a file of annotation information, create a
#' SCtkExperiment object.
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
#' @param create_logcounts If TRUE, create a log2(counts+1) normalized assay
#' and include it in the object. The default is TRUE
#'
#' @return a SCtkExperiment object
#' @export createSCE
#' @examples
#' \dontrun{
#' GSE60361_sce <- createSCE(countfile = "/path/to/input_counts.txt",
#'                           annotfile = "/path/to/input_annots.txt")
#'}
createSCE <- function(countfile=NULL, annotfile=NULL, featurefile=NULL,
                      inputdataframes=FALSE, create_logcounts=TRUE){
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
    annotin <- DataFrame(annotin)
  }
  if (is.null(featurefile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
    featurein <- DataFrame(featurein)
  }
  newassay <- SCtkExperiment(assays = list(counts = as.matrix(countsin)),
                                     colData = annotin,
                                     rowData = featurein)
  if (create_logcounts){
    assay(newassay, "logcounts") <- log2(assay(newassay, "counts") + 1)
  }
  return(newassay)
}

#' Filter Genes and Samples from a Single Cell Object
#'
#' @param insceset Input single cell object, required
#' @param use_assay Indicate which assay to use for filtering. Default is
#' "counts"
#' @param deletesamples List of samples to delete from the object.
#' @param remove_noexpress Remove genes that have no expression across all
#' samples. The default is true
#' @param remove_bottom Fraction of low expression genes to remove from the
#' single cell object. This occurs after remove_noexpress. The default is 0.50.
#' @param minimum_detect_genes Minimum number of genes with at least 1
#' count to include a sample in the single cell object. The default is 1700.
#' @param filter_spike Apply filtering to Spike in controls (indicated by
#' isSpike).
#' The default is TRUE.
#'
#' @return The filtered single cell object.
#' @export filterSCData
#'
#' @examples
#' data("GSE60361_subset_sce")
#' GSE60361_subset_sce <- filterSCData(GSE60361_subset_sce,
#'                                     deletesamples="X1772063061_G11")
filterSCData <- function(insceset, use_assay="counts", deletesamples=NULL,
                         remove_noexpress=TRUE, remove_bottom=0.5,
                         minimum_detect_genes=1700, filter_spike=TRUE){
  #delete specified samples
  insceset <- insceset[, !(colnames(insceset) %in% deletesamples)]

  if (filter_spike){
    nkeeprows <- ceiling((1 - remove_bottom) * as.numeric(nrow(insceset)))
    tokeeprow <- order(rowSums(assay(insceset, use_assay)), decreasing = TRUE)[1:nkeeprows]
  } else {
    nkeeprows <- ceiling((1 - remove_bottom) * as.numeric(nrow(insceset))) - sum(isSpike(insceset))
    tokeeprow <- order(rowSums(assay(insceset, use_assay)), decreasing = TRUE)
    tokeeprow <- setdiff(tokeeprow, which(isSpike(insceset)))
    tokeeprow <- tokeeprow[1:nkeeprows]
    tokeeprow <- c(tokeeprow, which(isSpike(insceset)))
  }
  tokeepcol <- colSums(assay(insceset, use_assay) != 0) >= minimum_detect_genes
  insceset <- insceset[tokeeprow, tokeepcol]

  #remove genes with no expression
  if (remove_noexpress){
    if (filter_spike){
      insceset <- insceset[rowSums(assay(insceset, use_assay)) != 0, ]
    } else {
      insceset <- insceset[(rowSums(assay(insceset, use_assay)) != 0 | isSpike(insceset)), ]
    }
  }

  return(insceset)
}
