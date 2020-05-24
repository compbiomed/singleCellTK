#' Create a SingleCellExperiment object from files
#'
#' Creates a SingleCellExperiment object from a counts file in various formats.
#' and a file of annotation information, .
#'
#' @param assayFile The path to a  file in .mtx, .txt, .csv, .tab, or .tsv format.
#' @param annotFile The path to a text file that contains columns of annotation
#' information for each sample in the assayFile. This file should have the same
#' number of rows as there are columns in the assayFile. If multiple samples are
#' represented in these files, this should be denoted by a column called \code{'sample'}
#' within the \code{annotFile}.
#' @param featureFile The path to a text file that contains columns of
#' annotation information for each gene in the count matrix. This file should
#' have the same genes in the same order as assayFile. This is optional.
#' @param assayName The name of the assay that you are uploading. The default
#' is "counts".
#' @param inputDataFrames If TRUE, assayFile and annotFile are read as data
#' frames instead of file paths. The default is FALSE.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @return a SingleCellExperiment object
#' @export
importFromFiles <- function(assayFile, annotFile = NULL, featureFile = NULL,
                            assayName = "counts", inputDataFrames = FALSE,
                            class = c("Matrix", "matrix"), delayedArray = FALSE){
  
  if (inputDataFrames){
    countsin <- assayFile
    annotin <- annotFile
    featurein <- featureFile
  } else{
    countsin <- readSingleCellMatrix(assayFile, class = class, delayedArray = delayedArray)
    if (!is.null(annotFile)){
      annotin <- utils::read.table(annotFile, sep = "\t", header = TRUE,
                                   row.names = 1)
    }
    if (!is.null(featureFile)){
      featurein <- utils::read.table(featureFile, sep = "\t", header = TRUE,
                                     row.names = 1)
    }
  }
  if (is.null(annotFile)){
    annotin <- data.frame(row.names = colnames(countsin))
    annotin$Sample <- rownames(annotin)
    annotin <- S4Vectors::DataFrame(annotin)
  }
  if (is.null(featureFile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
    featurein <- S4Vectors::DataFrame(featurein)
  }
  if (nrow(annotin) != ncol(countsin)){
    stop("Different number of samples in input matrix and annotations: annot: ",
         nrow(annotin), ", counts: ", ncol(countsin))
  }
  if (nrow(featurein) != nrow(countsin)){
    stop("Different number of samples in input matrix and feature annotation",
         nrow(featurein), ", counts: ", nrow(countsin))
  }
  if (any(rownames(annotin) != colnames(countsin))){
    stop("Sample names in input matrix and annotation do not match!\nExample: ",
         rownames(annotin)[rownames(annotin) != colnames(countsin)][1], " vs. ",
         colnames(countsin)[rownames(annotin) != colnames(countsin)][1])
  }
  if (any(rownames(featurein) != rownames(countsin))){
    stop("Sample names in input matrix and feature annotation do not match!")
  }
  assaylist <- list()
  assaylist[[assayName]] <- methods::as(countsin, "dgCMatrix")
  newassay <- SingleCellExperiment::SingleCellExperiment(assays = assaylist,
                                                         colData = annotin,
                                                         rowData = featurein)
  
  if(is.null(newassay$sample)) {
    newassay$sample <- "sample"
  }
  
  return(newassay)
}
