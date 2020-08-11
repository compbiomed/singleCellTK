.checkGzip <- function(path, gzipped){
    if (gzipped == "auto") {
      ext <- tools::file_ext(path)
      if (ext == "gz") {
            path <- gzfile(path)
        }
    } else if (isTRUE(gzipped)) {
        path <- gzfile(path)
    }

    return(path)
}

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
#' @param annotFileHeader Whether there's a header (colnames) in the cell annotation file. Default is FALSE  
#' @param annotFileRowName Which column is used as the rownames for the cell annotation file. Default is 1 (first column). 
#' @param annotFileSep Separater used for the cell annotation file. Default is "\\t". 
#' @param featureHeader Whether there's a header (colnames) in the feature annotation file. Default is FALSE  
#' @param featureRowName Which column is used as the rownames for the feature annotation file. Default is 1 (first column).
#' @param featureSep Separater used for the feature annotation file. Default is "\\t".
#' @param gzipped Whether the input file is gzipped. Default is "auto" and it will automatically detect whether the file is gzipped. Other options is TRUE or FALSE. 
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @return a SingleCellExperiment object
#' @export

importFromFiles <- function(assayFile, annotFile = NULL, featureFile = NULL,
                            assayName = "counts", inputDataFrames = FALSE,
                            class = c("Matrix", "matrix"), delayedArray = FALSE,
                            annotFileHeader = FALSE, annotFileRowName = 1, 
                            annotFileSep = "\t", featureHeader = FALSE,
                            featureRowName = 1, featureSep = "\t", gzipped = "auto"
                            ){
  
  class <- match.arg(class)
  
  if (inputDataFrames){
    countsin <- assayFile
    annotin <- annotFile
    featurein <- featureFile
  } else{
    countsin <- readSingleCellMatrix(assayFile, class = class, delayedArray = delayedArray)
    if (!is.null(annotFile)){
      annotFile <- .checkGzip(annotFile, gzipped = gzipped)
      annotin <- utils::read.table(annotFile, sep = annotFileSep, header = annotFileHeader,
                                   row.names = annotFileRowName, stringsAsFactors = FALSE)
    }
    if (!is.null(featureFile)){
      featureFile <- .checkGzip(featureFile, gzipped = gzipped)
      featurein <- utils::read.table(featureFile, sep = featureSep, header = featureHeader,
                                     row.names = featureRowName, stringsAsFactors = FALSE)
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
  if (is.null(rownames(countsin))){
    rownames(countsin) <- rownames(featurein)
  }
  if (is.null(colnames(countsin))){
    colnames(countsin) <- rownames(annotin)
  }
  #assaylist[[assayName]] <- methods::as(countsin, "dgCMatrix")
  assaylist[[assayName]] <- .convertToMatrix(countsin)

  newassay <- SingleCellExperiment::SingleCellExperiment(assays = assaylist,
                                                         colData = annotin,
                                                         rowData = featurein)
  
  if(is.null(newassay$sample)) {
    newassay$sample <- "sample"
  }
  
  return(newassay)
}
