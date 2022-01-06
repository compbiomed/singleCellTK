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
#' @details Creates a \linkS4class{SingleCellExperiment} object from a counts 
#' file in various formats, and files of cell and feature annotation.
#' @param assayFile The path to a file in .mtx, .txt, .csv, .tab, or .tsv 
#' format.
#' @param annotFile The path to a text file that contains columns of annotation
#' information for each cell in the \code{assayFile}. This file should have the 
#' same number of rows as there are columns in the \code{assayFile}. If multiple
#' samples are represented in the dataset, this should be denoted by a column 
#' called \code{'sample'} within the \code{annotFile}.
#' @param featureFile The path to a text file that contains columns of
#' annotation information for each gene in the count matrix. This file should
#' have the same genes in the same order as \code{assayFile}. This is optional.
#' @param assayName The name of the assay that you are uploading. The default
#' is \code{"counts"}.
#' @param inputDataFrames If \code{TRUE}, \code{assayFile}, \code{annotFile} and
#' \code{featureFile} should be \code{data.frames} object (or its inheritance) 
#' instead of file paths. The default is \code{FALSE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of \code{"Matrix"} (as returned by
#'  \link{readMM} function), or \code{"matrix"} (as returned by
#'  \link[base]{matrix} function). Default \code{"Matrix"}.
#' @param annotFileHeader Whether there's a header (colnames) in the cell 
#' annotation file. Default is \code{FALSE}.
#' @param annotFileRowName Which column is used as the rownames for the cell 
#' annotation file. This should match to the colnames of the \code{assayFile}. 
#' Default is \code{1} (first column).
#' @param annotFileSep Separater used for the cell annotation file. Default is 
#' \code{"\\t"}.
#' @param featureHeader Whether there's a header (colnames) in the feature 
#' annotation file. Default is \code{FALSE}.
#' @param featureRowName Which column is used as the rownames for the feature 
#' annotation file. This should match to the rownames of the \code{assayFile}. 
#' Default is \code{1}. (first column).
#' @param featureSep Separater used for the feature annotation file. Default is 
#' \code{"\\t"}.
#' @param gzipped Whether the input file is gzipped. Default is \code{"auto"} 
#' and it will automatically detect whether the file is gzipped. Other options 
#' are \code{TRUE} or \code{FALSE}.
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link{DelayedArray} object or not. Default \code{FALSE}.
#' @param rowNamesDedup Boolean. Whether to deduplicate rownames. Default 
#'  \code{TRUE}.
#' @return a \linkS4class{SingleCellExperiment} object
#' @export

importFromFiles <- function(assayFile, annotFile = NULL, featureFile = NULL,
                            assayName = "counts", inputDataFrames = FALSE,
                            class = c("Matrix", "matrix"), delayedArray = FALSE,
                            annotFileHeader = FALSE, annotFileRowName = 1,
                            annotFileSep = "\t", featureHeader = FALSE,
                            featureRowName = 1, featureSep = "\t", 
                            gzipped = "auto", rowNamesDedup = TRUE){

  class <- match.arg(class)

  if (inputDataFrames){
    countsin <- assayFile
    annotin <- annotFile
    featurein <- featureFile
  } else{
    countsin <- readSingleCellMatrix(assayFile, class = class, 
                                     delayedArray = delayedArray)
    if (!is.null(annotFile)){
      annotFile <- .checkGzip(annotFile, gzipped = gzipped)
      annotin <- utils::read.table(annotFile, sep = annotFileSep, 
                                   header = annotFileHeader,
                                   row.names = annotFileRowName, 
                                   stringsAsFactors = FALSE)
    }
    if (!is.null(featureFile)){
      featureFile <- .checkGzip(featureFile, gzipped = gzipped)
      featurein <- utils::read.table(featureFile, sep = featureSep, 
                                     header = featureHeader,
                                     row.names = featureRowName, 
                                     stringsAsFactors = FALSE)
    }
  }
  if (is.null(annotFile)){
    annotin <- data.frame(row.names = colnames(countsin))
    annotin <- S4Vectors::DataFrame(annotin)
  }
  if (is.null(featureFile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
    featurein <- S4Vectors::DataFrame(featurein)
  }

  if (nrow(annotin) != ncol(countsin)){
    stop("Different number of cells in input matrix and annotations: annot: ",
         nrow(annotin), ", counts: ", ncol(countsin))
  }
  if (nrow(featurein) != nrow(countsin)){
    stop("Different number of features in input matrix and feature annotation",
         nrow(featurein), ", counts: ", nrow(countsin))
  }
  if (any(rownames(annotin) != colnames(countsin))){
    stop("Cell names in input matrix and annotation do not match!\nExample: ",
         rownames(annotin)[rownames(annotin) != colnames(countsin)][1], " vs. ",
         colnames(countsin)[rownames(annotin) != colnames(countsin)][1])
  }
  if (any(rownames(featurein) != rownames(countsin))){
    stop("Feature names in input matrix and feature annotation do not match!")
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

  if (isTRUE(rowNamesDedup)) {
    if (any(duplicated(rownames(newassay)))) {
      message("Duplicated gene names found, adding '-1', '-2', ",
              "... suffix to them.")
    }
    newassay <- dedupRowNames(newassay)
  }
  
  return(newassay)
}
