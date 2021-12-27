#' @name importAlevin
#' @rdname importAlevin
#' @title Construct SCE object from Salmon-Alevin output
#' @param alevinDir Character. The output directory of salmon-Alevin pipeline.
#'  It should contain subfolder named 'alevin', which contains the count data
#'  which is stored
#'  in 'quants_mat.gz'. Default \code{NULL}.
#' @param sampleName Character. A user-defined sample name for the sample to be
#'  imported. The 'sampleName' will be appended to the begining of cell
#'  barcodes. Default is 'sample'.
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link{DelayedArray} object or not. Default \code{FALSE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "Matrix" (as returned by
#'  \link{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @param rowNamesDedup Boolean. Whether to deduplicate rownames. Default 
#'  \code{TRUE}.
#' @return A \code{SingleCellExperiment} object containing the count
#'  matrix, the feature annotations, and the cell annotation
#'  (which includes QC metrics stored in 'featureDump.txt').
#' @import fishpond
#' @export

importAlevin <- function(
  alevinDir = NULL,
  sampleName = 'sample',
  delayedArray = FALSE,
  class = c("Matrix", "matrix"),
  rowNamesDedup = TRUE) {

  class <- match.arg(class)

  matFile <- file.path(alevinDir, "alevin/quants_mat.gz")
  ### require package fishpond
  ma <- tximport::tximport(files = matFile, type = "alevin")
  if (!'counts' %in% names(ma)) {
    stop("RNA count matrix not found in the alevin output!")
  }
  mat <- ma$counts

  if (class == "Matrix") {
    mat <- .convertToMatrix(mat)
  } else if (class == "matrix") {
    mat <- base::as.matrix(mat)
  }
  
  if (delayedArray) {
    mat <- DelayedArray::DelayedArray(mat)
  }
  
  if (isTRUE(rowNamesDedup)) {
    if (any(duplicated(rownames(mat)))) {
      message("Duplicated gene names found, adding '-1', '-2', ",
              "... suffix to them.")
    }
    mat <- dedupRowNames(mat)
  }
  
  genes <- rownames(mat)
  cb <- .readBarcodes(file.path(alevinDir, 'alevin/featureDump.txt'),
                      header = 'auto',
                      colname = "cell_barcode",
                      colClasses = "character")
  coln <- paste(sampleName, cb[[1]], sep = "_")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat))
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(
    'feature_name' = genes,
    row.names = genes)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
    cb,
    column_name = coln,
    sample = sampleName,
    row.names = coln)
  
  return(sce)
}
