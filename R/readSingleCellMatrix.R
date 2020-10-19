
#' @importFrom tools file_ext
.getFileExt <- function(file) {
    ext1 <- tools::file_ext(file)
    if (!(ext1 %in% c("mtx", "txt", "csv", "tab", "tsv",
        "npz", "gz", "bz2", "xz"))) {

        warning("Unknown extension ", ext1, ". Treat as text file.")
        return(c("unknown"))
    }

    # if compressed
    if (ext1 %in% c("gz", "bz2", "xz")) {
        bn <- substr(file, start = 1, stop = nchar(file) - nchar(ext1) - 1)
        ext2 <- tools::file_ext(bn)
        if (!(ext2 %in% c("mtx", "txt", "csv", "tab", "npz"))) {
            warning("Unknown extension ", ext2, ". Treat as text file.")
            return(c(ext1, "unknown"))
        }
        return(c(ext1, ext2))
    } else {
        return(ext1)
    }
}


#' @name readSingleCellMatrix
#' @title Read single cell expression matrix
#' @description Automatically detact the format of the input file and read
#'  the file.
#' @param file Path to input file. Supported file endings include .mtx, .txt,
#'  .csv, .tab, .tsv, .npz, and their corresponding \code{gzip},
#'  \code{bzip2}, or \code{xz} compressed extensions (\code{*.gz},
#'  \code{*.bz2}, or \code{*.xz}).
#' @param class Character. Class of matrix. One of "Matrix" or "matrix". Specifying "Matrix"
#'  will convert to a sparse format which should be used
#'  for datasets with large numbers of cells.  Default "Matrix".
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @param colIndexLocation Character. For Optimus output, the path to the
#'  barcode index .npy file. Used only if \code{file} has .npz extension.
#'  Default \code{NULL}.
#' @param rowIndexLocation Character. For Optimus output, The path to the
#'  feature (gene) index .npy file. Used only if \code{file} has .npz extension.
#'  Default \code{NULL}.
#' @return A \link[DelayedArray]{DelayedArray} object or matrix.
#' @examples
#' mat <- readSingleCellMatrix(system.file("extdata/hgmm_1k_v3_20x20/outs/",
#'     "filtered_feature_bc_matrix/matrix.mtx.gz", package = "singleCellTK"))
#' @importFrom reticulate import
#' @export
readSingleCellMatrix <- function(file,
    class = c("Matrix", "matrix"),
    delayedArray = TRUE,
    colIndexLocation = NULL,
    rowIndexLocation = NULL) {

    class <- match.arg(class)
    ext <- .getFileExt(file)

    if (ext[1] == "gz") {
        file <- gzfile(file)
    } else if (ext[1] == "bz2") {
        file <- bzfile(file)
    } else if (ext[1] == "xz") {
        file <- xzfile(file)
    }

    ext2 <- data.table::last(ext)

    if (ext2 %in% c("txt", "csv", "tab", "tsv", "unknown")) {
        dt <- data.table::fread(file)
        if (class == "Matrix") {
            mat <- Matrix::Matrix(dt[, -1])
            rownames(mat) <- dt[[1]]
        } else if (class == "matrix") {
            mat <- as.matrix(dt[, -1])
            rownames(mat) <- dt[[1]]
        }
    } else if (ext2 == "npz") {
        mat <- .readMatrixNpz(matrixLocation = file,
            colIndexLocation = colIndexLocation,
            rowIndexLocation = rowIndexLocation,
            class = class)
    } else if (ext2 == "mtx") {
        mat <- .readMatrixMM(path = file,
            gzipped = FALSE,
            class = class,
            delayedArray = delayedArray)
    }
    return(mat)
}
