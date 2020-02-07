.readBarcodes <- function(path,
    header = FALSE,
    colname = "cell_barcode",
    colClasses = "character") {

    res <- data.table::fread(path, header = header, colClasses = colClasses)
    if (ncol(res) == 1) {
        colnames(res) <- colname
    } else {
        warning("'barcodes' file contains >1 columns!",
            " The column names are kept as is.")
    }
    return(res)
}


.readFeatures <- function(path,
    header = FALSE,
    colnames = c("feature_ID", "feature_name", "feature_type"),
    colClasses = "character") {

    res <- data.table::fread(path, header = header)
    if (ncol(res) == 1) {
        colnames(res) <- colnames[1]
    } else if (ncol(res) == 2) {
        colnames(res) <- colnames[1:2]
    } else if (ncol(res) == 3) {
        colnames(res) <- colnames
    } else {
        warning("'features' file contains >3 columns!",
            " The column names are kept as is.")
    }
    return(res)
}

#' @importFrom tools file_ext
.readMatrixMM <- function(path, gzipped, class, delayedArray) {

    if (gzipped == "auto") {
        ext <- tools::file_ext(path)
        if (ext == "gz") {
            path <- gzfile(path)
        }
    } else if (isTRUE(gzipped)) {
        path <- gzfile(path)
    }

    mat <- Matrix::readMM(path)

    if (class == "matrix") {
        mat <- as.matrix(mat)
    }

    if (isTRUE(delayedArray)) {
        mat <- DelayedArray::DelayedArray(mat)
    }
    return(mat)
}


# dir <- "outs/filtered_feature_bc_matrix/"
.constructSCEFromCellRangerOutputs <- function(dir,
    sampleName,
    matrixFileName,
    featuresFileName,
    barcodesFileName,
    gzipped,
    class,
    delayedArray) {

    cb <- .readBarcodes(file.path(dir, barcodesFileName))
    fe <- .readFeatures(file.path(dir, featuresFileName))
    ma <- .readMatrixMM(file.path(dir, matrixFileName),
        gzipped = gzipped,
        class = class,
        delayedArray = delayedArray)

    coln <- paste(sampleName, cb[[1]], sep = "_")
    rownames(ma) <- fe[[1]]

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = ma))
    SummarizedExperiment::rowData(sce) <- fe
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cb,
        column_name = coln,
        sample = sampleName,
        row.names = coln)

    return(sce)
}


.getOutputFolderPath <- function(samplePath, cellRangerOuts) {
    path <- file.path(samplePath, cellRangerOuts)
    return(path)
}


.checkArgsImportCellRanger <- function(cellRangerDirs,
    sampleDirs,
    sampleNames,
    cellRangerOuts,
    matrixFileNames,
    featuresFileNames,
    barcodesFileNames,
    gzipped) {

    if (!gzipped %in% c("auto", TRUE, FALSE)) {
        stop("Invalid 'gzipped' argument! Should be one of 'auto',",
            " TRUE, or FALSE")
    }

    if (is.null(cellRangerDirs)) {
        if (is.null(sampleDirs)) {
            stop("'sampleDirs' can not be NULL if 'cellRangerDirs' is NULL!")
        }

        for (i in seq_along(sampleDirs)) {
            if (!dir.exists(sampleDirs[i])) {
                stop("Sample folder ", samples[i], " does not exist!")
            }
        }

        sampleLength <- length(sampleDirs)

        if (!is.null(sampleNames)) {
            if (length(sampleNames) != sampleLength) {
                stop("'sampleDirs' and 'sampleNames' have unequal lengths!")
            }
        }

        if (length(cellRangerOuts) != 1) {
            if (length(cellRangerOuts) != sampleLength) {
                stop("'sampleDirs' and 'cellRangerOuts' have unequal lengths!")
            }
        }

        if (length(matrixFileNames) != 1) {
            if (length(matrixFileNames) != sampleLength) {
                stop("'sampleDirs' and 'matrixFileNames' have unequal lengths!")
            }
        }

        if (length(featuresFileNames) != 1) {
            if (length(featuresFileNames) != sampleLength) {
                stop("'sampleDirs' and 'featuresFileNames'",
                    " have unequal lengths!")
            }
        }

        if (length(barcodesFileNames) != 1) {
            if (length(barcodesFileNames) != sampleLength) {
                stop("'sampleDirs' and 'barcodesFileNames'",
                    " have unequal lengths!")
            }
        }

        if (gzipped != "auto") {
            if (length(gzipped) != sampleLength) {
                stop("'sampleDirs' and 'gzipped' have unequal lengths!")
            }
        }

    } else {
        if (is.null(sampleDirs)) {
            for (i in seq_along(cellRangerDirs)) {
                if (length(list.dirs(cellRangerDirs[i],
                    recursive = FALSE)) == 0) {
                    warning("Invalid folder. Skipping cellRangerDir ",
                        cellRangerDirs[i])
                }
            }

            sampleLength <- length(unlist(lapply(cellRangerDirs,
                list.dirs, recursive = FALSE)))

            if (!is.null(sampleNames)) {
                if (sampleLength != length(sampleNames)) {
                    stop("The length of 'sampleNames' does not match length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }

            if (length(cellRangerOuts) != 1) {
                if (sampleLength != length(cellRangerOuts)) {
                    stop("The length of 'cellRangerOuts' does not match",
                        " length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }

            if (length(matrixFileNames) != 1) {
                if (sampleLength != length(matrixFileNames)) {
                    stop("The length of 'matrixFileNames' does not match",
                        " length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }

            if (length(featuresFileNames) != 1) {
                if (sampleLength != length(featuresFileNames)) {
                    stop("The length of 'featuresFileNames' does not match",
                        " length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }

            if (length(barcodesFileNames) != 1) {
                if (sampleLength != length(barcodesFileNames)) {
                    stop("The length of 'barcodesFileNames' does not match",
                        " length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }

            if (gzipped != "auto") {
                if (sampleLength != length(gzipped)) {
                    stop("The length of 'gzipped' does not match",
                        " length of",
                        " subdirectories in 'cellRangerDirs'!")
                }
            }
        } else {
            if (length(sampleDirs) != length(cellRangerDirs)) {
                stop("'sampleDirs' and 'cellRangerDirs' have unequal lengths!")
            } else {
                for (i in seq_along(cellRangerDirs)) {
                    paths <- file.path(cellRangerDirs[i], sampleDirs[[i]])
                    for (j in seq_along(paths)) {
                        if (!dir.exists(paths[j])) {
                            stop("Sample folder does not exist!\n",
                                paths[j])
                        }
                    }
                }
            }
            # analogous to length(unlist(sampleDirs))
            sampleLength <- sum(vapply(sampleDirs, length, integer(1)))

            if (!is.null(sampleNames)) {
                if (length(sampleNames) != sampleLength) {
                    stop("'sampleNames' and 'unlist(sampleDirs)' have unequal",
                        " lengths!")
                }
            }

            if (length(cellRangerOuts) != 1) {
                if (length(cellRangerOuts) != sampleLength) {
                    stop("'cellRangerOuts' and 'unlist(sampleDirs)'",
                        " have unequal lengths!")
                }
            }

            if (length(matrixFileNames) != 1) {
                if (length(matrixFileNames) != sampleLength) {
                    stop("'matrixFileNames' and 'unlist(sampleDirs)'",
                        " have unequal lengths!")
                }
            }

            if (length(featuresFileNames) != 1) {
                if (length(featuresFileNames) != sampleLength) {
                    stop("'featuresFileNames' and 'unlist(sampleDirs)'",
                        " have unequal lengths!")
                }
            }

            if (length(barcodesFileNames) != 1) {
                if (length(barcodesFileNames) != sampleLength) {
                    stop("'barcodesFileNames' and 'unlist(sampleDirs)'",
                        " have unequal lengths!")
                }
            }

            if (gzipped != "auto") {
                if (length(gzipped) != sampleLength) {
                    stop("'gzipped' and 'unlist(sampleDirs)'",
                        " have unequal lengths!")
                }
            }
        }
    }
}


.getSamplesPaths <- function(cellRangerDirs, samplesDirs) {
    if (is.null(cellRangerDirs)) {
        res <- samplesDirs
    } else {
        if (is.null(samplesDirs)) {
            res <- list.dirs(cellRangerDirs, recursive = FALSE)
        } else {
            res <- vector("list", length = length(cellRangerDirs))
            for (i in seq_along(cellRangerDirs)) {
                res[[i]] <- file.path(cellRangerDirs[i], samplesDirs[[i]])
            }
            res <- unlist(res)
        }
    }
    return(res)
}


.getSampleNames <- function(samplePaths) {
    res <- basename(samplePaths)
    return(res)
}


.getVectorized <- function(arg, len) {
    if (length(arg) != 1) {
        arg <- rep(arg, len)
    }
    return(arg)
}

# main function
.importCellRanger <- function(
    cellRangerDirs,
    sampleDirs,
    sampleNames,
    cellRangerOuts,
    matrixFileNames,
    featuresFileNames,
    barcodesFileNames,
    gzipped,
    class,
    delayedArray) {

    .checkArgsImportCellRanger(cellRangerDirs,
        sampleDirs,
        sampleNames,
        cellRangerOuts,
        matrixFileNames,
        featuresFileNames,
        barcodesFileNames,
        gzipped)

    samplePaths <- .getSamplesPaths(cellRangerDirs, sampleDirs)

    res <- vector("list", length = length(samplePaths))

    cellRangerOuts <- .getVectorized(cellRangerOuts, length(samplePaths))
    matrixFileNames <- .getVectorized(matrixFileNames, length(samplePaths))
    featuresFileNames <- .getVectorized(featuresFileNames, length(samplePaths))
    barcodesFileNames <- .getVectorized(barcodesFileNames, length(samplePaths))
    gzipped <- .getVectorized(gzipped, length(samplePaths))

    if (is.null(sampleNames)) {
        for (i in seq_along(samplePaths)) {
            dir <- .getOutputFolderPath(samplePaths[i], cellRangerOuts[i])
            scei <- .constructSCEFromCellRangerOutputs(dir,
                sample = .getSampleNames(samplePaths[i]),
                matrixFileName = matrixFileNames[i],
                featuresFileName = featuresFileNames[i],
                barcodesFileName = barcodesFileNames[i],
                gzipped = gzipped[i],
                class = class,
                delayedArray = delayedArray)
            res[[i]] <- scei
        }
    } else {
        for (i in seq_along(samplePaths)) {
            dir <- .getOutputFolderPath(samplePaths[i], cellRangerOuts[i])
            scei <- .constructSCEFromCellRangerOutputs(dir,
                sample = sampleNames[i],
                matrixFileName = matrixFileNames[i],
                featuresFileName = featuresFileNames[i],
                barcodesFileName = barcodesFileNames[i],
                gzipped = gzipped[i],
                class = class,
                delayedArray = delayedArray)
            res[[i]] <- scei
        }
    }

    sce <- do.call(BiocGenerics::cbind, res)
    return(sce)
}


#' @name importCellRanger
#' @rdname importCellRanger
#' @title Construct SCE object from Cell Ranger output
#' @description Read the filtered barcodes, features, and matrices for all
#'  samples from (preferably a single run of) Cell Ranger output. Import and
#'  combine them
#'  as one big \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param cellRangerDirs The root directories where Cell Ranger was run. These
#'  folders should contain sample specific folders. Default \code{NULL},
#'  meaning the paths for each sample will be specified in \emph{samples}
#'  argument.
#' @param sampleDirs Default \code{NULL}. Can be one of
#' \itemize{
#'   \item \code{NULL}. All samples within \code{cellRangerDirs} will be
#'    imported. The order of samples will be first determined by the order of
#'    \code{cellRangerDirs} and then by \link[base]{list.dirs}. This is only
#'    for the case where \code{cellRangerDirs} is specified.
#'   \item A list of vectors containing the folder names for samples to import.
#'    Each vector in
#'    the list corresponds to samples from one of \code{cellRangerDirs}.
#'    These names are the same as the folder names under \code{cellRangerDirs}.
#'    This is only for the case where \code{cellRangerDirs} is specified.
#'   \item A vector of folder paths for the samples to import. This is only for
#'    the case where \code{cellRangerDirs} is \code{NULL}.
#' }
#' The cells in the final SCE object will be ordered in the same order of
#'  \code{sampleDirs}.
#' @param sampleNames A vector of user-defined sample names for the samples
#'  to be
#'  imported. Must have the same length as \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}. Default
#'  \code{NULL}, in which case the folder names will be used as sample names.
#' @param cellRangerOuts Character vector. The intermediate
#'  paths to filtered or raw cell barcode, feature, and matrix files
#'  for each sample. Must have length 1 or the same length as
#'  \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}.
#'  Reference genome names might need to be
#'  appended for CellRanger version below 3.0.0 if reads were mapped to
#'  multiple genomes when running Cell Ranger pipeline.
#' @param matrixFileNames Character vector. Filenames for the Market Exchange
#'  Format (MEX) sparse matrix files (matrix.mtx or matrix.mtx.gz files).
#'  Must have length 1 or the same
#'  length as \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}.
#' @param featuresFileNames Character vector. Filenames for the feature
#'  annotation files. They are usually named \emph{features.tsv.gz} or
#'  \emph{genes.tsv}. Must have length 1 or the same
#'  length as \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}.
#' @param barcodesFileNames Character vector. Filename for the cell barcode
#'  list files. They are usually named \emph{barcodes.tsv.gz} or
#'  \emph{barcodes.tsv}. Must have length 1 or the same
#'  length as \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}.
#' @param gzipped \code{TRUE} if the Cell Ranger output files
#'  (barcodes.tsv, features.tsv, and matrix.mtx) were
#'  gzip compressed. \code{FALSE} otherwise. This is true after Cell Ranger
#'  3.0.0 update. Default \code{"auto"} which automatically detects if the
#'  files are gzip compressed. If not \code{"auto"}, \code{gzipped} must have
#'  length 1 or the same
#'  length as \code{length(unlist(sampleDirs))} if
#'  \code{sampleDirs} is not \code{NULL}. Otherwise, make sure the length and
#'  order match the output of
#'  \code{unlist(lapply(cellRangerDirs, list.dirs, recursive = FALSE))}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @details
#'  \code{importCellRangerV2} imports output from Cell Ranger V2.
#'  \code{importCellRangerV2Sample} imports output from one sample from Cell
#'  Ranger V2.
#'  \code{importCellRangerV3} imports output from Cell Ranger V3.
#'  \code{importCellRangerV3} imports output from one sample from Cell Ranger
#'  V3.
#'  Some implicit
#'  assumptions which match the output structure of Cell Ranger V2 & V3
#'  are made in these 4 functions including \code{cellRangerOuts},
#'  \code{matrixFileName}, \code{featuresFileName}, \code{barcodesFileName},
#'  and \code{gzipped}.
#'  Alternatively, user can call \code{importCellRanger} to explicitly
#'  specify these arguments.
#' @return A \code{SingleCellExperiment} object containing the combined count
#'  matrix, the feature annotations, and the cell annotation.
#' @examples
#' # Example #1
#' # The following filtered feature, cell, and matrix files were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/
#' # 3.0.0/hgmm_1k_v3
#' # The top 10 hg19 & mm10 genes are included in this example.
#' # Only the first 20 cells are included.
#' sce <- importCellRanger(
#'     cellRangerDirs = system.file("extdata/", package = "singleCellTK"),
#'     sampleDirs = "hgmm_1k_v3_20x20",
#'     sampleNames = "hgmm1kv3")
#' @export
importCellRanger <- function(
    cellRangerDirs = NULL,
    sampleDirs = NULL,
    sampleNames = NULL,
    cellRangerOuts = "outs/filtered_feature_bc_matrix/",
    matrixFileNames = "matrix.mtx.gz",
    featuresFileNames = "features.tsv.gz",
    barcodesFileNames = "barcodes.tsv.gz",
    gzipped = "auto",
    class = c("Matrix", "matrix"),
    delayedArray = TRUE) {

    class <- match.arg(class)

    .importCellRanger(cellRangerDirs = cellRangerDirs,
        sampleDirs = sampleDirs,
        sampleNames = sampleNames,
        cellRangerOuts = cellRangerOuts,
        matrixFileNames = matrixFileNames,
        featuresFileNames = featuresFileNames,
        barcodesFileNames = barcodesFileNames,
        gzipped = gzipped,
        class = class,
        delayedArray = delayedArray)
}


#' @rdname importCellRanger
#' @param dataType The type of output to import. Whether to import the
#'  filtered or the raw data. Can be one of "filtered" or
#'  "raw". Default "filtered". When \code{dataType} is specified,
#'  \code{cellRangerDirs} will be ignored.
#' @export
importCellRangerV2 <- function(
    cellRangerDirs = NULL,
    sampleDirs = NULL,
    sampleNames = NULL,
    dataType = c("filtered", "raw"),
    class = c("Matrix", "matrix"),
    delayedArray = TRUE) {

    class <- match.arg(class)
    dataType <- match.arg(dataType)

    if (dataType == "filtered") {
        .importCellRanger(cellRangerDirs = cellRangerDirs,
            sampleDirs = sampleDirs,
            sampleNames = sampleNames,
            cellRangerOuts = "outs/filtered_gene_bc_matrices/",
            matrixFileNames = "matrix.mtx",
            featuresFileNames = "genes.tsv",
            barcodesFileNames = "barcodes.tsv",
            gzipped = FALSE,
            class = class,
            delayedArray = delayedArray)
    } else if (dataType == "raw") {
        .importCellRanger(cellRangerDirs = cellRangerDirs,
            sampleDirs = sampleDirs,
            sampleNames = sampleNames,
            cellRangerOuts = "outs/raw_gene_bc_matrices/",
            matrixFileNames = "matrix.mtx",
            featuresFileNames = "genes.tsv",
            barcodesFileNames = "barcodes.tsv",
            gzipped = FALSE,
            class = class,
            delayedArray = delayedArray)
    }
}


#' @rdname importCellRanger
#' @export
importCellRangerV2Sample <- function(
    sampleDir = NULL,
    sampleName = NULL,
    class = c("Matrix", "matrix"),
    delayedArray = TRUE) {

    class <- match.arg(class)

    .importCellRanger(cellRangerDirs = NULL,
        sampleDirs = sampleDir,
        sampleNames = sampleName,
        cellRangerOuts = "",
        matrixFileNames = "matrix.mtx",
        featuresFileNames = "genes.tsv",
        barcodesFileNames = "barcodes.tsv",
        gzipped = FALSE,
        class = class,
        delayedArray = delayedArray)
}


#' @rdname importCellRanger
#' @examples
#' sce <- importCellRangerV3(
#'     cellRangerDirs = system.file("extdata", package = "singleCellTK"),
#'     sampleDirs = "hgmm_1k_v3_20x20",
#'     sampleNames = "hgmm1kv3")
#' @export
importCellRangerV3 <- function(
    cellRangerDirs = NULL,
    sampleDirs = NULL,
    sampleNames = NULL,
    dataType = c("filtered", "raw"),
    class = c("Matrix", "matrix"),
    delayedArray = delayedArray) {

    class <- match.arg(class)
    dataType <- match.arg(dataType)

    if (dataType == "filtered") {
        .importCellRanger(cellRangerDirs = cellRangerDirs,
            sampleDirs = sampleDirs,
            sampleNames = sampleNames,
            cellRangerOuts = "outs/filtered_feature_bc_matrix/",
            matrixFileNames = "matrix.mtx.gz",
            featuresFileNames = "features.tsv.gz",
            barcodesFileNames = "barcodes.tsv.gz",
            gzipped = TRUE,
            class = class,
            delayedArray = delayedArray)
    } else if (dataType == "raw") {
        .importCellRanger(cellRangerDirs = cellRangerDirs,
            sampleDirs = sampleDirs,
            sampleNames = sampleNames,
            cellRangerOuts = "outs/raw_feature_bc_matrix/",
            matrixFileNames = "matrix.mtx.gz",
            featuresFileNames = "features.tsv.gz",
            barcodesFileNames = "barcodes.tsv.gz",
            gzipped = TRUE,
            class = class,
            delayedArray = delayedArray)
    }
}


#' @rdname importCellRanger
#' @examples
#' sce <- importCellRangerV3Sample(
#'     sampleDir = system.file("extdata/hgmm_1k_v3_20x20/outs/",
#'         "filtered_feature_bc_matrix", package = "singleCellTK"),
#'     sampleName = "hgmm1kv3")
#' @export
importCellRangerV3Sample <- function(
    sampleDir = NULL,
    sampleName = NULL,
    class = c("Matrix", "matrix"),
    delayedArray = TRUE) {

    class <- match.arg(class)

    .importCellRanger(cellRangerDirs = NULL,
        sampleDirs = sampleDir,
        sampleNames = sampleName,
        cellRangerOuts = "",
        matrixFileNames = "matrix.mtx.gz",
        featuresFileNames = "features.tsv.gz",
        barcodesFileNames = "barcodes.tsv.gz",
        gzipped = TRUE,
        class = class,
        delayedArray = delayedArray)
}
