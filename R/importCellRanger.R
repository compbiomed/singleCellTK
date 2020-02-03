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


.readMatrixMM <- function(path, gzipped, class) {
    if (isTRUE(gzipped)) {
        path <- gzfile(path)
    }

    res <- Matrix::readMM(path)

    if (class == "Matrix") {
        return(res)
    } else if (class == "DelayedArray") {
        res <- DelayedArray::DelayedArray(res)
        return(res)
    } else if (class == "matrix") {
        res <- as.matrix(res)
        return(res)
    }
}


# dir <- "outs/filtered_feature_bc_matrix/"
.constructSCEFromCellRangerOutputs <- function(dir,
    sample,
    matrixFileName,
    featuresFileName,
    barcodesFileName,
    gzipped,
    class) {

    cb <- .readBarcodes(file.path(dir, barcodesFileName))
    fe <- .readFeatures(file.path(dir, featuresFileName))
    ma <- .readMatrixMM(file.path(dir, matrixFileName),
        gzipped = gzipped,
        class = class)

    coln <- paste(sample, cb[[1]], sep = "_")
    rownames(ma) <- fe[[1]]

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = ma))
    SummarizedExperiment::rowData(sce) <- fe
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cb,
        column_name = coln,
        sample = sample,
        row.names = coln)

    return(sce)
}


.getOutputFolderPath <- function(sample, cellRangerOuts) {
    path <- file.path(sample, cellRangerOuts)
    return(path)
}


.checkArgsImportCellRanger <- function(cellRangerDirs, samples, class) {
    if (is.null(cellRangerDirs)) {
        if (is.null(samples)) {
            stop("samples can not be NULL if cellRangerDirs is NULL!")
        }
        for (i in seq_along(samples)) {
            if (!dir.exists(samples[i])) {
                stop("Sample folder does not exist!\n", samples[i])
            }
        }
    } else {
        if (is.null(samples)) {
            for (i in seq_along(cellRangerDirs)) {
                if (length(list.dirs(cellRangerDirs[i],
                    recursive = FALSE)) == 0) {
                    warning("Empty folder. Skipping cellRangerDir ",
                        cellRangerDirs[i])
                }
            }
        } else {
            if (!(length(samples) == length(cellRangerDirs))) {
                stop("Length of samples is not equal to length of ",
                    "cellRangerDirs!")
            } else {
                for (i in seq_along(cellRangerDirs)) {
                    paths <- file.path(cellRangerDirs[i], samples[[i]])
                    for (j in seq_along(paths)) {
                        if (!dir.exists(paths[j])) {
                            stop("Sample folder does not exist!\n",
                                paths[j])
                        }
                    }
                }
            }
        }
    }

    if (!(class %in% c("DelayedArray", "Matrix", "matrix"))) {
        stop("Invalid 'class' argument!")
    }
}


.getSamplesPaths <- function(cellRangerDirs, samples) {
    if (is.null(cellRangerDirs)) {
        res <- samples
    } else {
        if (is.null(samples)) {
            res <- list.dirs(cellRangerDirs, recursive = FALSE)
        } else {
            res <- vector("list", length = length(cellRangerDirs))
            for (i in seq_along(cellRangerDirs)) {
                res[[i]] <- file.path(cellRangerDirs[i], samples[[i]])
            }
            res <- unlist(res)
        }
    }
    return(res)
}


.getSampleNames <- function(samples) {
    res <- basename(samples)
    return(res)
}


# main function
.importCellRanger <- function(
    cellRangerDirs,
    samples,
    cellRangerOuts,
    matrixFileName,
    featuresFileName,
    barcodesFileName,
    gzipped,
    class) {

    .checkArgsImportCellRanger(cellRangerDirs, samples, class)
    samplePaths <- .getSamplesPaths(cellRangerDirs, samples)

    res <- vector("list", length = length(samples))

    for (i in seq_along(samples)) {
        dir <- .getOutputFolderPath(samplePaths[i], cellRangerOuts)
        scei <- .constructSCEFromCellRangerOutputs(dir,
            sample = .getSampleNames(samples[i]),
            matrixFileName = matrixFileName,
            featuresFileName = featuresFileName,
            barcodesFileName = barcodesFileName,
            gzipped = gzipped,
            class = class)
        res[[i]] <- scei
    }

    sce <- do.call(SingleCellExperiment::cbind, res)
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
#' @param samples Default \code{NULL}. Can be one of
#' \itemize{
#'   \item \code{NULL}. All samples within \emph{cellRangerDirs} will be
#'    imported. The order of samples will be first determined by the order of
#'    \emph{cellRangerDirs} and then by \link[base]{list.files}. This is only
#'    for the case where \emph{cellRangerDirs} is specified.
#'   \item A list of vectors containing sample names to import. Each vector in
#'    the list corresponds to samples from one of \emph{cellRangerDirs}.
#'    These names are the same as the folder names under \emph{cellRangerDirs}.
#'    This is only for the case where \emph{cellRangerDirs} is specified.
#'   \item A vector of folder paths for the samples to import. This is only for
#'    the case where \emph{cellRangerDirs} is \code{NULL}.
#' }
#' The cells in the final SCE object will be ordered in the same order of
#' samples.
#' @param cellRangerOuts Character. It is the intermediate
#'  path to filtered or raw feature count file saved in sparse matrix format
#'  for each of \emph{samples}.
#'  Reference genome names might need to be
#'  appended for CellRanger version below 3.0.0 if reads were mapped to
#'  multiple genomes when running Cell Ranger pipeline.
#' @param matrixFileName Filename for the Market Exchange Format (MEX) sparse
#'  matrix file (.mtx file).
#' @param featuresFileName Filename for the feature annotation file. It can be
#'  \emph{features.tsv.gz} or \emph{genes.tsv}.
#' @param barcodesFileName Filename for the cell barcode list file.
#' @param gzipped Boolean. \code{TRUE} if the Cell Ranger output files
#'  (barcodes.tsv, features.tsv, and matrix.mtx) were
#'  gzip compressed. \code{FALSE} otherwise. This is true after Cell Ranger
#'  3.0.0 update. Default \code{TRUE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "DelayedArray" (as returned by
#'  \link[DelayedArray]{DelayedArray} function), "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "DelayedArray".
#' @details
#'  \code{importCellRangerV2} imports output from Cell Ranger V2.
#'  \code{importCellRangerV3} imports output from Cell Ranger V3. Some implicit
#'  assumptions are made in these two functions including \code{cellRangerOuts},
#'  \code{matrixFileName}, \code{featuresFileName}, and \code{barcodesFileName}.
#'  Alternatively, user can call \code{importCellRanger} to explicitly
#'  specify \code{cellRangerOuts}, \code{matrixFileName},
#'  \code{featuresFileName}, and \code{barcodesFileName}.
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
#'     cellRangerDirs = system.file("extdata", package = "singleCellTK"),
#'     samples = "hgmm_1k_v3_20x20")
#' @export
importCellRanger <- function(
    cellRangerDirs = NULL,
    samples = NULL,
    cellRangerOuts = "outs/filtered_feature_bc_matrix/",
    matrixFileName = "matrix.mtx.gz",
    featuresFileName = "features.tsv.gz",
    barcodesFileName = "barcodes.tsv.gz",
    gzipped = TRUE,
    class = "DelayedArray") {

    .importCellRanger(cellRangerDirs = cellRangerDirs,
        samples = samples,
        cellRangerOuts = cellRangerOuts,
        matrixFileName = matrixFileName,
        featuresFileName = featuresFileName,
        barcodesFileName = barcodesFileName,
        gzipped = gzipped,
        class = class)
}


#' @rdname importCellRanger
#' @export
importCellRangerV2 <- function(
    cellRangerDirs = NULL,
    samples = NULL,
    matrixFileName = "matrix.mtx",
    featuresFileName = "genes.tsv",
    barcodesFileName = "barcodes.tsv",
    gzipped = TRUE,
    class = "DelayedArray") {

    .importCellRanger(cellRangerDirs = cellRangerDirs,
        samples = samples,
        cellRangerOuts = "outs/filtered_gene_bc_matrices/",
        matrixFileName = matrixFileName,
        featuresFileName = featuresFileName,
        barcodesFileName = barcodesFileName,
        gzipped = gzipped,
        class = class)
}


#' @rdname importCellRanger
#' @examples
#' sce <- importCellRangerV3(
#'     cellRangerDirs = system.file("extdata", package = "singleCellTK"),
#'     samples = "hgmm_1k_v3_20x20")
#' @export
importCellRangerV3 <- function(
    cellRangerDirs = NULL,
    samples = NULL,
    matrixFileName = "matrix.mtx.gz",
    featuresFileName = "features.tsv.gz",
    barcodesFileName = "barcodes.tsv.gz",
    gzipped = TRUE,
    class = "DelayedArray") {

    .importCellRanger(cellRangerDirs = cellRangerDirs,
        samples = samples,
        cellRangerOuts = "outs/filtered_feature_bc_matrix/",
        matrixFileName = matrixFileName,
        featuresFileName = featuresFileName,
        barcodesFileName = barcodesFileName,
        gzipped = gzipped,
        class = class)
}
