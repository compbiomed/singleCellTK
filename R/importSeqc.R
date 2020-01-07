# .readBarcodes <- function(path,
#     header = FALSE,
#     colname = "cell_barcode", 
#     removeFirstCol = TRUE) {

#     res <- data.table::fread(path, header = header)

#     if (ncol(res) == 1) {
#         colnames(res) <- colname
#     } else {
#         if (ncol(res) == 2) {
#             if (removeFirstCol) {
#                 message("First column of barcode file was row index and it was removed.")
#                 res <- res[, -1]
#                 colnames(res) <- colname
#             }
#         } else if (ncol(res) > 2) {
#             warning("'barcodes' file contains >2 columns!",
#             " The column names are kept as is. ")
#         }
#     }
#     return(res)
# }

# .readFeatures <- function(path,
#     header = FALSE,
#     colname = "feature_name",
#     removeFirstCol = TRUE) {

#     res <- data.table::fread(path, header = header)
#     if (ncol(res) == 1) {
#         colnames(res) <- colname
#     } else {
#         if (ncol(res) == 2) {
#             if (removeFirstCol) {
#                 message("First column of gene file was row index and it was removed.")
#                 res <- res[, -1]
#                 colnames(res) <- colname
#             }
#         } else if (ncol(res) > 2) {
#             warning("'barcodes' file contains >2 columns!",
#             " The column names are kept as is. ")
#         }
#     }
#     return(res)
# }


# .readMatrixMM <- function(path, gzipped = FALSE, class = "DelayedArray") {
#     if (isTRUE(gzipped)) {
#         path <- gzfile(path)
#     }

#     res <- Matrix::readMM(path)
#     res <- t(res)
#     if (class == "Matrix") {
#         return(res)
#     } else if (class == "DelayedArray") {
#         res <- DelayedArray::DelayedArray(res)
#         return(res)
#     } else if (class == "matrix") {
#         res <- as.matrix(res)
#         return(res)
#     }
# }

.constructSCEFromSeqcOutputs <- function(
    sampleName,
    matrix,
    features,
    barcodes) {

    coln <- paste(sampleName, barcodes[[1]], sep = "_")
    rownames(matrix) <- features[[1]]

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix))
    SummarizedExperiment::rowData(sce) <- features
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(barcodes,
        column_name = coln,
        sample = sampleName,
        row.names = coln)

    return(sce)
}

## We don't need .getOutputFolderPath for seqc output. 

.checkArgsImportSeqc <- function(SeqcDirs, samples, class, prefix) {
    if (is.null(SeqcDirs)) {
        if (is.null(samples)) {
            stop("samples can not be NULL if SeqcDirs is NULL!")
        }
        for (i in seq_along(samples)) {
            if (!dir.exists(samples[i])) {
                stop("Sample folder does not exist!\n", samples[i])
            }
        }
    } else {
        if (is.null(samples)) {
            for (i in seq_along(SeqcDirs)) {
                if (length(list.dirs(SeqcDirs[i],
                    recursive = FALSE)) == 0) {
                    warning("Empty folder. Skipping SeqcDirs ",
                        SeqcDirs[i])
                }
            }
        } else {
            if (!(length(samples) == length(SeqcDirs))) {
                stop("Length of samples is not equal to length of ",
                    "SeqcDirs!")
            } else {
                for (i in seq_along(SeqcDirs)) {
                    paths <- file.path(SeqcDirs[i], samples[i])
                    ## why need a for loop below
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
        stop("Invalid 'class' argument! ", "Only accept 'DelayedArray', 'Matric' or 'matrix'")
    }

    if (is.null(prefix)) {
        stop("prefix of output files could not be null ")
    }

}

# .getSamplesPaths <- function(SeqcDirs, samples){
#     if (is.null(SeqcDirs)){
#         res <- samples
#     } else {
#         if (is.null(samples)){
#             ## We assume there are only sample directories udner SeqcDirs
#             res <- list.dirs(SeqcDirs, recursive = FALSE)
#         } else {
#             res <- vector("list", length = length(SeqcDirs))
#             for (i in seq_along(SeqcDirs)) {
#                 res[[i]] <- file.path(SeqcDirs[i], samples[i])
#             }
#             res <- unlist(res)
#         }
#     }
#     return(res)
# }

# .getSampleNames <- function(samplesDir) {
#     res <- basename(samplesDir)
#     return(res)
# }

.unionGeneMatrix <- function(geneUnion, matrix){
    missGene <- geneUnion[!geneUnion %in% rownames(matrix)]
    missMat <- Matrix::Matrix(0, nrow = length(missGene), ncol = ncol(matrix),
        dimnames = list(missGene, NULL))

    mat <- rbind(matrix,missMat)
    if (anyDuplicated(rownames(mat))) {
        mat <- mat[!duplicated(rownames(mat)), ]
        warning('Duplicated genes exist in count matrix. Filtered duplicated genes.')
    }
    return(mat)
}

.getGeneUnion <- function(geneList){
    gene <- geneList
    for (i in seq_along(geneList)){
        gene[[i]] <- geneList[[i]][[1]]
    }
    
    geneUnion <- base::Reduce(union, gene)
    return(geneUnion)
}

.importSeqc <- function(
    SeqcDirs,
    samples,
    prefix, 
    gzipped,
    class,
    cbNotFirstCol,
    feNotFirstCol,
    combinedSample) {

    .checkArgsImportSeqc(SeqcDirs, samples, class, prefix)
    sampleDirs <- .getSamplesPaths(SeqcDirs, samples)

    res <- vector("list", length = length(sampleDirs))
    cb <- vector("list", length = length(sampleDirs))
    fe <- vector("list", length = length(sampleDirs))
    mat <- vector("list", length = length(sampleDirs))

    for (i in seq_along(sampleDirs)) {
        matrixFile <- paste(prefix[i], 'sparse_molecule_counts.mtx', sep = "_")
        featuresFile <- paste(prefix[i], 'sparse_counts_genes.csv', sep = "_")
        barcodesFile <- paste(prefix[i], 'sparse_counts_barcodes.csv', sep = "_")

        cb[[i]] <- .readBarcodes(file.path(sampleDirs[i], barcodesFile))
        if (isTRUE(cbNotFirstCol)) {
            message("First column of barcode file was row index and it was removed.")
            cb[[i]] <- cb[[i]][, -1]            
        }

        fe[[i]] <- .readFeatures(file.path(sampleDirs[i], featuresFile))
        if (isTRUE(feNotFirstCol)) {
            message("First column of gene file was row index and it was removed.")
            fe[[i]] <- fe[[i]][, -1]            
        }

        mat[[i]] <- .readMatrixMM(file.path(sampleDirs[i], matrixFile), 
            gzipped = gzipped, class = 'Matrix')
        mat[[i]] <- t(mat[[i]])
        rownames(mat[[i]]) <- fe[[i]][[1]]
    }

    if (isTRUE(combinedSample) & length(sampleDirs) > 1) {
        geneUnion <- .getGeneUnion(fe)
        for (i in seq_along(sampleDirs)) {
            matrix <- .unionGeneMatrix(geneUnion = geneUnion, matrix = mat[[i]])
            matrix <- matrix[geneUnion, ]
            feature <- S4Vectors::DataFrame('feature_name' = rownames(matrix))

            scei <- .constructSCEFromSeqcOutputs(
                sampleName = .getSampleNames(sampleDirs[i]),
                matrix = matrix,
                features = feature,
                barcodes = cb[[i]])
            res[[i]] <- scei

        }
        sce <- do.call(BiocGenerics::cbind, res)
        return(sce)

    } else {
        for (i in seq_along(sampleDirs)) {
            if (class == 'DelayedArray') {
                mat[[i]] <- DelayedArray::DelayedArray(mat[[i]])
            } else if (class == 'matrix') {
                mat[[i]] <- as.matrix(mat[[i]])
            }

            scei <- .constructSCEFromSeqcOutputs(
                sampleName = .getSampleNames(sampleDirs[i]),
                matrix = mat[[i]],
                features = fe[[i]],
                barcodes = cb[[i]])
            res[[i]] <- scei
        }
        if (length(sampleDirs) == 1){
            return(res[[1]])
        } else {
            return(res) 
        }
    }
}

#' @name importSeqc
#' @rdname importSeqc
#' @title Construct SCE object from seqc output
#' @description Read the filtered barcodes, features, and matrices for all
#'  samples from (preferably a single run of) seqc output. Import and
#'  combine them as one big \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param SeqcDirs The root directories where seqc was run. These
#'  folders should contain sample specific folders. Default \code{NULL},
#'  meaning the paths for each sample will be specified in \emph{samples}
#'  argument.
#' @param samples Default \code{NULL}. Can be one of
#' \itemize{
#'   \item \code{NULL}. All samples within \emph{SeqcDirs} will be
#'    imported. The order of samples will be first determined by the order of
#'    \emph{SeqcDirs} and then by \link[base]{list.files}. This is only
#'    for the case where \emph{SeqcDirs} is specified.
#'   \item A vector containing sample names to import.
#'    These names are the same as the folder names under \emph{SeqcDirs}.
#'    This is only for the case where \emph{SeqcDirs} is specified.
#'   \item A vector of folder paths for the samples to import. This is only for
#'    the case where \emph{SeqcDirs} is \code{NULL}.
#' }
#' The cells in the final SCE object will be ordered in the same order of
#' samples.
#' @param prefix A vector containing the prefix of file names within each sample directory. 
#' It cannot be null and the vector should have the same length as \emph{samples}.
#' @param gzipped Boolean. \code{TRUE} if the seqc output files
#' (sparse_counts_barcode.csv, sparse_counts_genes.csv, and sparse_molecule_counts.mtx) 
#' were gzip compressed. \code{FALSE} otherwise. Default seqc outputs are not gzipped.
#' Default \code{FALSE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#' object. Can be one of "DelayedArray" (as returned by
#' \link[DelayedArray]{DelayedArray} function), "Matrix" (as returned by
#' \link[Matrix]{readMM} function), or "matrix" (as returned by
#' \link[base]{matrix} function). Default "DelayedArray".
#' @param feNotFirstCol Boolean. \code{TRUE} if first column of sparse_counts_genes.csv
#' is row index and it will be removed. \code{FALSE} the first column will be kept. 
#' @param cbNotFirstCol Boolean. \code{TRUE} if first column of sparse_counts_barcode.csv
#' is row index and it will be removed. \code{FALSE} the first column will be kept. 
#' @param combinedSample Boolean. If \code{TRUE}, \code{importSeqc} returns a 
#' \code{SingleCellExperiment} object containing the combined count matrix, feature annotations 
#' and the cell annotations. If \code{FALSE}, \code{importSeqc} returns a list containing multiple 
#' \code{SingleCellExperiment} objects. Each \code{SingleCellExperiment} contains count matrix
#' , feature anotations and cell annotations for each sample. 
#' @details
#' \code{importSeqc} imports output from seqc.
#' The default sparse_counts_barcode.csv or sparse_counts_genes.csv from seqc output
#' contains two columns. The first column is row index and the second column is cell-barcode 
#' or gene symbol. \code{importSeqc} will remove first column. Alternatively, user can call 
#' \code{cbNotFirstCol} or \code{feNotFirstCol} as FALSE to keep the first column
#' of these files.  
#' When \code{combinedSample} is TRUE, \code{importSeqc} will combined count matrix
#' with genes detected in at least one sample. 
#' @return A \code{SingleCellExperiment} object containing the combined count
#'  matrix, the feature annotations, and the cell annotation.
#' @examples
#' # Example #1
#' # The following filtered feature, cell, and matrix files were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/
#' # 3.0.0/pbmc_1k_v3
#' # The top 50 hg38 genes are included in this example.
#' # Only the top 50 cells are included.
#' sce <- importSeqc(
#'     SeqcDirs = system.file("extdata", package = "singleCellTK"),
#'     samples = "pbmc_1k_50x50", 
#'     prefix = 'pbmc_1k', 
#'     combinedSample = FALSE)
#' @export

importSeqc <- function(
    SeqcDirs = NULL,
    samples = NULL,
    prefix = NULL, 
    gzipped = FALSE,
    class = "DelayedArray", 
    cbNotFirstCol = TRUE,
    feNotFirstCol = TRUE,
    combinedSample = TRUE) {

    .importSeqc(SeqcDirs = SeqcDirs,
        samples = samples,
        prefix = prefix, 
        gzipped = gzipped,
        class = class, 
        cbNotFirstCol = cbNotFirstCol,
        feNotFirstCol = feNotFirstCol, 
        combinedSample = combinedSample)
}
