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
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
        cell_barcode = as.character(barcodes[[1]]),
        column_name = coln,
        sample = sampleName,
        row.names = coln)

    return(sce)
}


.unionGeneMatrix <- function(geneUnion, matrix){
    missGene <- geneUnion[!geneUnion %in% rownames(matrix)]
    missMat <- Matrix::Matrix(0, nrow = length(missGene), ncol = ncol(matrix),
        dimnames = list(missGene, NULL))

    matb <- methods::as(matrix, "dgCMatrix")
    rownames(matb) <- rownames(matrix)

    mat <- rbind(matb, missMat)
    if (anyDuplicated(rownames(mat))) {
        mat <- mat[!duplicated(rownames(mat)), ]
        warning("Duplicated genes exist in count matrix. Filtered",
            " duplicated genes.")
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


.readBarcodesSEQC <- function(path) {
    res <- data.table::fread(path, header = FALSE, sep=",", colClasses = "character")
    res <- res[,-1,drop = FALSE]
    colnames(res) <- "cell_barcode"
    return(res)
}


.readFeaturesSEQC <- function(path) {
    res <- data.table::fread(path, header = FALSE, sep=",", colClasses = "character")
    res <- res[,-1,drop = FALSE]
    colnames(res) <- "feature_name"
    return(res)
}


.importSEQC <- function(
    seqcDirs,
    samples,
    prefix,
    gzipped,
    class,
    delayedArray,
    cbNotFirstCol,
    feNotFirstCol,
    combinedSample) {

    if (length(seqcDirs) != length(samples)) {
        stop("'seqcDirs' and 'samples' have unequal lengths!")
    }

    if (length(seqcDirs) != length(prefix)) {
        stop("'seqcDirs' and 'prefix' have unequal lengths!")
    }

    res <- vector("list", length = length(seqcDirs))
    cb <- vector("list", length = length(seqcDirs))
    fe <- vector("list", length = length(seqcDirs))
    mat <- vector("list", length = length(seqcDirs))

    for (i in seq_along(seqcDirs)) {
        dir <- seqcDirs[i]
        matrixFile <- paste(prefix[i], 'sparse_molecule_counts.mtx', sep = "_")
        featuresFile <- paste(prefix[i], 'sparse_counts_genes.csv', sep = "_")
        barcodesFile <- paste(prefix[i], 'sparse_counts_barcodes.csv',
            sep = "_")

        cb[[i]] <- .readBarcodesSEQC(file.path(dir, barcodesFile))
        fe[[i]] <- .readFeaturesSEQC(file.path(dir, featuresFile))

        mat[[i]] <- .readMatrixMM(file.path(dir, matrixFile),
            gzipped = gzipped, class = class, delayedArray = delayedArray)
        mat[[i]] <- t(mat[[i]])
        rownames(mat[[i]]) <- fe[[i]][[1]]
    }

    if (isTRUE(combinedSample) & length(seqcDirs) > 1) {
        geneUnion <- .getGeneUnion(fe)
        for (i in seq_along(seqcDirs)) {
            matrix <- .unionGeneMatrix(geneUnion = geneUnion, matrix = mat[[i]])
            matrix <- matrix[geneUnion, ]
            feature <- S4Vectors::DataFrame('feature_name' = rownames(matrix))

            scei <- .constructSCEFromSeqcOutputs(
                sampleName = samples[i],
                matrix = matrix,
                features = feature,
                barcodes = cb[[i]])
            res[[i]] <- scei

        }
        sce <- do.call(SingleCellExperiment::cbind, res)
        return(sce)

    } else {
        for (i in seq_along(seqcDirs)) {
            scei <- .constructSCEFromSeqcOutputs(
                sampleName = samples[i],
                matrix = mat[[i]],
                features = fe[[i]],
                barcodes = cb[[i]])
            res[[i]] <- scei
        }
        if (length(seqcDirs) == 1) {
            return(res[[1]])
        } else {
            return(res)
        }
    }
}


#' @name importSEQC
#' @rdname importSEQC
#' @title Construct SCE object from seqc output
#' @description Read the filtered barcodes, features, and matrices for all
#'  samples from (preferably a single run of) seqc output. Import and
#'  combine them as one big \link[SingleCellExperiment]{SingleCellExperiment}
#'  object.
#' @param seqcDirs A vector of paths to seqc output files. Each sample
#'  should have its own path. For example: \code{./pbmc_1k_50x50}.
#'  Must have the same length as \code{samples}.
#' @param samples A vector of user-defined sample names for the samples to be
#'  imported. Must have the same length as \code{seqcDirs}.
#' @param prefix A vector containing the prefix of file names within each
#'  sample directory. It cannot be null and the vector should have the same
#'  length as \emph{samples}.
#' @param gzipped Boolean. \code{TRUE} if the seqc output files
#'  (sparse_counts_barcode.csv, sparse_counts_genes.csv, and
#'  sparse_molecule_counts.mtx)
#'  were gzip compressed. \code{FALSE} otherwise. Default seqc outputs are
#'  not gzipped.
#' Default \code{FALSE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @param feNotFirstCol Boolean. \code{TRUE} if first column of
#'  sparse_counts_genes.csv
#' is row index and it will be removed. \code{FALSE} the first column will
#'  be kept.
#' @param cbNotFirstCol Boolean. \code{TRUE} if first column of
#'  sparse_counts_barcode.csv
#' is row index and it will be removed. \code{FALSE} the first column will
#'  be kept.
#' @param combinedSample Boolean. If \code{TRUE}, \code{importSEQC} returns a
#' \code{SingleCellExperiment} object containing the combined count matrix,
#'  feature annotations
#'  and the cell annotations. If \code{FALSE}, \code{importSEQC} returns a
#'  list containing multiple
#'  \code{SingleCellExperiment} objects. Each \code{SingleCellExperiment}
#'  contains count matrix, feature annotations and cell annotations for
#'  each sample.
#' @details
#' \code{importSEQC} imports output from seqc.
#'  The default sparse_counts_barcode.csv or sparse_counts_genes.csv from
#'  seqc output
#'  contains two columns. The first column is row index and the second column
#'  is cell-barcode
#'  or gene symbol. \code{importSEQC} will remove first column. Alternatively,
#'  user can call
#'  \code{cbNotFirstCol} or \code{feNotFirstCol} as FALSE to keep the first
#'  column of these files.
#'  When \code{combinedSample} is TRUE, \code{importSEQC} will combined count
#'  matrix with genes detected in at least one sample.
#' @return A \code{SingleCellExperiment} object containing the combined count
#'  matrix, the feature annotations, and the cell annotation.
#' @examples
#' # Example #1
#' # The following filtered feature, cell, and matrix files were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/
#' # 3.0.0/pbmc_1k_v3
#' # The top 50 hg38 genes are included in this example.
#' # Only the top 50 cells are included.
#' sce <- importSEQC(
#'     seqcDirs = system.file("extdata/pbmc_1k_50x50", package = "singleCellTK"),
#'     samples = "pbmc_1k_50x50",
#'     prefix = "pbmc_1k",
#'     combinedSample = FALSE)
#' @export
importSEQC <- function(
    seqcDirs = NULL,
    samples = NULL,
    prefix = NULL,
    gzipped = FALSE,
    class = c("Matrix", "matrix"),
    delayedArray = TRUE,
    cbNotFirstCol = TRUE,
    feNotFirstCol = TRUE,
    combinedSample = TRUE) {

    class <- match.arg(class)

    .importSEQC(seqcDirs = seqcDirs,
        samples = samples,
        prefix = prefix,
        gzipped = gzipped,
        class = class,
        delayedArray = delayedArray,
        cbNotFirstCol = cbNotFirstCol,
        feNotFirstCol = feNotFirstCol,
        combinedSample = combinedSample)
}
