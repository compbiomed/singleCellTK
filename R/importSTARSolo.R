
# dir <- "outs/filtered_feature_bc_matrix/"
.constructSCEFromSTARsoloOutputs <- function(dir,
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


# main function
.importSTARsolo <- function(STARsoloDirs,
    samples,
    STARsoloOuts,
    matrixFileName,
    featuresFileName,
    barcodesFileName,
    gzipped,
    class) {


    if (length(STARsoloDirs) != length(samples)) {
        stop("'STARsoloDirs' and 'samples' have unequal lengths!")
    }

    res <- vector("list", length = length(samples))

    for (i in seq_along(samples)) {
        dir <- file.path(BUStoolsDirs[i], STARsoloOuts)
        scei <- .constructSCEFromSTARsoloOutputs(dir,
            sample = samples[i],
            matrixFileName = matrixFileName,
            featuresFileName = featuresFileName,
            barcodesFileName = barcodesFileName,
            gzipped = gzipped,
            class = class)
        res[[i]] <- scei
    }

    sce <- do.call(BiocGenerics::cbind, res)
    return(sce)
}


#' @name importSTARsolo
#' @rdname importSTARsolo
#' @title Construct SCE object from STARsolo output
#' @description Read the barcodes, features (genes), and matrix from STARsolo
#'  output. Import them
#'  as one \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param STARsoloDirs A vector of root directories of STARsolo output files.
#'  The paths should be something like this:
#'  \bold{/PATH/TO/\emph{prefix}Solo.out}. For example: \code{./Solo.out}.
#'  Each sample should have its own path. Must have the same length as
#'  \code{samples}.
#' @param samples A vector of user-defined sample names for the sample to be
#'  imported. Must have the same length as \code{BUStoolsDirs}.
#' @param STARsoloOuts Character. It is the intermediate
#'  path to filtered or raw feature count file saved in sparse matrix format
#'  for each of \emph{samples}. Default "Gene/filtered"  which works for STAR
#'  2.7.3a.
#' @param matrixFileName Filename for the Market Exchange Format (MEX) sparse
#'  matrix file (.mtx file).
#' @param featuresFileName Filename for the feature annotation file.
#' @param barcodesFileName Filename for the cell barcode list file.
#' @param gzipped Boolean. \code{TRUE} if the STARsolo output files
#'  (barcodes.tsv, features.tsv, and matrix.mtx) were
#'  gzip compressed. \code{FALSE} otherwise. This is \code{FALSE} in STAR
#'  2.7.3a. Default \code{FALSE}.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "DelayedArray" (as returned by
#'  \link[DelayedArray]{DelayedArray} function), "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @return A \code{SingleCellExperiment} object containing the count
#'  matrix, the gene annotation, and the cell annotation.
#' @examples
#' # Example #1
#' # FASTQ files were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' # /pbmc_1k_v3
#' # They were concatenated as follows:
#' # cat pbmc_1k_v3_S1_L001_R1_001.fastq.gz pbmc_1k_v3_S1_L002_R1_001.fastq.gz >
#' # pbmc_1k_v3_R1.fastq.gz
#' # cat pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_S1_L002_R2_001.fastq.gz >
#' # pbmc_1k_v3_R2.fastq.gz
#' # The following STARsolo command generates the filtered feature, cell, and
#' # matrix files
#' # STAR \
#' #   --genomeDir ./index \
#' #   --readFilesIn ./pbmc_1k_v3_R2.fastq.gz \
#' #                 ./pbmc_1k_v3_R1.fastq.gz \
#' #   --readFilesCommand zcat \
#' #   --outSAMtype BAM Unsorted \
#' #   --outBAMcompression -1 \
#' #   --soloType CB_UMI_Simple \
#' #   --soloCBwhitelist ./737K-august-2016.txt \
#' #   --soloUMIlen 12
#'
#' # The top 20 genes and the first 20 cells are included in this example.
#' sce <- importSTARsolo(
#'   STARsoloDirs = system.file("extdata/STARsolo_PBMC_1k_v3_20x20",
#'     package = "singleCellTK"),
#'   samples = "PBMC_1k_v3_20x20")
#' @export
importSTARsolo <- function(
    STARsoloDirs,
    samples,
    STARsoloOuts = "Gene/filtered",
    matrixFileName = "matrix.mtx",
    featuresFileName = "features.tsv",
    barcodesFileName = "barcodes.tsv",
    gzipped = FALSE,
    class = "Matrix") {

    .importSTARsolo(
        STARsoloDirs = STARsoloDirs,
        samples = samples,
        STARsoloOuts = STARsoloOuts,
        matrixFileName = matrixFileName,
        featuresFileName = featuresFileName,
        barcodesFileName = barcodesFileName,
        gzipped = gzipped,
        class = class)
}
