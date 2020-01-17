.readMatrix <- function(path) {

  res <- readRDS(path)
  matrix <- t(as.matrix(res))

  if (class == "Matrix") {
    return(matrix)
  } else if (class == "DelayedArray") {
    res <- DelayedArray::DelayedArray(res)
    return(matrix)
  }
}


.readMetrics <- function(path) {

  res <- fread(path)
  res <- as.dataframe(res)
  return(res)

  if (ncol(res) == 1) {
    stop("There are no Cell or Gene Metrics!")
  }

}

.readEmptyDrops<-function(path){
  EmptyDrops<-read.csv(path)
}


.constructSCEFromOptimusOutputs <- function(matrixLocation,
  CellMetricsLocation,
  EmptyDropsLocation,
  GeneMetricsLocation) {

  c_me <- .readMetrics(file.path(CellMetricsLocation))
  cb <- c_me[,1]
  c_me <- c_me[,-1]
  emptydrops <- .readEmptyDrops(file.path(EmptyDropsLocation))
  c_me_full <- cbind(c_me, emptydrops)
  g_me <- .readMetrics(file.path(GeneMetricsLocation))
  g_me <- g_me[,g_me %in% rownames(ma)]
  gi <- g_me[,1]
  g_me <- g_me[,-1]
  ma <- .readMatrix(file.path(matrixLocation))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ma))
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(g_me,row.names = gi)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(c_me_full,row.names = cb)

  return(sce)
}


.checkArgsImportOptimus <- function(OptimusDirs) {

  if (!dir.exists(OptimusDirs)) {
      stop("OptimusDirs is NULL!")
  }
}


.importOptimus <- function(
  OptimusDirs,
  matrixLocation,
  CellMetricsLocation,
  EmptyDropsLocation,
  GeneMetricsLocation) {

  .checkArgsImportOptimus(OptimusDirs)
  sce <- .constructSCEFromOptimusOutputs(OptimusDirs,
    matrixLocation = matrixLocation,
    CellMetricsLocation = CellMetricsLocation,
    EmptyDropsLocation = EmptyDropsLocation,
    GeneMetricsLocation = GeneMetricsLocation)
   return(sce)
}


#' @name importOptimus
#' @rdname importOptimus
#' @title Construct SCE object from Optimus output
#' @description Read the barcodes, features (genes), and matrix from STARsolo
#'  output. Import them
#'  as one \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param STARsoloDirs A vector of root directories of STARsolo output files.
#'  The paths should be something like this:
#'  \bold{/PATH/TO/\emph{prefix}Solo.out}. For example: \code{./Solo.out}.
#'  Each sample should have its own path. Must have the same length as
#'  \code{samples}.
#' @param samples A vector of user-defined sample names for the sample to be
#'  imported. Must have the same length as \code{STARsoloDirs}.
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
importOptimus <- function() {

}
