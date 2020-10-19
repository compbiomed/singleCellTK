
#' @importFrom reticulate import
.readMatrixNpz <- function(matrixLocation,
  colIndexLocation,
  rowIndexLocation,
  class,
  delayedArray) {

  ## Now importing these functions in 'reticulate_setup.R' file
  #  sparse <- reticulate::import("scipy.sparse")
  #  numpy <- reticulate::import("numpy")
  if (!reticulate::py_module_available(module = "scipy.sparse")) {
    stop("Error!", "Cannot find python module 'scipy.sparse', please install Conda and run sctkPythonInstallConda()
         or run sctkPythonInstallVirtualEnv(). If one of these have been previously run to install the modules,
         make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(), respectively, if R has been
         restarted since the module installation. Alternatively, scipy can be installed on the local machine
         with pip (e.g. pip install scipy) and then the 'use_python()' function from the 'reticulate' package
         can be used to select the correct Python environment.")
  }
  if (!reticulate::py_module_available(module = "numpy")) {
    stop("Error!", "Cannot find python module 'numpy', please install Conda and run sctkPythonInstallConda()
         or run sctkPythonInstallVirtualEnv(). If one of these have been previously run to install the modules,
         make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(), respectively, if R has been
         restarted since the module installation. Alternatively, numpy can be installed on the local machine
         with pip (e.g. pip install numpy) and then the 'use_python()' function from the 'reticulate' package
         can be used to select the correct Python environment.")
  }

  error <- try({
    mat <- sparse$load_npz(matrixLocation)
    colIndex <- as.vector(numpy$load(colIndexLocation, allow_pickle = TRUE))
    rowIndex <- as.vector(numpy$load(rowIndexLocation, allow_pickle = TRUE))
    colnames(mat) <- colIndex
    rownames(mat) <- rowIndex
    mat <- t(mat)

    ## Convert to "dgCMatrix"
    newM <- Matrix::Matrix(mat[,1], nrow=nrow(mat))
    newM <- methods::as(newM, "dgCMatrix")
    breaks <- seq(2, ncol(mat), by=1000)
    if(length(breaks) > 2) {
      for(i in seq(2, length(breaks))) {
        ix <- seq(breaks[i-1], (breaks[i]-1))
        newM <- cbind(newM, mat[,ix])
      }
      ix <- seq(utils::tail(breaks, n = 1), ncol(mat))
      newM <- cbind(newM, mat[,ix])
    } else {
      ix <- seq(2, ncol(mat))
      newM <- cbind(newM, mat[,ix])
    }

    colnames(newM) <- colnames(mat)
    rownames(newM) <- rownames(mat)
    mat <- newM
  }, silent = TRUE)

  if(inherits(error, "try-error")) {
    stop(paste0("importOptimus did not complete successfully. SCE could not be generated. Error given during the import process: \n\n", error))
  }

  if (class == "matrix") {
    mat <- as.matrix(mat)
  }

  if (isTRUE(delayedArray)) {
    mat <- DelayedArray::DelayedArray(mat)
  }
  return(mat)
}


.readMetricsOptimus <- function(path) {
  metrics <- data.table::fread(path)
  return(metrics)
}


.readEmptyDrops <- function(path) {
  emptyDrops <- data.table::fread(path)
  colnames(emptyDrops) <- paste0("dropletUtils_emptyDrops_",
    colnames(emptyDrops))
  return(emptyDrops)
}


.combineColData <- function(colnames, cellMetrics, emptyDrops) {
  cd <- data.table::data.table(CellId = colnames)
  cd <- merge(cd,
    cellMetrics,
    by.x = "CellId",
    by.y = "V1",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE)

  if (!is.null(emptyDrops)) {
    cd <- merge(cd,
      emptyDrops,
      by.x = "CellId",
      by.y = "dropletUtils_emptyDrops_CellId",
      all.x = TRUE,
      all.y = FALSE,
      sort = FALSE)
  }

  return(cd)
}


.combineRowData <- function(rownames, geneMetrics) {
  rd <- data.table::data.table(feature_ID = rownames)
  rd <- merge(rd,
    geneMetrics,
    by.x = "feature_ID",
    by.y = "V1",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE)
  return(rd)
}


.constructSCEFromOptimusOutputs <- function(dir,
  sample,
  matrixLocation,
  colIndexLocation,
  rowIndexLocation,
  cellMetricsLocation,
  geneMetricsLocation,
  emptyDropsLocation,
  class,
  delayedArray) {

  mat <- .readMatrixNpz(file.path(dir, matrixLocation),
    file.path(dir, colIndexLocation),
    file.path(dir, rowIndexLocation),
    class,
    delayedArray)

  cellMetrics <- .readMetricsOptimus(file.path(dir, cellMetricsLocation))
  geneMetrics <- .readMetricsOptimus(file.path(dir, geneMetricsLocation))

  if (!is.null(geneMetricsLocation)) {
    emptyDrops <- .readEmptyDrops(file.path(dir, emptyDropsLocation))
  }

  cd <- .combineColData(colnames(mat), cellMetrics, emptyDrops)
  rd <- .combineRowData(rownames(mat), geneMetrics)

  coln <- paste(sample, colnames(mat), sep = "_")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat))
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd,
    column_name = coln,
    sample = sample,
    row.names = coln)
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(rd)

  return(sce)
}


.checkArgsImportOptimus <- function(OptimusDirs,
  samples) {

  if (length(OptimusDirs) != length(samples)) {
    stop("'OptimusDirs' and 'samples' have unequal lengths!")
  }
}


.importOptimus <- function(OptimusDirs,
  samples,
  matrixLocation,
  colIndexLocation,
  rowIndexLocation,
  cellMetricsLocation,
  geneMetricsLocation,
  emptyDropsLocation,
  class,
  delayedArray) {

  .checkArgsImportOptimus(OptimusDirs, samples)

  res <- vector("list", length = length(samples))

  for (i in seq_along(samples)) {
    scei <- .constructSCEFromOptimusOutputs(OptimusDirs[i],
      sample = samples[i],
      matrixLocation = matrixLocation,
      colIndexLocation = colIndexLocation,
      rowIndexLocation = rowIndexLocation,
      cellMetricsLocation = cellMetricsLocation,
      geneMetricsLocation = geneMetricsLocation,
      emptyDropsLocation = emptyDropsLocation,
      class = class,
      delayedArray = delayedArray)
    res[[i]] <- scei
  }

  sce <- do.call(SingleCellExperiment::cbind, res)
  return(sce)
}


#' @name importOptimus
#' @rdname importOptimus
#' @title Construct SCE object from Optimus output
#' @description Read the barcodes, features (genes), and matrices from Optimus
#'  outputs. Import them
#'  as one \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param OptimusDirs A vector of root directories of Optimus output files.
#'  The paths should be something like this:
#'  \code{/PATH/TO/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530}.
#'  Each entry in \code{OptimusDirs} is considered a sample and should have
#'  its own path. Must have the same length as \code{samples}.
#' @param samples A vector of user-defined sample names for the sample to be
#'  imported. Must have the same length as \code{OptimusDirs}.
#' @param matrixLocation Character. It is the intermediate
#'  path to the filtered count maxtrix file saved in sparse matrix format
#'  (\code{.npz}). Default
#'  \code{call-MergeCountFiles/sparse_counts.npz} which works for
#'  optimus_v1.4.0.
#' @param colIndexLocation Character. The intermediate path to the barcode
#'  index file. Default \code{call-MergeCountFiles/sparse_counts_col_index.npy}.
#' @param rowIndexLocation Character. The intermediate path to the feature
#'  (gene) index file. Default
#'  \code{call-MergeCountFiles/sparse_counts_row_index.npy}.
#' @param cellMetricsLocation Character. It is the intermediate
#'  path to the cell metrics file (\code{merged-cell-metrics.csv.gz}). Default
#'  \code{call-MergeCellMetrics/merged-cell-metrics.csv.gz} which works for
#'  optimus_v1.4.0.
#' @param geneMetricsLocation Character. It is the intermediate
#'  path to the feature (gene) metrics file (\code{merged-gene-metrics.csv.gz}).
#'  Default \code{call-MergeGeneMetrics/merged-gene-metrics.csv.gz} which works
#'  for optimus_v1.4.0.
#' @param emptyDropsLocation Character. It is the intermediate
#'  path to \link[DropletUtils]{emptyDrops} metrics file
#'  (\code{empty_drops_result.csv}).
#'  Default \code{call-RunEmptyDrops/empty_drops_result.csv} which works for
#'  optimus_v1.4.0.
#' @param class Character. The class of the expression matrix stored in the SCE
#'  object. Can be one of "Matrix" (as returned by
#'  \link[Matrix]{readMM} function), or "matrix" (as returned by
#'  \link[base]{matrix} function). Default "Matrix".
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  containing the count
#'  matrix, the gene annotation, and the cell annotation.
#' @examples
#' file.path <- system.file("extdata/Optimus_20x1000",
#'   package = "singleCellTK")
#' \dontrun{
#' sce <- importOptimus(OptimusDirs = file.path,
#'   samples = "Optimus_20x1000")
#' }
#' @export
importOptimus <- function(OptimusDirs,
  samples,
  matrixLocation = "call-MergeCountFiles/sparse_counts.npz",
  colIndexLocation = "call-MergeCountFiles/sparse_counts_col_index.npy",
  rowIndexLocation = "call-MergeCountFiles/sparse_counts_row_index.npy",
  cellMetricsLocation = "call-MergeCellMetrics/merged-cell-metrics.csv.gz",
  geneMetricsLocation = "call-MergeGeneMetrics/merged-gene-metrics.csv.gz",
  emptyDropsLocation = "call-RunEmptyDrops/empty_drops_result.csv",
  class = c("Matrix", "matrix"),
  delayedArray = TRUE) {

  class <- match.arg(class)

  .importOptimus(OptimusDirs = OptimusDirs,
    samples = samples,
    matrixLocation = matrixLocation,
    colIndexLocation = colIndexLocation,
    rowIndexLocation = rowIndexLocation,
    cellMetricsLocation = cellMetricsLocation,
    geneMetricsLocation = geneMetricsLocation,
    emptyDropsLocation = emptyDropsLocation,
    class = class,
    delayedArray = delayedArray)

}
