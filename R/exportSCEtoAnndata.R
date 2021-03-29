#' @title Export a \link[SingleCellExperiment]{SingleCellExperiment} R object as
#' Python annData object
#' @description Writes all assays, colData, rowData, reducedDims, and altExps objects in a
#' \link[SingleCellExperiment]{SingleCellExperiment} to a Python annData object in the .h5ad format
#' All parameters of Anndata.write_h5ad function (https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.write_h5ad.html)
#' are available as parameters to this export function and set to defaults. Defaults can be
#' overridden at function call.
#' @param sce \link[SingleCellExperiment]{SingleCellExperiment} R object to be
#'  exported.
#' @param useAssay Character. The name of assay of
#' interests that will be set as the primary matrix of the output AnnData.
#' Default \code{"counts"}.
#' @param outputDir Path to the directory where .h5ad outputs will be written. Default is the current working directory.
#' @param prefix Prefix to use for the name of the output file. Default \code{"sample"}.
#' @param overwrite Boolean. Default \code{TRUE}.
#' @param compression If output file compression is required, this variable accepts
#' 'gzip' or 'lzf' as inputs. Default \code{None}.
#' @param compressionOpts Integer. Sets the compression level
#' @param forceDense Default \code{False} Write sparse data as a dense matrix.
#' Refer \code{anndata.write_h5ad} documentation for details. Default \code{NULL}.
#' @return Generates a Python anndata object containing data from \code{inSCE}.
#' @examples
#' data(sce_chcl, package = "scds")
#' \dontrun{
#' exportSCEtoAnnData(sce=sce_chcl, compression="gzip")
#' }
#' @export
exportSCEtoAnnData <- function(sce,
                                useAssay = 'counts',
                                outputDir = "./",
                                prefix = "sample",
                                overwrite = TRUE,
                                compression = c('None','lzf','gzip'),
                                compressionOpts = NULL,
                                forceDense = c('False','True')){
  compression <- match.arg(compression)
  forceDense <- match.arg(forceDense)
  if (compression == 'None'){
    compression <- NULL
  }

  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scrublet can be installed on the local machine",
            "with pip (e.g. pip install --user scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(sce)}

  AssayName <- SummarizedExperiment::assayNames(sce)
  for (assay in AssayName){
    if (!methods::is(SummarizedExperiment::assay(sce, assay), 'dgCMatrix')) {
      SummarizedExperiment::assay(sce, assay) <- .convertToMatrix(SummarizedExperiment::assay(sce, assay))
    }
  }


  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  annData <- .sce2adata(sce,useAssay)
  fileName <- paste0(prefix,".h5ad")
  filePath <- file.path(outputDir,fileName)

  if (file.exists(filePath) && !isTRUE(overwrite)) {
    stop(paste0(path, " already exists. Change 'outputDir' or set 'overwrite' to TRUE."))
    }

  annData$write_h5ad(filePath,
                     compression = compression,
                     compression_opts = compressionOpts,
                     force_dense = forceDense)
}
