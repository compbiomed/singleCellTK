#' @title Export a \link[SingleCellExperiment]{SingleCellExperiment} R object as 
#' Python Anndata object 
#' @description Writes all assays, colData, rowData, reducedDims, and altExps objects in a
#' \link[SingleCellExperiment]{SingleCellExperiment} to a Python Anndata object in the .h5ad format
#' All parameters of Anndata.write_h5ad function (https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.write_h5ad.html)
#' are available as parameters to this export function and set to defaults. Defaults can be
#' overridden at function call.
#' @param sce \link[SingleCellExperiment]{SingleCellExperiment} R object to be
#'  exported.
#' @param outputDir Path to the directory where .h5ad outputs will be written
#' @param overwrite Boolean. Default \code{TRUE}.
#' @param compression Default \code{NULL}.If output file compression is required, this variable accepts
#' 'gzip' or 'lzf' as inputs.
#' #' @param compression_opts Integer. Default \code{NULL} Sets the compression level
#' @param as_dense Default \code{NULL} Sparse arrays in AnnData object to write as dense. 
#' Refer anndata.write_h5ad documentation for details
#' @param force_dense Default \code{NULL} Write sparse data as a dense matrix.
#' Refer anndata.write_h5ad documentation for details
#' @examples
#' data(sce_chcl, package = "scds")
#' exportSCEtoAnnData(sce=sce_chcl, compression="gzip")
#'
#' @export

exportSCEtoAnnData <- function(sce, 
                                mainAssay='counts',
                                outputDir="./",
                                sample = "sample",
                                overwrite=TRUE,
                                compression= NULL,
                                compression_opts = NULL,
                                as_dense="()",
                                force_dense= NULL){
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
  #sc <- reticulate::import("scanpy")
  anndata <- .sce2adata(sce)
  filename <- paste0(sample,".h5ad")
  filepath <- file.path(outputDir,filename)
  
  if (file.exists(filepath) && !isTRUE(overwrite)) {
    stop(paste0(path, " already exists. Change 'outputDir' or set 'overwrite' to TRUE."))
    }
  
  anndata$write_h5ad(filepath,
                     compression = compression, 
                     compression_opts = compression_opts,
                     as_dense=as_dense,
                     force_dense = force_dense)
}
