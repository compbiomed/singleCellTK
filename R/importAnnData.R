.importAnnDataSample <- function(sampleDir = './',
                                 sampleName = 'sample',
                                 delayedArray = TRUE){

  anndata_file <- file.path(sampleDir, paste0(sampleName,'.h5ad',sep=''))
  if (!file.exists(anndata_file)){
    stop("AnnData file not found at specified location. Please check path provided and/or filename.")
  }
  anndata <- ad$read_h5ad(anndata_file)

  counts_matrix <- t((reticulate::py_to_r(anndata$X)))
  if (isTRUE(delayedArray)) {
    counts_matrix <- DelayedArray::DelayedArray(counts_matrix)
  }

  sce_rowdata <- S4Vectors::DataFrame(reticulate::py_to_r(anndata$var))
  sce_coldata <- S4Vectors::DataFrame(reticulate::py_to_r(anndata$obs))
  sce <- SingleCellExperiment(list(counts = counts_matrix),
                              rowData = sce_rowdata,
                              colData = sce_coldata)
  colnames(sce) <- paste0(sampleName,"_",colnames(sce))

  multi_Assay <- reticulate::py_to_r(anndata$layers$as_dict())
  for(assay_name in names(multi_Assay)){
    tryCatch({
      SummarizedExperiment::assay(sce, assay_name, withDimnames = FALSE) <- t(reticulate::py_to_r(multi_Assay[[assay_name]]))
      base::dimnames(SummarizedExperiment::assay(sce, assay_name)) <- base::dimnames(SummarizedExperiment::assay(sce, "counts"))
    }, error = function(x){
      error_message <- paste0("Warning: unable to add '",assay_name,"' from .layers AnnData slot to SCE assay. Skipping. ")
      message(error_message)
    })
  }
  
  multidim_observations <- reticulate::py_to_r(anndata$obsm_keys())
  for(obsm_name in multidim_observations){
    tryCatch({
      SingleCellExperiment::reducedDims(sce)[[obsm_name]] <- reticulate::py_to_r(anndata$obsm[obsm_name])
    }, error = function(x){
      error_message <- paste0("Warning: unable to add '",obsm_name,"' from .obsm AnnData slot to SCE metadata. Skipping. ")
      message(error_message)
    })
  }

  unstructured_data <- reticulate::py_to_r(anndata$uns_keys())
  for(uns_name in unstructured_data){
    tryCatch({
      sce@metadata[[sampleName]]$annData[[uns_name]] <- reticulate::py_to_r(anndata$uns[uns_name])
    }, error = function(x){
      error_message <- paste0("Warning: unable to add unstructured data (.uns slot): '",uns_name,"' to SCE metadata. Skipping. ")
      message(error_message)
    })

  }

  return(sce)

}

#' @name importAnnData
#' @rdname importAnnData
#' @title Create a SingleCellExperiment Object from Python AnnData .h5ad files
#' @description This function reads in one or more Python AnnData files in the .h5ad format
#' and returns a single \link[SingleCellExperiment]{SingleCellExperiment} object containing all the
#' AnnData samples by concatenating their counts matrices and related information slots.
#' @param sampleDirs Folder containing the .h5ad file. Can be one of -
#' \itemize{
#'   \item Default \code{current working directory}.
#'   \item Full path to the directory containing the .h5ad file.
#'   E.g \code{sampleDirs = '/path/to/sample'}
#'   \item A vector of folder paths for the samples to import.
#'   E.g. \code{sampleDirs = c('/path/to/sample1', '/path/to/sample2','/path/to/sample3')}
#'   importAnnData will return a single SCE object containing all the samples
#'   with the sample name appended to each colname in colData
#' }
#' @param sampleNames The prefix/name of the .h5ad file without the .h5ad extension
#' e.g. if 'sample.h5ad' is the filename, pass \code{sampleNames = 'sample'}.
#' Can be one of -
#' \itemize{
#'   \item Default \code{sample}.
#'   \item A vector of samples to import. Length of vector must be equal to length of sampleDirs vector
#'   E.g. \code{sampleDirs = c('sample1', 'sample2','sample3')}
#'   importAnnData will return a single SCE object containing all the samples
#'   with the sample name appended to each colname in colData
#' }
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object. Default \code{TRUE}.
#' @details
#' \code{importAnnData} converts scRNA-seq data in the AnnData format to the
#' \code{SingleCellExperiment} object. The .X slot in AnnData is transposed to the features x cells
#' format and becomes the 'counts' matrix in the assay slot. The .vars AnnData slot becomes the SCE rowData
#' and the .obs AnnData slot becomes the SCE colData. Multidimensional data in the .obsm AnnData slot is
#' ported over to the SCE reducedDims slot. Additionally, unstructured data in the .uns AnnData slot is
#' available through the SCE metadata slot.
#' There are 2 currently known minor issues -
#' Anndata python module depends on another python module h5pyto read hd5 format files.
#' If there are errors reading the .h5ad files, such as "ValueError: invalid shape in fixed-type tuple."
#' the user will need to do downgrade h5py by running \code{pip3 install --user h5py==2.9.0}
#' Additionally there might be errors in converting some python objects in the unstructured data slots.
#' There are no known R solutions at present. Refer \url{https://github.com/rstudio/reticulate/issues/209}
#' @return A \code{SingleCellExperiment} object.
#' @examples
#' file.path <- system.file("extdata/annData_pbmc_3k", package = "singleCellTK")
#' \dontrun{
#' sce <- importAnnData(sampleDirs = file.path,
#'                      sampleNames = 'pbmc3k_20by20')
#' }
#' @export
importAnnData <- function(sampleDirs = NULL,
                          sampleNames = NULL,
                          delayedArray = TRUE) {

  if (length(sampleDirs)!=length(sampleNames)){
    stop("Number of sampleDirs must be equal to number of SampleNames. Please provide sample names for all input directories")
  }

  res <- vector("list", length = length(sampleDirs))

  for (i in seq_along(sampleDirs)){
    scei <- .importAnnDataSample(sampleDir = sampleDirs[[i]],
                                 sampleName = sampleNames[[i]],
                                 delayedArray = delayedArray)
    res[[i]] <- scei
  }
  sce <- do.call(SingleCellExperiment::cbind, res)
  return(sce)
}













