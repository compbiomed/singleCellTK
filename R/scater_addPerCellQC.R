#' @title Wrapper for calculating QC metrics with scater.
#' @description A wrapper function for \link[celda]{decontX}. Identify
#'  potential contamination from experimental factors such as ambient RNA.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object. 
#' @param assayName  A string specifying which assay in the SCE to use. Default 
#' 'counts'.
#' @param ... Additional arguments to pass to \link[celda]{decontX}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  cell QC metrics added to the \link[SummarizedExperiment]{colData} slot. 
#' @examples
#' data(emptyDropsSceExample, package = "singleCellTK")
#' mito.ix = grep("^MT-", rowData(emptyDropsSceExample)$feature_name)
#' geneSet <- list("Mito"=rownames(emptyDropsSceExample)[mito.ix])
#' sce <- runPerCellQC(emptyDropsSceExample, geneSet = geneSet)
#' @export
runPerCellQC <- function(sce,
    assayName = "counts",
    geneSets = NULL,
    ...
) {

  message(paste0(date(), " ... Running 'perCellQCMetrics'"))    
  
  ## TODO: Add error checking of gene sets
  
  sce <- addPerCellQC(x = sce, exprs_values = assayName, subsets = geneSets, ...)
  
  return(sce)
}

