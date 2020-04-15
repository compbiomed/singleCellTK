#' @title Detecting contamination with DecontX.
#' @description A wrapper function for \link[celda]{decontX}. Identify
#'  potential contamination from experimental factors such as ambient RNA.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' Default NULL. \link[celda]{decontX} will be run on cells from each
#' sample separately.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#' 'counts'.
#' @param ... Additional arguments to pass to \link[celda]{decontX}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  'decontX_Contamination' and 'decontX_Clusters' added to the
#'  \link[SummarizedExperiment]{colData} slot. Additionally, the
#' decontaminated counts will be added as an assay called 'decontXCounts'.
#' @examples
#' \dontrun{
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runDecontX(emptyDropsSceExample)
#' }
#' @export
runDecontX <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    ...
) {
  argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  }

  message(paste0(date(), " ... Running 'DecontX'"))

  inSCE <- celda::decontX(x = inSCE, batch = sample, assayName = useAssay, ...)
  argsList = argsList[!names(argsList) %in% ("...")]
  inSCE@metadata$runDecontX <- argsList[-1]
  inSCE@metadata$runDecontX$packageVersion <- packageDescription("celda")$Version

  return(inSCE)
}

