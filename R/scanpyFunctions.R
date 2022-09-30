# Helper/Wrapper Functions ---

#' .updateAssaySCEFromScanpy
#' Update/Modify/Add an assay in the provided SingleCellExperiment object from
#' an AnnData object
#' @param inSCE Input SingleCellExperiment object
#' @param scanpyObject Input annData object
#' @param assaySlotSCE Selected assay to update in the input
#' SingleCellExperiment object
#' @param scanpyAssaySlot Selected assay from annData object. Default
#' \code{"X"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  data from annData object appended to the \link{assay} slot.
#' @importFrom SummarizedExperiment assay<-
#' @noRd
.updateAssaySCEFromScanpy <- function(inSCE,
                                      scanpyObject,
                                      assaySlotSCE,
                                      scanpyAssaySlot = "X") {
  assay(inSCE, assaySlotSCE) <- NULL
  temp.matrix <- t(scanpyObject[[scanpyAssaySlot]])
  colnames(temp.matrix) <- colnames(inSCE)
  rownames(temp.matrix) <- rownames(inSCE)
  assay(inSCE, assaySlotSCE) <- temp.matrix
  
  return(inSCE)
}

#' runScanpyNormalizeData
#' Wrapper for NormalizeData() function from scanpy library
#' Normalizes the sce object according to the input parameters provided.
#' @param inSCE (sce) object to normalize
#' @param useAssay Assay containing raw counts to use for normalization.
#' @param normAssayName Name of new assay containing normalized data. Default
#' \code{scanpyNormData}.
#' @param normalizationMethod selected normalization method. Default
#' \code{"LogNormalize"}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
runScanpyNormalizeData <- function(inSCE,
                                   useAssay,
                                   normAssayName = "scanpyNormData",
                                   normalizationMethod = "LogNormalize") {
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  if(normalizationMethod == "LogNormalize"){
    # Total-count normalize (library-size correct) to 10,000 reads/cell
    sc$pp$normalize_per_cell(scanpyObject, counts_per_cell_after = 1e4)
    # log transform the data.
    sc$pp$log1p(scanpyObject)
  }
  else{
    scanpyObject <- sc$pp$normalize_total(scanpyObject, 
                                          target_sum=1e4, 
                                          inplace = FALSE)
  }
  
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, normAssayName)
  inSCE@metadata$scanpy$obj <- scanpyObject
  inSCE@metadata$scanpy$normAssay <- normAssayName
  inSCE <- expSetDataTag(inSCE = inSCE,
                         assayType = "normalized",
                         assays = normAssayName)
  return(inSCE)
}


