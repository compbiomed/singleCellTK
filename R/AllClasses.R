#' A lightweight S4 extension to the SingleCellExperiment class to store
#' additional information.
#'
#' @slot pcaVariances The percent variation contained in each PCA dimension
#'
#' @param value The DataFrame of pcaVariances
#'
#' @return A SingleCellExperiment like object with an addition pcaVariances
#' slot.
#' @exportClass SCtkExperiment
#' @examples
#' data("mouseBrainSubsetSCE")
#' counts_mat <- assay(mouseBrainSubsetSCE, "counts")
#' sample_annot <- colData(mouseBrainSubsetSCE)
#' row_annot <- rowData(mouseBrainSubsetSCE)
#' newSCE <- SCtkExperiment(assays=list(counts=counts_mat),
#'                          colData=sample_annot,
#'                          rowData=row_annot)
#' newSCE <- getPCA(newSCE, useAssay = "counts")
#' #View the percent variation of the PCA
#' pcaVariances(newSCE)
#'
setClass("SCtkExperiment",
         slots = c(pcaVariances = "DataFrame"),
         contains = "SingleCellExperiment")

#' Create a SCtkExperiment
#'
#' @param ... SingleCellExperiment and SummarizedExperiment components
#' @param pcaVariances The percent variation contained in each PCA dimension
#'
#' @return A SingleCellExperiment like object with an addition pcaVariances
#' slot.
#'
#' @import SingleCellExperiment SummarizedExperiment
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' counts_mat <- assay(mouseBrainSubsetSCE, "counts")
#' sample_annot <- colData(mouseBrainSubsetSCE)
#' row_annot <- rowData(mouseBrainSubsetSCE)
#' newSCE <- SCtkExperiment(assays=list(counts=counts_mat),
#'                          colData=sample_annot,
#'                          rowData=row_annot)
#' newSCE <- getPCA(newSCE, useAssay = "counts")
#' #View the percent variation of the PCA
#' pcaVariances(newSCE)
#'
SCtkExperiment <- function(..., pcaVariances = S4Vectors::DataFrame()) {
  sce <- SingleCellExperiment::SingleCellExperiment(...)
  out <- methods::new("SCtkExperiment", sce,
                      pcaVariances = S4Vectors::DataFrame())
  return(out)
}
