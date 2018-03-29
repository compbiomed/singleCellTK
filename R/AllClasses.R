#' A lightweight S4 extension to the SingleCellExperiment class to store
#' additional information.
#'
#' @slot pca_variances The percent variation contained in each PCA dimension
#'
#' @param value The DataFrame of pca_variances
#'
#' @return A SingleCellExperiment like object with an addition pca_variances
#' slot.
#' @exportClass SCtkExperiment
#' @examples
#' data("mouse_brain_subset_sce")
#' counts_mat <- assay(mouse_brain_subset_sce, "counts")
#' sample_annot <- colData(mouse_brain_subset_sce)
#' row_annot <- rowData(mouse_brain_subset_sce)
#' newSCE <- SCtkExperiment(assays=list(counts=counts_mat),
#'                          colData=sample_annot,
#'                          rowData=row_annot)
#' newSCE <- getPCA(newSCE, use_assay = "counts")
#' #View the percent variation of the PCA
#' pca_variances(newSCE)
#'
setClass("SCtkExperiment",
         slots = c(pca_variances = "DataFrame"),
         contains = "SingleCellExperiment")

#' Create a SCtkExperiment
#'
#' @param ... SingleCellExperiment and SummarizedExperiment components
#' @param pca_variances The percent variation contained in each PCA dimension
#'
#' @return A SingleCellExperiment like object with an addition pca_variances
#' slot.
#'
#' @import SingleCellExperiment SummarizedExperiment
#'
#' @export
#' @examples
#' data("mouse_brain_subset_sce")
#' counts_mat <- assay(mouse_brain_subset_sce, "counts")
#' sample_annot <- colData(mouse_brain_subset_sce)
#' row_annot <- rowData(mouse_brain_subset_sce)
#' newSCE <- SCtkExperiment(assays=list(counts=counts_mat),
#'                          colData=sample_annot,
#'                          rowData=row_annot)
#' newSCE <- getPCA(newSCE, use_assay = "counts")
#' #View the percent variation of the PCA
#' pca_variances(newSCE)
#'
SCtkExperiment <- function(..., pca_variances = S4Vectors::DataFrame()) {
  sce <- SingleCellExperiment::SingleCellExperiment(...)
  out <- methods::new("SCtkExperiment", sce, pca_variances = S4Vectors::DataFrame())
  return(out)
}
