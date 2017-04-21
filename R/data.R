#' Example Single Cell RNA-Seq data in SCESet Object, GSE60361 subest
#'
#' A subset of 30 samples from a single cell RNA-Seq experiment from Zeisel, et
#' al. Science 2015. The data was produced from cells from the mouse
#' somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were
#' identified as oligodendrocytes and 15 of the cell were identified as
#' microglia.
#'
#' @name GSE60361_subset
#' @docType data
#' @format List of two data frames, with counts and annotations. Use them as
#' input to createSCESet()
#' @source DOI: 10.1126/science.aaa1934
#' @keywords datasets
#' @examples
#' library(scater)
#' data("GSE60361_subset")
#' GSE60361_SCESet <- createSCESet(countfile = GSE60361_subset$counts,
#'                                 annotfile = GSE60361_subset$annot,
#'                                 inputdataframes = TRUE)
"GSE60361_subset"

#' Example Single Cell RNA-Seq data in SCESet Object, GSE73121
#'
#' 117 Single-cell transcriptome profiling for metastatic renal cell carcinoma
#' patient-derived cells
#'
#' @name GSE73121
#' @docType data
#' @format List of two data frames, with counts and annotations. Use them as
#' input to createSCESet()
#' @source DOI: 10.1186/s13059-016-0945-9
#' @keywords datasets
#' @examples
#' library(scater)
#' data("GSE73121")
#' GSE73121_SCESet <- createSCESet(countfile = GSE73121$counts,
#'                                 annotfile = GSE73121$annot,
#'                                 inputdataframes = TRUE)
"GSE73121"

#' Example Single Cell RNA-Seq data in SCESet Object, GSE36552
#'
#' 86 Single-cell transcriptome profiling from embryonic stem cells derived from
#' donated human pre-implatation embryos.
#'
#' @name GSE36552
#' @docType data
#' @format List of two data frames, with counts and annotations. Use them as
#' input to createSCESet()
#' @source DOI: 10.1038/nsmb.2660
#' @keywords datasets
#' @examples
#' library(scater)
#' data("GSE36552")
#' GSE36552_SCESet <- createSCESet(countfile = GSE36552$counts,
#'                                 annotfile = GSE36552$annot,
#'                                 inputdataframes = TRUE)
"GSE36552"
