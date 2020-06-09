#' Example Single Cell RNA-Seq data in SCtkExperiment Object, GSE60361
#' subset
#'
#' A subset of 30 samples from a single cell RNA-Seq experiment from Zeisel, et
#' al. Science 2015. The data was produced from cells from the mouse
#' somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were
#' identified as oligodendrocytes and 15 of the cell were identified as
#' microglia.
#'
#' @name mouseBrainSubsetSCE
#' @docType data
#' @format SCtkExperiment
#' @source DOI: 10.1126/science.aaa1934
#' @keywords datasets
#' @examples
#' data("mouseBrainSubsetSCE")
"mouseBrainSubsetSCE"

#' Example Single Cell RNA-Seq data in SingleCellExperiment Object, 
#' subset of 10x public dataset
#' https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
#' A subset of 390 barcodes and top 200 genes were included in this example.
#' Within 390 barcodes, 195 barcodes are empty droplet, 150 barcodes are cell barcode
#' and 45 barcodes are doublets predicted by scrublet and doubletFinder package. 
#' This example only serves as a proof of concept and a tutoriol on how to
#' run the functions in this package. The results should not be
#' used for drawing scientific conclusions.

#' @name sce
#' @docType data
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @keywords datasets
#' @examples
#' data("sceQCExample")
"sce"

#' Example PBMC_1k_v3_33538x20 SingleCellExperiment Object
#'
#' The following unfiltered PBMC_1k_v3 data were downloaded from
#' https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' /pbmc_1k_v3
#' Only the top 10 cells with most counts and the last 10 cells with non-zero
#' counts are included in this example.
#' This example only serves as a proof of concept and a tutoriol on how to
#' run the functions in this package. The results should not be
#' used for drawing scientific conclusions.
#' @examples
#' data("emptyDropsSceExample", package = "singleCellTK")
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object.
"emptyDropsSceExample"

#' Example Single Cell RNA-Seq data in SingleCellExperiment object, with
#' different batches annotated
#'
#' Two batches of pancreas scRNAseq dataset are combined with their original
#' counts. Cell types and batches are annotated in `colData(sceBatches)`.
#' Two batches came from Wang, et al., 2016, annotated as `'w'`; and Xin, et
#' al., 2016, annotated as `'x'`. Two common cell types, `'alpha'` and 
#' `'beta'`, that could be found in both original studies with relatively 
#' large population were kept for cleaner demonstration.
#'
#' @name sceBatches
#' @docType data
#' @format SingleCellExperiment
#' @source DOI: 10.2337/db16-0405 and 10.1016/j.cmet.2016.08.018
#' @keywords datasets
#' @examples
#' data('sceBatches')
"sceBatches"

#' Stably Expressed Gene (SEG) list obect, with SEG sets for human and mouse.
#' 
#' The two gene sets came from dataset called `segList` of package `scMerge`.
#' @name SEG
#' @docType data
#' @format list, with two entries `"human"` and `"mouse"`, each is a charactor
#' array.
#' @source `data('segList', package='scMerge')``
#' @keywords datasets
#' @examples 
#' data('SEG')
#' humanSEG <- SEG$human
"SEG"