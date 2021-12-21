#' Example Single Cell RNA-Seq data in SingleCellExperiment Object, GSE60361
#' subset
#'
#' A subset of 30 cells from a single cell RNA-Seq experiment from Zeisel, et
#' al. Science 2015. The data was produced from cells from the mouse
#' somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were
#' identified as oligodendrocytes and 15 of the cell were identified as
#' microglia.
#'
#' @name mouseBrainSubsetSCE
#' @docType data
#' @format SingleCellExperiment
#' @source DOI: 10.1126/science.aaa1934
#' @keywords datasets
#' @usage data("mouseBrainSubsetSCE")
#' @examples
#' data("mouseBrainSubsetSCE")
"mouseBrainSubsetSCE"

#' Example Single Cell RNA-Seq data in SingleCellExperiment Object,
#' subset of 10x public dataset
#' https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
#' A subset of 390 barcodes and top 200 genes were included in this example.
#' Within 390 barcodes, 195 barcodes are empty droplet, 150 barcodes are cell
#' barcode and 45 barcodes are doublets predicted by scrublet and doubletFinder
#' package. This example only serves as a proof of concept and a tutoriol on how
#' to run the functions in this package. The results should not be used for
#' drawing scientific conclusions.

#' @name sce
#' @docType data
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @keywords datasets
#' @usage data("scExample")
#' @examples
#' data("scExample")
"sce"

#' Example Single Cell RNA-Seq data in SingleCellExperiment object, with
#' different batches annotated
#'
#' Two batches of pancreas scRNAseq dataset are combined with their original
#' counts. Cell types and batches are annotated in `colData(sceBatches)`.
#' Two batches came from Wang, et al., 2016, annotated as `'w'`; and Xin, et
#' al., 2016, annotated as `'x'`. Two common cell types, `'alpha'` and
#' `'beta'`, that could be found in both original studies with relatively
#' large population were kept for cleaner demonstration.
#' @usage data('sceBatches')
"sceBatches"

#' Stably Expressed Gene (SEG) list obect, with SEG sets for human and mouse.
#'
#' The two gene sets came from dataset called `segList` of package `scMerge`.
#' @name SEG
#' @docType data
#' @format list, with two entries \code{"human"} and \code{"mouse"}, each is a
#' charactor vector.
#' @source \code{data('segList', package='scMerge')}
#' @keywords datasets
#' @usage data('SEG')
#' @examples
#' data('SEG')
#' humanSEG <- SEG$human
"SEG"

#' MSigDB gene get Category table
#'
#' A table of gene set categories that can be download from MSigDB. The
#' categories and descriptions can be found here:
#' https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp. The IDs in the
#' first column can be used to retrieve the gene sets for these categories
#' using the \link{importGeneSetsFromMSigDB} function.

#' @name msigdb_table
#' @docType data
#' @format A data.frame.
#' @keywords datasets
#' @usage data("msigdb_table")
#' @examples
#' data("msigdb_table")
"msigdb_table"

#' List of mitochondrial genes of multiple reference
#' 
#' A list of gene set that contains mitochondrial genes of multiple reference
#' (hg38, hg19, mm10 and mm9). It contains multiple types of gene identifier:
#' gene symbol, entrez ID, ensemble ID and ensemble transcript ID. It's used 
#' for the function 'importMitoGeneSet'. 

#' @name MitoGenes
#' @docType data
#' @format A list
#' @keywords datasets
#' @usage data("MitoGenes")
#' @examples
#' data("MitoGenes")
"MitoGenes"