#' @name importExampleData
#' @title Retrieve example datasets
#' @description Retrieves published example datasets stored in
#' \link[SingleCellExperiment]{SingleCellExperiment} using the
#' \link[scRNAseq]{scRNAseq} and
#' \link[TENxPBMCData]{TENxPBMCData} packages. See 'Details' for a
#' list of available datasets.
#' @param dataset Character. Name of the dataset to retrieve.
#' @param class Character. The class of the expression matrix stored in the SCE
#'   object. Can be one of \code{"Matrix"} or \code{"matrix"}. \code{"Matrix"}
#'   will store the data as a sparse matrix from package \link{Matrix} while
#'   \code{"matrix"} will store the data in a standard matrix. Default
#'   \code{"Matrix"}.
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'   \link[DelayedArray]{DelayedArray} object or not. Default \code{FALSE}.
#' @details See the list below for the available datasets and their
#' descriptions.
#' \describe{
#' \item{"fluidigm_pollen"}{Retrieved with
#' \code{\link[scRNAseq]{ReprocessedFluidigmData}}. Returns a dataset of 65
#'  human neural cells from Pollen et al. (2014), each sequenced at high and low
#'  coverage (SRA accession SRP041736).}
#' \item{"allen_tasic"}{Retrieved with
#' \code{\link[scRNAseq]{ReprocessedAllenData}}. Returns a dataset of 379 mouse
#' brain cells from Tasic et al. (2016).}
#' \item{"pbmc3k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  2,700 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' \item{"pbmc4k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  4,340 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' \item{"pbmc6k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  5,419 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' \item{"pbmc8k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  8,381 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' \item{"pbmc33k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  33,148 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' \item{"pbmc68k"}{Retrieved with \code{\link[TENxPBMCData]{TENxPBMCData}}.
#'  68,579 peripheral blood mononuclear cells (PBMCs) from 10X Genomics.}
#' }
#' @author Joshua D. Campbell, David Jenkins
#' @return The specified \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @examples
#' sce <- importExampleData("pbmc3k")
#' @export
#' @importFrom SummarizedExperiment colData rowData colData<- assay assays
importExampleData <- function(dataset, class = c("Matrix", "matrix"),
                              delayedArray = FALSE) {
  class <- match.arg(class)

  scRNAseqDatasets <- c("fluidigm_pollen", "allen_tasic")
  tenxPbmcDatasets <- c("pbmc3k", "pbmc4k", "pbmc6k", "pbmc8k", "pbmc33k", "pbmc68k")

  if(dataset %in% scRNAseqDatasets) {
    if(!("scRNAseq" %in% rownames(utils::installed.packages()))) {
      stop(paste0("Package 'scRNAseq' is not installed. Please install to load dataset '", dataset, "'."))
    }
    if (dataset == "fluidigm_pollen") {
      temp <- scRNAseq::ReprocessedFluidigmData()
      temp$sample <- paste0(colData(temp)$Biological_Condition, "_", colData(temp)$Coverage_Type)
    } else if (dataset == "allen_tasic") {
      temp <- scRNAseq::ReprocessedAllenData()
      temp$sample <- paste0(colData(temp)$driver_1_s, "_", colData(temp)$dissection_s)
    }
  } else if (dataset %in% tenxPbmcDatasets) {
    if(!("TENxPBMCData" %in% rownames(utils::installed.packages()))) {
      stop(paste0("Package 'TENxPBMCData' is not installed. Please install to load dataset '", dataset, "'."))
    }
    temp <- TENxPBMCData::TENxPBMCData(dataset = dataset)
    colnames(temp) <- paste(temp$Sample, temp$Barcode, sep="_")
    rownames(temp) <- rowData(temp)$Symbol_TENx
    colData(temp)$sample <- colData(temp)$Sample
  } else {
    stop("'dataset' must be one of: ", paste(c(scRNAseqDatasets, tenxPbmcDatasets), collapse = ","))
  }

  # Convert to sparseMatrix or regular matrix
  for(i in seq_along(names(assays(temp)))) {
    if (class == "matrix") {
      assay(temp, i) <- as.matrix(temp)
    } else if(class == "Matrix") {
      assay(temp, i) <- .convertToMatrix(assay(temp, i))
    }

    if (isTRUE(delayedArray)) {
      assay(temp, i) <- DelayedArray::DelayedArray(assay(temp, i))
    }
  }
  return(temp)
}


