.runEmptyDrops <- function(barcode.matrix, ...) {

  if (class(barcode.matrix) != "dgCMatrix") {
    barcode.matrix <- as(barcode.matrix, "dgCMatrix")
  }

  result <- DropletUtils::emptyDrops(m = barcode.matrix, ...)
  colnames(result) <- paste0("dropletUtils_emptyDrops_", colnames(result))

  return(result)
}


#' @title Identify empty droplets using \link[DropletUtils]{emptyDrops}.
#' @description Run \link[DropletUtils]{emptyDrops} on the count matrix in the
#'  provided \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain a raw counts matrix before empty droplets have been removed.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample separately.
#'  If NULL, then all cells will be processed together. Default NULL.
#' @param assayName  A string specifying which assay in the SCE to use.
#' @param ... Additional arguments to pass to \link[DropletUtils]{emptyDrops}.
#'  matrix.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{emptyDrops} output table appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{emptyDrops_total}, \emph{emptyDrops_logprob},
#'  \emph{emptyDrops_pvalue}, \emph{emptyDrops_limited}, \emph{emptyDrops_fdr}.
#'  Please refer to the documentation of \link[DropletUtils]{emptyDrops} for
#'  details.
#' @examples
#' # The following unfiltered PBMC_1k_v3 data were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' # /pbmc_1k_v3
#' # Only the top 10 cells with most counts and the last 10 cells with non-zero
#' # counts are included in this example.
#' # This example only serves as an proof of concept and a tutoriol on how to
#' # run the function. The results should not be
#' # used for drawing scientific conclusions.
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runEmptyDrops(sce = emptyDropsSceExample)
#' @import DropletUtils
#' @export
runEmptyDrops <- function(sce,
    sample = NULL,
    assayName = "counts",
    ...
) {
  if(!is.null(sample)) {
    if(length(sample) != ncol(sce)) {
      stop("'sample' must be the same length as the number of columns in 'sce'")
    }
  } else {
    sample = rep(1, ncol(sce))
  }

  message(date(), " ... Running 'emptyDrops'")

  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(row.names = colnames(sce),
                dropletUtils_emptyDrops_total = integer(ncol(sce)),
                dropletUtils_emptyDrops_logprob = numeric(ncol(sce)),
                dropletUtils_emptyDrops_pvalue = numeric(ncol(sce)),
                dropletUtils_emptyDrops_limited = logical(ncol(sce)),
                dropletUtils_emptyDrops_fdr = numeric(ncol(sce)))

  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
    sceSample <- sce[, sceSampleInd]

    mat <- SummarizedExperiment::assay(sceSample, i = assayName)
    result <- .runEmptyDrops(barcode.matrix = mat, ...)

    output[sceSampleInd, ] <- result
  }

  colData(sce) = cbind(colData(sce), output)

  return(sce)
}
