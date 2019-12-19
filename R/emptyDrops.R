
.runEmptyDrops <- function(sce,
    sampleColname,
    ...,
    assayType) {

    samples <- unique(SummarizedExperiment::colData(sce)[[sampleColname]])

    for (i in seq_len(length(samples))) {
        sceSampleInd <- which(SummarizedExperiment::colData(sce)
            [[sampleColname]] == samples[i])
        sceSample <- sce[, sceSampleInd]

        if ("DelayedMatrix" %in% class(SummarizedExperiment::assay(sce,
            i = assayType))) {
            SummarizedExperiment::assay(sce, i = assayType) <-
                as.matrix(SummarizedExperiment::assay(sce, i = assayType))
        }

        output <- DropletUtils::emptyDrops(m = SummarizedExperiment::assay(sce,
            i = assayType), ...)
        colnames(output) <- paste0("emptyDrops_", colnames(output))

        if (i == 1) {
            cd <- S4Vectors::DataFrame(row.names = colnames(sce),
                emptyDrops_total = integer(ncol(sce)),
                emptyDrops_logprob = numeric(ncol(sce)),
                emptyDrops_pvalue = numeric(ncol(sce)),
                emptyDrops_limited = logical(ncol(sce)),
                emptyDrops_fdr = numeric(ncol(sce)))
            cd[sceSampleInd, ] <- output
        } else {
            cd[sceSampleInd, ] <- output
        }
    }

    SummarizedExperiment::colData(sce) <-
        cbind(SummarizedExperiment::colData(sce), cd)
    return(sce)
}


#' @title Identify empty droplets using \link[DropletUtils]{emptyDrops}.
#' @description Run \link[DropletUtils]{emptyDrops} on the count matrix in the
#'  provided \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment. This function assumes
#'  the provided
#'  \link[SingleCellExperiment]{SingleCellExperiment} object contains cells
#'  from only one sample and any neccesary subsetting of cells happened
#'  beforehand.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain an unfiltered raw counts matrix.
#' @param sampleColname Character. The column name which specifies the sample
#'  origin in the \link[SummarizedExperiment]{colData} of the provided
#'  \code{sce} object. Default "sample".
#' @param ... Additional arguments to pass to \link[DropletUtils]{emptyDrops}.
#' @param assayType  A string specifying which assay values contain the count
#'  matrix.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{emptyDrops} output table appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
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
    sampleColname = "sample",
    ...,
    assayType = "counts") {

    sce <- .runEmptyDrops(sce = sce,
        sampleColname = sampleColname,
        ...,
        assayType = assayType)
    return(sce)
}
