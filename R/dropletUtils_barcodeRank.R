
.runBarcodeRankDrops <- function(barcode.matrix, ...) {

  if (class(barcode.matrix) != "dgCMatrix") {
    barcode.matrix <- as(barcode.matrix, "dgCMatrix")
  }

  output <- DropletUtils::barcodeRanks(m = barcode.matrix, ...)

  knee.ix <- as.integer(output@listData$total >= metadata(output)$knee)
  inflection.ix <- as.integer(output@listData$total >= metadata(output)$inflection)

  result <- cbind(knee.ix, inflection.ix)
  colnames(result) <- c("dropletUtils_BarcodeRank_Knee",
                        "dropletUtils_BarcodeRank_Knee")

  return(result)
}


#' @title Identify empty droplets using \link[DropletUtils]{barcodeRanks}.
#' @description Run \link[DropletUtils]{barcodeRanks} on a count matrix 
#'  provided in a\link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment. 
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain a raw counts matrix before empty droplets have been removed.
#' @param sample Character vector. Indicates which sample each cell belongs to
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample separately.
#'  If NULL, then all cells will be processed together. Default NULL.
#' @param ... Additional arguments to pass to \link[DropletUtils]{barcodeRanks}.
#' @param assayName  A string specifying which assay in the SCE to use.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{barcodeRanks} output table appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
#'  \emph{dropletUtils_BarcodeRank_Knee} and \emph{dropletUtils_BarcodeRank_Knee}
#'  Please refer to the documentation of \link[DropletUtils]{barcodeRanks} for
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
#' sce <- runBarcodeRankDrops(sce = emptyDropsSceExample)
#' @export
runBarcodeRankDrops <- function(sce,
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
  
  message(paste0(date(), " ... Running 'barcodeRanks'"))    
  
  ## Define result matrix for all samples  
  output <- S4Vectors::DataFrame(row.names = colnames(sce),
            dropletUtils_BarcodeRank_Knee = integer(ncol(sce)),
            dropletUtils_BarcodeRank_Inflection = integer(ncol(sce)))
            
  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
    sceSample <- sce[, sceSampleInd]

    mat <- SummarizedExperiment::assay(sceSample, i = assayName)  
    result <- .runBarcodeRankDrops(barcode.matrix = mat, ...)
      
    output[sceSampleInd, ] <- result
  }
  
  colData(sce) = cbind(colData(sce), output)
  
  return(sce)
}
