#' @title Find doublets using \code{scrublet}.
#' @description A wrapper function that calls \code{scrub_doublets} from python
#'  module \code{scrublet}. Simulates doublets from the observed data and uses
#'  a k-nearest-neighbor classifier to calculate a continuous
#'  \code{scrublet_score} (between 0 and 1) for each transcriptome. The score
#'  is automatically thresholded to generate \code{scrublet_call}, a boolean
#'  array that is \code{TRUE} for predicted doublets and \code{FALSE} otherwise.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  Scrublet will be run on cells from each sample separately. If NULL, then
#'  all cells will be processed together. Default \code{NULL}.
#' @param assayName  A string specifying which assay in the SCE to use. Default 'counts'.
#' @param seed Seed for the random number generator. Default 12345.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \code{scrub_doublets} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{scrublet_score} and \emph{scrublet_call}.
#' @examples
#' \dontrun{
#' data(sce_chcl, package = "scds")
#' sce <- runScrublet(sce_chcl)
#' }
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
runScrublet <- function(sce,
  sample = NULL,
  assayName = "counts",
  seed = 12345) {

  if (!reticulate::py_module_available(module = "scrublet")) {
    warning("Cannot find python module 'scrublet', please install through pip (e.g. pip install scrublet)
            or use 'use_python()' to select correct Python environment.")
    return(sce)
  }

  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }

  if (!is.null(sample)) {
    if (length(sample) != ncol(sce)) {
      stop("'sample' must be the same length as the number of",
        " columns in 'sce'")
    }
  } else {
    sample = rep(1, ncol(sce))
  }

  message(paste0(date(), " ... Running 'scrublet'"))

  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(row.names = colnames(sce),
    scrublet_score = numeric(ncol(sce)),
    scrublet_call = logical(ncol(sce)))

  ## Loop through each sample and run scrublet
  error <- try({
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
      sceSampleInd <- sample == samples[i]
      sceSample <- sce[, sceSampleInd]

      mat <- SummarizedExperiment::assay(sceSample, i = assayName)

      if (class(mat) != "dgCMatrix") {
        mat <- methods::as(mat, "dgCMatrix")
      }

      scr <- scrublet$Scrublet(t(mat))
      result <- scr$scrub_doublets()

      output[sceSampleInd, "scrublet_score"] <- result[[1]]
      output[sceSampleInd, "scrublet_call"] <- result[[2]]
    }

    colData(sce) = cbind(colData(sce), output)
  }, silent = TRUE)

  if(inherits(error, "try-error")) {
    warning("Scrublet did not complete successfully. Returning SCE without making any changes.")
  }

  return(sce)
}

