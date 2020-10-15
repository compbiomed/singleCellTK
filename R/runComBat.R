#' Apply ComBat batch effect correction method to SingleCellExperiment object
#'
#' The ComBat batch adjustment approach assumes that batch effects represent
#' non-biological but systematic shifts in the mean or variability of genomic
#' features for all samples within a processing batch. It uses either parametric
#' or non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param par.prior A logical scalar. TRUE indicates parametric adjustments
#' will be used, FALSE indicates non-parametric adjustments will be used.
#' Default \code{TRUE}.
#' @param covariates List of other column names in colData to be added to the
#' ComBat model as covariates. Default \code{NULL}.
#' @param mean.only If TRUE ComBat only corrects the mean of the batch effect.
#' Default \code{FALSE}.
#' @param ref.batch If given, will use the selected batch as a reference for
#' batch adjustment. Default \code{NULL}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"ComBat"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' \dontrun{
#' # parametric adjustment
#' sceCorr <- runComBat(sceBatches)
#' # non-parametric adjustment, mean-only version
#' sceCorr <- runComBat(sceBatches, par.prior=FALSE, mean.only=TRUE)
#' # reference-batch version, with covariates
#' sceCorr <- runComBat(sceBatches, covariates = "cell_type", ref.batch = 'w')
#' }
#' @export
runComBat <- function(inSCE, useAssay = "logcounts", batch = 'batch',
                      par.prior = TRUE, covariates = NULL,
                      mean.only = FALSE, ref.batch = NULL,
                      assayName = "ComBat") {
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("\"inSCE\" should be a SingleCellExperiment Object.")
  }
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
    stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
  }
  if(any(!c(batch, covariates) %in% names(SummarizedExperiment::colData(inSCE)))){
    anns <- c(batch, covariates)
    notFound <- which(!anns %in% names(SummarizedExperiment::colData(inSCE)))
    notFound <- anns[notFound]
    stop("\"annotation\" name:", paste(notFound, collapse = ', '), "not found")
  }
  #prepare model matrix
  mod <- NULL
  if (!is.null(covariates)){
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SummarizedExperiment::colData(inSCE)))
  }

  resassay <-
    sva::ComBat(dat = SummarizedExperiment::assay(inSCE, useAssay),
                batch = SummarizedExperiment::colData(inSCE)[[batch]],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)

  SummarizedExperiment::assay(inSCE, assayName) <- resassay
  return(inSCE)
}
