#' runComBat
#'
#' Run ComBat batch correction method on a SCtkExperiment object
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object. Required
#' @param batch A character scalar, for the name of a column in
#' \code{\link[SummarizedExperiment]{colData}} to use as the batch variable.
#' Default \code{"batch"}.
#' @param useAssay The assay to correct. Default \code{"logcounts"}.
#' @param par.prior A logical scalar. TRUE indicates parametric adjustments
#' will be used, FALSE indicates non-parametric adjustments will be used.
#' Default \code{TRUE}.
#' @param covariates List of other column names in colData to be added to the
#' ComBat model as covariates. Default \code{NULL}.
#' @param mean.only If TRUE ComBat only corrects the mean of the batch effect.
#' Default \code{FALSE}.
#' @param ref.batch If given, will use the selected batch as a reference for
#' batch adjustment. Default \code{NULL}.
#' @param assayName A character scalar for the new assay name to save the
#' corrected expression. Default \code{"ComBat"}
#'
#' @return Input SCTK object with corrected assay updated at \code{assay(inSCE, assayName)}.
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'
#'   #subset for testing
#'   dat <- bladderEset[1:50,]
#'   dat <- as(as(dat, "SummarizedExperiment"), "SCtkExperiment")
#'
#'   # parametric adjustment
#'   dat <- runComBat(inSCE = dat, useAssay = "exprs",
#'                    batch = "batch", covariates = NULL,
#'                    assayName = "parametric_combat")
#'
#'   # non-parametric adjustment, mean-only version
#'   dat <- runComBat(inSCE = dat, useAssay = "exprs",
#'                    batch = "batch", par.prior = FALSE,
#'                    mean.only = TRUE, covariates = NULL,
#'                    assayName = "nonparametric_combat_meanonly")
#'
#'   # reference-batch version, with covariates
#'   dat <- runComBat(inSCE = dat, useAssay = "exprs",
#'                    batch = "batch", covariates = "cancer",
#'                    ref.batch = 3, assayName = "refbatch_combat_wcov")
#' }
#' @export
runComBat <- function(inSCE, batch = 'batch', useAssay = "logcounts",
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
    stop(paste("\"annotation\" name:", batch, "not found"))
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
                batch = SummarizedExperiment::colData(inSCE)[, batch],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)

  assay(inSCE, assayName) <- resassay
  return(inSCE)
}
