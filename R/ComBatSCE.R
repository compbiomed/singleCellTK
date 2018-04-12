#' ComBatSCE
#'
#' Run ComBat on a SCtkExperiment object
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param batch The name of a column in colData to use as the batch variable.
#' Required
#' @param useAssay The assay to use for ComBat. The default is "logcounts"
#' @param par.prior TRUE indicates parametric adjustments will be used, FALSE
#' indicates non-parametric adjustments will be used. Accepted parameters:
#' "Parametric" or "Non-parametric"
#' @param covariates List of other column names in colData to be added to the
#' ComBat model as covariates
#' @param mean.only If TRUE ComBat only corrects the mean of the batch effect
#' @param ref.batch If given, will use the selected batch as a reference for
#' batch adjustment.
#'
#' @return ComBat matrix based on inputs. You can save this matrix into the
#' SCtkExperiment with assay()
#' @export
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'
#'   #subset for testing
#'   dat <- bladderEset[1:50,]
#'   dat <- as(as(dat, "SummarizedExperiment"), "SCtkExperiment")
#'   mod <- stats::model.matrix(~as.factor(cancer), data = colData(dat))
#'
#'   # parametric adjustment
#'   combat_edata1 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", covariates = NULL)
#'   assay(dat, "parametric_combat") <- combat_edata1
#'
#'   # non-parametric adjustment, mean-only version
#'   combat_edata2 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", par.prior = "Non-parametric",
#'                              mean.only = TRUE, covariates = NULL)
#'   assay(dat, "nonparametric_combat_meanonly") <- combat_edata2
#'
#'   # reference-batch version, with covariates
#'   combat_edata3 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", covariates = "cancer",
#'                              ref.batch = 3)
#'   assay(dat, "refbatch_combat_wcov") <- combat_edata3
#'   assays(dat)
#' }
#'
ComBatSCE <- function(inSCE, batch, useAssay="logcounts",
                      par.prior="Parametric", covariates=NULL, mean.only=FALSE,
                      ref.batch=NULL){

  #prepare model matrix
  mod <- NULL
  if (length(covariates) > 0){
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCE)[, covariates,
                                                               drop = FALSE]))
  }

  #prepare parametric
  if (par.prior == "Parametric"){
    par.prior <- TRUE
  } else {
    par.prior <- FALSE
  }

  resassay <-
    sva::ComBat(dat = SummarizedExperiment::assay(inSCE, useAssay),
                batch = SingleCellExperiment::colData(inSCE)[, batch],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)
  return(resassay)
}
