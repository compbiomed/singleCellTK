#' ComBat_SCE
#'
#' Run ComBat on a SCtkExperiment object
#'
#' @param SCEdata SCtkExperiment object. Required
#' @param batch The name of a column in colData to use as the batch variable.
#' Required
#' @param use_assay The assay to use for ComBat. The default is "logcounts"
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
#'
#' @export
ComBat_SCE <- function(SCEdata, batch, use_assay="logcounts",
                       par.prior="Parametric", covariates=NULL, mean.only=FALSE,
                       ref.batch=NULL){

  #prepare model matrix
  mod <- NULL
  if (length(covariates) > 0){
    mod <- model.matrix(as.formula(paste0("~", paste0(covariates,
                                                      collapse = "+"))),
                        data = data.frame(colData(SCEdata)[, covariates, drop = FALSE]))
  }

  #prepare parametric
  if (par.prior == "Parametric"){
    par.prior <- TRUE
  } else {
    par.prior <- FALSE
  }

  resassay <-
    sva::ComBat(dat = assay(SCEdata, use_assay),
                batch = colData(SCEdata)[, batch],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)
  return(resassay)
}
