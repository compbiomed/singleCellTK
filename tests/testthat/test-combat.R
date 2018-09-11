context("ComBat")

test_that("ComBat", {
  library(bladderbatch)
  data(bladderdata)
  #subset for testing
  dat <- bladderEset[1:50, ]
  dat <- as(as(dat, "SummarizedExperiment"), "SCtkExperiment")
  mod <- stats::model.matrix(~as.factor(cancer), data = colData(dat))
  # parametric adjustment
  expect_is(ComBatSCE(inSCE = dat, useAssay = "exprs", batch = "batch",
                      covariates = NULL), "matrix")

  # non-parametric adjustment, mean-only version
  expect_is(ComBatSCE(inSCE = dat, useAssay = "exprs", batch = "batch",
                      par.prior = "Non-parametric", mean.only = TRUE,
                      covariates = NULL), "matrix")

  expect_is(ComBatSCE(inSCE = dat, useAssay = "exprs", batch = "batch",
                      covariates = "cancer", ref.batch = 3), "matrix")

  expect_error(ComBatSCE(inSCE = dat, useAssay = "exprs", batch = "batch",
                         covariates = "cancer", par.prior = "NON-parametric",
                         ref.batch = 3), "Invalid option given to par.prior. Accepted values are Parametric and Non-parametric.")
})
