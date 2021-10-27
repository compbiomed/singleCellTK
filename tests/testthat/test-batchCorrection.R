# Batch Correction Functions
library(singleCellTK)
context("Testing Batch Correction functions")
data(sceBatches, package = "singleCellTK")

sceBatches <- scaterlogNormCounts(sceBatches, "logcounts")

test_that(desc = "Testing Limma Batch Correction", {
  sceBatches <- runLimmaBC(inSCE = sceBatches)
  testthat::expect_true("LIMMA" %in% assayNames(sceBatches))

  # Also plotting function at this point
  p <- plotBatchCorrCompare(sceBatches, "LIMMA")
  testthat::expect_is(p, "gtable")
})

test_that(desc = "Testing ComBat_seq without covariate", {
  sceBatches <- runComBatSeq(inSCE = sceBatches, assayName = "CBS1")

  testthat::expect_true("CBS1" %in% assayNames(sceBatches))
})

test_that(desc = "Testing ComBat_seq with covariate", {
  sceBatches <- runComBatSeq(inSCE = sceBatches, assayName = "CBS2",
                             covariates = "cell_type")

  testthat::expect_true("CBS2" %in% assayNames(sceBatches))
})

test_that(desc = "Testing MNN", {
  sceBatches <- runMNNCorrect(inSCE = sceBatches)

  testthat::expect_true("MNN" %in% assayNames(sceBatches))
})

if (isTRUE(py_available(initialize = FALSE))) {
  test_that(desc = "Testing BBKNN", {
    sceBatches <- runBBKNN(inSCE = sceBatches)
    testthat::expect_true("BBKNN" %in% reducedDimNames(sceBatches))
  })
  test_that(desc = "Testing SCANORAMA", {
    sceBatches <- runSCANORAMA(inSCE = sceBatches)
    testthat::expect_true("SCANORAMA" %in% assayNames(sceBatches))
  })
}
