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

test_that(desc = "Testing Harmony", {
  testthat::skip_if_not_installed("harmony")
  sceBatches <- scaterPCA(sceBatches, useAssay = "logcounts",
                          reducedDimName = "PCA", nComponents = 10,
                          useFeatureSubset = NULL)

  # from a reducedDim
  sceRD <- runHarmony(inSCE = sceBatches, useReducedDim = "PCA",
                      nComponents = 10, verbose = FALSE)
  testthat::expect_true("HARMONY" %in% reducedDimNames(sceRD))
  hRD <- SingleCellExperiment::reducedDim(sceRD, "HARMONY")
  testthat::expect_equal(nrow(hRD), ncol(sceBatches))
  testthat::expect_equal(ncol(hRD), 10)
  testthat::expect_true(all(is.finite(hRD)))
  testthat::expect_equal(metadata(sceRD)$batchCorr$HARMONY$method, "harmony")

  # from a full-size assay: PCA is computed internally for harmony >= 1.0.0.
  # The assay input is expected to warn; pin that warning rather than
  # suppressing all of them.
  testthat::expect_warning(
    sceAssay <- runHarmony(inSCE = sceBatches, useAssay = "logcounts",
                           nComponents = 10, verbose = FALSE),
    "recommended")
  testthat::expect_true("HARMONY" %in% reducedDimNames(sceAssay))
  hAssay <- SingleCellExperiment::reducedDim(sceAssay, "HARMONY")
  testthat::expect_equal(nrow(hAssay), ncol(sceBatches))
  testthat::expect_equal(ncol(hAssay), 10)
  testthat::expect_true(all(is.finite(hAssay)))

  # same seed gives the same result
  sceRD2 <- runHarmony(inSCE = sceBatches, useReducedDim = "PCA",
                       nComponents = 10, verbose = FALSE)
  testthat::expect_equal(
    hRD, SingleCellExperiment::reducedDim(sceRD2, "HARMONY"))
})

if (isTRUE(reticulate::py_available(initialize = FALSE))) {
  if(reticulate::py_module_available(module = "bbknn")){
    test_that(desc = "Testing BBKNN", {
      sceBatches <- runBBKNN(inSCE = sceBatches)
      testthat::expect_true("BBKNN" %in% reducedDimNames(sceBatches))
    })
  }
  if(reticulate::py_module_available(module = "scanorama")){
    test_that(desc = "Testing SCANORAMA", {
      sceBatches <- runSCANORAMA(inSCE = sceBatches)
      testthat::expect_true("SCANORAMA" %in% assayNames(sceBatches))
    })
  }
}
