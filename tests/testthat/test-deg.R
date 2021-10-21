# DEG Functions
library(singleCellTK)
context("Testing DEG functions")
data(sceBatches, package = "singleCellTK")

sceBatches <- scaterlogNormCounts(sceBatches, "logcounts")
sceBatches <- subsetSCECols(sceBatches, colData = "batch == 'w'")

test_that(desc = "Testing Limma DE", {
  sceBatches <- runLimmaDE(inSCE = sceBatches,
                           class = "cell_type",
                           classGroup1 = "alpha", classGroup2 = "beta",
                           groupName1 = "a", groupName2 = "b",
                           analysisName = "aVSbLimma")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbLimma" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing MAST DE", {
  sceBatches <- runMAST(inSCE = sceBatches,
                        class = "cell_type",
                        classGroup1 = "alpha", classGroup2 = "beta",
                        groupName1 = "a", groupName2 = "b",
                        analysisName = "aVSbMAST")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbMAST" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing MAST DE", {
  sceBatches <- runDESeq2(inSCE = sceBatches,
                          class = "cell_type",
                          classGroup1 = "alpha", classGroup2 = "beta",
                          groupName1 = "a", groupName2 = "b",
                          analysisName = "aVSbDESeq2")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbDESeq2" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing MAST DE", {
  sceBatches <- runANOVA(inSCE = sceBatches,
                         class = "cell_type",
                         classGroup1 = "alpha", classGroup2 = "beta",
                         groupName1 = "a", groupName2 = "b",
                         analysisName = "aVSbANOVA")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbANOVA" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing MAST DE", {
  sceBatches <- runWilcox(inSCE = sceBatches,
                          class = "cell_type",
                          classGroup1 = "alpha", classGroup2 = "beta",
                          groupName1 = "a", groupName2 = "b",
                          analysisName = "aVSbWilcox")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbWilcox" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing findMarker", {
  sceBatches <- findMarkerDiffExp(inSCE = sceBatches,
                                  cluster = "cell_type")

  testthat::expect_true("findMarker" %in% names(metadata(sceBatches)))
  testthat::expect_true(nrow(metadata(sceBatches)$findMarker) > 0)
})
