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

test_that(desc = "Testing DESeq2 DE", {
  sceBatches <- runDESeq2(inSCE = sceBatches,
                          class = "cell_type",
                          classGroup1 = "alpha", classGroup2 = "beta",
                          groupName1 = "a", groupName2 = "b",
                          analysisName = "aVSbDESeq2")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbDESeq2" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing ANOVA DE", {
  sceBatches <- runANOVA(inSCE = sceBatches,
                         class = "cell_type",
                         classGroup1 = "alpha", classGroup2 = "beta",
                         groupName1 = "a", groupName2 = "b",
                         analysisName = "aVSbANOVA")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbANOVA" %in% names(metadata(sceBatches)$diffExp))
})

test_that(desc = "Testing Wilcoxon DE", {
  sceBatches <- runWilcox(inSCE = sceBatches,
                          class = "cell_type",
                          classGroup1 = "alpha", classGroup2 = "beta",
                          groupName1 = "a", groupName2 = "b",
                          analysisName = "aVSbWilcox")

  testthat::expect_true("diffExp" %in% names(metadata(sceBatches)))
  testthat::expect_true("aVSbWilcox" %in% names(metadata(sceBatches)$diffExp))

  # Also Plotting functions at this point
  vlcn <- plotDEGVolcano(sceBatches, "aVSbWilcox")
  testthat::expect_is(vlcn, "ggplot")
  
  hm <- plotDEGHeatmap(sceBatches, "aVSbWilcox", 
                       minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL)
  testthat::expect_is(hm, "ggplot")

  pR <- plotDEGRegression(sceBatches, "aVSbWilcox")
  testthat::expect_is(pR, "ggplot")

  pV <- plotDEGViolin(sceBatches, "aVSbWilcox")
  testthat::expect_is(pV, "ggplot")
})

test_that(desc = "Testing findMarker", {
  sceBatches <- findMarkerDiffExp(inSCE = sceBatches,
                                  cluster = "cell_type")

  testthat::expect_true("findMarker" %in% names(metadata(sceBatches)))

  topTable <- findMarkerTopTable(sceBatches)
  testthat::expect_is(topTable, "data.frame")
  testthat::expect_named(topTable, c("Gene", "Log2_FC", "Pvalue", "FDR",
                                     "cell_type", "clusterExprPerc",
                                     "ControlExprPerc", "clusterAveExpr"))
  testthat::expect_gt(nrow(topTable), 0)

  hmFM <- plotMarkerDiffExp(sceBatches)
  testthat::expect_is(hmFM, "ggplot")
})
