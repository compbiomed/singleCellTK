# DEG Functions
library(singleCellTK)
context("Testing DEG functions")
data(sceBatches, package = "singleCellTK")

logcounts(sceBatches) <- log1p(counts(sceBatches))
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
  testthat::expect_is(hm, "Heatmap")

  pR <- plotDEGRegression(sceBatches, "aVSbWilcox")
  testthat::expect_is(pR, "ggplot")

  pV <- plotDEGViolin(sceBatches, "aVSbWilcox")
  testthat::expect_is(pV, "ggplot")
})

test_that(desc = "Testing findMarker", {
  sceBatches <- runFindMarker(inSCE = sceBatches,
                              cluster = "cell_type")
  testthat::expect_true("findMarker" %in% names(metadata(sceBatches)))

  topTable <- getFindMarkerTopTable(sceBatches, log2fcThreshold = 1,
                                    fdrThreshold = 0.05, minClustExprPerc = 0.7,
                                    maxCtrlExprPerc = 0.4, minMeanExpr = 1,
                                    topN = 10)
  testthat::expect_is(topTable, "data.frame")
  testthat::expect_named(topTable, c("Gene", "Log2_FC", "Pvalue", "FDR",
                                     "cell_type", "clusterExprPerc",
                                     "ControlExprPerc", "clusterAveExpr"))
  testthat::expect_gt(nrow(topTable), 0)

  hmFM <- plotFindMarkerHeatmap(sceBatches)
  testthat::expect_is(hmFM, "Heatmap")
})

test_that(desc = "Internal Wilcoxon test matches stats::wilcox.test", {
  set.seed(1)
  mat <- matrix(rpois(30 * 25, 2), nrow = 30,
                dimnames = list(paste0("g", 1:30), NULL))
  mat[1, ] <- 3
  mat[2, ] <- c(rep(0, 20), 1:5)
  ix1 <- 1:12
  ix2 <- 13:25
  res <- .wilcoxTest(mat, ix1, ix2, chunkSize = 7)
  expected <- apply(mat, 1, function(v) {
    p <- suppressWarnings(stats::wilcox.test(v[ix1], v[ix2], exact = FALSE,
                                             correct = TRUE)$p.value)
    if (is.na(p)) 1 else p
  })
  expect_equal(res$p.value, unname(expected))
  expect_equal(res$FDR, stats::p.adjust(unname(expected), method = "BH"))
  expect_equal(rownames(res), rownames(mat))
})

test_that(desc = "runWilcox runs without deprecations", {
  expect_no_warning(
    runWilcox(inSCE = sceBatches, class = "cell_type",
              classGroup1 = "alpha", classGroup2 = "beta",
              groupName1 = "a", groupName2 = "b",
              analysisName = "aVSbWilcox2"),
    class = "deprecatedWarning"
  )
})
