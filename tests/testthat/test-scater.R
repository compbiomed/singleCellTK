# Scater Functions
library(singleCellTK)
context("Testing scater functions")
data(scExample, package = "singleCellTK")
# Remove zero colsums cells required for scater functions
zeroCols <- which(colSums(assay(sce, "counts")) == 0)
sce <- sce[, -zeroCols]

test_that(desc = "Testing scaterCPM", {
  sce <- scaterCPM(sce)

  testthat::expect_true("ScaterCPMCounts" %in% assayNames(sce))
})

test_that(desc = "Testing scaterLogNormCounts & scaterPCA", {
  sce <- scaterlogNormCounts(sce, useAssay = "counts", assayName = "logNormCounts")
  sce <- scaterPCA(sce, useAssay = "logNormCounts", useFeatureSubset = NULL)

  testthat::expect_true("logNormCounts" %in% assayNames(sce))
  testthat::expect_true("PCA" %in% reducedDimNames(sce))
})

test_that(desc = "scaterlogNormCounts matches library-size normalization", {
  expect_no_warning(
    res <- scaterlogNormCounts(sce, useAssay = "counts",
                               assayName = "logcounts"),
    class = "deprecatedWarning"
  )
  counts <- as.matrix(assay(sce, "counts"))
  libSize <- colSums(counts)
  sf <- libSize / mean(libSize)
  expected <- log2(t(t(counts) / sf) + 1)
  expect_equal(as.matrix(assay(res, "logcounts")), expected)
  expect_equal(sizeFactors(res), stats::setNames(sf, colnames(sce)))
  expect_s4_class(assay(res, "logcounts"), "dgCMatrix")
})

test_that(desc = "scaterlogNormCounts centers existing size factors", {
  sceSF <- sce
  sizeFactors(sceSF) <- seq(1, 3, length.out = ncol(sceSF))
  res <- suppressWarnings(scaterlogNormCounts(sceSF, assayName = "logcounts"))
  sf <- sizeFactors(sceSF) / mean(sizeFactors(sceSF))
  counts <- as.matrix(assay(sce, "counts"))
  expect_equal(as.matrix(assay(res, "logcounts")),
               log2(t(t(counts) / sf) + 1))
})

test_that(desc = "scaterlogNormCounts stops on cells without counts", {
  sceZero <- sce
  counts(sceZero)[, 1] <- 0
  expect_error(suppressWarnings(scaterlogNormCounts(sceZero)),
               "size factors should be positive")
})
