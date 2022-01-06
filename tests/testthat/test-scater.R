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
  sce <- scaterPCA(sce, useAssay = "logNormCounts")

  testthat::expect_true("logNormCounts" %in% assayNames(sce))
  testthat::expect_true("PCA" %in% reducedDimNames(sce))
})
