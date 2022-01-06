# computeZScore.R
library(singleCellTK)
context("Testing computeZScore.R")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing computeZScore.R", {
  scaledMatrix <- computeZScore(assay(sce, "counts"))
  
  testthat::expect_equal(nrow(scaledMatrix), nrow(sce))
  testthat::expect_equal(ncol(scaledMatrix), ncol(sce))
})
