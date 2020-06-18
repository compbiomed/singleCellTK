context("misc functions")

data(scExample, package = "singleCellTK")
sce <- sce[, colData(sce)$type != 'EmptyDroplet']

test_that("summarizeSCE", {
  ta <- summarizeSCE(sce, sample = NULL)
  expect_is(ta, "data.frame")
})
